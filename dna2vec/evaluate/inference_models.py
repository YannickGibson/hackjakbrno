import torch
from typing import Literal
from transformers import PreTrainedModel


class EvalModel():
    def __init__(self, 
                 tokenizer, 
                 model, 
                 pooling, 
                 device):
        
        self.tokenizer = tokenizer
        #self.tokenizer.enable_truncation(max_length = 1024)
        self.model = model
        
        self.pooling = pooling
        self.pooling = self.pooling.to(device)
        
        self.model = self.model.to(device)
        self.device = device
        
        self.model.eval()
    
    def encode(self, x: list):
        #print("test")
        with torch.no_grad():
            input_data = self.tokenizer.tokenize(x).to_torch()
            #print(type(input_data))
            #print(input_data["input_ids"].shape)
            #print(input_data["attention_mask"])
            last_hidden_state = self.model(input_ids = input_data["input_ids"].to(self.device), 
                                           attention_mask = input_data["attention_mask"].to(self.device))
            #print(last_hidden_state)
            y = self.pooling(last_hidden_state, attention_mask=input_data["attention_mask"].to(self.device))
            return torch.nn.functional.normalize(y.squeeze(), dim=0).detach().cpu().numpy()



class Baseline():
    def __init__(self, 
                 option: Literal["dna2bert-max", 
                                 "dna2bert-mean", 
                                 "nucleotide-transformer-mean", 
                                 "nucleotide-transformer-max",
                                 "hyena-dna-mean",
                                 "hyena-dna-max"], 
                 device: str,
                 cache_dir: str = '/mnt/SSD2/pholur/dna2vec/baselines'):
        
        import os
        os.environ['TRANSFORMERS_CACHE'] = cache_dir
        self.device = device
        from transformers import AutoTokenizer, AutoModelForMaskedLM, AutoModel
        
        if option == "dna2bert-mean" or option == "dna2bert-max":
            self.tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
            self.model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
            self.model.to(device)
            self.model.eval()

        elif option == "nucleotide-transformer-max" or option == "nucleotide-transformer-mean":
            self.tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-500m-human-ref")
            self.model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-500m-human-ref")
            self.model.to(device)
            self.model.eval()
        
        elif option == "hyena-dna-mean" or option == "hyena-dna-max":
            
            from inference_model_helpers import CharacterTokenizer, HyenaDNAPreTrainedModel
            self.model = HyenaDNAPreTrainedModel.from_pretrained(
                './checkpoints',
                "hyenadna-small-32k-seqlen",
                download=True,
                config=None,
                device=device,
                use_head=False,
                n_classes=2,
            )
            
            max_length = 32_000
            self.tokenizer = CharacterTokenizer(
                            characters=['A', 'C', 'G', 'T', 'N'],  # add DNA characters, N is uncertain
                            model_max_length=max_length + 2,  # to account for special tokens, like EOS
                            add_special_tokens=False,  # we handle special tokens elsewhere
                            padding_side='left', # since HyenaDNA is causal, we pad on the left
                        )
            self.model.to(device)
            self.model.eval()
        
        else:
            raise NotImplementedError("Baseline is not implemented.")
        
        self.option = option
        
        
    def encode(self, x: list):
        with torch.no_grad():
            if self.option == "nucleotide-transformer-mean":
                tokens_ids = self.tokenizer.batch_encode_plus(x, return_tensors="pt", padding=True, truncation=True)["input_ids"]
                # Compute the embeddings
                attention_mask = tokens_ids != self.tokenizer.pad_token_id
                torch_outs = self.model(
                    tokens_ids.to(self.device),
                    attention_mask=attention_mask.to(self.device),
                    encoder_attention_mask=attention_mask.to(self.device),
                    output_hidden_states=True
                )
                
                embeddings = torch_outs['hidden_states'][-1]
                mean_sequence_embeddings = torch.sum(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)/\
                        torch.sum(attention_mask.to(self.device), axis=-1).unsqueeze(-1)
                return torch.nn.functional.normalize(mean_sequence_embeddings.squeeze(), dim=0).detach().cpu().numpy()
            
            
            elif self.option == "nucleotide-transformer-max":
                tokens_ids = self.tokenizer.batch_encode_plus(x, return_tensors="pt", padding=True, truncation=True)["input_ids"]
                # Compute the embeddings
                attention_mask = tokens_ids != self.tokenizer.pad_token_id
                torch_outs = self.model(
                    tokens_ids.to(self.device),
                    attention_mask=attention_mask.to(self.device),
                    encoder_attention_mask=attention_mask.to(self.device),
                    output_hidden_states=True
                )
                
                embeddings = torch_outs['hidden_states'][-1]
                max_sequence_embeddings = torch.max(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)
                return torch.nn.functional.normalize(max_sequence_embeddings.values.squeeze(), dim=0).detach().cpu().numpy()
            
            
            elif self.option == "dna2bert-mean":
                tokens_ids = self.tokenizer(x, return_tensors = 'pt', padding=True)["input_ids"]
                attention_mask = tokens_ids != self.tokenizer.pad_token_id
                
                torch_outs = self.model(
                    tokens_ids.to(self.device),
                    attention_mask=attention_mask.to(self.device),
                    encoder_attention_mask=attention_mask.to(self.device),
                    output_hidden_states=True,
                )

                embeddings = torch_outs[0]
                mean_sequence_embeddings = torch.sum(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)/\
                        torch.sum(attention_mask.to(self.device), axis=-1).unsqueeze(-1)
                return torch.nn.functional.normalize(mean_sequence_embeddings.squeeze(), dim=0).detach().cpu().numpy()
            
            elif self.option == "dna2bert-max":
                tokens_ids = self.tokenizer(x, return_tensors = 'pt', padding=True)["input_ids"]
                attention_mask = tokens_ids != self.tokenizer.pad_token_id
                
                torch_outs = self.model(
                    tokens_ids.to(self.device),
                    attention_mask=attention_mask.to(self.device),
                    encoder_attention_mask=attention_mask.to(self.device),
                    output_hidden_states=True,
                )

                embeddings = torch_outs[0]
                max_sequence_embeddings = torch.max(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)
                return torch.nn.functional.normalize(max_sequence_embeddings.values.squeeze(), dim=0).detach().cpu().numpy()
            
            elif self.option == "hyena-dna-mean":
                tok_seq = self.tokenizer(x)
                tok_seq = tok_seq["input_ids"]  # grab ids
                
                attention_mask = torch.zeros((len(tok_seq), 1252)) + 1
                for i,s in enumerate(tok_seq):
                    if len(s) < 1252: # VERY MUCH AN EDGE CASE - 1 in a million chance. Due to chromosome cuts.
                        attention_mask[i][len(s):1252] = 0
                        s.extend([4]*(1252 - len(s)))

                # place on device, convert to tensor
                tok_seq = torch.LongTensor(tok_seq) #.unsqueeze(0)  # unsqueeze for batch dim
                tok_seq = tok_seq.to(self.device)

                # prep model and forward
                with torch.inference_mode():
                    embeddings = self.model(tok_seq)
                
                mean_sequence_embeddings = torch.sum(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)/\
                        torch.sum(attention_mask.to(self.device), axis=-1).unsqueeze(-1)
                return torch.nn.functional.normalize(mean_sequence_embeddings.squeeze(), dim=0).detach().cpu().numpy()
                
            elif self.option == "hyena-dna-max":
                tok_seq = self.tokenizer(x)
                tok_seq = tok_seq["input_ids"]  # grab ids

                attention_mask = torch.zeros((len(tok_seq), 1252)) + 1
                for i,s in enumerate(tok_seq):
                    if len(s) < 1252: # VERY MUCH AN EDGE CASE - 1 in a million chance. Due to chromosome cuts.
                        attention_mask[i][len(s):1252] = 0
                        s.extend([4]*(1252 - len(s)))

                # place on device, convert to tensor
                tok_seq = torch.LongTensor(tok_seq) #.unsqueeze(0)  # unsqueeze for batch dim
                tok_seq = tok_seq.to(self.device)

                # prep model and forward
                with torch.inference_mode():
                    embeddings = self.model(tok_seq)
                
                max_sequence_embeddings = torch.max(
                    attention_mask.to(self.device).unsqueeze(-1)*embeddings, axis=-2)
                return torch.nn.functional.normalize(max_sequence_embeddings.values.squeeze(), dim=0).detach().cpu().numpy()


            else:
                raise NotImplementedError("Baseline Option undefined")


    def get_sentence_embedding_dimension(self):
        
        if self.option == "dna2bert-max" or self.option == "dna2bert-mean":
            return 768
        elif self.option == "nucleotide-transformer-max" or self.option == "nucleotide-transformer-mean":
            return 1280
        elif self.option == "hyena-dna-max" or self.option == "hyena-dna-mean":
            return 256
        else:
            raise NotImplementedError("Dimension unspecified")