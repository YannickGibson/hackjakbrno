import torch
from transformers import AutoTokenizer, AutoModel
from transformers.models.bert.configuration_bert import BertConfig

class DNABERT2:
    def __init__(self):
        self.device = "cuda:0" if torch.cuda.is_available() else "cpu"
        self.tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
        self.config = BertConfig.from_pretrained("zhihan1996/DNABERT-2-117M")
        self.model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True, config=self.config)
        self.model.to(self.device)
    def get_embedding(self, dna_sequence, method="mean"):
        """
            Returns embedded sequences in shape [batch, embedding_size]
            Padding is broken, better results if batch_size=1
        """
        #check if its a batch or a single sequence            
        tokens = self.tokenizer(dna_sequence, return_tensors = 'pt', padding=True)["input_ids"].to(self.device)
        embeddings = self.model(tokens)[0]
        if method == "mean":
            return torch.mean(embeddings, dim=1)
        elif method == "max":
            return torch.max(embeddings, dim=1)
        elif method == "min":
            return torch.min(embeddings, dim=1)
        else:
            raise ValueError(f"Embedding method is invalid ({method})")
