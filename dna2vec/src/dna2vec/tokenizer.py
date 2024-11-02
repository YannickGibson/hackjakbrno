"""
DNA Tokenizer template
- MetaSpace enforces preprocessing of substrings prior to sub-token detection.
https://colab.research.google.com/github/huggingface/notebooks/blob/master/examples/tokenizer_training.ipynb
- 

"""
from abc import ABC, abstractmethod
from typing import Dict, List

import torch
from pydantic import BaseModel
from tokenizers import Encoding, Tokenizer, decoders, models, pre_tokenizers, processors
from tokenizers.trainers import BpeTrainer


class TokenizationOutput(BaseModel):
    encodings: List[Encoding]

    class Config:
        arbitrary_types_allowed = True

    def to_torch(self) -> Dict[str, torch.Tensor]:
        return {
            "input_ids": torch.tensor([x.ids for x in self.encodings]),
            "attention_mask": torch.tensor([x.attention_mask for x in self.encodings]),
        }


class AbstractTokenizer(ABC):
    @abstractmethod
    def train_from_generator(self, generator: List[str]):
        pass

    @abstractmethod
    def save(self, filename: str):
        pass

    @abstractmethod
    def load(self, filename: str):
        pass

    @abstractmethod
    def tokenize(self, sequence: str) -> List[str]:
        pass


class BPTokenizer(AbstractTokenizer):
    def __init__(self, vocab_size: int = 20, min_frequency: int = 2):
        self.tokenizer = Tokenizer(models.BPE())
        self.tokenizer.pre_tokenizer = pre_tokenizers.Metaspace()
        self.tokenizer.decoder = decoders.Metaspace()

        self.tokenizer.post_processor = processors.TemplateProcessing(
            single="[CLS] $A [SEP]",
            special_tokens=[
                ("[CLS]", vocab_size + 1),
                ("[SEP]", vocab_size + 2),
            ],
        )
        self.min_frequency = min_frequency
        self.tokenizer.enable_padding(pad_id=vocab_size + 3, pad_token="[PAD]")
        self.max_length = 1024
        ######
        #self.tokenizer.enable_truncation(max_length = 1024)

    @property
    def vocab_size(self):
        return self.tokenizer.get_vocab_size() + 3  # +3 for [PAD], [CLS], [SEP]

    def train_from_generator(self, generator: List[str]):
        self.tokenizer.train_from_iterator(
            generator,
            trainer=BpeTrainer(
                min_frequency=self.min_frequency, vocab_size=self.vocab_size
            ),
        )

    def save(self, filename: str):
        self.tokenizer.save(filename)

    @classmethod
    def load(cls, filename: str) -> "BPTokenizer":
        print(filename)
        tokenizer = cls()
        tokenizer.tokenizer = Tokenizer.from_file(filename)
        cls.min_frequency = None
        return tokenizer

    def tokenize(self, sequences: List[str]) -> TokenizationOutput:
        #self.tokenizer.enable_truncation(max_length=self.max_length)
        
        # Encode the sequences with truncation applied
        encoded = self.tokenizer.encode_batch(sequences)
        return TokenizationOutput(encodings=encoded)


def generate_sequence(file_path, sequence_length, max_sequences=None):
    with open(file_path, "r") as file:
        buffer = ""
        num_sequences = 0
        for line in file:
            buffer += line.strip()
            while len(buffer) >= sequence_length:
                yield buffer[:sequence_length]
                buffer = buffer[sequence_length:]
                num_sequences += 1
                if max_sequences is not None and num_sequences >= max_sequences:
                    return


if __name__ == "__main__":
    # Create an instance of the CustomTokenizer
    tokenizer = BPTokenizer()

    # Train the tokenizer from a generator
    # The generator needs to be sufficient.
    generator = generate_sequence("NC_000002.12.txt", 1000, 500)
    tokenizer.train_from_generator(generator)

    # Save the tokenizer to a file
    tokenizer.save("dna_tokenizer_20.json")
