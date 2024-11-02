from pathlib import Path
from typing import Dict, List, Literal, Tuple, Union

import torch
from torch.utils.data import IterableDataset


class FastaSamplerDataset(IterableDataset):
    def __init__(
        self,
        range_mean: float,
        range_std: float,
        subsequence_range_mean: float,
        subsequence_range_std: float,
        fasta_file: Path,
        sampling_strategy: Literal["local", "subsequence", "random_subsequence"] = "random_subsequence",
    ):
        super().__init__()
        self.range_mean = range_mean
        self.range_std = range_std
        self.fasta_file = fasta_file
        self.subsequence_range_mean = subsequence_range_mean
        self.subsequence_range_std = subsequence_range_std

        # load in text file
        with open(self.fasta_file, "r") as f:
            self.text = f.read()

        self.len_text = len(self.text)
        self.sampling_strategy = sampling_strategy

    def iter_local_sequence(self):
        """
        Randomly samples two sequence from the fasta file which constitute a positive
        sample.

        1) sample a random length L_1, and random sequence index i_1 as well as a direction [left, right]
        2) sample sequence x_1 from the range [i_1, i_1 +/- L_1]
        2) From the range of the first sequence sample a index
        3) the sample a direction (left or right) for the second sequence
        4) Sample a random length for the second sequence
        5) then sample the second sequence from the range of the first sequence
        """

        while True:
            i_1 = torch.randint(0, self.len_text, (1,))
            L_1 = torch.normal(self.range_mean, self.range_std, (1,)).int()
            # sample a direction
            direction_1 = torch.randint(0, 2, (1,))
            # sample a sequence
            if direction_1 == 0:
                x_1 = self.text[i_1 : i_1 + L_1]
            else:
                x_1 = self.text[i_1 - L_1 : i_1]

            # sample a second index
            i_2 = torch.randint(0, self.len_text, (1,))
            direction_2 = torch.randint(0, 2, (1,))
            L_2 = torch.normal(self.range_mean, self.range_std, (1,)).int()
            # sample a sequence
            if direction_2 == 0:
                x_2 = self.text[i_2 : i_2 + L_2]
            else:
                x_2 = self.text[i_2 - L_2 : i_2]

            yield x_1, x_2

    def iter_subsequence(self):
        """
        Randomly sampled a sequence from the fasta file and then samples and then samples a subsequence from that sequence
        """

        while True:
            i_1 = torch.randint(0, self.len_text, (1,))
            L_1 = torch.normal(self.range_mean, self.range_std, (1,)).int()
            direction = torch.randint(0, 2, (1,))
            if direction == 0:
                x_1 = self.text[i_1 : i_1 + L_1]
            else:
                x_1 = self.text[i_1 - L_1 : i_1]

            # sample a second sequence from the first sequence
            i_2 = torch.randint(0, int(L_1), (1,))
            L_2 = torch.normal(self.range_mean, self.range_std, (1,)).int()
            # choose the direction which leads to the longest sequence
            if i_2 + L_2 > L_1:
                direction = 0
            else:
                direction = 1
            if direction == 0:
                x_2 = x_1[i_2 : i_2 + L_2]
            else:
                x_2 = x_1[i_2 - L_2 : i_2]

            yield x_1, x_2



    def iter_random_subsequence(self, subsequence_mean_length: int=200, subsequence_std_length: int=20):
        """
        differs from iter_subsequence in that the second sequence does not try to capture the longest sequence
        but rather samples a random subsequence from the first sequence
        """

        while True:
            L_1 = torch.normal(self.range_mean, self.range_std, (1,)).int()
            i_1 = torch.randint(0, int(self.len_text-L_1), (1,))

            # sample lengh of second sequence
            L_2 = torch.normal(subsequence_mean_length, subsequence_std_length, (1,)).int()
            # sample the start of the second sequence from the first sequence [0, L_1 - L_2]
            i_2 = torch.randint(0, int(L_1 - L_2), (1,))

            x_1 = self.text[i_1 : i_1 + L_1]
            x_2 = x_1[i_2 : i_2 + L_2]
            yield x_1, x_2
            



    def __iter__(self):
        if self.sampling_strategy == "local":
            return self.iter_local_sequence()
        elif self.sampling_strategy == "subsequence":
            return self.iter_subsequence()
        elif self.sampling_strategy == "random_subsequence":
            return self.iter_random_subsequence(self.subsequence_range_mean, self.subsequence_range_std)
        else:
            raise ValueError(
                f"Sampling strategy {self.sampling_strategy} not implemented"
            )


class FastaUniformSampler(IterableDataset):
    def __init__(
        self,
        range_min: int,
        range_max: int,
        subsequence_range_min: int,
        subsequence_range_max: int,
        fasta_file: Union[Path, List[Path]],
        sampling_strategy: Literal["random_subsequence", "random_subsequence_uppercase"] = "random_subsequence",
    ):
        super().__init__()
        self.range_min = range_min
        self.range_max = range_max
        self.subsequence_range_min = subsequence_range_min
        self.subsequence_range_max = subsequence_range_max

        if isinstance(fasta_file, Path):
            fasta_file = [fasta_file]
        self.fasta_file = fasta_file

        # load in text file and concatenate
        files = []
        self.text = ""
        for path in self.fasta_file:
            with open(path, "r") as f:
                self.text += f.read()

        self.len_text = len(self.text)
        self.sampling_strategy = sampling_strategy


    def iter_random_subsequence(self):
        while True:
            L_1 = torch.randint(self.range_min, self.range_max, (1,)).int()
            i_1 = torch.randint(low=0, high=int(self.len_text) - int(L_1), size= (1,))
            x_1 = self.text[i_1 : i_1 + L_1]

            # sample lengh of second sequence
            L_2 = torch.randint(self.subsequence_range_min, self.subsequence_range_max, (1,)).int()
            # sample the start of the second sequence from the first sequence [0, L_1 - L_2]
            i_2 = torch.randint(0, int(L_1 - L_2), (1,))
            x_2 = x_1[i_2 : i_2 + L_2]
            yield x_1, x_2

    def iter_random_subsequence_uppercase(self):
        for x_1, x_2 in self.iter_random_subsequence():
            yield x_1.upper(), x_2.upper()

    def __iter__(self):
        if self.sampling_strategy == "random_subsequence":
            return self.iter_random_subsequence()
        elif self.sampling_strategy == "random_subsequence_uppercase":
            return self.iter_random_subsequence_uppercase()
        else:
            raise ValueError(
                f"Sampling strategy {self.sampling_strategy} not implemented"
            )
        



def collate_fn(
    batch, tokenizer
) -> Tuple[Dict[str, torch.Tensor], Dict[str, torch.Tensor]]:
    """
    collate to max batch size and output a dictionary with two elements
    ids = matrix of shape (batch_size, max_sequence_length)
    attention_mask = matrix of shape (batch_size, max_sequence_length)
    """
    x_1, x_2 = list(zip(*batch))
    x_1 = tokenizer.tokenize(x_1)
    x_2 = tokenizer.tokenize(x_2)

    return x_1.to_torch(), x_2.to_torch()
