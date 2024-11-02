import sys
sys.path.append("../src/")
sys.append("../evaluate/")

import torch
from tqdm import tqdm

from pinecone_store import EvalModel
from helpers import pick_random_lines, pick_from_special_gene_list, pick_from_chimp_2a_2b, pick_from_chromosome3

checkpoints = {}

checkpoints["trained"] = {}
checkpoints["trained"]["model"] = torch.load("/mnt/SSD5/pholur/checkpoints/checkpoint_smooth-rock-29.pt")

checkpoints["init"] = {}
checkpoints["init"]["model"] = torch.load("/mnt/SSD5/pholur/checkpoints/checkpoint_initalialized.pt")

print("Checkpoints loaded.")

from dna2vec.trainer import ContrastiveTrainer
from dna2vec.model import model_from_config
# The path for the tokenizer isn't relative.

embedders = {}
device = "cuda:0"


for baseline in checkpoints:
    config = checkpoints[baseline]["model"]["config"]
    config.model_config.tokenizer_path = "/home/pholur/dna2vec/src/model/tokenizers/dna_tokenizer_10k.json"
    encoder, pooling, tokenizer = model_from_config(config.model_config)
    
    encoder.load_state_dict(checkpoints[baseline]["model"]["model"])
    encoder.eval()
    
    embedders[baseline] = EvalModel(
        tokenizer, encoder, pooling, device
    )



def encode(sequences: list):
    import numpy as np
    train_mat = np.zeros((len(sequences), 384))
    init_mat = np.zeros((len(sequences), 384))
    labels = []
    for i, sequence in tqdm(enumerate(sequences)):
        for baseline, model in embedders.items():
            if baseline == "trained":
                train_mat[i,:] = model.encode([sequence[0]])
            elif baseline == "init":
                init_mat[i,:] = model.encode([sequence[0]])
            else:
                raise NotImplementedError("Baseline not defined.")
        labels.append(sequence[1])
    return train_mat, init_mat, labels



def main(samples: int, mode: str, sequences_prior: int = -1):
    import numpy as np
    
    if mode == "train + random":
        sequences = pick_random_lines(
            path = "/home/pholur/dna2vec/tests/data/subsequences_sample_train.txt",
            k = samples,
        )
        
    elif mode == "train + sequenced":
        sequences = pick_random_lines(
            path = "/home/pholur/dna2vec/tests/data/subsequences_sample_train.txt",
            k = samples,
            mode = "sequenced",
        )
        
    elif mode == "test + random":
        sequences = pick_random_lines(
            path = "/home/pholur/dna2vec/tests/data/subsequences_sample_test.txt",
            k = samples,
        )
        
    elif mode == "train + subsequenced":
        sequences = pick_random_lines(
            path = "/home/pholur/dna2vec/tests/data/subsequences_sample_train.txt",
            k = samples,
            mode = "subsequenced",
            sequences_prior=sequences_prior,
        )
    
    elif mode == "select genes":
        sequences = pick_from_special_gene_list(
            gene_path = "/home/pholur/dna2vec/tests/data/ch2_genes.csv",
            full_path = "/home/pholur/dna2vec/tests/data/NC_000002.12.txt",
            samples = samples
        )
    
    elif mode == "select chimp 2a 2b":
        sequences = pick_from_chimp_2a_2b(
            path_2a = "/home/pholur/dna2vec/tests/data/chimp2a2b/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa",
            path_2b = "/home/pholur/dna2vec/tests/data/chimp2a2b/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa",
            samples = samples
        )
    
    elif mode == "from chromosome 3":
        sequences = pick_from_chromosome3(
            path = "/home/pholur/dna2vec/tests/data/NC_000003.fasta",
            samples = samples,
        )
    
    else:
        raise NotImplementedError("Mode is not defined.")
    
    mat_train, mat_init, labels = encode(sequences)
    np.savez_compressed(f"test_cache/clustering/cluster_{mode}_{samples}_{sequences_prior}.npz", train = mat_train, init = mat_init, labels = labels)






if __name__ == "__main__":
    # main(samples=2000, mode="train + random")
    # main(samples=1000, mode="test + random")
    # main(samples=5000, mode="train + sequenced")
    # main(samples=5000, mode="train + subsequenced", sequences_prior=5)
    # main(samples=5000, mode='select genes')
    # main(samples=6000, mode='select chimp 2a 2b')
    main(samples=6000, mode='from chromosome 3')


