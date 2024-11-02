import random
import yaml
import numpy as np

from pinecone_store import PineconeStore
from aligners.smith_waterman import bwamem_align, bwamem_align_parallel

from alignment_metrics import calculate_smith_waterman_distance
from collections import defaultdict

# import sys
# sys.path.append("../src/")
from dna2vec.model import model_from_config

from collections import defaultdict
import time


with open("configs/data_recipes.yaml", 'r') as stream:
    data_recipes = yaml.safe_load(stream)
    
with open("configs/model_checkpoints.yaml", 'r') as stream:
    checkpoints = yaml.safe_load(stream)
    
with open("configs/raw.yaml", 'r') as stream:
    raw_fasta_files = yaml.safe_load(stream)

from tqdm import tqdm





import scipy.stats as stats

def clopper_pearson_interval(successes, trials, confidence_level=0.95):
    alpha = 1 - confidence_level
    
    lower_bound = stats.beta.ppf(alpha / 2, successes, trials - successes + 1)
    upper_bound = stats.beta.ppf(1 - alpha / 2, successes + 1, trials - successes)
    
    return lower_bound, upper_bound





import time
# def main_align(store, 
#                 queries, 
#                 indices, 
#                 top_k,
#                 exactness=0,
#                 distance_bound=0,
#                 flex=False,
#                 batch_size=64):
    
    
#     finer_flag = np.zeros((len(queries), 1))
#     returned = store.query_batch(queries, indices, top_k=top_k)
    
#     for i,returned_unit in enumerate(returned): # can potentially be shuffled but that's okay
#         returned_unit_matched = returned_unit["matches"]
#         trained_positions = [sample["metadata"]["position"] for sample in returned_unit_matched]
#         metadata_set = [sample["metadata"]["metadata"] for sample in returned_unit_matched]
        
#         all_candidate_strings = [sample["metadata"]["text"] for sample in returned_unit_matched]
        
#         identified_sub_indices, identified_indices, _, timer, smallest_distance = bwamem_align_parallel(
#                                                                                             all_candidate_strings, 
#                                                                                             trained_positions, 
#                                                                                             metadata_set,
#                                                                                             returned_unit["query"])
        
#         series = [int(tup[0]) + int(tup[1]) for tup in zip(identified_sub_indices, identified_indices)]

#         if flex:
#             if is_within_range_of_any_element(returned_unit["index"], series, exactness) or \
#                     abs(smallest_distance + 2*len(returned_unit["query"])) < distance_bound + 1: # the 1 here helps with instabilities
#                 finer_flag[i,0] = 1
#             else:
#                 finer_flag[i,0] = 0
#         else:
#             if (returned_unit["index"] in series) or abs(smallest_distance + 2*len(returned_unit["query"])) < 1: # exact SW match
#                 finer_flag[i,0] = 1
#             else:
#                 finer_flag[i,0] = 0

#     return finer_flag



def main_align(store, 
                queries, 
                indices, 
                top_k,
                exactness=0,
                distance_bound=0,
                flex=False,
                batch_size=64,
                match=True,
                distributed=False,
                per_k=0,
                namespaces=None,
                namespace_dict=None,
                dictionary_of_values=None):

    if not distributed:
        
        num_queries = len(queries)
        finer_flag = np.zeros((num_queries, 1))

        for batch_start in tqdm(range(0, num_queries, batch_size)):
            batch_end = min(batch_start + batch_size, num_queries)
            batch_queries = queries[batch_start:batch_end]
            batch_indices = indices[batch_start:batch_end]
            
            returned = store.query_batch(batch_queries, batch_indices, top_k=top_k) #1

            for i, returned_unit in enumerate(returned):
                
                returned_unit_matched = returned_unit["matches"]
                trained_positions = [sample["metadata"]["position"] for sample in returned_unit_matched]
                metadata_set = [sample["metadata"]["metadata"] for sample in returned_unit_matched]
                all_candidate_strings = [sample["metadata"]["text"] for sample in returned_unit_matched]
                
                identified_sub_indices, identified_indices, meta_retrieve, timer, smallest_distance = bwamem_align_parallel(
                    all_candidate_strings, 
                    trained_positions, 
                    metadata_set,
                    returned_unit["query"])
                
                series = [int(tup[0]) + int(tup[1]) for tup in zip(identified_sub_indices, identified_indices)]
                if flex:
                    if (is_within_range_of_any_element(returned_unit["index"], series, exactness) and match) or \
                        abs(smallest_distance + 2 * len(returned_unit["query"])) < distance_bound + 1:
                        finer_flag[batch_start + i, 0] = 1
                    else:
                        finer_flag[batch_start + i, 0] = 0
                else:
                    if ((returned_unit["index"] in series) and match) or \
                        abs(smallest_distance + 2 * len(returned_unit["query"])) < 1:
                            finer_flag[batch_start + i, 0] = 1
                    else:
                        finer_flag[batch_start + i, 0] = 0

        return finer_flag[:num_queries, 0]
    
    else: # recall - hotstart
        
        num_queries = len(queries)
        finer_flag = np.zeros((num_queries, 1))

        for batch_start in tqdm(range(0, num_queries, batch_size)):
            batch_end = min(batch_start + batch_size, num_queries)
            batch_queries = queries[batch_start:batch_end]
            batch_indices = indices[batch_start:batch_end]
            
            # TODO: Does not complete the search.
            returned = store.query_batch(batch_queries, batch_indices, hotstart_list=namespaces[batch_start:batch_end], meta_dict=namespace_dict, prioritize=True, top_k=per_k) #1

            for i, returned_unit in enumerate(returned):
                
                returned_unit_matched = returned_unit["matches"]
                trained_positions = [sample["metadata"]["position"] for sample in returned_unit_matched]
                metadata_set = [sample["metadata"]["metadata"] for sample in returned_unit_matched]
                all_candidate_strings = [sample["metadata"]["text"] for sample in returned_unit_matched]
                
                if dictionary_of_values is not None:
                    dictionary_of_values[returned_unit["query"]].append(all_candidate_strings)
                
                identified_sub_indices, identified_indices, meta_retrieve, timer, smallest_distance = bwamem_align_parallel(
                    all_candidate_strings, 
                    trained_positions, 
                    metadata_set,
                    returned_unit["query"])
                
                series = [int(tup[0]) + int(tup[1]) for tup in zip(identified_sub_indices, identified_indices)]
                
                if flex:
                    if (is_within_range_of_any_element(returned_unit["index"], series, exactness) and match) or \
                        abs(smallest_distance + 2 * len(returned_unit["query"])) < distance_bound + 1:
                        finer_flag[batch_start + i, 0] = 1
                    else:
                        finer_flag[batch_start + i, 0] = 0
                else:
                    if ((returned_unit["index"] in series) and match) or \
                        abs(smallest_distance + 2 * len(returned_unit["query"])) < 1:
                            finer_flag[batch_start + i, 0] = 1
                    else:
                        finer_flag[batch_start + i, 0] = 0

                if dictionary_of_values is not None:
                        dictionary_of_values[returned_unit["query"]].append(bool(finer_flag[batch_start + i, 0]))
                        
        if dictionary_of_values is not None:
            return [finer_flag[:num_queries, 0], dictionary_of_values]
        
        return finer_flag[:num_queries, 0]






def read_fasta_chromosomes(file_path):
    """
    Generator function to read a FASTA file and extract each chromosome sequence with its header.

    Parameters:
        file_path (str): Path to the FASTA file.

    Yields:
        tuple: A tuple containing the chromosome header and sequence data.
    """
    with open(file_path, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            if line.startswith('>'):
                
                # If the line starts with '>', it is a chromosome header
                if header is not None:
                    with open("test_cache/logs/headers", "a+") as f:
                        f.write(header)
                        f.write("\n")
                    yield (header, sequence)
                    
                    
                header = line  # Extract the header without '>'
                sequence = ''
                
            else:
                sequence += line
        # Yield the last chromosome entry in the file
        if header is not None:
            with open("test_cache/logs/headers", "a+") as f:
                f.write(header)
                f.write("\n")
            yield (header, sequence)



def is_within_range_of_any_element(X: int, Y: list, Z: int):
    for element in Y:
        if abs(X - element) <= Z:
            return True
    return False



def initialize_pinecone(checkpoint_queue: list[str], 
                        data_queue: list[str], 
                        device:str
):
    
    import torch
    
    for alias in checkpoint_queue:
        
        if alias in checkpoints and checkpoints[alias] != "Baseline":
            received = torch.load(checkpoints[alias], map_location="cpu")
            config = received["config"]
            config.model_config.tokenizer_path = checkpoints["tokenizer"]
            encoder, pooling, tokenizer = model_from_config(config.model_config)
            encoder.load_state_dict(received["model"])
            encoder.eval()
            model_params = {
                                        "tokenizer": tokenizer,
                                        "model": encoder,
                                        "pooling": pooling,
                            }
            baseline = False
            baseline_name = None
            
        elif alias in checkpoints and checkpoints[alias] == "Baseline":
            model_params = None
            baseline = True
            baseline_name = alias
        
        for data_alias in data_queue:
            config = str("config-" + alias + "-" + data_alias).lower()
            store = PineconeStore(
                                    device = torch.device(device),
                                    index_name = str("config-" + alias + "-" + data_alias.replace(",","-")).lower(),
                                    metric = "cosine",
                                    model_params = model_params,
                                    baseline_name=baseline_name,
                                    baseline=baseline
                                )
            
            yield store, data_alias, config


def sample_subsequence(string: str, 
                       min_length: int = 150, 
                       max_length: int = 350
):
    
    subseq_length = random.randint(min_length, max_length)
    # Generate a random starting index
    start_index = random.randint(0, len(string) - subseq_length)
    # Extract the subsequence
    subsequence = string[start_index : start_index + subseq_length]
    return subsequence


def pick_random_lines(path:str = "/home/pholur/dna2vec/tests/data/subsequences_sample_train.txt",
                      k = 100,
                      mode: str = "random",
                      sequences_prior: int = 0,
):
    import pickle
    with open(path, "rb") as file:
        units = pickle.load(file)
    
    if mode == "random":
        random_lines = random.sample(units, k)
        return random_lines
    
    elif mode == "sequenced":
        random_index = random.sample(range(len(units) - k), 1)[0]
        random_lines = units[random_index:random_index+k]
        return random_lines
    
    elif mode == "subsequenced":

        if k % sequences_prior != 0:
            print("k does not divide perfetly by sequences_prior resulting in fewer samples.")
        
        random_lines = []
        real_k = k // sequences_prior
        random_indices = random.sample(range(len(units)), sequences_prior)
        for random_index in random_indices:
            line = units[random_index]
            random_lines.append(line)

            for _ in range(real_k):
                random_lines.append((sample_subsequence(line["text"]), str(random_index)))
        return random_lines

    else:
        raise NotImplementedError("Mode not defined.")


def pick_from_special_gene_list(
            gene_path = "/home/pholur/dna2vec/tests/data/ch2_genes.csv",
            full_path = "/home/pholur/dna2vec/tests/data/NC_000002.12.txt",
            samples = 5000,
):
    import pandas as pd
    df = pd.read_csv(gene_path)
    
    with open(full_path, "r") as f:
        characters = f.read()
    
    sequences = []
    samples_per = samples // df.shape[0]
    
    for _,row in df.iterrows():
        
        label = row["name"]
        indices = row["indices"].split(";")
        big_sequence = characters[int(indices[0]): int(indices[1])]
        
        big_sequence = big_sequence[len(big_sequence) // 2  - 500: len(big_sequence) // 2  + 500]
        t = 0
        sequences.append((big_sequence, label + "_anchor"))

        for _ in range(samples_per):
            sequences.append((sample_subsequence(big_sequence), label))
    
    return sequences


def pick_from_chimp_2a_2b(path_2a: str, 
                          path_2b: str, 
                          samples: int, 
                          per_window: int = 1000
):
    
    per_region = samples // 2
    number_of_cuts = per_region // per_window
    
    with open(path_2a, "r") as f:
        chromosome_2a = f.read()
        
    with open(path_2b, "r") as f:
        chromosome_2b = f.read()
    
    def return_sequences_from_chimp(text_sequence: str, 
                                    label: str
    ):
        random_lines = []
        random_indices = random.sample(range(len(text_sequence) - per_window), number_of_cuts)
        for random_index in random_indices:
            full_sequence = text_sequence[random_index:random_index + per_window]
            for _ in range(per_window):
                random_lines.append((sample_subsequence(full_sequence), label))
        return random_lines
    
    full_sequences = return_sequences_from_chimp(chromosome_2a, "chimp_2a")
    full_sequences.extend(return_sequences_from_chimp(chromosome_2b, "chimp_2b"))
    
    return full_sequences


def pick_from_chromosome3(path, samples, per_window = 1000):
    number_of_cuts = samples // per_window
    
    with open(path, "r") as f:
        chromosome_3 = f.read()

    random_lines = []
    random_indices = random.sample(range(len(chromosome_3) - per_window), number_of_cuts)
    for random_index in random_indices:
        full_sequence = chromosome_3[random_index:random_index + per_window]
        random_lines.append((full_sequence, "anchor_" + str(random_index)))

        for _ in range(per_window):
            random_lines.append((sample_subsequence(full_sequence), str(random_index)))
            
    return random_lines




