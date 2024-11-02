import sys
sys.path.append("../src/")
sys.path.append("../evaluate/")

import torch
from tqdm import tqdm
import yaml
import argparse
from pinecone_store import PineconeStore
from dna2vec.model import model_from_config
from helpers import pick_random_lines, initialize_pinecone, data_recipes, checkpoints, main_align
import numpy as np
from typing import Literal
import random
from aligners.smith_waterman import bwamem_align_parallel, bwamem_align, bwamem_align_parallel_process

parser = argparse.ArgumentParser(description="Permute evaluate")
parser.add_argument('--recipes', type=str)
parser.add_argument('--checkpoints', type=str)
parser.add_argument('--generalize', type=int)
parser.add_argument('--test_k', type=int)
parser.add_argument("--topk", type=str)
parser.add_argument("--device", type=str)



def identify_random_subsequence(query, length):
    if length > len(query):
        return query  # Return None if the requested length is longer than the query
    start_index = random.randint(0, len(query) - length)  # Generate a random starting index within the valid range
    end_index = start_index + length  # Calculate the end index
    return start_index, query[start_index:end_index]




def main(paths:list,
         store: PineconeStore = None,
         config: str = None,
         test_k: int = 200,
         top_k: int = 50,
         generalize: int = 5,
):
    
    per_samples = test_k // len(paths)
    test_lines = []
    
    for path in paths:
        test_lines.extend(pick_random_lines(path, per_samples))
    
    #finer_flags = np.zeros((test_k, len(test_lines[0]["text"])))   
    finer_flags = np.zeros((test_k, test_k))   

    start = 0
    for k in tqdm(range(1, (len(test_lines) // generalize) + 1)):
        
        sub_indices = []
        queries = []
        indices = []
        
        for i,line in enumerate(test_lines):
            
            query = line["text"]
            index = line["position"]
            
            sub_index, substring = identify_random_subsequence(
                                    query, min(generalize*k, len(query)))
            sub_indices.append(sub_index)
            queries.append(substring)
            indices.append(index)
            
            
        full_indices = [int(ind_base) + int(ind_fine) for ind_base, ind_fine in zip(indices, sub_indices)]
        finer_flag = main_align(store, 
                            queries, 
                            full_indices, 
                            top_k)
        
        finer_flags[:,start:start+generalize] = finer_flag.reshape(-1,1)
        start += generalize

                
        np.savez_compressed(f"test_cache/permute/run_{config}_{test_k}_{top_k}_{generalize}.npz", finer = finer_flags)
                



if __name__ == "__main__":

    args = parser.parse_args()
    
    data_queue = args.recipes.split(";")
    checkpoint_queue = args.checkpoints.split(";")
    
    for (store, data_alias, config) in initialize_pinecone(checkpoint_queue, data_queue, args.device):
        list_of_data_sources = []
        sources = data_alias.split(",")
        for source in sources:
            if source in data_recipes:
                list_of_data_sources.append(data_recipes[source])
            else:
                list_of_data_sources.append(source)
        
        for topk in args.topk.split(";"):       
            main(list_of_data_sources, 
                store, 
                config, 
                top_k=int(topk), 
                test_k = args.test_k, 
                generalize=args.generalize)
        # main(list_of_data_sources, store, config, top_k=50, edit_mode=args.mode, test_k = 1000, generalize=25)
