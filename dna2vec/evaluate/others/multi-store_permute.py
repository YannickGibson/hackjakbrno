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
         stores: PineconeStore = None,
         config: str = None,
         test_k: int = 200,
         top_k: int = 50,
         generalize: int = 5,
):
    
    per_samples = test_k // len(paths)
    test_lines = {}
    
    for path in paths:
        test_lines.append(pick_random_lines(path, per_samples))
    
    for i,sub_test_lines in enumerate(test_lines):
        
        sub_indices = []
        queries = []
        indices = []
        
        for _,line in enumerate(sub_test_lines):
            
            query = line["text"]
            index = line["position"] # these positions are relative to that data source
            
            sub_index, substring = identify_random_subsequence(
                                    query, 250) #min(generalize*k, len(query)))
            sub_indices.append(sub_index)
            queries.append(substring)
            indices.append(index)
            
        all_matches = None    
        full_indices = [int(ind_base) + int(ind_fine) for ind_base, ind_fine in zip(indices, sub_indices)]
        
        for j,store in enumerate(stores):
            finer_flag = main_align(store, 
                                queries, 
                                full_indices, 
                                top_k,
                                match=i==j)
            
            if all_matches is None:
                all_matches = finer_flag
            else:
                all_matches += finer_flag

        print("Custom perf: ", np.count_nonzero(all_matches) / len(all_matches))
        exit()



if __name__ == "__main__":

    args = parser.parse_args()
    
    data_alias = args.recipes
    checkpoint_queue = args.checkpoints.split(";")
    list_of_data_sources = []
    data_queue = args.recipes.split(";")
    
    stores = [store for (store, _, _) in initialize_pinecone(checkpoint_queue, data_queue, args.device)]
    
    sources = data_queue
    for source in sources:
        if source in data_recipes:
            list_of_data_sources.append(data_recipes[source])
        else:
            list_of_data_sources.append(source)
        
    for topk in args.topk.split(";"):       
        main(list_of_data_sources, 
            stores, 
            None, 
            top_k=int(topk), 
            test_k = args.test_k, 
            generalize=args.generalize)
