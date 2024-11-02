import os
os.environ["DNA2VEC_CACHE_DIR"] = "/mnt/SSD2/pholur/dna2vec"
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "2,3"

from helpers import raw_fasta_files
import numpy as np
from tqdm import tqdm

from helpers import initialize_pinecone, is_within_range_of_any_element, main_align


from pathlib import Path
from datetime import datetime
import logging
import random
from itertools import product

from dna2vec.config_schema import DatasetConfigSchemaUniformSampling
from dna2vec.simulate import simulate_mapped_reads

from pinecone_store import PineconeStore

grid = {
    "read_length": [250], #[150, 300, 500],
    "insertion_rate": [0.0001],
    "deletion_rate" : [0.0001],
    "qq": [(60,90)], # https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf
    "topk": [50],
    "distance_bound": [25],
    "exactness": [10]
}


###
import argparse
parser = argparse.ArgumentParser(description="Grid Search")
parser.add_argument('--recipe', type=str)
parser.add_argument('--checkpoints', type=str)
parser.add_argument('--topk', type=int)
parser.add_argument('--test', type=int)
parser.add_argument('--system', type=str)
parser.add_argument('--device', type=str)
###




if __name__ == "__main__":
    
    args = parser.parse_args()
    fasta_file_paths = []
    
    data_queue = args.recipe.split(";")
    for data_alias in data_queue:
        fasta_file_paths.append(raw_fasta_files[data_alias])

    checkpoint_queue = args.checkpoints.split(";")
    stores = [store for (store, _, _) in initialize_pinecone(checkpoint_queue, data_queue, args.device)]

    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d_%H_%M_%S")
    
    all_reads = []
    for fasta_file_path in fasta_file_paths:
        for read_length, insertion_rate, deletion_rate, quality, \
                topk, distance_bound, exactness in \
                product(grid["read_length"], grid["insertion_rate"], 
                        grid["deletion_rate"], grid["qq"], grid["topk"], 
                        grid["distance_bound"], grid["exactness"]):
                    
            mapped_reads = simulate_mapped_reads(
                n_reads_pr_amplicon=args.test,
                read_length=read_length,
                insertion_rate=insertion_rate,
                deletion_rate=deletion_rate,
                reference_genome=fasta_file_path,
                sequencing_system=args.system,
                quality = quality
            )
            
            all_reads.append(mapped_reads)
            break

    all_successes = 0
    count = 0
    
    for i,mapped_reads in enumerate(all_reads):
    
        perf_read = 0
        perf_true = 0

        queries = []
        small_indices = []
        start_indices = []


        for sample in tqdm(mapped_reads):
            queries.append(sample.read.query_sequence)
            small_indices.append(int(sample.read.reference_start))
            start_indices.append(int(sample.seq_offset))

        ground_truth = [index_main + inter_fine for index_main, inter_fine in zip(start_indices, small_indices)]    
        total_perf = None
        
        for j,store in enumerate(stores):  
            results = main_align(store, queries, ground_truth, grid["topk"][0], 
                                 exactness=grid["exactness"][0], distance_bound=grid["distance_bound"][0], 
                                 flex = True, match=i==j)
            if total_perf is None:
                total_perf = results
            else:
                total_perf += results
        
        successes = np.count_nonzero(results)
        all_successes += successes
        count += len(results)
        
    print(all_successes / count)
    exit()

