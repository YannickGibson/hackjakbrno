from helpers import raw_fasta_files
import numpy as np
from tqdm import tqdm

from helpers import initialize_pinecone, is_within_range_of_any_element
from aligners.smith_waterman import bwamem_align, bwamem_align_parallel

from pathlib import Path
from datetime import datetime
import logging
import random
from itertools import product

from dna2vec.config_schema import DatasetConfigSchemaUniformSampling
from dna2vec.simulate import simulate_mapped_reads

from pinecone_store import PineconeStore

import os
os.environ["DNA2VEC_CACHE_DIR"] = "/mnt/SSD2/pholur/dna2vec"

grid = {
    "read_length": [250], #[150, 300, 500],
    "insertion_rate": [0.0, 0.0001, 0.01],
    "deletion_rate" : [0.0, 0.0001, 0.01],
    "qq": [(60,90), (30,60)], # https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf
    "topk": [5, 25, 50],
    "distance_bound": [0, 10, 20],
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




def evaluate(store: PineconeStore, 
             query: str, 
             top_k: int
):
    returned = store.query([query], top_k=top_k)["matches"]
    trained_positions = [sample["metadata"]["position"] for sample in returned]
    metadata_set = [sample["metadata"]["metadata"] for sample in returned]
    all_candidate_strings = [sample["metadata"]["text"] for sample in returned]
    
    # Parallelized BWA mem - FAST
    identified_sub_indices, identified_indices, _, timer, smallest_distance = bwamem_align_parallel(
                                                                                all_candidate_strings, 
                                                                                trained_positions, 
                                                                                metadata_set,
                                                                                query)

    return [(int(full) + int(fine))
            for (full, fine) in zip(identified_indices, identified_sub_indices)], \
                    timer, smallest_distance


if __name__ == "__main__":
    
    args = parser.parse_args()
    fasta_file_path = raw_fasta_files[args.recipe]
    
    data_queue = args.recipe.split(";")
    checkpoint_queue = args.checkpoints.split(";")
    
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d_%H_%M_%S")

    logging.basicConfig(filename = Path(os.environ["DNA2VEC_CACHE_DIR"]) / "Logs" / f"log_{formatted_date}", 
                        level=logging.INFO)
    
    logging.info("Parameters:")
    
    f = open(Path(os.environ["DNA2VEC_CACHE_DIR"]) / "Results" / f"result_{formatted_date}.csv", "w+")
        
    for arg in vars(args):
        f.write(f"# {arg}: {getattr(args, arg)}\nqqq")
        
    f.write("Quality,Read length,Insertion rate,Deletion rate,TopK,Accuracy\n")
    f.flush()
    
    for arg in vars(args):
        logging.info(f"{arg}: {getattr(args, arg)}")
        
    for (store, _, _) in initialize_pinecone(checkpoint_queue, data_queue, args.device):
        
        for read_length, insertion_rate, deletion_rate, quality, topk in \
            product(grid["read_length"], grid["insertion_rate"], grid["deletion_rate"], grid["qq"], grid["topk"]):

            perf_read = 0
            perf_true = 0
            count = 0

            mapped_reads = simulate_mapped_reads(
                n_reads_pr_amplicon=args.test,
                read_length=read_length,
                insertion_rate=insertion_rate,
                deletion_rate=deletion_rate,
                reference_genome=fasta_file_path,
                sequencing_system=args.system,
                quality = quality
            )
            
            for sample in tqdm(mapped_reads):
                query = sample.read.query_sequence
                beginning = sample.read.reference_start

                matches, timer, smallest_distance = evaluate(store, query, topk)
                full_start = beginning + sample.seq_offset
                if is_within_range_of_any_element(full_start, matches, args.exactness) or \
                    abs(smallest_distance + 2*len(query)) < args.distance_bound + 1: # the 1 here helps with instabilities
                    perf_read += 1
                else:
                    logging.info(f"############# Error: \nQuery: {query}\nStart: {beginning}\nOriginal: \
{sample.reference}\nMatches: {matches}\n{str(quality).replace(',',';')},{read_length},{insertion_rate},\
{deletion_rate},{topk},{perf_read/(count+0.0001)} \n #############") 
                count += 1
            
            f.write(f"{str(quality).replace(',',';')},{read_length},{insertion_rate},{deletion_rate},{topk},{perf_read/count}\n")
            f.flush()
                        
    f.close()
                        
                    