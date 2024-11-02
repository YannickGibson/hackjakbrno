import time
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


def calculate_smith_waterman_distance(  string1,
                                        string2,
                                        match_score = 2, # do not change
                                        mismatch_penalty = -1,
                                        open_gap_penalty = -0.5,
                                        continue_gap_penalty = -0.1,
                                        debug = False):
    start = time.time()
    # Perform the Smith-Waterman alignment
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = open_gap_penalty
    aligner.extend_gap_score = continue_gap_penalty
    aligner.mismatch_score = mismatch_penalty
    aligner.match_score = match_score
    aligner.mode = 'local'
    alignments = aligner.align(
        string1, string2)
    # alignments = pairwise2.align.localms(
    #     string1, string2, match_score, mismatch_penalty, 
    #     open_gap_penalty, continue_gap_penalty
    # )
    if debug:
        print("Best match:")
        print(alignments[0])
    
    if alignments == [] or alignments is None: # exponentially large
        return {"elapsed time": time.time() - start, 
            "distance": float("inf"), 
            "begins":[]}

    try:
        if len(alignments) < 1:
            return {"elapsed time": time.time() - start, 
            "distance": float("inf"), 
            "begins":[]}
    except OverflowError: # too many matches
        pass
        
    alignment = alignments[0]
    alignment_score = None

    alignment_score = alignment.score
    # Calculate the Smith-Waterman distance
    smith_waterman_distance = -alignment_score
    begins = []
    begins.append(alignment.aligned[0][0][0]) # target, first match, first index TODO what if many matches - rare at this length
    
    
    return {"elapsed time": time.time() - start, 
            "distance": smith_waterman_distance, 
            "begins":list(set(begins))}














def bwamem_align(all_candidate_strings: list[str], 
                 trained_positions: list[str], 
                 metadata_set: list[str], 
                 substring: list[str]
):
    '''
    Currently implemented sequentially
    '''
    
    total_time = 0
    refined_results = defaultdict(list)
    
    for long_string, train_pos, metadata in zip(
        all_candidate_strings, trained_positions, metadata_set):
        
        returned_object = calculate_smith_waterman_distance(long_string, substring)
        total_time += returned_object["elapsed time"]

        for starting_sub_index in returned_object["begins"]:
            refined_results[returned_object["distance"]].append(
                (starting_sub_index, train_pos, metadata)
            )
    
    try:
        smallest_key = min(refined_results.keys())
        
    except ValueError:
        return [], [], [], total_time
    
    smallest_values = refined_results[smallest_key]
    
    identified_sub_indices = []
    identified_indices = []
    metadata_indices = []
    
    for term in smallest_values:
        identified_sub_indices.append(term[0])
        identified_indices.append(term[1])
        metadata_indices.append(term[2])      
          
    return identified_sub_indices, identified_indices, metadata_indices, total_time, smallest_key


def process_single_string(args:tuple):
    long_string, substring, train_pos, metadata = args
    returned_object = calculate_smith_waterman_distance(long_string, substring)
    return returned_object["distance"], returned_object["begins"], train_pos, metadata, returned_object["elapsed time"]


def bwamem_align_parallel(all_candidate_strings: list[str], 
                          trained_positions: list[str], 
                          metadata_set: list[str], 
                          substring: str, 
                          max_workers: int=50
):
    
    total_time = time.time()
    refined_results = defaultdict(list)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_single_string, 
                                    zip(all_candidate_strings, 
                                        [substring]*len(all_candidate_strings), 
                                        trained_positions, 
                                        metadata_set)))

    for distance, begins, train_pos, metadata, _ in results:
        for starting_sub_index in begins:
            refined_results[distance].append(
                (starting_sub_index, train_pos, metadata)
            )
    
    try:
        smallest_key = min(refined_results.keys())
    except ValueError:
        return [], [], [], time.time() - total_time
    

        
    smallest_values = refined_results[smallest_key]
    
    identified_sub_indices = []
    identified_indices = []
    metadata_indices = []
    
    for term in smallest_values:
        identified_sub_indices.append(term[0])
        identified_indices.append(term[1])
        metadata_indices.append(term[2])      
          
    return  identified_sub_indices, \
            identified_indices, \
            metadata_indices, \
            time.time() - total_time, \
            smallest_key




def process_batch(batch):
    return [process_single_string(args) for args in batch]

def bwamem_align_parallel_process(all_candidate_strings: list[str], 
                          trained_positions: list[str], 
                          metadata_set: list[str], 
                          substring: str, 
                          batch_size: int = 5,
                          max_workers: int = 10
):
    total_time = time.time()
    refined_results = defaultdict(list)

    batches = [
        batch for batch in zip(
            all_candidate_strings, 
            [substring] * len(all_candidate_strings), 
            trained_positions, 
            metadata_set
        )
    ]
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        batch_results = list(executor.map(process_batch, [batches[i:i+batch_size] for i in range(0, len(batches), batch_size)]))

    for batch_result in batch_results:
        for results in batch_result:
            distance, begins, train_pos, metadata, _ = results
            for starting_sub_index in begins:
                refined_results[distance].append(
                    (starting_sub_index, train_pos, metadata)
                )
    
    if not refined_results:
        return [], [], [], time.time() - total_time
        
    smallest_key = min(refined_results.keys())
    smallest_values = refined_results[smallest_key]
    
    identified_sub_indices = []
    identified_indices = []
    metadata_indices = []
    
    for term in smallest_values:
        identified_sub_indices.append(term[0])
        identified_indices.append(term[1])
        metadata_indices.append(term[2])      
          
    return  identified_sub_indices, \
            identified_indices, \
            metadata_indices, \
            time.time() - total_time, \
            smallest_key