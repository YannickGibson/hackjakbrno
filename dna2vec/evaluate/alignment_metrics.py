import time
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align

def calculate_smith_waterman_distance(  string1,
                                        string2,
                                        match_score = 2,
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
    begins.append(alignment.aligned[0][0][0]) # target, first match, first index TODO what if many matches
    
    
    return {"elapsed time": time.time() - start, 
            "distance": smith_waterman_distance, 
            "begins":list(set(begins))}


