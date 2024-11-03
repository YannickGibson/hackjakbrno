from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def is_match_pairwise2(gene: str, sequence: str, verbose_if_matched: bool = True, verbose: bool = False) -> tuple[bool, float]:
    SCORE_PERCENT_THRESHOLD = 0.385

    # Find alignment
    match_bonus = 2
    mismatch_penalty = -3
    gap_extending_penalty = -.3
    gap_penalty = -2
    _alignments = pairwise2.align.globalms(sequence, gene, match_bonus, mismatch_penalty, gap_penalty, gap_extending_penalty, one_alignment_only=True)

    if len(_alignments) == 0:
        return False, 0

    alignment_score = _alignments[0].score
    max_score = (len(gene) * 2)
    score_percent = alignment_score / max_score

    # Evaluate result
    if len(_alignments) > 0 and score_percent > SCORE_PERCENT_THRESHOLD:
        if verbose_if_matched:
            print(format_alignment(*_alignments[0]))
            pass
        if verbose:
            print(f"Score Percentage: {score_percent:.3f}")
        return True, score_percent
    elif verbose:
        print(f"Score Percentage: {score_percent:.3f}")
    return False, score_percent


import time
from Bio import Align
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


def is_match_waterman(gene: str, sequence: str, verbose_if_matched: bool = True, verbose: bool = False) -> tuple[bool, float]:
    # Define threshold
    SCORE_PERCENT_THRESHOLD = 0.92

    # Define arguments
    match_bonus = 2
    mismatch_penalty = -3
    gap_extending_penalty = -.3
    gap_penalty = -2

    # Inicialize aligner
    aligner = Align.PairwiseAligner()
    aligner.match_score = match_bonus
    aligner.mismatch_score = mismatch_penalty
    aligner.extend_gap_score = gap_extending_penalty
    aligner.open_gap_score = gap_penalty
    aligner.mode = 'global'

    # Perform alignment
    alignments = aligner.align(gene, sequence)
    
    if alignments == [] or alignments is None: # exponentially large
        return False, 0

    try:
        if len(alignments) < 1:
            return False, 0
    except OverflowError: # too many matches
        pass
        
    # Extract score
    alignment = alignments[0]
    alignment_score = alignment.score
    max_score = (len(gene) * 2)
    score_percent = alignment_score / max_score

    if score_percent > SCORE_PERCENT_THRESHOLD:
        if verbose_if_matched:
            print(alignments[0])
        if verbose:
            print(f"Score Percentage: {score_percent:.3f}")
        return True, score_percent
    
    if verbose:
        print(f"Score Percentage: {score_percent:.3f}")
    return False, score_percent
