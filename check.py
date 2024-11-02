from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def is_match_pairwise2(gene: str, sequence: str, verbose_if_matched=True, verbose=False) -> tuple[bool, float]:
    SCORE_PERCENT_THRESHOLD = 0.4

    # Find alignment
    match_bonus = 2
    mismatch_penalty = -3
    gap_extending_penalty = -.3
    gap_penalty = -2
    _alignments = pairwise2.align.localms(sequence, gene, match_bonus, mismatch_penalty, gap_penalty, gap_extending_penalty, one_alignment_only=True)

    if len(_alignments) == 0:
        return False, 0

    score = _alignments[0].score
    max_score = (len(gene) * 3)
    score_percent = score / max_score
    print(score, max_score, score_percent)

    # Evaluate result
    if len(_alignments) > 0 and score > score_percent > SCORE_PERCENT_THRESHOLD:
        if verbose_if_matched:
            print(format_alignment(*_alignments[0]))
        print(f"Score Percentage: {score_percent:.3f}")
        return True, score
    elif verbose:
        print(f"Score Percentage: {score_percent:.3f}")
    return False, score
