from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def is_match_pairwise2(gene: str, sequence: str, verbose_if_matched=True, verbose=False) -> tuple[bool, float]:
    THRESHOLD = 1000

    # Find alignment
    match_bonus = 2
    mismatch_penalty = -3
    gap_extending_penalty = -.3
    gap_penalty = -2
    _alignments = pairwise2.align.localms(sequence, gene, match_bonus, mismatch_penalty, gap_penalty, gap_extending_penalty, one_alignment_only=True)

    # Evaluate result
    if len(_alignments) > 0 and _alignments[0].score > THRESHOLD:
        if verbose_if_matched:
            print(format_alignment(*_alignments[0]))
        elif verbose:
            print(f"Score: {_alignments[0].score}")
        return True, _alignments[0].score
    elif verbose:
        print(f"Score: {_alignments[0].score}")
    return False, _alignments[0].score
