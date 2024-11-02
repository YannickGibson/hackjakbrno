from pathlib import Path
from typing import List

import pytest
from pysam.libcalignedsegment import AlignedSegment

from dna2vec.simulate import (
    ReadAndReference,
    load_simulated_reads_from_disk,
    map_reads_to_reference,
    simulate_mapped_reads,
    simulate_reads_to_disk,
)


@pytest.fixture
def simulated_path() -> Path:
    project_path = Path(__file__).parent.parent
    simulated_file = project_path / "tests" / "data" / "sample_simulated_reads"

    assert simulated_file.with_suffix(".fq").exists()
    assert simulated_file.with_suffix(".sam").exists()
    assert simulated_file.with_suffix(".aln").exists()

    return simulated_file


@pytest.fixture
def reads(simulated_path) -> List[AlignedSegment]:
    return load_simulated_reads_from_disk(simulated_path)


def test_reads(reads):
    for read in reads:
        assert read.query_sequence is not None
        assert read.query_alignment_sequence is not None
        assert read.query_sequence == read.query_alignment_sequence


def test_simulate_reads_to_disk(simulated_path: Path):
    n_reads = 10
    read_length = 100

    out_path = simulate_reads_to_disk(
        n_reads_pr_amplicon=n_reads,
        read_length=read_length,
        output_path=simulated_path,
        insertion_rate=0.1,
        deletion_rate=0.1,
    )

    assert out_path == simulated_path
    assert out_path.with_suffix(".fq").exists()
    assert out_path.with_suffix(".sam").exists()
    assert out_path.with_suffix(".aln").exists()

    # test that they can be read in
    reads = load_simulated_reads_from_disk(out_path)

    assert isinstance(reads, list)
    assert (
        len(reads) > n_reads
    ), f"there should be at least {n_reads} reads per amplicon"
    assert isinstance(reads[0], AlignedSegment)
    assert (
        len(reads[0].query_sequence) == read_length
    ), f"read length should be {read_length}"


def test_load_reads_from_disk(simulated_path: Path):
    reads = load_simulated_reads_from_disk(simulated_path)
    read = reads[0]

    assert isinstance(read, AlignedSegment)


def is_approximately_the_same(seq1: str, seq2: str, max_prob_diff: float = 0.1):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    assert len(seq1) == len(seq2)
    n_diff = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            n_diff += 1

    prob_diff = n_diff / len(seq1)
    assert (
        prob_diff < max_prob_diff
    ), f"sequences should be approximately the same, but are different by {prob_diff}"


def test_map_reads_to_reference(reads: List[AlignedSegment]):
    mapped_reads = map_reads_to_reference(reads)

    assert isinstance(mapped_reads, list)
    for m_read in mapped_reads:
        assert isinstance(m_read, ReadAndReference)
        assert m_read.is_mapped(), "all reads should be mapped"
        read, ref = m_read.get_pair_as_string()
        is_approximately_the_same(read, ref)  # type: ignore


def test_simulate_reads():
    mapped_reads = simulate_mapped_reads(
        n_reads_pr_amplicon=10,
        read_length=100,
    )

    assert isinstance(mapped_reads, list)
    assert len(mapped_reads) > 10

    m_read = mapped_reads[0]
    assert isinstance(m_read, ReadAndReference)
    assert m_read.id is not None