from pathlib import Path

from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqRecord import SeqRecord

from dna2vec.utils import (
    cfg_to_wandb_dict,
    download_human_reference_genome,
    get_config_from_path,
    load_human_reference_genome,
)


def test_cfg_to_wandb_dict():
    cfg = get_config_from_path(
        Path(__file__).parents[1] / "configs" / "default_config.py"
    )
    print(cfg)
    wandb_dict = cfg_to_wandb_dict(cfg)
    # check that it is json serializable
    import json

    print(wandb_dict)
    json.dumps(wandb_dict)


def test_get_human_reference_genome():
    path = download_human_reference_genome(use_uncertified_ssl=True)
    assert isinstance(path, Path)
    assert path.exists()
    assert path.is_file()
    assert path.suffix == ".fna"
    assert path.name == "GRCh38_latest_genomic.fna"


def test_load_human_reference_genome():
    ref = load_human_reference_genome()

    assert isinstance(ref, FastaIterator)
    seq = next(ref)
    assert isinstance(seq, SeqRecord)
