import importlib
import importlib.util
import json
import os
import subprocess
import urllib.request
from pathlib import Path
from typing import Literal, Optional

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from pydantic import BaseModel
from tqdm import tqdm

from dna2vec.config_schema import ConfigSchema

CACHE_DIR = Path.home() / ".cache" / "dna2vec"


def download_human_reference_genome(
    reference: Literal["GRCh38", "CHM13"] = "CHM13",
    force: bool = False, use_uncertified_ssl: bool = False
) -> Path:
    """
    Download the human reference genome to the cache directory. If the file already
    exists, it will not be downloaded again unless force is set to True.

    Args:
        reference: Which reference genome to download. Currently supported are
            "GRCh38" and "CHM13".
        force: If True, the file will be downloaded even if it already exists.
        use_uncertified_ssl: If True, use uncertified ssl when downloading the file.
    """
    cache_dir = get_cache_dir()
    cache_dir.mkdir(exist_ok=True, parents=True)

    urls = {
        "GRCh38":
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz",
        "CHM13":
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    }
    url = urls[reference]
    output = cache_dir / "GRCh38_latest_genomic.fna.gz"
    unpacked = cache_dir / "GRCh38_latest_genomic.fna"

    if not unpacked.exists() or force:
        download_url(url, output, use_uncertified_ssl=use_uncertified_ssl)
        subprocess.run(["gzip", "-d", str(output)], check=True)

    return cache_dir / "GRCh38_latest_genomic.fna"


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path, use_uncertified_ssl: bool = False):
    if use_uncertified_ssl:
        import ssl

        ssl._create_default_https_context = ssl._create_unverified_context

    with DownloadProgressBar(
        unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]
    ) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def get_cache_dir() -> Path:
    """
    Get the cache directory for dna2vec. Can be overridden by setting the environment
    variable DNA2VEC_CACHE_DIR.
    """
    cache_dir = os.environ.get("DNA2VEC_CACHE_DIR")
    if cache_dir is not None:
        return Path(cache_dir)
    return CACHE_DIR


def get_config_from_path(config_path: Path) -> ConfigSchema:
    spec = importlib.util.spec_from_file_location(config_path.stem, config_path)
    module = importlib.util.module_from_spec(spec)  # type: ignore
    spec.loader.exec_module(module)  # type: ignore

    CONFIG = module.CONFIG
    return CONFIG


def cfg_to_wandb_dict(cfg: BaseModel) -> dict:
    """
    Convert a ConfigSchema instance to a dictionary that can be logged to wandb by
    recursively unnesting the dict structure and joining then with a dot.

    E.g. {"model_config": {"embedding_dim": 384}} -> {"model_config.embedding_dim": 384}
    if it is a callable, the name of the callable is used instead of the value

    Args:
        cfg: ConfigSchema
    """
    cfg_dict = cfg.dict()

    def _cfg_to_wandb_dict(cfg_dict):
        for key, value in cfg_dict.items():
            if isinstance(value, dict):
                for k, v in _cfg_to_wandb_dict(value):
                    yield f"{key}.{k}", v
            elif callable(value):
                try:
                    name = value.__name__
                except:
                    name = value.__class__.__name__

                yield key, name
            else:
                # check if value if json serializable
                try:
                    json.dumps(value)
                except TypeError:
                    # if not, convert to string
                    value = str(value)
                yield key, value

    return dict(_cfg_to_wandb_dict(cfg_dict))


def load_human_reference_genome(path: Optional[Path] = None) -> FastaIterator:
    """
    Load the human reference genome from the given path. If no path is given, the

    Args:
        path: Path to the human reference genome. If None, the genome will be
            downloaded to the cache directory and loaded from there.

    Returns:
        An iterator of sequence records.
    """
    # if path is None:
    #     path = download_human_reference_genome()
    print(path)
    return SeqIO.parse(path, "fasta")
