"""
Base configurations
"""

from functools import partial
from pathlib import Path
from typing import Literal, Optional, Type

import torch
from pydantic import BaseModel
from torch import nn
from torch.optim.lr_scheduler import LRScheduler, OneCycleLR
from torch.utils.data import Dataset

from dna2vec.dataset import FastaSamplerDataset
from dna2vec.model import AveragePooler, SinusoidalPositionalEncoding
from dna2vec.similarity import SimilarityWithTemperature

scheduler = partial(OneCycleLR, 
                       max_lr = 1e-4, # Upper learning rate boundaries in the cycle for each parameter group
                       anneal_strategy = 'cos')

project_path = Path(__file__).parent.parent.parent
tokenizer_path = (
    project_path / "src" / "model" / "tokenizers" / "dna_tokenizer_10k.json"
)


class ModelConfigSchema(BaseModel):
    embedding_dim: int = 1024 # changed from 384
    dim_feedforward: int = 1536
    vocab_size: Optional[int] = None  # derived from tokenizer
    num_heads: int = 12
    num_layers: int = 6
    dropout: float = 0.1
    activation: Literal["relu", "gelu"] = "gelu"
    pos_embedding: Type[nn.Module] = SinusoidalPositionalEncoding
    pooling: nn.Module = AveragePooler()
    max_position_embeddings: int = 2048 # changed from 1024
    tokenizer_path: Path = tokenizer_path
    model_path: Optional[Path] = None # where to load the model from

    class Config:
        arbitrary_types_allowed = True


class OptimizerConfigSchema(BaseModel):
    lr: float = 0.0001
    betas: tuple[float, float] = (0.9, 0.999)
    weight_decay: float = 0.01
    eps = 1e-8

class SchedulerConfigSchema(BaseModel):
    max_lr: float = 1e-4
    anneal_strategy: Literal["cos", "linear", "polynomial", "constant"] = "cos"
    total_steps: Optional[int] =None # derived from training config




class TrainingConfigSchema(BaseModel):
    batch_size: int = 64
    optimizer: Type[torch.optim.Optimizer] = torch.optim.Adam
    optimizer_config: OptimizerConfigSchema = OptimizerConfigSchema()
    similarity: Type[nn.Module] = SimilarityWithTemperature
    temperature: float = 0.05  # default derived from SimCSE
    loss: nn.Module = nn.CrossEntropyLoss()
    accumulation_steps: int = 1
    max_grad_norm: float = 1.0
    scheduler: Type[LRScheduler] = OneCycleLR
    scheduler_config: SchedulerConfigSchema = SchedulerConfigSchema()
    save_path: Path = project_path / "models"

    max_steps: int = 1000
    log_interval: int = 100
    device: torch.device = torch.device("cpu")

    class Config:
        arbitrary_types_allowed = True


class BaseDatasetConfigSchema(BaseModel):
    dataset: Type[Dataset]
    fasta_file: Path


class DatasetConfigSchemaUniformSampling(BaseDatasetConfigSchema):
    fasta_file: list
    range_min = 1000
    range_max = 2000
    subsequence_range_min = 100
    subsequence_range_max = 500
    sampling_strategy: str = "random_subsequence"


class DatasetConfigSchema(BaseDatasetConfigSchema):
    dataset: Type[Dataset] = FastaSamplerDataset
    fasta_file: Path = project_path / "tests" / "data" / "NC_000002.12.txt"
    range_mean: float = 3000 # 3000, changed from 1000
    range_std: float = 3500 # ?, changed from 100
    subsequence_range_mean: float = 1000 # 1000, 200
    subsequence_range_std: float = 200 # spocitame 945642
    sampling_strategy: str = "random_subsequence"


class ConfigSchema(BaseModel):
    model_config: ModelConfigSchema = ModelConfigSchema()
    training_config: TrainingConfigSchema = TrainingConfigSchema()
    dataset_config: BaseDatasetConfigSchema = DatasetConfigSchema()


