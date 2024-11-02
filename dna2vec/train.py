from pathlib import Path

import torch

from dna2vec.config_schema import (ConfigSchema,
                                   DatasetConfigSchemaUniformSampling,
                                   SchedulerConfigSchema, TrainingConfigSchema)
from dna2vec.dataset import FastaUniformSampler
from dna2vec.main import main

device = torch.device("cuda")
CONFIG = ConfigSchema(
    training_config=TrainingConfigSchema(
        max_steps=100_000,
        batch_size=256,
        device=device,
        log_interval=50,
        accumulation_steps=4,
        scheduler_config=SchedulerConfigSchema(
            max_lr=1e-4,
        ),
    ),
    dataset_config=DatasetConfigSchemaUniformSampling(
        fasta_file = [Path("ncbi_dataset/ncbi_dataset/data/GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna")], 
        range_min = 800,
        range_max = 9500,
        subsequence_range_min = 800,
        subsequence_range_max = 1200, 
        dataset=FastaUniformSampler,
        sampling_strategy="random_subsequence_uppercase",
        ),
)


main(CONFIG, watch_watch=True)
