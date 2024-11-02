from dna2vec.config_schema import (
    ConfigSchema,
    TrainingConfigSchema,
    DatasetConfigSchema,
    SchedulerConfigSchema,
)
from dna2vec.main import main
from pathlib import Path
import warnings
import torch
from dna2vec.trainer import ContrastiveTrainer


def test_main():
    test_folder = Path(__file__).parent
    CONFIG = ConfigSchema(
        training_config=TrainingConfigSchema(
            max_steps=2,
            batch_size=2,
            log_interval=1,
            save_path=test_folder / "test_models",
            scheduler_config=SchedulerConfigSchema(total_steps=4),

        ),
        dataset_config=DatasetConfigSchema(
            range_mean=200,
            range_std=20,
            subsequence_range_mean=50,
            subsequence_range_std=5,
        ),
    )
    main(CONFIG, wandb_mode="dryrun")

    # test that we can load from checkpoint
    save_path = test_folder / "test_models" / "checkpoint.pt"

    trainer = ContrastiveTrainer.load_from_disk(str(save_path))

    assert trainer.config.training_config.max_steps == CONFIG.training_config.max_steps, "config not loaded correctly"
    trainer.train(max_steps=2)


def test_main_with_accelerator():
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available():  # type: ignore
        device = torch.device("mps")
    else:
        warnings.warn("No accelerators found, skipping test")
        return

    CONFIG = ConfigSchema(
        training_config=TrainingConfigSchema(max_steps=2, batch_size=2, device=device),
        dataset_config=DatasetConfigSchema(
            range_mean=200,
            range_std=20,
            subsequence_range_mean=50,
            subsequence_range_std=5,
        ),
    )
    main(CONFIG, wandb_mode="dryrun")
