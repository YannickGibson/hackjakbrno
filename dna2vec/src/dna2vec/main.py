"""
cli interface for training the contrastive self-supervised sequence model
"""


from functools import partial
from pathlib import Path
import sched
from typing import Optional

import typer
from torch.utils.data import DataLoader

import wandb
from dna2vec.config_schema import ConfigSchema
from dna2vec.dataset import collate_fn
from dna2vec.model import model_from_config
from dna2vec.trainer import ContrastiveTrainer
from dna2vec.utils import cfg_to_wandb_dict, get_config_from_path


def main(config: ConfigSchema, wandb_mode: str = "online", watch_watch: bool = True):
    """
    Args:
        config: The configuration object
    """
    model_cfg = config.model_config
    training_cfg = config.training_config
    dataset_cfg = config.dataset_config

    # MODEL: Model and tokenizer
    model, pooler, tokenizer = model_from_config(model_cfg)
    _collate_fn = partial(collate_fn, tokenizer=tokenizer)

    # DATASET: Dataset creation
    dataset_kwargs = dataset_cfg.dict()
    dataset_fn = dataset_kwargs.pop("dataset")
    dataset = dataset_fn(**dataset_kwargs)

    # TRAINING: Optimizer, data loader and trainer
    optimizer_cfg = training_cfg.optimizer_config
    optimizer = training_cfg.optimizer(model.parameters(), **optimizer_cfg.dict())
    if training_cfg.scheduler_config.total_steps is None:
        training_cfg.scheduler_config.total_steps = int(training_cfg.max_steps / training_cfg.accumulation_steps)
    scheduler = training_cfg.scheduler(optimizer, **training_cfg.scheduler_config.dict())
    dataloader = DataLoader(
        dataset, batch_size=training_cfg.batch_size, collate_fn=_collate_fn
    )
    sim = training_cfg.similarity(temperature=training_cfg.temperature)

    trainer = ContrastiveTrainer(
        encoder=model,
        pooling=pooler,
        loss=training_cfg.loss,
        optimizer=optimizer,
        train_dataloader=dataloader,
        device=training_cfg.device,
        similarity=sim,
        scheduler=scheduler, 
        config=config,
        tokenizer=tokenizer,
    )

    # log config to wandb
    wandb.init(project="dna2vec", config=cfg_to_wandb_dict(config), mode=wandb_mode)
    if watch_watch:
        wandb.watch(model, log="all", log_freq=1, log_graph=True) # just for debugging

    # TRAINING: Training loop
    trainer.train(
        max_steps=training_cfg.max_steps,
        log_interval=training_cfg.log_interval,
    )


def main_cli(config_path: Optional[str] = None):
    """
    The main cli interface for training the model

    Args:
        cfg: Path to the configuration file, should be a .py file containing a Config instance
    """
    if config_path is None:
        cfg = Path(__file__).parents[2] / "configs" / "default_config.py"
    else:
        cfg = Path(config_path)

    config = get_config_from_path(cfg)
    main(config)


if __name__ == "__main__":
    typer.run(main_cli)
