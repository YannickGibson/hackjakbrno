"""
The trainer module contains the Trainer class, which is responsible for training the model contrastive learning
"""

from functools import partial
from pathlib import Path
from typing import Optional

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.optim.lr_scheduler import LRScheduler
import wandb
from dna2vec.config_schema import ConfigSchema
from dna2vec.dataset import collate_fn
from dna2vec.model import model_from_config

from dna2vec.tokenizer import BPTokenizer

class ContrastiveTrainer:
    def __init__(
        self,
        encoder: nn.Module,
        pooling: nn.Module,
        similarity: nn.Module,
        loss: nn.Module,
        optimizer: torch.optim.Optimizer,
        train_dataloader: DataLoader,
        scheduler: LRScheduler,
        device: torch.device,
        config: ConfigSchema,
        tokenizer: BPTokenizer,
    ):
        self.encoder = encoder
        self.loss = loss
        self.train_dataloader = train_dataloader
        self.optimizer = optimizer
        self.device = device
        self.similarity = similarity
        self.pooling = pooling
        self.scheduler = scheduler
        self.config = config
        self.tokenizer = tokenizer
        self.best_loss = float("inf")

        self.training_config = config.training_config


    def model_to_device(self, device: Optional[torch.device] = None) -> None:
        """
        Move the model to the specified device

        Args:
            device: Device to move the model to
        """
        if device is not None:
            self.device = device
        self.encoder.to(self.device)

    def dict_to_device(self, d: dict):
        for k, v in d.items():
            if isinstance(v, torch.Tensor):
                d[k] = v.to(self.device)

    def train(
        self,
        max_steps: Optional[int] = None,
        log_interval: int = 100,
    ) -> None:
        """
        Train the model for the specified number of steps

        Args:
            steps: Number of steps to train the model for
            log_interval: Number of steps after which to log the training loss
        """


        self.model_to_device()
        self.encoder.train()
        for step, (x_1, x_2) in enumerate(self.train_dataloader):
            self.dict_to_device(x_1)
            self.dict_to_device(x_2)

            #print(f"Current step {step}")
            if max_steps is not None and step >= max_steps:
                break


            last_hidden_x_1 = self.encoder(**x_1) # long sequence
            last_hidden_x_2 = self.encoder(**x_2) # subsequence
            y_1 = self.pooling(last_hidden_x_1, attention_mask=x_1["attention_mask"])
            y_2 = self.pooling(last_hidden_x_2, attention_mask=x_2["attention_mask"])

            # Calculate similarity
            # y_1 # names: [batch, embedding]
            # y_2 # names: [batch, embedding]
            sim = self.similarity(y_1.unsqueeze(1), y_2.unsqueeze(0)) # outer-product

            labels = torch.arange(sim.size(0)).long().to(self.device)

            loss = self.loss(sim, labels)

            loss = loss / self.training_config.accumulation_steps
            loss.backward()
            if (step + 1) % self.training_config.accumulation_steps == 0:
                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(self.encoder.parameters(), self.training_config.max_grad_norm) # type: ignore

                self.optimizer.step()
                self.optimizer.zero_grad()
                self.scheduler.step()


            #print(f"{step%log_interval}")
            if step % log_interval == 0:
                current_lr = self.optimizer.param_groups[0]["lr"]
                wandb.log({"loss": loss, "step": step, "lr": current_lr})

            # save the model
            if loss < self.best_loss:
                self.best_loss = loss
                self.save_to_disk()


            # trying to resolve the CUDA out of memory error
            if step % 1000 == 0:
                torch.cuda.empty_cache()

            # delete the tensors: https://discuss.pytorch.org/t/gpu-memory-consumption-increases-while-training/2770/3?u=nagabhushansn95
            del y_1, y_2, sim, labels, loss
            for k, v in x_1.items():
                del v
            for k, v in x_2.items():
                del v
            del last_hidden_x_1, last_hidden_x_2, x_1, x_2

    def save_to_disk(self, path: Optional[str] = None):
        if path is None:
            save_path = self.config.training_config.save_path
        else:
            save_path = Path(path)
        save_path.mkdir(parents=True, exist_ok=True)

        # record model state train/eval
        state = self.encoder.training

        save_dict = {
            "model": self.encoder.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "scheduler": self.scheduler.state_dict(),
            "config": self.config,
        }

        torch.save(save_dict, save_path / "checkpoint.pt")

    @staticmethod
    def load_from_disk(path: str) -> "ContrastiveTrainer":
    
        checkpoint = torch.load(path)
        config = checkpoint["config"]
        encoder, pooling, tokenizer = model_from_config(config.model_config)
        encoder.load_state_dict(checkpoint["model"])
        optimizer_state = checkpoint["optimizer"]
        scheduler_state = checkpoint["scheduler"]

        dataset_kwargs = config.dataset_config.dict()
        dataset_fn = dataset_kwargs.pop("dataset")
        dataset = dataset_fn(**dataset_kwargs)

        _collate_fn = partial(collate_fn, tokenizer=tokenizer)

        train_cfg = config.training_config
        dataloader = DataLoader(
            dataset, batch_size=train_cfg.batch_size, collate_fn=_collate_fn
        )

        sim = train_cfg.similarity(temperature=train_cfg.temperature)

        # recreate optimizer and scheduler states
        opt_kwargs = train_cfg.optimizer_config.dict()
        optimizer_ = train_cfg.optimizer(encoder.parameters(), **opt_kwargs)
        optimizer_.load_state_dict(optimizer_state)
        scheduler_ = train_cfg.scheduler(optimizer_, **train_cfg.scheduler_config.dict())
        scheduler_.load_state_dict(scheduler_state)

        trainer = ContrastiveTrainer(
            encoder=encoder,
            pooling=pooling,
            loss=train_cfg.loss,
            optimizer=optimizer_,
            scheduler=scheduler_,
            device=train_cfg.device,
            train_dataloader=dataloader,
            similarity=sim,
            config=config,
            tokenizer=tokenizer,
        )
        
        return trainer


