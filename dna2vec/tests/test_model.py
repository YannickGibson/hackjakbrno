import warnings

import pytest
import torch

from dna2vec.model import Encoder, LearnedPositionalEncoding


def test_model():
    # test constructor
    encoder = Encoder(
        vocab_size=5,
        embedding_dim=10,
        max_position_embeddings=10,
        num_heads=2,
        num_layers=2,
    )

    assert encoder.embedding.weight.shape == (5, 10)

    # forward
    x = torch.tensor([[1, 2, 3, 4, 4], [1, 2, 3, 4, 0]])
    emb = encoder(x)
    # assert emb.names == ("batch", "sequence", "embedding")
    assert emb.shape == (2, 5, 10)


def test_learned_positional_embedding():
    # test positional embedding
    pos_emb = LearnedPositionalEncoding
    encoder = Encoder(
        vocab_size=5,
        embedding_dim=10,
        max_position_embeddings=10,
        num_heads=2,
        num_layers=2,
        pos_embedding=pos_emb,
    )

    # forward
    x = torch.tensor([[1, 2, 3, 4, 4], [1, 2, 3, 4, 0]])
    emb = encoder(x)
    # assert emb.names == ("batch", "sequence", "embedding")
    assert emb.shape == (2, 5, 10)
