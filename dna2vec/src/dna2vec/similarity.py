from torch import nn


class SimilarityWithTemperature(nn.Module):
    """
    Dot product or cosine similarity

    derived from:
    https://github.com/princeton-nlp/SimCSE/blob/13361d0e29da1691e313a94f003e2ed1cfa97fef/simcse/models.py#L35
    """

    def __init__(self, temperature):
        super().__init__()
        self.temp = temperature
        self.cos = nn.CosineSimilarity(dim=-1)

    def forward(self, x, y):
        return self.cos(x, y) / self.temp
