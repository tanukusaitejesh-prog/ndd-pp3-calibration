"""A deterministic stand-in for evo2.Evo2 (v0.6.0 call signatures), for CPU tests of evo2_features.py.

Real API mirrored:
  Evo2(name).tokenizer.tokenize(str) -> list[int] (byte values)
  Evo2(name).model(input_ids)        -> (logits, inference_params)          # StripedHyena forward
  Evo2(name)(input_ids, return_embeddings=True, layer_names=[...])
                                     -> ((logits, inference_params), {layer_name: tensor (B, L, D)})
"""
import torch
import torch.nn as nn


class _Tok:
    def tokenize(self, s):
        return list(s.encode())


class _Net(nn.Module):
    def __init__(self, D=4096, V=512):
        super().__init__()
        g = torch.Generator().manual_seed(0)
        self.emb = nn.Parameter(torch.randn(V, 8, generator=g))
        self.conv = nn.Parameter(torch.randn(8, 8, 3, generator=g) * 0.3)
        self.out = nn.Parameter(torch.randn(8, V, generator=g))
        self.proj = nn.Parameter(torch.randn(8, D, generator=g) * 0.1)

    def hidden(self, ids):
        x = self.emb[ids].transpose(1, 2)                                  # (B, 8, L)
        x = torch.nn.functional.pad(x, (2, 0))                              # causal
        return torch.tanh(torch.nn.functional.conv1d(x, self.conv)).transpose(1, 2)  # (B, L, 8)

    def forward(self, ids):
        h = self.hidden(ids)
        return (h @ self.out, None)


class Evo2:
    def __init__(self, name="mock"):
        self.tokenizer = _Tok()
        self.model = _Net()

    def __call__(self, ids, return_embeddings=False, layer_names=None):
        out = self.model(ids)
        if not return_embeddings:
            return out, None
        h = self.model.hidden(ids) @ self.model.proj
        return out, {name: h for name in layer_names}
