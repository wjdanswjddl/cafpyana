"""Shared WireMod walk helpers (RSS, universe merge, CV merge, shard split)."""
from __future__ import annotations

import os
from typing import Optional, Sequence

import numpy as np


def rss_gb() -> float:
    try:
        with open(f"/proc/{os.getpid()}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1]) / (1024.0 * 1024.0)
    except Exception:
        pass
    return float("nan")


def shard_files(files: Sequence[str], n_shards: int, shard_id: int) -> list:
    """Files owned by ``shard_id`` when workers split ``i % n_shards``."""
    files = list(files)
    if n_shards <= 1:
        return files
    return [f for i, f in enumerate(files) if i % n_shards == shard_id]


def merge_univ_products(acc: dict, chunk: dict) -> dict:
    """Add chunk universe hists into the accumulator."""
    if not acc:
        return {
            "by_universe": {
                u: {
                    "hists_cut": {k: np.asarray(v, dtype=float).copy() for k, v in p["hists_cut"].items()},
                    "hists_final": {k: np.asarray(v, dtype=float).copy() for k, v in p["hists_final"].items()},
                }
                for u, p in chunk["by_universe"].items()
            },
            "pot": float(chunk["pot"]),
            "cut_var_names": list(chunk["cut_var_names"]),
            "final_var_names": list(chunk["final_var_names"]),
            "universes": list(chunk["universes"]),
        }
    acc["pot"] = float(acc["pot"]) + float(chunk["pot"])
    for u, payload in chunk["by_universe"].items():
        dst = acc["by_universe"].setdefault(
            u,
            {
                "hists_cut": {k: np.zeros_like(v, dtype=float) for k, v in payload["hists_cut"].items()},
                "hists_final": {k: np.zeros_like(v, dtype=float) for k, v in payload["hists_final"].items()},
            },
        )
        for key in ("hists_cut", "hists_final"):
            for var, hist in payload[key].items():
                dst[key][var] = np.asarray(dst[key].get(var, 0.0), dtype=float) + np.asarray(
                    hist, dtype=float
                )
    return acc


def merge_cv_products(acc: dict, chunk: dict) -> dict:
    """Add a matched-CV hist chunk (no calorimetry universes)."""
    if not acc:
        return {
            "hists_cut": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_cut"].items()},
            "hists_final": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_final"].items()},
            "pot": float(chunk["pot"]),
            "cut_var_names": list(chunk["cut_var_names"]),
            "final_var_names": list(chunk["final_var_names"]),
        }
    acc["pot"] = float(acc["pot"]) + float(chunk["pot"])
    for key in ("hists_cut", "hists_final"):
        for var, hist in chunk[key].items():
            acc[key][var] = np.asarray(acc[key].get(var, 0.0), dtype=float) + np.asarray(hist, dtype=float)
    return acc
