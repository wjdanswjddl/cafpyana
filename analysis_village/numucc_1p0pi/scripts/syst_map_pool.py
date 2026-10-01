"""Forked worker pool shared by the systematic map dispatchers.

Each dispatcher still builds its own job list and log lines. This module only
runs ``imap_unordered`` under ``fork`` with the same ``maxtasksperchild`` the
individual pools used.
"""
from __future__ import annotations

import multiprocessing as mp
from typing import Any, Iterable, Iterator, Optional


def fork_imap(
    worker,
    payloads: Iterable[Any],
    *,
    processes: int,
    maxtasksperchild: Optional[int] = 32,
) -> Iterator[Any]:
    """Yield worker results. The pool is closed when iteration finishes."""
    ctx = mp.get_context("fork")
    pool_kw = {"processes": processes}
    if maxtasksperchild is not None:
        pool_kw["maxtasksperchild"] = maxtasksperchild
    pool = ctx.Pool(**pool_kw)
    try:
        yield from pool.imap_unordered(worker, list(payloads), chunksize=1)
    finally:
        pool.terminate()
        pool.join()
