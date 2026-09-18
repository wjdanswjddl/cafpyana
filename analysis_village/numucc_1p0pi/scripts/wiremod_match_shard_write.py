#!/usr/bin/env python3
"""Shard WireMod matched writes across exclusive file partitions.

Each worker owns inputs where ``index % n_shards == shard_id``. Combined with
``skip_existing``, this prevents duplicate writes when the sequential writer is
stopped first. Does not start hist / Product A/B.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.scripts.dent_match_common_events import (
    list_df_files,
    matched_out_path,
    save_matched_sel_all_files,
)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--matched-out-dir", required=True, help="…/matched/<variation>")
    p.add_argument("--common-keys-pkl", required=True)
    p.add_argument("--variation-name", required=True)
    p.add_argument("--shard-id", type=int, required=True)
    p.add_argument("--n-shards", type=int, required=True)
    p.add_argument("--filename-str", default="sel_all")
    p.add_argument("--matched-suffix", default="_matched")
    p.add_argument(
        "--drop-partials",
        action="store_true",
        help="Remove matched outputs that lack a .source sidecar (interrupted writes).",
    )
    args = p.parse_args(argv)

    if not (0 <= args.shard_id < args.n_shards):
        raise SystemExit(f"shard-id must be in [0, {args.n_shards})")

    keys_path = Path(args.common_keys_pkl)
    with open(keys_path, "rb") as fh:
        common_keys = pickle.load(fh)
    print(
        f"[shard {args.shard_id}/{args.n_shards}] loaded {len(common_keys)} keys "
        f"from {keys_path}",
        flush=True,
    )

    files = list_df_files(args.raw_dir, args.filename_str)
    files = sorted(files)
    mine = [f for i, f in enumerate(files) if i % args.n_shards == args.shard_id]
    print(f"[shard {args.shard_id}] {len(mine)}/{len(files)} assigned inputs", flush=True)

    out_dir = args.matched_out_dir
    os.makedirs(out_dir, exist_ok=True)

    if args.drop_partials:
        n_drop = 0
        for fpath in mine:
            out_path = matched_out_path(fpath, suffix=args.matched_suffix, out_dir=out_dir)
            if os.path.isfile(out_path) and not os.path.isfile(out_path + ".source"):
                os.remove(out_path)
                n_drop += 1
                print(f"[shard {args.shard_id}] dropped partial {out_path}", flush=True)
        print(f"[shard {args.shard_id}] dropped {n_drop} partials", flush=True)

    summary = save_matched_sel_all_files(
        mine,
        common_keys,
        matched_suffix=args.matched_suffix,
        variation_name=f"{args.variation_name}_s{args.shard_id}",
        out_dir=out_dir,
        skip_existing=True,
    )
    n_written = 0 if summary is None or len(summary) == 0 else len(summary)
    print(
        f"[shard {args.shard_id}] done written={n_written} "
        f"rss_hint=/proc/{os.getpid()}/status",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
