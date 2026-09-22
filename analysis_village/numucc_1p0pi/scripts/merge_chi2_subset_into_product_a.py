#!/usr/bin/env python
"""Merge chi2 track-subset packs into Product A **without overwriting existing keys**.

Adds only ``CHI2_TRACK_SUBSET_SLUGS`` from campaign outputs into::

    PRL/systematics/productA_sel_all/{Flux,G4,Cosmics,Detector,GENIE}/

Updates manifests with an append-only ``chi2_track_subset`` section.
"""
from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
from datetime import datetime, timezone
from os import path
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.syst_pipeline_walker import CHI2_TRACK_SUBSET_SLUGS

PRL_A = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/systematics/productA_sel_all"
)


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _backup(p: Path) -> Path:
    bak = p.with_suffix(p.suffix + f".bak_before_chi2subset_{datetime.now():%Y%m%d_%H%M%S}")
    shutil.copy2(p, bak)
    return bak


def _merge_npz_var_blocks(
    dest: Path,
    src_by_var: Mapping[str, Any],
    *,
    dry_run: bool,
) -> list[str]:
    """Merge top-level var→pack entries (Flux/G4/Cosmics style)."""
    z = dict(np.load(dest, allow_pickle=True))
    added = []
    for slug in CHI2_TRACK_SUBSET_SLUGS:
        if slug not in src_by_var:
            continue
        if slug in z:
            print(f"[skip] {dest.name} already has {slug}", flush=True)
            continue
        z[slug] = src_by_var[slug]
        added.append(slug)
    if not added:
        print(f"[merge] nothing new for {dest}", flush=True)
        return []
    if dry_run:
        print(f"[dry-run] would add to {dest}: {added}", flush=True)
        return added
    _backup(dest)
    np.savez_compressed(dest, **z)
    print(f"[merge] added {added} → {dest}", flush=True)
    return added


def _merge_flux_g4(family: str, src_npz: Path, *, dry_run: bool) -> list[str]:
    dest = PRL_A / family / (
        "flux_syst_dict.npz" if family == "Flux" else "g4_syst_dict.npz"
    )
    zsrc = np.load(src_npz, allow_pickle=True)
    # Support both flat {slug: pack} and nested {family: {slug: pack}}
    if family in zsrc.files:
        block = zsrc[family].item()
        src_by_var = {k: block[k] for k in CHI2_TRACK_SUBSET_SLUGS if k in block}
        # also by_knob if present
        bk_key = f"{family}_by_knob"
        dest_z = dict(np.load(dest, allow_pickle=True))
        added = []
        for slug in CHI2_TRACK_SUBSET_SLUGS:
            if slug not in block:
                continue
            if slug in dest_z:
                print(f"[skip] {dest.name} already has {slug}", flush=True)
                continue
            dest_z[slug] = {
                ("flux" if family == "Flux" else "G4"): block[slug]
                if "cov_frac" in block[slug]
                else block[slug],
            }
            # normalize: Product A Flux stores {slug: {flux: pack, flux_by_knob?: ...}}
            pack = block[slug]
            if isinstance(pack, dict) and "cov_frac" in pack:
                cell = {("flux" if family == "Flux" else "G4"): pack}
            else:
                cell = pack
            # Prefer matching existing Product A layout
            existing_sample = next(iter(dest_z.values()))
            if isinstance(existing_sample, np.ndarray):
                existing_sample = existing_sample.item()
            comb = "flux" if family == "Flux" else "G4"
            if comb in existing_sample:
                cell = {comb: pack if "cov_frac" in pack else pack.get(comb, pack)}
                if bk_key.replace("_by_knob", "") :  # noqa
                    pass
            dest_z[slug] = cell
            added.append(slug)
        # Fix properly by inspecting existing layout once
        return _merge_flux_g4_layout(dest, src_npz, family, dry_run=dry_run)
    # flat files: each key is a var
    src_by_var = {k: zsrc[k] for k in CHI2_TRACK_SUBSET_SLUGS if k in zsrc.files}
    return _merge_npz_var_blocks(dest, src_by_var, dry_run=dry_run)


def _merge_flux_g4_layout(
    dest: Path, src_npz: Path, family: str, *, dry_run: bool
) -> list[str]:
    comb = "flux" if family == "Flux" else "G4"
    by_knob_key = "flux_by_knob" if family == "Flux" else "G4_by_knob"
    zsrc = np.load(src_npz, allow_pickle=True)
    dest_z = dict(np.load(dest, allow_pickle=True))

    # Discover src layout
    if family in zsrc.files:
        fam_block = zsrc[family].item()
        bk_block = zsrc[f"{family}_by_knob"].item() if f"{family}_by_knob" in zsrc.files else {}
    else:
        # histcounts productA layout: each var is {flux: pack, flux_by_knob: {...}}
        fam_block = {}
        bk_block = {}
        for slug in CHI2_TRACK_SUBSET_SLUGS:
            if slug not in zsrc.files:
                continue
            cell = zsrc[slug].item()
            if comb in cell:
                fam_block[slug] = cell[comb]
                if by_knob_key in cell:
                    bk_block[slug] = cell[by_knob_key]
            elif "cov_frac" in cell:
                fam_block[slug] = cell

    added = []
    for slug in CHI2_TRACK_SUBSET_SLUGS:
        if slug not in fam_block:
            print(f"[warn] {family} source missing {slug}", flush=True)
            continue
        if slug in dest_z:
            print(f"[skip] {dest.name} already has {slug}", flush=True)
            continue
        cell = {comb: fam_block[slug]}
        if slug in bk_block:
            cell[by_knob_key] = bk_block[slug]
        dest_z[slug] = cell
        added.append(slug)

    if not added:
        return []
    if dry_run:
        print(f"[dry-run] would add to {dest}: {added}", flush=True)
        return added
    _backup(dest)
    # np.savez needs arrays; store dicts via object arrays like existing
    out = {}
    for k, v in dest_z.items():
        out[k] = v
    np.savez_compressed(dest, **out)
    print(f"[merge] added {added} → {dest}", flush=True)
    _append_manifest(
        dest.parent / "manifest.json",
        family,
        added,
        source=str(src_npz),
    )
    return added


def _append_manifest(
    mani_path: Path,
    component: str,
    added: Sequence[str],
    *,
    source: str,
) -> None:
    if not mani_path.is_file():
        data: Dict[str, Any] = {}
    else:
        data = json.loads(mani_path.read_text())
    section = data.setdefault("chi2_track_subset", {})
    section["updated"] = _now()
    section["source"] = source
    section["added_slugs"] = list(added)
    section["note"] = (
        "Append-only merge of len>50 / not_mu avg-χ² @ 2prong-vtxdist; "
        "existing all-track packs unchanged."
    )
    # also update status maps if present
    status = data.get("status_by_var")
    if isinstance(status, dict):
        for s in added:
            status[s] = "filled_chi2_track_subset"
    nb = data.get("target_nbins")
    if isinstance(nb, dict):
        for s in added:
            nb[s] = 60
    data["n_vars"] = int(data.get("n_vars", 0)) + len(
        [s for s in added if s not in (data.get("_counted_subset") or [])]
    )
    mani_path.write_text(json.dumps(data, indent=2) + "\n")
    print(f"[manifest] updated {mani_path}", flush=True)


def merge_cosmics(src_npz: Path, *, dry_run: bool) -> list[str]:
    dest = PRL_A / "Cosmics" / "cosmics_syst_dict.npz"
    zsrc = np.load(src_npz, allow_pickle=True)
    src = {k: zsrc[k] for k in CHI2_TRACK_SUBSET_SLUGS if k in zsrc.files}
    added = _merge_npz_var_blocks(dest, src, dry_run=dry_run)
    if added and not dry_run:
        mani = PRL_A / "Cosmics" / "cosmics_covariance_manifest.json"
        _append_manifest(mani, "Cosmics", added, source=str(src_npz))
        # keep n_vars in sync
        if mani.is_file():
            d = json.loads(mani.read_text())
            vars_ = d.get("variables")
            if isinstance(vars_, list):
                for s in added:
                    if s not in vars_:
                        vars_.append(s)
                d["variables"] = vars_
                d["n_vars"] = len(vars_)
                mani.write_text(json.dumps(d, indent=2) + "\n")
    return added


def merge_detector(src_npz: Path, *, dry_run: bool) -> list[str]:
    dest = PRL_A / "Detector" / "detector_sel_syst_dict.npz"
    zsrc = np.load(src_npz, allow_pickle=True)
    src_det = zsrc["detector"].item()
    dest_z = dict(np.load(dest, allow_pickle=True))
    dest_det = dest_z["detector"].item()
    added = []
    for slug in CHI2_TRACK_SUBSET_SLUGS:
        if slug not in src_det:
            print(f"[warn] Detector source missing {slug}", flush=True)
            continue
        if slug in dest_det:
            print(f"[skip] Detector already has {slug}", flush=True)
            continue
        dest_det[slug] = src_det[slug]
        added.append(slug)
        # per-knob blocks
        for key in dest_z:
            if not key.startswith("detector-"):
                continue
            src_knob = zsrc[key].item() if key in zsrc.files else {}
            if slug in src_knob:
                blk = dest_z[key].item()
                blk[slug] = src_knob[slug]
                dest_z[key] = blk
        if "detector_by_wiremod" in dest_z.files if False else "detector_by_wiremod" in dest_z:
            by = dest_z["detector_by_wiremod"].item()
            src_by = (
                zsrc["detector_by_wiremod"].item()
                if "detector_by_wiremod" in zsrc.files
                else {}
            )
            if slug in src_by:
                by[slug] = src_by[slug]
            dest_z["detector_by_wiremod"] = by

    if not added:
        return []
    dest_z["detector"] = dest_det
    if dry_run:
        print(f"[dry-run] would add Detector: {added}", flush=True)
        return added
    _backup(dest)
    np.savez_compressed(dest, **dest_z)
    print(f"[merge] added {added} → {dest}", flush=True)
    _append_manifest(
        PRL_A / "Detector" / "manifest.json",
        "Detector",
        added,
        source=str(src_npz),
    )
    mani = PRL_A / "Detector" / "manifest.json"
    if mani.is_file():
        d = json.loads(mani.read_text())
        vars_ = d.get("variables")
        if isinstance(vars_, list):
            for s in added:
                if s not in vars_:
                    vars_.append(s)
            d["variables"] = sorted(vars_)
            d["n_vars"] = len(vars_)
            mani.write_text(json.dumps(d, indent=2) + "\n")
    return added


def merge_genie(src_cov_pkl: Path, src_slim_npz: Optional[Path], *, dry_run: bool) -> list[str]:
    dest_pkl = PRL_A / "GENIE" / "cov_mat_dict.pkl"
    with open(src_cov_pkl, "rb") as f:
        src = pickle.load(f)
    with open(dest_pkl, "rb") as f:
        dest = pickle.load(f)
    added = []
    for slug in CHI2_TRACK_SUBSET_SLUGS:
        if slug not in src:
            print(f"[warn] GENIE source missing {slug}", flush=True)
            continue
        if slug in dest:
            print(f"[skip] GENIE cov_mat already has {slug}", flush=True)
            continue
        row = src[slug]
        # Prefer genie_rate / GENIE_slim_v3_rate → store as genie_rate only (match Product A)
        if "genie_rate" in row:
            dest[slug] = {"genie_rate": row["genie_rate"]}
        elif "GENIE_slim_v3_rate" in row:
            dest[slug] = {"genie_rate": row["GENIE_slim_v3_rate"]}
        elif "rate" in row and isinstance(row["rate"], dict) and "cov_frac" in row["rate"]:
            dest[slug] = {"genie_rate": row["rate"]["cov_frac"]}
        else:
            # FSI pack style may nest under GENIE_slim_v3
            print(f"[warn] GENIE {slug} unexpected keys {list(row)[:8]}", flush=True)
            continue
        added.append(slug)

    if added and not dry_run:
        _backup(dest_pkl)
        with open(dest_pkl, "wb") as f:
            pickle.dump(dest, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"[merge] added {added} → {dest_pkl}", flush=True)

        if src_slim_npz and src_slim_npz.is_file():
            dest_slim = PRL_A / "GENIE" / "genie_slim_v3.npz"
            zs = np.load(src_slim_npz, allow_pickle=True)
            zd = dict(np.load(dest_slim, allow_pickle=True))
            syst = zd["syst"].item()
            pack = syst.setdefault("GENIE_slim_v3", {})
            src_syst = zs["syst"].item()
            src_pack = src_syst.get("GENIE_slim_v3") or src_syst.get("syst", {}).get(
                "GENIE_slim_v3"
            )
            if src_pack is None and "GENIE_slim_v3" in src_syst:
                src_pack = src_syst["GENIE_slim_v3"]
            # also accept flat by-var under syst
            for slug in added:
                if src_pack and slug in src_pack:
                    pack[slug] = src_pack[slug]
                elif slug in src_syst:
                    pack[slug] = src_syst[slug]
            zd["syst"] = {"GENIE_slim_v3": pack}
            _backup(dest_slim)
            np.savez_compressed(dest_slim, **zd)
            print(f"[merge] updated slim_v3 for {added}", flush=True)

        _append_manifest(
            PRL_A / "GENIE" / "manifest.json",
            "GENIE",
            added,
            source=str(src_cov_pkl),
        )
    elif dry_run:
        print(f"[dry-run] would add GENIE: {added}", flush=True)
    return added


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flux-npz", type=Path, default=None)
    p.add_argument("--g4-npz", type=Path, default=None)
    p.add_argument("--cosmics-npz", type=Path, default=None)
    p.add_argument("--detector-npz", type=Path, default=None)
    p.add_argument("--genie-cov-pkl", type=Path, default=None)
    p.add_argument("--genie-slim-npz", type=Path, default=None)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.flux_npz:
        _merge_flux_g4_layout(
            PRL_A / "Flux" / "flux_syst_dict.npz",
            args.flux_npz,
            "Flux",
            dry_run=args.dry_run,
        )
    if args.g4_npz:
        _merge_flux_g4_layout(
            PRL_A / "G4" / "g4_syst_dict.npz",
            args.g4_npz,
            "G4",
            dry_run=args.dry_run,
        )
    if args.cosmics_npz:
        merge_cosmics(args.cosmics_npz, dry_run=args.dry_run)
    if args.detector_npz:
        merge_detector(args.detector_npz, dry_run=args.dry_run)
    if args.genie_cov_pkl:
        merge_genie(args.genie_cov_pkl, args.genie_slim_npz, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
