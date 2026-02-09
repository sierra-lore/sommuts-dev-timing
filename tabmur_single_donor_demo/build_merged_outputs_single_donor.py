#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

SUFFIX_GENO = ".single_cell_genotype.tsv"
SUFFIX_SITES = ".SitesPerCell.tsv"
BASE_DIR = Path(__file__).resolve().parent
GENO_KEEP_COLS = [
    "#CHROM",
    "Start",
    "End",
    "REF",
    "ALT_expected",
    "CB",
    "Cell_type_observed",
    "Base_observed",
    "Num_reads",
]


def norm_cb(s):
    return s.astype(str).str.strip().str.upper()


def parse_celltype_from_geno(name):
    if not name.endswith(SUFFIX_GENO):
        return None
    return name[: -len(SUFFIX_GENO)]


def parse_celltype_from_sites(donor, name):
    if not name.endswith(SUFFIX_SITES):
        return None
    base = name[: -len(SUFFIX_SITES)]
    prefix = f"{donor}_"
    if not base.startswith(prefix):
        return None
    parts = base.rsplit(".", 1)
    if len(parts) != 2:
        return None
    return parts[1]


def read_genotype(fp):
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    if "CB" not in df.columns:
        return None
    df["CB"] = norm_cb(df["CB"])
    return df


def read_sites_per_cell(fp):
    df = pd.read_csv(fp, sep=",", low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    if "CB" not in df.columns:
        return None
    if "SitesPerCell" not in df.columns:
        return pd.DataFrame(columns=["CB", "SitesPerCell"])
    df["CB"] = norm_cb(df["CB"])
    return df[["CB", "SitesPerCell"]]


def collect_geno_files_for_donor(root, donor):
    base = root / donor
    ct_map = {}
    if not base.is_dir():
        return ct_map
    for lib_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        for gfp in lib_dir.glob(f"*{SUFFIX_GENO}"):
            ct = parse_celltype_from_geno(gfp.name)
            if ct:
                ct_map.setdefault(ct, []).append(gfp)
    return ct_map


def collect_sites_files_for_donor(root, donor):
    base = root / donor
    ct_map = {}
    if not base.is_dir():
        return ct_map
    for lib_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        for sfp in lib_dir.glob(f"*{SUFFIX_SITES}"):
            ct = parse_celltype_from_sites(donor, sfp.name)
            if ct:
                ct_map.setdefault(ct, []).append(sfp)
    return ct_map


def merge_and_write_for_donor(donor, geno_root, sites_root, dest_root):
    geno_map = collect_geno_files_for_donor(geno_root, donor)
    sites_map = collect_sites_files_for_donor(sites_root, donor)

    if not geno_map and not sites_map:
        raise SystemExit(f"No files found for donor {donor}")

    all_cts = sorted(set(geno_map.keys()) | set(sites_map.keys()))
    dest_geno_dir = dest_root / donor
    dest_geno_dir.mkdir(parents=True, exist_ok=True)

    geno_written = 0
    sites_written = 0

    for ct in all_cts:
        gfiles = geno_map.get(ct, [])
        sfiles = sites_map.get(ct, [])

        if gfiles:
            gdfs = [read_genotype(fp) for fp in gfiles]
            gdfs = [df for df in gdfs if df is not None]
            if gdfs:
                g_all = pd.concat(gdfs, ignore_index=True)
                if all(c in g_all.columns for c in GENO_KEEP_COLS):
                    g_all = g_all[GENO_KEEP_COLS]
                    out_g = dest_geno_dir / f"{ct}.single_cell_genotype.tsv"
                    g_all.to_csv(out_g, sep="\t", index=False)
                    geno_written += 1

        if sfiles:
            sdfs = [read_sites_per_cell(fp) for fp in sfiles]
            sdfs = [df for df in sdfs if df is not None]
            if sdfs:
                s_all = pd.concat(sdfs, ignore_index=True)
                s_all = s_all.groupby("CB", as_index=False)["SitesPerCell"].max()
                out_s = dest_geno_dir / f"{donor}.{ct}.SitesPerCell.tsv"
                s_all.to_csv(out_s, sep="\t", index=False)
                sites_written += 1

    if geno_written == 0 and sites_written == 0:
        raise SystemExit(f"No merged files written for donor {donor}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor", default="1-M-62")
    parser.add_argument(
        "--geno-root",
        type=Path,
        default=Path("/home/ubuntu/sommuts/TabMur/outputs/1mo_outputs/SingleCellAlleles"),
    )
    parser.add_argument(
        "--sites-root",
        type=Path,
        default=Path("/home/ubuntu/sommuts/TabMur/outputs/1mo_outputs/UniqueCellCallableSites"),
    )
    parser.add_argument(
        "--dest-root",
        type=Path,
        default=BASE_DIR / "demo_data" / "merged_outputs_demo",
    )
    args = parser.parse_args()

    if not args.geno_root.exists():
        raise SystemExit(f"geno-root not found: {args.geno_root}")
    if not args.sites_root.exists():
        raise SystemExit(f"sites-root not found: {args.sites_root}")

    args.dest_root.mkdir(parents=True, exist_ok=True)
    merge_and_write_for_donor(args.donor, args.geno_root, args.sites_root, args.dest_root)


if __name__ == "__main__":
    main()
