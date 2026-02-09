#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

GENO_SUFFIX = ".single_cell_genotype.tsv"
BASE_DIR = Path(__file__).resolve().parent
GENO_SUB_COLS = [
    "#CHROM",
    "Start",
    "End",
    "REF",
    "ALT_expected",
    "CB",
    "Cell_type_observed",
    "Num_reads",
]
GENO_OUT_COLS = [
    "#CHROM",
    "Start",
    "End",
    "REF",
    "ALT_expected",
    "donor",
    "CB",
    "Cell_type_observed",
    "Num_reads",
]


def norm_cb(s):
    return s.astype(str).str.strip().str.upper()


def mode_one(s):
    s = s.dropna()
    vc = s.value_counts(dropna=False)
    return vc.index[0] if not vc.empty else pd.NA


def load_genotypes(outputs_root, donor_filter=None):
    rows = []
    for donor_dir in sorted(p for p in outputs_root.iterdir() if p.is_dir()):
        donor = donor_dir.name
        if donor_filter and donor not in donor_filter:
            continue
        for fp in sorted(donor_dir.glob(f"*{GENO_SUFFIX}")):
            df = pd.read_csv(fp, sep="\t", low_memory=False)
            df["CB"] = norm_cb(df["CB"])
            mask = df["ALT_expected"].astype(str) == df["Base_observed"].astype(str)
            sub = df.loc[mask, GENO_SUB_COLS].copy()
            if sub.empty:
                continue
            sub["donor"] = donor
            rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=GENO_OUT_COLS)
    out = pd.concat(rows, ignore_index=True)
    return out[GENO_OUT_COLS]


def add_tissue(cb_table, barcode_path):
    bar = pd.read_csv(barcode_path)
    bar = bar.rename(columns={"mouse.id": "donor"})
    bar["CB"] = norm_cb(bar["CB"])
    cb_table["CB"] = norm_cb(cb_table["CB"])
    collapsed = bar.groupby(["donor", "CB"], as_index=False)["tissue"].agg(mode_one)
    return cb_table.merge(collapsed, on=["donor", "CB"], how="left")


def add_germ_layer(cb_table, germ_map_path):
    gmap = pd.read_csv(germ_map_path)
    if "germ_layer" in cb_table.columns:
        cb_table = cb_table.drop(columns=["germ_layer"])
    gmap = gmap.rename(columns={"cell_type": "Cell_type_observed", "tissue_type": "tissue"})
    return cb_table.merge(
        gmap[["Cell_type_observed", "tissue", "germ_layer"]],
        on=["Cell_type_observed", "tissue"],
        how="left",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outputs-root",
        type=Path,
        default=BASE_DIR / "demo_data" / "merged_outputs_demo",
    )
    parser.add_argument("--donor", default="1-M-62")
    parser.add_argument(
        "--cellbarcodes",
        type=Path,
        default=BASE_DIR / "demo_data" / "cellbarcodes_with_tissue.csv",
    )
    parser.add_argument(
        "--germ-map",
        type=Path,
        default=BASE_DIR / "demo_data" / "cell_type_to_germ_layer_map.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=BASE_DIR / "demo_data" / "cell_metadata_table_1-M-62.csv",
    )
    args = parser.parse_args()

    if not args.outputs_root.exists():
        raise SystemExit(f"outputs-root not found: {args.outputs_root}")

    cb_table = load_genotypes(args.outputs_root, donor_filter={args.donor})
    if cb_table.empty:
        raise SystemExit("No ALT matches found in genotype files.")

    cb_table = add_tissue(cb_table, args.cellbarcodes)
    cb_table = add_germ_layer(cb_table, args.germ_map)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cb_table.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
