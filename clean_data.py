#!/usr/bin/env python3
"""
Clean and filter the raw SynLet datasets for our 3 target organisms:
  - Human: Homo sapiens (Taxon ID: 9606)
  - Yeast: Saccharomyces cerevisiae (S288c) (Taxon ID: 559292)
  - Fruit fly (DROME): Drosophila melanogaster (Taxon ID: 7227)

Cleaning performed (raw files are left untouched):
  1. BioGRID:
     - Keep only rows where BOTH interactors belong to {9606, 559292, 7227}
       (including cross-species pairs among these 3 organisms).
     - Remove self-interactions (Interactor A == Interactor B in same organism).
     - Remove 100%-empty columns: Modification, Tags.
     - Remove 100%-identical/constant columns: Experimental System,
       Experimental System Type, Source Database.
  2. G-HumanEssential:
     - Copied untouched (already Human-only).
  3. OMA (HUMAN-YEAST & HUMAN-DROME):
     - Remove rows missing EntrezGene_1 or EntrezGene_2 so every ortholog pair
       can be linked via Entrez Gene ID.
"""

import csv
import shutil
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "cleaned"
OUT_DIR.mkdir(exist_ok=True)

TARGET_ORGS = {"9606", "559292", "7227"}

BIOGRID_DROP_COLS = {
    "Modification",             # 100% empty ('-')
    "Tags",                     # 100% empty ('-')
    "Experimental System",      # 100% identical ('Synthetic Lethality')
    "Experimental System Type", # 100% identical ('genetic')
    "Source Database",          # 100% identical ('BIOGRID')
}


def clean_biogrid():
    src = BASE_DIR / "BIOGRID-SYSTEM-Synthetic_Lethality-5.0.263.tab3.txt"
    dst = OUT_DIR / "BIOGRID-SYSTEM-Synthetic_Lethality-5.0.263.filtered.tab3.txt"

    total = 0
    kept = 0
    dropped_other_org = 0
    dropped_self_loop = 0

    with src.open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        out_fields = [c for c in reader.fieldnames if c not in BIOGRID_DROP_COLS]

        with dst.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()

            for row in reader:
                total += 1
                oa = row["Organism ID Interactor A"]
                ob = row["Organism ID Interactor B"]

                if oa not in TARGET_ORGS or ob not in TARGET_ORGS:
                    dropped_other_org += 1
                    continue

                if (
                    oa == ob
                    and row["Entrez Gene Interactor A"] == row["Entrez Gene Interactor B"]
                ):
                    dropped_self_loop += 1
                    continue

                writer.writerow({col: row[col] for col in out_fields})
                kept += 1

    print(
        f"[BioGRID] total={total}, kept={kept}, "
        f"dropped_other_org={dropped_other_org}, dropped_self_loop={dropped_self_loop}, "
        f"columns={len(out_fields)} (dropped {len(BIOGRID_DROP_COLS)})"
    )


def clean_human_essential():
    src = BASE_DIR / "G-HumanEssential.tsv"
    dst = OUT_DIR / "G-HumanEssential.tsv"
    shutil.copy2(src, dst)
    print(f"[HumanEssential] copied untouched -> {dst.name}")


def clean_oma(filename: str):
    src = BASE_DIR / filename
    dst = OUT_DIR / filename

    total = 0
    kept = 0
    dropped_missing_entrez = 0

    with src.open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin)
        out_fields = reader.fieldnames

        with dst.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, lineterminator="\n")
            writer.writeheader()

            for row in reader:
                total += 1
                if not row["EntrezGene_1"].strip() or not row["EntrezGene_2"].strip():
                    dropped_missing_entrez += 1
                    continue
                writer.writerow(row)
                kept += 1

    print(
        f"[OMA {filename}] total={total}, kept={kept}, "
        f"dropped_missing_entrez={dropped_missing_entrez}"
    )


if __name__ == "__main__":
    for old in OUT_DIR.iterdir():
        if old.is_file():
            old.unlink()

    clean_biogrid()
    clean_human_essential()
    clean_oma("OMA-Pairs_HUMAN-YEAST_2026-11-16_150430.csv")
    clean_oma("OMA-Pairs_HUMAN-DROME_2026-11-16_150443.csv")
    print(f"\nCleaned files written to: {OUT_DIR}")
