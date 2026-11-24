#!/usr/bin/env python3
"""
Build the SynLet SQLite database (synlet.db) from the filtered datasets in cleaned/.
"""

import csv
import sqlite3
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
CLEAN_DIR = BASE_DIR / "cleaned"
SCHEMA_PATH = BASE_DIR / "sql" / "schema.sql"
DB_PATH = BASE_DIR / "synlet.db"

ORGANISMS = [
    (9606, "Homo sapiens", "Human"),
    (559292, "Saccharomyces cerevisiae (S288c)", "Yeast"),
    (7227, "Drosophila melanogaster", "Fruit fly"),
]


def dash_to_none(value: str):
    val = value.strip()
    return None if not val or val == "-" else val


def build_database():
    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")

    with SCHEMA_PATH.open("r", encoding="utf-8") as f:
        conn.executescript(f.read())

    # 1. Seed organisms
    conn.executemany(
        "INSERT INTO organism (organism_id, scientific_name, common_name) VALUES (?, ?, ?)",
        ORGANISMS,
    )

    # 2. Load G-HumanEssential.tsv
    ess_path = CLEAN_DIR / "G-HumanEssential.tsv"
    genes_ess = []
    pubs = {}
    ess_sources = []

    with ess_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            entrez_id = int(row["entrez_gene_id"] if "entrez_gene_id" in row else row["Gene ID"])
            ess_val = row["Essentiality (determined from multiple datasets)"].strip()
            pmids = [p.strip() for p in row["Data Source"].split(",") if p.strip()]
            genes_ess.append((entrez_id, 9606, ess_val, len(pmids)))
            for pmid in pmids:
                pubs.setdefault(pmid, None)
                ess_sources.append((entrez_id, pmid))

    conn.executemany(
        """
        INSERT INTO gene (entrez_id, organism_id, essentiality, essentiality_source_count)
        VALUES (?, ?, ?, ?)
        """,
        genes_ess,
    )
    conn.executemany(
        "INSERT OR IGNORE INTO publication (publication_id, author) VALUES (?, ?)",
        pubs.items(),
    )
    conn.executemany(
        "INSERT OR IGNORE INTO gene_essentiality_source (entrez_id, publication_id) VALUES (?, ?)",
        ess_sources,
    )

    # 3. Load BioGRID filtered dataset
    bg_path = CLEAN_DIR / "BIOGRID-SYSTEM-Synthetic_Lethality-5.0.263.filtered.tab3.txt"
    gene_upserts = {}
    aliases = set()
    bg_pubs = {}
    interactions = []

    alias_columns = [
        ("SYNONYM", "Synonyms Interactor"),
        ("SWISS-PROT", "SWISS-PROT Accessions Interactor"),
        ("TREMBL", "TREMBL Accessions Interactor"),
        ("REFSEQ", "REFSEQ Accessions Interactor"),
    ]

    with bg_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for side in ("A", "B"):
                entrez_id = int(row[f"Entrez Gene Interactor {side}"])
                org_id = int(row[f"Organism ID Interactor {side}"])
                biogrid_id = int(row[f"BioGRID ID Interactor {side}"])
                symbol = dash_to_none(row[f"Official Symbol Interactor {side}"])
                sys_name = dash_to_none(row[f"Systematic Name Interactor {side}"])

                gene_upserts[entrez_id] = (entrez_id, org_id, biogrid_id, symbol, sys_name)

                for alias_type, col_prefix in alias_columns:
                    raw = dash_to_none(row[f"{col_prefix} {side}"])
                    if raw:
                        for token in raw.split("|"):
                            t = token.strip()
                            if t and t != "-":
                                aliases.add((entrez_id, alias_type, t))

            pub_id = row["Publication Source"].strip()
            author = dash_to_none(row["Author"])
            if pub_id not in bg_pubs or bg_pubs[pub_id] is None:
                bg_pubs[pub_id] = author

            interactions.append(
                (
                    int(row["#BioGRID Interaction ID"]),
                    int(row["Entrez Gene Interactor A"]),
                    int(row["Entrez Gene Interactor B"]),
                    pub_id,
                    row["Throughput"].strip(),
                    dash_to_none(row["Score"]),
                    dash_to_none(row["Qualifications"]),
                    dash_to_none(row["Ontology Term Names"]),
                )
            )

    conn.executemany(
        """
        INSERT INTO gene (entrez_id, organism_id, biogrid_id, official_symbol, systematic_name)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(entrez_id) DO UPDATE SET
            biogrid_id = excluded.biogrid_id,
            official_symbol = COALESCE(gene.official_symbol, excluded.official_symbol),
            systematic_name = COALESCE(gene.systematic_name, excluded.systematic_name)
        """,
        gene_upserts.values(),
    )
    conn.executemany(
        "INSERT OR IGNORE INTO gene_alias (entrez_id, alias_type, alias_value) VALUES (?, ?, ?)",
        aliases,
    )
    conn.executemany(
        """
        INSERT INTO publication (publication_id, author) VALUES (?, ?)
        ON CONFLICT(publication_id) DO UPDATE SET
            author = COALESCE(publication.author, excluded.author)
        """,
        bg_pubs.items(),
    )
    conn.executemany(
        """
        INSERT INTO synthetic_lethality (
            biogrid_interaction_id, gene_a_entrez_id, gene_b_entrez_id,
            publication_id, throughput, score, qualifications, ontology_term_names
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        interactions,
    )

    # 4. Load OMA ortholog pairs (Human-Yeast & Human-DROME)
    oma_files = [
        ("OMA-Pairs_HUMAN-YEAST_2026-11-16_150430.csv", 559292),
        ("OMA-Pairs_HUMAN-DROME_2026-11-16_150443.csv", 7227),
    ]
    oma_genes = {}
    oma_map = {}
    orthology_rows = []

    for filename, model_org_id in oma_files:
        path = CLEAN_DIR / filename
        with path.open("r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                oma1 = row["OMAid_1"].strip()
                e1 = int(row["EntrezGene_1"].strip())
                oma2 = row["OMAid_2"].strip()
                e2 = int(row["EntrezGene_2"].strip())
                rel_type = row["RelType"].strip()
                oma_grp_raw = row["OmaGroup"].strip()
                oma_group = int(oma_grp_raw) if oma_grp_raw else None

                oma_genes[e1] = 9606
                oma_genes[e2] = model_org_id
                oma_map[oma1] = e1
                oma_map[oma2] = e2
                orthology_rows.append((oma1, oma2, rel_type, oma_group))

    conn.executemany(
        "INSERT OR IGNORE INTO gene (entrez_id, organism_id) VALUES (?, ?)",
        oma_genes.items(),
    )
    conn.executemany(
        "INSERT OR IGNORE INTO gene_oma (oma_id, entrez_id) VALUES (?, ?)",
        oma_map.items(),
    )
    conn.executemany(
        "INSERT INTO orthology (oma_id_1, oma_id_2, rel_type, oma_group) VALUES (?, ?, ?, ?)",
        orthology_rows,
    )

    conn.commit()

    # Verify foreign keys
    fk_errors = conn.execute("PRAGMA foreign_key_check;").fetchall()
    assert len(fk_errors) == 0, f"Foreign key check failed: {fk_errors[:5]}"

    # Print summary table
    tables = [
        "organism",
        "gene",
        "gene_oma",
        "gene_alias",
        "publication",
        "gene_essentiality_source",
        "orthology",
        "synthetic_lethality",
    ]
    print(f"Database created: {DB_PATH}")
    print("Foreign key check: OK (0 violations)\n")
    print(f"{'Table':<28} {'Rows':>10}")
    print("-" * 40)
    for t in tables:
        cnt = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"{t:<28} {cnt:>10}")

    conn.close()


if __name__ == "__main__":
    build_database()
