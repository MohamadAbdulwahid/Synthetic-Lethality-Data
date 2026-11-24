#!/usr/bin/env python3
"""
SynLet CLI Search Tool
Look up ranked synthetic lethal partner genes for a mutated human gene using
Entrez Gene ID, BioGRID ID, or OMA ID.
"""

import argparse
import sqlite3
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DB = BASE_DIR / "synlet.db"

LOOKUP_SQL = """
WITH evidence AS (
    SELECT
        g2.entrez_id AS partner_entrez_id,
        'Direct Human' AS path_type,
        50 AS directness_pts,
        9606 AS evidence_org_id,
        sl.biogrid_interaction_id,
        sl.publication_id,
        sl.throughput,
        0 AS rel_pts,
        0 AS has_oma_group
    FROM v_sl_undirected sl
    JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
    WHERE sl.gene_1_entrez_id = :input_entrez_id
      AND g2.organism_id = 9606
      AND g2.entrez_id <> :input_entrez_id

    UNION ALL

    SELECT
        o2.human_entrez_id AS partner_entrez_id,
        'Cross-Species' AS path_type,
        30 AS directness_pts,
        gm.organism_id AS evidence_org_id,
        sl.biogrid_interaction_id,
        sl.publication_id,
        sl.throughput,
        o2.rel_pts,
        o2.has_oma_group
    FROM v_sl_undirected sl
    JOIN gene gm ON gm.entrez_id = sl.gene_2_entrez_id
    JOIN v_gene_orthology o2 ON o2.model_entrez_id = gm.entrez_id
    WHERE sl.gene_1_entrez_id = :input_entrez_id
      AND gm.organism_id <> 9606
      AND o2.human_entrez_id <> :input_entrez_id

    UNION ALL

    SELECT
        g2.entrez_id AS partner_entrez_id,
        'Ortholog' AS path_type,
        15 AS directness_pts,
        o1.model_organism_id AS evidence_org_id,
        sl.biogrid_interaction_id,
        sl.publication_id,
        sl.throughput,
        o1.rel_pts,
        o1.has_oma_group
    FROM v_gene_orthology o1
    JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
    JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
    WHERE o1.human_entrez_id = :input_entrez_id
      AND g2.organism_id = 9606
      AND g2.entrez_id <> :input_entrez_id

    UNION ALL

    SELECT
        o2.human_entrez_id AS partner_entrez_id,
        'Ortholog' AS path_type,
        15 AS directness_pts,
        o1.model_organism_id AS evidence_org_id,
        sl.biogrid_interaction_id,
        sl.publication_id,
        sl.throughput,
        MIN(o1.rel_pts, o2.rel_pts) AS rel_pts,
        MIN(o1.has_oma_group, o2.has_oma_group) AS has_oma_group
    FROM v_gene_orthology o1
    JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
    JOIN v_gene_orthology o2 ON o2.model_entrez_id = sl.gene_2_entrez_id
    WHERE o1.human_entrez_id = :input_entrez_id
      AND o2.human_entrez_id <> :input_entrez_id
),
scored AS (
    SELECT
        e.partner_entrez_id,
        g.biogrid_id,
        COALESCE(g.official_symbol, '-') AS official_symbol,
        g.essentiality,
        CASE g.essentiality
            WHEN 'Non-essential' THEN 1
            WHEN 'Unknown' THEN 2
            ELSE 3
        END AS essentiality_rank,
        GROUP_CONCAT(DISTINCT e.path_type) AS evidence_paths,
        GROUP_CONCAT(DISTINCT org.common_name) AS evidence_organisms,
        COUNT(DISTINCT e.publication_id) AS study_count,
        COUNT(DISTINCT e.biogrid_interaction_id) AS interaction_count,
        CASE
            WHEN MAX(CASE WHEN e.throughput LIKE '%Low Throughput%' THEN 1 ELSE 0 END) = 1
             AND MAX(CASE WHEN e.throughput LIKE '%High Throughput%' THEN 1 ELSE 0 END) = 1
                THEN 'Both'
            WHEN MAX(CASE WHEN e.throughput LIKE '%Low Throughput%' THEN 1 ELSE 0 END) = 1
                THEN 'Low'
            ELSE 'High'
        END AS throughput_summary,
        MAX(e.directness_pts) AS directness_pts,
        MIN(COUNT(DISTINCT e.publication_id) * 10, 50) AS study_pts,
        CASE WHEN COUNT(DISTINCT e.evidence_org_id) > 1 THEN 20 ELSE 0 END AS multi_species_pts,
        CASE
            WHEN MAX(CASE WHEN e.throughput LIKE '%Low Throughput%' THEN 1 ELSE 0 END) = 1
             AND MAX(CASE WHEN e.throughput LIKE '%High Throughput%' THEN 1 ELSE 0 END) = 1
                THEN 15
            WHEN MAX(CASE WHEN e.throughput LIKE '%Low Throughput%' THEN 1 ELSE 0 END) = 1
                THEN 10
            ELSE 5
        END AS throughput_pts,
        MAX(CASE WHEN e.rel_pts > 0 THEN e.rel_pts + (3 * e.has_oma_group) ELSE 0 END) AS orthology_pts
    FROM evidence e
    JOIN gene g ON g.entrez_id = e.partner_entrez_id
    JOIN organism org ON org.organism_id = e.evidence_org_id
    GROUP BY e.partner_entrez_id, g.biogrid_id, g.official_symbol, g.essentiality
)
SELECT
    partner_entrez_id,
    biogrid_id,
    official_symbol,
    essentiality,
    (directness_pts + study_pts + multi_species_pts + throughput_pts + orthology_pts) AS total_score,
    evidence_paths,
    evidence_organisms,
    study_count,
    throughput_summary
FROM scored
ORDER BY
    essentiality_rank ASC,
    total_score DESC,
    study_count DESC,
    partner_entrez_id ASC
LIMIT :limit;
"""

EXPLAIN_SQL = """
SELECT
    'Direct Human' AS path_type,
    'Human' AS evidence_org,
    COALESCE(g1.official_symbol, CAST(g1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(g2.official_symbol, CAST(g2.entrez_id AS TEXT)) AS interactor_2,
    '-' AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    COALESCE(p.author, '-') AS author,
    sl.throughput,
    COALESCE(sl.score, '-') AS raw_score,
    COALESCE(sl.qualifications, '-') AS qualifications
FROM v_sl_undirected sl
JOIN gene g1 ON g1.entrez_id = sl.gene_1_entrez_id
JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE sl.gene_1_entrez_id = :input_entrez_id
  AND sl.gene_2_entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Cross-Species' AS path_type,
    org.common_name AS evidence_org,
    COALESCE(g1.official_symbol, CAST(g1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(gm.official_symbol, CAST(gm.entrez_id AS TEXT)) AS interactor_2,
    o2.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    COALESCE(p.author, '-') AS author,
    sl.throughput,
    COALESCE(sl.score, '-') AS raw_score,
    COALESCE(sl.qualifications, '-') AS qualifications
FROM v_sl_undirected sl
JOIN gene g1 ON g1.entrez_id = sl.gene_1_entrez_id
JOIN gene gm ON gm.entrez_id = sl.gene_2_entrez_id
JOIN organism org ON org.organism_id = gm.organism_id
JOIN v_gene_orthology o2 ON o2.model_entrez_id = gm.entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE sl.gene_1_entrez_id = :input_entrez_id
  AND gm.organism_id <> 9606
  AND o2.human_entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Ortholog' AS path_type,
    org.common_name AS evidence_org,
    COALESCE(gm1.official_symbol, CAST(gm1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(g2.official_symbol, CAST(g2.entrez_id AS TEXT)) AS interactor_2,
    o1.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    COALESCE(p.author, '-') AS author,
    sl.throughput,
    COALESCE(sl.score, '-') AS raw_score,
    COALESCE(sl.qualifications, '-') AS qualifications
FROM v_gene_orthology o1
JOIN gene gm1 ON gm1.entrez_id = o1.model_entrez_id
JOIN organism org ON org.organism_id = o1.model_organism_id
JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE o1.human_entrez_id = :input_entrez_id
  AND g2.entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Ortholog' AS path_type,
    org.common_name AS evidence_org,
    COALESCE(gm1.official_symbol, CAST(gm1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(gm2.official_symbol, CAST(gm2.entrez_id AS TEXT)) AS interactor_2,
    o1.best_rel_type || ' -> ' || o2.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    COALESCE(p.author, '-') AS author,
    sl.throughput,
    COALESCE(sl.score, '-') AS raw_score,
    COALESCE(sl.qualifications, '-') AS qualifications
FROM v_gene_orthology o1
JOIN gene gm1 ON gm1.entrez_id = o1.model_entrez_id
JOIN organism org ON org.organism_id = o1.model_organism_id
JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
JOIN gene gm2 ON gm2.entrez_id = sl.gene_2_entrez_id
JOIN v_gene_orthology o2 ON o2.model_entrez_id = gm2.entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE o1.human_entrez_id = :input_entrez_id
  AND o2.human_entrez_id = :partner_entrez_id;
"""


def resolve_gene(conn: sqlite3.Connection, id_type: str, id_value: str):
    base_select = """
        SELECT g.entrez_id, g.biogrid_id, COALESCE(g.official_symbol, '-') AS symbol,
               g.essentiality, o.organism_id, o.common_name, o.scientific_name
        FROM gene g
        JOIN organism o ON o.organism_id = g.organism_id
    """
    if id_type == "entrez":
        try:
            val = int(id_value)
        except ValueError:
            sys.exit(f"Error: Entrez ID must be an integer, got '{id_value}'.")
        return conn.execute(base_select + " WHERE g.entrez_id = ?", (val,)).fetchone()
    if id_type == "biogrid":
        try:
            val = int(id_value)
        except ValueError:
            sys.exit(f"Error: BioGRID ID must be an integer, got '{id_value}'.")
        return conn.execute(base_select + " WHERE g.biogrid_id = ?", (val,)).fetchone()
    if id_type == "oma":
        return conn.execute(
            base_select + " JOIN gene_oma go ON go.entrez_id = g.entrez_id WHERE UPPER(go.oma_id) = UPPER(?)",
            (id_value.strip(),),
        ).fetchone()
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Search synthetic lethal partner genes for a mutated human gene."
    )
    parser.add_argument("--type", required=True, choices=["entrez", "biogrid", "oma"], help="Input ID type")
    parser.add_argument("--id", required=True, help="Gene identifier value")
    parser.add_argument("--limit", type=int, default=20, help="Max results to display (default: 20)")
    parser.add_argument("--explain", type=int, default=None, help="Show detailed evidence for a partner Entrez ID")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="Path to SQLite database (default: synlet.db)")
    args = parser.parse_args()

    db_path = Path(args.db)
    if not db_path.exists():
        sys.exit(f"Error: Database not found at {db_path}. Run build_db.py first.")

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    gene = resolve_gene(conn, args.type, args.id)
    if gene is None:
        sys.exit(f"No gene found for {args.type} ID '{args.id}'.")

    print(
        f"Input Gene: {gene['symbol']} | Entrez: {gene['entrez_id']} | "
        f"BioGRID: {gene['biogrid_id'] or '-'} | Organism: {gene['common_name']} ({gene['organism_id']}) | "
        f"Essentiality: {gene['essentiality']}"
    )

    if gene["organism_id"] != 9606:
        print(
            f"\nWarning: Input gene belongs to {gene['common_name']} ({gene['scientific_name']}), "
            "not Human (9606). The 3-step lookup expects a Human gene."
        )
        sys.exit(1)

    if args.explain is not None:
        partner = resolve_gene(conn, "entrez", str(args.explain))
        p_label = f"{partner['symbol']} ({args.explain})" if partner else str(args.explain)
        rows = conn.execute(
            EXPLAIN_SQL, {"input_entrez_id": gene["entrez_id"], "partner_entrez_id": args.explain}
        ).fetchall()
        print(f"\nDetailed Evidence for Partner {p_label} ({len(rows)} interaction record(s)):")
        if not rows:
            print("  No evidence rows found between these genes.")
        for r in rows:
            print(
                f"  - [{r['path_type']} | {r['evidence_org']}] {r['interactor_1']} <-> {r['interactor_2']} "
                f"(OMA: {r['orthology_chain']}) | BioGRID ID: {r['biogrid_interaction_id']} | "
                f"{r['publication_id']} ({r['author']}) | {r['throughput']} | Score: {r['raw_score']}"
            )
            if r["qualifications"] != "-":
                print(f"    Qualifications: {r['qualifications']}")
            print(f"    Link: https://thebiogrid.org/interaction/{r['biogrid_interaction_id']}")
        return

    results = conn.execute(
        LOOKUP_SQL, {"input_entrez_id": gene["entrez_id"], "limit": args.limit}
    ).fetchall()

    print(f"\nFound {len(results)} ranked partner gene(s) (showing up to {args.limit}):\n")
    if not results:
        print("No synthetic lethal partners found in the database.")
        return

    header = (
        f"{'#':<4} {'Entrez':<10} {'BioGRID':<10} {'Symbol':<12} {'Essentiality':<15} "
        f"{'Score':<6} {'Studies':<8} {'TP':<6} {'Paths':<28} {'BioGRID Link'}"
    )
    print(header)
    print("-" * len(header) + "-" * 22)
    for idx, r in enumerate(results, start=1):
        bg_id = str(r["biogrid_id"]) if r["biogrid_id"] else "-"
        bg_link = f"https://thebiogrid.org/{bg_id}/" if r["biogrid_id"] else "-"
        print(
            f"{idx:<4} {r['partner_entrez_id']:<10} {bg_id:<10} {r['official_symbol']:<12} "
            f"{r['essentiality']:<15} {r['total_score']:<6} {r['study_count']:<8} "
            f"{r['throughput_summary']:<6} {r['evidence_paths']:<28} {bg_link}"
        )


if __name__ == "__main__":
    main()
