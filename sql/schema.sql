PRAGMA foreign_keys = ON;

-- 1. Organisms (Human, Yeast, Fruit fly)
CREATE TABLE organism (
    organism_id     INTEGER PRIMARY KEY,
    scientific_name TEXT NOT NULL,
    common_name     TEXT NOT NULL
);

-- 2. Canonical Genes (keyed by Entrez Gene ID)
CREATE TABLE gene (
    entrez_id                 INTEGER PRIMARY KEY,
    organism_id               INTEGER NOT NULL REFERENCES organism(organism_id),
    biogrid_id                INTEGER UNIQUE,
    official_symbol           TEXT,
    systematic_name           TEXT,
    essentiality              TEXT NOT NULL DEFAULT 'Unknown'
                                  CHECK (essentiality IN ('Non-essential', 'Unknown', 'Essential')),
    essentiality_source_count INTEGER NOT NULL DEFAULT 0
);

-- 3. OMA Protein/Isoform IDs (1:N from gene)
CREATE TABLE gene_oma (
    oma_id    TEXT PRIMARY KEY,
    entrez_id INTEGER NOT NULL REFERENCES gene(entrez_id)
);

-- 4. Gene Aliases & External Accessions (preserved from BioGRID)
CREATE TABLE gene_alias (
    entrez_id   INTEGER NOT NULL REFERENCES gene(entrez_id),
    alias_type  TEXT NOT NULL CHECK (alias_type IN ('SYNONYM', 'SWISS-PROT', 'TREMBL', 'REFSEQ')),
    alias_value TEXT NOT NULL,
    PRIMARY KEY (entrez_id, alias_type, alias_value)
);

-- 5. Publications (BioGRID studies + G-HumanEssential PMIDs)
CREATE TABLE publication (
    publication_id TEXT PRIMARY KEY,
    author         TEXT
);

-- 6. Supporting PMIDs for Human Gene Essentiality (M:N)
CREATE TABLE gene_essentiality_source (
    entrez_id      INTEGER NOT NULL REFERENCES gene(entrez_id),
    publication_id TEXT NOT NULL REFERENCES publication(publication_id),
    PRIMARY KEY (entrez_id, publication_id)
);

-- 7. Orthology Pairs (Human OMA ID <-> Yeast/DROME OMA ID)
CREATE TABLE orthology (
    oma_id_1  TEXT NOT NULL REFERENCES gene_oma(oma_id),
    oma_id_2  TEXT NOT NULL REFERENCES gene_oma(oma_id),
    rel_type  TEXT NOT NULL CHECK (rel_type IN ('1:1', '1:n', 'm:1', 'm:n')),
    oma_group INTEGER,
    PRIMARY KEY (oma_id_1, oma_id_2)
);

-- 8. Synthetic Lethality Interactions (BioGRID experimental records)
CREATE TABLE synthetic_lethality (
    biogrid_interaction_id INTEGER PRIMARY KEY,
    gene_a_entrez_id       INTEGER NOT NULL REFERENCES gene(entrez_id),
    gene_b_entrez_id       INTEGER NOT NULL REFERENCES gene(entrez_id),
    publication_id         TEXT NOT NULL REFERENCES publication(publication_id),
    throughput             TEXT NOT NULL,
    score                  TEXT,
    qualifications         TEXT,
    ontology_term_names    TEXT,
    CHECK (gene_a_entrez_id <> gene_b_entrez_id)
);

-- Indexes for fast ID resolution and graph traversal
CREATE INDEX idx_gene_organism   ON gene(organism_id);
CREATE INDEX idx_gene_symbol     ON gene(official_symbol);
CREATE INDEX idx_gene_oma_entrez ON gene_oma(entrez_id);
CREATE INDEX idx_orthology_oma2  ON orthology(oma_id_2);
CREATE INDEX idx_sl_gene_a       ON synthetic_lethality(gene_a_entrez_id);
CREATE INDEX idx_sl_gene_b       ON synthetic_lethality(gene_b_entrez_id);
CREATE INDEX idx_sl_pub          ON synthetic_lethality(publication_id);

-- Deduplicated Gene-to-Gene Orthology View (prevents isoform row multiplication)
CREATE VIEW v_gene_orthology AS
WITH ranked AS (
    SELECT
        go1.entrez_id AS human_entrez_id,
        go2.entrez_id AS model_entrez_id,
        g2.organism_id AS model_organism_id,
        o.rel_type,
        CASE o.rel_type
            WHEN '1:1' THEN 1
            WHEN '1:n' THEN 2
            WHEN 'm:1' THEN 2
            ELSE 3
        END AS rel_rank,
        CASE WHEN o.oma_group IS NOT NULL THEN 1 ELSE 0 END AS has_oma_group
    FROM orthology o
    JOIN gene_oma go1 ON go1.oma_id = o.oma_id_1
    JOIN gene_oma go2 ON go2.oma_id = o.oma_id_2
    JOIN gene g2 ON g2.entrez_id = go2.entrez_id
)
SELECT
    human_entrez_id,
    model_entrez_id,
    model_organism_id,
    CASE MIN(rel_rank)
        WHEN 1 THEN '1:1'
        WHEN 2 THEN MIN(CASE WHEN rel_rank = 2 THEN rel_type END)
        ELSE 'm:n'
    END AS best_rel_type,
    CASE MIN(rel_rank)
        WHEN 1 THEN 10
        WHEN 2 THEN 5
        ELSE 2
    END AS rel_pts,
    MAX(has_oma_group) AS has_oma_group
FROM ranked
GROUP BY human_entrez_id, model_entrez_id, model_organism_id;

-- Bidirectional Synthetic Lethality View (makes all SL lookups symmetric)
CREATE VIEW v_sl_undirected AS
SELECT
    biogrid_interaction_id,
    gene_a_entrez_id AS gene_1_entrez_id,
    gene_b_entrez_id AS gene_2_entrez_id,
    publication_id,
    throughput,
    score,
    qualifications,
    ontology_term_names
FROM synthetic_lethality
UNION ALL
SELECT
    biogrid_interaction_id,
    gene_b_entrez_id AS gene_1_entrez_id,
    gene_a_entrez_id AS gene_2_entrez_id,
    publication_id,
    throughput,
    score,
    qualifications,
    ontology_term_names
FROM synthetic_lethality;
