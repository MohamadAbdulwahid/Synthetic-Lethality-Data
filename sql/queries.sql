-- =============================================================================
-- SynLet Project - Reference SQL Queries
-- =============================================================================

-- -----------------------------------------------------------------------------
-- 1. ID RESOLUTION (Search by Entrez ID, BioGRID ID, or OMA ID + ID Type)
-- -----------------------------------------------------------------------------

-- 1a. Resolve by Entrez ID
SELECT g.entrez_id, g.biogrid_id, g.official_symbol, g.essentiality,
       o.organism_id, o.common_name, o.scientific_name
FROM gene g
JOIN organism o ON o.organism_id = g.organism_id
WHERE g.entrez_id = :input_entrez_id;

-- 1b. Resolve by BioGRID ID
SELECT g.entrez_id, g.biogrid_id, g.official_symbol, g.essentiality,
       o.organism_id, o.common_name, o.scientific_name
FROM gene g
JOIN organism o ON o.organism_id = g.organism_id
WHERE g.biogrid_id = :input_biogrid_id;

-- 1c. Resolve by OMA ID
SELECT g.entrez_id, g.biogrid_id, g.official_symbol, g.essentiality,
       o.organism_id, o.common_name, o.scientific_name, go.oma_id
FROM gene_oma go
JOIN gene g ON g.entrez_id = go.entrez_id
JOIN organism o ON o.organism_id = g.organism_id
WHERE UPPER(go.oma_id) = UPPER(:input_oma_id);


-- -----------------------------------------------------------------------------
-- 2. MAIN 3-STEP SYNTHETIC LETHALITY LOOKUP & RANKING QUERY
--    Parameter: :input_entrez_id (a Human gene Entrez ID, organism_id = 9606)
--
--    Paths combined:
--      1. Direct Human:   H1 --(SL)-- H2
--      2. Cross-Species:  H1 --(SL)-- Model M2 --(OMA)-- H2
--      3. Ortholog:       H1 --(OMA)-- Model M1 --(SL)-- (H2 or M2 --(OMA)-- H2)
--
--    Sort order:
--      Primary:   Essentiality tier (Non-essential -> Unknown -> Essential)
--      Secondary: Confidence score DESC
-- -----------------------------------------------------------------------------
WITH evidence AS (
    -- Path 1: Direct Human <-> Human
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

    -- Path 2: Cross-Species Human H1 <-> Model M2 -> Human H2
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

    -- Path 3a: Human H1 -> Model M1 <-> Human H2 (Cross-species via ortholog)
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

    -- Path 3b: Human H1 -> Model M1 <-> Model M2 -> Human H2
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
    directness_pts,
    study_pts,
    multi_species_pts,
    throughput_pts,
    orthology_pts,
    evidence_paths,
    evidence_organisms,
    study_count,
    interaction_count,
    throughput_summary
FROM scored
ORDER BY
    essentiality_rank ASC,
    total_score DESC,
    study_count DESC,
    partner_entrez_id ASC;


-- -----------------------------------------------------------------------------
-- 3. EVIDENCE DETAIL QUERY FOR A SPECIFIC PAIR (:input_entrez_id, :partner_entrez_id)
-- -----------------------------------------------------------------------------
SELECT
    'Direct Human' AS path_type,
    9606 AS evidence_org_id,
    g1.official_symbol AS interactor_1,
    g2.official_symbol AS interactor_2,
    NULL AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    p.author,
    sl.throughput,
    sl.score,
    sl.qualifications,
    'https://thebiogrid.org/interaction/' || sl.biogrid_interaction_id AS biogrid_url
FROM v_sl_undirected sl
JOIN gene g1 ON g1.entrez_id = sl.gene_1_entrez_id
JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE sl.gene_1_entrez_id = :input_entrez_id
  AND sl.gene_2_entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Cross-Species' AS path_type,
    gm.organism_id AS evidence_org_id,
    g1.official_symbol AS interactor_1,
    COALESCE(gm.official_symbol, CAST(gm.entrez_id AS TEXT)) AS interactor_2,
    o2.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    p.author,
    sl.throughput,
    sl.score,
    sl.qualifications,
    'https://thebiogrid.org/interaction/' || sl.biogrid_interaction_id AS biogrid_url
FROM v_sl_undirected sl
JOIN gene g1 ON g1.entrez_id = sl.gene_1_entrez_id
JOIN gene gm ON gm.entrez_id = sl.gene_2_entrez_id
JOIN v_gene_orthology o2 ON o2.model_entrez_id = gm.entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE sl.gene_1_entrez_id = :input_entrez_id
  AND gm.organism_id <> 9606
  AND o2.human_entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Ortholog' AS path_type,
    o1.model_organism_id AS evidence_org_id,
    COALESCE(gm1.official_symbol, CAST(gm1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(g2.official_symbol, CAST(g2.entrez_id AS TEXT)) AS interactor_2,
    o1.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    p.author,
    sl.throughput,
    sl.score,
    sl.qualifications,
    'https://thebiogrid.org/interaction/' || sl.biogrid_interaction_id AS biogrid_url
FROM v_gene_orthology o1
JOIN gene gm1 ON gm1.entrez_id = o1.model_entrez_id
JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
JOIN gene g2 ON g2.entrez_id = sl.gene_2_entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE o1.human_entrez_id = :input_entrez_id
  AND g2.entrez_id = :partner_entrez_id

UNION ALL

SELECT
    'Ortholog' AS path_type,
    o1.model_organism_id AS evidence_org_id,
    COALESCE(gm1.official_symbol, CAST(gm1.entrez_id AS TEXT)) AS interactor_1,
    COALESCE(gm2.official_symbol, CAST(gm2.entrez_id AS TEXT)) AS interactor_2,
    o1.best_rel_type || ' -> ' || o2.best_rel_type AS orthology_chain,
    sl.biogrid_interaction_id,
    sl.publication_id,
    p.author,
    sl.throughput,
    sl.score,
    sl.qualifications,
    'https://thebiogrid.org/interaction/' || sl.biogrid_interaction_id AS biogrid_url
FROM v_gene_orthology o1
JOIN gene gm1 ON gm1.entrez_id = o1.model_entrez_id
JOIN v_sl_undirected sl ON sl.gene_1_entrez_id = o1.model_entrez_id
JOIN gene gm2 ON gm2.entrez_id = sl.gene_2_entrez_id
JOIN v_gene_orthology o2 ON o2.model_entrez_id = gm2.entrez_id
JOIN publication p ON p.publication_id = sl.publication_id
WHERE o1.human_entrez_id = :input_entrez_id
  AND o2.human_entrez_id = :partner_entrez_id;


-- -----------------------------------------------------------------------------
-- 4. DATABASE SUMMARY STATS
-- -----------------------------------------------------------------------------
SELECT o.organism_id, o.common_name, o.scientific_name,
       COUNT(g.entrez_id) AS total_genes,
       SUM(CASE WHEN g.essentiality = 'Non-essential' THEN 1 ELSE 0 END) AS non_essential,
       SUM(CASE WHEN g.essentiality = 'Unknown' THEN 1 ELSE 0 END) AS unknown_ess,
       SUM(CASE WHEN g.essentiality = 'Essential' THEN 1 ELSE 0 END) AS essential
FROM organism o
LEFT JOIN gene g ON g.organism_id = o.organism_id
GROUP BY o.organism_id, o.common_name, o.scientific_name
ORDER BY o.organism_id;
