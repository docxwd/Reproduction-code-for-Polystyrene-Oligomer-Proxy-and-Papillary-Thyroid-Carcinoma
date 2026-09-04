# Data inputs

Only compact analysis inputs and reference outputs are included. Large raw downloads and complete expression matrices are not bundled.

| Input | Role and source |
|---|---|
| processed/target_source_membership_144.csv | CTD/TargetNet/SwissTargetPrediction processed source membership |
| processed/thyroid_cancer_genes_2139.csv | GeneCards/OMIM processed disease-gene union |
| processed/genecards_top3000.csv | Ranked GeneCards export subset for sensitivity checks |
| processed/shared_targets_75.csv | Shared-target identifiers and provenance |
| processed/discovery_7gene_expression.csv | The 7 × 94 GSE33630 matrix used for feature selection; copied without numeric modification |
| processed/ml_resampling_plan.json | Sample-labelled outer and inner partitions for the 50-fold fixed-panel resampling |
| processed/external_core_expression.csv | GSE60542 expression input used for individual-gene validation; copied without numeric modification |
| processed/immune_core_expression.csv | Three core genes across 94 GSE33630 samples, extracted from the saved expression matrix after trim/uppercase symbol grouping and arithmetic-mean aggregation |
| processed/core_genes_no_header.csv | Three-gene list for ROC scripts |
| processed/candidate_7genes.txt | Candidate list for single-cell input preparation |
| processed/final_core_genes_BCL2_BAX_FN1.csv | Core list for single-cell input preparation |
| processed/ihc_hscore_deidentified.csv | Recorded paired H-score measurements with deidentified Pair_ID labels |
| processed/BAX_Hscore_deidentified.csv and corresponding BCL2/FN1 files | Marker-specific inputs to IHC plotting entry points |
| Other processed HPA/clinical/cohort tables | Source identifiers and deidentified clinical context accompanying the manuscript |
| derived_outputs/cibersort_retained_samples.csv | Saved 42 × 22 CIBERSORT fractions with sample groups; input to downstream statistics, not a new deconvolution |
| Other derived_outputs tables | Reference numerical outputs for comparison, including the GEO2R-derived gene summary |

Expression sample names carry the public GEO accession and `_con`/`_tre` (or `_treat`) group suffix. The compact matrices are processed analytical inputs, not raw microarray intensities.

GEO records: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60542

Reference numerical tables are separate from newly generated outputs. Script 13 reports which correlations are estimable. Missing/undefined correlations are not filled with zeros or P=1.

No identifying clinical originals or Pair_ID linkage key are included. Large raw-data backups remain separate from this code package. `manifest.tsv` describes input/output roles; `docs/checksums_sha256.txt` covers release file integrity.
