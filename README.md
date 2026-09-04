# Supplementary Code 1

R 4.3.3 code and compact data for the polystyrene-oligomer/PTC study.

## Run

Install the packages in `environment/requirements-r.txt` into an R 4.3.3 environment. No script installs software. From the package directory:

```bash
Rscript run_all.R --include-go
```

Alternatively, in R 4.3.3/RStudio set the working directory to this package and run `source("run_all.R")`. This runs the compact workflows without GO; use the command above to include GO. No Python installation is required.

Outputs are created under `reproduced_outputs/complete_run/`. Existing runs are not overwritten; select another directory with `--output reproduced_outputs/run_02`. Generated figures and logs are not distributed with the code.

## Workflow

| Entry | Analysis |
|---|---|
| 01 | Target intersection and source/GeneCards-threshold sensitivity |
| 02–05 | CTD filtering, STRING PPI and GO/KEGG enrichment |
| 06 / 06b | GEO2R DEG extraction / intersection with 75 shared genes |
| 07 / 08 | LASSO and SVM-RFE feature selection |
| 08b | Repeated nested resampling of the seven-gene panel (S13) |
| 09 / 10 | Discovery and external individual-gene ROC |
| 11 | CIBERSORT with external expression/signature inputs |
| 12 | Immune group statistics (S16) and separate Figure 8A–C panels |
| 13 | Gene–immune correlations (S17) and separate Figure 9A–B panels |
| 14–16 | Single-cell analyses with external matrices |
| 17 / 17b / 18–20 | IHC plotting functions / paired statistics / marker-specific figures |

The compact run includes 01, 06b, cached PPI plotting (03b), 07–10, 12–13 and 17b. GO (04) is optional. Reference data are under `derived_outputs/` and are never overwritten by a run. The target chain is 144 targets, 2,139 disease genes, 75 shared genes, seven candidates and three core genes (BCL2/BAX/FN1).

### Individual analyses

```bash
Rscript scripts/01_target_screening_and_sensitivity.R
Rscript scripts/06b_candidate_intersection.R
Rscript scripts/08b_nested_resampling.R
Rscript scripts/12_immune_group_comparison.R
Rscript scripts/13_gene_immune_correlations.R
Rscript scripts/17b_IHC_paired_statistics.R
```

Scripts 12 and 13 compute statistics from fractions/expression and then plot the results. They do not substitute saved P values for calculation. Their arguments are output directory, figure directory and retained-fraction CSV; script 13 accepts the core-expression CSV as a fourth argument. All have package-relative defaults.

## Statistical settings

S13 uses the fixed seven-gene panel and the sample-labelled plan in `data/processed/ml_resampling_plan.json`: five outer folds repeated ten times, with five inner folds per outer training set. Scaling uses training-set means and population standard deviations. L1 logistic regression minimizes mean log loss plus the L1 penalty, including a unit bias term, with C from 10^-3 to 10^3 (13 values). RFE ranks squared coefficients from an L2-penalized squared-hinge linear model (C=1); feature count is selected by inner AUC. The final linear SVM uses C=2^-4 to 2^4. Objectives and convergence checks are explicit in script 08b. This is a conditional fixed-panel assessment; the 50 outer folds are dependent resamples.

ROC CIs use 1,000 bootstrap resamples without a fixed statistical seed. AUCs, expression means and Welch P values are deterministic for the same inputs; CI endpoints can fluctuate. Verification records these differences. Plot-only jitter seeds do not fix the statistical bootstrap.

Immune group comparisons use BH over 22 tests. Gene–fraction correlations use 66 estimable pooled tests and 63 estimable PTC-only tests out of 66 planned pairs. Three constant-input pairs remain undefined. Script 13 exports test counts and estimability. The reference-compatible column `PTC_only_BH_FDR_66_tests` identifies the planned grid; its q values are calculated over 63 finite P values.

Figure 8C is a descriptive Pearson matrix without significance stars. Figure 9A uses a descriptive Spearman cell–cell matrix; gene–fraction links and Figure 9B stars use pooled BH q values.

## External data

Large raw datasets are not included. Sources and formats are described in `data/README.md` and `data/source_files/README.md`. To extract DEGs from an external GEO2R export:

```bash
Rscript run_all.R --output reproduced_outputs/with_geo2r --geo2r-table external_data/GSE33630_GEO2R_top_table.tsv
```

Script 06 filters the GEO2R table; it does not refit raw-microarray models. Large single-cell and upstream CIBERSORT analyses are separate from the compact run. The CIBERSORT signature must be supplied under its applicable terms. Live database/annotation updates can change network and enrichment outputs.

## Checks

```bash
Rscript tests/check_package.R
Rscript tests/verify_reproduced_results.R reproduced_outputs/complete_run
Rscript tests/test_ihc_statistics.R
```

The run checks numerical outputs automatically. `--skip-s13` can omit resampling for a shorter run; the verification report explicitly states that S13 was not checked.

Clinical files contain deidentified pair labels only. Do not distribute linkage keys, raw downloads or execution logs. The MIT license applies to author-owned code; third-party materials retain their source terms. Insert the repository address into `CITATION.cff` before publication.
