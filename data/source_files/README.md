# External inputs

Place external files outside the tracked release, for example under `external_data/`. These files are not bundled. A GEO accession identifies the study; the exact input required by each script is listed below.

## Bulk expression and GEO2R

- Script 06: tab-separated GEO2R export for GSE33630, with ID, adj.P.Val, P.Value, t, B, logFC, Gene.symbol and Gene.title. Supply it with `run_all.R --geo2r-table`. The discovery dataset has 45 normal and 49 PTC samples. Script 06 extracts DEGs at adj.P.Val<0.05 and |logFC|>1.
- Scripts 07–10: compact matrices and the core-gene list are included; `run_all.R` prepares separate output work directories.
- Script 11: full gene-by-sample expression CSV (first column gene symbol; remaining columns public sample IDs with group suffixes), plus an appropriately licensed reference signature as a tab-separated gene-by-cell-type matrix. Use the expression scale appropriate to the original deconvolution workflow. The script accepts analysis directory, expression filename, and signature filename as its first three arguments.
- Scripts 12–13: the included retained-fraction matrix and core-expression matrix are sufficient for the group/correlation statistics and Figure 8/9 panels. A full expression matrix is not required for these downstream analyses.

Bulk study pages:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60542

## Single-cell inputs for scripts 14–16

Provide one external project root containing the following subdirectories expected by the scripts:

```text
28.单细胞分析_GSE191288_2026-06-12/
  raw/h5_files/       GEO-supplied per-sample H5 matrices
  inputs/
    candidate_7genes.txt
    final_core_genes_BCL2_BAX_FN1.csv
29.单细胞分析_GSE184362_2026-06-12/
  raw/mtx_files/      matching *_matrix.mtx.gz, *_features.tsv.gz, *_barcodes.tsv.gz
  inputs/
    candidate_7genes.txt
    final_core_genes_BCL2_BAX_FN1.csv
```

Copy the two small gene lists from `data/processed/` into each `inputs/` directory. Preserve the GEO sample prefixes because the scripts parse group/patient metadata from filenames. Run scripts 15 and 14 with that external project root, then script 16 with the same root. Script 16 reads the gene-by-cell-type summaries produced by the two analyses.

GSE184362 uses the full retained cells for expression summaries and at most 4,000 cells per sample for its UMAP/clustering subset. GSE191288 uses all retained cells. These scripts are separate from the compact default reproduction run.

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191288

## Database-source workflows

Script 02 requires the CTD Polystyrenes interaction export (D011137) to apply the distinct-PubMed support rule. TargetNet/SwissTargetPrediction and GeneCards/OMIM raw export preparation precedes the packaged processed membership tables. Preserve the recorded query and thresholds in `config/parameters.yml`.

Script 03 queries STRING using the shared-target input and script 05 queries KEGG. Network/enrichment reference snapshots are included; live service outputs may change. Cached PPI plotting needs no external data and runs through script 03b.

Third-party reference signatures, database exports and software components retain their own source terms; verify permission before redistributing external files.
