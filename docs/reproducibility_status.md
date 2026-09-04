# Reproducibility

Use R 4.3.3 and the listed dependencies. Fixed-panel resampling uses the supplied sample-labelled fold plan. The compact workflow recomputes statistics from processed inputs, and its R checks compare results with reference tables.

S13 reference values correspond to script 08b. Mean outer AUCs are 0.9699506173 (L1 logistic regression) and 0.9661234568 (linear SVM-RFE). BCL2/BAX/FN1 are selected 50/49/50 times by L1 and 42/44/45 times by SVM-RFE.

ROC bootstrap CI endpoints are stochastic; verification requires agreement of AUC, expression means and P values and records CI differences. Undefined constant-input immune correlations remain missing.

Cached STRING, GO and KEGG tables are included as snapshots. The compact run does not represent an end-to-end raw-microarray, upstream CIBERSORT, live KEGG or large single-cell rerun.
