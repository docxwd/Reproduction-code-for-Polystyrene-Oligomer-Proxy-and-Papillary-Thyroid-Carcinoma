## =========================================================
## Clean SVM-RFE + Export in your original format
## Input: geneexp.csv (rows=genes, cols=samples, colnames end with _con/_tre)
## Output (same style as original script):
##   feature_svm.txt, errors.jpeg, accuracy.jpeg, combined_plots.jpeg, SVM-RFE.gene.txt
## =========================================================

## ---- packages ----
pkgs <- c("e1071", "caret", "pROC")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0L) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       ". See environment/requirements-r.txt.")
}
library(e1071)
library(caret)
library(pROC)

set.seed(20260127)
cli_args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
if (!dir.exists(work_dir)) stop("Analysis directory not found: ", work_dir)
setwd(work_dir)

## ---------------------- Plot functions (publication style) ----------------------
draw_label_box <- function(x, y, label, cex = 0.72) {
  label_lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  box_width <- max(strwidth(label_lines, cex = cex)) + 0.18
  box_height <- length(label_lines) * strheight("0", cex = cex) * 1.55
  rect(
    x - box_width / 2, y - box_height / 2,
    x + box_width / 2, y + box_height / 2,
    col = grDevices::adjustcolor("white", alpha.f = 0.94),
    border = "#D6DCE3",
    lwd = 0.8
  )
  text(x, y, labels = label, cex = cex, col = "#263238", font = 2)
}

plot_metric_curve <- function(values, optimum = c("min", "max"),
                              xlab = "Number of features",
                              ylab = "CV error",
                              main = "Cross-validation error",
                              line_col = "#3B73B9",
                              best_col = "#E15759",
                              best_name = "Best k",
                              metric_name = "Error") {
  optimum <- match.arg(optimum)
  x <- which(!is.na(values))
  y <- values[x]
  y_range <- range(y, na.rm = TRUE)
  y_pad <- diff(y_range) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) {
    y_pad <- max(abs(y_range), 1) * 0.04
  }
  ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

  oldPar <- par(
    mar = c(4.4, 4.7, 3.1, 1.0),
    mgp = c(2.7, 0.75, 0),
    tck = -0.018,
    xaxs = "i",
    yaxs = "i"
  )
  on.exit(par(oldPar))

  plot(
    x, y,
    type = "n",
    xlim = c(min(x) - 0.35, max(x) + 0.35),
    ylim = ylim,
    xaxt = "n",
    yaxt = "n",
    xlab = xlab,
    ylab = ylab,
    main = main,
    bty = "n",
    cex.lab = 1.02,
    cex.axis = 0.88,
    cex.main = 1.08,
    font.lab = 2,
    font.main = 2,
    col.lab = "#222222",
    col.main = "#222222"
  )

  abline(h = pretty(ylim, n = 5), col = "#E8ECF1", lwd = 0.9)
  abline(v = x, col = "#F3F5F7", lwd = 0.75)
  lines(x, y, col = line_col, lwd = 2.4)
  points(x, y, pch = 21, bg = "white", col = line_col, lwd = 1.35, cex = 1.05)

  best_x <- if (optimum == "min") which.min(values) else which.max(values)
  best_y <- values[best_x]
  segments(best_x, ylim[1], best_x, best_y, col = best_col, lty = 3, lwd = 1.1)
  points(best_x, best_y, pch = 21, bg = best_col, col = "#222222", lwd = 1.1, cex = 1.45)

  label_x <- best_x + if (best_x <= stats::median(x)) 0.9 else -0.9
  label_x <- min(max(label_x, min(x) + 0.85), max(x) - 0.85)
  label_y <- best_y + if (optimum == "min") diff(ylim) * 0.19 else -diff(ylim) * 0.19
  label_y <- min(max(label_y, ylim[1] + diff(ylim) * 0.18), ylim[2] - diff(ylim) * 0.18)
  draw_label_box(
    label_x,
    label_y,
    sprintf("%s = %d\n%s = %.3f", best_name, best_x, metric_name, best_y)
  )

  y_ticks <- pretty(ylim, n = 5)
  y_ticks <- y_ticks[y_ticks >= ylim[1] & y_ticks <= ylim[2]]
  axis(1, at = x, labels = x, lwd = 0.9, col = "#222222", col.axis = "#222222")
  axis(
    2,
    at = y_ticks,
    labels = formatC(y_ticks, format = "f", digits = 3),
    las = 1,
    lwd = 0.9,
    col = "#222222",
    col.axis = "#222222"
  )
  box(bty = "l", lwd = 0.9, col = "#222222")
}

PlotErrors <- function(errors, errors2 = NULL, no.info = 0.5,
                       ylim = range(c(errors, errors2), na.rm = TRUE),
                       xlab = "Number of features", ylab = "CV error") {
  plot_metric_curve(
    errors,
    optimum = "min",
    xlab = xlab,
    ylab = ylab,
    main = "SVM-RFE Error",
    line_col = "#3B73B9",
    best_col = "#E15759",
    metric_name = "Error"
  )
}

PlotAccuracy <- function(acc, acc2 = NULL, no.info = 0.5,
                         ylim = range(c(acc, acc2), na.rm = TRUE),
                         xlab = "Number of features", ylab = "CV accuracy") {
  plot_metric_curve(
    acc,
    optimum = "max",
    xlab = xlab,
    ylab = ylab,
    main = "SVM-RFE Accuracy",
    line_col = "#59A14F",
    best_col = "#F28E2B",
    metric_name = "Accuracy"
  )
}

## ---------------------- 1) Read & reshape ----------------------
inputFile <- if (length(cli_args) >= 2L) cli_args[[2L]] else "geneexp.csv"
raw <- read.csv(inputFile, header = TRUE, check.names = FALSE)

stopifnot("ID" %in% colnames(raw))
genes <- raw$ID
expr  <- raw[, -1, drop = FALSE]
colnames(expr) <- trimws(colnames(expr))

group <- ifelse(grepl("_con$", colnames(expr)), "con",
                ifelse(grepl("_tre$", colnames(expr)), "tre", NA))
if (any(is.na(group))) stop("存在列名不以 _con/_tre 结尾，无法分组：\n", paste(colnames(expr)[is.na(group)], collapse = ", "))
group <- factor(group, levels = c("con", "tre"))

data <- as.data.frame(t(as.matrix(expr)))
colnames(data) <- genes
data <- cbind(group = group, data)
data$group <- factor(data$group, levels = c("con", "tre"))

p <- ncol(data) - 1
cat("Samples:", nrow(data), " Features:", p, "\n")

## ---------------------- 2) SVM-RFE (clean) ----------------------
get_linear_svm_weights <- function(model) {
  w <- t(model$coefs) %*% model$SV
  as.numeric(w)
}

svm_rfe_rank_one <- function(train_df, cost = 1) {
  feats <- setdiff(colnames(train_df), "group")
  alive <- feats
  rank_pos <- rep(NA_integer_, length(feats)); names(rank_pos) <- feats
  
  while (length(alive) > 1) {
    m <- svm(x = train_df[, alive, drop = FALSE],
             y = train_df$group,
             kernel = "linear", scale = TRUE, cost = cost,
             type = "C-classification")
    w <- get_linear_svm_weights(m); names(w) <- alive
    drop_feat <- names(which.min(abs(w)))
    rank_pos[drop_feat] <- length(alive)
    alive <- setdiff(alive, drop_feat)
  }
  rank_pos[alive] <- 1L
  list(rank = rank_pos)
}

svm_rfe_cv_rank <- function(dat, folds = 5, cost = 1) {
  idx <- createFolds(dat$group, k = folds, returnTrain = TRUE)
  feats <- setdiff(colnames(dat), "group")
  rank_sum <- setNames(rep(0, length(feats)), feats)
  
  for (i in seq_along(idx)) {
    tr <- idx[[i]]
    r <- svm_rfe_rank_one(dat[tr, , drop = FALSE], cost = cost)
    rank_sum <- rank_sum + r$rank
  }
  rank_sum
}

folds_used <- 5
rank_sum <- svm_rfe_cv_rank(data, folds = folds_used, cost = 1)

final_order <- names(sort(rank_sum, decreasing = FALSE))
# 输出成你原脚本的列名/风格：FeatureName / FeatureID / AvgRank
topFeatures <- data.frame(
  FeatureName = final_order,
  FeatureID   = match(final_order, colnames(data)[-1]),  # 1..p
  AvgRank     = as.numeric(rank_sum[final_order]) / folds_used,
  stringsAsFactors = FALSE
)

write.table(topFeatures, file = "feature_svm.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

## ---------------------- 3) CV performance sweep (k = 1..p) ----------------------
eval_topk_cv <- function(dat, ordered_feats, folds = 5) {
  idx <- createFolds(dat$group, k = folds, returnTrain = TRUE)
  p <- length(ordered_feats)
  
  out <- data.frame(k = 1:p, Accuracy = NA_real_, AUC = NA_real_)
  
  for (k in 1:p) {
    feats_k <- ordered_feats[1:k]
    acc_vec <- c(); auc_vec <- c()
    
    for (i in seq_along(idx)) {
      tr <- idx[[i]]
      te <- setdiff(seq_len(nrow(dat)), tr)
      
      # 小网格调参（你特征只有1~5个，不需要巨网格）
      tune_obj <- tune(
        svm,
        train.x = dat[tr, feats_k, drop = FALSE],
        train.y = dat$group[tr],
        kernel = "linear", scale = TRUE,
        ranges = list(cost = 2^(-2:4)),
        tunecontrol = tune.control(sampling = "fix")
      )
      best_cost <- tune_obj$best.parameters$cost
      
      m <- svm(x = dat[tr, feats_k, drop = FALSE],
               y = dat$group[tr],
               kernel = "linear", scale = TRUE,
               cost = best_cost, probability = TRUE)
      
      pred <- predict(m, dat[te, feats_k, drop = FALSE], probability = TRUE)
      acc_vec <- c(acc_vec, mean(pred == dat$group[te]))
      
      prob <- attr(pred, "probabilities")[, "tre"]
      roc_obj <- pROC::roc(dat$group[te], prob, levels = c("con", "tre"), quiet = TRUE)
      auc_vec <- c(auc_vec, as.numeric(pROC::auc(roc_obj)))
    }
    
    out$Accuracy[k] <- mean(acc_vec)
    out$AUC[k] <- mean(auc_vec)
  }
  out
}

perf <- eval_topk_cv(data, topFeatures$FeatureName, folds = folds_used)
errors <- 1 - perf$Accuracy
no.info <- min(prop.table(table(data$group)))  # 和原脚本一致思路：无信息错误率 :contentReference[oaicite:4]{index=4}

## ---------------------- 4) Plot exports (same filenames as original) ----------------------
jpeg(filename = "errors.jpeg", width = 5, height = 5, units = "in", res = 600, quality = 95, bg = "white")
PlotErrors(errors, no.info = no.info)
dev.off()

jpeg(filename = "accuracy.jpeg", width = 5, height = 5, units = "in", res = 600, quality = 95, bg = "white")
PlotAccuracy(1 - errors, no.info = no.info)
dev.off()

jpeg(filename = "combined_plots.jpeg", width = 10, height = 5, units = "in", res = 600, quality = 95, bg = "white")
par(mfrow = c(1, 2))
PlotErrors(errors, no.info = no.info)
PlotAccuracy(1 - errors, no.info = no.info)
dev.off()

## ---------------------- 5) Pick optimal k & export gene list ----------------------
optimalFeatureCount <- which.min(errors)
featureGenes <- topFeatures[1:optimalFeatureCount, "FeatureName", drop = FALSE]

write.table(featureGenes, file = "SVM-RFE.gene.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Done! Saved: feature_svm.txt, errors.jpeg, accuracy.jpeg, combined_plots.jpeg, SVM-RFE.gene.txt\n")



