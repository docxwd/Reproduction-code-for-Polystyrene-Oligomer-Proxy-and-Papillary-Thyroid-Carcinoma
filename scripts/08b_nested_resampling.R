# Run with R 4.3.3; paths are relative to this package.
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")))
}
if (as.character(getRversion()) != "4.3.3") stop("Use R 4.3.3")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else tail(source_files[nzchar(source_files)], 1)
if (!length(this_file)) stop("Run this file with Rscript or source().")
root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts/lib/workflow_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)

# S13: repeated nested resampling of the fixed seven-gene panel.
# Folds are supplied with sample identifiers; no Python runtime is needed.
# L1-logistic: mean log-loss + sum(abs(w))/(C*n), with a penalized unit bias.
# RFE ranking: sum(w^2)/2 + sum(squared hinge loss), with a penalized unit bias.
# Final classifier: linear C-SVM (e1071); selection/tuning are within outer training sets.
require_packages(c("jsonlite", "e1071"))
suppressPackageStartupMessages(library(e1071))
out <- safe_output(if (length(args)) args[1] else file.path(root, "reproduced_outputs/nested_ml"), root)
d <- read_data(file.path(root, "data/processed/discovery_7gene_expression.csv"))
genes <- d$ID
X <- t(as.matrix(d[-1])); colnames(X) <- genes
y <- as.integer(endsWith(rownames(X), "_tre"))
plan <- jsonlite::fromJSON(file.path(root, "data/processed/ml_resampling_plan.json"), simplifyVector = FALSE)
stopifnot(identical(rownames(X),unlist(plan$sample_ids)),length(genes)==7L,
          !anyDuplicated(genes),all(is.finite(X)),sum(y==0)==45L,sum(y==1)==49L)
folds <- plan$folds
stopifnot(length(folds)==50L)
for (s in folds) {
  tr <- unlist(s$train); te <- unlist(s$test)
  stopifnot(!anyDuplicated(c(tr,te)),setequal(c(tr,te),seq_len(nrow(X))),length(s$inner)==5L)
  for (f in s$inner) {
    inner_tr <- unlist(f$train); inner_te <- unlist(f$test)
    stopifnot(!anyDuplicated(c(inner_tr,inner_te)),
              setequal(c(inner_tr,inner_te),seq_along(tr)))
  }
}
solver_checks <- list()
auc <- function(y,s) (sum(rank(s)[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))
scaled <- function(a,b) {
  center <- colMeans(a); sd_pop <- sqrt(colMeans(sweep(a,2,center)^2))
  sd_pop[sd_pop==0] <- 1
  list(train=sweep(sweep(a,2,center),2,sd_pop,"/"),test=sweep(sweep(b,2,center),2,sd_pop,"/"))
}
l1_fit <- function(x,y,C) {
  z <- cbind(x,1); p <- ncol(z); lambda <- 1/(C*nrow(z))
  objective <- function(a) {
    eta <- drop(z %*% (a[seq_len(p)]-a[p+seq_len(p)]))
    mean(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta)+lambda*sum(a)
  }
  gradient <- function(a) {
    w <- a[seq_len(p)]-a[p+seq_len(p)]
    g <- drop(crossprod(z,plogis(z %*% w)-y))/nrow(z)
    c(g+lambda,-g+lambda)
  }
  fit <- optim(rep(0,2*p),objective,gradient,method="L-BFGS-B",lower=rep(0,2*p),
               control=list(maxit=10000,factr=1e4,pgtol=1e-8))
  w <- fit$par[seq_len(p)]-fit$par[p+seq_len(p)]
  g <- drop(crossprod(z,plogis(z%*%w)-y))/nrow(z)
  kkt <- max(ifelse(abs(w)>1e-10,abs(g+lambda*sign(w)),pmax(abs(g)-lambda,0)))
  if (fit$convergence != 0L || !is.finite(kkt) || kkt > 1e-4) stop("L1 optimizer convergence check failed")
  solver_checks[[length(solver_checks)+1L]] <<- c(convergence=fit$convergence,kkt=kkt)
  w
}
linear_rank_fit <- function(x,y) {
  z <- cbind(x,1); sign_y <- 2*y-1
  fn <- function(w) sum(w*w)/2+sum(pmax(1-sign_y*drop(z%*%w),0)^2)
  gr <- function(w) w-2*drop(crossprod(z,sign_y*pmax(1-sign_y*drop(z%*%w),0)))
  fit <- optim(rep(0,ncol(z)),fn,gr,method="BFGS",control=list(maxit=20000,reltol=1e-12))
  if (fit$convergence!=0) stop("RFE optimizer did not converge")
  fit$par
}
rank_path <- function(a,ya,b=NULL,yb=NULL,stop_k=1L) {
  alive <- seq_len(ncol(a)); scores <- rep(NA_real_,ncol(a))
  repeat {
    scale <- scaled(a[,alive,drop=FALSE], if(is.null(b)) a[,alive,drop=FALSE] else b[,alive,drop=FALSE])
    w <- linear_rank_fit(scale$train,ya)
    if(!is.null(b)) scores[length(alive)] <- auc(yb,drop(cbind(scale$test,1)%*%w))
    if(length(alive)==stop_k) break
    # Equal importances retain original column order for tie-breaking.
    remove <- order(w[seq_along(alive)]^2,seq_along(alive))[1]
    alive <- alive[-remove]
  }
  list(scores=scores,selected=alive)
}
svm_score <- function(a,ya,b,C) {
  scale <- scaled(a,b)
  model <- svm(x=scale$train,y=factor(ya,levels=c(0,1)),kernel="linear",
               type="C-classification",cost=C,scale=FALSE,tolerance=0.001)
  pred <- predict(model,scale$test,decision.values=TRUE)
  score <- drop(attr(pred,"decision.values"))
  if(model$levels[model$labels[1]]!="1") score <- -score
  score
}
c_grid <- 10^seq(-3,3,length.out=13); svm_grid <- 2^(-4:4)
rows <- list()
for(i in seq_along(folds)) {
  split <- folds[[i]]; tr <- unlist(split$train); te <- unlist(split$test)
  a <- X[tr,,drop=FALSE]; b <- X[te,,drop=FALSE]; ya <- y[tr]; yb <- y[te]
  l1_cv <- matrix(NA_real_,5,length(c_grid))
  rfe_cv <- matrix(NA_real_,5,ncol(a))
  for(j in seq_along(split$inner)) {
    it <- unlist(split$inner[[j]]$train); iv <- unlist(split$inner[[j]]$test)
    sc <- scaled(a[it,,drop=FALSE],a[iv,,drop=FALSE])
    for(k in seq_along(c_grid)) {
      w <- l1_fit(sc$train,ya[it],c_grid[k])
      l1_cv[j,k] <- auc(ya[iv],plogis(drop(cbind(sc$test,1)%*%w)))
    }
    rfe_cv[j,] <- rank_path(a[it,,drop=FALSE],ya[it],a[iv,,drop=FALSE],ya[iv])$scores
  }
  best_C <- c_grid[which.max(colMeans(l1_cv))]
  sc <- scaled(a,b); w <- l1_fit(sc$train,ya,best_C)
  kept_l1 <- genes[abs(w[seq_along(genes)])>1e-10]
  l1_auc <- auc(yb,plogis(drop(cbind(sc$test,1)%*%w)))
  n_keep <- which.max(colMeans(rfe_cv))
  kept <- rank_path(a,ya,stop_k=n_keep)$selected
  svm_cv <- matrix(NA_real_,5,length(svm_grid))
  for(j in seq_along(split$inner)) {
    it <- unlist(split$inner[[j]]$train); iv <- unlist(split$inner[[j]]$test)
    for(k in seq_along(svm_grid)) svm_cv[j,k] <- auc(ya[iv],svm_score(a[it,kept,drop=FALSE],ya[it],a[iv,kept,drop=FALSE],svm_grid[k]))
  }
  svm_C <- svm_grid[which.max(colMeans(svm_cv))]
  svm_auc <- auc(yb,svm_score(a[,kept,drop=FALSE],ya,b[,kept,drop=FALSE],svm_C))
  rows[[i]] <- data.frame(Split=i,LASSO_selected=paste(kept_l1,collapse="; "),LASSO_best_C=best_C,
                          LASSO_outer_AUC=l1_auc,SVM_RFE_selected=paste(genes[kept],collapse="; "),
                          SVM_RFE_n_features=n_keep,SVM_best_C=svm_C,SVM_outer_AUC=svm_auc)
  write.csv(do.call(rbind,rows),file.path(out,"ml_nested_resampling_split_results.csv"),row.names=FALSE)
  cat("Nested R split",i,"/50 completed\n"); flush.console()
}

# Summarize selection counts and dependent outer-fold performance.
split_results <- do.call(rbind,rows)
published <- read_data(file.path(root,"derived_outputs/original_ML_feature_selection.csv"))
as_bool <- function(x) tolower(as.character(x)) %in% c("true","yes","1")
counts <- function(x,g) sum(vapply(strsplit(x,"; ",fixed=TRUE),function(z) g %in% z,logical(1)))
selection <- do.call(rbind,lapply(genes,function(g) {
  l1 <- counts(split_results$LASSO_selected,g)
  svm <- counts(split_results$SVM_RFE_selected,g)
  p <- published[published$Gene==g,,drop=FALSE]
  stopifnot(nrow(p)==1L)
  data.frame(Gene=g, Published_LASSO_selected=as_bool(p$LASSO_selected),
             Published_SVM_RFE_selected=as_bool(p$SVM_RFE_selected),
             Nested_LASSO_selection_count=l1,Nested_LASSO_selection_frequency=l1/nrow(split_results),
             Nested_SVM_RFE_selection_count=svm,Nested_SVM_RFE_selection_frequency=svm/nrow(split_results),
             Both_nested_frequencies_ge_0_50=l1/nrow(split_results)>=0.5 & svm/nrow(split_results)>=0.5)
}))
selection <- selection[order(-selection$Both_nested_frequencies_ge_0_50,
                            -selection$Nested_LASSO_selection_frequency,
                            -selection$Nested_SVM_RFE_selection_frequency),,drop=FALSE]
note <- "Conditional evaluation of a seven-gene panel selected before resampling; scaling, selection and tuning use outer training sets. Outer folds are dependent resamples."
performance_row <- function(method, a) data.frame(Method=method,Outer_split_count=length(a),
  Mean_outer_AUC=mean(a),SD_outer_AUC=sd(a),Median_outer_AUC=median(a),
  Min_outer_AUC=min(a),Max_outer_AUC=max(a),Evaluation_note=note)
performance <- rbind(performance_row("L1 logistic regression",split_results$LASSO_outer_AUC),
                     performance_row("Linear SVM-RFE",split_results$SVM_outer_AUC))
write_data(selection,file.path(out,"ml_nested_resampling_selection_stability.csv"))
write_data(performance,file.path(out,"ml_nested_resampling_performance.csv"))
diagnostics <- do.call(rbind,solver_checks)
write_data(data.frame(L1_fits=nrow(diagnostics),Convergence_failures=sum(diagnostics[,"convergence"]!=0),
                      Max_KKT_residual=max(diagnostics[,"kkt"])),file.path(out,"solver_convergence.csv"))
summary <- list(sample_count=length(y),control_count=sum(y==0),tumor_count=sum(y==1),
                feature_count=length(genes),outer_split_count=nrow(split_results),
                genes_with_both_selection_frequencies_ge_0_50=selection$Gene[selection$Both_nested_frequencies_ge_0_50],
                mean_outer_auc_lasso=mean(split_results$LASSO_outer_AUC),
                mean_outer_auc_svm=mean(split_results$SVM_outer_AUC))
jsonlite::write_json(summary,file.path(out,"nested_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(performance[,1:4])

