# Generate S25, including correction across all three markers within each test.
if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
if (as.character(getRversion()) != "4.3.3") stop("Use R 4.3.3")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else file.path(root,"reproduced_outputs/ihc_statistics")
out <- normalizePath(out,winslash="/",mustWork=FALSE)
for (p in c("data","derived_outputs","scripts","tests","figures")) {
 protected <- paste0(root,"/",p)
 if (out==protected || startsWith(out,paste0(protected,"/"))) stop("Output overlaps protected input")
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dat <- read.csv(file.path(root,"data/processed/ihc_hscore_deidentified.csv"))
stopifnot(nrow(dat)==45L,setequal(unique(dat$Marker),c("BAX","BCL2","FN1")))
result <- do.call(rbind,lapply(c("BAX","BCL2","FN1"),function(gene) {
 d <- dat[dat$Marker==gene,]; a <- d$Con_H_score; b <- d$Tumor_H_score; delta <- b-a
 stopifnot(nrow(d)==15L,!anyDuplicated(d$Pair_ID), all(is.finite(c(a,b))),
  all(c(a,b)>=0 & c(a,b)<=300),setequal(d$Pair_ID,dat$Pair_ID))
 t <- t.test(b,a,paired=TRUE); w <- wilcox.test(b,a,paired=TRUE,exact=TRUE,conf.int=TRUE)
 nonzero <- delta[delta!=0]; ranks <- rank(abs(nonzero))
 rb <- (sum(ranks[nonzero>0])-sum(ranks[nonzero<0]))/sum(ranks)
 data.frame(Gene=gene,Biological_pairs_n=length(a),Adjacent_mean=mean(a),Adjacent_SD=sd(a),
  Tumor_mean=mean(b),Tumor_SD=sd(b),Mean_difference_Tumor_minus_Adjacent=mean(delta),
  Mean_difference_CI_low=t$conf.int[1],Mean_difference_CI_high=t$conf.int[2],
  Cohen_dz=round(mean(delta)/sd(delta),5),Wilcoxon_pseudomedian_difference=unname(w$estimate),
  Wilcoxon_CI_low=w$conf.int[1],Wilcoxon_CI_high=w$conf.int[2],Matched_rank_biserial=round(rb,6),
  Shapiro_P_paired_differences=shapiro.test(delta)$p.value,Paired_t_P=t$p.value,
  Exact_Wilcoxon_P=w$p.value,Direction_in_tumor=ifelse(mean(delta)>0,"Higher","Lower"))
}))
result$Paired_t_Bonferroni_3 <- p.adjust(result$Paired_t_P,method="bonferroni")
result$Exact_Wilcoxon_Bonferroni_3 <- p.adjust(result$Exact_Wilcoxon_P,method="bonferroni")
result <- result[c("Gene","Biological_pairs_n","Adjacent_mean","Adjacent_SD","Tumor_mean","Tumor_SD",
 "Mean_difference_Tumor_minus_Adjacent","Mean_difference_CI_low","Mean_difference_CI_high","Cohen_dz",
 "Wilcoxon_pseudomedian_difference","Wilcoxon_CI_low","Wilcoxon_CI_high","Matched_rank_biserial",
 "Shapiro_P_paired_differences","Paired_t_P","Paired_t_Bonferroni_3","Exact_Wilcoxon_P",
 "Exact_Wilcoxon_Bonferroni_3","Direction_in_tumor")]
write.csv(result,file.path(out,"ihc_paired_statistics.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(result[,c("Gene","Paired_t_Bonferroni_3","Exact_Wilcoxon_Bonferroni_3")])
