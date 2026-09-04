## STRING PPI network topology analysis for PS-MPs and thyroid cancer project
## Main network: submitted 75 proxy-associated/thyroid-cancer shared targets
## Highlight nodes: 7 DEG-compound-disease candidate genes and final core genes BAX/BCL2/FN1
## R version used during analysis: R 4.3.3
## Date: 2026-06-12

options(stringsAsFactors = FALSE)
options(timeout = 600)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggrepel)
  library(openxlsx)
  library(httr)
})

cli_args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "scripts/03_STRING_PPI_topology.R"
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input_gene_file <- normalizePath(
  if (length(cli_args) >= 1L) cli_args[[1L]] else file.path(package_root, "data", "processed", "shared_targets_75.csv"),
  winslash = "/",
  mustWork = TRUE
)
out_dir <- normalizePath(
  if (length(cli_args) >= 2L) cli_args[[2L]] else file.path(package_root, "reproduced_outputs", "STRING_submitted_75"),
  winslash = "/",
  mustWork = FALSE
)
result_dir <- file.path(out_dir, "results")
figure_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(log_dir, "STRING_PPI_topology_session.log"), split = TRUE)
cat("STRING PPI topology analysis\n")
cat("Input gene list:", basename(input_gene_file), "\n")
cat("Output directory name:", basename(out_dir), "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:", R.version.string, "\n\n")

read_gene_list <- function(path) {
  if (!file.exists(path)) {
    stop("Input gene file not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv")) {
    sep <- ifelse(ext == "csv", ",", "\t")
    x <- fread(path, sep = sep, header = TRUE)
    if ("Gene" %in% names(x)) {
      genes <- x[["Gene"]]
    } else {
      genes <- unlist(x[[1]], use.names = FALSE)
    }
  } else {
    genes <- readLines(path, warn = FALSE)
  }
  genes <- trimws(genes)
  genes <- genes[genes != "" & toupper(genes) != "GENE"]
  sort(unique(toupper(genes)))
}

safe_download_tsv <- function(url, destfile) {
  ok <- FALSE
  last_error <- NULL
  for (i in seq_len(3)) {
    tryCatch({
      resp <- httr::GET(
        url,
        httr::timeout(180),
        httr::user_agent("PS_MPs_thyroid_cancer_PPI_topology_R4.3.3")
      )
      if (httr::http_error(resp)) {
        stop("HTTP ", httr::status_code(resp), ": ", httr::content(resp, as = "text", encoding = "UTF-8"))
      }
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      if (!nzchar(txt)) {
        stop("Empty response from STRING API")
      }
      writeLines(txt, destfile, useBytes = TRUE)
      ok <<- file.exists(destfile) && file.info(destfile)$size > 0
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      Sys.sleep(2)
    })
    if (ok) break
  }
  if (!ok) {
    stop("Failed to download STRING API result: ", url, "\nLast error: ", last_error)
  }
  fread(destfile, sep = "\t", header = TRUE, quote = "")
}

string_post_tsv <- function(endpoint, body, destfile) {
  ok <- FALSE
  last_error <- NULL
  url <- paste0("https://string-db.org/api/tsv/", endpoint)
  for (i in seq_len(3)) {
    attempt <- tryCatch({
      resp <- httr::POST(
        url,
        body = body,
        encode = "form",
        httr::timeout(180),
        httr::user_agent("PS_MPs_thyroid_cancer_PPI_topology_R4.3.3")
      )
      if (httr::http_error(resp)) {
        stop("HTTP ", httr::status_code(resp), ": ", httr::content(resp, as = "text", encoding = "UTF-8"))
      }
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      if (!nzchar(txt)) {
        stop("Empty response from STRING API endpoint: ", endpoint)
      }
      writeLines(txt, destfile, useBytes = TRUE)
      list(ok = file.exists(destfile) && isTRUE(file.info(destfile)$size > 0), error = NULL)
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (!nzchar(msg)) {
        msg <- paste("Unknown error in STRING request attempt", i, "for endpoint", endpoint)
      }
      list(ok = FALSE, error = msg)
    })
    ok <- isTRUE(attempt$ok)
    if (!ok) {
      last_error <- attempt$error
      if (!nzchar(last_error)) {
        last_error <- paste("Unknown error in STRING request attempt", i, "for endpoint", endpoint)
      }
      Sys.sleep(2)
    }
    if (ok) break
  }
  if (!ok) {
    stop("Failed to download STRING API result from endpoint: ", endpoint, "\nLast error: ", last_error)
  }
  fread(destfile, sep = "\t", header = TRUE, quote = "")
}

string_network <- function(genes, label, required_score = 700) {
  genes <- sort(unique(toupper(genes)))
  identifiers <- paste(genes, collapse = "\r")
  caller <- "PS_MPs_thyroid_cancer_PPI_topology_R4.3.3"

  mapping <- string_post_tsv(
    "get_string_ids",
    body = list(
      identifiers = identifiers,
      species = 9606,
      limit = 1,
      echo_query = 1,
      caller_identity = caller
    ),
    file.path(result_dir, paste0(label, "_STRING_mapping.tsv"))
  )
  raw_edges <- string_post_tsv(
    "network",
    body = list(
      identifiers = identifiers,
      species = 9606,
      required_score = required_score,
      network_type = "functional",
      caller_identity = caller
    ),
    file.path(result_dir, paste0(label, "_STRING_network_required_score_", required_score, ".tsv"))
  )

  mapping[, preferredName := toupper(preferredName)]
  mapping[, queryItem := toupper(queryItem)]
  unmatched <- setdiff(genes, mapping$preferredName)
  if (length(unmatched) > 0) {
    fwrite(data.table(unmatched_gene = unmatched),
           file.path(result_dir, paste0(label, "_unmatched_input_genes.csv")))
  }

  if (nrow(raw_edges) == 0) {
    vertices <- data.table(name = genes)
    edges <- data.table(from = character(), to = character(), combined_score = numeric())
  } else {
    edges <- data.table(
      from = toupper(raw_edges$preferredName_A),
      to = toupper(raw_edges$preferredName_B),
      stringId_A = raw_edges$stringId_A,
      stringId_B = raw_edges$stringId_B,
      combined_score = as.numeric(raw_edges$score),
      nscore = as.numeric(raw_edges$nscore),
      fscore = as.numeric(raw_edges$fscore),
      pscore = as.numeric(raw_edges$pscore),
      ascore = as.numeric(raw_edges$ascore),
      escore = as.numeric(raw_edges$escore),
      dscore = as.numeric(raw_edges$dscore),
      tscore = as.numeric(raw_edges$tscore)
    )
    edges <- unique(edges[from != to])
    vertices <- data.table(name = sort(unique(c(genes, edges$from, edges$to))))
  }

  fwrite(edges, file.path(result_dir, paste0(label, "_PPI_edges.csv")))
  fwrite(vertices, file.path(result_dir, paste0(label, "_PPI_nodes.csv")))

  list(label = label, genes = genes, mapping = mapping, edges = edges, vertices = vertices)
}

factorial_safe <- function(x) {
  ifelse(x <= 1, 1, gamma(x + 1))
}

mcc_score <- function(g) {
  cliques <- max_cliques(g, min = 2)
  scores <- setNames(rep(0, vcount(g)), V(g)$name)
  if (length(cliques) == 0) return(scores)
  for (clq in cliques) {
    members <- V(g)$name[clq]
    contribution <- factorial_safe(length(members) - 1)
    scores[members] <- scores[members] + contribution
  }
  scores
}

neighbor_connectivity <- function(g) {
  vals <- setNames(rep(0, vcount(g)), V(g)$name)
  for (v in V(g)) {
    nb <- neighbors(g, v)
    if (length(nb) < 2) {
      vals[V(g)$name[v]] <- length(nb)
    } else {
      subg <- induced_subgraph(g, nb)
      vals[V(g)$name[v]] <- max(components(subg)$csize)
    }
  }
  vals
}

topology_metrics <- function(network_obj) {
  label <- network_obj$label
  edges <- copy(network_obj$edges)
  vertices <- copy(network_obj$vertices)

  if (nrow(edges) > 0) {
    g <- graph_from_data_frame(
      d = edges[, .(from, to, combined_score)],
      directed = FALSE,
      vertices = vertices
    )
    E(g)$weight <- edges$combined_score
  } else {
    g <- make_empty_graph(n = nrow(vertices), directed = FALSE)
    V(g)$name <- vertices$name
  }

  comp <- components(g)
  deg <- degree(g, mode = "all")
  strength_score <- if (ecount(g) > 0) strength(g, weights = E(g)$weight) else deg
  bet <- if (ecount(g) > 0) betweenness(g, directed = FALSE, weights = NULL, normalized = TRUE) else deg * 0
  clo <- if (ecount(g) > 0) closeness(g, mode = "all", weights = NULL, normalized = TRUE) else deg * 0
  eig <- if (ecount(g) > 0) eigen_centrality(g, directed = FALSE, weights = E(g)$weight)$vector else deg * 0
  pr <- if (ecount(g) > 0) page_rank(g, directed = FALSE, weights = E(g)$weight)$vector else rep(1 / max(vcount(g), 1), vcount(g))
  mcc <- if (ecount(g) > 0) mcc_score(g) else deg * 0
  mnc <- if (ecount(g) > 0) neighbor_connectivity(g) else deg * 0
  clustering <- transitivity(g, type = "local", isolates = "zero")

  metrics <- data.table(
    gene = V(g)$name,
    degree = as.numeric(deg),
    weighted_degree = as.numeric(strength_score),
    betweenness = as.numeric(bet),
    closeness = as.numeric(clo),
    eigenvector = as.numeric(eig),
    pagerank = as.numeric(pr),
    MCC = as.numeric(mcc),
    MNC = as.numeric(mnc),
    clustering_coefficient = as.numeric(clustering),
    component_id = comp$membership,
    component_size = comp$csize[comp$membership]
  )

  rank_cols <- c("degree", "weighted_degree", "betweenness", "closeness", "eigenvector", "pagerank", "MCC", "MNC")
  for (col in rank_cols) {
    metrics[, paste0("rank_", col) := frank(-get(col), ties.method = "min")]
  }
  metrics[, mean_topology_rank := rowMeans(.SD), .SDcols = paste0("rank_", rank_cols)]
  setorder(metrics, mean_topology_rank, -degree, -MCC, gene)

  fwrite(metrics, file.path(result_dir, paste0(label, "_topology_metrics.csv")))

  wb <- createWorkbook()
  addWorksheet(wb, "topology_metrics")
  writeData(wb, "topology_metrics", metrics)
  addWorksheet(wb, "PPI_edges")
  writeData(wb, "PPI_edges", edges)
  addWorksheet(wb, "STRING_mapping")
  writeData(wb, "STRING_mapping", network_obj$mapping)
  saveWorkbook(wb, file.path(result_dir, paste0(label, "_PPI_topology_results.xlsx")), overwrite = TRUE)

  list(label = label, graph = g, metrics = metrics, edges = edges)
}

node_category <- function(genes, candidate_genes, core_genes) {
  category <- rep("Other common target", length(genes))
  category[genes %in% candidate_genes] <- "7-gene candidate"
  category[genes %in% core_genes] <- "Final core gene"
  factor(category, levels = c("Other common target", "7-gene candidate", "Final core gene"))
}

plot_network <- function(topo, candidate_genes, core_genes = c("BAX", "BCL2", "FN1")) {
  label <- topo$label
  g <- topo$graph
  metrics <- topo$metrics
  if (vcount(g) == 0) return(invisible(NULL))

  node_df <- metrics[, .(name = gene, degree, MCC, mean_topology_rank)]
  V(g)$degree <- node_df$degree[match(V(g)$name, node_df$name)]
  V(g)$MCC <- node_df$MCC[match(V(g)$name, node_df$name)]
  V(g)$rank <- node_df$mean_topology_rank[match(V(g)$name, node_df$name)]
  V(g)$category <- node_category(V(g)$name, candidate_genes, core_genes)
  V(g)$top10 <- V(g)$name %in% head(metrics$gene, min(10, nrow(metrics)))
  V(g)$candidate <- V(g)$name %in% candidate_genes
  V(g)$core <- V(g)$name %in% core_genes
  if (ecount(g) > 0 && is.null(E(g)$weight)) {
    E(g)$weight <- 0.7
  }

  set.seed(20260612)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = ifelse(is.na(weight), 0.7, weight)),
                   color = "grey78", alpha = 0.28, show.legend = FALSE) +
    geom_node_point(aes(size = degree, fill = category),
                    shape = 21, color = "white", stroke = 0.6) +
    geom_node_text(
      aes(label = ifelse(core | candidate | top10, name, "")),
      repel = TRUE,
      size = 3.2,
      family = "sans",
      color = "grey15"
    ) +
    scale_fill_manual(values = c(
      "Other common target" = "#4C78A8",
      "7-gene candidate" = "#F28E2B",
      "Final core gene" = "#E45756"
    )) +
    scale_edge_width(range = c(0.12, 0.85)) +
    scale_size_continuous(range = c(3, 11)) +
    labs(
      title = paste0(label, " STRING PPI network"),
      subtitle = "Node size: degree; orange: 7 candidate genes; red: final core genes BAX, BCL2, FN1"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
      legend.position = "bottom",
      legend.title = element_blank()
    )

  ggsave(file.path(figure_dir, paste0(label, "_STRING_PPI_network.png")), p, width = 8.5, height = 7, dpi = 320)
  ggsave(file.path(figure_dir, paste0(label, "_STRING_PPI_network.pdf")), p, width = 8.5, height = 7)
}

plot_topology_bars <- function(topo, candidate_genes, core_genes) {
  label <- topo$label
  metrics <- topo$metrics
  top_n <- min(20, nrow(metrics))
  dat <- copy(metrics[seq_len(top_n)])
  dat[, gene := factor(gene, levels = rev(gene))]
  dat[, category := node_category(as.character(gene), candidate_genes, core_genes)]

  category_colors <- c(
    "Other common target" = "#59A14F",
    "7-gene candidate" = "#F28E2B",
    "Final core gene" = "#E15759"
  )

  p1 <- ggplot(dat, aes(x = gene, y = degree, fill = category)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = category_colors) +
    labs(title = paste0(label, " top genes by integrated topology rank"), x = NULL, y = "Degree") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold"))

  p2 <- ggplot(dat, aes(x = gene, y = MCC, fill = category)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = category_colors) +
    labs(x = NULL, y = "MCC") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())

  ggsave(file.path(figure_dir, paste0(label, "_topology_degree_bar.png")), p1, width = 7.5, height = 5.8, dpi = 320)
  ggsave(file.path(figure_dir, paste0(label, "_topology_MCC_bar.png")), p2, width = 7.5, height = 5.8, dpi = 320)
  ggsave(file.path(figure_dir, paste0(label, "_topology_degree_bar.pdf")), p1, width = 7.5, height = 5.8)
  ggsave(file.path(figure_dir, paste0(label, "_topology_MCC_bar.pdf")), p2, width = 7.5, height = 5.8)
}

plot_metric_heatmap <- function(topo) {
  label <- topo$label
  metrics <- topo$metrics
  top_n <- min(25, nrow(metrics))
  dat <- metrics[seq_len(top_n), .(gene, degree, weighted_degree, betweenness, closeness, eigenvector, pagerank, MCC, MNC)]
  long <- melt(dat, id.vars = "gene", variable.name = "metric", value.name = "value")
  long[, z := as.numeric(scale(value)), by = metric]
  long[is.na(z), z := 0]
  long[, gene := factor(gene, levels = rev(dat$gene))]

  p <- ggplot(long, aes(x = metric, y = gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756", midpoint = 0) +
    labs(title = paste0(label, " topology metrics heatmap"), x = NULL, y = NULL, fill = "Z-score") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(figure_dir, paste0(label, "_topology_metrics_heatmap.png")), p, width = 8.5, height = 6.5, dpi = 320)
  ggsave(file.path(figure_dir, paste0(label, "_topology_metrics_heatmap.pdf")), p, width = 8.5, height = 6.5)
}

make_summary <- function(results, candidate_genes, core_genes) {
  rows <- rbindlist(lapply(results, function(x) {
    g <- x$graph
    data.table(
      network = x$label,
      input_genes = length(x$metrics$gene),
      nodes_in_network = vcount(g),
      edges = ecount(g),
      connected_components = components(g)$no,
      largest_component_size = max(components(g)$csize),
      top_by_integrated_rank = paste(head(x$metrics$gene, min(10, nrow(x$metrics))), collapse = ";"),
      candidate_gene_ranks = paste(
        x$metrics[gene %in% candidate_genes,
                  paste0(gene, "(rank=", round(mean_topology_rank, 2), ",degree=", degree, ",MCC=", MCC, ")")],
        collapse = "; "
      ),
      core_gene_ranks = paste(
        x$metrics[gene %in% core_genes,
                  paste0(gene, "(rank=", round(mean_topology_rank, 2), ",degree=", degree, ",MCC=", MCC, ")")],
        collapse = "; "
      )
    )
  }))
  fwrite(rows, file.path(result_dir, "PPI_network_summary.csv"))
  rows
}

candidate_genes <- c("BAX", "BCL2", "CDKN1A", "COL1A1", "CXCL8", "FN1", "TGFB1")
main_genes <- read_gene_list(input_gene_file)
core_genes <- c("BAX", "BCL2", "FN1")

cat("Main PPI network gene set size:", length(main_genes), "\n")
cat("Candidate genes highlighted within main network:", paste(candidate_genes, collapse = ", "), "\n")
cat("Final core genes:", paste(core_genes, collapse = ", "), "\n\n")

main_network <- string_network(main_genes, "main_75_common_targets", required_score = 700)

main_topo <- topology_metrics(main_network)

plot_network(main_topo, candidate_genes, core_genes)
plot_topology_bars(main_topo, candidate_genes, core_genes)
plot_metric_heatmap(main_topo)

candidate_focus <- main_topo$metrics[gene %in% candidate_genes]
fwrite(candidate_focus, file.path(result_dir, "candidate_7genes_topology_within_main_75_network.csv"))

summary_table <- make_summary(list(main_topo), candidate_genes, core_genes)
print(summary_table)

cat("\nTop genes in the submitted 75-gene PPI network by integrated rank:\n")
print(main_topo$metrics[1:min(20, .N)])

cat("\nCandidate 7 genes within the submitted 75-gene network:\n")
print(candidate_focus)

cat("\nSession info:\n")
print(sessionInfo())
sink()
