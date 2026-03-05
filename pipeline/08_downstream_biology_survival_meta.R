# 08_downstream_biology_survival_meta.R
# Downstream analyses after selecting compact piRNA features:
# 1) Pearson correlation to mRNAs + Reactome enrichment
# 2) piRNA-gene-pathway network
# 3) stratified ROC + T-score boxplots with Mann-Whitney tests
# 4) KM curves for each biomarker and composite score
# 5) multivariable Cox
# 6) cross-dataset expression heatmap + meta-analysis forest + correlation heatmap
# 7) KEGG/GO enrichment figures

source("pipeline/00_config.R")

get_signature_features <- function(model_artifact_path = file.path(config$output$models_dir, "best_model_artifact.rds")) {
  art <- readRDS(model_artifact_path)
  unique(art$features)
}

prepare_long_expression <- function(df, id_col, label_col, signature_features) {
  keep <- c(id_col, label_col, signature_features)
  keep <- keep[keep %in% colnames(df)]
  x <- df[, keep, drop = FALSE]
  long <- tidyr::pivot_longer(x, cols = tidyselect::all_of(signature_features), names_to = "feature", values_to = "expr")
  long
}

pearson_correlated_genes <- function(expr_df, signature_features, min_abs_cor = 0.3, p_cut = 0.05) {
  # expr_df is expected to include piRNA + mRNA expression columns in the same samples.
  all_features <- colnames(expr_df)
  genes <- setdiff(all_features, signature_features)
  out <- list()

  for (p in signature_features) {
    if (!(p %in% all_features)) next
    rows <- lapply(genes, function(g) {
      test <- suppressWarnings(stats::cor.test(expr_df[[p]], expr_df[[g]], method = "pearson"))
      data.frame(piRNA = p, gene = g, cor = unname(test$estimate), p_value = test$p.value)
    })
    tbl <- dplyr::bind_rows(rows)
    tbl$FDR <- p.adjust(tbl$p_value, method = "BH")
    tbl <- tbl %>% dplyr::filter(abs(cor) >= min_abs_cor, FDR <= p_cut)
    out[[p]] <- tbl
  }

  dplyr::bind_rows(out)
}

reactome_enrichment <- function(cor_tbl) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    warning("ReactomePA not installed; skipping Reactome enrichment")
    return(NULL)
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning("clusterProfiler not installed; skipping Reactome enrichment")
    return(NULL)
  }

  genes <- unique(cor_tbl$gene)
  enr <- ReactomePA::enrichPathway(gene = genes, pAdjustMethod = "BH", readable = TRUE)
  as.data.frame(enr)
}

plot_pirna_gene_pathway_network <- function(cor_tbl, reactome_tbl, out_file) {
  if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
    warning("ggraph/igraph not installed; skipping network plot")
    return(invisible(NULL))
  }

  if (is.null(reactome_tbl) || nrow(reactome_tbl) == 0 || nrow(cor_tbl) == 0) return(invisible(NULL))

  top_pathways <- reactome_tbl %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 10)
  path2gene <- top_pathways %>%
    dplyr::transmute(pathway = Description, gene = strsplit(geneID, "/")) %>%
    tidyr::unnest(gene)

  edges_pg <- cor_tbl %>% dplyr::distinct(piRNA, gene) %>% dplyr::rename(from = piRNA, to = gene)
  edges_gp <- path2gene %>% dplyr::rename(from = gene, to = pathway)
  edges <- dplyr::bind_rows(edges_pg, edges_gp)

  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  p <- ggraph::ggraph(graph, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.25) +
    ggraph::geom_node_point(ggplot2::aes(color = ifelse(name %in% unique(cor_tbl$piRNA), "piRNA", ifelse(name %in% unique(path2gene$pathway), "Pathway", "Gene"))), size = 3) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 2.7) +
    ggplot2::theme_void() +
    ggplot2::labs(title = "piRNA-gene-pathway network")

  ggplot2::ggsave(out_file, p, width = 10, height = 8, dpi = 300)
}

stratified_roc_and_box <- function(df, label_col, prob_col, score_col, strata_cols, out_prefix = "strat") {
  for (scol in strata_cols) {
    if (!(scol %in% colnames(df))) next

    # ROC by subgroup
    auc_rows <- list()
    roc_curve <- list()
    for (g in na.omit(unique(df[[scol]]))) {
      tmp <- df[df[[scol]] == g, , drop = FALSE]
      if (length(unique(tmp[[label_col]])) < 2) next
      r <- pROC::roc(tmp[[label_col]], tmp[[prob_col]], levels = rev(levels(tmp[[label_col]])), quiet = TRUE)
      ci <- pROC::ci.auc(r)
      auc_rows[[as.character(g)]] <- data.frame(group = g, AUC = as.numeric(pROC::auc(r)), CI_low = ci[1], CI_high = ci[3], n = nrow(tmp))
      roc_curve[[as.character(g)]] <- data.frame(group = g, fpr = 1 - r$specificities, tpr = r$sensitivities)
    }

    auc_tbl <- dplyr::bind_rows(auc_rows)
    curve_tbl <- dplyr::bind_rows(roc_curve)
    data.table::fwrite(auc_tbl, file.path(config$output$report_dir, paste0(out_prefix, "_", scol, "_auc_ci.csv")))

    if (nrow(curve_tbl) > 0) {
      p_roc <- ggplot2::ggplot(curve_tbl, ggplot2::aes(fpr, tpr, color = group)) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste0("ROC by ", scol), x = "FPR", y = "TPR")
      ggplot2::ggsave(file.path(config$output$report_dir, paste0(out_prefix, "_", scol, "_roc.png")), p_roc, width = 7, height = 5, dpi = 300)
    }

    # T-score boxplot + Mann-Whitney per stratum level
    pvals <- lapply(na.omit(unique(df[[scol]])), function(g) {
      tmp <- df[df[[scol]] == g, , drop = FALSE]
      if (length(unique(tmp[[label_col]])) < 2) return(NULL)
      w <- suppressWarnings(stats::wilcox.test(tmp[[score_col]] ~ tmp[[label_col]]))
      data.frame(group = g, p_value = w$p.value, n = nrow(tmp))
    })
    pval_tbl <- dplyr::bind_rows(pvals)
    data.table::fwrite(pval_tbl, file.path(config$output$report_dir, paste0(out_prefix, "_", scol, "_mann_whitney.csv")))

    p_box <- ggplot2::ggplot(df, ggplot2::aes_string(x = scol, y = score_col, fill = label_col)) +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.75), outlier.alpha = 0.45) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste0("T-score tumor vs control by ", scol), x = scol, y = "T-score")
    ggplot2::ggsave(file.path(config$output$report_dir, paste0(out_prefix, "_", scol, "_tscore_box.png")), p_box, width = 9, height = 5, dpi = 300)
  }
}

run_km_univariable <- function(df, time_col, event_col, biomarkers, score_bin_col) {
  if (!requireNamespace("survminer", quietly = TRUE)) {
    warning("survminer not installed; KM plots skipped")
    return(invisible(NULL))
  }

  for (b in biomarkers) {
    if (!(b %in% colnames(df))) next
    cut <- stats::median(df[[b]], na.rm = TRUE)
    grp <- factor(ifelse(df[[b]] >= cut, "High", "Low"))
    sfit <- survival::survfit(survival::Surv(df[[time_col]], df[[event_col]]) ~ grp)
    p <- survminer::ggsurvplot(sfit, data = df, risk.table = TRUE, pval = TRUE, conf.int = FALSE, palette = c("#2b8cbe", "#de2d26"), title = paste0("KM: ", b))
    ggplot2::ggsave(file.path(config$output$report_dir, paste0("km_", b, ".png")), p$plot, width = 7, height = 5, dpi = 300)
  }

  if (score_bin_col %in% colnames(df)) {
    sfit2 <- survival::survfit(survival::Surv(df[[time_col]], df[[event_col]]) ~ df[[score_bin_col]])
    p2 <- survminer::ggsurvplot(sfit2, data = df, risk.table = TRUE, pval = TRUE, conf.int = FALSE, palette = c("#3182bd", "#ffd92f"), title = "KM: composite risk score")
    ggplot2::ggsave(file.path(config$output$report_dir, "km_composite_score.png"), p2$plot, width = 7, height = 5, dpi = 300)
  }
}

run_multivariable_cox <- function(df, time_col, event_col, covariates) {
  covariates <- covariates[covariates %in% colnames(df)]
  if (length(covariates) == 0) return(NULL)
  cdf <- df[, c(time_col, event_col, covariates), drop = FALSE]
  cdf <- cdf[stats::complete.cases(cdf), , drop = FALSE]
  f <- stats::as.formula(paste0("survival::Surv(", time_col, ",", event_col, ") ~ ", paste(covariates, collapse = " + ")))
  fit <- survival::coxph(f, data = cdf)
  tbl <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  data.table::fwrite(tbl, file.path(config$output$report_dir, "multivariable_cox_results.csv"))
  tbl
}

plot_cross_dataset_heatmap <- function(df, signature_features, dataset_col, label_col, out_file) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed; heatmap skipped")
    return(invisible(NULL))
  }
  use_cols <- c(dataset_col, label_col, signature_features)
  use_cols <- use_cols[use_cols %in% colnames(df)]
  mat <- as.matrix(df[, signature_features[signature_features %in% colnames(df)], drop = FALSE])
  ann <- df[, c(dataset_col, label_col), drop = FALSE]
  rownames(ann) <- seq_len(nrow(ann))
  pheatmap::pheatmap(mat, annotation_row = ann, scale = "row", show_rownames = FALSE, filename = out_file, width = 9, height = 8)
}

run_smd_meta <- function(df, feature, label_col, dataset_col) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    warning("metafor not installed; meta-analysis skipped")
    return(NULL)
  }
  dsets <- na.omit(unique(df[[dataset_col]]))
  rows <- lapply(dsets, function(ds) {
    tmp <- df[df[[dataset_col]] == ds, , drop = FALSE]
    if (length(unique(tmp[[label_col]])) < 2) return(NULL)
    es <- metafor::escalc(measure = "SMD", m1i = mean(tmp[tmp[[label_col]] == config$input$positive_class, feature], na.rm = TRUE),
                          sd1i = stats::sd(tmp[tmp[[label_col]] == config$input$positive_class, feature], na.rm = TRUE),
                          n1i = sum(tmp[[label_col]] == config$input$positive_class),
                          m2i = mean(tmp[tmp[[label_col]] == config$input$negative_class, feature], na.rm = TRUE),
                          sd2i = stats::sd(tmp[tmp[[label_col]] == config$input$negative_class, feature], na.rm = TRUE),
                          n2i = sum(tmp[[label_col]] == config$input$negative_class))
    data.frame(dataset = ds, yi = es$yi, vi = es$vi)
  })
  est <- dplyr::bind_rows(rows)
  if (nrow(est) < 2) return(NULL)
  fit <- metafor::rma(yi, vi, data = est, method = "REML")
  png(file.path(config$output$report_dir, paste0("meta_forest_", feature, ".png")), width = 950, height = 650)
  metafor::forest(fit, slab = est$dataset, xlab = "SMD (tumor vs normal)")
  dev.off()
  est
}

if (interactive()) {
  message("Template script: provide project-specific inputs before full downstream run.")
}
