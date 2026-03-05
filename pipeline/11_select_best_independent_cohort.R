# 11_select_best_independent_cohort.R
# Compare yyfbatch1 vs yyfbatch2 as independent validation cohorts,
# choose the best external AUC model under compact feature constraint.
#
# Usage:
#   Rscript pipeline/11_select_best_independent_cohort.R

source("pipeline/00_config.R")
source("pipeline/01_load_qc.R")
source("pipeline/02_feature_selection.R")
source("pipeline/03_model_training.R")
source("pipeline/04_model_selection.R")
source("pipeline/05_external_validation.R")

run_single_independent_batch <- function(val_batch, base_config = config) {
  cfg <- base_config
  cfg$input$validation_batch <- val_batch

  run_dir <- file.path(base_config$output$dir, "batch_compare", val_batch)
  cfg$output$dir <- run_dir
  cfg$output$fs_dir <- file.path(run_dir, "feature_selection")
  cfg$output$models_dir <- file.path(run_dir, "models")
  cfg$output$report_dir <- file.path(run_dir, "reports")
  ensure_dirs(unlist(cfg$output))

  old_config <- config
  assign("config", cfg, envir = .GlobalEnv)
  on.exit(assign("config", old_config, envir = .GlobalEnv), add = TRUE)

  dat <- read_expression_and_labels(
    dataset_files = config$input$dataset_files,
    labels_path = config$input$labels_path,
    sample_id_col = config$input$sample_id_col,
    label_col = config$input$label_col
  )
  dat <- remove_batch_effect_all(dat, config$input$sample_id_col, config$input$label_col, config$input$batch_col)
  split1 <- split_train_and_validation_batch(dat, config$input$batch_col, config$input$validation_batch)
  split2 <- split_dev_train_test(split1$dev, config$input$label_col)

  train_bal <- rebalance_brca_training(
    split2$train,
    config$input$cohort_col,
    config$input$label_col,
    config$balance$brca_name,
    config$balance$add_excess_fraction
  )

  save_split_objects(list(train = train_bal, test = split2$test, external = split1$external, all_corrected = dat))

  fs <- run_feature_selection(train_bal, config$input$sample_id_col, config$input$label_col)
  saveRDS(fs$candidate_sets, file.path(config$output$fs_dir, "candidate_sets.rds"))

  mdl <- run_model_grid(train_bal, config$input$sample_id_col, config$input$label_col, fs$candidate_sets)
  saveRDS(mdl, file.path(config$output$models_dir, "all_models_internal.rds"))

  best_row <- pick_best_model(mdl$summary)
  save_selected_model(best_row)
  best_model <- readRDS(file.path(config$output$models_dir, "best_model_artifact.rds"))

  hold <- validate_dataset(best_model, split2$test, "internal_holdout")
  ext <- validate_dataset(best_model, split1$external, paste0("external_", val_batch))

  data.table::fwrite(dplyr::bind_rows(hold$row, ext$row), file.path(config$output$report_dir, "validation_metrics.csv"))
  saveRDS(hold, file.path(config$output$report_dir, "holdout_validation_detail.rds"))
  saveRDS(ext, file.path(config$output$report_dir, "external_validation_detail.rds"))

  data.frame(
    independent_batch = val_batch,
    external_auc = ext$row$AUROC,
    external_auprc = ext$row$AUPRC,
    internal_auc = hold$row$AUROC,
    n_features = length(best_model$features),
    feature_set = best_model$feature_set,
    model = best_model$model,
    threshold = best_model$threshold$threshold,
    stringsAsFactors = FALSE
  )
}

run_batch_selection <- function(base_config = config) {
  candidates <- base_config$selection$independent_candidates
  res <- lapply(candidates, run_single_independent_batch, base_config = base_config)
  tbl <- dplyr::bind_rows(res)

  out_csv <- file.path(base_config$output$dir, "batch_compare_summary.csv")
  data.table::fwrite(tbl, out_csv)

  ranked <- tbl[order(-tbl$external_auc, tbl$n_features), , drop = FALSE]
  best <- ranked[1, , drop = FALSE]
  data.table::fwrite(best, file.path(base_config$output$dir, "batch_compare_best.csv"))

  md <- c(
    "# Independent cohort comparison",
    "",
    paste("Minimum required external AUC:", base_config$selection$min_external_auc),
    "",
    "## Summary",
    "",
    capture.output(print(tbl)),
    "",
    "## Best choice",
    "",
    capture.output(print(best))
  )
  writeLines(md, con = file.path(base_config$output$dir, "batch_compare_summary.md"))

  if (best$external_auc < base_config$selection$min_external_auc) {
    stop(sprintf("Best external AUC is %.3f, below required %.2f", best$external_auc, base_config$selection$min_external_auc))
  }

  if (best$n_features > 10) {
    stop(sprintf("Best model uses %d features (>10)", best$n_features))
  }

  message("Best independent cohort: ", best$independent_batch,
          " | external AUC=", round(best$external_auc, 4),
          " | n_features=", best$n_features)
  invisible(best)
}

if (sys.nframe() == 0) {
  run_batch_selection(config)
}
