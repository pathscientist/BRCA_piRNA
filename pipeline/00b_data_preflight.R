# 00b_data_preflight.R
# One-command preflight checker for real-data readiness before running the pipeline.
# Usage:
#   Rscript pipeline/00b_data_preflight.R

source("pipeline/00_config.R")

resolve_dataset_path <- function(p) {
  candidates <- unique(c(
    p,
    file.path("processed_results", p),
    file.path(getwd(), p),
    file.path(getwd(), "processed_results", p)
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

validate_single_dataset <- function(path, cfg) {
  resolved <- resolve_dataset_path(path)
  out <- list(
    dataset = path,
    resolved_path = resolved,
    exists = !is.na(resolved),
    n_rows = NA_integer_,
    n_cols = NA_integer_,
    missing_required = character(0),
    missing_recommended = character(0),
    notes = character(0)
  )

  if (is.na(resolved)) {
    out$notes <- c(out$notes, "File not found")
    return(out)
  }

  df <- tryCatch(data.table::fread(resolved, data.table = FALSE), error = function(e) e)
  if (inherits(df, "error")) {
    out$notes <- c(out$notes, paste("Read error:", df$message))
    return(out)
  }

  out$n_rows <- nrow(df)
  out$n_cols <- ncol(df)

  required <- c(cfg$input$sample_id_col, cfg$input$label_col)
  recommended <- c(cfg$input$batch_col, cfg$input$cohort_col, "age", "stage", "subtype", "OS_time", "OS_event")

  out$missing_required <- setdiff(required, colnames(df))
  out$missing_recommended <- setdiff(recommended, colnames(df))

  # detect expression columns as all columns not in known metadata set
  known_meta <- unique(c(required, recommended, cfg$input$optional_covariates))
  expr_cols <- setdiff(colnames(df), known_meta)
  if (length(expr_cols) == 0) {
    out$notes <- c(out$notes, "No expression columns detected after excluding metadata columns")
  }

  # label sanity
  if (cfg$input$label_col %in% colnames(df)) {
    lv <- unique(as.character(df[[cfg$input$label_col]]))
    if (!all(c(cfg$input$positive_class, cfg$input$negative_class) %in% lv)) {
      out$notes <- c(out$notes, paste0("Label values do not include both '", cfg$input$negative_class, "' and '", cfg$input$positive_class, "'"))
    }
  }

  out
}

validate_labels_mapping <- function(cfg) {
  p <- cfg$input$labels_path
  if (is.null(p) || identical(p, "")) {
    return(list(enabled = FALSE, exists = FALSE, missing_cols = character(0), note = "labels_path is NULL/empty (expected when labels are in each dataset)"))
  }

  rp <- resolve_dataset_path(p)
  if (is.na(rp)) {
    return(list(enabled = TRUE, exists = FALSE, missing_cols = character(0), note = "labels mapping file configured but not found"))
  }

  lbl <- tryCatch(data.table::fread(rp, data.table = FALSE), error = function(e) e)
  if (inherits(lbl, "error")) {
    return(list(enabled = TRUE, exists = TRUE, missing_cols = character(0), note = paste("labels mapping read error:", lbl$message)))
  }

  req <- c(cfg$input$sample_id_col, cfg$input$label_col)
  miss <- setdiff(req, colnames(lbl))
  list(enabled = TRUE, exists = TRUE, missing_cols = miss, note = ifelse(length(miss) == 0, "OK", "Missing required columns in labels mapping"))
}

write_preflight_reports <- function(results, labels_check, out_dir) {
  rows <- lapply(results, function(x) {
    data.frame(
      dataset = x$dataset,
      resolved_path = ifelse(is.na(x$resolved_path), "", x$resolved_path),
      exists = x$exists,
      n_rows = x$n_rows,
      n_cols = x$n_cols,
      missing_required = paste(x$missing_required, collapse = ";"),
      missing_recommended = paste(x$missing_recommended, collapse = ";"),
      notes = paste(x$notes, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  tbl <- dplyr::bind_rows(rows)

  data.table::fwrite(tbl, file.path(out_dir, "data_preflight_report.csv"))

  md <- c(
    "# Data Preflight Report",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Dataset checks",
    "",
    "|dataset|exists|n_rows|n_cols|missing_required|missing_recommended|notes|",
    "|---|---:|---:|---:|---|---|---|"
  )

  for (i in seq_len(nrow(tbl))) {
    md <- c(md, paste0("|", tbl$dataset[i], "|", tbl$exists[i], "|", tbl$n_rows[i], "|", tbl$n_cols[i], "|", tbl$missing_required[i], "|", tbl$missing_recommended[i], "|", gsub("\\|", "/", tbl$notes[i]), "|"))
  }

  md <- c(md,
          "",
          "## labels_path check",
          "",
          paste("- enabled:", labels_check$enabled),
          paste("- exists:", labels_check$exists),
          paste("- missing_cols:", paste(labels_check$missing_cols, collapse = ",")),
          paste("- note:", labels_check$note))

  writeLines(md, con = file.path(out_dir, "data_preflight_report.md"))
}

run_preflight <- function(cfg = config) {
  ensure_dirs(unlist(cfg$output))
  results <- lapply(cfg$input$dataset_files, validate_single_dataset, cfg = cfg)
  labels_check <- validate_labels_mapping(cfg)

  write_preflight_reports(results, labels_check, cfg$output$report_dir)

  hard_fail <- any(vapply(results, function(x) !x$exists || length(x$missing_required) > 0, logical(1))) ||
    (labels_check$enabled && (!labels_check$exists || length(labels_check$missing_cols) > 0))

  if (hard_fail) {
    message("Preflight FAILED. See artifacts/reports/data_preflight_report.md")
    return(invisible(FALSE))
  }

  message("Preflight PASSED. See artifacts/reports/data_preflight_report.md")
  invisible(TRUE)
}

if (sys.nframe() == 0) {
  ok <- run_preflight(config)
  if (!isTRUE(ok)) quit(save = "no", status = 1)
}
