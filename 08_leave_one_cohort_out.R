################################################################################
#                                                                              #
#   08 — Leave-One-Cohort-Out (LOCO) Validation                                #
#                                                                              #
#   For each cohort in processed_results/: the other cohorts are the           #
#   discovery set, the held-out cohort is the independent validation set.     #
#   With all 7 cohorts present this is 6 discovery -> 1 validation, 7 folds.   #
#   Missing cohort files are skipped (3 folds with the current data).          #
#                                                                              #
#   No leakage: every step that looks at labels uses discovery cohorts only.   #
#     1. Common piRNAs, log2(TPM+1), per-cohort z-score (label-free).          #
#     2. Screening: expressed (>0) in >= 50% of samples of EVERY discovery     #
#        cohort; Welch t-test per cohort -> signed z; same direction in all    #
#        discovery cohorts; ranked by Stouffer combined z; top K_FEATURES.     #
#     3. Model: L2-penalised logistic regression (class-balanced) on pooled   #
#        discovery cohorts. Threshold = Youden's J on discovery.               #
#     4. Discovery estimate: 5-fold CV with screening redone in every fold.   #
#     5. Validation: locked model + threshold on the held-out cohort.          #
#                                                                              #
#   Produces: results/loco/loco_results.csv                                    #
#             results/loco/Fig_LOCO_validation.png / .pdf                      #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  08: Leave-one-cohort-out validation\n")
cat("=================================================================\n\n")

for (pkg in c("ggplot2", "dplyr", "pROC", "glmnet", "cowplot"))
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(pROC); library(glmnet); library(cowplot)
})

SEED <- 2024; set.seed(SEED)
K_FEATURES <- 10
MIN_EXPR   <- 0.5
OUT <- "results/loco"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COHORTS  <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049",
              "yyfbatch1", "yyfbatch2")
SPECIMEN <- c(yyfbatch1 = "plasma", yyfbatch2 = "plasma")   # all others: tissue
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7a5cc7")

# ==============================================================================
# 1. LOAD
# ==============================================================================
cat("========== STEP 1: Load cohorts ==========\n")
load_cohort <- function(nm) {
  f <- file.path("processed_results", paste0(nm, c("_processed.csv", "_processed_1.csv")))
  f <- f[file.exists(f)]
  if (length(f) == 0) { cat("  missing:", nm, "\n"); return(NULL) }
  df <- read.csv(f[1], check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- make.unique(as.character(df[[1]])); df[[1]] <- NULL
  g <- as.character(df$Group); df$Group <- NULL
  y <- ifelse(g %in% c("Tumor", "Cancer", "cancer", "tumor"), 1,
       ifelse(g %in% c("Normal", "Benign", "benign", "normal", "control", "Control",
                       "healthy", "Healthy"), 0, NA))
  keep <- !is.na(y)
  X <- as.matrix(df[keep, , drop = FALSE]); storage.mode(X) <- "numeric"
  X[is.na(X) | X < 0] <- 0
  list(X = log2(X + 1), y = y[keep])
}
cohorts <- Filter(Negate(is.null), setNames(lapply(COHORTS, load_cohort), COHORTS))
common  <- Reduce(intersect, lapply(cohorts, function(co) colnames(co$X)))
cohorts <- lapply(cohorts, function(co) { co$X <- co$X[, common, drop = FALSE]; co })
cat("  Cohorts:", paste(names(cohorts), collapse = ", "),
    "| common piRNAs:", length(common), "\n")
if (length(cohorts) < 3) stop("Need at least 3 cohorts for leave-one-cohort-out.")

zscore <- function(X) {
  s <- apply(X, 2, sd); s[s == 0 | is.na(s)] <- 1
  sweep(sweep(X, 2, colMeans(X)), 2, s, "/")
}

# ==============================================================================
# 2. HELPERS: screening, fitting
# ==============================================================================
screen <- function(train) {
  feats <- common[Reduce(`&`, lapply(train, function(co) colMeans(co$X > 0) >= MIN_EXPR))]
  if (length(feats) == 0) return(list(feats = character(0), n_expr = 0, n_cons = 0))
  Z <- sapply(train, function(co) {
    X <- co$X[, feats, drop = FALSE]
    apply(X, 2, function(v) {
      a <- v[co$y == 1]; b <- v[co$y == 0]
      if (sd(a) == 0 && sd(b) == 0) return(0)
      tt <- tryCatch(t.test(a, b), error = function(e) NULL)
      if (is.null(tt)) return(0)
      sign(tt$statistic) * qnorm(max(tt$p.value, 1e-300) / 2, lower.tail = FALSE)
    })
  })
  Z <- matrix(Z, nrow = length(feats), dimnames = list(feats, NULL))
  consistent <- apply(Z, 1, function(z) all(z != 0) && length(unique(sign(z))) == 1)
  stouffer <- abs(rowSums(Z)) / sqrt(ncol(Z))
  score <- ifelse(consistent, stouffer, 0)
  score <- sort(score[score > 0], decreasing = TRUE)
  list(feats = head(names(score), K_FEATURES), n_expr = length(feats), n_cons = sum(consistent))
}

pool <- function(train, feats) {
  list(X = do.call(rbind, lapply(train, function(co) zscore(co$X[, feats, drop = FALSE]))),
       y = unlist(lapply(train, `[[`, "y")))
}

fit_model <- function(train, feats) {
  d <- pool(train, feats)
  w <- ifelse(d$y == 1, 0.5 / mean(d$y == 1), 0.5 / mean(d$y == 0))
  # alpha = 0 (ridge), fixed lambda — mirrors C = 0.1 L2 logistic regression
  glmnet(d$X, d$y, family = "binomial", alpha = 0, lambda = 1 / (0.1 * nrow(d$X)),
         weights = w, standardize = FALSE)
}
predict_prob <- function(m, X) as.numeric(predict(m, X, type = "response"))

auc_ci <- function(y, p) {
  r <- roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  as.numeric(ci.auc(r, method = "delong"))[c(2, 1, 3)]   # AUC, low, high
}

# ==============================================================================
# 3. LOCO LOOP
# ==============================================================================
cat("\n========== STEP 2: Leave-one-cohort-out ==========\n")
rows <- list(); roc_df <- list()
for (held in names(cohorts)) {
  train <- cohorts[setdiff(names(cohorts), held)]
  sc <- screen(train)
  if (length(sc$feats) == 0) { cat("  ", held, ": no consistent features\n"); next }

  # discovery estimate: 5-fold CV, screening redone inside each fold
  set.seed(SEED)
  folds <- lapply(train, function(co) {
    f <- integer(length(co$y))
    for (cl in 0:1) { i <- which(co$y == cl); f[i] <- sample(rep(1:5, length.out = length(i))) }
    f
  })
  oof_y <- c(); oof_p <- c()
  for (k in 1:5) {
    tr <- Map(function(co, f) list(X = co$X[f != k, , drop = FALSE], y = co$y[f != k]), train, folds)
    te <- Map(function(co, f) list(X = co$X[f == k, , drop = FALSE], y = co$y[f == k]), train, folds)
    fk <- screen(tr)$feats
    if (length(fk) == 0) next
    m <- fit_model(tr, fk)
    for (co in te) {
      oof_p <- c(oof_p, predict_prob(m, zscore(co$X[, fk, drop = FALSE])))
      oof_y <- c(oof_y, co$y)
    }
  }
  cv <- auc_ci(oof_y, oof_p)

  # final model on all discovery cohorts, locked threshold
  m  <- fit_model(train, sc$feats)
  d  <- pool(train, sc$feats)
  rt <- roc(d$y, predict_prob(m, d$X), levels = c(0, 1), direction = "<", quiet = TRUE)
  cut <- coords(rt, "best", best.method = "youden", ret = "threshold")$threshold[1]

  v  <- cohorts[[held]]
  pv <- predict_prob(m, zscore(v$X[, sc$feats, drop = FALSE]))
  va <- auc_ci(v$y, pv)
  pred <- as.integer(pv >= cut)
  sens <- mean(pred[v$y == 1] == 1); spec <- mean(pred[v$y == 0] == 0)

  rv <- roc(v$y, pv, levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_df[[held]] <- data.frame(cohort = held, fpr = 1 - rv$specificities,
                               tpr = rv$sensitivities)
  rows[[held]] <- data.frame(
    validation = held, specimen = ifelse(held %in% names(SPECIMEN), SPECIMEN[held], "tissue"),
    discovery = paste(names(train), collapse = " + "),
    n_val = length(v$y), n_val_tumor = sum(v$y == 1), n_val_normal = sum(v$y == 0),
    n_expressed = sc$n_expr, n_consistent = sc$n_cons, n_features = length(sc$feats),
    features = paste(sc$feats, collapse = ";"),
    cv_auc = cv[1], cv_lo = cv[2], cv_hi = cv[3],
    val_auc = va[1], val_lo = va[2], val_hi = va[3],
    threshold = cut, sensitivity = sens, specificity = spec, stringsAsFactors = FALSE)
  cat(sprintf("  %-12s discovery CV AUC %.3f | validation AUC %.3f (%.3f-%.3f)\n",
              held, cv[1], va[1], va[2], va[3]))
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(OUT, "loco_results.csv"), row.names = FALSE)
cat("  Saved:", file.path(OUT, "loco_results.csv"), "\n")

# ==============================================================================
# 4. FIGURE
# ==============================================================================
cat("\n========== STEP 3: Figure ==========\n")
roc_all <- do.call(rbind, roc_df)
roc_all$label <- with(res[match(roc_all$cohort, res$validation), ],
                      sprintf("%s (n=%d): AUC %.2f (%.2f-%.2f)",
                              validation, n_val, val_auc, val_lo, val_hi))
cols <- setNames(PAL[seq_len(nrow(res))], unique(roc_all$label))

pA <- ggplot(roc_all, aes(fpr, tpr, color = label)) +
  geom_abline(linetype = "dashed", color = "grey70") +
  geom_step(direction = "vh", linewidth = 1) +
  scale_color_manual(values = cols, name = NULL) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "ROC in each held-out validation cohort") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = c(0.62, 0.16),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", color = NA))

forest <- rbind(
  data.frame(res[, c("validation", "specimen", "discovery", "n_val_tumor", "n_val_normal")],
             set = "Discovery (5-fold CV)", auc = res$cv_auc, lo = res$cv_lo, hi = res$cv_hi),
  data.frame(res[, c("validation", "specimen", "discovery", "n_val_tumor", "n_val_normal")],
             set = "Independent validation", auc = res$val_auc, lo = res$val_lo, hi = res$val_hi))
forest$row <- sprintf("%s (%s; %dT/%dN)\ndiscovery: %s", forest$validation, forest$specimen,
                      forest$n_val_tumor, forest$n_val_normal, forest$discovery)
forest$row <- factor(forest$row, levels = rev(unique(forest$row)))

pB <- ggplot(forest, aes(auc, row, color = set, shape = set)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_text(data = subset(forest, set == "Independent validation"),
            aes(x = 1.04, label = sprintf("%.2f", auc)), color = "grey15", size = 3.2,
            position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_color_manual(values = c("Discovery (5-fold CV)" = "#6b6a64",
                                "Independent validation" = "#2a78d6"), name = NULL) +
  scale_shape_manual(values = c("Discovery (5-fold CV)" = 16,
                                "Independent validation" = 18), name = NULL) +
  scale_x_continuous(limits = c(0.2, 1.08), breaks = seq(0.2, 1, 0.1)) +
  labs(x = "AUC (95% DeLong CI)", y = NULL, title = "Discovery vs held-out AUC per fold") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        axis.text.y = element_text(size = 8.5), panel.grid.minor = element_blank())

title <- ggdraw() + draw_label(sprintf(
  "Leave-one-cohort-out validation (%d cohorts; top-%d direction-consistent piRNAs, logistic regression)",
  nrow(res), K_FEATURES), fontface = "bold", size = 13)
fig <- plot_grid(title, plot_grid(pA, pB, ncol = 2, labels = c("A", "B"),
                                  rel_widths = c(1, 1.25)),
                 ncol = 1, rel_heights = c(0.06, 1))
h <- max(6.5, 1.3 * nrow(res) + 3)
ggsave(file.path(OUT, "Fig_LOCO_validation.png"), fig, width = 15, height = h, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_LOCO_validation.pdf"), fig, width = 15, height = h, bg = "white")
cat("  Saved:", file.path(OUT, "Fig_LOCO_validation.png / .pdf"), "\n")

cat(sprintf("\n=== 08 complete. Runtime: %.1f min ===\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
