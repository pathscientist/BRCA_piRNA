# ==============================================================================
# A) 分三阶段：Discovery / Hold-out / Independent(yyfbatch1)
# ==============================================================================

val_name  <- "yyfbatch1"
test_prop <- 0.30  # Hold-out 比例（从非独立池里抽）

# 确保 Tumor 是 positive class（便于 caret::twoClassSummary）
combat_df_all$Group <- factor(combat_df_all$Group, levels = c("Tumor", "Normal"))

valid_data <- combat_df_all[combat_df_all$Batch == val_name, ]
pool_data  <- combat_df_all[combat_df_all$Batch != val_name, ]

set.seed(2026)
idx_tr <- createDataPartition(pool_data$Group, p = 1 - test_prop, list = FALSE)
discovery_data <- pool_data[idx_tr, ]
test_data      <- pool_data[-idx_tr, ]

cat("Discovery n=", nrow(discovery_data),
    "| Hold-out n=", nrow(test_data),
    "| Independent(", val_name, ") n=", nrow(valid_data), "\n")

# ==============================================================================
# B) 特征筛选（只用 Discovery） + 训练模型（只用 Discovery）
# ==============================================================================

set.seed(123)
rf_fs  <- randomForest(x = discovery_data[, gene_cols], y = discovery_data$Group, ntree = 200)
rf_imp <- importance(rf_fs)
top_feats <- rownames(rf_imp)[order(rf_imp[, "MeanDecreaseGini"], decreasing = TRUE)][1:5]
cat("Top 5:", paste(top_feats, collapse = ", "), "\n")

fitControl <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"     # 用于画“Training(CV)”ROC/PRC
)

set.seed(888)
model <- train(
  Group ~ ., data = discovery_data[, c("Group", top_feats)],
  method = "rf", metric = "ROC",
  trControl = fitControl, ntree = 500
)

# ==============================================================================
# C) 统一：算 ROC/AUC + PRC/AUPRC + bootstrap CI
# ==============================================================================

calc_metrics <- function(y_true01, y_score, n_boot = 1000, seed = 42) {
  roc_obj <- roc(y_true01, y_score, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  pr_obj <- PRROC::pr.curve(
    scores.class0 = y_score[y_true01 == 1],
    scores.class1 = y_score[y_true01 == 0],
    curve = TRUE
  )
  auprc_val <- pr_obj$auc.integral
  
  set.seed(seed)
  boot_auc <- boot_auprc <- numeric(n_boot)
  for (i in 1:n_boot) {
    idx <- sample(seq_along(y_true01), length(y_true01), replace = TRUE)
    yt <- y_true01[idx]; ys <- y_score[idx]
    if (length(unique(yt)) < 2) { boot_auc[i] <- NA; boot_auprc[i] <- NA; next }
    
    boot_auc[i] <- as.numeric(auc(roc(yt, ys, direction = "<", quiet = TRUE)))
    pr_tmp <- PRROC::pr.curve(
      scores.class0 = ys[yt == 1],
      scores.class1 = ys[yt == 0],
      curve = FALSE
    )
    boot_auprc[i] <- pr_tmp$auc.integral
  }
  boot_auc <- na.omit(boot_auc); boot_auprc <- na.omit(boot_auprc)
  
  list(
    roc_obj = roc_obj,
    auc = auc_val,
    auc_ci = quantile(boot_auc, c(0.025, 0.975), na.rm = TRUE),
    pr_obj = pr_obj,
    auprc = auprc_val,
    auprc_ci = quantile(boot_auprc, c(0.025, 0.975), na.rm = TRUE)
  )
}

plot_roc <- function(roc_obj, label, title) {
  df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
  df <- df[order(df$fpr), ]
  ggplot(df, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_path(linewidth = 1.1) +
    annotate("label", x = 0.95, y = 0.05, label = label,
             hjust = 1, vjust = 0, size = 4.2, fill = "white", alpha = 0.85, label.size = NA) +
    scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(title = title, x = "False Positive Rate (1 - Specificity)", y = "Sensitivity") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 13, hjust = 0.5))
}

plot_prc <- function(pr_obj, prevalence, label, title) {
  df <- data.frame(recall = pr_obj$curve[, 1], precision = pr_obj$curve[, 2])
  ggplot(df, aes(recall, precision)) +
    geom_hline(yintercept = prevalence, linetype = "dashed") +
    geom_path(linewidth = 1.1) +
    annotate("label", x = 0.95, y = 0.05, label = label,
             hjust = 1, vjust = 0, size = 4.2, fill = "white", alpha = 0.85, label.size = NA) +
    scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(title = title, x = "Recall", y = "Precision") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 13, hjust = 0.5))
}

dir.create("plots", showWarnings = FALSE)

# ==============================================================================
# D) 三阶段 ROC/PRC：Training 用 CV out-of-fold；Testing/Independent 用外推预测
# ==============================================================================

# ---- Training (CV out-of-fold) ----
pred_cv <- model$pred
# caret 的 pred 可能包含多组调参结果，筛 bestTune
bt <- model$bestTune
for (nm in names(bt)) pred_cv <- pred_cv[pred_cv[[nm]] == bt[[nm]], ]

y_true_tr <- ifelse(pred_cv$obs == "Tumor", 1, 0)
y_score_tr <- pred_cv$Tumor
mt_tr <- calc_metrics(y_true_tr, y_score_tr, n_boot = 1000, seed = 101)

lab_tr_roc <- sprintf("AUC = %.2f (95%% CI %.2f-%.2f)", mt_tr$auc, mt_tr$auc_ci[1], mt_tr$auc_ci[2])
lab_tr_prc <- sprintf("AUPRC = %.2f (95%% CI %.2f-%.2f)", mt_tr$auprc, mt_tr$auprc_ci[1], mt_tr$auprc_ci[2])

p_roc_tr <- plot_roc(mt_tr$roc_obj, lab_tr_roc, "ROC - Discovery (Training, CV)")
p_prc_tr <- plot_prc(mt_tr$pr_obj, mean(y_true_tr), lab_tr_prc, "PRC - Discovery (Training, CV)")

ggsave("plots/ROC_training_CV.png", p_roc_tr, width = 6.5, height = 5, dpi = 300)
ggsave("plots/PRC_training_CV.png", p_prc_tr, width = 6.5, height = 5, dpi = 300)

# ---- Hold-out (Testing) ----
prob_te <- predict(model, test_data[, top_feats], type = "prob")$Tumor
y_true_te <- ifelse(test_data$Group == "Tumor", 1, 0)
mt_te <- calc_metrics(y_true_te, prob_te, n_boot = 1000, seed = 202)

lab_te_roc <- sprintf("AUC = %.2f (95%% CI %.2f-%.2f)", mt_te$auc, mt_te$auc_ci[1], mt_te$auc_ci[2])
lab_te_prc <- sprintf("AUPRC = %.2f (95%% CI %.2f-%.2f)", mt_te$auprc, mt_te$auprc_ci[1], mt_te$auprc_ci[2])

p_roc_te <- plot_roc(mt_te$roc_obj, lab_te_roc, "ROC - Hold-out (Testing)")
p_prc_te <- plot_prc(mt_te$pr_obj, mean(y_true_te), lab_te_prc, "PRC - Hold-out (Testing)")

ggsave("plots/ROC_holdout.png", p_roc_te, width = 6.5, height = 5, dpi = 300)
ggsave("plots/PRC_holdout.png", p_prc_te, width = 6.5, height = 5, dpi = 300)

# ---- Independent Validation (yyfbatch1) ----
prob_iv <- predict(model, valid_data[, top_feats], type = "prob")$Tumor
y_true_iv <- ifelse(valid_data$Group == "Tumor", 1, 0)
mt_iv <- calc_metrics(y_true_iv, prob_iv, n_boot = 1000, seed = 303)

lab_iv_roc <- sprintf("AUC = %.2f (95%% CI %.2f-%.2f)", mt_iv$auc, mt_iv$auc_ci[1], mt_iv$auc_ci[2])
lab_iv_prc <- sprintf("AUPRC = %.2f (95%% CI %.2f-%.2f)", mt_iv$auprc, mt_iv$auprc_ci[1], mt_iv$auprc_ci[2])

p_roc_iv <- plot_roc(mt_iv$roc_obj, lab_iv_roc, paste0("ROC - Independent (", val_name, ")"))
p_prc_iv <- plot_prc(mt_iv$pr_obj, mean(y_true_iv), lab_iv_prc, paste0("PRC - Independent (", val_name, ")"))

ggsave("plots/ROC_independent.png", p_roc_iv, width = 6.5, height = 5, dpi = 300)
ggsave("plots/PRC_independent.png", p_prc_iv, width = 6.5, height = 5, dpi = 300)

print(p_roc_tr); print(p_prc_tr)
print(p_roc_te); print(p_prc_te)
print(p_roc_iv); print(p_prc_iv)

# ==============================================================================
# E) 给每个样本打分：T-score = P(Tumor)*100，并画 3 张 boxplot（按 Batch 分 Tumor/Normal）
# ==============================================================================

score_phase <- function(df, phase) {
  p <- predict(model, df[, top_feats], type = "prob")$Tumor
  data.frame(
    Sample = rownames(df),
    Phase  = phase,
    Batch  = df$Batch,
    Group  = df$Group,
    T_score = 100 * p,
    stringsAsFactors = FALSE
  )
}

score_df <- bind_rows(
  score_phase(discovery_data, "Discovery"),
  score_phase(test_data,      "Hold-out"),
  score_phase(valid_data,     "Independent")
)

plot_box_phase <- function(df, phase_name) {
  df$Batch <- factor(df$Batch, levels = sort(unique(df$Batch)))
  ggplot(df, aes(x = Batch, y = T_score, fill = Group)) +
    geom_boxplot(outlier.size = 0.6, width = 0.7) +
    labs(
      title = paste0("T-scores by Dataset (", phase_name, ")\nTumor vs Normal"),
      x = "Dataset (Batch)",
      y = "T-score = P(Tumor) × 100"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title  = element_text(size = 14, face = "bold"),
      plot.title  = element_text(size = 13, hjust = 0.5)
    )
}

p_box_dis <- plot_box_phase(filter(score_df, Phase == "Discovery"), "Discovery")
p_box_hol <- plot_box_phase(filter(score_df, Phase == "Hold-out"), "Hold-out")
p_box_ind <- plot_box_phase(filter(score_df, Phase == "Independent"), paste0("Independent (", val_name, ")"))

ggsave("plots/Tscore_boxplot_discovery.png", p_box_dis, width = 10.5, height = 5.5, dpi = 300)
ggsave("plots/Tscore_boxplot_holdout.png",   p_box_hol, width = 10.5, height = 5.5, dpi = 300)
ggsave("plots/Tscore_boxplot_independent.png", p_box_ind, width = 10.5, height = 5.5, dpi = 300)

print(p_box_dis)
print(p_box_hol)
print(p_box_ind)

cat("Saved to ./plots/\n")


















# ==============================================================================
# F) Confusion matrices (Training / Hold-out / Independent)
# ==============================================================================

dir.create("plots", showWarnings = FALSE)
thr <- 0.5  # diagnosis threshold on P(Tumor)

plot_cm <- function(truth, prob_tumor, cohort_name, thr = 0.5, file_prefix = NULL) {
  # truth: factor with levels c("Tumor","Normal")
  stopifnot(all(levels(truth) == c("Tumor","Normal")))
  
  pred <- factor(ifelse(prob_tumor >= thr, "Tumor", "Normal"), levels = c("Tumor","Normal"))
  
  cm <- caret::confusionMatrix(data = pred, reference = truth, positive = "Tumor")
  print(cohort_name)
  print(cm)
  
  # build plotting df: Reference (rows) x Prediction (cols)
  tab <- as.data.frame(cm$table)
  colnames(tab) <- c("Prediction", "Reference", "N")  # caret table dims: data (pred) vs reference
  tab <- tab %>%
    group_by(Reference) %>%
    mutate(RowPct = 100 * N / sum(N)) %>%
    ungroup()
  
  p <- ggplot(tab, aes(x = Prediction, y = Reference, fill = N)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(N, "\n(", sprintf("%.1f", RowPct), "%)")), size = 5) +
    scale_y_discrete(limits = rev(levels(truth))) +
    labs(
      title = paste0("Confusion Matrix - ", cohort_name, "\nThreshold: P(Tumor) ≥ ", thr),
      x = "Predicted", y = "True (Reference)", fill = "Count"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 13, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 13, hjust = 0.5)
    )
  
  if (!is.null(file_prefix)) {
    ggsave(paste0("plots/", file_prefix, ".png"), p, width = 6.2, height = 5.2, dpi = 300)
    write.csv(cm$table, paste0("plots/", file_prefix, "_table.csv"), row.names = TRUE)
  }
  
  return(list(cm = cm, plot = p))
}

# ---- Training: use CV out-of-fold predictions (no leakage from in-sample fit) ----
# pred_cv assumed already filtered to bestTune, and has columns: obs, Tumor
truth_tr <- factor(pred_cv$obs, levels = c("Tumor","Normal"))
prob_tr  <- pred_cv$Tumor
res_cm_tr <- plot_cm(truth_tr, prob_tr, cohort_name = "Discovery (Training, CV out-of-fold)",
                     thr = thr, file_prefix = "CM_training_CV")

# ---- Hold-out ----
truth_te <- factor(test_data$Group, levels = c("Tumor","Normal"))
# prob_te assumed already computed: predict(model, test_data[, top_feats], type="prob")$Tumor
res_cm_te <- plot_cm(truth_te, prob_te, cohort_name = "Hold-out (Testing)",
                     thr = thr, file_prefix = "CM_holdout")

# ---- Independent ----
truth_iv <- factor(valid_data$Group, levels = c("Tumor","Normal"))
# prob_iv assumed already computed: predict(model, valid_data[, top_feats], type="prob")$Tumor
res_cm_iv <- plot_cm(truth_iv, prob_iv, cohort_name = paste0("Independent (", val_name, ")"),
                     thr = thr, file_prefix = "CM_independent")

# Print plots to device (optional)
print(res_cm_tr$plot)
print(res_cm_te$plot)
print(res_cm_iv$plot)

cat(">>> Confusion matrices saved to ./plots/ (PNG + raw table CSV)\n")








# 打印选中的 piRNA 名字
print(top_feats)

# 或者把它们保存到一个 CSV 文件里，方便写文章用
write.csv(top_feats, file = "Final_5_piRNA_Signature.csv", row.names = FALSE)

































































library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(ggplot2)
library(RColorBrewer) # 用于配色

# ==============================================================================
# 1. 数据准备 & 合并指定数据集
# ==============================================================================
# 读取所有文件
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))

# Log2 处理
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# --- 仅筛选 BRCA, yyfbatch1, yyfbatch2 ---
target_sets <- c("BRCA1", "yyfbatch1", "yyfbatch2") # 确保名字和你文件名一致
if("BRCA1" %in% names(clean_list)) {
  # 稍微平衡一下 BRCA1 (防止它数量太大淹没其他两个)
  set.seed(123)
  df <- clean_list[["BRCA1"]]
  idx_norm <- which(df$Group == "Normal")
  idx_tum <- sample(which(df$Group == "Tumor"), length(idx_norm) * 2) # 取 Normal 的2倍
  clean_list[["BRCA1"]] <- df[c(idx_norm, idx_tum), ]
}

selected_list <- clean_list[names(clean_list) %in% target_sets]

# --- 合并 & 去除批次效应 (Harmonization) ---
cat(">>> 正在合并数据集: ", paste(names(selected_list), collapse=", "), "...\n")

expr_matrices <- lapply(selected_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(selected_list), function(n) rep(n, nrow(selected_list[[n]]))))
group_vec <- unlist(lapply(selected_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组为大宽表
merged_df <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                        t(combat_expr), check.names = FALSE)
merged_df$Batch <- batch_vec
merged_df$SampleID <- rownames(merged_df) # 确保有 SampleID 用于匹配临床信息

# Z-score 标准化
gene_cols <- colnames(merged_df)[!colnames(merged_df) %in% c("Group", "Batch", "SampleID")]
merged_df[, gene_cols] <- scale(merged_df[, gene_cols])

# ==============================================================================
# 2. 准备/加载 临床信息 (关键步骤)
# ==============================================================================

# *** 请注意：这里我生成了模拟数据 ***
# *** 你需要读取你真实的 clinical.csv 文件并合并进来 ***

set.seed(42)
n_samples <- nrow(merged_df)
# 模拟临床数据
mock_clinical <- data.frame(
  SampleID = merged_df$SampleID,
  Age = sample(30:80, n_samples, replace = T),
  # 简化分期: Stage I, II, III, IV
  Stage = sample(c("Stage I", "Stage II", "Stage III", "Stage IV"), n_samples, replace = T, prob=c(0.3, 0.3, 0.2, 0.2)),
  # 简化分型
  Subtype = sample(c("Luminal A", "Luminal B", "HER2+", "Basal-like"), n_samples, replace = T)
)

# 如果是 Normal 样本，通常没有分期/分型，设为 NA
is_normal <- merged_df$Group == "Normal"
mock_clinical$Stage[is_normal] <- NA
mock_clinical$Subtype[is_normal] <- NA

# 合并到主数据表
final_df <- merge(merged_df, mock_clinical, by = "SampleID")

# 处理年龄：转换为二分类 (Age < 60 vs Age >= 60)
final_df$Age_Group <- ifelse(final_df$Age >= 60, "Age >= 60", "Age < 60")

# ==============================================================================
# 3. 训练模型并获取 T-Score (使用 5-fold CV)
# ==============================================================================
# 为了在整个合并队列上画图，我们使用交叉验证预测值，避免过拟合
cat(">>> 正在计算全队列预测分数 (5-fold CV)...\n")

# 特征筛选 (在全集上做一次快速筛选，或者使用你确定的5个)
set.seed(123)
rf_fs <- randomForest(x = final_df[, gene_cols], y = final_df$Group, ntree=100)
rf_imp <- importance(rf_fs)
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:5]
cat("Selected Features:", paste(top_feats, collapse=", "), "\n")

# 定义 CV
folds <- createFolds(final_df$Group, k = 5, list = TRUE)
final_df$T_Score <- NA # 初始化分数列

for(i in 1:5) {
  # 训练集索引
  idx_train <- unlist(folds[-i])
  idx_test  <- unlist(folds[i])
  
  train_data <- final_df[idx_train, c("Group", top_feats)]
  test_data  <- final_df[idx_test, c("Group", top_feats)]
  
  set.seed(123)
  model <- randomForest(Group ~ ., data = train_data, ntree = 300)
  
  # 预测测试集概率
  preds <- predict(model, test_data, type = "prob")[, "Tumor"]
  final_df$T_Score[idx_test] <- preds
}

# ==============================================================================
# 4. 亚组绘图通用函数
# ==============================================================================
plot_subgroup_roc <- function(df, group_col, title_text) {
  
  # 1. 提取非 NA 的子集 (排除 Normal 中没有分期/分型的情况，但也需要保留 Normal 作为对照!)
  # 策略：亚组分析通常是 "特定亚组肿瘤 vs 所有正常样本" 或者 "特定亚组肿瘤 vs 匹配的正常样本"
  # 这里采用：特定亚组肿瘤 (Group=Tumor & Subtype=X) vs 所有正常样本 (Group=Normal)
  
  df_normal <- df[df$Group == "Normal", ]
  df_tumor  <- df[df$Group == "Tumor", ]
  
  # 获取该临床特征下的所有类别 (去除 NA)
  subgroups <- unique(na.omit(df_tumor[[group_col]]))
  subgroups <- sort(subgroups)
  
  # 准备绘图数据
  plot_data <- data.frame()
  color_map <- c()
  auc_labels <- c()
  
  for(grp in subgroups) {
    # 选出该亚组的肿瘤
    sub_tumor <- df_tumor[df_tumor[[group_col]] == grp & !is.na(df_tumor[[group_col]]), ]
    
    # 如果该亚组样本太少(<5)，跳过
    if(nrow(sub_tumor) < 5) next
    
    # 合并：该亚组肿瘤 + 所有正常样本
    sub_df <- rbind(sub_tumor, df_normal)
    
    # 计算 ROC
    roc_obj <- roc(ifelse(sub_df$Group=="Tumor", 1, 0), sub_df$T_Score, direction="<", quiet=T)
    auc_val <- round(as.numeric(auc(roc_obj)), 3)
    
    # 提取坐标
    tmp_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      Subgroup = grp
    )
    plot_data <- rbind(plot_data, tmp_df)
    
    # 保存 AUC 标签
    auc_labels <- c(auc_labels, paste0(grp, " (AUC = ", auc_val, ")"))
  }
  
  # 绘图
  ggplot(plot_data, aes(x = fpr, y = tpr, color = Subgroup)) +
    geom_path(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    scale_color_brewer(palette = "Set1", labels = auc_labels) + # 自动配色
    labs(title = title_text,
         x = "1 - Specificity",
         y = "Sensitivity",
         color = "Subgroup AUC") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = c(0.7, 0.25), # 图例放在右下角内部
      legend.background = element_rect(fill = "white", color = "black"),
      panel.grid = element_blank()
    )
}

# ==============================================================================
# 5. 执行绘图
# ==============================================================================

# A. 按年龄分组 (Age Group)
# 逻辑：Age < 60 Tumor vs Normal; Age >= 60 Tumor vs Normal
# 注意：这里需要把 Normal 也按年龄分吗？
# 通常做法是：Tumor按年龄分，Normal使用全部Normal作为共有背景，或者Normal也按年龄匹配。
# 为了简化且更有鲁棒性，这里比较：特定年龄段的肿瘤 vs 该年龄段的正常人
# 重新定义函数逻辑有点复杂，我们简化为：Tumor (Subgroup) vs All Normal (Control)
p1 <- plot_subgroup_roc(final_df, "Age_Group", "ROC Stratified by Age")

# B. 按分期分组 (Stage)
p2 <- plot_subgroup_roc(final_df, "Stage", "ROC Stratified by Cancer Stage")

# C. 按分型分组 (Subtype)
p3 <- plot_subgroup_roc(final_df, "Subtype", "ROC Stratified by Subtype")

# 打印
print(p1)
print(p2)
print(p3)

# 保存图片
ggsave("ROC_Subgroup_Age.pdf", p1, width = 6, height = 6)
ggsave("ROC_Subgroup_Stage.pdf", p2, width = 6, height = 6)
ggsave("ROC_Subgroup_Subtype.pdf", p3, width = 6, height = 6)




library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# ==============================================================================
# 0. 准备工作 (确保 final_df 存在)
# ==============================================================================
# 假设你已经运行了上一步的代码，内存中已经有了 'final_df'
# final_df 包含: Group, T_Score, Age_Group, Stage, Subtype 等列

# 如果没有 final_df，请先运行上一段代码生成它。
# 这里做一个简单的检查
if(!exists("final_df")) stop("请先运行上一段代码以生成 final_df！")

# ==============================================================================
# 1. 定义通用绘图函数 (自动匹配 Normal 对照)
# ==============================================================================
plot_subgroup_boxplot <- function(df, subgroup_col, title_text) {
  
  # --- A. 数据重构 (关键步骤) ---
  # 提取 Tumor 和 Normal
  df_tumor  <- df[df$Group == "Tumor", ]
  df_normal <- df[df$Group == "Normal", ]
  
  # 获取所有有效的亚组类别 (去除 NA)
  valid_subgroups <- sort(unique(na.omit(df_tumor[[subgroup_col]])))
  
  # 构建绘图数据：
  # 对于每个亚组，我们取该亚组的 Tumor + *所有的* Normal
  plot_data <- data.frame()
  
  for(grp in valid_subgroups) {
    # 1. 取出该亚组的 Tumor
    sub_tumor <- df_tumor[df_tumor[[subgroup_col]] == grp, ]
    
    # 如果样本太少，跳过
    if(nrow(sub_tumor) < 3) next
    
    # 2. 取出所有 Normal (赋予当前的亚组标签，以便画在同一栏)
    # 注意：这里我们给 Normal 强行打上 grp 标签，只是为了让它在图上显示在该 grp 的旁边
    tmp_normal <- df_normal
    tmp_normal[[subgroup_col]] <- grp 
    
    # 3. 合并
    combined <- rbind(sub_tumor, tmp_normal)
    plot_data <- rbind(plot_data, combined)
  }
  
  # 确保亚组列是因子 (保持顺序)
  plot_data[[subgroup_col]] <- factor(plot_data[[subgroup_col]], levels = valid_subgroups)
  
  # --- B. 绘图 ---
  # 定义颜色: Normal=灰色, Tumor=红色
  my_colors <- c("#868686FF", "#CD534CFF")
  
  p <- ggplot(plot_data, aes_string(x = subgroup_col, y = "T_Score", fill = "Group")) +
    stat_boxplot(geom ='errorbar', width = 0.2, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, position = position_dodge(0.8)) +
    # 添加抖动散点
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                size = 0.5, alpha = 0.3, color="black") +
    # 添加统计检验 (Tumor vs Normal in each subgroup)
    stat_compare_means(aes(group = Group), method = "wilcox.test", 
                       label = "p.signif", label.y = 1.05, size = 5) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
    labs(title = title_text, x = "", y = "Predicted T-Score") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 15, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid = element_blank()
    )
  
  return(p)
}

# ==============================================================================
# 2. 执行绘图
# ==============================================================================

# A. 年龄分组 (Age Group)
p_box_age <- plot_subgroup_boxplot(final_df, "Age_Group", "T-Score by Age Group")

# B. 分期分组 (Stage)
# 确保 Stage 顺序正确
valid_stages <- c("Stage I", "Stage II", "Stage III", "Stage IV")
# 如果数据里只有部分分期，取交集
available_stages <- intersect(valid_stages, unique(final_df$Stage))
final_df$Stage <- factor(final_df$Stage, levels = available_stages)

p_box_stage <- plot_subgroup_boxplot(final_df, "Stage", "T-Score by Cancer Stage")

# C. 分型分组 (Subtype)
p_box_subtype <- plot_subgroup_boxplot(final_df, "Subtype", "T-Score by Subtype")

# ==============================================================================
# 3. 输出与保存
# ==============================================================================
print(p_box_age)
print(p_box_stage)
print(p_box_subtype)

# 保存
ggsave("Boxplot_Age.pdf", p_box_age, width = 6, height = 5)
ggsave("Boxplot_Stage.pdf", p_box_stage, width = 8, height = 5)
ggsave("Boxplot_Subtype.pdf", p_box_subtype, width = 8, height = 5)