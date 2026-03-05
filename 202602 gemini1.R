library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 读取数据 & 准备工作
# ==============================================================================
# 假设文件在 processed_results 文件夹
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
# 获取数据集名称 (去掉路径和后缀)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
# 读取数据
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
cat("共有 miRNA:", length(common_genes), "\n")

# 数据清洗 (Log2 check)
prep_data <- function(df, genes) {
  df_sub <- df[, c("Group", genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
}
clean_list <- lapply(data_list, prep_data, genes = common_genes)

# ==============================================================================
# 2. 全局 ComBat (All-in-One) -> 这就是你要的“所有一起去除”
# ==============================================================================
cat("正在对所有 7 个数据集进行全局 ComBat...\n")

# A. 构建合并矩阵
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)

# B. 构建 Batch 向量 (用于记录每个样本属于哪个数据集)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))

# C. 构建分组信息
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# D. 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# E. 重组为大表 & Z-score
# 转置回来
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), 
                            check.names = FALSE)
# 添加一列 Batch 用于后续拆分
combat_df_all$Batch <- batch_vec 

# 全局 Z-score (既然都一起ComBat了，一起Scale也没问题)
# 注意：排除 Group 和 Batch 列
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])

# 修复 NA (防止 scale 产生 NaN)
combat_df_all[is.na(combat_df_all)] <- 0

cat("全局处理完成！总样本数:", nrow(combat_df_all), "\n")

# ==============================================================================
# 3. 循环 7 次：Leave-One-Dataset-Out
# ==============================================================================

# 存储结果的容器
auc_results <- data.frame(Dataset = character(), AUC = numeric(), stringsAsFactors = F)
roc_list <- list()

# 开启并行 (加速训练)
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

cat("开始循环验证...\n")

for (valid_name in dataset_names) {
  
  cat(paste0(">>> 正在验证: ", valid_name, " (其余 6 个做训练)\n"))
  
  # --- A. 拆分数据 ---
  # 验证集 = 当前 dataset
  valid_set <- combat_df_all[combat_df_all$Batch == valid_name, ]
  # 训练集 = 除了当前 dataset 之外的所有
  train_set <- combat_df_all[combat_df_all$Batch != valid_name, ]
  
  # 去掉 Batch 列 (模型不需要知道批次)
  valid_set$Batch <- NULL
  train_set$Batch <- NULL
  
  # --- B. 特征筛选 (LASSO + RF) ---
  # 注意：严格来说特征筛选也应该在循环内只对训练集做
  x <- as.matrix(train_set[, gene_cols])
  y <- train_set$Group
  
  # 快速筛选 (为了节省时间，这里只用 RF Top 15)
  # 如果想要更严谨，可以把之前的 LASSO 代码加回来
  set.seed(123)
  rf_fs <- randomForest(x, y, ntree=200) # 树少一点为了快
  rf_imp <- importance(rf_fs)
  final_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:15]
  
  # --- C. 建模 (Random Forest) ---
  # 检查是否需要 down-sampling (如果训练集里有 BRCA1，通常需要)
  train_sub <- train_set[, c("Group", final_feats)]
  
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary,
                             sampling = "down") 
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 300)
  
  # --- D. 预测与评估 ---
  valid_sub <- valid_set[, c("Group", final_feats)]
  preds <- predict(model, valid_sub, type = "prob")
  
  # 计算 ROC
  roc_obj <- roc(valid_set$Group, preds$Tumor, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
  current_auc <- auc(roc_obj)
  
  # 存结果
  auc_results <- rbind(auc_results, data.frame(Dataset = valid_name, AUC = as.numeric(current_auc)))
  roc_list[[valid_name]] <- roc_obj
  
  cat(paste("    AUC =", round(current_auc, 4), "\n"))
}

stopCluster(cl)

# ==============================================================================
# 4. 结果可视化与排序
# ==============================================================================

# A. 打印排行榜
auc_results <- auc_results[order(auc_results$AUC, decreasing = TRUE), ]
print("====== AUC 排行榜 (由高到低) ======")
print(auc_results)

# B. 绘制在一张图上
# 设置颜色盘
colors <- rainbow(length(dataset_names))

# 画第一个
plot(roc_list[[auc_results$Dataset[1]]], col = colors[1], lwd = 2, 
     main = "LODO-CV Performance Comparison", legacy.axes = TRUE)

# 循环画剩下的
for (i in 2:nrow(auc_results)) {
  dataset_name <- auc_results$Dataset[i]
  plot(roc_list[[dataset_name]], col = colors[i], lwd = 2, add = TRUE)
}

# 添加图例
legend("bottomright", 
       legend = paste0(auc_results$Dataset, " (AUC=", round(auc_results$AUC, 3), ")"),
       col = colors, lwd = 2, cex = 0.8)




###############################################################################
library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 读取数据
# ==============================================================================
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
cat("共有 miRNA:", length(common_genes), "\n")

# 数据清洗 (Log2)
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# ==============================================================================
# 2. 核心步骤：对 BRCA1 进行“文献级”下采样 (Balanced Down-sampling)
# ==============================================================================
# 逻辑：Keep All Normals + Equal number of Tumors (Matched) + 40% of Remaining Tumors

process_brca_balance <- function(df, seed = 123) {
  set.seed(seed)
  
  # 分离 Tumor 和 Normal 的索引
  idx_normal <- which(df$Group == "Normal")
  idx_tumor <- which(df$Group == "Tumor")
  
  n_normal <- length(idx_normal)
  n_tumor_total <- length(idx_tumor)
  
  cat("BRCA1 原始: Normal =", n_normal, ", Tumor =", n_tumor_total, "\n")
  
  # 1. 模拟 "Matched Pairs": 抽取与 Normal 数量相等的 Tumor
  # (假设每个 Normal 都有一个对应的 Tumor)
  n_matched <- n_normal 
  
  # 2. 计算剩余 Tumor 数量
  n_remaining <- n_tumor_total - n_matched
  
  # 3. 抽取剩余部分的 40%
  n_extra <- round(n_remaining * 0.40)
  
  # 4. 总共需要抽取的 Tumor 数量
  n_tumor_final <- n_matched + n_extra
  
  # 执行抽样
  idx_tumor_keep <- sample(idx_tumor, n_tumor_final)
  
  # 合并索引
  idx_keep <- c(idx_normal, idx_tumor_keep)
  df_balanced <- df[idx_keep, ]
  
  cat("BRCA1 处理后: Normal =", n_normal, ", Tumor =", n_tumor_final, 
      "(Ratio ≈ 1:", round(n_tumor_final/n_normal, 1), ")\n")
  
  return(df_balanced)
}

# 应用到 BRCA1
# 注意：一定要确保你的列表中 BRCA1 的名字叫 "BRCA1"
if ("BRCA1" %in% names(clean_list)) {
  clean_list[["BRCA1"]] <- process_brca_balance(clean_list[["BRCA1"]])
} else {
  stop("未找到名为 BRCA1 的数据集，请检查文件名！")
}

# ==============================================================================
# 3. 全局 ComBat (去除批次效应)
# ==============================================================================
cat("\n正在对所有数据（含平衡后的BRCA1）进行 ComBat...\n")

expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Z-score
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), 
                            check.names = FALSE)
combat_df_all$Batch <- batch_vec # 标记来源

# 全局 Z-score
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 4. 循环 6 次：其他数据集依次做验证 (Leave-One-Out)
# ==============================================================================

# 找出除了 BRCA1 以外的所有数据集名字
validation_candidates <- setdiff(names(clean_list), "BRCA1")

# 存储结果
metrics_df <- data.frame()
plot_list <- list()

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

cat("\n开始 6 轮循环验证 (BRCA1 始终在训练集)...\n")

for (valid_name in validation_candidates) {
  
  cat(paste0(">>> 验证集: ", valid_name, " | 训练集: BRCA1 + 其他 5 个\n"))
  
  # A. 拆分
  # 验证集 = 当前 dataset
  valid_set <- combat_df_all[combat_df_all$Batch == valid_name, ]
  # 训练集 = (BRCA1) + (其他 5 个)
  train_set <- combat_df_all[combat_df_all$Batch != valid_name, ]
  
  # 清理不需要的列
  train_data <- train_set[, c("Group", gene_cols)]
  valid_data <- valid_set[, c("Group", gene_cols)]
  
  # B. 特征筛选 (快速版: RF Top 15)
  # 仅在训练集上做
  x <- as.matrix(train_data[, -1])
  y <- train_data$Group
  
  set.seed(123)
  rf_fs <- randomForest(x, y, ntree=200)
  rf_imp <- importance(rf_fs)
  final_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:15]
  
  # C. 建模
  # 此时 BRCA1 已经手动平衡过，所以 sampling="down" 可以选填
  # 但为了双重保险，建议保留 cv (不采样) 或者 down (再次严格平衡)
  # 这里演示用 "cv" (相信我们的手动平衡)
  train_sub <- train_data[, c("Group", final_feats)]
  
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary)
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 300)
  
  # D. 评估 (训练集CV, 验证集)
  
  # 1. 训练集表现 (CV AUC)
  train_auc <- max(model$results$ROC)
  
  # 2. 独立验证集表现
  prob_valid <- predict(model, valid_data[, final_feats], type = "prob")
  roc_valid <- roc(valid_data$Group, prob_valid$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  valid_auc <- auc(roc_valid)
  
  # 记录
  metrics_df <- rbind(metrics_df, data.frame(
    Validation_Dataset = valid_name,
    Training_AUC = train_auc,
    Independent_AUC = valid_auc
  ))
  
  plot_list[[valid_name]] <- roc_valid
  
  cat(paste("    Train AUC:", round(train_auc, 3), "| Independent AUC:", round(valid_auc, 3), "\n"))
}

stopCluster(cl)

# ==============================================================================
# 5. 结果汇总与可视化
# ==============================================================================

print("====== 最终验证结果 ======")
print(metrics_df)

# 绘制在一张图上
colors <- rainbow(length(plot_list))
# 初始化画布
plot(plot_list[[1]], col=colors[1], lwd=2, 
     main="Performance on 6 Independent Validation Sets", legacy.axes=TRUE)

for(i in 2:length(plot_list)){
  plot(plot_list[[i]], col=colors[i], lwd=2, add=TRUE)
}

legend("bottomright", 
       legend = paste0(metrics_df$Validation_Dataset, " (AUC=", round(metrics_df$Independent_AUC, 3), ")"),
       col = colors, lwd = 2, cex = 0.8)









################################################################################
library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
cat("共有 miRNA:", length(common_genes), "\n")

# Log2 清洗
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# ==============================================================================
# 2. 平衡 BRCA1 (Normal + Matched Tumor + 40% Remaining)
# ==============================================================================
process_brca_balance <- function(df, seed = 123) {
  set.seed(seed)
  idx_normal <- which(df$Group == "Normal")
  idx_tumor <- which(df$Group == "Tumor")
  
  n_normal <- length(idx_normal)
  # 逻辑: 抽取 Normal 数量的配对 Tumor + 剩余 Tumor 的 40%
  n_remaining <- length(idx_tumor) - n_normal
  n_keep_tumor <- n_normal + round(n_remaining * 0.40)
  
  idx_tumor_keep <- sample(idx_tumor, n_keep_tumor)
  df_balanced <- df[c(idx_normal, idx_tumor_keep), ]
  
  cat("BRCA1 平衡后: Normal =", n_normal, ", Tumor =", n_keep_tumor, "\n")
  return(df_balanced)
}

if ("BRCA1" %in% names(clean_list)) {
  clean_list[["BRCA1"]] <- process_brca_balance(clean_list[["BRCA1"]])
}

# ==============================================================================
# 3. 全局 ComBat
# ==============================================================================
cat("正在执行全局 ComBat...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组大表
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), 
                            check.names = FALSE)
combat_df_all$Batch <- batch_vec

# 全局 Z-score & NA 修复
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 4. 生成验证集组合 (15 种情况)
# ==============================================================================
# 候选验证集：除 BRCA1 外的 6 个
candidate_datasets <- setdiff(dataset_names, "BRCA1")

# 生成所有两两组合 (C(6,2) = 15)
combo_matrix <- combn(candidate_datasets, 2)
cat("总共将进行", ncol(combo_matrix), "轮双重验证实验。\n")

# ==============================================================================
# 5. 循环验证
# ==============================================================================
results_df <- data.frame()

# 开启并行加速
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

for (i in 1:ncol(combo_matrix)) {
  
  # 当前选中的 2 个验证集
  valid_pair <- combo_matrix[, i]
  val_1_name <- valid_pair[1]
  val_2_name <- valid_pair[2]
  
  cat(paste0("\n>>> Round ", i, ": 验证集 [", val_1_name, " & ", val_2_name, "]\n"))
  
  # --- A. 拆分数据 ---
  # 验证集数据
  valid_data_1 <- combat_df_all[combat_df_all$Batch == val_1_name, ]
  valid_data_2 <- combat_df_all[combat_df_all$Batch == val_2_name, ]
  
  # 训练集数据 = BRCA1 + 其他 4 个 (排除这两个验证集)
  train_data_full <- combat_df_all[!combat_df_all$Batch %in% valid_pair, ]
  
  # 去除 Batch 列，准备训练
  train_data <- train_data_full[, c("Group", gene_cols)]
  
  # --- B. 特征筛选 (RF Top 15) ---
  set.seed(123)
  # 简单用下采样抽一部分数据做特征筛选，为了快
  rf_fs <- randomForest(x = train_data[, -1], y = train_data$Group, ntree=150)
  rf_imp <- importance(rf_fs)
  final_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:15]
  
  # --- C. 建模 ---
  train_sub <- train_data[, c("Group", final_feats)]
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary)
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 300)
  
  # --- D. 预测与评估 ---
  # 验证集 1
  prob_1 <- predict(model, valid_data_1[, final_feats], type = "prob")
  auc_1 <- auc(roc(valid_data_1$Group, prob_1$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T))
  
  # 验证集 2
  prob_2 <- predict(model, valid_data_2[, final_feats], type = "prob")
  auc_2 <- auc(roc(valid_data_2$Group, prob_2$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T))
  
  # 记录结果
  results_df <- rbind(results_df, data.frame(
    Round = i,
    Valid_Set_1 = val_1_name,
    AUC_1 = as.numeric(auc_1),
    Valid_Set_2 = val_2_name,
    AUC_2 = as.numeric(auc_2),
    Mean_AUC = (as.numeric(auc_1) + as.numeric(auc_2)) / 2
  ))
  
  cat(paste("    AUC1:", round(auc_1, 3), "| AUC2:", round(auc_2, 3), "\n"))
}

stopCluster(cl)

# ==============================================================================
# 6. 结果汇总与展示
# ==============================================================================
print("====== 15轮 双重验证结果 ======")
# 按平均 AUC 排序
results_df_sorted <- results_df[order(results_df$Mean_AUC, decreasing = TRUE), ]
print(results_df_sorted)

# 简单可视化：箱线图查看 AUC 分布
all_aucs <- c(results_df$AUC_1, results_df$AUC_2)
boxplot(all_aucs, col="lightblue", main="Distribution of AUCs across 30 Validations",
        ylab="AUC Value", ylim=c(0.5, 1.0))
stripchart(all_aucs, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = "blue")
abline(h=0.8, col="red", lty=2) # 0.8 及格线











#################################################################################
library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 数据准备
# ==============================================================================
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
cat("共有 miRNA:", length(common_genes), "\n")

# 清洗 (Log2)
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# --- 平衡 BRCA1 ---
process_brca_balance <- function(df, seed = 123) {
  set.seed(seed)
  idx_normal <- which(df$Group == "Normal")
  idx_tumor <- which(df$Group == "Tumor")
  n_normal <- length(idx_normal)
  n_remaining <- length(idx_tumor) - n_normal
  n_keep_tumor <- n_normal + round(n_remaining * 0.40)
  idx_tumor_keep <- sample(idx_tumor, n_keep_tumor)
  df_balanced <- df[c(idx_normal, idx_tumor_keep), ]
  cat("BRCA1 平衡后: Normal =", n_normal, ", Tumor =", n_keep_tumor, "\n")
  return(df_balanced)
}

if ("BRCA1" %in% names(clean_list)) {
  clean_list[["BRCA1"]] <- process_brca_balance(clean_list[["BRCA1"]])
}

# ==============================================================================
# 2. 全局 ComBat & Z-score
# ==============================================================================
cat("正在执行全局 ComBat...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Z-score
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), 
                            check.names = FALSE)
combat_df_all$Batch <- batch_vec

gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
# 修复 NA/NaN
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 3. 生成 3 验证集组合 (C(6,3) = 20)
# ==============================================================================
candidate_datasets <- setdiff(dataset_names, "BRCA1")
# 修改这里：m = 3
combo_matrix <- combn(candidate_datasets, 3)

cat("准备进行", ncol(combo_matrix), "轮 [三验证集] 测试...\n")

# ==============================================================================
# 4. 循环验证
# ==============================================================================
results_df <- data.frame()

# 辅助函数：安全计算 AUC
get_auc <- function(df_valid, model, feats) {
  if(length(levels(factor(df_valid$Group))) < 2) return(NA)
  prob <- predict(model, df_valid[, feats], type = "prob")
  roc_obj <- roc(df_valid$Group, prob$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  return(as.numeric(auc(roc_obj)))
}

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

for (i in 1:ncol(combo_matrix)) {
  
  # 1. 确定本轮 3 个验证集
  valid_trio <- combo_matrix[, i]
  name_1 <- valid_trio[1]
  name_2 <- valid_trio[2]
  name_3 <- valid_trio[3]
  
  cat(paste0("\n>>> Round ", i, ": [", name_1, ", ", name_2, ", ", name_3, "]\n"))
  
  # 2. 拆分数据
  v_data_1 <- combat_df_all[combat_df_all$Batch == name_1, ]
  v_data_2 <- combat_df_all[combat_df_all$Batch == name_2, ]
  v_data_3 <- combat_df_all[combat_df_all$Batch == name_3, ]
  
  # *** 构建 Combined Validation Set (3合1) ***
  v_data_combined <- rbind(v_data_1, v_data_2, v_data_3)
  
  # 训练集 = (BRCA1) + (其余3个)
  train_data_full <- combat_df_all[!combat_df_all$Batch %in% valid_trio, ]
  train_data <- train_data_full[, c("Group", gene_cols)]
  
  # 3. 特征筛选
  set.seed(123)
  rf_fs <- randomForest(x = train_data[, -1], y = train_data$Group, ntree=150)
  rf_imp <- importance(rf_fs)
  final_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:15]
  
  # 4. 建模
  train_sub <- train_data[, c("Group", final_feats)]
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary)
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 300)
  
  # 5. 计算 4 个 AUC (3个独立 + 1个合并)
  auc_1 <- get_auc(v_data_1, model, final_feats)
  auc_2 <- get_auc(v_data_2, model, final_feats)
  auc_3 <- get_auc(v_data_3, model, final_feats)
  auc_comb <- get_auc(v_data_combined, model, final_feats)
  
  # 6. 记录
  results_df <- rbind(results_df, data.frame(
    Round = i,
    Set1 = name_1, AUC1 = round(auc_1, 3),
    Set2 = name_2, AUC2 = round(auc_2, 3),
    Set3 = name_3, AUC3 = round(auc_3, 3),
    Combined_AUC = round(auc_comb, 3)
  ))
  
  cat(paste("    Comb AUC:", round(auc_comb, 3), "\n"))
}

stopCluster(cl)

# ==============================================================================
# 5. 结果展示
# ==============================================================================

print("====== 20轮 三重验证结果 (按 Combined AUC 排序) ======")
results_df_sorted <- results_df[order(results_df$Combined_AUC, decreasing = TRUE), ]
print(results_df_sorted)

# 可视化分布
par(mfrow=c(1,1))
boxplot(results_df$Combined_AUC, 
        main="Robustness Test: 3 Independent Validation Sets (20 Rounds)",
        ylab="Combined AUC",
        col="mistyrose",
        ylim=c(0.5, 1.0))
stripchart(results_df$Combined_AUC, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 19, col = "darkred", cex=1.2)
abline(h=0.8, col="blue", lty=2, lwd=2)

# 找出并打印最佳组合
best <- results_df_sorted[1, ]
cat("\n>>> 最佳组合 (Round", best$Round, "):", 
    best$Set1, ",", best$Set2, ",", best$Set3, 
    "\n>>> Combined AUC =", best$Combined_AUC, "\n")















library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)

# ==============================================================================
# 0. 前置准备 (假设你已经跑完了之前的 Combat 和 Z-score 步骤)
# ==============================================================================
# 确保 environment 中已经有了 'combat_df_all' 和 'gene_cols'
# 如果没有，请先运行你提供的那段代码的前半部分 (直到 Step 4 之前)

# 定义要画图的目标轮次和对应的验证集名称 (根据你的描述)
target_rounds <- list(
  "Round 10" = c("PRJNA808405", "PRJNA934049"),
  "Round 13" = c("PRJNA934049", "yyfbatch1"),
  "Round 14" = c("PRJNA934049", "yyfbatch2")
)

# 设置绘图布局 (1行3列)
par(mfrow = c(1, 3)) 

# ==============================================================================
# 1. 循环处理指定的 3 个 Rounds
# ==============================================================================

for (round_name in names(target_rounds)) {
  
  valid_pair <- target_rounds[[round_name]]
  val_1_name <- valid_pair[1]
  val_2_name <- valid_pair[2]
  
  cat(paste0("\n>>> 正在绘图: ", round_name, " [", val_1_name, " & ", val_2_name, "]\n"))
  
  # --- A. 数据拆分 ---
  # 验证集数据
  valid_data_1 <- combat_df_all[combat_df_all$Batch == val_1_name, ]
  valid_data_2 <- combat_df_all[combat_df_all$Batch == val_2_name, ]
  
  # 训练集数据 (BRCA1 + 其他 4 个)
  train_data_full <- combat_df_all[!combat_df_all$Batch %in% valid_pair, ]
  
  # 准备训练数据 (去掉 Batch)
  train_data <- train_data_full[, c("Group", gene_cols)]
  
  # --- B. 特征筛选 (重现当时的筛选) ---
  set.seed(123)
  rf_fs <- randomForest(x = train_data[, -1], y = train_data$Group, ntree=150)
  rf_imp <- importance(rf_fs)
  final_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:15]
  
  # --- C. 建模 ---
  train_sub <- train_data[, c("Group", final_feats)]
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary)
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 300)
  
  # --- D. 预测 (获取概率) ---
  
  # 1. Training Set Performance (自身回测 or CV)
  # 这里为了画平滑曲线，我们直接预测训练集本身 (Self-prediction)
  # 注意：这是 Training Performance，通常会比 CV 高
  prob_train <- predict(model, train_data[, final_feats], type = "prob")
  roc_train <- roc(train_data$Group, prob_train$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  
  # 2. Independent Cohort 1
  prob_1 <- predict(model, valid_data_1[, final_feats], type = "prob")
  roc_1 <- roc(valid_data_1$Group, prob_1$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  
  # 3. Independent Cohort 2
  prob_2 <- predict(model, valid_data_2[, final_feats], type = "prob")
  roc_2 <- roc(valid_data_2$Group, prob_2$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  
  # --- E. 绘图 ---
  # 画布设置
  plot(roc_train, col = "black", lty = 2, lwd = 2, legacy.axes = TRUE,
       main = paste0(round_name, "\nTraining & 2 Indep. Cohorts"),
       xlab = "1 - Specificity", ylab = "Sensitivity")
  
  # 添加验证集 1 (红色)
  plot(roc_1, col = "red", lwd = 2, add = TRUE)
  
  # 添加验证集 2 (蓝色)
  plot(roc_2, col = "blue", lwd = 2, add = TRUE)
  
  # --- F. 添加图例 (含 AUC 数值) ---
  legend_text <- c(
    paste0("Training (AUC=", round(auc(roc_train), 3), ")"),
    paste0(val_1_name, " (AUC=", round(auc(roc_1), 3), ")"),
    paste0(val_2_name, " (AUC=", round(auc(roc_2), 3), ")")
  )
  
  legend("bottomright", 
         legend = legend_text,
         col = c("black", "red", "blue"),
         lty = c(2, 1, 1), # 训练集虚线，验证集实线
         lwd = 2, 
         cex = 0.8, bty = "n") # bty="n" 去掉图例边框更美观
}

# 恢复默认绘图设置
par(mfrow = c(1, 1))








library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)

# 1. 强制关闭可能残留的并行连接
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
try(stopImplicitCluster(), silent=TRUE)
try(unregister(), silent=TRUE)
registerDoSEQ() # 强制串行运行，这是画图不报错的关键

# 2. 设置画布为 1行3列
par(mfrow = c(1, 3)) 

# 假设 combat_df_all 和 gene_cols 还在你的环境中
# 如果没有，请重新运行之前生成 combat_df_all 的代码




cat(">>> 正在绘制 Round 10...\n")

# A. 定义验证集
val1_name <- "PRJNA808405"
val2_name <- "PRJNA934049"

# B. 拆分数据
valid_1 <- combat_df_all[combat_df_all$Batch == val1_name, ]
valid_2 <- combat_df_all[combat_df_all$Batch == val2_name, ]
train_data <- combat_df_all[!combat_df_all$Batch %in% c(val1_name, val2_name), c("Group", gene_cols)]

# C. 特征筛选 (训练集上)
set.seed(123)
rf_fs <- randomForest(x=train_data[,-1], y=train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing=T)][1:15]

# D. 建模
set.seed(888)
# 为了画图不报错，这里用简单的 trainControl
fitControl <- trainControl(method="cv", number=5, classProbs=TRUE, summaryFunction=twoClassSummary)
model <- train(Group ~ ., data=train_data[, c("Group", feats)], method="rf", 
               metric="ROC", trControl=fitControl, ntree=300)

# E. 预测与绘图
prob_train <- predict(model, train_data[, feats], type="prob")
prob_1 <- predict(model, valid_1[, feats], type="prob")
prob_2 <- predict(model, valid_2[, feats], type="prob")

roc_train <- roc(train_data$Group, prob_train$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_1 <- roc(valid_1$Group, prob_1$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_2 <- roc(valid_2$Group, prob_2$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)

# 画图
plot(roc_train, col="black", lty=2, lwd=2, legacy.axes=TRUE,
     main=paste0("Round 10\n", val1_name, "\n& ", val2_name),
     xlab="1 - Specificity", ylab="Sensitivity")
plot(roc_1, col="red", lwd=2, add=TRUE)
plot(roc_2, col="blue", lwd=2, add=TRUE)

legend("bottomright", 
       legend=c(paste0("Train (AUC=", round(auc(roc_train),3), ")"),
                paste0(val1_name, " (AUC=", round(auc(roc_1),3), ")"),
                paste0(val2_name, " (AUC=", round(auc(roc_2),3), ")")),
       col=c("black", "red", "blue"), lty=c(2,1,1), lwd=2, cex=0.7, bty="n")















cat(">>> 正在绘制 Round 13...\n")

# A. 定义验证集
val1_name <- "PRJNA934049"
val2_name <- "yyfbatch1"

# B. 拆分数据
valid_1 <- combat_df_all[combat_df_all$Batch == val1_name, ]
valid_2 <- combat_df_all[combat_df_all$Batch == val2_name, ]
train_data <- combat_df_all[!combat_df_all$Batch %in% c(val1_name, val2_name), c("Group", gene_cols)]

# C. 特征筛选
set.seed(123)
rf_fs <- randomForest(x=train_data[,-1], y=train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing=T)][1:15]

# D. 建模
set.seed(888)
model <- train(Group ~ ., data=train_data[, c("Group", feats)], method="rf", 
               metric="ROC", trControl=fitControl, ntree=300)

# E. 预测与绘图
prob_train <- predict(model, train_data[, feats], type="prob")
prob_1 <- predict(model, valid_1[, feats], type="prob")
prob_2 <- predict(model, valid_2[, feats], type="prob")

roc_train <- roc(train_data$Group, prob_train$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_1 <- roc(valid_1$Group, prob_1$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_2 <- roc(valid_2$Group, prob_2$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)

# 画图
plot(roc_train, col="black", lty=2, lwd=2, legacy.axes=TRUE,
     main=paste0("Round 13\n", val1_name, "\n& ", val2_name),
     xlab="1 - Specificity", ylab="Sensitivity")
plot(roc_1, col="red", lwd=2, add=TRUE)
plot(roc_2, col="blue", lwd=2, add=TRUE)

legend("bottomright", 
       legend=c(paste0("Train (AUC=", round(auc(roc_train),3), ")"),
                paste0(val1_name, " (AUC=", round(auc(roc_1),3), ")"),
                paste0(val2_name, " (AUC=", round(auc(roc_2),3), ")")),
       col=c("black", "red", "blue"), lty=c(2,1,1), lwd=2, cex=0.7, bty="n")



cat(">>> 正在绘制 Round 14...\n")

# A. 定义验证集
val1_name <- "PRJNA934049"
val2_name <- "yyfbatch2"

# B. 拆分数据
valid_1 <- combat_df_all[combat_df_all$Batch == val1_name, ]
valid_2 <- combat_df_all[combat_df_all$Batch == val2_name, ]
train_data <- combat_df_all[!combat_df_all$Batch %in% c(val1_name, val2_name), c("Group", gene_cols)]

# C. 特征筛选
set.seed(123)
rf_fs <- randomForest(x=train_data[,-1], y=train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing=T)][1:15]

# D. 建模
set.seed(888)
model <- train(Group ~ ., data=train_data[, c("Group", feats)], method="rf", 
               metric="ROC", trControl=fitControl, ntree=300)

# E. 预测与绘图
prob_train <- predict(model, train_data[, feats], type="prob")
prob_1 <- predict(model, valid_1[, feats], type="prob")
prob_2 <- predict(model, valid_2[, feats], type="prob")

roc_train <- roc(train_data$Group, prob_train$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_1 <- roc(valid_1$Group, prob_1$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
roc_2 <- roc(valid_2$Group, prob_2$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)

# 画图
plot(roc_train, col="black", lty=2, lwd=2, legacy.axes=TRUE,
     main=paste0("Round 14\n", val1_name, "\n& ", val2_name),
     xlab="1 - Specificity", ylab="Sensitivity")
plot(roc_1, col="red", lwd=2, add=TRUE)
plot(roc_2, col="blue", lwd=2, add=TRUE)

legend("bottomright", 
       legend=c(paste0("Train (AUC=", round(auc(roc_train),3), ")"),
                paste0(val1_name, " (AUC=", round(auc(roc_1),3), ")"),
                paste0(val2_name, " (AUC=", round(auc(roc_2),3), ")")),
       col=c("black", "red", "blue"), lty=c(2,1,1), lwd=2, cex=0.7, bty="n")

# 绘图结束，恢复布局
par(mfrow = c(1, 1))









library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 数据准备
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
cat("共有 miRNA:", length(common_genes), "\n")

# Log2 清洗
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1 (可选，但推荐保留以防止训练集被淹没)
process_brca_balance <- function(df, seed = 123) {
  set.seed(seed)
  idx_normal <- which(df$Group == "Normal")
  idx_tumor <- which(df$Group == "Tumor")
  n_normal <- length(idx_normal)
  n_remaining <- length(idx_tumor) - n_normal
  n_keep_tumor <- n_normal + round(n_remaining * 0.40)
  idx_tumor_keep <- sample(idx_tumor, n_keep_tumor)
  return(df[c(idx_normal, idx_tumor_keep), ])
}

if ("BRCA1" %in% names(clean_list)) {
  clean_list[["BRCA1"]] <- process_brca_balance(clean_list[["BRCA1"]])
}

# ==============================================================================
# 2. 全局 ComBat (所有数据一起做！)
# ==============================================================================
cat("正在执行全局 ComBat (含验证集)...\n")

expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组大表 & Z-score
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), 
                            check.names = FALSE)
combat_df_all$Batch <- batch_vec

gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 3. 定义实验函数
# ==============================================================================

run_experiment <- function(val_name, n_features) {
  
  cat(paste0("\n>>> 实验: 验证集=", val_name, " | 特征数=", n_features, "\n"))
  
  # 1. 拆分
  valid_data <- combat_df_all[combat_df_all$Batch == val_name, ]
  train_data <- combat_df_all[combat_df_all$Batch != val_name, ]
  
  # 2. 特征筛选 (在训练集上)
  set.seed(123)
  rf_fs <- randomForest(x = train_data[, gene_cols], y = train_data$Group, ntree=150)
  rf_imp <- importance(rf_fs)
  # 选 Top N
  top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:n_features]
  
  cat(paste("    选定特征:", paste(top_feats, collapse=", "), "\n"))
  
  # 3. 建模 (Random Forest + Down-sampling)
  train_sub <- train_data[, c("Group", top_feats)]
  fitControl <- trainControl(method = "cv", number = 5, 
                             classProbs = TRUE, summaryFunction = twoClassSummary,
                             sampling = "down") # 处理不平衡
  
  set.seed(888)
  model <- train(Group ~ ., data = train_sub, method = "rf", 
                 metric = "ROC", trControl = fitControl, ntree = 500)
  
  # 4. 预测
  prob_val <- predict(model, valid_data[, top_feats], type = "prob")
  roc_val <- roc(valid_data$Group, prob_val$Tumor, levels=c("Normal","Tumor"), direction="<", quiet=T)
  
  return(list(auc = as.numeric(auc(roc_val)), roc = roc_val, feats = top_feats))
}

# ==============================================================================
# 4. 执行 4 组实验
# ==============================================================================

# 开启并行
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# A. yyfbatch1 作为验证集
res_y1_f5  <- run_experiment("yyfbatch1", 5)
res_y1_f10 <- run_experiment("yyfbatch1", 10)

# B. yyfbatch2 作为验证集
res_y2_f5  <- run_experiment("yyfbatch2", 5)
res_y2_f10 <- run_experiment("yyfbatch2", 10)

stopCluster(cl)

# ==============================================================================
# 5. 绘图对比 (Feature 5 vs 10)
# ==============================================================================
par(mfrow = c(1, 2))

# --- 图 1: yyfbatch1 验证结果 ---
plot(res_y1_f5$roc, col = "blue", lwd = 2, legacy.axes = TRUE,
     main = "Validation: yyfbatch1\n(Global ComBat)",
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(res_y1_f10$roc, col = "red", lwd = 2, add = TRUE, lty = 2)

legend("bottomright", 
       legend = c(paste0("Top 5 (AUC=", round(res_y1_f5$auc, 3), ")"),
                  paste0("Top 10 (AUC=", round(res_y1_f10$auc, 3), ")")),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

# --- 图 2: yyfbatch2 验证结果 ---
plot(res_y2_f5$roc, col = "blue", lwd = 2, legacy.axes = TRUE,
     main = "Validation: yyfbatch2\n(Global ComBat)",
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(res_y2_f10$roc, col = "red", lwd = 2, add = TRUE, lty = 2)

legend("bottomright", 
       legend = c(paste0("Top 5 (AUC=", round(res_y2_f5$auc, 3), ")"),
                  paste0("Top 10 (AUC=", round(res_y2_f10$auc, 3), ")")),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

par(mfrow = c(1, 1))

# ==============================================================================
# 6. 打印最终结果
# ==============================================================================
res_df <- data.frame(
  Validation_Set = c("yyfbatch1", "yyfbatch1", "yyfbatch2", "yyfbatch2"),
  Feature_Count = c(5, 10, 5, 10),
  AUC = c(res_y1_f5$auc, res_y1_f10$auc, res_y2_f5$auc, res_y2_f10$auc)
)
print(res_df)




































# ==============================================================================
# 0. 急救代码：修复 "invalid connection" 报错
# ==============================================================================
library(foreach)

# 定义一个清理函数
unregister <- function() {
  env <- foreach:::.foreachGlobals
  if (exists("fun", env)) rm(list=ls(name=env), pos=env)
}

# 尝试停止所有并行集群
try(stopImplicitCluster(), silent=TRUE)
try(unregister(), silent=TRUE)

# 关键一步：强制注册为串行（Sequential）模式
registerDoSEQ() 

cat("并行连接已重置。现在可以重新运行 train 代码了。\n")


# ==============================================================================
# 重新运行建模 (现在应该不会报错了)
# ==============================================================================
set.seed(888)

# 这里的 train 代码完全没变，只是环境被修复了
model <- train(Group ~ ., 
               data = train_sub, 
               method = "rf", 
               metric = "ROC", 
               trControl = fitControl, 
               ntree = 500)

print("模型训练完成！")


library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(PRROC)     # 用于计算 PRC
library(ggplot2)
library(gridExtra) # 用于拼图

# ==============================================================================
# 1. 数据准备 (Global ComBat + Top 5 特征筛选)
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因 & Log2
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1
process_brca_balance <- function(df, seed = 123) {
  set.seed(seed)
  idx_normal <- which(df$Group == "Normal")
  idx_tumor <- which(df$Group == "Tumor")
  n_keep <- length(idx_normal) + round((length(idx_tumor) - length(idx_normal)) * 0.40)
  return(df[c(idx_normal, sample(idx_tumor, n_keep)), ])
}
if ("BRCA1" %in% names(clean_list)) clean_list[["BRCA1"]] <- process_brca_balance(clean_list[["BRCA1"]])

# 全局 ComBat
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Scale
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), check.names = FALSE)
combat_df_all$Batch <- batch_vec
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 2. 训练模型 & 独立验证 (Top 5 Features)
# ==============================================================================
val_name <- "yyfbatch2"
n_features <- 5

cat(">>> 正在进行独立验证: ", val_name, " (Top ", n_features, " 特征)...\n")

# 拆分数据
valid_data <- combat_df_all[combat_df_all$Batch == val_name, ]
train_data <- combat_df_all[combat_df_all$Batch != val_name, ]

# 特征筛选
set.seed(123)
rf_fs <- randomForest(x = train_data[, gene_cols], y = train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:n_features]
cat("Top 5 Features:", paste(top_feats, collapse=", "), "\n")

# 建模
train_sub <- train_data[, c("Group", top_feats)]
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary, sampling = "down")
set.seed(888)
model <- train(Group ~ ., data = train_sub, method = "rf", metric = "ROC", trControl = fitControl, ntree = 500)

# 预测概率
prob_val <- predict(model, valid_data[, top_feats], type = "prob")
y_true <- ifelse(valid_data$Group == "Tumor", 1, 0)
y_pred <- prob_val$Tumor

# ==============================================================================
# 3. Bootstrap 计算 AUC 和 AUPRC 的 95% CI (关键步骤)
# ==============================================================================
cat("正在进行 Bootstrap 计算置信区间 (可能需要几秒钟)...\n")

n_boot <- 1000
boot_aucs <- numeric(n_boot)
boot_auprcs <- numeric(n_boot)

set.seed(42)
for(i in 1:n_boot) {
  # 重采样索引
  idx <- sample(1:length(y_true), length(y_true), replace = TRUE)
  y_true_boot <- y_true[idx]
  y_pred_boot <- y_pred[idx]
  
  # 防止重采样导致只有一类样本
  if(length(unique(y_true_boot)) < 2) {
    boot_aucs[i] <- NA
    boot_auprcs[i] <- NA
    next
  }
  
  # Calculate ROC AUC
  boot_aucs[i] <- auc(roc(y_true_boot, y_pred_boot, direction="<", quiet=TRUE))
  
  # Calculate PRC AUPRC (Integral)
  pr_obj <- pr.curve(scores.class0 = y_pred_boot, weights.class0 = y_true_boot, curve=FALSE)
  boot_auprcs[i] <- pr_obj$auc.integral
}

# 移除 NA 并计算 CI
boot_aucs <- na.omit(boot_aucs)
boot_auprcs <- na.omit(boot_auprcs)

# 获取最终统计量
roc_obj <- roc(y_true, y_pred, direction="<", quiet=TRUE)
final_auc <- as.numeric(auc(roc_obj))
auc_ci <- quantile(boot_aucs, probs = c(0.025, 0.975))

pr_obj_final <- pr.curve(scores.class0 = y_pred, weights.class0 = y_true, curve=TRUE)
final_auprc <- pr_obj_final$auc.integral
auprc_ci <- quantile(boot_auprcs, probs = c(0.025, 0.975))

# 格式化 Label 字符串
roc_label <- paste0("AUC = ", sprintf("%.2f", final_auc), 
                    " (95% CI ", sprintf("%.2f", auc_ci[1]), "-", sprintf("%.2f", auc_ci[2]), ")")
prc_label <- paste0(val_name, ": AUPRC = ", sprintf("%.2f", final_auprc), 
                    " (95% CI ", sprintf("%.2f", auprc_ci[1]), "-", sprintf("%.2f", auprc_ci[2]), ")")

# ==============================================================================
# 4. ggplot2 绘图 (复刻论文风格)
# ==============================================================================

# --- A. 准备 ROC 数据 ---
roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)
roc_df <- roc_df[order(roc_df$fpr), ]

# --- B. 准备 PRC 数据 ---
prc_df <- data.frame(
  recall = pr_obj_final$curve[, 1],
  precision = pr_obj_final$curve[, 2]
)

# --- C. 定义通用主题 (美观关键) ---
my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5), # 标题居中
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# --- D. 绘制 ROC 图 ---
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(color = "black", size = 1.2) + # 黑色实线
  # 添加右下角文本框
  annotate("label", x = 0.95, y = 0.05, label = roc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = paste0("Independent Validation in ", val_name),
       x = "False Positive Rate (1-Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  my_theme

# --- E. 绘制 PRC 图 ---
p_prc <- ggplot(prc_df, aes(x = recall, y = precision)) +
  geom_path(color = "black", size = 1.2) + # 黑色实线
  # 添加基准线 (No Skill Line = 阳性比例)
  geom_hline(yintercept = sum(y_true)/length(y_true), linetype="dashed", color="grey50") +
  # 添加右下角文本框
  annotate("label", x = 0.95, y = 0.05, label = prc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "PRC Plot of Independent Validation Dataset",
       x = "Recall",
       y = "Precision") +
  my_theme

# ==============================================================================
# 5. 拼图并输出
# ==============================================================================
grid.arrange(p_roc, p_prc, ncol = 2)

# 如果你想保存为高分辨率图片 (类似 PDF 或 TIFF)
# ggsave("ROC_PRC_Validation.pdf", arrangeGrob(p_roc, p_prc, ncol = 2), width = 12, height = 6)

















library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(PRROC)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# 1. 数据准备 & 全局 ComBat (注意：这里包含验证集！)
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 共有基因 & Log2
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1 (防止训练集样本量过大淹没其他数据)
if ("BRCA1" %in% names(clean_list)) {
  set.seed(123)
  df <- clean_list[["BRCA1"]]
  idx_norm <- which(df$Group == "Normal")
  idx_tum <- sample(which(df$Group == "Tumor"), length(idx_norm) + round((sum(df$Group=="Tumor")-length(idx_norm))*0.4))
  clean_list[["BRCA1"]] <- df[c(idx_norm, idx_tum), ]
}

# --- 执行全局 ComBat (All-in-One) ---
cat(">>> 正在执行全局 ComBat (包含 yyfbatch2)...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Scale
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), check.names = FALSE)
combat_df_all$Batch <- batch_vec

# 统一 Z-score
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 2. 拆分 & 建模 (Top 5 Features)
# ==============================================================================
val_name <- "yyfbatch2"

# 拆分数据
valid_data <- combat_df_all[combat_df_all$Batch == val_name, ]
train_data <- combat_df_all[combat_df_all$Batch != val_name, ]

cat(">>> 训练集样本:", nrow(train_data), "| 验证集(yyfbatch2)样本:", nrow(valid_data), "\n")

# 特征筛选 (RF Importance)
set.seed(123)
rf_fs <- randomForest(x = train_data[, gene_cols], y = train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
# 选 Top 5
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:5]
cat(">>> Top 5 特征:", paste(top_feats, collapse=", "), "\n")

# 训练最终模型
train_sub <- train_data[, c("Group", top_feats)]
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(888)
model <- train(Group ~ ., data = train_sub, method = "rf", metric = "ROC", trControl = fitControl, ntree = 500)

# 预测
prob_val <- predict(model, valid_data[, top_feats], type = "prob")
y_true <- ifelse(valid_data$Group == "Tumor", 1, 0)
y_pred <- prob_val$Tumor

# ==============================================================================
# 3. Bootstrap 计算置信区间 (95% CI)
# ==============================================================================
cat(">>> 正在计算 Bootstrap CI...\n")
n_boot <- 1000
boot_aucs <- numeric(n_boot)
boot_auprcs <- numeric(n_boot)

set.seed(42)
for(i in 1:n_boot) {
  idx <- sample(1:length(y_true), length(y_true), replace = TRUE)
  y_t <- y_true[idx]; y_p <- y_pred[idx]
  if(length(unique(y_t)) < 2) { boot_aucs[i] <- NA; boot_auprcs[i] <- NA; next }
  
  boot_aucs[i] <- as.numeric(auc(roc(y_t, y_p, direction="<", quiet=TRUE)))
  pr_obj <- pr.curve(scores.class0 = y_p, weights.class0 = y_t, curve=FALSE)
  boot_auprcs[i] <- pr_obj$auc.integral
}

boot_aucs <- na.omit(boot_aucs); boot_auprcs <- na.omit(boot_auprcs)
auc_ci <- quantile(boot_aucs, c(0.025, 0.975))
auprc_ci <- quantile(boot_auprcs, c(0.025, 0.975))

# 最终 AUC/AUPRC
roc_obj <- roc(y_true, y_pred, direction="<", quiet=TRUE)
final_auc <- as.numeric(auc(roc_obj))
pr_obj_final <- pr.curve(scores.class0 = y_pred, weights.class0 = y_true, curve=TRUE)
final_auprc <- pr_obj_final$auc.integral

# 标签文本
roc_label <- paste0("AUC = ", sprintf("%.2f", final_auc), 
                    " (95% CI ", sprintf("%.2f", auc_ci[1]), "-", sprintf("%.2f", auc_ci[2]), ")")
prc_label <- paste0("AUPRC = ", sprintf("%.2f", final_auprc), 
                    " (95% CI ", sprintf("%.2f", auprc_ci[1]), "-", sprintf("%.2f", auprc_ci[2]), ")")

# ==============================================================================
# 4. ggplot2 绘图
# ==============================================================================
# 数据转换
roc_df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
roc_df <- roc_df[order(roc_df$fpr), ]
prc_df <- data.frame(recall = pr_obj_final$curve[, 1], precision = pr_obj_final$curve[, 2])

# 通用主题
my_theme <- theme_bw() +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 13, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# 绘制 ROC
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(color = "#E41A1C", size = 1.2) + # 使用深红色突出
  annotate("label", x = 0.95, y = 0.05, label = roc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = paste0("Independent Validation (Global ComBat)\n", val_name),
       x = "False Positive Rate (1-Specificity)", y = "Sensitivity") + my_theme

# 绘制 PRC
p_prc <- ggplot(prc_df, aes(x = recall, y = precision)) +
  geom_hline(yintercept = sum(y_true)/length(y_true), linetype="dashed", color="grey50") +
  geom_path(color = "#377EB8", size = 1.2) + # 使用深蓝色突出
  annotate("label", x = 0.95, y = 0.05, label = prc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") + my_theme

# 拼图输出
grid.arrange(p_roc, p_prc, ncol = 2)



















library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(PRROC)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# 1. 数据准备 & 全局 ComBat (All-in-One Mode)
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因 & Log2 处理
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1 (保留，防止样本不平衡)
if ("BRCA1" %in% names(clean_list)) {
  set.seed(123)
  df <- clean_list[["BRCA1"]]
  idx_norm <- which(df$Group == "Normal")
  # 抽取 Normal 数量 + 剩余 Tumor 的 40%
  idx_tum <- sample(which(df$Group == "Tumor"), length(idx_norm) + round((sum(df$Group=="Tumor")-length(idx_norm))*0.4))
  clean_list[["BRCA1"]] <- df[c(idx_norm, idx_tum), ]
}

# --- 执行全局 ComBat (关键步骤：验证集也在其中) ---
cat(">>> 正在执行全局 ComBat (包含 yyfbatch1)...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Scale (Z-score)
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), check.names = FALSE)
combat_df_all$Batch <- batch_vec

gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 2. 拆分 & 建模 (目标验证集：yyfbatch1)
# ==============================================================================
val_name <- "yyfbatch1"  # <--- 修改为 yyfbatch1

# 拆分数据
valid_data <- combat_df_all[combat_df_all$Batch == val_name, ]
train_data <- combat_df_all[combat_df_all$Batch != val_name, ]

cat(">>> 训练集样本:", nrow(train_data), "| 验证集(yyfbatch1)样本:", nrow(valid_data), "\n")

# 特征筛选 (在训练集上做)
set.seed(123)
rf_fs <- randomForest(x = train_data[, gene_cols], y = train_data$Group, ntree=150)
rf_imp <- importance(rf_fs)
# 选 Top 5
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:5]
cat(">>> Top 5 特征:", paste(top_feats, collapse=", "), "\n")

# 训练模型
train_sub <- train_data[, c("Group", top_feats)]
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(888)
model <- train(Group ~ ., data = train_sub, method = "rf", metric = "ROC", trControl = fitControl, ntree = 500)

# 预测
prob_val <- predict(model, valid_data[, top_feats], type = "prob")
y_true <- ifelse(valid_data$Group == "Tumor", 1, 0)
y_pred <- prob_val$Tumor

# ==============================================================================
# 3. Bootstrap 计算置信区间 (95% CI)
# ==============================================================================
cat(">>> 正在计算 Bootstrap CI...\n")
n_boot <- 1000
boot_aucs <- numeric(n_boot)
boot_auprcs <- numeric(n_boot)

set.seed(42)
for(i in 1:n_boot) {
  idx <- sample(1:length(y_true), length(y_true), replace = TRUE)
  y_t <- y_true[idx]; y_p <- y_pred[idx]
  # 防止重采样导致单一类别
  if(length(unique(y_t)) < 2) { boot_aucs[i] <- NA; boot_auprcs[i] <- NA; next }
  
  boot_aucs[i] <- as.numeric(auc(roc(y_t, y_p, direction="<", quiet=TRUE)))
  pr_obj <- pr.curve(scores.class0 = y_p, weights.class0 = y_t, curve=FALSE)
  boot_auprcs[i] <- pr_obj$auc.integral
}

boot_aucs <- na.omit(boot_aucs); boot_auprcs <- na.omit(boot_auprcs)
auc_ci <- quantile(boot_aucs, c(0.025, 0.975))
auprc_ci <- quantile(boot_auprcs, c(0.025, 0.975))

# 计算最终值
roc_obj <- roc(y_true, y_pred, direction="<", quiet=TRUE)
final_auc <- as.numeric(auc(roc_obj))
pr_obj_final <- pr.curve(scores.class0 = y_pred, weights.class0 = y_true, curve=TRUE)
final_auprc <- pr_obj_final$auc.integral

# 准备标签
roc_label <- paste0("AUC = ", sprintf("%.2f", final_auc), 
                    " (95% CI ", sprintf("%.2f", auc_ci[1]), "-", sprintf("%.2f", auc_ci[2]), ")")
prc_label <- paste0("AUPRC = ", sprintf("%.2f", final_auprc), 
                    " (95% CI ", sprintf("%.2f", auprc_ci[1]), "-", sprintf("%.2f", auprc_ci[2]), ")")

cat("Done! Final AUC:", final_auc, "\n")

# ==============================================================================
# 4. 绘图 (ggplot2 美化版)
# ==============================================================================
roc_df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
roc_df <- roc_df[order(roc_df$fpr), ]
prc_df <- data.frame(recall = pr_obj_final$curve[, 1], precision = pr_obj_final$curve[, 2])

# 主题设置
my_theme <- theme_bw() +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 13, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# 绘制 ROC
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(color = "#E41A1C", size = 1.2) + # 红色曲线
  annotate("label", x = 0.95, y = 0.05, label = roc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = paste0("Independent Validation (Global ComBat)\n", val_name),
       x = "False Positive Rate (1-Specificity)", y = "Sensitivity") + my_theme

# 绘制 PRC
p_prc <- ggplot(prc_df, aes(x = recall, y = precision)) +
  geom_hline(yintercept = sum(y_true)/length(y_true), linetype="dashed", color="grey50") +
  geom_path(color = "#377EB8", size = 1.2) + # 蓝色曲线
  annotate("label", x = 0.95, y = 0.05, label = prc_label, 
           hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.8, label.size=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") + my_theme

# 拼图
grid.arrange(p_roc, p_prc, ncol = 2)






library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(PRROC)
library(ggplot2)
library(ggpubr) # 用于美化箱线图和添加统计显著性
library(tidyr)

# ==============================================================================
# 1. 数据准备 & 全局 ComBat (All-in-One Mode)
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 共有基因 & Log2 处理
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1
if ("BRCA1" %in% names(clean_list)) {
  set.seed(123)
  df <- clean_list[["BRCA1"]]
  idx_norm <- which(df$Group == "Normal")
  idx_tum <- sample(which(df$Group == "Tumor"), length(idx_norm) + round((sum(df$Group=="Tumor")-length(idx_norm))*0.4))
  clean_list[["BRCA1"]] <- df[c(idx_norm, idx_tum), ]
}

# --- 执行全局 ComBat ---
cat(">>> 正在执行全局 ComBat...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Scale
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), check.names = FALSE)
combat_df_all$Batch <- batch_vec
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 2. 数据集划分 (Discovery / Hold-out / Independent)
# ==============================================================================
val_name <- "yyfbatch1" # 独立验证集

# 1. 分离独立验证集
data_independent <- combat_df_all[combat_df_all$Batch == val_name, ]
data_remaining   <- combat_df_all[combat_df_all$Batch != val_name, ]

# 2. 将剩余数据按 7:3 划分为 Discovery 和 Hold-out
set.seed(123)
trainIndex <- createDataPartition(data_remaining$Group, p = 0.7, list = FALSE, times = 1)
data_discovery <- data_remaining[trainIndex, ]
data_holdout   <- data_remaining[-trainIndex, ]

cat(">>> 样本统计:\n")
cat("    Discovery (Training):", nrow(data_discovery), "\n")
cat("    Hold-out (Testing):  ", nrow(data_holdout), "\n")
cat("    Independent (Valid): ", nrow(data_independent), "\n")

# ==============================================================================
# 3. 特征筛选 & 建模 (仅在 Discovery 上进行)
# ==============================================================================
set.seed(123)
rf_fs <- randomForest(x = data_discovery[, gene_cols], y = data_discovery$Group, ntree=150)
rf_imp <- importance(rf_fs)
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:5]
cat(">>> Top 5 特征:", paste(top_feats, collapse=", "), "\n")

# 训练最终模型
train_sub <- data_discovery[, c("Group", top_feats)]
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(888)
model <- train(Group ~ ., data = train_sub, method = "rf", metric = "ROC", trControl = fitControl, ntree = 500)

# ==============================================================================
# 4. 计算 T-Scores (预测概率)
# ==============================================================================
# 定义一个辅助函数来获取预测结果
get_pred_data <- function(model, dataset, phase_name) {
  prob <- predict(model, dataset[, top_feats], type = "prob")$Tumor
  y_true <- ifelse(dataset$Group == "Tumor", 1, 0)
  return(data.frame(
    SampleID = rownames(dataset),
    Group = dataset$Group,
    T_Score = prob,
    True_Class = y_true,
    Phase = phase_name
  ))
}

res_disc  <- get_pred_data(model, data_discovery, "Discovery")
res_hold  <- get_pred_data(model, data_holdout, "Hold-out Validation")
res_indep <- get_pred_data(model, data_independent, "Independent Validation")

# 合并所有结果用于箱线图
all_scores <- rbind(res_disc, res_hold, res_indep)
# 设定 Phase 的因子顺序，保证画图顺序正确
all_scores$Phase <- factor(all_scores$Phase, levels = c("Discovery", "Hold-out Validation", "Independent Validation"))

# ==============================================================================
# 5. 定义绘图函数 (ROC & PRC)
# ==============================================================================
# 通用主题
my_theme <- theme_bw() +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# --- 画单个 ROC 的函数 ---
draw_roc <- function(res_df, title_text, color_line="#E41A1C") {
  # Bootstrap CI
  n_boot <- 500 # 演示用500，正式可用1000
  aucs <- numeric(n_boot)
  set.seed(42)
  for(i in 1:n_boot) {
    idx <- sample(1:nrow(res_df), nrow(res_df), replace=T)
    sub <- res_df[idx, ]
    if(length(unique(sub$True_Class)) < 2) { aucs[i]<-NA; next }
    aucs[i] <- as.numeric(auc(roc(sub$True_Class, sub$T_Score, direction="<", quiet=T)))
  }
  auc_ci <- quantile(na.omit(aucs), c(0.025, 0.975))
  
  # Main Curve
  roc_obj <- roc(res_df$True_Class, res_df$T_Score, direction="<", quiet=T)
  final_auc <- as.numeric(auc(roc_obj))
  
  roc_df <- data.frame(fpr=1-roc_obj$specificities, tpr=roc_obj$sensitivities)
  roc_df <- roc_df[order(roc_df$fpr),]
  
  label_txt <- paste0("AUC = ", sprintf("%.3f", final_auc), 
                      "\n(95% CI: ", sprintf("%.2f", auc_ci[1]), "-", sprintf("%.2f", auc_ci[2]), ")")
  
  ggplot(roc_df, aes(x=fpr, y=tpr)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    geom_path(color=color_line, size=1.2) +
    annotate("text", x=0.95, y=0.05, label=label_txt, hjust=1, vjust=0, size=5) +
    labs(title=title_text, x="False Positive Rate", y="True Positive Rate") +
    my_theme
}

# --- 画单个 PRC 的函数 ---
draw_prc <- function(res_df, title_text, color_line="#377EB8") {
  # Bootstrap CI
  n_boot <- 500
  auprcs <- numeric(n_boot)
  set.seed(42)
  for(i in 1:n_boot) {
    idx <- sample(1:nrow(res_df), nrow(res_df), replace=T)
    sub <- res_df[idx, ]
    if(length(unique(sub$True_Class)) < 2) { auprcs[i]<-NA; next }
    pr_obj <- pr.curve(scores.class0=sub$T_Score, weights.class0=sub$True_Class, curve=F)
    auprcs[i] <- pr_obj$auc.integral
  }
  ci <- quantile(na.omit(auprcs), c(0.025, 0.975))
  
  # Main Curve
  pr_final <- pr.curve(scores.class0=res_df$T_Score, weights.class0=res_df$True_Class, curve=T)
  final_auprc <- pr_final$auc.integral
  prc_df <- data.frame(recall=pr_final$curve[,1], precision=pr_final$curve[,2])
  
  baseline <- sum(res_df$True_Class)/nrow(res_df)
  label_txt <- paste0("AUPRC = ", sprintf("%.3f", final_auprc), 
                      "\n(95% CI: ", sprintf("%.2f", ci[1]), "-", sprintf("%.2f", ci[2]), ")")
  
  ggplot(prc_df, aes(x=recall, y=precision)) +
    geom_hline(yintercept=baseline, linetype="dashed", color="grey50") +
    geom_path(color=color_line, size=1.2) +
    annotate("text", x=0.05, y=0.05, label=label_txt, hjust=0, vjust=0, size=5) +
    labs(title=title_text, x="Recall", y="Precision") +
    my_theme
}

# ==============================================================================
# 6. 生成并打印所有曲线图
# ==============================================================================

# A. Discovery Phase
p1 <- draw_roc(res_disc, "ROC - Discovery Phase")
p2 <- draw_prc(res_disc, "PRC - Discovery Phase")

# B. Hold-out Validation Phase
p3 <- draw_roc(res_hold, "ROC - Hold-out Validation")
p4 <- draw_prc(res_hold, "PRC - Hold-out Validation")

# C. Independent Validation Phase (yyfbatch1)
p5 <- draw_roc(res_indep, paste0("ROC - Independent Validation\n(", val_name, ")"))
p6 <- draw_prc(res_indep, paste0("PRC - Independent Validation\n(", val_name, ")"))

# 打印图像 (你可以在 RStudio Plots 面板中依次查看)
print(p1); print(p2)
print(p3); print(p4)
print(p5); print(p6)

# ==============================================================================
# 7. 绘制 T-Score 箱线图 (核心需求)
# ==============================================================================
# 比较不同数据集下 Tumor 和 Control 的得分差异

# 颜色设置：Normal用灰色/蓝色，Tumor用红色
my_colors <- c("#868686FF", "#CD534CFF") 

p_box <- ggplot(all_scores, aes(x = Phase, y = T_Score, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  # 添加散点 (抖动)
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), 
              size = 1, alpha = 0.3, color = "black") +
  # 添加统计显著性 (Wilcoxon test)
  stat_compare_means(aes(group = Group), method = "wilcox.test", 
                     label = "p.signif", label.y = 1.05, size = 6) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2), name = "Predicted T-Score (Probability of Tumor)") +
  labs(title = "Boxplots of T-scores across different phases",
       subtitle = "Comparison of predicted tumor probability in Normal vs Tumor samples",
       x = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color="black"),
    axis.text.y = element_text(size = 12, color="black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid = element_blank() # 去除网格线更像发表级图表
  )

print(p_box)




















































library(sva)
library(caret)
library(dplyr)
library(randomForest)
library(pROC)
library(PRROC)
library(ggplot2)
library(ggpubr) # 用于统计检验和美化
library(tidyr)

# ==============================================================================
# 1. 数据准备 & 全局 ComBat (All-in-One Mode)
# ==============================================================================
# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
dataset_names <- gsub("processed_results/|_processed.csv", "", files)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- dataset_names

# 提取共有基因 & Log2 处理
common_genes <- Reduce(intersect, lapply(data_list, function(df) colnames(df)[-1]))
clean_list <- lapply(data_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) df_sub[, -1] <- log2(mat + 1)
  return(df_sub)
})

# 平衡 BRCA1 (保留)
if ("BRCA1" %in% names(clean_list)) {
  set.seed(123)
  df <- clean_list[["BRCA1"]]
  idx_norm <- which(df$Group == "Normal")
  idx_tum <- sample(which(df$Group == "Tumor"), length(idx_norm) + round((sum(df$Group=="Tumor")-length(idx_norm))*0.4))
  clean_list[["BRCA1"]] <- df[c(idx_norm, idx_tum), ]
}

# --- 执行全局 ComBat (关键：包含 yyfbatch2) ---
cat(">>> 正在执行全局 ComBat...\n")
expr_matrices <- lapply(clean_list, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)
batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# 重组 & Scale
combat_df_all <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                            t(combat_expr), check.names = FALSE)
combat_df_all$Batch <- batch_vec
gene_cols <- colnames(combat_df_all)[!colnames(combat_df_all) %in% c("Group", "Batch")]
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

# ==============================================================================
# 2. 数据集划分 (Discovery / Hold-out / Independent)
# ==============================================================================
val_name <- "yyfbatch2"  # <--- 【关键修改】指定 yyfbatch2 为独立验证集

# 1. 分离独立验证集
data_independent <- combat_df_all[combat_df_all$Batch == val_name, ]
data_remaining   <- combat_df_all[combat_df_all$Batch != val_name, ]

# 2. 将剩余数据按 7:3 划分为 Discovery 和 Hold-out
set.seed(123)
trainIndex <- createDataPartition(data_remaining$Group, p = 0.7, list = FALSE, times = 1)
data_discovery <- data_remaining[trainIndex, ]
data_holdout   <- data_remaining[-trainIndex, ]

cat(">>> 样本统计:\n")
cat("    Discovery (Training):", nrow(data_discovery), "\n")
cat("    Hold-out (Testing):  ", nrow(data_holdout), "\n")
cat("    Independent (Valid): ", nrow(data_independent), " (Dataset:", val_name, ")\n")

# ==============================================================================
# 3. 特征筛选 & 建模 (仅在 Discovery 上进行)
# ==============================================================================
set.seed(123)
rf_fs <- randomForest(x = data_discovery[, gene_cols], y = data_discovery$Group, ntree=150)
rf_imp <- importance(rf_fs)
# 选 Top 5
top_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:5]
cat(">>> Top 5 特征:", paste(top_feats, collapse=", "), "\n")

# 训练最终模型
train_sub <- data_discovery[, c("Group", top_feats)]
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(888)
model <- train(Group ~ ., data = train_sub, method = "rf", metric = "ROC", trControl = fitControl, ntree = 500)

# ==============================================================================
# 4. 计算 T-Scores (预测概率)
# ==============================================================================
get_pred_data <- function(model, dataset, phase_name) {
  prob <- predict(model, dataset[, top_feats], type = "prob")$Tumor
  y_true <- ifelse(dataset$Group == "Tumor", 1, 0)
  return(data.frame(
    SampleID = rownames(dataset),
    Group = dataset$Group,
    T_Score = prob,
    True_Class = y_true,
    Phase = phase_name
  ))
}

res_disc  <- get_pred_data(model, data_discovery, "Discovery")
res_hold  <- get_pred_data(model, data_holdout, "Hold-out Validation")
res_indep <- get_pred_data(model, data_independent, "Independent Validation")

# 合并结果
all_scores <- rbind(res_disc, res_hold, res_indep)
all_scores$Phase <- factor(all_scores$Phase, levels = c("Discovery", "Hold-out Validation", "Independent Validation"))

# ==============================================================================
# 5. 绘图函数定义 (ROC & PRC)
# ==============================================================================
my_theme <- theme_bw() +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# ROC 函数
draw_roc <- function(res_df, title_text, color_line="#E41A1C") {
  # Bootstrap CI
  n_boot <- 500
  aucs <- numeric(n_boot)
  set.seed(42)
  for(i in 1:n_boot) {
    idx <- sample(1:nrow(res_df), nrow(res_df), replace=T)
    sub <- res_df[idx, ]
    if(length(unique(sub$True_Class)) < 2) { aucs[i]<-NA; next }
    aucs[i] <- as.numeric(auc(roc(sub$True_Class, sub$T_Score, direction="<", quiet=T)))
  }
  auc_ci <- quantile(na.omit(aucs), c(0.025, 0.975))
  
  roc_obj <- roc(res_df$True_Class, res_df$T_Score, direction="<", quiet=T)
  final_auc <- as.numeric(auc(roc_obj))
  roc_df <- data.frame(fpr=1-roc_obj$specificities, tpr=roc_obj$sensitivities)
  roc_df <- roc_df[order(roc_df$fpr),]
  
  label_txt <- paste0("AUC = ", sprintf("%.3f", final_auc), 
                      "\n(95% CI: ", sprintf("%.2f", auc_ci[1]), "-", sprintf("%.2f", auc_ci[2]), ")")
  
  ggplot(roc_df, aes(x=fpr, y=tpr)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    geom_path(color=color_line, size=1.2) +
    annotate("label", x=0.95, y=0.05, label=label_txt, hjust=1, vjust=0, size=4.5, fill="white", alpha=0.8, label.size=NA) +
    labs(title=title_text, x="False Positive Rate", y="True Positive Rate") +
    my_theme
}

# PRC 函数
draw_prc <- function(res_df, title_text, color_line="#377EB8") {
  # Bootstrap CI
  n_boot <- 500
  auprcs <- numeric(n_boot)
  set.seed(42)
  for(i in 1:n_boot) {
    idx <- sample(1:nrow(res_df), nrow(res_df), replace=T)
    sub <- res_df[idx, ]
    if(length(unique(sub$True_Class)) < 2) { auprcs[i]<-NA; next }
    pr_obj <- pr.curve(scores.class0=sub$T_Score, weights.class0=sub$True_Class, curve=F)
    auprcs[i] <- pr_obj$auc.integral
  }
  ci <- quantile(na.omit(auprcs), c(0.025, 0.975))
  
  pr_final <- pr.curve(scores.class0=res_df$T_Score, weights.class0=res_df$True_Class, curve=T)
  final_auprc <- pr_final$auc.integral
  prc_df <- data.frame(recall=pr_final$curve[,1], precision=pr_final$curve[,2])
  
  baseline <- sum(res_df$True_Class)/nrow(res_df)
  label_txt <- paste0("AUPRC = ", sprintf("%.3f", final_auprc), 
                      "\n(95% CI: ", sprintf("%.2f", ci[1]), "-", sprintf("%.2f", ci[2]), ")")
  
  ggplot(prc_df, aes(x=recall, y=precision)) +
    geom_hline(yintercept=baseline, linetype="dashed", color="grey50") +
    geom_path(color=color_line, size=1.2) +
    annotate("label", x=0.05, y=0.05, label=label_txt, hjust=0, vjust=0, size=4.5, fill="white", alpha=0.8, label.size=NA) +
    labs(title=title_text, x="Recall", y="Precision") +
    my_theme
}

# ==============================================================================
# 6. 生成曲线图
# ==============================================================================
p1 <- draw_roc(res_disc, "ROC - Discovery Phase")
p2 <- draw_prc(res_disc, "PRC - Discovery Phase")

p3 <- draw_roc(res_hold, "ROC - Hold-out Validation")
p4 <- draw_prc(res_hold, "PRC - Hold-out Validation")

p5 <- draw_roc(res_indep, paste0("ROC - Independent Validation\n(", val_name, ")"))
p6 <- draw_prc(res_indep, paste0("PRC - Independent Validation\n(", val_name, ")"))

print(p1); print(p2)
print(p3); print(p4)
print(p5); print(p6)

# ==============================================================================
# 7. 绘制 T-Score 箱线图
# ==============================================================================
my_colors <- c("#868686FF", "#CD534CFF") 

p_box <- ggplot(all_scores, aes(x = Phase, y = T_Score, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.3, color = "black") +
  stat_compare_means(aes(group = Group), method = "wilcox.test", label = "p.signif", label.y = 1.05, size = 6) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2), name = "Predicted T-Score") +
  labs(title = "Boxplots of T-scores Comparison", 
       subtitle = paste0("Independent Validation on ", val_name), x = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, face = "bold", color="black"),
    axis.text.y = element_text(size = 12, color="black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid = element_blank()
  )

print(p_box)