library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 读取数据 (假设之前的 CSV 文件都在 processed_results 文件夹下)
# ==============================================================================
# 确保你的工作目录设置正确
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
# 根据文件名读取数据到列表
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- gsub("processed_results/|_processed.csv", "", files)

# 定义分组
pool_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049", "yyfbatch2")
valid_name <- "yyfbatch1"

# 提取列表
pool_list <- data_list[pool_names]
valid_list <- data_list[valid_name]

# ==============================================================================
# 2. 寻找共有基因 (Common Genes)
# ==============================================================================
all_dfs <- c(pool_list, valid_list)
common_genes <- Reduce(intersect, lapply(all_dfs, function(df) colnames(df)[-1])) # -1 排除 Group

cat("共有基因数量:", length(common_genes), "\n")

# 数据清洗辅助函数 (取共有基因 + Log2)
prep_data <- function(df, genes) {
  df_sub <- df[, c("Group", genes)]
  # 简单的 Log2 转换 (如果最大值 > 50)
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) {
    df_sub[, -1] <- log2(mat + 1)
  }
  return(df_sub)
}

pool_clean <- lapply(pool_list, prep_data, genes = common_genes)
valid_clean <- lapply(valid_list, prep_data, genes = common_genes)[[1]] # 提取出 dataframe

# ==============================================================================
# 3. 对 6 个数据集进行 ComBat 去批次
# ==============================================================================
cat("正在对 6 个训练数据集进行合并和 ComBat...\n")

# A. 构建合并矩阵 (ComBat 需要 基因在行，样本在列)
# 我们需要先 cbind 所有表达矩阵
expr_matrices <- lapply(pool_clean, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)

# B. 构建 Batch 向量 (告诉 ComBat 哪个样本属于哪个数据集)
batch_vec <- unlist(lapply(names(pool_clean), function(n) rep(n, nrow(pool_clean[[n]]))))

# C. 构建分组信息 (保护生物学差异不被去除)
# 必须保证合并后的顺序和 expr_matrices 一致
group_vec <- unlist(lapply(pool_clean, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# D. 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# E. 转置回 (行=样本，列=基因) 并加上 Group
combat_df <- as.data.frame(t(combat_expr))
combat_df$Group <- factor(group_vec, levels = c("Normal", "Tumor"))

# 将 Group 放回第一列
combat_df <- combat_df[, c("Group", colnames(combat_df)[-1])]

cat("ComBat 完成！数据维度:", dim(combat_df)[1], "行,", dim(combat_df)[2], "列\n")




# ==============================================================================
# 修正后的步骤 3 & 4：ComBat 与 数据重组
# ==============================================================================

# ... (前面的 ComBat 计算过程 combat_expr 不变，从这里开始替换) ...

# D. 运行 ComBat (假设你已经运行完了这一步，得到了 combat_expr)
# combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# --- 修正开始 ---

# 1. 将 ComBat 后的矩阵转置，并直接与 Group 拼接
# 这样 Group 自动就在第一列，不需要后续调整顺序，也不会出错
combat_df_final <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                              t(combat_expr), 
                              check.names = FALSE) # 防止基因名里的 "-" 被变成 "."

# 检查一下结构，确保全是数字 (除了第一列)
# str(combat_df_final[, 1:5]) 

# 2. Z-score 标准化
# 注意：scale 函数返回的是矩阵，我们需要把它塞回 dataframe
# combat_df_final[, -1] 取出除 Group 外的所有列
combat_scaled_mat <- scale(combat_df_final[, -1])

# 3. 重新组装为最终的 Pool 数据
pool_final <- data.frame(Group = combat_df_final$Group, 
                         combat_scaled_mat,
                         check.names = FALSE)

# --- 修正结束 ---

cat("修正完成！数据维度:", dim(pool_final)[1], "行,", dim(pool_final)[2], "列\n")
cat("前5列的数据类型:", paste(sapply(pool_final[,1:5], class), collapse=", "), "\n") 
# 输出应该是: factor, numeric, numeric, numeric, numeric




# B. 对 yyfbatch1 单独做 Scale (模拟临床单中心标准化)
valid_scaled_mat <- scale(valid_clean[, -1])
valid_final <- data.frame(Group = as.factor(valid_clean$Group), valid_scaled_mat)
valid_final$Group <- factor(valid_final$Group, levels = c("Normal", "Tumor"))

# ==============================================================================
# 5. 70% / 30% 拆分
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(pool_final$Group, p = 0.7, list = FALSE, times = 1)

train_set <- pool_final[trainIndex, ]
internal_test_set <- pool_final[-trainIndex, ]

cat("训练集 (ComBat+Scaled):", nrow(train_set), "\n")
cat("内验集 (ComBat+Scaled):", nrow(internal_test_set), "\n")
cat("独立外验集 (yyfbatch1 Scaled):", nrow(valid_final), "\n")




# ==============================================================================
# 6. 特征筛选 (仅在 70% 训练集上)
# ==============================================================================
cat("开始特征筛选...\n")
x <- as.matrix(train_set[, -1])
y <- train_set$Group

# 简单演示：LASSO + Random Forest 投票
# (你可以把之前复杂的7种方法放回来，这里为了代码简洁只写两种)

# LASSO
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)
coefs <- coef(cv_lasso, s = "lambda.min")
lasso_feats <- rownames(coefs)[coefs[,1] != 0][-1]

# Random Forest
rf_fs <- randomForest(x, y, ntree=300)
rf_imp <- importance(rf_fs)
rf_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:20]

# 取并集 (或者交集，取决于你想筛多少)
final_miRNA_set <- unique(c(lasso_feats, rf_feats))
# 如果太多，只取前10个最重要的
if(length(final_miRNA_set) > 10) final_miRNA_set <- final_miRNA_set[1:10]

cat("最终筛选特征:", paste(final_miRNA_set, collapse=", "), "\n")

# ==============================================================================
# 7. 训练最佳模型 (解决样本不平衡)
# ==============================================================================
train_subset <- train_set[, c("Group", final_miRNA_set)]

# 开启多核
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

fitControl <- trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           sampling = "down") # <--- 关键：解决BRCA1不平衡

cat("开始训练 Random Forest (带 Down-sampling)...\n")
set.seed(888)
final_model <- train(Group ~ ., 
                     data = train_subset, 
                     method = "rf", 
                     metric = "ROC", 
                     trControl = fitControl,
                     ntree = 500)

stopCluster(cl)
print(final_model)







# ==============================================================================
# 8. 验证与绘图
# ==============================================================================

# 定义评估函数
eval_model <- function(model, data, features, label) {
  data_sub <- data[, c("Group", features)]
  probs <- predict(model, data_sub, type = "prob")
  roc_res <- roc(data_sub$Group, probs$Tumor, levels=c("Normal", "Tumor"), direction="<")
  return(list(roc=roc_res, auc=auc(roc_res)))
}

# A. 内部验证 (30% Split)
res_internal <- eval_model(final_model, internal_test_set, final_miRNA_set, "Internal")

# B. 外部独立验证 (yyfbatch1)
res_external <- eval_model(final_model, valid_final, final_miRNA_set, "External")

# 打印结果
cat("\n----------------------------------\n")
cat("内部验证集 (30%) AUC:", res_internal$auc, "\n")
cat("独立验证集 (yyfbatch1) AUC:", res_external$auc, "\n")
cat("----------------------------------\n")

# 绘制 ROC
plot(res_internal$roc, col="blue", lwd=2, main="ROC Curves Performance")
plot(res_external$roc, col="red", lwd=2, add=TRUE)
legend("bottomright", 
       legend=c(paste0("Internal Valid (AUC=", round(res_internal$auc, 3), ")"),
                paste0("External yyfbatch1 (AUC=", round(res_external$auc, 3), ")")),
       col=c("blue", "red"), lwd=2)

# 输出混淆矩阵查看特异性
cat("\nyyfbatch1 混淆矩阵:\n")
pred_raw <- predict(final_model, valid_final[, c("Group", final_miRNA_set)], type="raw")
print(confusionMatrix(pred_raw, valid_final$Group, positive="Tumor"))













library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 准备数据名单 (Swap: yyfbatch1 进训练池, yyfbatch2 做验证)
# ==============================================================================
# 训练池：BRCA1 + 4个PRJNA + yyfbatch1 (共6个)
pool_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049", "yyfbatch1")

# 独立验证集：yyfbatch2
valid_name <- "yyfbatch2"

# 读取数据 (假设在 processed_results 目录下)
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- gsub("processed_results/|_processed.csv", "", files)

# 提取列表
pool_list <- data_list[pool_names]
valid_list <- data_list[valid_name]

# ==============================================================================
# 2. 数据清洗与提取共有基因
# ==============================================================================
# 寻找所有 7 个数据集共有的基因
all_dfs <- c(pool_list, valid_list)
common_genes <- Reduce(intersect, lapply(all_dfs, function(df) colnames(df)[-1])) # -1 排除 Group

cat("共有 miRNA 数量:", length(common_genes), "\n")

# 清洗函数 (取共有基因 + Log2处理)
prep_data <- function(df, genes) {
  df_sub <- df[, c("Group", genes)]
  # 检查是否需要 Log2 (简单阈值判断)
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) {
    df_sub[, -1] <- log2(mat + 1)
  }
  return(df_sub)
}

pool_clean <- lapply(pool_list, prep_data, genes = common_genes)
valid_clean <- lapply(valid_list, prep_data, genes = common_genes)[[1]]

# ==============================================================================
# 3. 对 6 个训练数据集进行 ComBat 去批次 (关键步骤)
# ==============================================================================
cat("正在对 6 个训练数据集 (含yyfbatch1) 进行 ComBat...\n")

# A. 构建合并矩阵 (ComBat 要求: 行=基因, 列=样本)
expr_matrices <- lapply(pool_clean, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)

# B. 构建 Batch 向量
batch_vec <- unlist(lapply(names(pool_clean), function(n) rep(n, nrow(pool_clean[[n]]))))

# C. 构建分组信息 (保护生物学差异)
group_vec <- unlist(lapply(pool_clean, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# D. 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# ==============================================================================
# 4. 数据重组与 Z-score 标准化 (修复之前的报错)
# ==============================================================================

# A. 重组训练池数据 (ComBat后)
# 关键修复: 直接在 data.frame 构建时放入 Group，并使用 check.names=FALSE 保护 miRNA 名字
combat_df_final <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                              t(combat_expr), 
                              check.names = FALSE)

# Scale (Z-score) 训练池
# 注意：combat_df_final[,-1] 确保只取数值列进行计算
combat_scaled_mat <- scale(combat_df_final[, -1])
pool_final <- data.frame(Group = combat_df_final$Group, 
                         combat_scaled_mat, 
                         check.names = FALSE)

# B. 处理独立验证集 (yyfbatch2)
# 验证集不做 ComBat，但必须做 Z-score 以对齐量纲
valid_scaled_mat <- scale(valid_clean[, -1])
valid_final <- data.frame(Group = factor(valid_clean$Group, levels = c("Normal", "Tumor")), 
                          valid_scaled_mat, 
                          check.names = FALSE)

# ==============================================================================
# 5. 70% 训练 / 30% 内验 拆分
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(pool_final$Group, p = 0.7, list = FALSE, times = 1)

train_set <- pool_final[trainIndex, ]
internal_test_set <- pool_final[-trainIndex, ]

cat("训练集 (ComBat+Scaled):", nrow(train_set), "\n")
cat("内验集 (ComBat+Scaled):", nrow(internal_test_set), "\n")
cat("外验集 (yyfbatch2 Scaled):", nrow(valid_final), "\n")

# ==============================================================================
# 6. 特征筛选 (仅在训练集上)
# ==============================================================================
cat("开始特征筛选 (LASSO + RF)...\n")
x <- as.matrix(train_set[, -1])
y <- train_set$Group

# 1. LASSO
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)
coefs <- coef(cv_lasso, s = "lambda.min")
lasso_feats <- rownames(coefs)[coefs[,1] != 0][-1]

# 2. Random Forest Importance
rf_fs <- randomForest(x, y, ntree=300)
rf_imp <- importance(rf_fs)
rf_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:20]

# 取并集
final_miRNA_set <- unique(c(lasso_feats, rf_feats))
# 如果数量太多，限制在前15个
if(length(final_miRNA_set) > 15) final_miRNA_set <- final_miRNA_set[1:15]

cat("最终特征 Panel:", paste(final_miRNA_set, collapse=", "), "\n")

# ==============================================================================
# 7. 建模 (Random Forest + Down-sampling)
# ==============================================================================
train_subset <- train_set[, c("Group", final_miRNA_set)]

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# 使用 Down-sampling 解决 BRCA1 带来的不平衡
fitControl <- trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           sampling = "down") 

cat("开始训练模型...\n")
set.seed(888)
final_model <- train(Group ~ ., 
                     data = train_subset, 
                     method = "rf", 
                     metric = "ROC", 
                     trControl = fitControl,
                     ntree = 500)

stopCluster(cl)

# ==============================================================================
# 8. 最终验证与画图
# ==============================================================================
# 定义绘图函数
plot_roc_comparison <- function(model, internal_data, external_data, features) {
  
  # 内部验证
  prob_int <- predict(model, internal_data[, c("Group", features)], type = "prob")
  roc_int <- roc(internal_data$Group, prob_int$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 外部验证 (yyfbatch2)
  prob_ext <- predict(model, external_data[, c("Group", features)], type = "prob")
  roc_ext <- roc(external_data$Group, prob_ext$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 绘图
  plot(roc_int, col="blue", lwd=2, main="ROC Analysis: yyfbatch2 as Validation", legacy.axes=TRUE)
  plot(roc_ext, col="red", lwd=2, add=TRUE)
  
  legend("bottomright", 
         legend=c(paste0("Internal Valid (AUC=", round(auc(roc_int), 3), ")"),
                  paste0("External yyfbatch2 (AUC=", round(auc(roc_ext), 3), ")")),
         col=c("blue", "red"), lwd=2)
  
  # 输出混淆矩阵 (查看 yyfbatch2 的特异性)
  cat("\n=== yyfbatch2 Confusion Matrix ===\n")
  pred_raw <- predict(model, external_data[, c("Group", features)], type="raw")
  print(confusionMatrix(pred_raw, external_data$Group, positive="Tumor"))
}

# 执行绘图
plot_roc_comparison(final_model, internal_test_set, valid_final, final_miRNA_set)














library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 准备数据名单 (Swap: BRCA1 做验证，其他 6 个做训练)
# ==============================================================================
# 训练池：4个PRJNA + 2个yyfbatch (共6个)
pool_names <- c("PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049", "yyfbatch1", "yyfbatch2")

# 独立验证集：BRCA1 (TCGA)
valid_name <- "BRCA1"

# 读取数据
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- gsub("processed_results/|_processed.csv", "", files)

pool_list <- data_list[pool_names]
valid_list <- data_list[valid_name]

# ==============================================================================
# 2. 数据清洗与提取共有基因
# ==============================================================================
all_dfs <- c(pool_list, valid_list)
common_genes <- Reduce(intersect, lapply(all_dfs, function(df) colnames(df)[-1]))

cat("共有 miRNA 数量:", length(common_genes), "\n")

# 清洗函数
prep_data <- function(df, genes) {
  df_sub <- df[, c("Group", genes)]
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) {
    df_sub[, -1] <- log2(mat + 1)
  }
  return(df_sub)
}

pool_clean <- lapply(pool_list, prep_data, genes = common_genes)
valid_clean <- lapply(valid_list, prep_data, genes = common_genes)[[1]]

# 检查训练池中是否包含 Normal (防止 ComBat 报错)
total_normal <- sum(sapply(pool_clean, function(df) sum(df$Group == "Normal")))
if(total_normal < 2) stop("警告：你的 6 个训练数据集中 Normal 样本太少，无法进行分类训练！")

# ==============================================================================
# 3. 对 6 个训练数据集进行 ComBat 去批次
# ==============================================================================
cat("正在对 6 个小数据集进行合并和 ComBat...\n")

# A. 构建合并矩阵
expr_matrices <- lapply(pool_clean, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)

# B. 构建 Batch 向量
batch_vec <- unlist(lapply(names(pool_clean), function(n) rep(n, nrow(pool_clean[[n]]))))

# C. 构建分组信息
group_vec <- unlist(lapply(pool_clean, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# D. 运行 ComBat
# 注意：如果某个小数据集只有 Tumor 没有 Normal，这里可能会有 Warning，通常不影响运行
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# ==============================================================================
# 4. 数据重组与 Z-score 标准化
# ==============================================================================

# A. 处理训练池 (ComBat后)
# 修复列名问题 check.names = FALSE
combat_df_final <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                              t(combat_expr), 
                              check.names = FALSE)

# Scale (Z-score)
combat_scaled_mat <- scale(combat_df_final[, -1])
pool_final <- data.frame(Group = combat_df_final$Group, 
                         combat_scaled_mat, 
                         check.names = FALSE)

# B. 处理独立验证集 (BRCA1)
# BRCA1 不做 ComBat，只做 Z-score
valid_scaled_mat <- scale(valid_clean[, -1])
valid_final <- data.frame(Group = factor(valid_clean$Group, levels = c("Normal", "Tumor")), 
                          valid_scaled_mat, 
                          check.names = FALSE)

# ==============================================================================
# 5. 70% 训练 / 30% 内验 拆分
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(pool_final$Group, p = 0.7, list = FALSE, times = 1)

train_set <- pool_final[trainIndex, ]
internal_test_set <- pool_final[-trainIndex, ]

cat("训练集样本数:", nrow(train_set), "\n")
cat("独立验证集(BRCA1)样本数:", nrow(valid_final), "\n")

# ==============================================================================
# 6. 特征筛选 (在小样本集合上筛选)
# ==============================================================================
cat("开始特征筛选...\n")
x <- as.matrix(train_set[, -1])
y <- train_set$Group

# LASSO
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)
coefs <- coef(cv_lasso, s = "lambda.min")
lasso_feats <- rownames(coefs)[coefs[,1] != 0][-1]

# Random Forest
rf_fs <- randomForest(x, y, ntree=300)
rf_imp <- importance(rf_fs)
rf_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:20]

final_miRNA_set <- unique(c(lasso_feats, rf_feats))
# 限制特征数量
if(length(final_miRNA_set) > 15) final_miRNA_set <- final_miRNA_set[1:15]

cat("最终特征 Panel:", paste(final_miRNA_set, collapse=", "), "\n")

# ==============================================================================
# 7. 建模 (Random Forest)
# ==============================================================================
# 注意：现在的训练集(PRJNA+yyf)可能不像BRCA1那样有极端的10:1不平衡
# 我们可以先检查一下比例，决定是否需要 down-sampling
print(table(train_set$Group))

train_subset <- train_set[, c("Group", final_miRNA_set)]

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# 依然保持 down-sampling 是个保险的策略，或者改成 'cv' 不做采样
fitControl <- trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           sampling = "down") 

cat("开始训练模型...\n")
set.seed(888)
final_model <- train(Group ~ ., 
                     data = train_subset, 
                     method = "rf", 
                     metric = "ROC", 
                     trControl = fitControl,
                     ntree = 500)

stopCluster(cl)

# ==============================================================================
# 8. 最终验证 (在 BRCA1 上大考)
# ==============================================================================

# 定义绘图函数
plot_roc_comparison <- function(model, internal_data, external_data, features) {
  
  # 内部验证
  prob_int <- predict(model, internal_data[, c("Group", features)], type = "prob")
  roc_int <- roc(internal_data$Group, prob_int$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 外部验证 (BRCA1)
  prob_ext <- predict(model, external_data[, c("Group", features)], type = "prob")
  roc_ext <- roc(external_data$Group, prob_ext$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 绘图
  plot(roc_int, col="blue", lwd=2, main="Testing Model Generalization on BRCA1", legacy.axes=TRUE)
  plot(roc_ext, col="red", lwd=2, add=TRUE)
  
  legend("bottomright", 
         legend=c(paste0("Internal Valid (AUC=", round(auc(roc_int), 3), ")"),
                  paste0("External BRCA1 (AUC=", round(auc(roc_ext), 3), ")")),
         col=c("blue", "red"), lwd=2)
  
  # 输出混淆矩阵 (BRCA1)
  cat("\n=== BRCA1 Confusion Matrix ===\n")
  pred_raw <- predict(model, external_data[, c("Group", features)], type="raw")
  print(confusionMatrix(pred_raw, external_data$Group, positive="Tumor"))
}

# 执行绘图
plot_roc_comparison(final_model, internal_test_set, valid_final, final_miRNA_set)










library(sva)
library(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 准备数据名单 (Target: PRJNA808405 做验证)
# ==============================================================================
# 训练池：BRCA1 + 3个PRJNA + 2个yyfbatch (共6个)
pool_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA934049", "yyfbatch1", "yyfbatch2")

# 独立验证集：PRJNA808405
valid_name <- "PRJNA808405"

# 读取数据 (假设在 processed_results 目录下)
files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
data_list <- lapply(files, function(x) read.csv(x, row.names = 1, stringsAsFactors = FALSE))
names(data_list) <- gsub("processed_results/|_processed.csv", "", files)

# 提取列表
pool_list <- data_list[pool_names]
valid_list <- data_list[valid_name]

# ==============================================================================
# 2. 数据清洗与提取共有基因
# ==============================================================================
# 寻找所有 7 个数据集共有的基因
all_dfs <- c(pool_list, valid_list)
common_genes <- Reduce(intersect, lapply(all_dfs, function(df) colnames(df)[-1])) # -1 排除 Group

cat("共有 miRNA 数量:", length(common_genes), "\n")

# 清洗函数 (取共有基因 + Log2处理)
prep_data <- function(df, genes) {
  df_sub <- df[, c("Group", genes)]
  # 检查是否需要 Log2 (简单阈值判断)
  mat <- df_sub[, -1]
  if(max(mat, na.rm=TRUE) > 50) {
    df_sub[, -1] <- log2(mat + 1)
  }
  return(df_sub)
}

pool_clean <- lapply(pool_list, prep_data, genes = common_genes)
valid_clean <- lapply(valid_list, prep_data, genes = common_genes)[[1]]

# ==============================================================================
# 3. 对 6 个训练数据集进行 ComBat 去批次
# ==============================================================================
cat("正在对 6 个训练数据集进行 ComBat...\n")

# A. 构建合并矩阵 (ComBat 要求: 行=基因, 列=样本)
expr_matrices <- lapply(pool_clean, function(df) t(df[, -1]))
combined_expr <- do.call(cbind, expr_matrices)

# B. 构建 Batch 向量
batch_vec <- unlist(lapply(names(pool_clean), function(n) rep(n, nrow(pool_clean[[n]]))))

# C. 构建分组信息 (保护生物学差异)
group_vec <- unlist(lapply(pool_clean, function(df) df$Group))
mod <- model.matrix(~as.factor(group_vec))

# D. 运行 ComBat
combat_expr <- ComBat(dat = combined_expr, batch = batch_vec, mod = mod, par.prior = TRUE)

# ==============================================================================
# 4. 数据重组与 Z-score 标准化 (关键修正版)
# ==============================================================================

# A. 重组训练池数据 (ComBat后)
# 关键: 直接在 data.frame 构建时放入 Group，并使用 check.names=FALSE 保护 miRNA 名字
# 这样可以避免 "must be numeric" 报错和列名错乱
combat_df_final <- data.frame(Group = factor(group_vec, levels = c("Normal", "Tumor")), 
                              t(combat_expr), 
                              check.names = FALSE)

# Scale (Z-score) 训练池
# 注意：combat_df_final[,-1] 确保只取数值列进行计算
combat_scaled_mat <- scale(combat_df_final[, -1])
pool_final <- data.frame(Group = combat_df_final$Group, 
                         combat_scaled_mat, 
                         check.names = FALSE)

# B. 处理独立验证集 (PRJNA808405)
# 验证集不做 ComBat，但必须做 Z-score 以对齐量纲
valid_scaled_mat <- scale(valid_clean[, -1])
valid_final <- data.frame(Group = factor(valid_clean$Group, levels = c("Normal", "Tumor")), 
                          valid_scaled_mat, 
                          check.names = FALSE)

# ==============================================================================
# 5. 70% 训练 / 30% 内验 拆分
# ==============================================================================
set.seed(123)
# 基于 Group 进行分层抽样，保证训练集和内验集里 Tumor/Normal 比例一致
trainIndex <- createDataPartition(pool_final$Group, p = 0.7, list = FALSE, times = 1)

train_set <- pool_final[trainIndex, ]
internal_test_set <- pool_final[-trainIndex, ]

cat("训练集 (ComBat+Scaled):", nrow(train_set), "\n")
cat("内验集 (ComBat+Scaled):", nrow(internal_test_set), "\n")
cat("外验集 (PRJNA808405 Scaled):", nrow(valid_final), "\n")

# ==============================================================================
# 6. 特征筛选 (仅在 70% 训练集上)
# ==============================================================================
cat("开始特征筛选 (LASSO + RF)...\n")
x <- as.matrix(train_set[, -1])
y <- train_set$Group

# 1. LASSO
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)
coefs <- coef(cv_lasso, s = "lambda.min")
lasso_feats <- rownames(coefs)[coefs[,1] != 0][-1]

# 2. Random Forest Importance
rf_fs <- randomForest(x, y, ntree=300)
rf_imp <- importance(rf_fs)
rf_feats <- rownames(rf_imp)[order(rf_imp[,"MeanDecreaseGini"], decreasing = T)][1:20]

# 取并集
final_miRNA_set <- unique(c(lasso_feats, rf_feats))
# 限制特征数量 (防止过拟合)
if(length(final_miRNA_set) > 15) final_miRNA_set <- final_miRNA_set[1:15]

cat("最终特征 Panel:", paste(final_miRNA_set, collapse=", "), "\n")

# ==============================================================================
# 7. 建模 (Random Forest + Down-sampling)
# ==============================================================================
train_subset <- train_set[, c("Group", final_miRNA_set)]

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# 使用 Down-sampling 解决 BRCA1 (在训练池中) 带来的不平衡
fitControl <- trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           sampling = "down") 

cat("开始训练模型...\n")
set.seed(888)
final_model <- train(Group ~ ., 
                     data = train_subset, 
                     method = "rf", 
                     metric = "ROC", 
                     trControl = fitControl,
                     ntree = 500)

stopCluster(cl)

# ==============================================================================
# 8. 最终验证与画图
# ==============================================================================
# 定义绘图函数
plot_roc_comparison <- function(model, internal_data, external_data, features) {
  
  # 内部验证
  prob_int <- predict(model, internal_data[, c("Group", features)], type = "prob")
  roc_int <- roc(internal_data$Group, prob_int$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 外部验证 (PRJNA808405)
  prob_ext <- predict(model, external_data[, c("Group", features)], type = "prob")
  roc_ext <- roc(external_data$Group, prob_ext$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # 绘图
  plot(roc_int, col="blue", lwd=2, main="Validation on PRJNA808405", legacy.axes=TRUE)
  plot(roc_ext, col="red", lwd=2, add=TRUE)
  
  legend("bottomright", 
         legend=c(paste0("Internal Valid (AUC=", round(auc(roc_int), 3), ")"),
                  paste0("External PRJNA808405 (AUC=", round(auc(roc_ext), 3), ")")),
         col=c("blue", "red"), lwd=2)
  
  # 输出混淆矩阵 (查看外验特异性)
  cat("\n=== PRJNA808405 Confusion Matrix ===\n")
  pred_raw <- predict(model, external_data[, c("Group", features)], type="raw")
  print(confusionMatrix(pred_raw, external_data$Group, positive="Tumor"))
}

# 执行绘图
plot_roc_comparison(final_model, internal_test_set, valid_final, final_miRNA_set)


# ==============================================================================
# 修复步骤：处理 Scale 产生的 NA/NaN
# ==============================================================================

# 定义一个简单的修复函数
fix_na_inf <- function(df) {
  # 将所有 NA, NaN, Inf 替换为 0
  df[is.na(df)] <- 0
  df[sapply(df, is.infinite)] <- 0
  return(df)
}

# 1. 修复 内部验证集
internal_test_set <- fix_na_inf(internal_test_set)

# 2. 修复 外部验证集 (PRJNA808405)
valid_final <- fix_na_inf(valid_final)

cat("数据修复完成：已将 Z-score 产生的 NaN 替换为 0。\n")

# ==============================================================================
# 重新运行绘图函数
# ==============================================================================

plot_roc_comparison <- function(model, internal_data, external_data, features) {
  
  # 确保特征名一致性 (防止 "-" 被 R 自动变成 ".")
  # 这一步是为了防止 caret 训练时改了名字导致匹配不上
  model_features <- model$finalModel$xNames
  
  # 如果模型里的名字有变化 (比如变成了点号)，我们需要尝试匹配
  # 这里假设名字保持一致，如果有问题，通常是 caret 自动 sanitized 了
  
  # --- 内部验证 ---
  # 提取数据
  data_int <- internal_data[, c("Group", features)]
  # 再次确保无 NA (双重保险)
  data_int[is.na(data_int)] <- 0 
  
  prob_int <- predict(model, data_int, type = "prob")
  roc_int <- roc(internal_data$Group, prob_int$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # --- 外部验证 ---
  data_ext <- external_data[, c("Group", features)]
  # 再次确保无 NA
  data_ext[is.na(data_ext)] <- 0
  
  prob_ext <- predict(model, data_ext, type = "prob")
  roc_ext <- roc(external_data$Group, prob_ext$Tumor, levels=c("Normal", "Tumor"), direction="<")
  
  # --- 绘图 ---
  plot(roc_int, col="blue", lwd=2, main="Validation on PRJNA808405", legacy.axes=TRUE)
  plot(roc_ext, col="red", lwd=2, add=TRUE)
  
  legend("bottomright", 
         legend=c(paste0("Internal Valid (AUC=", round(auc(roc_int), 3), ")"),
                  paste0("External PRJNA808405 (AUC=", round(auc(roc_ext), 3), ")")),
         col=c("blue", "red"), lwd=2)
  
  cat("\n=== PRJNA808405 Confusion Matrix ===\n")
  pred_raw <- predict(model, data_ext, type="raw")
  print(confusionMatrix(pred_raw, external_data$Group, positive="Tumor"))
}

# 执行绘图
plot_roc_comparison(final_model, internal_test_set, valid_final, final_miRNA_set)