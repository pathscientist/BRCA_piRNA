setwd("C:/HP/PhD/Projects/GEO BRCA/brca_dataset_1101tpm")

BRCA1 <- read.table("BRCA1.txt", header = TRUE, row.names = 1, sep = "\t")
#plasma <- read.table("plasma.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA294226 <- read.table("PRJNA294226.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA482141 <- read.table("PRJNA482141.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA808405 <- read.table("PRJNA808405.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA934049 <- read.table("PRJNA934049.txt", header = TRUE, row.names = 1, sep = "\t")
#tissue <- read.table("tissue.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch1 <- read.table("yyfbatch1.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch2 <- read.table("yyfbatch2.txt", header = TRUE, row.names = 1, sep = "\t")

#input the pheno data
BRCA1_pheno <- read.table("BRCA1_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA294226_pheno <- read.table("PRJNA294226_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA482141_pheno <- read.table("PRJNA482141_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA808405_pheno <- read.table("PRJNA808405_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA934049_pheno <- read.table("PRJNA934049_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch1_pheno <- read.table("yyfbatch1_pheno.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch2_pheno <- read.table("yyfbatch2_pheno.txt", header = TRUE, row.names = 1, sep = "\t")

#########################

# 1. 准备数据列表 (和之前一样)
expr_list <- list(
  BRCA1 = BRCA1, 
  PRJNA294226 = PRJNA294226, 
  PRJNA482141 = PRJNA482141,
  PRJNA808405 = PRJNA808405,
  PRJNA934049 = PRJNA934049,
  yyfbatch1 = yyfbatch1,
  yyfbatch2 = yyfbatch2
)

pheno_list <- list(
  BRCA1 = BRCA1_pheno,
  PRJNA294226 = PRJNA294226_pheno,
  PRJNA482141 = PRJNA482141_pheno,
  PRJNA808405 = PRJNA808405_pheno,
  PRJNA934049 = PRJNA934049_pheno,
  yyfbatch1 = yyfbatch1_pheno,
  yyfbatch2 = yyfbatch2_pheno
)

# 2. 定义包含清洗和筛选功能的处理函数
process_and_clean <- function(expr_data, pheno_data, dataset_name) {
  
  # --- A. 转置 ---
  expr_t <- t(expr_data)
  expr_df <- as.data.frame(expr_t)
  
  # --- B. 合并 Group ---
  # 取交集样本，防止报错
  common_samples <- intersect(rownames(expr_df), rownames(pheno_data))
  
  if(length(common_samples) == 0) {
    warning(paste("警告:", dataset_name, "没有匹配的样本名，请检查行名格式！"))
    return(NULL)
  }
  
  expr_df <- expr_df[common_samples, , drop = FALSE]
  pheno_subset <- pheno_data[common_samples, , drop = FALSE]
  
  # 赋予初始 Group
  expr_df$Group <- pheno_subset$Group
  
  # --- C. 值替换 (Benign->Normal, Cancer->Tumor) ---
  # 使用 gsub 做字符替换，或者直接全等匹配替换
  expr_df$Group[expr_df$Group == "Benign"] <- "Normal"
  expr_df$Group[expr_df$Group == "Cancer"] <- "Tumor"
  
  # --- D. 筛选 (只保留 Tumor 和 Normal) ---
  # 找出要保留的行
  keep_rows <- expr_df$Group %in% c("Tumor", "Normal")
  
  # 打印过滤信息，让你知道删了多少数据
  removed_count <- sum(!keep_rows)
  if(removed_count > 0) {
    cat(paste0("[", dataset_name, "] 过滤掉了 ", removed_count, " 个非 Tumor/Normal 的样本\n"))
  }
  
  expr_df <- expr_df[keep_rows, , drop = FALSE]
  
  # --- E. 整理列顺序 (Group 在第一列) ---
  genes <- setdiff(colnames(expr_df), "Group")
  expr_df <- expr_df[, c("Group", genes)]
  
  return(expr_df)
}

# 3. 批量处理
# 使用 Map 并传递名字，方便打印日志
processed_list <- Map(process_and_clean, expr_list, pheno_list, names(expr_list))

# 4. 批量保存为 CSV
# 创建一个文件夹存放结果，保持目录整洁
if(!dir.exists("processed_results")) dir.create("processed_results")

for (name in names(processed_list)) {
  data_to_save <- processed_list[[name]]
  
  if (!is.null(data_to_save) && nrow(data_to_save) > 0) {
    # 构建文件名
    file_name <- paste0("processed_results/", name, "_processed.csv")
    
    # 保存 CSV
    # row.names = TRUE 会把样本名作为第一列保留
    write.csv(data_to_save, file = file_name, row.names = TRUE)
    
    cat(paste("成功保存:", file_name, "\n"))
  }
}


#################################

reformat_dataset <- function(expr, pheno, group_col = "Group") {
  ## expr: 行=gene, 列=sample
  ## pheno: 行=sample, 至少包含一列group信息（比如 "Group"）
  
  # 1. 只保留两边都有的样本，且顺序一致
  common_samples <- intersect(colnames(expr), rownames(pheno))
  expr  <- expr[, common_samples, drop = FALSE]
  pheno <- pheno[common_samples, , drop = FALSE]
  
  # 2. 转置表达矩阵：现在行 = sample, 列 = gene
  expr_t <- t(expr)
  
  # 3. 变成 data.frame，并把 group 放到第一列
  df <- data.frame(
    group = pheno[[group_col]],
    expr_t,
    check.names = FALSE  # 不要改基因名
  )
  rownames(df) <- common_samples
  
  return(df)
}

BRCA1_ready      <- reformat_dataset(BRCA1,      BRCA1_pheno,      group_col = "Group")
PRJNA294226_ready<- reformat_dataset(PRJNA294226,PRJNA294226_pheno,group_col = "Group")
PRJNA482141_ready<- reformat_dataset(PRJNA482141,PRJNA482141_pheno,group_col = "Group")
PRJNA808405_ready<- reformat_dataset(PRJNA808405,PRJNA808405_pheno,group_col = "Group")
PRJNA934049_ready<- reformat_dataset(PRJNA934049,PRJNA934049_pheno,group_col = "Group")
yyfbatch1_ready  <- reformat_dataset(yyfbatch1,  yyfbatch1_pheno,  group_col = "Group")
yyfbatch2_ready  <- reformat_dataset(yyfbatch2,  yyfbatch2_pheno,  group_col = "Group")


dat_list <- list(
  BRCA1_ready       = BRCA1_ready,
  PRJNA294226_ready = PRJNA294226_ready,
  PRJNA482141_ready = PRJNA482141_ready,
  PRJNA808405_ready = PRJNA808405_ready,
  PRJNA934049_ready = PRJNA934049_ready,
  yyfbatch1_ready   = yyfbatch1_ready,
  yyfbatch2_ready   = yyfbatch2_ready
)

for (nm in names(dat_list)) {
  write.csv(dat_list[[nm]],
            file = paste0(nm, ".csv"),
            row.names = TRUE)
}


recode_group <- function(df, group_col = "group") {
  g <- df[[group_col]]
  
  # 文本替换
  g[g == "Benign"] <- "Normal"
  g[g == "Cancer"] <- "Tumor"
  
  # 只允许 Tumor / Normal
  allowed <- c("Tumor", "Normal")
  if (!all(unique(g) %in% allowed)) {
    stop("group 中仍有不是 Tumor/Normal 的取值：",
         paste(setdiff(unique(g), allowed), collapse = ", "))
  }
  
  # 重新写回去，并可选转成因子
  df[[group_col]] <- factor(g, levels = allowed)
  df
}

BRCA1_ready       <- recode_group(BRCA1_ready)
PRJNA294226_ready <- recode_group(PRJNA294226_ready)
PRJNA482141_ready <- recode_group(PRJNA482141_ready)
PRJNA808405_ready <- recode_group(PRJNA808405_ready)
PRJNA934049_ready <- recode_group(PRJNA934049_ready)
yyfbatch1_ready   <- recode_group(yyfbatch1_ready)
yyfbatch2_ready   <- recode_group(yyfbatch2_ready)

write.csv(BRCA1_ready,       "BRCA1_ready.csv",       row.names = TRUE)
write.csv(PRJNA294226_ready, "PRJNA294226_ready.csv", row.names = TRUE)
write.csv(PRJNA482141_ready, "PRJNA482141_ready.csv", row.names = TRUE)
write.csv(PRJNA808405_ready, "PRJNA808405_ready.csv", row.names = TRUE)
write.csv(PRJNA934049_ready, "PRJNA934049_ready.csv", row.names = TRUE)
write.csv(yyfbatch1_ready,   "yyfbatch1_ready.csv",   row.names = TRUE)
write.csv(yyfbatch2_ready,   "yyfbatch2_ready.csv",   row.names = TRUE)

#########################################################################

library(caret)
library(glmnet)
library(Boruta)
library(randomForest)
library(xgboost)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 数据准备与分类
# ==============================================================================

# 假设 processed_list 包含了所有处理好的数据 (第一列 Group, 后面是基因)
# 确保所有数据的列名(Gene)是一致的

# 定义名单
pool_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049", "yyfbatch2")
external_name <- "yyfbatch1"

# 提取数据对象
pool_list <- processed_list[pool_names]
external_list <- processed_list[external_name]

# --- 寻找共有基因 ---
# 必须保证训练集和外部验证集有相同的基因
all_dfs <- c(pool_list, external_list)
common_genes <- Reduce(intersect, lapply(all_dfs, function(df) colnames(df)[-1]))

cat("所有 7 个数据集共有的 miRNA 数量:", length(common_genes), "\n")

# --- 数据标准化辅助函数 ---
clean_and_log <- function(df, genes) {
  # 提取共有基因
  sub_df <- df[, c("Group", genes)]
  # 简单的 log2 转换 (如果数据还没取log)
  mat <- sub_df[, -1]
  if(max(mat, na.rm=T) > 50) {
    sub_df[, -1] <- log2(mat + 1)
  }
  return(sub_df)
}

# 应用清洗
pool_clean_list <- lapply(pool_list, clean_and_log, genes = common_genes)
external_clean_list <- lapply(external_list, clean_and_log, genes = common_genes)

# --- 合并大数据池 (Pool A) ---
big_pool_data <- do.call(rbind, pool_clean_list)
# 准备外部验证集 (Pool B)
external_validation_data <- external_clean_list[[1]]

# 确保 Group 是因子，且 Tumor 是 Positive
big_pool_data$Group <- factor(big_pool_data$Group, levels = c("Normal", "Tumor"))
external_validation_data$Group <- factor(external_validation_data$Group, levels = c("Normal", "Tumor"))

# ==============================================================================
# 2. 70% 训练 / 30% 内部验证 拆分
# ==============================================================================
set.seed(123) # 保证结果可复现

# 在 big_pool_data 内部进行拆分
trainIndex <- createDataPartition(big_pool_data$Group, p = 0.7, list = FALSE, times = 1)

train_set <- big_pool_data[trainIndex, ]       # 真正的训练集 (用于特征筛选 + 建模)
internal_test_set <- big_pool_data[-trainIndex, ] # 内部验证集 (用于初步评估)

cat("Training Set (70%):", nrow(train_set), "samples\n")
cat("Internal Test Set (30%):", nrow(internal_test_set), "samples\n")
cat("External Validation (yyfbatch1):", nrow(external_validation_data), "samples\n")

# 开启多核并行
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# ==============================================================================
# 3. 特征筛选 (Feature Selection) - 仅在 train_set 上操作！
# ==============================================================================
cat(">>> 开始在 Training Set 上进行特征筛选...\n")

x_train <- as.matrix(train_set[, -1])
y_train <- train_set$Group

features_list <- list()

# 1. LASSO
cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
features_list$LASSO <- rownames(coef(cv_lasso, s = "lambda.min"))[coef(cv_lasso, s = "lambda.min")[,1]!= 0][-1]

# 2. Random Forest Importance
rf_fs <- randomForest(x_train, y_train, ntree = 300)
imp_rf <- importance(rf_fs)
features_list$RF <- rownames(imp_rf)[order(imp_rf[,"MeanDecreaseGini"], decreasing = T)][1:50]

# 3. Wilcoxon Test (P < 0.0001)
pvals <- apply(x_train, 2, function(col) wilcox.test(col ~ y_train)$p.value)
features_list$Wilcoxon <- names(pvals)[pvals < 1e-4]

# 4. XGBoost Gain
dtrain_fs <- xgb.DMatrix(data = x_train, label = ifelse(y_train=="Tumor", 1, 0))
xgb_fs <- xgb.train(params = list(objective="binary:logistic"), data=dtrain_fs, nrounds=50, verbose=0)
xgb_imp <- xgb.importance(model = xgb_fs)
features_list$XGBoost <- xgb_imp$Feature[1:min(50, nrow(xgb_imp))]

# --- 投票选出最佳 miRNA ---
all_votes <- unlist(features_list)
vote_counts <- sort(table(all_votes), decreasing = T)

# 阈值：至少被 3 种方法选中
final_miRNA_set <- names(vote_counts)[vote_counts >= 3]

# 如果太少，强制取前 5 个
if(length(final_miRNA_set) < 5) final_miRNA_set <- names(vote_counts)[1:5]

cat("最终筛选出的 miRNA (Panel):", length(final_miRNA_set), "个\n")
print(final_miRNA_set)

# ==============================================================================
# 4. 模型训练 (Model Training) - 基于 train_set
# ==============================================================================
cat(">>> 开始训练诊断模型...\n")

# 只保留筛选出的特征
train_subset <- train_set[, c("Group", final_miRNA_set)]

fitControl <- trainControl(method = "cv", number = 5, # 5折交叉验证
                           classProbs = TRUE, summaryFunction = twoClassSummary)

model_list <- list()

# 训练 RF, SVM, GLM (Logistic), XGB
set.seed(888)
model_list$RF <- train(Group ~ ., data = train_subset, method = "rf", metric = "ROC", trControl = fitControl)
set.seed(888)
model_list$SVM <- train(Group ~ ., data = train_subset, method = "svmRadial", metric = "ROC", trControl = fitControl)
set.seed(888)
model_list$GLM <- train(Group ~ ., data = train_subset, method = "glmnet", metric = "ROC", trControl = fitControl)
set.seed(888)
model_list$XGB <- train(Group ~ ., data = train_subset, method = "xgbTree", metric = "ROC", trControl = fitControl, verbosity=0)

# 选择最佳模型
res <- resamples(model_list)
best_model_name <- names(which.max(summary(res)$statistics$ROC[,"Mean"]))
final_model <- model_list[[best_model_name]]

cat("最佳模型是:", best_model_name, "\n")

# ==============================================================================
# 5. 验证阶段 (Validation Phase)
# ==============================================================================

# 定义一个绘图和评估的函数
evaluate_performance <- function(model, data, title_text) {
  # 提取特征列
  data_sub <- data[, c("Group", final_miRNA_set)]
  
  # 预测概率
  pred_prob <- predict(model, newdata = data_sub, type = "prob")
  # 预测类别
  pred_raw <- predict(model, newdata = data_sub, type = "raw")
  
  # 计算 ROC
  roc_obj <- roc(response = data_sub$Group, predictor = pred_prob$Tumor, 
                 levels = c("Normal", "Tumor"), direction = "<")
  auc_val <- auc(roc_obj)
  
  # 打印混淆矩阵
  cat("\n------------------------------------------------\n")
  cat("验证集: ", title_text, "\n")
  cat("AUC: ", auc_val, "\n")
  print(confusionMatrix(pred_raw, data_sub$Group, positive = "Tumor"))
  
  # 返回 ROC 对象用于绘图
  return(list(roc = roc_obj, auc = auc_val))
}

# --- A. 内部验证 (30% Split) ---
res_internal <- evaluate_performance(final_model, internal_test_set, "内部测试集 (30% Split)")

# --- B. 外部独立验证 (yyfbatch1) ---
res_external <- evaluate_performance(final_model, external_validation_data, "外部独立验证集 (yyfbatch1)")

# ==============================================================================
# 6. 绘图 (将两个 ROC 画在一起)
# ==============================================================================
plot(res_internal$roc, col = "blue", lwd = 2, main = "Diagnosis Model ROC Curves")
plot(res_external$roc, col = "red", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(paste0("Internal Valid (AUC=", round(res_internal$auc, 3), ")"),
                  paste0("External yyfbatch1 (AUC=", round(res_external$auc, 3), ")")),
       col = c("blue", "red"), lwd = 2)

stopCluster(cl)

