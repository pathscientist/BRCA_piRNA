## =============================
## 0) Packages
## =============================
# install.packages(c("randomForest", "pROC"))
library(randomForest)
library(pROC)

## =============================
## 1) Functions
## =============================

# expr: rows=gene, cols=sample
# pheno: rows=sample
reformat_dataset <- function(expr, pheno, group_col = "Group") {
  common_samples <- intersect(colnames(expr), rownames(pheno))
  expr  <- expr[, common_samples, drop = FALSE]
  pheno <- pheno[common_samples, , drop = FALSE]
  
  expr_t <- t(expr)
  
  df <- data.frame(
    group = pheno[[group_col]],
    expr_t,
    check.names = FALSE
  )
  rownames(df) <- common_samples
  df
}

# Benign->Normal, Cancer->Tumor, only Tumor/Normal allowed
recode_group <- function(df, group_col = "group") {
  g <- as.character(df[[group_col]])
  
  g[g == "Benign"] <- "Normal"
  g[g == "Cancer"] <- "Tumor"
  
  allowed <- c("Tumor", "Normal")
  bad <- setdiff(unique(g), allowed)
  if (length(bad) > 0) {
    stop("group 中仍有不是 Tumor/Normal 的取值：",
         paste(bad, collapse = ", "))
  }
  
  df[[group_col]] <- factor(g, levels = allowed)
  df
}

get_genes <- function(df) setdiff(colnames(df), "group")

# cap each cohort per class to avoid BRCA1 dominating
cap_cohort <- function(df, n_per_class = 120, seed = 1) {
  set.seed(seed)
  df$group <- factor(df$group, levels = c("Tumor","Normal"))
  idx <- unlist(lapply(levels(df$group), function(g) {
    ids <- which(df$group == g)
    if (length(ids) <= n_per_class) return(ids)
    sample(ids, n_per_class)
  }))
  df[idx, , drop = FALSE]
}

# align to a gene set
align_genes <- function(df, genes) {
  df[, c("group", genes), drop = FALSE]
}

# stratified 70/30 split on pooled internal
stratified_split <- function(df, train_ratio = 0.7, seed = 123) {
  set.seed(seed)
  df$group <- factor(df$group, levels = c("Tumor","Normal"))
  train_idx <- unlist(lapply(levels(df$group), function(g) {
    ids <- which(df$group == g)
    sample(ids, size = floor(train_ratio * length(ids)))
  }))
  list(
    train = df[train_idx, , drop = FALSE],
    test  = df[-train_idx, , drop = FALSE]
  )
}

# train RF with balanced sampling
fit_rf_balanced <- function(train_df, ntree = 1000, seed = 123) {
  set.seed(seed)
  genes <- get_genes(train_df)
  x <- train_df[, genes, drop = FALSE]
  y <- train_df$group
  
  # balanced per tree
  min_n <- min(table(y))
  
  randomForest(
    x = x, y = y,
    ntree = ntree,
    strata = y,
    sampsize = rep(min_n, length(levels(y)))
  )
}

# predict AUC
auc_on <- function(model, df) {
  genes <- intersect(colnames(df), colnames(model$importance))
  # 上面这一句可能太苛刻，改成以训练时用的基因为准更稳：
  # 但 randomForest对象不直接存gene列表，这里用 df 的gene列即可
  genes <- get_genes(df)
  
  prob <- predict(model, df[, genes, drop=FALSE], type = "prob")[, "Tumor"]
  roc_obj <- roc(df$group, prob, levels = c("Normal","Tumor"), direction = "<", quiet = TRUE)
  as.numeric(auc(roc_obj))
}

## =============================
## 2) Your 7 datasets
## =============================

expr_list <- list(
  BRCA1       = BRCA1,
  PRJNA294226 = PRJNA294226,
  PRJNA482141 = PRJNA482141,
  PRJNA808405 = PRJNA808405,
  PRJNA934049 = PRJNA934049,
  yyfbatch1   = yyfbatch1,
  yyfbatch2   = yyfbatch2
)

pheno_list <- list(
  BRCA1       = BRCA1_pheno,
  PRJNA294226 = PRJNA294226_pheno,
  PRJNA482141 = PRJNA482141_pheno,
  PRJNA808405 = PRJNA808405_pheno,
  PRJNA934049 = PRJNA934049_pheno,
  yyfbatch1   = yyfbatch1_pheno,
  yyfbatch2   = yyfbatch2_pheno
)

## =============================
## 3) Reformat + recode
## =============================

ready_list <- list()
for (nm in names(expr_list)) {
  df <- reformat_dataset(expr_list[[nm]], pheno_list[[nm]], group_col = "Group")
  df <- recode_group(df, "group")
  ready_list[[nm]] <- df
}

## =============================
## 4) Define internal/external
## =============================

external_name <- "yyfbatch1"
internal_names <- c("BRCA1","PRJNA294226","PRJNA482141",
                    "PRJNA808405","PRJNA934049","yyfbatch2")

external_val_raw <- ready_list[[external_name]]
internal_list_raw <- ready_list[internal_names]

## =============================
## 5) Per-cohort cap (recommended)
## =============================

internal_list_cap <- lapply(internal_list_raw, cap_cohort, n_per_class = 120, seed = 1)

## =============================
## 6) Use common genes across INTERNAL only
##     (avoid letting external influence feature space)
## =============================

common_genes_internal <- Reduce(intersect, lapply(internal_list_cap, get_genes))

internal_list_cap <- lapply(internal_list_cap, align_genes, genes = common_genes_internal)
external_val <- align_genes(external_val_raw, genes = common_genes_internal)

## =============================
## 7) LOCO internal evaluation
## =============================

loco_auc <- c()

for (test_cohort in internal_names) {
  train_cohorts <- setdiff(internal_names, test_cohort)
  
  train_df <- do.call(rbind, internal_list_cap[train_cohorts])
  test_df  <- internal_list_cap[[test_cohort]]
  
  # Fit balanced RF
  rf <- fit_rf_balanced(train_df, ntree = 1000, seed = 123)
  
  # AUC on left-out cohort
  prob <- predict(rf, test_df[, common_genes_internal, drop=FALSE], type="prob")[, "Tumor"]
  roc_obj <- roc(test_df$group, prob, levels = c("Normal","Tumor"), direction = "<", quiet = TRUE)
  
  loco_auc[test_cohort] <- as.numeric(auc(roc_obj))
}

print(loco_auc)
cat("LOCO mean AUC:", mean(loco_auc), "\n\n")

## =============================
## 8) Pooled internal 70/30 AFTER cap
## =============================

internal_all <- do.call(rbind, internal_list_cap)
spl <- stratified_split(internal_all, train_ratio = 0.7, seed = 123)

internal_train <- spl$train
internal_test  <- spl$test

rf_final <- fit_rf_balanced(internal_train, ntree = 1000, seed = 123)

# Internal pooled AUC
prob_int <- predict(rf_final, internal_test[, common_genes_internal, drop=FALSE], type="prob")[, "Tumor"]
roc_int  <- roc(internal_test$group, prob_int, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)
auc_int  <- as.numeric(auc(roc_int))

# External AUC
prob_ext <- predict(rf_final, external_val[, common_genes_internal, drop=FALSE], type="prob")[, "Tumor"]
roc_ext  <- roc(external_val$group, prob_ext, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)
auc_ext  <- as.numeric(auc(roc_ext))

cat("Internal pooled 70/30 AUC:", auc_int, "\n")
cat("External (yyfbatch1) AUC:", auc_ext, "\n\n")

## =============================
## 9) Optional: save splits
## =============================

write.csv(internal_train, "internal_train_cap70.csv", row.names = TRUE)
write.csv(internal_test,  "internal_test_cap30.csv",  row.names = TRUE)
write.csv(external_val,   "yyfbatch1_external_validation_aligned.csv", row.names = TRUE)



library(sva)

# 1. 分离数据
training_datasets <- list(BRCA1, PRJNA294226, PRJNA482141, PRJNA808405, PRJNA934049, yyfbatch2)
validation_dataset <- yyfbatch1 # 把它踢出 ComBat 的圈子

# 2. 仅在训练集内部做 ComBat
# (先合并这6个数据，建立 batch 向量)
train_combined <- do.call(cbind, lapply(training_datasets, t)) # 注意行列转换
batch_vec <- c(...) # 对应这6个数据集的批次
mod <- model.matrix(~as.factor(Group), data=pheno_train)

# 只对训练集去批次！
train_combat <- ComBat(dat=train_combined, batch=batch_vec, mod=mod)

# 3. 验证集保持原样 (只做 log2)
valid_data <- log2(validation_dataset + 1)

# 4. 训练与预测...


# 定义一个函数：对每个数据集单独做 Z-score
# (x - mean) / sd
scale_per_dataset <- function(df) {
  # 假设第一列是 Group，后面是基因
  genes <- df[, -1]
  # 对每一列(基因)进行 scale
  genes_scaled <- scale(genes) 
  
  # 重新组合
  return(data.frame(Group = df$Group, genes_scaled))
}

# 1. 对 list 中的每一个数据集，单独执行标准化
# 注意：一定要在合并之前做！
processed_list_scaled <- lapply(processed_list, scale_per_dataset)

# 2. 之后再进行合并训练集、提取验证集的操作
# ... (后续训练流程不变)


# 在 trainControl 中加入 sampling = "down"
fitControl <- trainControl(method = "cv",
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           sampling = "down") # <--- 关键就是这一句！

# 然后正常训练
set.seed(888)
model_rf_balanced <- train(Group ~ ., data = train_subset, 
                           method = "rf", 
                           metric = "ROC", 
                           trControl = fitControl)


# 需要先安装额外的包，caret 会自动调用
# install.packages("DMwR") # 旧版 R
# 或者 install.packages("themis") 

fitControl <- trainControl(method = "cv",
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           sampling = "smote") # <--- 改成 smote

# 正常训练...

# 计算权重：Normal 很少，所以权重设大一点
model_rf_weighted <- train(Group ~ ., data = train_subset, 
                           method = "rf", 
                           metric = "ROC", 
                           trControl = fitControl,
                           # classwt 是 randomForest 特有的参数
                           # 假设 Normal 是第一个 level，Tumor 是第二个
                           classwt = c(10, 1))



library(caret)
library(dplyr)
library(pROC)
library(doParallel)

# ==============================================================================
# 1. 解决批次效应：单数据集独立标准化 (Z-score)
# ==============================================================================
# 原理：不强行把所有数据拉齐，而是把每个数据集都变成标准正态分布。
# 这样 BRCA1 的 "高表达" 和 yyfbatch1 的 "高表达" 在数值上就一致了。

# 定义标准化函数
normalize_per_dataset <- function(df) {
  # 假设第一列是 Group，后面是基因
  # 1. 提取数值矩阵
  mat <- df[, -1] 
  
  # 2. 对每个基因(列)进行 (x - mean)/sd 处理
  # scale 函数默认就是做 Z-score
  mat_scaled <- scale(mat)
  
  # 3. 重新拼合
  df_final <- data.frame(Group = df$Group, mat_scaled)
  return(df_final)
}

# 应用到你的 list 中 (假设 processed_list 包含你所有的7个原始数据)
# 这一步非常关键！一定要在合并之前做！
processed_list_norm <- lapply(processed_list, normalize_per_dataset)

# ==============================================================================
# 2. 数据划分：构建训练池和独立验证池
# ==============================================================================

# 定义名单
train_pool_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "PRJNA934049", "yyfbatch2")
valid_pool_name <- "yyfbatch1" # 你的独立验证集

# 提取并寻找共有基因
train_list <- processed_list_norm[train_pool_names]
valid_list <- processed_list_norm[valid_pool_name]

common_genes <- Reduce(intersect, lapply(c(train_list, valid_list), function(x) colnames(x)[-1]))
cat("共有基因数量:", length(common_genes), "\n")

# 合并训练集
train_data_all <- do.call(rbind, lapply(train_list, function(x) x[, c("Group", common_genes)]))
# 提取验证集
validation_data <- valid_list[[1]][, c("Group", common_genes)]

# 确保 Group 格式正确 (Normal在前, Tumor在后)
train_data_all$Group <- factor(train_data_all$Group, levels = c("Normal", "Tumor"))
validation_data$Group <- factor(validation_data$Group, levels = c("Normal", "Tumor"))

# ==============================================================================
# 3. 特征筛选 (仅在训练集上做，防止泄露)
# ==============================================================================
# 这里省略具体的 Feature Selection 代码，假设你已经跑完 LASSO/Boruta
# 并得到了最终的特征列表: final_miRNA_set
# 示例：
# final_miRNA_set <- c("hsa-miR-21-5p", "hsa-miR-155-5p", ...) 假设这是选出来的

# ！！！重要：为了演示，这里假设你已经有了 final_miRNA_set ！！！
# 请务必确保特征筛选时没有用到 validation_data

# ==============================================================================
# 4. 解决样本不平衡：Down-sampling 建模
# ==============================================================================
# 准备最终训练数据
train_subset <- train_data_all[, c("Group", final_miRNA_set)]

# 开启并行
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# *** 关键配置 ***
fitControl <- trainControl(method = "cv",
                           number = 10,       # 10折交叉验证
                           classProbs = TRUE, # 需要概率来算 AUC
                           summaryFunction = twoClassSummary,
                           
                           # >>> 核心修改：下采样 <<<
                           # 每次 CV 训练时，caret 会随机扔掉多余的 Tumor，
                           # 保持 Tumor 和 Normal 数量 1:1。
                           sampling = "down") 

cat("开始训练模型 (使用 Down-sampling 处理不平衡)...\n")

set.seed(123)
# 推荐使用 Random Forest (rf) 或 Elastic Net (glmnet)
final_model <- train(Group ~ ., 
                     data = train_subset, 
                     method = "rf",  # 随机森林
                     metric = "ROC", # 以优化 ROC 为目标
                     trControl = fitControl,
                     ntree = 500)    # 树的数量

print(final_model)

# ==============================================================================
# 5. 最终大考：独立验证
# ==============================================================================
valid_subset <- validation_data[, c("Group", final_miRNA_set)]

# 预测概率
pred_prob <- predict(final_model, newdata = valid_subset, type = "prob")
# 预测类别
pred_raw <- predict(final_model, newdata = valid_subset, type = "raw")

# 计算 AUC
roc_valid <- roc(validation_data$Group, pred_prob$Tumor, levels = c("Normal", "Tumor"), direction = "<")
auc_val <- auc(roc_valid)

cat("\n=========================================\n")
cat("独立验证集 (yyfbatch1) 最终 AUC:", auc_val, "\n")
cat("=========================================\n")

# 查看混淆矩阵 (重点看 Specificity 是否提高了)
print(confusionMatrix(pred_raw, validation_data$Group, positive = "Tumor"))

stopCluster(cl)




# ==============================================================================
# 1. 策略一：双重验证集划分
# ==============================================================================

# 训练集名单 (TCGA + 4个PRJNA + 1个私有)
train_names <- c("BRCA1", "PRJNA294226", "PRJNA482141", "PRJNA808405", "yyfbatch2") 
# 注：留下了 PRJNA934049 和 yyfbatch1 做验证

# 验证集 1 (公共)
valid_public_name <- "PRJNA934049"

# 验证集 2 (私有)
valid_private_name <- "yyfbatch1"

# ==============================================================================
# 2. 数据处理 (核心：独立 Z-score + 下采样)
# ==============================================================================
# 假设 processed_list 是你所有的原始数据

# A. 独立标准化 (Z-score) - 必须做！
normalized_list <- lapply(processed_list, function(df) {
  df[, -1] <- scale(df[, -1]) # 对基因列做 z-score
  return(df)
})

# B. 组装数据
# 提取共有基因
common_genes <- Reduce(intersect, lapply(normalized_list, function(x) colnames(x)[-1]))

# 训练集
train_data <- do.call(rbind, normalized_list[train_names])
train_data <- train_data[, c("Group", common_genes)]

# 验证集 1
valid_public <- normalized_list[[valid_public_name]][, c("Group", common_genes)]

# 验证集 2
valid_private <- normalized_list[[valid_private_name]][, c("Group", common_genes)]

# ==============================================================================
# 3. 训练模型 (带 Down-sampling)
# ==============================================================================
library(caret)
train_subset <- train_data[, c("Group", final_miRNA_set)] # 假设你已经选好特征

fitControl <- trainControl(method = "cv", number = 10, 
                           classProbs = TRUE, summaryFunction = twoClassSummary,
                           sampling = "down") # 解决 TCGA 不平衡

set.seed(123)
final_model <- train(Group ~ ., data = train_subset, method = "rf", 
                     metric = "ROC", trControl = fitControl)

# ==============================================================================
# 4. 双重验证与绘图
# ==============================================================================
library(pROC)

# 预测公共验证集
prob_pub <- predict(final_model, valid_public, type = "prob")
roc_pub <- roc(valid_public$Group, prob_pub$Tumor, levels=c("Normal","Tumor"), direction="<")

# 预测私有验证集
prob_pri <- predict(final_model, valid_private, type = "prob")
roc_pri <- roc(valid_private$Group, prob_pri$Tumor, levels=c("Normal","Tumor"), direction="<")

# 画图：一张图两条线
plot(roc_pub, col="blue", lwd=2, main="Double External Validation")
plot(roc_pri, col="red", lwd=2, add=TRUE)
legend("bottomright", 
       legend=c(paste0("Public Valid (AUC=", round(auc(roc_pub),3), ")"),
                paste0("Private Valid (AUC=", round(auc(roc_pri),3), ")")),
       col=c("blue", "red"), lwd=2)