# Advanced Machine Learning Pipeline for piRNA-based Cancer Classification
# with Beautiful ggplot2 Visualizations and Model Comparison

# Load required libraries
library(randomForest)
library(caret)
library(ggplot2)
library(dplyr)
library(pROC)
library(e1071)        # For SVM
library(glmnet)       # For Elastic Net
library(xgboost)      # For XGBoost
library(gridExtra)    # For arranging plots
library(RColorBrewer) # For color palettes
library(viridis)      # For color palettes
library(plotly)       # For interactive plots
library(corrplot)     # For correlation plots
library(reshape2)     # For data reshaping
library(scales)       # For pretty breaks

# Set ggplot2 theme
theme_set(theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "white"),
    strip.text = element_text(face = "bold")
  ))

# Custom color palette
colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928")

# Assume your data frame is named 'pirna_data'
# If you need to load data, uncomment one of these:
# pirna_data <- read.csv("your_data_file.csv", row.names = 1)
# pirna_data <- read.table("your_data_file.txt", header = TRUE, row.names = 1)

# Data preprocessing and exploration
print("Data dimensions:")
print(dim(pirna_data))

# Convert label to factor
pirna_data$label <- as.factor(pirna_data$label)

# Simplified data overview (removing the problematic section)
cat("Data overview:\n")
cat("- Number of samples:", nrow(pirna_data), "\n")
cat("- Number of piRNAs:", ncol(pirna_data) - 1, "\n")  # Subtract 1 for label column
cat("- Label column:", "label", "\n")
cat("- piRNA columns:", paste(colnames(pirna_data)[2:4], "...", colnames(pirna_data)[ncol(pirna_data)]), "\n")

# Label distribution plot
label_dist <- pirna_data %>%
  count(label) %>%
  mutate(percentage = n / sum(n) * 100)

p_label_dist <- ggplot(label_dist, aes(x = label, y = n, fill = label)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = colors[1:2]) +
  labs(title = "Sample Distribution by Label",
       subtitle = "Number and percentage of samples in each class",
       x = "Sample Type", y = "Count") +
  theme(legend.position = "none") +
  ylim(0, max(label_dist$n) * 1.2)

print(p_label_dist)

# Remove missing values if present
if(sum(is.na(pirna_data)) > 0) {
  pirna_data <- na.omit(pirna_data)
  cat("Removed rows with missing values\n")
}

# Split data into training and testing sets
set.seed(123)
train_indices <- createDataPartition(pirna_data$label, p = 0.8, list = FALSE)
train_data <- pirna_data[train_indices, ]
test_data <- pirna_data[-train_indices, ]

cat("Training set size:", nrow(train_data), "\n")
cat("Test set size:", nrow(test_data), "\n")

# Prepare data for modeling (separate features and labels)
train_x <- train_data[, -1]  # Remove label column
train_y <- train_data$label
test_x <- test_data[, -1]
test_y <- test_data$label

# Define cross-validation control
cv_control <- trainControl(
  method = "cv",
  number = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = TRUE,
  verboseIter = FALSE
)

# Initialize results storage
model_results <- list()
all_predictions <- data.frame()
roc_data <- list()

# ===== MODEL 1: RANDOM FOREST =====
cat("Training Random Forest...\n")
set.seed(456)
rf_model <- train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = cv_control,
  tuneGrid = expand.grid(mtry = c(2, 4, 6, 8, sqrt(ncol(train_x)))),
  ntree = 500,
  importance = TRUE,
  metric = "ROC"
)

rf_pred <- predict(rf_model, test_x)
rf_prob <- predict(rf_model, test_x, type = "prob")
model_results$RandomForest <- confusionMatrix(rf_pred, test_y)
roc_data$RandomForest <- roc(test_y, rf_prob[,2], quiet = TRUE)

# ===== MODEL 2: SUPPORT VECTOR MACHINE =====
cat("Training Support Vector Machine...\n")
set.seed(456)
svm_model <- train(
  x = train_x, y = train_y,
  method = "svmRadial",
  trControl = cv_control,
  tuneGrid = expand.grid(C = c(0.1, 1, 10, 100), sigma = c(0.01, 0.1, 1)),
  metric = "ROC",
  preProcess = c("center", "scale")
)

svm_pred <- predict(svm_model, test_x)
svm_prob <- predict(svm_model, test_x, type = "prob")
model_results$SVM <- confusionMatrix(svm_pred, test_y)
roc_data$SVM <- roc(test_y, svm_prob[,2], quiet = TRUE)

# ===== MODEL 3: ELASTIC NET =====
cat("Training Elastic Net...\n")
set.seed(456)
glmnet_model <- train(
  x = train_x, y = train_y,
  method = "glmnet",
  trControl = cv_control,
  tuneGrid = expand.grid(alpha = seq(0, 1, 0.2), lambda = seq(0.001, 0.1, length = 10)),
  metric = "ROC",
  preProcess = c("center", "scale")
)

glmnet_pred <- predict(glmnet_model, test_x)
glmnet_prob <- predict(glmnet_model, test_x, type = "prob")
model_results$ElasticNet <- confusionMatrix(glmnet_pred, test_y)
roc_data$ElasticNet <- roc(test_y, glmnet_prob[,2], quiet = TRUE)

# ===== MODEL 4: XGBOOST =====
cat("Training XGBoost...\n")
set.seed(456)
xgb_model <- train(
  x = train_x, y = train_y,
  method = "xgbTree",
  trControl = cv_control,
  tuneGrid = expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 6),
    eta = c(0.1, 0.3),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  ),
  metric = "ROC",
  verbosity = 0
)

xgb_pred <- predict(xgb_model, test_x)
xgb_prob <- predict(xgb_model, test_x, type = "prob")
model_results$XGBoost <- confusionMatrix(xgb_pred, test_y)
roc_data$XGBoost <- roc(test_y, xgb_prob[,2], quiet = TRUE)

# ===== MODEL 5: NAIVE BAYES =====
cat("Training Naive Bayes...\n")
set.seed(456)
nb_model <- train(
  x = train_x, y = train_y,
  method = "naive_bayes",
  trControl = cv_control,
  tuneGrid = expand.grid(laplace = c(0, 1, 2), usekernel = c(TRUE, FALSE), adjust = c(0.5, 1, 1.5)),
  metric = "ROC"
)

nb_pred <- predict(nb_model, test_x)
nb_prob <- predict(nb_model, test_x, type = "prob")
model_results$NaiveBayes <- confusionMatrix(nb_pred, test_y)
roc_data$NaiveBayes <- roc(test_y, nb_prob[,2], quiet = TRUE)

# ===== BEAUTIFUL VISUALIZATIONS =====

# 1. Model Performance Comparison
performance_data <- data.frame(
  Model = names(model_results),
  Accuracy = sapply(model_results, function(x) x$overall['Accuracy']),
  Sensitivity = sapply(model_results, function(x) x$byClass['Sensitivity']),
  Specificity = sapply(model_results, function(x) x$byClass['Specificity']),
  AUC = sapply(roc_data, function(x) as.numeric(x$auc))
)

# Reshape for plotting
performance_long <- performance_data %>%
  gather(Metric, Value, -Model) %>%
  mutate(Model = factor(Model, levels = performance_data$Model[order(-performance_data$AUC)]))

# Performance comparison plot
p_performance <- ggplot(performance_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_text(aes(label = round(Value, 3)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Model Performance Comparison",
       subtitle = "Accuracy, Sensitivity, Specificity, and AUC across different models",
       x = "Machine Learning Model", y = "Performance Score") +
  ylim(0, 1.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_performance)

# 2. ROC Curves Comparison
roc_df <- data.frame()
for(model_name in names(roc_data)) {
  roc_obj <- roc_data[[model_name]]
  temp_df <- data.frame(
    Model = model_name,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    AUC = as.numeric(roc_obj$auc)
  )
  roc_df <- rbind(roc_df, temp_df)
}

# Add AUC to model names for legend
roc_df$Model_AUC <- paste0(roc_df$Model, " (AUC = ", round(roc_df$AUC, 3), ")")

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model_AUC)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves Comparison",
       subtitle = "Receiver Operating Characteristic curves for all models",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Model") +
  coord_equal() +
  theme(legend.position = "right")

print(p_roc)

# 3. Feature Importance (Random Forest)
rf_importance <- varImp(rf_model, scale = TRUE)$importance
rf_importance$piRNA <- rownames(rf_importance)
rf_importance <- rf_importance[order(-rf_importance[,1]), ]
top_20_rf <- head(rf_importance, 20)

p_rf_importance <- ggplot(top_20_rf, aes(x = reorder(piRNA, top_20_rf[,1]), y = top_20_rf[,1])) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = round(top_20_rf[,1], 1)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Top 20 Most Important piRNAs (Random Forest)",
       subtitle = "Variable importance based on mean decrease in accuracy",
       x = "piRNA", y = "Variable Importance Score") +
  theme(axis.text.y = element_text(size = 9))

print(p_rf_importance)

# 4. Model Complexity vs Performance
complexity_data <- data.frame(
  Model = names(model_results),
  AUC = sapply(roc_data, function(x) as.numeric(x$auc)),
  Parameters = c(
    rf_model$bestTune$mtry,
    length(svm_model$bestTune),
    2,  # alpha and lambda for glmnet
    nrow(xgb_model$bestTune),
    length(nb_model$bestTune)
  )
)

p_complexity <- ggplot(complexity_data, aes(x = Parameters, y = AUC)) +
  geom_point(aes(color = Model), size = 4, alpha = 0.8) +
  geom_text(aes(label = Model), vjust = -1, hjust = 0.5, size = 3) +
  scale_color_manual(values = colors) +
  labs(title = "Model Complexity vs Performance",
       subtitle = "AUC score vs number of key hyperparameters",
       x = "Number of Key Parameters", y = "AUC Score") +
  ylim(min(complexity_data$AUC) - 0.02, max(complexity_data$AUC) + 0.02) +
  theme(legend.position = "none")

print(p_complexity)

# 5. Confusion Matrix Heatmaps
create_cm_plot <- function(cm, model_name) {
  cm_data <- as.data.frame(cm$table)
  ggplot(cm_data, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste("Confusion Matrix:", model_name),
         x = "Predicted", y = "Actual") +
    theme(legend.position = "none")
}

# Create confusion matrix plots for top 3 models
top_models <- names(sort(sapply(roc_data, function(x) as.numeric(x$auc)), decreasing = TRUE))[1:3]
cm_plots <- lapply(top_models, function(x) create_cm_plot(model_results[[x]], x))
cm_grid <- do.call(grid.arrange, c(cm_plots, ncol = 3))

# 6. Cross-validation Performance Distribution
cv_results <- resamples(list(
  RandomForest = rf_model,
  SVM = svm_model,
  ElasticNet = glmnet_model,
  XGBoost = xgb_model,
  NaiveBayes = nb_model
))

# Extract CV results and convert to long format for ggplot
cv_data <- cv_results$values
cv_long <- data.frame(
  Model = rep(c("RandomForest", "SVM", "ElasticNet", "XGBoost", "NaiveBayes"), each = 10),
  ROC = c(cv_data[,2], cv_data[,4], cv_data[,6], cv_data[,8], cv_data[,10])
)

p_cv_box <- ggplot(cv_long, aes(x = Model, y = ROC)) +
  geom_boxplot(aes(fill = Model), alpha = 0.7) +
  geom_jitter(aes(y = ROC), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = colors) +
  labs(title = "Cross-Validation Performance Distribution",
       subtitle = "10-fold CV AUC scores across all models",
       x = "Model", y = "AUC Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p_cv_box)

# ===== ENHANCED SUMMARY TABLE =====
summary_table <- performance_data %>%
  arrange(desc(AUC)) %>%
  mutate(
    Accuracy = round(Accuracy, 4),
    Sensitivity = round(Sensitivity, 4),
    Specificity = round(Specificity, 4),
    AUC = round(AUC, 4),
    PR_AUC = round(performance_data$PR_AUC, 4),
    Rank = row_number()
  ) %>%
  select(Rank, Model, AUC, PR_AUC, Accuracy, Sensitivity, Specificity)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("COMPREHENSIVE MODEL COMPARISON SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n")
print(summary_table)
cat(paste(rep("-", 80), collapse=""), "\n")

# Model ranking analysis
cat("MODEL RANKING ANALYSIS:\n")
cat("- ROC-AUC Best:", summary_table$Model[which.max(summary_table$AUC)], 
    "(", max(summary_table$AUC), ")\n")
cat("- PR-AUC Best:", summary_table$Model[which.max(summary_table$PR_AUC)], 
    "(", max(summary_table$PR_AUC), ")\n")
cat("- Accuracy Best:", summary_table$Model[which.max(summary_table$Accuracy)], 
    "(", max(summary_table$Accuracy), ")\n")

# Class balance information
cat("\nCLASS BALANCE INFORMATION:\n")
cat("- Positive class ratio:", round(baseline_precision, 3), "\n")
cat("- Negative class ratio:", round(1 - baseline_precision, 3), "\n")
if(baseline_precision < 0.3 || baseline_precision > 0.7) {
  cat("- Dataset is IMBALANCED - PR curves are more informative than ROC curves\n")
} else {
  cat("- Dataset is relatively BALANCED - both ROC and PR curves are informative\n")
}

cat(paste(rep("=", 80), collapse=""), "\n")model_name, "\n")
cat("AUC Score:", best_model_auc, "\n")
cat("Training samples:", nrow(train_data), "\n")
cat("Test samples:", nrow(test_data), "\n")
cat(paste(rep("=", 70), collapse=""), "\n")

# Save models (optional)
# saveRDS(list(rf = rf_model, svm = svm_model, glmnet = glmnet_model, 
#              xgb = xgb_model, nb = nb_model), "all_pirna_models.rds")