# ============================================================
# 07_ml_comparison.R
# Machine Learning: Random Forest & XGBoost
# Compare predictive performance against ARIMA
# Variable importance analysis
# ============================================================

library(tidyverse)
library(lubridate)
library(ranger)
library(xgboost)
library(caret)
library(vip)
library(ggplot2)
library(patchwork)

cat("=== Module 7: Machine Learning Comparison ===\n")

# ============================================================
# Load and prepare features
# ============================================================
pollution <- readRDS("data/pollution_clean.rds")
ts_results <- readRDS("data/ts_results.rds")

# Build a rich feature matrix from all sites
kc1_daily <- ts_results$kc1_daily

# Create lagged features (ML models can use many lags easily; ARIMA is constrained)
ml_data <- kc1_daily %>%
  rename(pm2.5 = pm2.5) %>%
  mutate(
    # Target: next-day PM2.5
    pm2.5_next = lead(pm2.5, 1),
    
    # Lag features
    pm2.5_lag1  = lag(pm2.5, 1),
    pm2.5_lag2  = lag(pm2.5, 2),
    pm2.5_lag3  = lag(pm2.5, 3),
    pm2.5_lag7  = lag(pm2.5, 7),
    pm2.5_lag14 = lag(pm2.5, 14),
    pm2.5_lag30 = lag(pm2.5, 30),
    
    # Rolling statistics
    pm2.5_roll7  = zoo::rollmean(pm2.5, 7,  fill = NA, align = "right"),
    pm2.5_roll14 = zoo::rollmean(pm2.5, 14, fill = NA, align = "right"),
    pm2.5_sd7    = zoo::rollapply(pm2.5, 7, sd, fill = NA, align = "right"),
    
    # Calendar features
    month      = lubridate::month(date),
    day_of_year = lubridate::yday(date),
    year       = lubridate::year(date),
    weekday    = lubridate::wday(date),
    is_weekend = as.integer(weekday %in% c(1, 7)),
    
    # Cyclical encoding (sin/cos) for calendar features
    month_sin  = sin(2 * pi * month / 12),
    month_cos  = cos(2 * pi * month / 12),
    doy_sin    = sin(2 * pi * day_of_year / 365),
    doy_cos    = cos(2 * pi * day_of_year / 365),
    
    # Trend
    t          = as.numeric(date - min(date))
  ) %>%
  filter(complete.cases(.))

cat("Feature matrix:", nrow(ml_data), "rows ×", ncol(ml_data), "columns\n")

# ============================================================
# SECTION A: Train/test split
# ============================================================

# Use last 90 days as test set (same as time series module for fair comparison)
n_test  <- 90
n_train <- nrow(ml_data) - n_test

train_data <- head(ml_data, n_train)
test_data  <- tail(ml_data, n_test)

feature_cols <- c(
  "pm2.5_lag1", "pm2.5_lag2", "pm2.5_lag3", "pm2.5_lag7", "pm2.5_lag14", "pm2.5_lag30",
  "pm2.5_roll7", "pm2.5_roll14", "pm2.5_sd7",
  "month_sin", "month_cos", "doy_sin", "doy_cos",
  "weekday", "is_weekend", "t"
)

X_train <- as.matrix(train_data[, feature_cols])
y_train <- train_data$pm2.5_next

X_test  <- as.matrix(test_data[, feature_cols])
y_test  <- test_data$pm2.5_next

cat("Training observations:", n_train, "\n")
cat("Test observations:    ", n_test, "\n")
cat("Features used:        ", length(feature_cols), "\n")

# ============================================================
# SECTION B: Random Forest
# ============================================================
cat("\nTraining Random Forest...\n")

set.seed(42)
rf_model <- ranger(
  pm2.5_next ~ .,
  data            = train_data[, c(feature_cols, "pm2.5_next")],
  num.trees       = 500,
  mtry            = floor(sqrt(length(feature_cols))),  # sqrt(p) features per split
  min.node.size   = 5,
  importance      = "impurity",     # Enable variable importance
  seed            = 42
)

rf_preds  <- predict(rf_model, data = test_data[, feature_cols])$predictions
rf_rmse   <- sqrt(mean((rf_preds - y_test)^2))
rf_mae    <- mean(abs(rf_preds - y_test))
rf_mape   <- mean(abs((rf_preds - y_test) / y_test)) * 100

cat("Random Forest results:\n")
cat("  RMSE:", round(rf_rmse, 3), "\n")
cat("  MAE: ", round(rf_mae,  3), "\n")
cat("  MAPE:", round(rf_mape, 2), "%\n")

# ============================================================
# SECTION C: XGBoost with cross-validation
# ============================================================
cat("\nTraining XGBoost with 5-fold cross-validation...\n")

# XGBoost needs a DMatrix
dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test,  label = y_test)

# Hyperparameters (tuned manually; could use grid search)
xgb_params <- list(
  objective        = "reg:squarederror",
  eval_metric      = "rmse",
  eta              = 0.05,       # learning rate
  max_depth        = 6,          # tree depth
  subsample        = 0.8,        # row subsampling (prevents overfitting)
  colsample_bytree = 0.8,        # column subsampling
  min_child_weight = 3,
  gamma            = 0.1
)

# Cross-validation to find optimal number of trees
set.seed(42)
xgb_cv <- xgb.cv(
  params   = xgb_params,
  data     = dtrain,
  nrounds  = 500,
  nfold    = 5,
  verbose  = 0,
  early_stopping_rounds = 30
)

best_nrounds <- xgb_cv$best_iteration
cat("  Best number of trees (CV):", best_nrounds, "\n")

# Final model
xgb_model <- xgb.train(
  params  = xgb_params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

xgb_preds <- predict(xgb_model, newdata = dtest)
xgb_rmse  <- sqrt(mean((xgb_preds - y_test)^2))
xgb_mae   <- mean(abs(xgb_preds - y_test))
xgb_mape  <- mean(abs((xgb_preds - y_test) / y_test)) * 100

cat("XGBoost results:\n")
cat("  RMSE:", round(xgb_rmse, 3), "\n")
cat("  MAE: ", round(xgb_mae,  3), "\n")
cat("  MAPE:", round(xgb_mape, 2), "%\n")

# ============================================================
# SECTION D: Comparison table (all models)
# ============================================================
ts_results_acc <- ts_results$acc_table

comparison_table <- data.frame(
  Model      = c("ARIMA", "ETS", "Random Forest", "XGBoost"),
  RMSE       = c(ts_results_acc$RMSE[1], ts_results_acc$RMSE[2], rf_rmse, xgb_rmse),
  MAE        = c(ts_results_acc$MAE[1],  ts_results_acc$MAE[2],  rf_mae,  xgb_mae),
  MAPE       = c(ts_results_acc$MAPE[1], ts_results_acc$MAPE[2], rf_mape, xgb_mape),
  Type       = c("Statistical", "Statistical", "ML", "ML")
) %>%
  arrange(RMSE) %>%
  mutate(
    RMSE = round(RMSE, 3),
    MAE  = round(MAE,  3),
    MAPE = round(MAPE, 2),
    Rank = 1:4
  )

cat("\n=== FULL MODEL COMPARISON ===\n")
print(comparison_table)

# ============================================================
# SECTION E: Visualise predictions vs actual
# ============================================================
cat("Creating prediction comparison plot...\n")

pred_comparison <- data.frame(
  date   = test_data$date,
  actual = y_test,
  rf     = rf_preds,
  xgb    = xgb_preds
)

p_pred_compare <- ggplot(pred_comparison, aes(x = date)) +
  geom_line(aes(y = actual), colour = "black",   linewidth = 0.8, alpha = 0.9) +
  geom_line(aes(y = rf),     colour = "#5B8DBE", linewidth = 0.8, linetype = "solid") +
  geom_line(aes(y = xgb),    colour = "#FF7043", linewidth = 0.8, linetype = "dashed") +
  labs(
    title    = "ML Model Predictions vs Actual PM2.5 (Test Set: 90 Days)",
    subtitle = "Black = actual; Blue = Random Forest; Orange = XGBoost",
    x = NULL, y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/25_ml_predictions.png", p_pred_compare, width = 12, height = 5, dpi = 150)

# Scatter: predicted vs actual
p_scatter_rf <- ggplot(pred_comparison, aes(x = actual, y = rf)) +
  geom_abline(colour = "red", linetype = "dashed") +
  geom_point(alpha = 0.4, colour = "#5B8DBE", size = 1.5) +
  labs(title = paste0("Random Forest\nRMSE=", round(rf_rmse, 2)),
       x = "Actual PM2.5", y = "Predicted PM2.5") +
  theme_minimal(base_size = 11)

p_scatter_xgb <- ggplot(pred_comparison, aes(x = actual, y = xgb)) +
  geom_abline(colour = "red", linetype = "dashed") +
  geom_point(alpha = 0.4, colour = "#FF7043", size = 1.5) +
  labs(title = paste0("XGBoost\nRMSE=", round(xgb_rmse, 2)),
       x = "Actual PM2.5", y = "Predicted PM2.5") +
  theme_minimal(base_size = 11)

p_scatter_combined <- p_scatter_rf | p_scatter_xgb
ggsave("outputs/26_ml_scatter.png", p_scatter_combined, width = 10, height = 5, dpi = 150)

# ============================================================
# SECTION F: Variable Importance
# ============================================================
cat("Computing variable importance...\n")

# Random Forest variable importance
rf_importance <- data.frame(
  variable   = names(rf_model$variable.importance),
  importance = rf_model$variable.importance
) %>%
  arrange(desc(importance)) %>%
  head(15)

p_rf_imp <- ggplot(rf_importance, aes(x = importance, y = reorder(variable, importance))) +
  geom_col(fill = "#5B8DBE", alpha = 0.85) +
  labs(
    title    = "Random Forest: Variable Importance (Top 15)",
    subtitle = "Most important predictors for next-day PM2.5",
    x = "Importance (Impurity)", y = NULL
  ) +
  theme_minimal(base_size = 12)

# XGBoost variable importance
xgb_importance_raw <- xgb.importance(
  feature_names = feature_cols,
  model         = xgb_model
) %>%
  as.data.frame() %>%
  head(15)

p_xgb_imp <- ggplot(xgb_importance_raw,
                    aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_col(fill = "#FF7043", alpha = 0.85) +
  labs(
    title    = "XGBoost: Variable Importance (Gain, Top 15)",
    subtitle = "Features that most reduce prediction error",
    x = "Gain", y = NULL
  ) +
  theme_minimal(base_size = 12)

p_importance_combined <- p_rf_imp | p_xgb_imp
ggsave("outputs/27_variable_importance.png", p_importance_combined, width = 13, height = 6, dpi = 150)

# ============================================================
# SECTION G: RMSE comparison bar chart
# ============================================================
p_comparison <- comparison_table %>%
  mutate(Model = factor(Model, levels = comparison_table$Model)) %>%
  ggplot(aes(x = Model, y = RMSE, fill = Type)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(RMSE, 2)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Statistical" = "#5B8DBE", "ML" = "#FF7043")) +
  labs(
    title    = "Model Comparison: RMSE on 90-Day Test Set",
    subtitle = "Lower = better. Statistical models may outperform ML with limited data.",
    x        = NULL, y = "RMSE (µg/m³)", fill = "Model type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("outputs/28_model_comparison.png", p_comparison, width = 8, height = 5, dpi = 150)

# ============================================================
# SECTION H: Save all results for Shiny dashboard
# ============================================================
ml_results <- list(
  rf_model          = rf_model,
  xgb_model         = xgb_model,
  comparison_table  = comparison_table,
  pred_comparison   = pred_comparison,
  rf_importance     = rf_importance,
  xgb_importance    = xgb_importance_raw,
  feature_cols      = feature_cols,
  ml_data           = ml_data
)

saveRDS(ml_results, "data/ml_results.rds")

cat("\n=== Final Model Comparison ===\n")
print(comparison_table)

cat("\n🏆 Best model (lowest RMSE):", comparison_table$Model[1], "\n")
cat("\nKey insight: ML models capture non-linear relationships and interactions\n")
cat("that ARIMA cannot, but interpretability is harder to achieve.\n")

cat("\n✅ Module 7 complete! All analysis modules done.\n")
cat("Next step: rmarkdown::render('reports/full_report.Rmd')\n")
cat("OR: shiny::runApp('shiny/app.R')\n")
