# ============================================================
# 03_glm_models.R
# Generalised Regression Models
# Poisson GLM → Negative Binomial → Mixed Effects
# Covers: GLMs, overdispersion, model diagnostics, interpretation
# ============================================================

library(tidyverse)
library(MASS)
library(lme4)
library(sandwich)
library(lmtest)
library(ggplot2)
library(patchwork)
library(broom)

cat("=== Module 3: Generalised Regression Models ===\n")

# Load data
pollution <- readRDS("data/pollution_clean.rds")
health    <- readRDS("data/health_data.rds")

# ============================================================
# SECTION A: Prepare modelling dataset
# ============================================================

# Merge pollution and health data
model_data <- health %>%
  left_join(
    pollution %>%
      group_by(date) %>%
      summarise(
        pm2.5_mean = mean(pm2.5, na.rm = TRUE),
        no2_mean   = mean(no2, na.rm = TRUE),
        temp_mean  = mean(temp, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "date"
  ) %>%
  mutate(
    year     = as.factor(lubridate::year(date)),
    month    = as.factor(lubridate::month(date)),
    weekday  = as.factor(lubridate::wday(date)),
    # Centred and scaled PM2.5 for interpretable coefficients
    pm25_scaled = scale(pm2.5_mean)[,1],
    # Lag: pollution today affects admissions tomorrow
    pm25_lag1 = lag(pm2.5_mean, 1),
    pm25_lag2 = lag(pm2.5_mean, 2),
    pm25_lag3 = lag(pm2.5_mean, 3)
  ) %>%
  filter(complete.cases(.))

cat("Model dataset:", nrow(model_data), "complete observations\n")

# ============================================================
# SECTION B: Model 1 — Simple Poisson GLM
# ============================================================
cat("\nFitting Model 1: Simple Poisson GLM...\n")

m1_poisson <- glm(
  respiratory_admissions ~ pm2.5_mean + temp_mean + month + weekday,
  family = poisson(link = "log"),
  data   = model_data
)

cat("Model 1 summary:\n")
cat("  AIC:", round(AIC(m1_poisson), 1), "\n")

# Check for overdispersion (key diagnostic for count data)
# If residual deviance >> residual df, we have overdispersion
dispersion_ratio <- m1_poisson$deviance / m1_poisson$df.residual
cat("  Dispersion ratio:", round(dispersion_ratio, 3),
    ifelse(dispersion_ratio > 1.5, "⚠ OVERDISPERSED", "✓ OK"), "\n")

# Pearson residuals vs fitted
poisson_diag <- data.frame(
  fitted    = fitted(m1_poisson),
  residuals = residuals(m1_poisson, type = "pearson"),
  observed  = model_data$respiratory_admissions
)

p_diag1 <- ggplot(poisson_diag, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "loess", colour = "#FF7043", se = FALSE) +
  labs(
    title    = "Poisson GLM: Residuals vs Fitted",
    subtitle = paste("Dispersion ratio:", round(dispersion_ratio, 2),
                     "— values >1.5 indicate overdispersion"),
    x = "Fitted values", y = "Pearson residuals"
  ) +
  theme_minimal(base_size = 12)

# ============================================================
# SECTION C: Model 2 — Negative Binomial (handles overdispersion)
# ============================================================
cat("Fitting Model 2: Negative Binomial GLM...\n")

m2_negbin <- glm.nb(
  respiratory_admissions ~ pm2.5_mean + temp_mean + month + weekday,
  data = model_data
)

cat("Model 2 (NB) summary:\n")
cat("  AIC:", round(AIC(m2_negbin), 1), "\n")
cat("  Theta (overdispersion parameter):", round(m2_negbin$theta, 3), "\n")

# ============================================================
# SECTION D: Model 3 — With distributed lag (multi-day effect)
# ============================================================
cat("Fitting Model 3: Negative Binomial with lag structure...\n")

m3_lag <- glm.nb(
  respiratory_admissions ~ pm2.5_mean + pm25_lag1 + pm25_lag2 + pm25_lag3 +
    temp_mean + month + weekday,
  data = model_data
)

cat("Model 3 (NB + lags) AIC:", round(AIC(m3_lag), 1), "\n")

# ============================================================
# SECTION E: Compare models using AIC
# ============================================================
model_comparison <- data.frame(
  Model       = c("Poisson (simple)", "Negative Binomial", "NB + lag effects"),
  AIC         = c(AIC(m1_poisson), AIC(m2_negbin), AIC(m3_lag)),
  Deviance    = c(deviance(m1_poisson), deviance(m2_negbin), deviance(m3_lag)),
  df_residual = c(df.residual(m1_poisson), df.residual(m2_negbin), df.residual(m3_lag))
) %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

cat("\nModel Comparison Table:\n")
print(model_comparison)

# ============================================================
# SECTION F: Interpret the best model
# ============================================================
cat("\nInterpreting best model (Negative Binomial + lags)...\n")

# Extract clean coefficients
coefs <- tidy(m3_lag, conf.int = TRUE) %>%
  filter(!str_detect(term, "^month|^weekday|Intercept")) %>%
  mutate(
    # Convert log-scale coefficients to Rate Ratios (RR)
    rr      = exp(estimate),
    rr_low  = exp(conf.low),
    rr_high = exp(conf.high),
    # For PM2.5: RR per 10 µg/m³ increase
    term_label = case_when(
      term == "pm2.5_mean" ~ "PM2.5 (same day)",
      term == "pm25_lag1"  ~ "PM2.5 (lag 1 day)",
      term == "pm25_lag2"  ~ "PM2.5 (lag 2 days)",
      term == "pm25_lag3"  ~ "PM2.5 (lag 3 days)",
      term == "temp_mean"  ~ "Temperature",
      TRUE                 ~ term
    )
  )

cat("\nRate Ratios (RR) for pollution and weather:\n")
print(coefs %>% select(term_label, rr, rr_low, rr_high, p.value))

# Forest plot of Rate Ratios
p_forest <- coefs %>%
  ggplot(aes(x = rr, y = reorder(term_label, rr))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "gray50") +
  geom_errorbarh(aes(xmin = rr_low, xmax = rr_high), height = 0.25, colour = "#5B8DBE") +
  geom_point(size = 3, colour = "#FF7043") +
  labs(
    title    = "Rate Ratios: Effect of Pollution on Respiratory Admissions",
    subtitle = "Values >1 increase risk; <1 decrease risk. Error bars = 95% CI.",
    x        = "Rate Ratio (RR)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/10_forest_plot_rr.png", p_forest, width = 9, height = 5, dpi = 150)

# ============================================================
# SECTION G: Prediction plot
# ============================================================
cat("Creating prediction plot...\n")

# Predict admissions across a range of PM2.5 values
pm25_range <- seq(2, 50, by = 0.5)

pred_data <- data.frame(
  pm2.5_mean = pm25_range,
  pm25_lag1  = mean(model_data$pm25_lag1, na.rm = TRUE),
  pm25_lag2  = mean(model_data$pm25_lag2, na.rm = TRUE),
  pm25_lag3  = mean(model_data$pm25_lag3, na.rm = TRUE),
  temp_mean  = mean(model_data$temp_mean, na.rm = TRUE),
  month      = factor(7, levels = levels(model_data$month)),
  weekday    = factor(3, levels = levels(model_data$weekday))
)

pred_data$predicted     <- predict(m3_lag, newdata = pred_data, type = "response")
pred_se                 <- predict(m3_lag, newdata = pred_data, type = "link", se.fit = TRUE)
pred_data$lower         <- exp(pred_se$fit - 1.96 * pred_se$se.fit)
pred_data$upper         <- exp(pred_se$fit + 1.96 * pred_se$se.fit)

p_pred <- ggplot(pred_data, aes(x = pm2.5_mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#5B8DBE", alpha = 0.3) +
  geom_line(aes(y = predicted), colour = "#5B8DBE", linewidth = 1.2) +
  geom_vline(xintercept = 15, linetype = "dashed", colour = "orange") +
  annotate("text", x = 16.5, y = min(pred_data$lower) + 2,
           label = "WHO guideline\n(15 µg/m³)", colour = "orange", size = 3, hjust = 0) +
  labs(
    title    = "Predicted Respiratory Admissions vs PM2.5 (Negative Binomial Model)",
    subtitle = "Shaded area = 95% confidence interval. Controlling for temperature, seasonality, day-of-week.",
    x        = "PM2.5 (µg/m³)",
    y        = "Predicted daily admissions"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/11_glm_prediction.png", p_pred, width = 9, height = 5, dpi = 150)

# ============================================================
# SECTION H: Residual diagnostics for best model
# ============================================================
nb_diag <- data.frame(
  fitted    = fitted(m3_lag),
  residuals = residuals(m3_lag, type = "pearson"),
  std_res   = rstandard(m3_lag),
  date      = model_data$date
)

p_diag2 <- ggplot(nb_diag, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "loess", colour = "#FF7043", se = FALSE) +
  labs(title = "Negative Binomial: Residuals vs Fitted",
       x = "Fitted values", y = "Pearson residuals") +
  theme_minimal(base_size = 12)

p_diag3 <- ggplot(nb_diag, aes(x = date, y = residuals)) +
  geom_line(alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "loess", colour = "#FF7043", se = FALSE) +
  labs(title = "Residuals over Time — Check for Temporal Autocorrelation",
       x = NULL, y = "Pearson residuals") +
  theme_minimal(base_size = 12)

p_combined_diag <- p_diag2 / p_diag3
ggsave("outputs/12_nb_diagnostics.png", p_combined_diag, width = 10, height = 8, dpi = 150)

# ============================================================
# Save model results
# ============================================================
glm_results <- list(
  m1_poisson       = m1_poisson,
  m2_negbin        = m2_negbin,
  m3_lag           = m3_lag,
  model_comparison = model_comparison,
  coefs            = coefs,
  model_data       = model_data
)

saveRDS(glm_results, "data/glm_results.rds")

cat("\n=== Model Summary ===\n")
cat("Best model: Negative Binomial with lag effects (lowest AIC)\n")
cat("Key finding: Each 10 µg/m³ increase in PM2.5 increases admissions by ~",
    round((exp(coef(m3_lag)["pm2.5_mean"] * 10) - 1) * 100, 1), "%\n")

cat("\n✅ Module 3 complete! Plots saved to outputs/\n")
cat("Next step: source('R/04_time_series.R')\n")
