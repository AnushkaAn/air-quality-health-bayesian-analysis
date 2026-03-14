# ============================================================
# 04_time_series.R
# Time Series Analysis: ARIMA, ETS, forecasting, STL decomposition
# ============================================================

library(tidyverse)
library(forecast)
library(tseries)
library(zoo)
library(lubridate)
library(ggplot2)
library(patchwork)

cat("=== Module 4: Time Series Analysis ===\n")

# Load data
pollution <- readRDS("data/pollution_clean.rds")

# ============================================================
# SECTION A: Create a clean daily time series (one site)
# ============================================================

# Use KC1 (Kings College) as our primary series — urban background
kc1_daily <- pollution %>%
  filter(site_code == "KC1") %>%
  group_by(date) %>%
  summarise(pm2.5 = mean(pm2.5, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# Fill any missing dates with linear interpolation
full_dates <- data.frame(date = seq(min(kc1_daily$date), max(kc1_daily$date), by = "day"))
kc1_daily  <- full_dates %>%
  left_join(kc1_daily, by = "date") %>%
  mutate(pm2.5 = zoo::na.approx(pm2.5, na.rm = FALSE))

cat("Time series length:", nrow(kc1_daily), "days\n")
cat("Missing after interpolation:", sum(is.na(kc1_daily$pm2.5)), "\n")

# Convert to ts object for forecast package (daily, frequency=365)
ts_pm25 <- ts(kc1_daily$pm2.5, start = c(2018, 1), frequency = 365)

# ============================================================
# SECTION B: STL Decomposition
# ============================================================
cat("Performing STL decomposition...\n")

# STL = Seasonal and Trend decomposition using Loess
stl_fit <- stl(ts_pm25, s.window = "periodic", robust = TRUE)

# Extract components
stl_df <- data.frame(
  date      = kc1_daily$date,
  observed  = kc1_daily$pm2.5,
  trend     = as.numeric(stl_fit$time.series[, "trend"]),
  seasonal  = as.numeric(stl_fit$time.series[, "seasonal"]),
  remainder = as.numeric(stl_fit$time.series[, "remainder"])
)

# Plot decomposition manually with ggplot2
p_stl1 <- ggplot(stl_df, aes(x = date, y = observed)) +
  geom_line(alpha = 0.5, colour = "#5B8DBE", linewidth = 0.4) +
  labs(title = "Observed", y = "PM2.5", x = NULL) + theme_minimal(base_size = 10)

p_stl2 <- ggplot(stl_df, aes(x = date, y = trend)) +
  geom_line(colour = "#FF7043", linewidth = 0.8) +
  labs(title = "Trend Component", y = NULL, x = NULL) + theme_minimal(base_size = 10)

p_stl3 <- ggplot(stl_df, aes(x = date, y = seasonal)) +
  geom_line(colour = "#4CAF50", linewidth = 0.5, alpha = 0.7) +
  labs(title = "Seasonal Component", y = NULL, x = NULL) + theme_minimal(base_size = 10)

p_stl4 <- ggplot(stl_df, aes(x = date, y = remainder)) +
  geom_line(colour = "gray40", linewidth = 0.4, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  labs(title = "Remainder", y = NULL, x = NULL) + theme_minimal(base_size = 10)

p_stl_combined <- p_stl1 / p_stl2 / p_stl3 / p_stl4 +
  plot_annotation(
    title    = "STL Decomposition of Daily PM2.5 (KC1 — Kings College London)",
    subtitle = "Decomposed into trend, annual seasonality, and random remainder"
  )

ggsave("outputs/13_stl_decomposition.png", p_stl_combined, width = 12, height = 10, dpi = 150)

# ============================================================
# SECTION C: Stationarity tests
# ============================================================
cat("Testing stationarity...\n")

# Augmented Dickey-Fuller test (H0: unit root = non-stationary)
adf_result <- adf.test(na.omit(ts_pm25))
cat("ADF test p-value:", round(adf_result$p.value, 4),
    ifelse(adf_result$p.value < 0.05, "✓ Stationary", "⚠ Non-stationary"), "\n")

# KPSS test (H0: stationary — opposite null)
kpss_result <- kpss.test(na.omit(ts_pm25))
cat("KPSS test p-value:", round(kpss_result$p.value, 4),
    ifelse(kpss_result$p.value > 0.05, "✓ Stationary", "⚠ Non-stationary"), "\n")

# ============================================================
# SECTION D: ACF / PACF plots to inform ARIMA order
# ============================================================
cat("Creating ACF/PACF plots...\n")

acf_vals  <- acf(na.omit(ts_pm25),  lag.max = 60, plot = FALSE)
pacf_vals <- pacf(na.omit(ts_pm25), lag.max = 60, plot = FALSE)

acf_df  <- data.frame(lag = acf_vals$lag[,,1][-1],  acf  = acf_vals$acf[,,1][-1])
pacf_df <- data.frame(lag = pacf_vals$lag[,,1],      pacf = pacf_vals$acf[,,1])

ci <- qnorm(0.975) / sqrt(length(na.omit(ts_pm25)))

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_col(fill = "#5B8DBE", width = 0.3) +
  geom_hline(yintercept = c(ci, -ci), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  labs(title = "ACF: Autocorrelation Function",
       subtitle = "Slow decay → AR process; blue bars beyond dashed lines = significant",
       x = "Lag (days)", y = "Autocorrelation") +
  theme_minimal(base_size = 12)

p_pacf <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
  geom_col(fill = "#FF7043", width = 0.3) +
  geom_hline(yintercept = c(ci, -ci), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  labs(title = "PACF: Partial Autocorrelation Function",
       subtitle = "Cutoff after lag p → AR(p) component",
       x = "Lag (days)", y = "Partial Autocorrelation") +
  theme_minimal(base_size = 12)

p_acf_combined <- p_acf / p_pacf
ggsave("outputs/14_acf_pacf.png", p_acf_combined, width = 10, height = 7, dpi = 150)

# ============================================================
# SECTION E: Fit ARIMA models
# ============================================================
cat("Fitting ARIMA models...\n")

# Use last 90 days as test set
n_test   <- 90
n_train  <- length(na.omit(ts_pm25)) - n_test
ts_train <- head(na.omit(ts_pm25), n_train)
ts_test  <- tail(na.omit(ts_pm25), n_test)

# Auto ARIMA: automatically selects best (p,d,q) by AIC
cat("  Running auto.arima()...\n")
arima_auto <- auto.arima(
  ts_train,
  stepwise   = TRUE,
  approximation = TRUE,
  ic = "aic"
)
cat("  Best ARIMA:", arimaorder(arima_auto)[1], arimaorder(arima_auto)[2], arimaorder(arima_auto)[3], "\n")
cat("  AIC:", round(AIC(arima_auto), 1), "\n")

# ETS (Exponential Smoothing State Space)
cat("  Fitting ETS model...\n")
ets_model <- ets(ts_train, ic = "aic")
cat("  Best ETS model:", ets_model$method, "\n")

# ============================================================
# SECTION F: Forecast and evaluate
# ============================================================
cat("Generating forecasts...\n")

# Forecast 90 days
arima_fc <- forecast(arima_auto, h = n_test, level = c(80, 95))
ets_fc   <- forecast(ets_model,  h = n_test, level = c(80, 95))

# Calculate accuracy metrics
arima_acc <- accuracy(arima_fc, ts_test)
ets_acc   <- accuracy(ets_fc,   ts_test)

acc_table <- data.frame(
  Model = c("ARIMA", "ETS"),
  RMSE  = c(arima_acc["Test set", "RMSE"],  ets_acc["Test set", "RMSE"]),
  MAE   = c(arima_acc["Test set", "MAE"],   ets_acc["Test set", "MAE"]),
  MAPE  = c(arima_acc["Test set", "MAPE"],  ets_acc["Test set", "MAPE"])
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

cat("\nForecast Accuracy (test set, 90 days):\n")
print(acc_table)

# ============================================================
# SECTION G: Forecast plot
# ============================================================

# Build a tidy forecast dataframe
test_dates <- tail(kc1_daily$date, n_test)

forecast_df <- bind_rows(
  data.frame(
    date   = test_dates,
    actual = as.numeric(ts_test),
    model  = "Actual"
  ),
  data.frame(
    date  = test_dates,
    mean  = as.numeric(arima_fc$mean),
    lo80  = as.numeric(arima_fc$lower[,1]),
    hi80  = as.numeric(arima_fc$upper[,1]),
    lo95  = as.numeric(arima_fc$lower[,2]),
    hi95  = as.numeric(arima_fc$upper[,2]),
    model = "ARIMA"
  ),
  data.frame(
    date  = test_dates,
    mean  = as.numeric(ets_fc$mean),
    lo80  = as.numeric(ets_fc$lower[,1]),
    hi80  = as.numeric(ets_fc$upper[,1]),
    lo95  = as.numeric(ets_fc$lower[,2]),
    hi95  = as.numeric(ets_fc$upper[,2]),
    model = "ETS"
  )
)

# Training history for context
train_df <- data.frame(
  date  = tail(kc1_daily$date, n_train + n_test)[1:n_train],
  pm2.5 = as.numeric(ts_train)
)

p_forecast <- ggplot() +
  geom_line(data = tail(train_df, 120), aes(x = date, y = pm2.5),
            colour = "gray50", alpha = 0.7) +
  geom_ribbon(data = filter(forecast_df, model == "ARIMA"),
              aes(x = date, ymin = lo95, ymax = hi95), fill = "#5B8DBE", alpha = 0.2) +
  geom_ribbon(data = filter(forecast_df, model == "ARIMA"),
              aes(x = date, ymin = lo80, ymax = hi80), fill = "#5B8DBE", alpha = 0.3) +
  geom_line(data = filter(forecast_df, model == "ARIMA"),
            aes(x = date, y = mean), colour = "#5B8DBE", linewidth = 1) +
  geom_ribbon(data = filter(forecast_df, model == "ETS"),
              aes(x = date, ymin = lo95, ymax = hi95), fill = "#FF7043", alpha = 0.15) +
  geom_line(data = filter(forecast_df, model == "ETS"),
            aes(x = date, y = mean), colour = "#FF7043", linewidth = 1, linetype = "dashed") +
  geom_line(data = filter(forecast_df, model == "Actual"),
            aes(x = date, y = actual), colour = "black", linewidth = 0.8) +
  annotate("text", x = test_dates[5], y = max(kc1_daily$pm2.5, na.rm=TRUE)*0.85,
           label = "Blue = ARIMA\nOrange = ETS\nBlack = Actual",
           size = 3, hjust = 0) +
  labs(
    title    = "90-Day PM2.5 Forecast: ARIMA vs ETS",
    subtitle = "Shaded = 80% and 95% prediction intervals. Black line = actual observed values.",
    x = NULL, y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/15_forecast.png", p_forecast, width = 12, height = 6, dpi = 150)

# ============================================================
# SECTION H: Rolling 7-day forecast (operational use case)
# ============================================================
cat("Computing rolling forecasts...\n")

# Show how we'd use this operationally: weekly forecasts
future_fc <- forecast(arima_auto, h = 30)

future_df <- data.frame(
  date = max(kc1_daily$date) + 1:30,
  mean = as.numeric(future_fc$mean),
  lo95 = as.numeric(future_fc$lower[,2]),
  hi95 = as.numeric(future_fc$upper[,2])
)

p_future <- ggplot() +
  geom_line(data = tail(kc1_daily, 90), aes(x = date, y = pm2.5),
            colour = "gray50") +
  geom_ribbon(data = future_df, aes(x = date, ymin = lo95, ymax = hi95),
              fill = "#5B8DBE", alpha = 0.3) +
  geom_line(data = future_df, aes(x = date, y = mean),
            colour = "#5B8DBE", linewidth = 1.2) +
  geom_vline(xintercept = max(kc1_daily$date),
             linetype = "dashed", colour = "red") +
  annotate("text", x = max(kc1_daily$date) + 2, y = 30,
           label = "Forecast →", size = 3.5, colour = "red") +
  labs(
    title    = "30-Day Operational PM2.5 Forecast (ARIMA)",
    subtitle = "Gray = recent history; Blue = forecast with 95% CI",
    x = NULL, y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/16_future_forecast.png", p_future, width = 10, height = 5, dpi = 150)

# Save results
ts_results <- list(
  arima_model  = arima_auto,
  ets_model    = ets_model,
  acc_table    = acc_table,
  future_fc    = future_df,
  stl_fit      = stl_df,
  kc1_daily    = kc1_daily
)

saveRDS(ts_results, "data/ts_results.rds")

cat("\n=== Time Series Summary ===\n")
cat("Best model:", ifelse(acc_table$RMSE[1] < acc_table$RMSE[2], "ARIMA", "ETS"), "(lower RMSE)\n")
cat("ARIMA RMSE:", round(acc_table$RMSE[1], 3), "\n")
cat("ETS RMSE:  ", round(acc_table$RMSE[2], 3), "\n")

cat("\n✅ Module 4 complete! Plots saved to outputs/\n")
cat("Next step: source('R/05_bayesian_model.R')\n")
