# ============================================================
# 06_causal_inference.R
# Causal Inference: Interrupted Time Series (ITS) Analysis
# Evaluating the causal effect of the London ULEZ expansion
# (25 October 2021) on PM2.5 and respiratory admissions
# ============================================================

library(tidyverse)
library(lubridate)
library(sandwich)
library(lmtest)
library(ggplot2)
library(patchwork)
library(broom)

cat("=== Module 6: Causal Inference ===\n")
cat("Evaluating the causal impact of the London ULEZ expansion (Oct 2021)\n\n")

# ============================================================
# THE CAUSAL QUESTION
# ============================================================
# Q: Did the London ULEZ expansion CAUSE a reduction in air pollution?
# 
# Problem: Pollution was already declining before ULEZ (confounding trend).
# Naive before/after comparison would overestimate the ULEZ effect.
#
# Solution: Interrupted Time Series (ITS) design
# - Model the pre-intervention TREND
# - Estimate the CHANGE in level AND slope at the intervention point
# - This separates "ULEZ effect" from "pre-existing trend"
#
# Bonus: Difference-in-Differences (DiD)
# - Compare inner-city sites (treated by ULEZ) vs outer London (control)
# - This removes any London-wide confounders (e.g. COVID, economic factors)

# ============================================================
# Load and prepare data
# ============================================================
pollution <- readRDS("data/pollution_clean.rds")
health    <- readRDS("data/health_data.rds")

ulez_date <- as.Date("2021-10-25")

# Create ITS dataset for KC1 (urban — treated by ULEZ)
its_data <- pollution %>%
  filter(site_code == "KC1", !is.na(pm2.5)) %>%
  group_by(date) %>%
  summarise(pm2.5 = mean(pm2.5, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    # Time index from start of series
    t               = as.numeric(date - min(date)),
    # Intervention indicator: 0 before ULEZ, 1 after
    post_ulez       = as.integer(date >= ulez_date),
    # Time since intervention (0 before ULEZ)
    t_post          = pmax(0, as.numeric(date - ulez_date)),
    # Exclude 2020 COVID period to avoid contamination
    covid_period    = as.integer(date >= as.Date("2020-03-01") & date <= as.Date("2021-06-30")),
    # Season dummies
    month           = lubridate::month(date),
    sin_annual      = sin(2 * pi * month / 12),
    cos_annual      = cos(2 * pi * month / 12)
  )

# For DiD: compare KC1 (inner London, ULEZ zone) vs HRL (Heathrow, outside ULEZ)
did_data <- pollution %>%
  filter(site_code %in% c("KC1", "HRL"), !is.na(pm2.5)) %>%
  group_by(date, site_code) %>%
  summarise(pm2.5 = mean(pm2.5, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    treated    = as.integer(site_code == "KC1"),   # KC1 = in ULEZ zone
    post_ulez  = as.integer(date >= ulez_date),
    # DiD interaction term: the causal effect
    treated_x_post = treated * post_ulez,
    month      = lubridate::month(date),
    sin_annual = sin(2 * pi * month / 12),
    cos_annual = cos(2 * pi * month / 12),
    covid      = as.integer(date >= as.Date("2020-03-01") & date <= as.Date("2021-06-30"))
  )

# ============================================================
# SECTION A: Interrupted Time Series (ITS) Model
# ============================================================
cat("Fitting Interrupted Time Series model...\n")

# ITS model: pm2.5 ~ t + post_ulez + t_post + seasonality
# t           = pre-existing trend (slope before intervention)
# post_ulez   = immediate level change at ULEZ date
# t_post      = change in slope after ULEZ

its_model <- lm(
  pm2.5 ~ t + post_ulez + t_post + sin_annual + cos_annual,
  data   = filter(its_data, covid_period == 0)
)

# Use Newey-West standard errors (robust to autocorrelation in time series)
its_se     <- NeweyWest(its_model, lag = 14, prewhite = FALSE)
its_robust <- coeftest(its_model, vcov = its_se)

cat("\nITS Model Results (Newey-West robust SE):\n")
print(its_robust)

# Extract key estimates
its_coefs <- tidy(its_robust) %>%
  mutate(
    interpretation = case_when(
      term == "t"           ~ "Pre-ULEZ trend (µg/m³ per day)",
      term == "post_ulez"   ~ "Immediate level change at ULEZ",
      term == "t_post"      ~ "Change in slope after ULEZ",
      term == "sin_annual"  ~ "Seasonal component (sin)",
      term == "cos_annual"  ~ "Seasonal component (cos)",
      term == "(Intercept)" ~ "Baseline PM2.5",
      TRUE ~ term
    )
  )

level_change  <- its_coefs$estimate[its_coefs$term == "post_ulez"]
slope_change  <- its_coefs$estimate[its_coefs$term == "t_post"]

cat("\nKey ITS findings:\n")
cat("  Immediate level change:", round(level_change, 3), "µg/m³\n")
cat("  Change in daily slope:", round(slope_change, 5), "µg/m³ per day\n")
cat("  Annualised slope change:", round(slope_change * 365, 3), "µg/m³ per year\n")

# ============================================================
# SECTION B: ITS Visualisation — the key plot
# ============================================================
cat("Creating ITS visualisation...\n")

# Generate counterfactual (what would have happened without ULEZ?)
its_plot_data <- its_data %>%
  filter(covid_period == 0) %>%
  mutate(
    # Fitted (with ULEZ effect)
    fitted = predict(its_model, newdata = .),
    # Counterfactual: set post_ulez = 0 and t_post = 0
    counterfactual = predict(its_model, newdata = mutate(., post_ulez = 0, t_post = 0))
  )

p_its <- ggplot(its_plot_data) +
  # Raw data
  geom_point(aes(x = date, y = pm2.5), alpha = 0.15, size = 0.6, colour = "gray40") +
  # Fitted model
  geom_line(aes(x = date, y = fitted), colour = "#5B8DBE", linewidth = 1) +
  # Counterfactual (no ULEZ)
  geom_line(aes(x = date, y = counterfactual),
            colour = "#FF7043", linewidth = 1, linetype = "dashed") +
  # ULEZ line
  geom_vline(xintercept = ulez_date, colour = "red", linewidth = 1) +
  annotate("text", x = ulez_date + 15, y = 42,
           label = "ULEZ expansion\n25 Oct 2021", colour = "red", size = 3.5, hjust = 0) +
  # Label counterfactual
  annotate("text", x = as.Date("2023-01-01"), y = 28,
           label = "Counterfactual\n(no ULEZ)", colour = "#FF7043", size = 3.5, hjust = 0) +
  annotate("text", x = as.Date("2023-01-01"), y = 20,
           label = "Observed\n(with ULEZ)", colour = "#5B8DBE", size = 3.5, hjust = 0) +
  labs(
    title    = "Interrupted Time Series: Effect of ULEZ on PM2.5 (KC1 — Kings College London)",
    subtitle = "Blue = fitted model; Orange dashed = counterfactual (predicted without ULEZ). Grey = raw daily data.",
    x        = NULL,
    y        = "PM2.5 (µg/m³)"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/22_its_main.png", p_its, width = 12, height = 6, dpi = 150)

# ============================================================
# SECTION C: Difference-in-Differences (DiD)
# ============================================================
cat("\nFitting Difference-in-Differences model...\n")

# DiD logic:
# Effect = (Treated_after - Treated_before) - (Control_after - Control_before)
# The control group (HRL) absorbs London-wide confounders (weather, COVID recovery, etc.)

did_model <- lm(
  pm2.5 ~ treated * post_ulez + sin_annual + cos_annual + covid,
  data = did_data
)

did_se     <- NeweyWest(did_model, lag = 14, prewhite = FALSE)
did_robust <- coeftest(did_model, vcov = did_se)

cat("DiD Model Results:\n")
print(did_robust)

did_effect <- coef(did_model)["treated:post_ulez"]
did_se_val <- sqrt(did_se["treated:post_ulez", "treated:post_ulez"])

cat("\nDiD Estimate: ULEZ reduced PM2.5 by",
    round(-did_effect, 2), "µg/m³\n")
cat("  95% CI:", round(-did_effect - 1.96*did_se_val, 2),
    "to", round(-did_effect + 1.96*did_se_val, 2), "µg/m³\n")

# ============================================================
# SECTION D: DiD visualisation — parallel trends check
# ============================================================
cat("Creating DiD parallel trends plot...\n")

# Monthly averages for cleaner visualisation
did_monthly <- did_data %>%
  filter(covid == 0) %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month, site_code, treated) %>%
  summarise(pm2.5 = mean(pm2.5, na.rm = TRUE), .groups = "drop")

p_did <- ggplot(did_monthly, aes(x = year_month, y = pm2.5, colour = site_code)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = ulez_date, colour = "red", linetype = "dashed") +
  annotate("text", x = ulez_date + 20, y = 26,
           label = "ULEZ", colour = "red", size = 3) +
  scale_colour_manual(
    values = c("KC1" = "#5B8DBE", "HRL" = "#FF7043"),
    labels = c("KC1" = "Kings College (inner, treated)", "HRL" = "Heathrow (outer, control)")
  ) +
  labs(
    title    = "Difference-in-Differences: Treated vs Control Site Monthly PM2.5",
    subtitle = "Before ULEZ: both sites trend similarly (parallel trends assumption). After: diverge.",
    x        = NULL, y = "Mean monthly PM2.5 (µg/m³)", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("outputs/23_did_parallel_trends.png", p_did, width = 11, height = 6, dpi = 150)

# ============================================================
# SECTION E: Health outcome ITS — did ULEZ reduce admissions?
# ============================================================
cat("\nFitting ITS model for health outcomes...\n")

health_its <- health %>%
  mutate(
    t             = as.numeric(date - min(date)),
    post_ulez     = as.integer(date >= ulez_date),
    t_post        = pmax(0, as.numeric(date - ulez_date)),
    month         = lubridate::month(date),
    sin_annual    = sin(2 * pi * month / 12),
    cos_annual    = cos(2 * pi * month / 12),
    covid         = as.integer(date >= as.Date("2020-03-01") & date <= as.Date("2021-06-30"))
  ) %>%
  filter(covid == 0)

health_its_model <- glm(
  respiratory_admissions ~ t + post_ulez + t_post + sin_annual + cos_annual,
  family = poisson(link = "log"),
  data   = health_its
)

health_its_se     <- NeweyWest(health_its_model, lag = 14, prewhite = FALSE)
health_its_robust <- coeftest(health_its_model, vcov = health_its_se)

health_ulez_rr <- exp(coef(health_its_model)["post_ulez"])
cat("ULEZ effect on admissions: RR =", round(health_ulez_rr, 3),
    "(", round((health_ulez_rr - 1) * 100, 1), "% change)\n")

# Plot health ITS
health_plot_data <- health_its %>%
  mutate(
    fitted        = predict(health_its_model, type = "response"),
    counterfactual = predict(health_its_model, 
                             newdata = mutate(., post_ulez = 0, t_post = 0),
                             type = "response")
  )

p_health_its <- ggplot(health_plot_data) +
  geom_point(aes(x = date, y = respiratory_admissions), alpha = 0.15, size = 0.6, colour = "gray40") +
  geom_line(aes(x = date, y = fitted),        colour = "#5B8DBE", linewidth = 1) +
  geom_line(aes(x = date, y = counterfactual), colour = "#FF7043", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = ulez_date, colour = "red") +
  annotate("text", x = ulez_date + 15, y = 75, label = "ULEZ", colour = "red", size = 3.5) +
  labs(
    title    = "ITS: Effect of ULEZ on Respiratory Hospital Admissions",
    subtitle = "Blue = observed; Orange dashed = counterfactual without ULEZ",
    x = NULL, y = "Daily respiratory admissions"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/24_health_its.png", p_health_its, width = 12, height = 5, dpi = 150)

# ============================================================
# SECTION F: Results summary table
# ============================================================
causal_summary <- data.frame(
  Analysis              = c("ITS (pollution)", "DiD (pollution)", "ITS (health)"),
  Outcome               = c("PM2.5 change (µg/m³)", "PM2.5 change (µg/m³)", "Admission change (RR)"),
  Estimate              = c(round(level_change, 3),
                             round(did_effect, 3),
                             round(health_ulez_rr, 3)),
  Interpretation        = c("Immediate reduction in PM2.5",
                             "Causal reduction (DiD removes confounders)",
                             "Rate ratio of admissions post-ULEZ")
)

cat("\nCausal Inference Summary:\n")
print(causal_summary)

# Save results
causal_results <- list(
  its_model      = its_model,
  did_model      = did_model,
  health_model   = health_its_model,
  causal_summary = causal_summary,
  its_plot_data  = its_plot_data
)
saveRDS(causal_results, "data/causal_results.rds")

cat("\n✅ Module 6 complete! Plots saved to outputs/\n")
cat("Next step: source('R/07_ml_comparison.R')\n")
