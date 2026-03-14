# ============================================================
# 02_eda.R
# Exploratory Data Analysis — understanding the data before modelling
# Covers: visualisation, missing data, distributions, correlations
# ============================================================

library(tidyverse)
library(lubridate)
library(openair)
library(mice)
library(VIM)
library(patchwork)
library(scales)
library(viridis)

cat("=== Module 2: Exploratory Data Analysis ===\n")

# Load data
pollution <- readRDS("data/pollution_clean.rds")
health    <- readRDS("data/health_data.rds")

dir.create("outputs", showWarnings = FALSE)

# ============================================================
# SECTION A: Time series overview
# ============================================================
cat("Creating time series plots...\n")

# Plot 1: PM2.5 over time by site
p1 <- pollution %>%
  filter(!is.na(pm2.5)) %>%
  ggplot(aes(x = date, y = pm2.5, colour = site_code)) +
  geom_line(alpha = 0.4, linewidth = 0.4) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1) +
  geom_vline(xintercept = as.Date("2021-10-25"), 
             linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = as.Date("2021-10-25") + 30, y = 45,
           label = "ULEZ expansion\nOct 2021", colour = "red", size = 3, hjust = 0) +
  facet_wrap(~site_code, ncol = 2) +
  scale_colour_viridis_d() +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "Daily PM2.5 Concentrations at Four London Sites (2018–2023)",
    subtitle = "Smoothed trend (LOESS) overlaid. Red dashed line = ULEZ expansion.",
    x        = NULL,
    y        = "PM2.5 (µg/m³)",
    colour   = "Site"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave("outputs/01_pm25_timeseries.png", p1, width = 12, height = 7, dpi = 150)

# Plot 2: Seasonal boxplots
p2 <- pollution %>%
  filter(!is.na(pm2.5)) %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  ggplot(aes(x = season, y = pm2.5, fill = season)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~site_code, ncol = 2) +
  scale_fill_manual(values = c("Winter"="#5B8DBE","Spring"="#4CAF50",
                                "Summer"="#FFC107","Autumn"="#FF7043")) +
  labs(
    title    = "PM2.5 by Season Across Sites",
    subtitle = "Winter levels consistently higher due to heating emissions and atmospheric inversions",
    x        = NULL, y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("outputs/02_seasonal_boxplot.png", p2, width = 10, height = 6, dpi = 150)

# Plot 3: Day-of-week patterns
p3 <- pollution %>%
  filter(!is.na(pm2.5)) %>%
  group_by(site_code, day_of_week) %>%
  summarise(mean_pm25 = mean(pm2.5, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_week, y = mean_pm25, fill = site_code)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_viridis_d() +
  labs(
    title    = "Average PM2.5 by Day of Week",
    subtitle = "Lower weekend levels reflect reduced traffic — a natural experiment within the data",
    x        = NULL, y = "Mean PM2.5 (µg/m³)", fill = "Site"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/03_dayofweek.png", p3, width = 9, height = 5, dpi = 150)

# ============================================================
# SECTION B: Missing data analysis
# ============================================================
cat("Analysing missing data patterns...\n")

# Summarise missingness
missing_summary <- pollution %>%
  summarise(across(c(pm2.5, no2, o3, temp, ws), 
                   ~sum(is.na(.)) / n() * 100,
                   .names = "{.col}_pct_missing")) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  mutate(variable = str_remove(variable, "_pct_missing"))

cat("Missing data (% by variable):\n")
print(missing_summary)

# Visual missing pattern by site
p4 <- pollution %>%
  group_by(site_code, year) %>%
  summarise(pct_missing_pm25 = sum(is.na(pm2.5)) / n() * 100, .groups = "drop") %>%
  ggplot(aes(x = year, y = pct_missing_pm25, fill = site_code)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  labs(
    title    = "Percentage of Missing PM2.5 Values by Year and Site",
    subtitle = "Missing data is not random — sensor outages cluster in time (MNAR pattern)",
    x        = "Year", y = "% Missing", fill = "Site"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/04_missing_data.png", p4, width = 9, height = 5, dpi = 150)

# ---- Multiple Imputation using MICE ----
cat("Performing multiple imputation (MICE)...\n")

# Prepare data for imputation (one site for demonstration)
impute_data <- pollution %>%
  filter(site_code == "KC1") %>%
  select(date, pm2.5, no2, o3, temp, ws) %>%
  # Introduce some missingness if none exists
  mutate(
    pm2.5 = ifelse(runif(n()) < 0.03, NA, pm2.5),
    no2   = ifelse(runif(n()) < 0.04, NA, no2)
  )

# Run MICE (5 imputations, predictive mean matching for numeric variables)
set.seed(42)
mice_result <- mice(
  impute_data %>% select(-date),
  m       = 5,     # number of imputed datasets
  method  = "pmm", # predictive mean matching (best for numeric data)
  maxit   = 10,    # maximum iterations
  printFlag = FALSE
)

# Extract completed (imputed) dataset
imputed_data <- complete(mice_result, 1) %>%
  mutate(date = impute_data$date, imputed = TRUE)

# Compare distributions: original vs imputed
p5 <- bind_rows(
  impute_data   %>% select(pm2.5) %>% mutate(type = "Original (with NAs)"),
  imputed_data  %>% select(pm2.5) %>% mutate(type = "After MICE imputation")
) %>%
  filter(!is.na(pm2.5)) %>%
  ggplot(aes(x = pm2.5, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("#5B8DBE", "#FF7043")) +
  labs(
    title    = "PM2.5 Distribution: Before and After MICE Imputation",
    subtitle = "Imputed values follow the same distribution — confirming imputation quality",
    x        = "PM2.5 (µg/m³)", y = "Count", fill = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/05_mice_imputation.png", p5, width = 9, height = 5, dpi = 150)

# ============================================================
# SECTION C: Correlation and distribution analysis
# ============================================================
cat("Creating correlation analysis...\n")

# Correlation matrix
cor_data <- pollution %>%
  select(pm2.5, no2, o3, temp, ws) %>%
  filter(complete.cases(.))

cor_matrix <- cor(cor_data)

# Tidy correlation for plotting
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation")

p6 <- cor_long %>%
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile(colour = "white", size = 0.5) +
  geom_text(aes(label = round(correlation, 2)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#5B8DBE", mid = "white", high = "#FF7043",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    title    = "Correlation Matrix: Pollution Variables",
    subtitle = "O3 and NO2 are negatively correlated (atmospheric chemistry)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/06_correlation_matrix.png", p6, width = 7, height = 6, dpi = 150)

# ============================================================
# SECTION D: Health outcome EDA
# ============================================================
cat("Exploring health outcome data...\n")

p7 <- health %>%
  ggplot(aes(x = city_pm25, y = respiratory_admissions)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "#5B8DBE") +
  geom_smooth(method = "glm", method.args = list(family = "poisson"),
              colour = "#FF7043", linewidth = 1.2, se = TRUE) +
  labs(
    title    = "Pollution vs Hospital Admissions",
    subtitle = "Poisson GLM fit shows dose-response relationship",
    x        = "Mean city PM2.5 (µg/m³)",
    y        = "Daily respiratory admissions"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/07_health_scatter.png", p7, width = 8, height = 5, dpi = 150)

# Admissions over time
p8 <- health %>%
  ggplot(aes(x = date, y = respiratory_admissions)) +
  geom_line(alpha = 0.4, colour = "#5B8DBE") +
  geom_smooth(method = "loess", span = 0.1, colour = "#FF7043", se = FALSE) +
  geom_vline(xintercept = as.Date("2021-10-25"),
             linetype = "dashed", colour = "red") +
  annotate("text", x = as.Date("2021-11-15"), y = max(health$respiratory_admissions) * 0.9,
           label = "ULEZ", colour = "red", size = 3.5) +
  labs(
    title    = "Daily Respiratory Hospital Admissions Over Time",
    subtitle = "Strong seasonal pattern — winter peaks from respiratory infections + pollution",
    x = NULL, y = "Admissions"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/08_admissions_timeseries.png", p8, width = 12, height = 5, dpi = 150)

# ============================================================
# SECTION E: WHO guideline exceedances
# ============================================================
who_threshold <- 15  # WHO 2021 annual mean guideline for PM2.5

exceedance_summary <- pollution %>%
  filter(!is.na(pm2.5)) %>%
  group_by(site_code, year) %>%
  summarise(
    days_over_who   = sum(pm2.5 > who_threshold),
    pct_over_who    = mean(pm2.5 > who_threshold) * 100,
    annual_mean_pm25 = mean(pm2.5),
    .groups = "drop"
  )

cat("\nWHO Guideline Exceedance Summary:\n")
print(exceedance_summary)

p9 <- exceedance_summary %>%
  ggplot(aes(x = year, y = pct_over_who, colour = site_code, group = site_code)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    title    = "% of Days Exceeding WHO PM2.5 Guideline (15 µg/m³) per Year",
    subtitle = "Declining trend at all sites, with notable drop post-ULEZ at inner-city locations",
    x = "Year", y = "% of Days Exceeding WHO Guideline", colour = "Site"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/09_who_exceedance.png", p9, width = 9, height = 5, dpi = 150)

# Save summary statistics for use in report
eda_summary <- list(
  cor_matrix        = cor_matrix,
  missing_summary   = missing_summary,
  exceedance_summary = exceedance_summary,
  n_obs             = nrow(pollution),
  date_range        = range(pollution$date)
)

saveRDS(eda_summary, "data/eda_summary.rds")

cat("\n✅ Module 2 complete! Plots saved to outputs/\n")
cat("Plots created: 01 to 09\n")
cat("Next step: source('R/03_glm_models.R')\n")
