# ============================================================
# 01_data_acquisition.R
# Download real UK air quality data using the openair package
# Data comes from the UK government's DEFRA AURN network
# ============================================================

library(openair)
library(tidyverse)
library(lubridate)

cat("=== Module 1: Data Acquisition ===\n")
cat("Downloading real UK air quality data from DEFRA AURN...\n")

# ---- 1. Download pollution data for multiple London sites ----
# We use multiple sites to show spatial variation (links to Design & Sampling)
# MY1  = Marylebone Road (roadside - high pollution)
# KC1  = Kings College (urban background)
# BT1  = Bloomsbury (central London background)
# HRL  = Harlington (outer London / Heathrow area)

sites <- c("MY1", "KC1", "BT1", "HRL")

cat("Downloading data for sites:", paste(sites, collapse = ", "), "\n")
cat("Time range: 2018-01-01 to 2023-12-31\n")

# Download each site separately for reliability
download_site <- function(site_code) {
  cat("  Fetching", site_code, "...\n")
  tryCatch({
    df <- importAURN(
      site   = site_code,
      year   = 2018:2023,
      pollutant = c("pm2.5", "no2", "o3", "nox", "ws", "wd", "temp"),
      meta   = TRUE
    )
    df$site_code <- site_code
    return(df)
  }, error = function(e) {
    cat("  Warning: Could not download", site_code, "- generating synthetic data\n")
    return(NULL)
  })
}

raw_list <- lapply(sites, download_site)

# Check what downloaded successfully
successful <- !sapply(raw_list, is.null)
cat(sum(successful), "of", length(sites), "sites downloaded successfully\n")

# ---- 2. If download fails (e.g. no internet), generate realistic synthetic data ----
generate_synthetic_pollution <- function(site_code, start = "2018-01-01", end = "2023-12-31") {
  
  cat("  Generating synthetic data for", site_code, "\n")
  
  dates <- seq(as.Date(start), as.Date(end), by = "day")
  n     <- length(dates)
  
  # Baseline levels differ by site type
  pm25_base <- switch(site_code,
    "MY1" = 18,   # Roadside — highest
    "KC1" = 12,
    "BT1" = 11,
    "HRL" = 9,    # Outer London — lowest
    12
  )
  
  set.seed(42 + which(sites == site_code))
  
  # Seasonal pattern: higher PM2.5 in winter (wood burning, inversions)
  day_of_year   <- as.numeric(format(dates, "%j"))
  seasonal      <- 5 * cos(2 * pi * (day_of_year - 15) / 365)
  
  # Weekly pattern: lower on weekends (less traffic)
  weekday       <- as.numeric(format(dates, "%u"))  # 1=Mon, 7=Sun
  weekly        <- ifelse(weekday >= 6, -2, 0)
  
  # Long-term downward trend (air quality improving over years)
  year_num      <- as.numeric(format(dates, "%Y")) - 2018
  trend         <- -0.5 * year_num
  
  # Random noise
  noise         <- rnorm(n, 0, 3)
  
  # Simulate ULEZ effect: ~15% reduction after Oct 2021 for inner London sites
  ulez_date     <- as.Date("2021-10-25")
  ulez_effect   <- ifelse(
    dates >= ulez_date & site_code %in% c("MY1", "KC1", "BT1"),
    -3,    # approximately 15-20% reduction
    0
  )
  
  pm25 <- pmax(1, pm25_base + seasonal + weekly + trend + ulez_effect + noise)
  
  # NO2 correlated with PM2.5 but noisier
  no2  <- pmax(5, pm25 * 2.8 + rnorm(n, 0, 8) + ifelse(dates >= ulez_date & site_code %in% c("MY1","KC1","BT1"), -5, 0))
  
  # O3 is inversely related to NO2 (chemistry!)
  o3   <- pmax(5, 45 - 0.3 * no2 + rnorm(n, 0, 6))
  
  # Temperature: realistic UK seasonal pattern
  temp <- 10 + 8 * sin(2 * pi * (day_of_year - 80) / 365) + rnorm(n, 0, 2)
  
  # Wind speed: log-normal
  ws   <- pmax(0.5, rlnorm(n, log(4), 0.5))
  
  data.frame(
    date      = dates,
    site      = site_code,
    site_code = site_code,
    pm2.5     = round(pm25, 2),
    no2       = round(no2, 2),
    o3        = round(o3, 2),
    nox       = round(no2 * 1.4, 2),
    temp      = round(temp, 1),
    ws        = round(ws, 1),
    wd        = sample(0:359, n, replace = TRUE),
    latitude  = switch(site_code, "MY1"=51.523, "KC1"=51.511, "BT1"=51.524, "HRL"=51.486, 51.5),
    longitude = switch(site_code, "MY1"=-0.154, "KC1"=-0.116, "BT1"=-0.126, "HRL"=-0.441, -0.1)
  )
}

# Build the final dataset: use real data where available, synthetic otherwise
final_list <- lapply(seq_along(sites), function(i) {
  if (successful[i]) {
    df <- raw_list[[i]]
    df$site_code <- sites[i]
    # Standardise column names from openair
    if ("date" %in% names(df)) {
      df$date <- as.Date(df$date)
    }
    return(df)
  } else {
    return(generate_synthetic_pollution(sites[i]))
  }
})

# Combine all sites
pollution_raw <- bind_rows(final_list)

# ---- 3. Standardise and clean column names ----
pollution_clean <- pollution_raw %>%
  select(
    date, site_code,
    pm2.5  = any_of(c("pm2.5", "PM2.5")),
    no2    = any_of(c("no2",   "NO2")),
    o3     = any_of(c("o3",    "O3")),
    temp   = any_of(c("temp",  "air_temperature", "temperature")),
    ws     = any_of(c("ws",    "wind_speed")),
    wd     = any_of(c("wd",    "wind_direction"))
  ) %>%
  # Ensure numeric
  mutate(across(c(pm2.5, no2, o3, temp, ws), as.numeric)) %>%
  # Remove physically impossible values
  mutate(
    pm2.5 = ifelse(pm2.5 < 0 | pm2.5 > 500, NA, pm2.5),
    no2   = ifelse(no2   < 0 | no2   > 500, NA, no2),
    o3    = ifelse(o3    < 0 | o3    > 300, NA, o3),
    temp  = ifelse(temp  < -20 | temp > 45, NA, temp)
  ) %>%
  arrange(site_code, date)

# ---- 4. Add derived variables useful for modelling ----
pollution_clean <- pollution_clean %>%
  mutate(
    year         = year(date),
    month        = month(date),
    day_of_week  = wday(date, label = TRUE, abbr = TRUE),
    is_weekend   = wday(date) %in% c(1, 7),
    season       = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5)  ~ "Spring",
      month %in% c(6, 7, 8)  ~ "Summer",
      TRUE                    ~ "Autumn"
    ),
    # ULEZ expanded in London on 25 October 2021
    post_ulez    = date >= as.Date("2021-10-25"),
    # Air quality index categories (WHO guidelines)
    pm25_category = cut(pm2.5,
      breaks = c(-Inf, 10, 20, 25, 50, Inf),
      labels = c("Good", "Moderate", "Unhealthy (sensitive)", "Unhealthy", "Very Unhealthy")
    )
  )

# ---- 5. Generate synthetic daily hospital respiratory admissions ----
# Based on published epidemiological dose-response relationships
# A 10 µg/m³ increase in PM2.5 increases admissions by ~1.5% (Bell et al. 2008)
cat("Generating synthetic health outcome data...\n")

set.seed(123)

# Average PM2.5 across all sites each day (city-wide exposure)
avg_pm25 <- pollution_clean %>%
  group_by(date) %>%
  summarise(city_pm25 = mean(pm2.5, na.rm = TRUE), .groups = "drop")

health_data <- avg_pm25 %>%
  mutate(
    year    = year(date),
    month   = month(date),
    
    # Baseline admissions: seasonal + day-of-week effects
    seasonal_effect = 15 + 10 * cos(2 * pi * (yday(date) - 15) / 365),
    
    # Pollution effect: log-linear, ~1.5% per 10 µg/m³
    pollution_rr    = exp(0.0015 * (city_pm25 - 12)),
    
    # Flu season spike (Jan-Feb)
    flu_spike       = ifelse(month %in% c(1, 2), rpois(n(), 5), 0),
    
    # Expected count
    mu              = (seasonal_effect * pollution_rr + flu_spike),
    
    # Observed: Poisson with overdispersion (negative binomial)
    respiratory_admissions = rnbinom(n(), mu = mu, size = 10),
    
    # ULEZ intervention column
    post_ulez = date >= as.Date("2021-10-25")
  ) %>%
  select(date, city_pm25, respiratory_admissions, seasonal_effect, post_ulez)

# ---- 6. Save processed data ----
dir.create("data", showWarnings = FALSE)

saveRDS(pollution_clean, "data/pollution_clean.rds")
saveRDS(health_data,     "data/health_data.rds")

write_csv(pollution_clean, "data/pollution_clean.csv")
write_csv(health_data,     "data/health_data.csv")

# ---- 7. Summary ----
cat("\n=== Data Summary ===\n")
cat("Pollution dataset:\n")
cat("  Rows:", nrow(pollution_clean), "\n")
cat("  Sites:", paste(unique(pollution_clean$site_code), collapse = ", "), "\n")
cat("  Date range:", format(min(pollution_clean$date)), "to", format(max(pollution_clean$date)), "\n")
cat("  Missing PM2.5:", sum(is.na(pollution_clean$pm2.5)), "days\n")

cat("\nHealth dataset:\n")
cat("  Rows:", nrow(health_data), "\n")
cat("  Mean admissions/day:", round(mean(health_data$respiratory_admissions), 1), "\n")
cat("  Correlation (PM2.5 ~ admissions):", round(cor(health_data$city_pm25, health_data$respiratory_admissions, use="complete.obs"), 3), "\n")

cat("\n✅ Module 1 complete! Data saved to data/ folder.\n")
cat("Next step: source('R/02_eda.R')\n")
