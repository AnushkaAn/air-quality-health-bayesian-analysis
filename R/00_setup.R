# ============================================================
# 00_setup.R
# Install and load all packages needed for the project
# Run this FIRST before anything else
# ============================================================

cat("=== Installing packages (this may take 5-10 mins first time) ===\n")

# --- List of all packages needed ---
packages <- c(
  # Data import & wrangling
  "tidyverse",      # ggplot2, dplyr, tidyr, readr etc.
  "lubridate",      # easy date handling
  "janitor",        # clean column names
  
  # Air quality data
  "openair",        # UK DEFRA pollution data (the key data source)
  
  # Missing data
  "mice",           # multiple imputation
  "VIM",            # visualise missing data
  
  # Statistical modelling
  "MASS",           # negative binomial regression
  "lme4",           # mixed effects models
  "sandwich",       # robust standard errors
  "lmtest",         # coefficient tests
  
  # Time series
  "forecast",       # ARIMA, ETS
  "tseries",        # unit root tests
  "zoo",            # rolling averages
  
  # Bayesian modelling
  "brms",           # Bayesian regression using Stan (main Bayesian package)
  "bayesplot",      # plot posterior distributions
  "loo",            # model comparison
  "tidybayes",      # tidy Bayesian outputs
  
  # Machine learning
  "ranger",         # fast random forests
  "xgboost",        # gradient boosting
  "caret",          # cross-validation helpers
  "vip",            # variable importance plots
  
  # Visualisation
  "ggplot2",        # (part of tidyverse, but explicit)
  "patchwork",      # combine multiple plots
  "scales",         # axis formatting
  "viridis",        # colour-blind friendly palettes
  "plotly",         # interactive plots
  "kableExtra",     # nice tables in reports
  
  # Reporting & app
  "rmarkdown",      # generate reports
  "knitr",          # R Markdown engine
  "shiny",          # interactive web app
  "shinydashboard", # dashboard layout for Shiny
  "DT",             # interactive tables in Shiny
  "rsconnect"       # deploy to shinyapps.io
)

# --- Install any that are missing ---
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]

if (length(to_install) > 0) {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, dependencies = TRUE)
} else {
  cat("All packages already installed!\n")
}

# --- Load all packages ---
cat("\nLoading packages...\n")
suppressPackageStartupMessages({
  lapply(packages, library, character.only = TRUE)
})

cat("\n✅ Setup complete! All packages loaded.\n")
cat("Next step: source('R/01_data_acquisition.R')\n")
