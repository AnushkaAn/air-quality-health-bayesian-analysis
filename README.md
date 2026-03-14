# 🏙️ Bayesian Urban Air Quality Forecasting & Health Risk Analysis

An end-to-end statistical analysis project in R covering Bayesian modelling, time series forecasting, causal inference, generalised regression models, and machine learning — built to demonstrate skills required for the MSc Statistics with Data Science at the University of Edinburgh.

---

## 📁 Project Structure

```
air-quality-project/
├── README.md
├── air_quality_project.Rproj
├── data/                          # Raw and processed data (auto-downloaded)
├── R/
│   ├── 00_setup.R                 # Install & load all packages
│   ├── 01_data_acquisition.R      # Download real UK pollution data
│   ├── 02_eda.R                   # Exploratory data analysis
│   ├── 03_glm_models.R            # Generalised regression models
│   ├── 04_time_series.R           # ARIMA & forecasting
│   ├── 05_bayesian_model.R        # Bayesian hierarchical model
│   ├── 06_causal_inference.R      # Interrupted time series (ULEZ)
│   └── 07_ml_comparison.R         # Random forest & XGBoost
├── reports/
│   ├── full_report.Rmd            # Complete R Markdown report
│   └── full_report.html           # (generated output)
└── shiny/
    └── app.R                      # Interactive Shiny dashboard
```

---

## 🚀 How to Run

### Step 1 — Clone and open the project
```r
# Open air_quality_project.Rproj in RStudio
```

### Step 2 — Install all packages
```r
source("R/00_setup.R")
```

### Step 3 — Run modules in order
```r
source("R/01_data_acquisition.R")   # Downloads real data (~2 min)
source("R/02_eda.R")
source("R/03_glm_models.R")
source("R/04_time_series.R")
source("R/05_bayesian_model.R")     # Takes 3–5 min (MCMC sampling)
source("R/06_causal_inference.R")
source("R/07_ml_comparison.R")
```

### Step 4 — Generate the full report
```r
rmarkdown::render("reports/full_report.Rmd")
```

### Step 5 — Launch the Shiny dashboard
```r
shiny::runApp("shiny/app.R")
```

---

## 📊 What This Project Covers

| Module | MSc Course Covered |
|---|---|
| Bayesian hierarchical model | Bayesian Data Analysis, Bayesian Theory |
| Poisson & negative binomial GLMs | Generalised Regression Models |
| ARIMA + ETS forecasting | Time Series |
| Interrupted time series (ULEZ) | Methods for Causal Inference |
| Random forest + XGBoost | Applied Machine Learning |
| Missing data imputation | Incomplete Data Analysis |
| tidyverse + ggplot2 + Shiny | Extended Statistical Programming |
| Spatial sampling strategy | Design and Sampling for Data Science |

---

## 🌍 Data Sources

- **Air quality**: `openair` R package — wraps UK DEFRA Automatic Urban and Rural Network (AURN) data
- **Intervention**: London ULEZ expansion (October 2021) as the causal event
- **Synthetic health outcomes**: Generated with realistic epidemiological parameters (Poisson-distributed, correlated with pollution levels)

---

## 🌐 Deployment

Deploy the Shiny dashboard for free at [shinyapps.io](https://www.shinyapps.io/):
```r
install.packages("rsconnect")
rsconnect::deployApp("shiny/")
```

---

## 📚 Key R Packages Used

`openair`, `tidyverse`, `ggplot2`, `forecast`, `brms`, `bayesplot`, `loo`, `ranger`, `xgboost`, `shiny`, `plotly`, `mice`, `sandwich`, `lmtest`

---

## 👤 Author

Built as a portfolio project.
