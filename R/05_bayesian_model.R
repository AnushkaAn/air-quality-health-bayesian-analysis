# ============================================================
# 05_bayesian_model.R
# Bayesian Hierarchical Model using brms (Bayesian Regression Models)
# This is the flagship statistical module — directly maps to
# Bayesian Data Analysis and Bayesian Theory courses at Edinburgh
# Note: MCMC sampling takes 3–5 minutes — normal and expected
# ============================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(loo)
library(tidybayes)
library(ggplot2)
library(patchwork)

cat("=== Module 5: Bayesian Hierarchical Model ===\n")
cat("Note: MCMC sampling will take 3-5 minutes. This is normal.\n\n")

# Load data
glm_results <- readRDS("data/glm_results.rds")
model_data  <- glm_results$model_data
pollution   <- readRDS("data/pollution_clean.rds")

# ============================================================
# SECTION A: Why Bayesian? — set the scene
# ============================================================
# Classical (frequentist) GLM gives a point estimate for each coefficient.
# Bayesian model gives us the FULL PROBABILITY DISTRIBUTION of each coefficient.
# This means we can say:
#   "There's a 94% probability that PM2.5 increases admissions by between 0.8% and 2.1%"
# instead of just:
#   "The coefficient is 0.0015, p < 0.05"
# Much more informative for public health decision-making.

# ============================================================
# SECTION B: Prior specification (Bayesian Theory)
# ============================================================
cat("Specifying prior distributions...\n")

# We use weakly informative priors based on epidemiological knowledge:
# - PM2.5 effect: we expect a small positive effect, centred on 0 with modest SD
#   Normal(0, 0.05) on log scale = we think RR is between 0.9 and 1.1 a priori
# - Intercept: Normal prior on log scale
# - Overdispersion: Half-Cauchy for the negative binomial shape

priors <- c(
  prior(normal(0, 0.05),  class = "b",         coef = "pm2.5_mean"),
  prior(normal(0, 0.05),  class = "b",         coef = "pm25_lag1"),
  prior(normal(0, 0.05),  class = "b",         coef = "pm25_lag2"),
  prior(normal(0, 0.05),  class = "b",         coef = "pm25_lag3"),
  prior(normal(0, 0.1),   class = "b",         coef = "temp_mean"),
  prior(normal(3, 0.5),   class = "Intercept"),
  prior(cauchy(0, 2.5),   class = "shape")
)

cat("Prior specification:\n")
cat("  PM2.5 effect:  Normal(0, 0.05) — weakly informative\n")
cat("  Intercept:     Normal(3, 0.5)  — baseline ~20 admissions/day\n")
cat("  Shape (NB):    Half-Cauchy(0, 2.5) — regularises overdispersion\n")

# ============================================================
# SECTION C: Fit Bayesian model with brms
# ============================================================
cat("\nFitting Bayesian model...\n")
cat("This runs 4 MCMC chains × 2000 iterations = 8000 total samples\n")

# Check if model already cached (saves time when re-running)
model_cache_path <- "data/bayesian_model.rds"

if (file.exists(model_cache_path)) {
  cat("Loading cached model (delete data/bayesian_model.rds to refit)...\n")
  bayes_model <- readRDS(model_cache_path)
} else {
  
  bayes_model <- brm(
    formula = respiratory_admissions ~ pm2.5_mean + pm25_lag1 + pm25_lag2 + pm25_lag3 +
                                       temp_mean + month + weekday,
    data    = model_data,
    family  = negbinomial(),    # Negative binomial for overdispersed counts
    prior   = priors,
    chains  = 4,                # 4 parallel MCMC chains
    iter    = 2000,             # 2000 iterations per chain (1000 warmup + 1000 sampling)
    warmup  = 1000,
    cores   = 4,                # Use 4 CPU cores in parallel
    seed    = 42,
    control = list(
      adapt_delta   = 0.95,     # Higher = safer but slower (reduces divergent transitions)
      max_treedepth = 12
    ),
    file    = "data/bayesian_model"   # Cache the fitted model
  )
  
}

cat("Model fitted!\n")

# ============================================================
# SECTION D: MCMC diagnostics — did the sampler converge?
# ============================================================
cat("\nChecking MCMC convergence...\n")

# R-hat: should be very close to 1.0 (< 1.01 is excellent)
# Bulk ESS and Tail ESS: should be > 400 (more = better)
model_summary <- summary(bayes_model)

cat("Convergence diagnostics (Rhat values — should all be ≤ 1.01):\n")
rhat_values <- rhat(bayes_model)
rhat_df     <- data.frame(
  parameter = names(rhat_values),
  rhat      = as.numeric(rhat_values)
)

pollution_params <- rhat_df %>%
  filter(str_detect(parameter, "pm2.5|lag|temp|Intercept|shape"))

print(pollution_params)

all_converged <- all(pollution_params$rhat < 1.01, na.rm = TRUE)
cat(ifelse(all_converged, "✓ All parameters converged!\n", "⚠ Some chains may not have converged\n"))

# ============================================================
# SECTION E: Posterior distribution plots
# ============================================================
cat("Plotting posterior distributions...\n")

# Extract posterior samples
posterior_draws <- as_draws_df(bayes_model)

# Parameters of interest
params_of_interest <- c("b_pm2.5_mean", "b_pm25_lag1", "b_pm25_lag2", "b_pm25_lag3", "b_temp_mean")

# Convert to Rate Ratios for interpretability (exponentiate the log-scale estimates)
rr_draws <- posterior_draws %>%
  select(all_of(params_of_interest)) %>%
  mutate(across(everything(), exp)) %>%  # Convert to Rate Ratios
  pivot_longer(everything(), names_to = "parameter", values_to = "rr") %>%
  mutate(parameter = case_when(
    parameter == "b_pm2.5_mean" ~ "PM2.5 (day 0)",
    parameter == "b_pm25_lag1"  ~ "PM2.5 (lag 1)",
    parameter == "b_pm25_lag2"  ~ "PM2.5 (lag 2)",
    parameter == "b_pm25_lag3"  ~ "PM2.5 (lag 3)",
    parameter == "b_temp_mean"  ~ "Temperature",
    TRUE ~ parameter
  ))

# Posterior density plots
p_posterior <- rr_draws %>%
  filter(str_detect(parameter, "PM2.5")) %>%
  ggplot(aes(x = rr, fill = parameter)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  facet_wrap(~parameter, scales = "free", ncol = 2) +
  scale_fill_viridis_d() +
  labs(
    title    = "Posterior Distributions of Rate Ratios (PM2.5 Effect)",
    subtitle = "Distribution to the RIGHT of 1 = increased admissions. Wider = more uncertainty.",
    x        = "Rate Ratio (RR)",
    y        = "Posterior density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("outputs/17_posterior_distributions.png", p_posterior, width = 10, height = 7, dpi = 150)

# ============================================================
# SECTION F: Credible intervals (Bayesian equivalent of confidence intervals)
# ============================================================
cat("Computing credible intervals...\n")

# 95% Highest Density Interval (HDI) — the shortest interval containing 95% of the posterior
cred_intervals <- rr_draws %>%
  group_by(parameter) %>%
  summarise(
    median_rr = median(rr),
    mean_rr   = mean(rr),
    lower_95  = quantile(rr, 0.025),
    upper_95  = quantile(rr, 0.975),
    lower_80  = quantile(rr, 0.10),
    upper_80  = quantile(rr, 0.90),
    prob_gt_1 = mean(rr > 1),          # P(RR > 1) — probability of harmful effect
    .groups   = "drop"
  )

cat("\nCredible Intervals and Posterior Probabilities:\n")
print(cred_intervals %>% mutate(across(where(is.numeric), ~round(., 4))))

# ---- Bayesian forest plot ----
p_bayes_forest <- ggplot(cred_intervals, aes(y = reorder(parameter, median_rr))) +
  geom_vline(xintercept = 1, colour = "gray50", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0, linewidth = 2,
                 colour = "#5B8DBE", alpha = 0.4) +
  geom_errorbarh(aes(xmin = lower_80, xmax = upper_80), height = 0, linewidth = 3,
                 colour = "#5B8DBE", alpha = 0.7) +
  geom_point(aes(x = median_rr), size = 4, colour = "#FF7043") +
  geom_text(aes(x = upper_95 + 0.002,
                label = paste0("P(RR>1)=", round(prob_gt_1*100, 1), "%")),
            hjust = 0, size = 3.5, colour = "gray40") +
  labs(
    title    = "Bayesian Rate Ratios: Effect on Respiratory Admissions",
    subtitle = "Thick bars = 80% credible interval; thin bars = 95% CI. Orange dot = posterior median.",
    x        = "Rate Ratio",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  xlim(NA, max(cred_intervals$upper_95) + 0.015)

ggsave("outputs/18_bayesian_forest.png", p_bayes_forest, width = 10, height = 5, dpi = 150)

# ============================================================
# SECTION G: Posterior predictive checks
# ============================================================
cat("Performing posterior predictive checks...\n")

# These check if the model can reproduce the data it was trained on
pp_check_plot <- pp_check(bayes_model, ndraws = 100, type = "dens_overlay")

ggsave("outputs/19_pp_check.png", pp_check_plot, width = 8, height = 5, dpi = 150)

# Zero-inflation check
pp_zeros <- pp_check(bayes_model, ndraws = 100, type = "stat", stat = "mean")
ggsave("outputs/20_pp_check_mean.png", pp_zeros, width = 7, height = 5, dpi = 150)

# ============================================================
# SECTION H: Model comparison using LOO (Leave-One-Out cross-validation)
# ============================================================
cat("Computing LOO-CV for model comparison...\n")

# Compare: model with lags vs model without lags
bayes_simple_path <- "data/bayesian_simple.rds"

if (file.exists(bayes_simple_path)) {
  bayes_simple <- readRDS(bayes_simple_path)
} else {
  bayes_simple <- brm(
    formula = respiratory_admissions ~ pm2.5_mean + temp_mean + month + weekday,
    data    = model_data,
    family  = negbinomial(),
    prior   = priors[1:4],   # subset of priors
    chains  = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.95),
    file    = "data/bayesian_simple"
  )
}

# LOO comparison (lower ELPD = better predictive performance)
loo_full   <- loo(bayes_model)
loo_simple <- loo(bayes_simple)

loo_comparison <- loo_compare(loo_full, loo_simple)
cat("\nLOO Model Comparison:\n")
print(loo_comparison)

# ============================================================
# SECTION I: Prior vs Posterior comparison (visualise Bayesian learning)
# ============================================================
cat("Comparing prior vs posterior (showing what data taught us)...\n")

# Prior: Normal(0, 0.05) on log scale → on RR scale
set.seed(42)
n_prior <- 4000
prior_rr <- exp(rnorm(n_prior, 0, 0.05))

posterior_rr <- filter(rr_draws, parameter == "PM2.5 (day 0)")$rr

p_prior_posterior <- bind_rows(
  data.frame(rr = prior_rr,     type = "Prior"),
  data.frame(rr = posterior_rr, type = "Posterior")
) %>%
  ggplot(aes(x = rr, fill = type)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "gray40") +
  scale_fill_manual(values = c("Prior" = "#5B8DBE", "Posterior" = "#FF7043")) +
  coord_cartesian(xlim = c(0.9, 1.15)) +
  labs(
    title    = "Prior vs Posterior: PM2.5 Rate Ratio",
    subtitle = "The posterior is the updated belief AFTER seeing data — the core Bayesian idea",
    x        = "Rate Ratio for same-day PM2.5", y = "Density", fill = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/21_prior_vs_posterior.png", p_prior_posterior, width = 9, height = 5, dpi = 150)

# ============================================================
# Save results
# ============================================================
bayes_results <- list(
  model         = bayes_model,
  cred_intervals = cred_intervals,
  rhat_df        = rhat_df,
  loo_comparison = loo_comparison
)

saveRDS(bayes_results, "data/bayes_results.rds")

cat("\n=== Bayesian Model Summary ===\n")
cat("Key finding: PM2.5 effect on admissions:\n")
pm25_ci <- filter(cred_intervals, parameter == "PM2.5 (day 0)")
cat("  Posterior median RR:", round(pm25_ci$median_rr, 4), "\n")
cat("  95% credible interval:", round(pm25_ci$lower_95, 4), "to", round(pm25_ci$upper_95, 4), "\n")
cat("  P(harmful effect) =", round(pm25_ci$prob_gt_1 * 100, 1), "%\n")

cat("\n✅ Module 5 complete! Plots saved to outputs/\n")
cat("Next step: source('R/06_causal_inference.R')\n")
