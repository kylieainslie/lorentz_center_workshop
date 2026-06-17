library(dplyr)
library(cmdstanr)
library(posterior)
library(ggplot2)

# Simulate event counts under KNOWN parameter values, using the REAL
# covariate structure (Ibar, person-days, vax, tvax_bin) produced by
# prepare_stan_data.R. Only Y (the outcome) is replaced with simulated
# counts. We then refit model.stan and check whether log_beta and the
# time-varying gamma(t) curve are recovered.
#
# Data-generating model (matches model.stan exactly):
#   log_mu = log_beta_true + log_Ibar + log_persdays + gamma_true[tvax_bin] * vax
#   Y ~ Poisson(exp(log_mu))

set.seed(42)

# ---- Load real covariate structure -----------------------------------------
# Running the real pipeline means log_Ibar, log_persdays, vax, and tvax_bin
# are identical to what the actual model fit sees; stan_data$Y is overwritten
# below with simulated counts.
source("analysis/group-D/prepare_stan_data.R")

N     <- stan_data$N
T_vax <- stan_data$T_vax

# ---- True parameters (data-generating model) -------------------------------

log_beta_true <- -4.4   # near the real-data fitted value (~-4.36) so simulated
                        # event counts are in a realistic regime

# Rise-then-wane VE curve: gamma = 0 at the moment of vaccination, becomes
# most negative (most protective) at t_peak days post-vaccination, then
# wanes back toward 0 over the remainder of the study.
t_peak <- 14
g_max  <- 1.2   # peak |gamma|, i.e. ~70% VE at peak (1 - exp(-1.2) = 0.70)

day_since_vax <- 0:(T_vax - 1)
gamma_true    <- -g_max * (day_since_vax / t_peak) * exp(1 - day_since_vax / t_peak)

# ---- Simulate event counts per cell -----------------------------------------

log_mu_true <- log_beta_true +
  stan_data$log_Ibar +
  stan_data$log_persdays +
  gamma_true[stan_data$tvax_bin] * stan_data$vax

Y_sim <- rpois(N, lambda = exp(log_mu_true))

cat("\nSimulated total events:", sum(Y_sim), "\n")
cat("Real total events:     ", sum(stan_data$Y), "\n")

stan_data_sim    <- stan_data
stan_data_sim$Y  <- as.integer(Y_sim)

# ---- Fit the real model to the simulated data -------------------------------

model <- cmdstan_model("analysis/group-D/model_noncentered.stan")

fit_sim <- model$sample(
  data            = stan_data_sim,
  seed            = 42,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 200,
  adapt_delta     = 0.99,
  max_treedepth   = 12
)

fit_sim$cmdstan_diagnose()

cat("\nlog_beta_true:    ", log_beta_true, "\n")
print(fit_sim$summary("log_beta"))

# ---- Compare recovered gamma(t) to truth -------------------------------------

gamma_summary <- fit_sim$summary("gamma") |>
  mutate(
    tvax_bin      = as.integer(sub("gamma\\[(\\d+)\\]", "\\1", variable)),
    day_since_vax = tvax_bin - 1L,
    gamma_true    = gamma_true[tvax_bin]
  )

p_recovery <- gamma_summary |>
  ggplot(aes(x = day_since_vax)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = median), colour = "steelblue", linewidth = 1) +
  geom_line(aes(y = gamma_true), colour = "firebrick", linetype = "dashed", linewidth = 1) +
  labs(
    title    = "Parameter recovery: gamma(t)",
    subtitle = "Blue = posterior median (90% CrI); red dashed = true simulated value",
    x        = "Days since vaccination",
    y        = "gamma (log hazard ratio)"
  ) +
  theme_bw()

ggsave("analysis/group-D/simulate_validate_recovery.png", p_recovery,
       width = 8, height = 5, dpi = 150)
print(p_recovery)

fit_sim$save_object("analysis/group-D/fit_sim.rds")
cat("Saved simulation fit to analysis/group-D/fit_sim.rds\n")
