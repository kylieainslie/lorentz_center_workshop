library(cmdstanr)
library(dplyr)
library(ggplot2)
library(posterior)   # for draws summaries

# ---- Load pre-prepared Stan data -----------------------------------------
# Run prepare_stan_data.R first if these files don't exist
stan_data        <- readRDS("analysis/group-D/stan_data.rds")
stan_data_simple <- readRDS("analysis/group-D/stan_data_simple.rds")

cat("Age-structured cells:  ", stan_data$N, "\n")
cat("No-age-structure cells:", stan_data_simple$N, "\n")
cat("T_vax bins:            ", stan_data$T_vax, "\n")

# ---- Compile Stan model (shared by both fits) ----------------------------
model <- cmdstan_model("analysis/group-D/model.stan")

# ---- Helper: extract VE summary from a fit object ------------------------
extract_ve <- function(fit, label) {
  fit$summary("VE") |>
    mutate(
      tvax_bin      = as.integer(sub("VE\\[(\\d+)\\]", "\\1", variable)),
      day_since_vax = tvax_bin - 1L,
      model         = label
    ) |>
    filter(tvax_bin > 1)   # drop unvaccinated bin
}

mcmc_args <- list(
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 100,
  iter_sampling   = 1000,
  refresh         = 200,
  adapt_delta     = 0.95
)

# ---- Fit age-structured model --------------------------------------------
cat("\n--- Fitting age-structured model ---\n")
fit_age <- do.call(model$sample, c(list(data = stan_data, seed = 42), mcmc_args))
fit_age$cmdstan_diagnose()
cat("\nAge-structured scalar parameters:\n")
print(fit_age$summary(c("log_beta", "sigma_gamma")))
fit_age$save_object("analysis/group-D/fit_age.rds")

# ---- Fit no-age-structure model ------------------------------------------
cat("\n--- Fitting no-age-structure model ---\n")
fit_simple <- do.call(model$sample, c(list(data = stan_data_simple, seed = 42), mcmc_args))
fit_simple$cmdstan_diagnose()
cat("\nNo-age-structure scalar parameters:\n")
print(fit_simple$summary(c("log_beta", "sigma_gamma")))
fit_simple$save_object("analysis/group-D/fit_simple.rds")

# ---- Combine VE summaries ------------------------------------------------
ve_combined <- bind_rows(
  extract_ve(fit_age,    "Age-structured FoI"),
  extract_ve(fit_simple, "Overall prevalence (no age structure)")
)

print(ve_combined |> select(model, day_since_vax, median, q5, q95))

# ---- Plot: VE comparison -------------------------------------------------
p_ve <- ggplot(ve_combined, aes(x = day_since_vax, y = median,
                                colour = model, fill = model)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Age-structured FoI"                   = "steelblue",
                                 "Overall prevalence (no age structure)" = "tomato")) +
  scale_fill_manual(values   = c("Age-structured FoI"                   = "steelblue",
                                 "Overall prevalence (no age structure)" = "tomato")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.5, 1), xlim = c(0, 45)) +
  labs(
    title    = "Time-varying vaccine effectiveness",
    subtitle = "Median and 90% credible interval",
    x        = "Days since vaccination",
    y        = "Vaccine effectiveness (1 - HR)",
    colour   = NULL, fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("analysis/group-D/ve_over_time.png", p_ve,
       width = 8, height = 5, dpi = 150)
print(p_ve)

cat("Plot saved to analysis/group-D/ve_over_time.png\n")
