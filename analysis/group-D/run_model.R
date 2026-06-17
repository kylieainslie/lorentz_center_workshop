library(cmdstanr)
library(dplyr)
library(ggplot2)
library(posterior)   # for draws summaries

# ---- Load pre-prepared Stan data -----------------------------------------
# Run prepare_stan_data.R first if these files don't exist
stan_data  <- readRDS("analysis/group-D/stan_data.rds")
stan_cells <- readRDS("analysis/group-D/stan_cells.rds")

cat("Data cells: ", stan_data$N, "\n")
cat("T_vax bins: ", stan_data$T_vax, "\n")

# ---- Compile Stan model --------------------------------------------------
model <- cmdstan_model("analysis/group-D/model.stan")

# ---- Run MCMC ------------------------------------------------------------
fit <- model$sample(
  data            = stan_data,
  seed            = 42,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 500,
  iter_sampling   = 1000,
  refresh         = 200,
  adapt_delta     = 0.95
)

# ---- Basic diagnostics ---------------------------------------------------
fit$cmdstan_diagnose()

print(fit$summary(c("log_beta", "sigma_gamma")))

# ---- Vaccine effectiveness over time -------------------------------------
ve_summary <- fit$summary("VE") |>
  mutate(
    tvax_bin = as.integer(sub("VE\\[(\\d+)\\]", "\\1", variable)),
    day_since_vax = tvax_bin - 1L   # 0 = unvaccinated bin, 1..90 = days since vax
  )

print(ve_summary |> select(variable, day_since_vax, median, q5, q95))

# ---- Plot VE over time ---------------------------------------------------
p_ve <- ve_summary |>
  filter(tvax_bin > 1) |>   # exclude unvaccinated bin
  ggplot(aes(x = day_since_vax, y = median)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.5, 1)) +
  labs(
    title    = "Time-varying vaccine effectiveness",
    subtitle = "Median and 90% credible interval",
    x        = "Days since vaccination",
    y        = "Vaccine effectiveness (1 - HR)"
  ) +
  theme_bw()

ggsave("analysis/group-D/ve_over_time.png", p_ve,
       width = 8, height = 5, dpi = 150)
print(p_ve)

# ---- Save fit for further analysis ---------------------------------------
fit$save_object("analysis/group-D/fit.rds")
cat("Model fit saved to analysis/group-D/fit.rds\n")
