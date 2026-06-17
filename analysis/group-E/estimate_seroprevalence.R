here::here("estimate_seroprevalence.R")

library(here)
library(mgcv)
library(ggplot2)

# =====================================================================
# 1. DATA PREP
# =====================================================================

sero <- readRDS(here::here("data", "sero_data.rds"))
demo <- readRDS(here::here("data", "demographic_data.rds"))
inc <- readRDS(here::here("data", "incidence_data_reduced.rds"))

max(unique(inc$time_of_reporting))

demo_summary <- demo |>
  dplyr::mutate(
    age_cat = cut(
      age,
      breaks = seq(0, 90, by = 10),
      labels = seq(10, 90, by = 10)
    )
  ) |>
  dplyr::group_by(age_cat) |>
  dplyr::summarise(count = dplyr::n())

inc_summary <- inc |>
  dplyr::left_join(demo, by = "ind_id") |>
  dplyr::select(ind_id, age, time_of_reporting) |>
  dplyr::mutate(
    age_cat = cut(
      age,
      breaks = seq(0, 90, by = 10),
      labels = seq(10, 90, by = 10)
    ),
    reporting_cat = cut(
      time_of_reporting,
      breaks = seq(0, 100, by = 10),
      labels = seq(10, 100, by = 10)
    )
  ) |>
  dplyr::group_by(reporting_cat, age_cat) |>
  dplyr::summarise(total_cases = dplyr::n()) |>
  dplyr::left_join(demo_summary, by = "age_cat") |>
  dplyr::group_by(reporting_cat, age_cat, total_cases) |>
  dplyr::summarise(prev = total_cases / count)

df <- dplyr::left_join(sero, demo, by = "ind_id") |>
  dplyr::select(ind_id, age, sex, log_igg, time_of_sampling) |>
  dplyr::mutate(
    igg = exp(log_igg),
    log_age = log(age),
    log_time_of_sampling = log(time_of_sampling)
  ) |>
  dplyr::arrange(time_of_sampling, age)

n_samples <- length(unique(df$time_of_sampling))
sampling_times <- unique(df$time_of_sampling)

# =====================================================================
# 2. EM ALGORITHM FUNCTION
# =====================================================================
fit_gam_mixture <- function(
  y,
  a,
  b,
  link = "logit",
  max_iter = 100,
  tol = 1e-5
) {
  n <- length(y)

  # Data-based Initialization using K-means
  km <- kmeans(y, centers = 2)
  mu1 <- km$centers[1, 1]
  mu2 <- km$centers[2, 1]
  sd1 <- sd(y[km$cluster == 1])
  sd2 <- sd(y[km$cluster == 2])

  # Start with uniform weights
  w1 <- rep(0.5, n)

  loglik_old <- -Inf

  for (iter in 1:max_iter) {
    # ------------------ E-STEP ------------------
    dens1 <- dnorm(y, mean = mu1, sd = sd1)
    dens2 <- dnorm(y, mean = mu2, sd = sd2)

    # Calculate responsibilities
    gamma1 <- (w1 * dens1) / (w1 * dens1 + (1 - w1) * dens2)
    gamma2 <- 1 - gamma1

    # Monitor convergence via Log-Likelihood
    loglik <- sum(log(w1 * dens1 + (1 - w1) * dens2))
    if (abs(loglik - loglik_old) < tol) {
      cat("Converged at iteration:", iter, "\n")
      break
    }
    loglik_old <- loglik

    # ------------------ M-STEP ------------------
    # Update Gaussian Means
    mu1 <- sum(gamma1 * y) / sum(gamma1)
    mu2 <- sum(gamma2 * y) / sum(gamma2)

    # Update Gaussian Standard Deviations
    sd1 <- sqrt(sum(gamma1 * (y - mu1)^2) / sum(gamma1))
    sd2 <- sqrt(sum(gamma2 * (y - mu2)^2) / sum(gamma2))

    # Update GAM for the dynamic weight
    # Note: Using quasi family with logit link handles fractional pseudo-probabilities smoothly
    gam_mod <- gam(
      gamma1 ~ s(a, b),
      family = quasi(link = link, variance = "mu(1-mu)"),
      method = "REML"
    )
    w1 <- predict(gam_mod, type = "response")
  }

  # Return final parameters and estimated GAM model
  return(list(
    mu = c(Comp1 = mu1, Comp2 = mu2),
    sd = c(Comp1 = sd1, Comp2 = sd2),
    gam_model = gam_mod,
    loglik = loglik
  ))
}

# =====================================================================
# 3. RUN MODEL
# =====================================================================
model_fit <- fit_gam_mixture(
  df$log_igg,
  df$log_age,
  df$log_time_of_sampling,
  link = "cloglog"
)

# Print results vs Truth
cat("\n--- Estimated Means ---\n")
print(model_fit$mu)

cat("\n--- Estimated SDs ---\n")
print(model_fit$sd)

# Summary of the fitted GAM gates
summary(model_fit$gam_model)

# =====================================================================
# 4. VISUALIZE RESULTS
# =====================================================================

hist(df$log_igg, breaks = 50, xlab = "log(igg)", probability = TRUE)
abline(v = model_fit$mu, col = c("blue", "red"), lwd = 2)

df[, c("seroprev", "seroprev_se")] <- predict(
  model_fit$gam_mod,
  type = "response",
  se.fit = TRUE
)
if (model_fit$mu[1] < model_fit$mu[2]) {
  df$seroprev <- 1 - df$seroprev
}

df$incidence <- NA

for (i in 1:n_samples) {
  newd <- data.frame(
    log_age = log(seq(1, 90, by = 0.1)),
    log_time_of_sampling = log(sampling_times[i])
  )
  X0 <- predict(model_fit$gam_mod, newdata = newd, type = "lpmatrix")
}

# inc_summary |>
#   dplyr::filter(reporting_cat %in% as.character(seq(10, 50, by = 10))) |>
# ggplot() +
#   geom_line(
#     data = df,
#     aes(age, seroprev, group = 1),
#     size = 1
#   ) +
#   geom_bar(
#     data = inc_summary,
#     aes(x = age_cat, y = prev, fill = reporting_cat),
#     stat = "identity",
#     position = "dodge"
#   ) +
#   labs(x = "Age Category", y = "Prevalence", fill = "Reporting Time") +
#   theme_minimal()

ggplot() +
  geom_bar(
    data = subset(
      inc_summary,
      reporting_cat %in% as.character(seq(10, 40, by = 10))
    ),
    aes(x = as.numeric(paste(age_cat)), y = prev, fill = reporting_cat),
    stat = "identity",
    position = "dodge"
  ) +
  # scale_y_continuous() +
  geom_line(
    data = df,
    aes(
      age,
      seroprev,
      group = time_of_sampling,
      color = factor(time_of_sampling)
    ),
    size = 1,
    show.legend = FALSE
  ) +
  labs(
    x = "Age",
    y = "Estimated seroprevalence",
    fill = "Time from baseline"
  ) +
  theme_minimal()

# with(
#   df,
#   plot(
#     age[time_of_sampling == sampling_times[1]],
#     seroprev[time_of_sampling == sampling_times[1]],
#     xlab = "Age",
#     ylab = "Estimated Seroprevalence",
#     type = "l",
#     ylim = range(seroprev)
#   )
# )
# for (i in 2:n_samples) {
#   with(
#     df,
#     lines(
#       age[time_of_sampling == sampling_times[i]],
#       seroprev[time_of_sampling == sampling_times[i]],
#       col = i
#     )
#   )
# }
# legend(
#   "topright",
#   legend = c("Sampling time", sampling_times),
#   col = c(NA, 1:n_samples),
#   lty = c(NA, rep(1, 4)),
#   cex = 0.8
# )

ggplot() +
  geom_bar(
    data = subset(
      inc_summary,
      reporting_cat %in% as.character(seq(20, 50, by = 10))
    ),
    aes(x = as.numeric(paste(age_cat)), y = total_cases, fill = reporting_cat),
    stat = "identity",
    position = "dodge"
  )
