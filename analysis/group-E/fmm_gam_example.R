here::here("fmm_gam_example.R")

library(here)
library(mgcv)

# =====================================================================
# 1. SIMULATE DATA
# =====================================================================
set.seed(42)
n <- 1200

# Covariates
a <- runif(n, -2, 2)
b <- runif(n, -2, 2)

# True varying logit weight function: w ~ s(a, b)
logit_w <- sin(a) + cos(b)
w_true <- 1 / (1 + exp(-logit_w))

# Latent component assignment (Component 1 or 2)
z <- rbinom(n, 1, w_true)

# True Gaussian parameters
mu1_true <- -2
sd1_true <- 0.8
mu2_true <- 3
sd2_true <- 1.2

# Generate observed response y
y <- ifelse(z == 1, rnorm(n, mu1_true, sd1_true), rnorm(n, mu2_true, sd2_true))
df <- data.frame(y = y, a = a, b = b)

# =====================================================================
# 2. EM ALGORITHM FUNCTION
# =====================================================================
fit_gam_mixture <- function(y, a, b, max_iter = 100, tol = 1e-5) {
  n <- length(y)

  # Smart Initialization using K-means
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
      family = quasi(link = "logit", variance = "mu(1-mu)"),
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
# 3. RUN AND EVALUATE
# =====================================================================
model_fit <- fit_gam_mixture(df$y, df$a, df$b)

# Print results vs Truth
cat("\n--- Estimated Means ---\n")
print(model_fit$mu)
cat("True Means:", mu1_true, ",", mu2_true, "\n")

cat("\n--- Estimated SDs ---\n")
print(model_fit$sd)
cat("True SDs:", sd1_true, ",", sd2_true, "\n")

# Summary of the fitted GAM gates
summary(model_fit$gam_model)
