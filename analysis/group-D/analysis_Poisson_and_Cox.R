source("analysis//group-D//data_wrangling.R")
table(table(cohort_tv_daily$ind_id))

# Add time difference and logI etc as offset variables
d <- cohort_tv_daily |>
  mutate(delta_t = tstop - tstart,
         logI = log(I),
         logIa1 = log(Ia1),
         logIa2 = log(Ia2),
         logIa3 = log(Ia3),
         logIa4 = log(Ia4),
         logIa5 = log(Ia5),
         logdelta_t = log(delta_t),
         tv = tstart - time_of_vaccination # time since vaccination (at start of interval)
  )

# We start off with the simplest Poisson model
# We estimate a single hazard ratio (vax_status),
# based on non-structured I(t), hence offset(logI) term,
# offset(logdelta_t) is not really necessary here, but would be if
# tstart - tstop intervals have unequal length
glm(event ~ vax_status + offset(logdelta_t) + offset(logI), family = "poisson", data = d)
# Cox model, not using SIR input
coxph(Surv(tstart, tstop, event) ~ vax_status, data = d)
# Adding offset(logI) doesn't do anything, other than changing the
# baseline hazard, which isn't shown in the output anyway
coxph(Surv(tstart, tstop, event) ~ vax_status + offset(logI), data = d)

# A straightforward way to estimate time-varying VE is through
# estimating piecewise constant hazard ratios on prespecified
# intervals
# Working with tv, which has still got negative values in the period
# before vaccination, set these to 0 to indicate the period before
# vaccination where vax_status = 0
d$tv[d$tv < 0] <- 0
d$tv[is.na(d$tv)] <- 0 # for those never vaccinated
d$period <- cut(d$tv, c(-Inf, 0, 4.99, 9.99, 19.99, 29.99, Inf))
table(d$period)

pois_pwc <- glm(event ~ vax_status : period + offset(logdelta_t) + offset(logI), family = "poisson", data = d)
summary(pois_pwc)
# I am a little unsure what to do with the vax_status:period(-Inf,0] term
# I am inclined to think that this still is part of a "baseline", no
# vaccination, so each of the next terms should be interpreted as log(HR)'s
betas <- pois_pwc$coefficients
VEs <- 1 - exp(betas[3:7])
VEs
# Very rough picture (should be much nicer with ggplot2 shaded areas)
# First extract SE's of coefficients
covmat <- summary(pois_pwc)$cov.unscaled
SEs <- sqrt(diag(covmat))
lower <- 1 - exp(betas[3:7] + qnorm(0.975) * SEs[3:7])
upper <- 1 - exp(betas[3:7] - qnorm(0.975) * SEs[3:7])
plot(c(0, 5), rep(VEs[1], 2), xlim = c(0, 45), ylim = c(-0.35, 1), type = "s",
     xlab = "Days since vaccination", ylab = "VE")
lines(c(0, 5), rep(lower[1], 2), lty = 3)
lines(c(0, 5), rep(upper[1], 2), lty = 3)
lines(c(5, 10), rep(VEs[2], 2))
lines(c(5, 10), rep(lower[2], 2), lty = 3)
lines(c(5, 10), rep(upper[2], 2), lty = 3)
lines(c(10, 20), rep(VEs[3], 2))
lines(c(10, 20), rep(lower[3], 2), lty = 3)
lines(c(10, 20), rep(upper[3], 2), lty = 3)
lines(c(20, 30), rep(VEs[4], 2))
lines(c(20, 30), rep(lower[4], 2), lty = 3)
lines(c(20, 30), rep(upper[4], 2), lty = 3)
lines(c(30, 45), rep(VEs[5], 2))
lines(c(30, 45), rep(lower[5], 2), lty = 3)
lines(c(30, 45), rep(upper[5], 2), lty = 3) # lower bound way less than 0
# Would like to try out piecewise linear estimates

# TODO: bring in the age structure, using logIa1 ... logIa5 and
