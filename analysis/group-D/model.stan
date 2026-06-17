// Cox proportional hazards model (Poisson approximation, person-day format)
//
// Hazard for individual i on calendar day t_c, time-since-vax t:
//   h(i) = exp(log_beta + log_Ibar[i] + gamma[t_bin[i]] * vax[i])
//
// gamma[t] is the log-hazard-ratio for vaccination, allowed to vary
// by day since vaccination via a random-walk prior.
// Negative gamma => protection; gamma = 0 => no effect.
//
// Aggregated likelihood: rows are (age_group × vax_status × day × tvax_bin)
// cells; Y = events, offset = log(Ibar) + log(n_persdays).

data {
  int<lower=1> N;                   // number of aggregated cells
  int<lower=1> K;                   // number of age groups
  int<lower=1> T_vax;               // number of time-since-vax bins (weeks)

  array[N] int<lower=0> Y;          // event count per cell
  vector[N] log_Ibar;               // pre-computed log force-of-infection offset
  vector[N] log_persdays;           // log(person-days at risk) per cell
  array[N] int<lower=0,upper=1> vax;          // vaccination indicator (0/1)
  array[N] int<lower=1,upper=T_vax> tvax_bin; // weekly bin since vax (1 = unvax/bin 1)
}

parameters {
  real log_beta;                    // log transmission scaling
  vector[T_vax] gamma;              // time-varying log-VE (hazard ratio); gamma < 0 = protection
  real<lower=0> sigma_gamma;        // random-walk SD
}

transformed parameters {
  vector[N] log_mu;
  for (n in 1:N) {
    log_mu[n] = log_beta + log_Ibar[n] + log_persdays[n]
                + gamma[tvax_bin[n]] * vax[n];
  }
}

model {
  // Priors
  log_beta    ~ normal(0, 2);
  sigma_gamma ~ normal(0, 0.5);

  // Random walk on gamma: first element anchored, rest follow random walk
  gamma[1] ~ normal(-1, 1);        // prior: ~63 % VE at t=0
  for (t in 2:T_vax)
    gamma[t] ~ normal(gamma[t-1], sigma_gamma);

  // Poisson likelihood (Poisson approximation to Cox)
  Y ~ poisson_log(log_mu);
}

generated quantities {
  // Vaccine effectiveness = 1 - exp(gamma)
  vector[T_vax] VE;
  for (t in 1:T_vax)
    VE[t] = 1 - exp(gamma[t]);
}
