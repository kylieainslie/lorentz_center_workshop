# =====================================================================
# ANALYSIS FOR EXCESS MORTALITY
# Author: Miguel Beynaerts
# Email: miguel.beynaerts@uhasselt.be
# Last modified: 17/06/2026 13:52:05
# License: MIT
#
# Description:
# Helper file for excess_mort.R
# ======================================================================

save_plot <- function(filename, plot_object, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

write_csv <- function(data, filename) {
  utils::write.csv(
    data,
    file = file.path(results_dir, filename),
    row.names = FALSE
  )
}

make_base_date <- function(anchor_year) {
  as.Date(sprintf("%s-01-01", anchor_year))
}

clean_strata_label <- function(x) {
  sub(".*=", "", x)
}

relabel_model_term <- function(term) {
  dplyr::case_when(
    term == "sexmale" ~ "Male vs female",
    term == "vaccinated_td" ~ "Vaccinated (time-dependent)",
    term == "school_memberYes" ~ "School member: Yes vs No",
    term == "workerYes" ~ "Worker: Yes vs No",
    term %in% c("age_band5-17", "age_band5.17") ~ "Age 5-17 vs 0-4",
    term %in% c("age_band18-39", "age_band18.39") ~ "Age 18-39 vs 0-4",
    term %in% c("age_band40-64", "age_band40.64") ~ "Age 40-64 vs 0-4",
    term %in% c("age_band65-85", "age_band65.85") ~ "Age 65-85 vs 0-4",
    grepl("^fu \\[", term) ~ paste("Baseline excess hazard", term),
    TRUE ~ term
  )
}

tidy_surv_object <- function(fit, label = NULL) {
  fit_summary <- summary(fit)

  if (length(fit_summary$time) == 0) {
    return(tibble::tibble())
  }

  strata_values <- if (is.null(fit_summary$strata)) {
    rep(label %||% "Overall", length(fit_summary$time))
  } else {
    as.character(fit_summary$strata)
  }

  lower <- if (is.null(fit_summary$lower)) {
    rep(NA_real_, length(fit_summary$time))
  } else {
    fit_summary$lower
  }

  upper <- if (is.null(fit_summary$upper)) {
    rep(NA_real_, length(fit_summary$time))
  } else {
    fit_summary$upper
  }

  base_rows <- tibble::tibble(
    time = 0,
    survival = 1,
    lower = 1,
    upper = 1,
    strata = unique(strata_values)
  )

  bind_rows(
    base_rows,
    tibble::tibble(
      time = fit_summary$time,
      survival = fit_summary$surv,
      lower = lower,
      upper = upper,
      strata = strata_values
    )
  ) |>
    distinct(strata, time, .keep_all = TRUE)
}

overall_curve_data <- function(data, poptab, anchor_year) {
  study_data <- data |>
    mutate(study_date = make_base_date(anchor_year))

  observed_fit <- survfit(Surv(followup_days, event) ~ 1, data = study_data)
  expected_fit <- survexp(
    Surv(followup_days, event) ~ 1,
    data = study_data,
    ratetable = poptab,
    rmap = list(age = age_days, year = study_date, sex = sex)
  )

  observed_summary <- summary(observed_fit, times = 1:90, extend = TRUE)
  expected_summary <- summary(expected_fit, times = 1:90, extend = TRUE)

  bind_rows(
    tibble::tibble(
      time = 0:90,
      measure = "Observed",
      survival = c(1, observed_summary$surv)
    ),
    tibble::tibble(
      time = 0:90,
      measure = "Expected",
      survival = c(1, expected_summary$surv)
    )
  ) |>
    mutate(cumulative_mortality = 1 - survival)
}

# Calculate the expected survival probability and hazard for each individual in the observed dataset
# based on the population lifetime table
add_expected_risk <- function(data, poptab, anchor_year) {
  study_data <- data |>
    mutate(study_date = make_base_date(anchor_year))

  study_data |>
    mutate(
      expected_survival = survival::survexp(
        followup_days ~ 1,
        data = study_data,
        method = "individual.s",
        ratetable = poptab,
        rmap = list(age = age_days, year = study_date, sex = sex)
      ),
      expected_hazard = survival::survexp(
        followup_days ~ 1,
        data = study_data,
        method = "individual.h",
        ratetable = poptab,
        rmap = list(age = age_days, year = study_date, sex = sex)
      ),
      expected_death_prob = 1 - expected_survival
    )
}

summarise_obs_expected <- function(data, group_vars) {
  data |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      n = n(),
      observed_deaths = sum(event),
      expected_deaths = sum(expected_death_prob),
      observed_mortality = mean(event),
      expected_mortality = mean(expected_death_prob),
      excess_per_1000 = (observed_mortality - expected_mortality) * 1000,
      oe_ratio = observed_deaths / expected_deaths,
      .groups = "drop"
    )
}

# Create start-stop times and covariate for vaccination status for a counting-process additive hazard model.
# This is extremely useful when modelling the time-dependent nature of vaccination
# Unvaccinated -> 1 row per individual with start = 0, stop = 90, event = 0/1, vaccinated_td = 0
# Vaccinated -> 2 rows per individual:
# ROW 1: start = 0, stop = time_of_vaccination, event = 0 (the individual lives long enough to be vaccinated),
# vaccinated_td = 0 (not vaccinated in this interval)
# ROW 2: start = time_of_vaccination, stop = 90, event = 0/1, vaccinated_td = 1 (since vaccinated in this interval)
# TODO: Could be extended to model lockdown effect after day 20 as well. Very interesting for final presentation!
make_split_dataset <- function(data, anchor_year) {
  # Subset of individuals who are not vaccinated during the epidemic
  never_or_late <- data |>
    filter(
      is.na(time_of_vaccination) |
        time_of_vaccination <= 0 |
        time_of_vaccination >= followup_days
    ) |>
    transmute(
      ind_id,
      age,
      age_band,
      sex,
      school_member,
      worker,
      start = 0,
      stop = followup_days,
      event,
      vaccinated_td = as.integer(
        !is.na(time_of_vaccination) & time_of_vaccination <= 0
      ),
      age_days = age * days_per_year
    )

  # Subset of individuals who are vaccinated during the epidemic (after day 45)
  pre_vaccination <- data |>
    filter(
      !is.na(time_of_vaccination),
      time_of_vaccination > 0,
      time_of_vaccination < followup_days
    )

  #
  before_vaccination <- pre_vaccination |>
    transmute(
      ind_id,
      age,
      age_band,
      sex,
      school_member,
      worker,
      start = 0,
      stop = time_of_vaccination,
      event = 0L,
      vaccinated_td = 0L,
      age_days = age * days_per_year
    )

  after_vaccination <- pre_vaccination |>
    transmute(
      ind_id,
      age,
      age_band,
      sex,
      school_member,
      worker,
      start = time_of_vaccination,
      stop = followup_days,
      event,
      vaccinated_td = 1L,
      age_days = age * days_per_year + time_of_vaccination
    )

  bind_rows(never_or_late, before_vaccination, after_vaccination) |>
    mutate(
      interval_days = stop - start,
      study_date = make_base_date(anchor_year) + start
    ) |>
    arrange(ind_id, start, stop)
}

# Calculate the expected survival probability and hazard for each individual in the observed dataset
# during the followup interval
# ie. What is the expected survival probability for an individual in the interval stop-start
add_expected_interval_risk <- function(data, poptab) {
  data |>
    mutate(
      expected_interval_survival = survexp(
        interval_days ~ 1,
        data = data,
        method = "individual.s",
        ratetable = poptab,
        rmap = list(age = age_days, year = study_date, sex = sex)
      ),
      expected_interval_death_prob = 1 - expected_interval_survival
    )
}

fit_rs_models <- function(split_data, poptab) {
  list(
    no_covariates = rsadd(
      Surv(start, stop, event) ~ 1,
      data = split_data,
      ratetable = poptab,
      int = baseline_breaks_years,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    minimal = rsadd(
      Surv(start, stop, event) ~ age_band + sex + vaccinated_td,
      data = split_data,
      ratetable = poptab,
      int = baseline_breaks_years,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    extended = rsadd(
      Surv(start, stop, event) ~ age_band +
        sex +
        vaccinated_td +
        school_member +
        worker,
      data = split_data,
      ratetable = poptab,
      int = baseline_breaks_years,
      rmap = list(age = age_days, year = study_date, sex = sex)
    )
  )
}

extract_model_coefficients <- function(model, model_name, anchor_year) {
  coefficient_table <- summary(model)$coefficients
  coefficients <- coefficient_table[, "Estimate"]
  standard_errors <- coefficient_table[, "Std. Error"]
  z_values <- coefficient_table[, "z value"]
  p_values <- coefficient_table[, "Pr(>|z|)"]

  tibble::tibble(
    model_name = model_name,
    anchor_year = anchor_year,
    term = names(coefficients),
    term_label = vapply(names(coefficients), relabel_model_term, character(1)),
    component = ifelse(
      grepl("^fu \\[", names(coefficients)),
      "baseline",
      "covariate"
    ),
    estimate = unname(coefficients),
    std_error = unname(standard_errors),
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error,
    hazard_ratio = exp(estimate),
    hazard_ratio_low = exp(conf_low),
    hazard_ratio_high = exp(conf_high),
    z_value = unname(z_values),
    p_value = unname(p_values)
  )
}

extract_model_overview <- function(model, model_name, anchor_year, split_data) {
  coefficient_table <- summary(model)$coefficients

  tibble::tibble(
    model_name = model_name,
    anchor_year = anchor_year,
    formula = paste(deparse(formula(model)), collapse = " "),
    n_rows = nrow(split_data),
    n_individuals = dplyr::n_distinct(split_data$ind_id),
    n_events = sum(split_data$event),
    all_coefficients_finite = all(is.finite(stats::coef(model))),
    all_confidence_limits_finite = all(is.finite(coefficient_table[, c(
      "Estimate",
      "Std. Error"
    )]))
  )
}

fit_rs_curves <- function(data, poptab, anchor_year) {
  study_data <- data |>
    mutate(study_date = make_base_date(anchor_year))

  list(
    overall = rs.surv(
      Surv(followup_days, event) ~ 1,
      data = study_data,
      ratetable = poptab,
      add.times = 1:90,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    age_band = rs.surv(
      Surv(followup_days, event) ~ age_band,
      data = study_data,
      ratetable = poptab,
      add.times = 1:90,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    sex = rs.surv(
      Surv(followup_days, event) ~ sex,
      data = study_data,
      ratetable = poptab,
      add.times = 1:90,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    school_member = rs.surv(
      Surv(followup_days, event) ~ school_member,
      data = study_data,
      ratetable = poptab,
      add.times = 1:90,
      rmap = list(age = age_days, year = study_date, sex = sex)
    ),
    worker = rs.surv(
      Surv(followup_days, event) ~ worker,
      data = study_data,
      ratetable = poptab,
      add.times = 1:90,
      rmap = list(age = age_days, year = study_date, sex = sex)
    )
  )
}

baseline_cumulative_hazard <- function(model, times) {
  stats::approx(
    x = model$times,
    y = model$Lambda0,
    xout = times,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

covariate_linear_predictor <- function(model, newdata) {
  rhs_terms <- delete.response(stats::terms(model))
  model_matrix <- stats::model.matrix(rhs_terms, data = newdata)
  coefficients <- stats::coef(model)
  covariate_coefficients <- coefficients[!grepl("^fu \\[", names(coefficients))]
  colnames(model_matrix) <- make.names(colnames(model_matrix))
  names(covariate_coefficients) <- make.names(names(covariate_coefficients))

  common_terms <- intersect(
    colnames(model_matrix),
    names(covariate_coefficients)
  )
  sum(
    model_matrix[, common_terms, drop = FALSE] *
      covariate_coefficients[common_terms]
  )
}

predict_profile_curve <- function(
  model,
  poptab,
  anchor_year,
  age_years,
  age_band,
  sex,
  school_member,
  worker,
  vaccination_day,
  profile_label,
  times = profile_times
) {
  base_times <- sort(unique(times))
  base_lambda <- baseline_cumulative_hazard(model, base_times)
  cumulative_excess_hazard <- numeric(length(base_times))

  for (i in 2:length(base_times)) {
    vaccinated_td <- as.integer(base_times[i - 1] >= vaccination_day)

    profile_row <- data.frame(
      age_band = factor(age_band, levels = age_band_labels),
      sex = factor(sex, levels = c("female", "male")),
      vaccinated_td = vaccinated_td,
      school_member = factor(school_member, levels = c("No", "Yes")),
      worker = factor(worker, levels = c("No", "Yes"))
    )

    linear_predictor <- covariate_linear_predictor(model, profile_row)
    baseline_increment <- base_lambda[i] - base_lambda[i - 1]

    cumulative_excess_hazard[i] <-
      cumulative_excess_hazard[i - 1] +
      baseline_increment * exp(linear_predictor)
  }

  expected_profile <- data.frame(
    followup_days = max(base_times),
    age_days = age_years * days_per_year,
    study_date = make_base_date(anchor_year),
    sex = factor(sex, levels = c("female", "male"))
  )

  expected_fit <- survexp(
    Surv(followup_days, 0) ~ 1,
    data = expected_profile,
    ratetable = poptab,
    method = "conditional",
    rmap = list(age = age_days, year = study_date, sex = sex)
  )

  expected_summary <- summary(
    expected_fit,
    times = base_times[-1],
    extend = TRUE
  )
  population_survival <- c(1, expected_summary$surv)
  cumulative_population_hazard <- -log(population_survival)
  overall_survival <- exp(
    -(cumulative_excess_hazard + cumulative_population_hazard)
  )
  cumulative_excess_mortality <- c(
    0,
    cumsum(overall_survival[-1] * diff(cumulative_excess_hazard))
  )

  tibble::tibble(
    profile = profile_label,
    time = base_times,
    cumulative_excess_hazard = cumulative_excess_hazard,
    cumulative_population_hazard = cumulative_population_hazard,
    net_survival = exp(-cumulative_excess_hazard),
    overall_survival = overall_survival,
    cumulative_excess_mortality = cumulative_excess_mortality,
    vaccination_pattern = ifelse(
      is.finite(vaccination_day),
      "Vaccinated at day 45",
      "Never vaccinated"
    )
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
