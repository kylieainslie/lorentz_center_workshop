# =====================================================================
# ANALYSIS FOR EXCESS MORTALITY
# Author: Miguel Beynaerts
# Email: miguel.beynaerts@uhasselt.be
# Last modified: 17/06/2026 14:12:47
# License: MIT
#
# Description:
# R script for the analysis of excess mortality using simulated IBM data
# during the workshop "Connecting Survival Analysis and Infectious
# Disease Modeling", organized at the Lorenz Centre in Leiden
# (15-19 June 2026)
# ======================================================================

here::i_am("excess_mort.R")

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(relsurv)

set_theme(theme_minimal(base_size = 12))

source(here("excess_mort_functions.R"))

days_per_year <- 365.241

# Years at which population mortality rates are used from the lifetime tables
# The primary year for analysis is 2020
anchor_years <- c(2019L, 2020L, 2023L)
primary_year <- 2020L

# Some work to make analysis easier
age_band_breaks <- c(0, 5, 18, 40, 65, 86)
age_band_labels <- c("0-4", "5-17", "18-39", "40-64", "65-85")
baseline_breaks_days <- c(0, 30, 60, 90)
baseline_breaks_years <- baseline_breaks_days / days_per_year
profile_times <- 1:90

figures_dir <- file.path(here("figures"))
results_dir <- file.path(here("results"))

# dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. PREPARE DATA ----
# =====================================================================

mort_data <- readRDS(here("data", "mortality_data.rds"))
demo_data <- readRDS(here("data", "demographic_data.rds"))
vax_data <- readRDS(here("data", "vaccination_data.rds"))

# Make population lifetime table readable for the relsurv package
poptab <- transrate.hmd(
  male = here("data", "mltper_1x1.txt"),
  female = here("data", "fltper_1x1.txt")
)

# Combining mortality, demographic and vaccination data
cohort_data <- demo_data |>
  left_join(mort_data, by = "ind_id") |>
  left_join(vax_data, by = "ind_id") |>
  transmute(
    ind_id,
    age,
    sex = factor(
      dplyr::recode(sex, Male = "male", Female = "female"),
      levels = c("female", "male"), # This is needed for relsurv!!
    ),
    hh_id,
    classroom_id,
    workplace_id,
    time_of_death,
    time_of_vaccination,
    followup_days = if_else(is.na(time_of_death), 90, pmin(time_of_death, 90)),
    event = if_else(is.na(time_of_death) | time_of_death > 90, 0L, 1L), # Censoring indicator for death before t=90
    vaccinated_ever = if_else(
      # Indicator for vaccination status
      !is.na(time_of_vaccination) & time_of_vaccination <= 90,
      1L,
      0L
    ),
    age_days = age * days_per_year, # Age should be in days (instead of years) for relsurv
    age_band = cut(
      age,
      breaks = age_band_breaks,
      labels = age_band_labels,
      right = FALSE,
      include.lowest = TRUE
    ),
    school_member = factor(
      # Indicator for schoolplace membership (Yes/No)
      if_else(is.na(classroom_id), "No", "Yes"),
      levels = c("No", "Yes")
    ),
    worker = factor(
      # Indicator for workplace membership (Yes/No)
      if_else(is.na(workplace_id), "No", "Yes"),
      levels = c("No", "Yes")
    )
  )

# stopifnot(nrow(cohort_data) == 200000)
# stopifnot(dplyr::n_distinct(cohort_data$ind_id) == 200000)
# stopifnot(sum(cohort_data$event) == 29447)
# stopifnot(sum(cohort_data$vaccinated_ever) == 55221)
# stopifnot(all(cohort_data$followup_days >= 0))

# See helper file for detailed explanation
# Basically, this function creates start and stop times based on vaccination status for each individual
# In this way, we can account for the fact that vaccination is not a baseline covariate,
# but rather a time-dependent covariate
primary_split <- make_split_dataset(cohort_data, primary_year)

# Some sanity checks on the created datasets, uncomment to do these

# stopifnot(all(primary_split$interval_days > 0))
# stopifnot(all(primary_split$start < primary_split$stop))
# stopifnot(all(
#   primary_split |>
#     group_by(ind_id) |>
#     summarise(events = sum(event), .groups = "drop") |>
#     pull(events) <=
#     1
# ))

# data_checks <- tibble::tribble(
#   ~check                             ,
#   ~value                             , ~passed                                                                                                    ,
#   "Unique individuals"               , dplyr::n_distinct(cohort_data$ind_id)                                                                      , dplyr::n_distinct(cohort_data$ind_id) == 200000                                                                 ,
#   "Deaths preserved"                 , sum(cohort_data$event)                                                                                     , sum(cohort_data$event) == 29447                                                                                 ,
#   "Vaccinations preserved"           , sum(cohort_data$vaccinated_ever)                                                                           , sum(cohort_data$vaccinated_ever) == 55221                                                                       ,
#   "Non-negative follow-up"           , min(cohort_data$followup_days)                                                                             , min(cohort_data$followup_days) >= 0                                                                             ,
#   "Positive split intervals"         , min(primary_split$interval_days)                                                                           , min(primary_split$interval_days) > 0                                                                            ,
#   "At most one event per individual" , max(primary_split |> group_by(ind_id) |> summarise(events = sum(event), .groups = "drop") |> pull(events)) , max(primary_split |> group_by(ind_id) |> summarise(events = sum(event), .groups = "drop") |> pull(events)) <= 1
# )

# write_csv(data_checks, "data_checks.csv")

# =====================================================================
# 2. EXPLORATORY ANALYSIS ----
# =====================================================================

# Population-baseline risk (expected risk) for each individual
cohort_expected_primary <- add_expected_risk(cohort_data, poptab, primary_year)
primary_split_expected <- add_expected_interval_risk(primary_split, poptab) # Expected survival during observed followup

# This is just observed overall survival vs expected overall survival for each age in the dataset
overall_curve <- overall_curve_data(cohort_data, poptab, primary_year)
write_csv(overall_curve, "overall_observed_vs_expected_curve.csv")

# we calculate the excess cumulative mortality for each time point 0 to 90 as
# observed cumul. mortality - expected cumul. mortality
overall_excess_curve <- overall_curve |>
  select(time, measure, cumulative_mortality) |>
  tidyr::pivot_wider(
    names_from = measure,
    values_from = cumulative_mortality
  ) |>
  mutate(excess_cumulative_mortality = Observed - Expected)

write_csv(overall_excess_curve, "overall_excess_curve.csv")

# The same as above, but stratified over risk factors
age_sex_summary <- summarise_obs_expected(
  cohort_expected_primary,
  c("age_band", "sex")
)
age_band_summary <- summarise_obs_expected(
  cohort_expected_primary,
  c("age_band")
)

risk_factor_summary <- bind_rows(
  summarise_obs_expected(cohort_expected_primary, c("age_band", "sex")) |>
    rename(level = sex) |>
    mutate(risk_factor = "Sex"),
  summarise_obs_expected(
    cohort_expected_primary,
    c("age_band", "school_member")
  ) |>
    rename(level = school_member) |>
    mutate(risk_factor = "School membership"),
  summarise_obs_expected(cohort_expected_primary, c("age_band", "worker")) |>
    rename(level = worker) |>
    mutate(risk_factor = "Workplace membership")
)

# Same as above, but now stratified on vaccination status
vaccination_eda <- primary_split_expected |>
  mutate(
    vaccination_state = if_else(
      vaccinated_td == 1L,
      "After vaccination",
      "Before vaccination"
    )
  ) |>
  group_by(age_band, vaccination_state) |>
  summarise(
    person_days = sum(interval_days),
    observed_deaths = sum(event),
    expected_deaths = sum(expected_interval_death_prob),
    observed_rate_per_1000_pd = observed_deaths / person_days * 1000,
    expected_rate_per_1000_pd = expected_deaths / person_days * 1000,
    .groups = "drop"
  )

write_csv(age_sex_summary, "age_sex_90day_summary.csv")
write_csv(age_band_summary, "age_band_90day_summary.csv")
write_csv(risk_factor_summary, "risk_factor_90day_summary.csv")
write_csv(vaccination_eda, "vaccination_eda_summary.csv")

overall_curve_plot <- ggplot(
  overall_curve,
  aes(x = time, y = cumulative_mortality, color = measure)
) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c("Observed" = "#b22222", "Expected" = "#1f77b4")
  ) +
  labs(
    title = "Observed vs expected cumulative mortality over 90 days",
    x = "Days of follow-up",
    y = "Cumulative mortality",
    color = NULL
  )

overall_excess_plot <- ggplot(
  overall_excess_curve,
  aes(x = time, y = excess_cumulative_mortality)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.1, color = "#7f0000") +
  labs(
    title = "Excess cumulative mortality over 90 days",
    subtitle = "Observed minus expected cumulative mortality",
    x = "Days of follow-up",
    y = "Excess cumulative mortality"
  )

age_sex_plot_data <- age_sex_summary |>
  select(age_band, sex, observed_mortality, expected_mortality) |>
  pivot_longer(
    cols = c(observed_mortality, expected_mortality),
    names_to = "measure",
    values_to = "mortality"
  ) |>
  mutate(
    measure = recode(
      measure,
      observed_mortality = "Observed",
      expected_mortality = "Expected"
    )
  )

age_sex_plot <- ggplot(
  age_sex_plot_data,
  aes(
    x = age_band,
    y = mortality,
    color = sex,
    group = interaction(sex, measure)
  )
) +
  geom_line(aes(linetype = measure), linewidth = 1) +
  geom_point(aes(shape = measure), size = 2.2) +
  scale_color_manual(values = c("female" = "#d95f02", "male" = "#1b9e77")) +
  labs(
    title = "Observed and expected 90-day mortality by age band and sex",
    x = "Age band",
    y = "90-day mortality",
    color = "Sex",
    linetype = NULL,
    shape = NULL
  )

age_band_excess_plot <- ggplot(
  age_band_summary,
  aes(x = age_band, y = excess_per_1000, fill = age_band)
) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("O/E %.1f", oe_ratio)),
    vjust = -0.4,
    size = 3.3
  ) +
  labs(
    title = "Absolute excess mortality by age band",
    x = "Age band",
    y = "Excess deaths per 1,000 people"
  )

risk_factor_plot <- ggplot(
  risk_factor_summary,
  aes(x = age_band, y = excess_per_1000, fill = level)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~risk_factor, scales = "free_y") +
  labs(
    title = "Excess mortality for key risk factors within age bands",
    x = "Age band",
    y = "Excess deaths per 1,000 people",
    fill = NULL
  )

vaccination_plot_data <- vaccination_eda |>
  select(
    age_band,
    vaccination_state,
    observed_rate_per_1000_pd,
    expected_rate_per_1000_pd
  ) |>
  pivot_longer(
    cols = c(observed_rate_per_1000_pd, expected_rate_per_1000_pd),
    names_to = "measure",
    values_to = "rate"
  ) |>
  mutate(
    measure = recode(
      measure,
      observed_rate_per_1000_pd = "Observed rate",
      expected_rate_per_1000_pd = "Expected rate"
    )
  )

vaccination_plot <- ggplot(
  vaccination_plot_data,
  aes(
    x = age_band,
    y = rate,
    color = vaccination_state,
    group = vaccination_state
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~measure, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Before vaccination" = "#8c510a",
      "After vaccination" = "#01665e"
    )
  ) +
  labs(
    title = "Observed and expected death rates by vaccination state and age band",
    x = "Age band",
    y = "Deaths per 1,000 person-days",
    color = NULL
  )

save_plot(
  "overall_observed_vs_expected.pdf",
  overall_curve_plot,
  width = 9,
  height = 5.5
)
save_plot(
  "overall_excess_mortality.pdf",
  overall_excess_plot,
  width = 9,
  height = 5.5
)
save_plot("age_sex_mortality.pdf", age_sex_plot, width = 9, height = 5.5)
save_plot(
  "age_band_excess.pdf",
  age_band_excess_plot,
  width = 8.5,
  height = 5.5
)
save_plot("risk_factor_excess.pdf", risk_factor_plot, width = 10, height = 6.5)
save_plot(
  "vaccination_observed_expected_rates.pdf",
  vaccination_plot,
  width = 10,
  height = 6
)

# =====================================================================
# 3. FORMAL RELATIVE SURVIVAL ANALYSIS ----
# =====================================================================

# Now we estimate the net survival (estimand of interest) using several models:
# 1. Firstly, we make a distinction between the year in which the pandemic takes place
# Since the code uses baseline population rates from only 1 year, we have to make a choice
# as to when the pandemic takes place. We do this for several years as a sensitivity analysis
# 2. We fit several models wrt the covariates included in the model (see fit_rs_models ):
# no_covariates = overall net survival
# minimal = adjustment for age class, sex and vaccination status
# extended = minimal + workplace membership + household membership

# Again, we create start-stop times for each individual based on vaccination status
# but now we do this for multiple study years
split_data_by_year <- lapply(anchor_years, function(year) {
  make_split_dataset(cohort_data, year)
})
names(split_data_by_year) <- as.character(anchor_years)

# Fit additive hazard models
models_by_year <- lapply(
  split_data_by_year,
  function(split_data) fit_rs_models(split_data, poptab)
)

# Extract hazard ratios from the models
model_coefficients <- bind_rows(lapply(
  names(models_by_year),
  function(year_name) {
    purrrless <- names(models_by_year[[year_name]])
    bind_rows(lapply(purrrless, function(model_name) {
      extract_model_coefficients(
        model = models_by_year[[year_name]][[model_name]],
        model_name = model_name,
        anchor_year = as.integer(year_name)
      )
    }))
  }
))

# This is basically a csv reference of the models fitted to the data
model_overview <- bind_rows(lapply(names(models_by_year), function(year_name) {
  purrrless <- names(models_by_year[[year_name]])
  bind_rows(lapply(purrrless, function(model_name) {
    extract_model_overview(
      model = models_by_year[[year_name]][[model_name]],
      model_name = model_name,
      anchor_year = as.integer(year_name),
      split_data = split_data_by_year[[year_name]]
    )
  }))
}))

write_csv(model_coefficients, "model_coefficients.csv")
write_csv(model_overview, "model_overview.csv")
saveRDS(
  models_by_year[[as.character(primary_year)]],
  file.path(results_dir, "primary_models_2020.rds")
)

# These are overall net survival curves stratified over all risk factors separately (much like KM-estimator)
rs_curves <- fit_rs_curves(cohort_data, poptab, primary_year)

# This is just data cleanup for readability
net_survival_overall <- tidy_surv_object(rs_curves$overall, label = "Overall")
net_survival_age_band <- tidy_surv_object(rs_curves$age_band) |>
  mutate(strata = clean_strata_label(strata))
net_survival_sex <- tidy_surv_object(rs_curves$sex) |>
  mutate(strata = clean_strata_label(strata), variable = "Sex")
net_survival_school <- tidy_surv_object(rs_curves$school_member) |>
  mutate(strata = clean_strata_label(strata), variable = "School membership")
net_survival_worker <- tidy_surv_object(rs_curves$worker) |>
  mutate(strata = clean_strata_label(strata), variable = "Workplace membership")

write_csv(net_survival_overall, "net_survival_overall.csv")
write_csv(net_survival_age_band, "net_survival_age_band.csv")
write_csv(net_survival_sex, "net_survival_sex.csv")
write_csv(net_survival_school, "net_survival_school.csv")
write_csv(net_survival_worker, "net_survival_worker.csv")

net_survival_age_plot <- ggplot(
  net_survival_age_band,
  aes(x = time, y = survival, color = strata)
) +
  geom_line(linewidth = 1) +
  labs(
    title = "Net survival by age band",
    x = "Days of follow-up",
    y = "Net survival",
    color = "Age band"
  )

net_survival_risk_plot <- bind_rows(
  net_survival_sex,
  net_survival_school,
  net_survival_worker
) |>
  ggplot(aes(x = time, y = survival, color = strata)) +
  geom_line(linewidth = 1) +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    title = "Net survival for fixed risk-factor strata",
    x = "Days of follow-up",
    y = "Net survival",
    color = NULL
  )

save_plot(
  "net_survival_age_band.pdf",
  net_survival_age_plot,
  width = 9,
  height = 5.5
)
save_plot(
  "net_survival_risk_factors.pdf",
  net_survival_risk_plot,
  width = 10,
  height = 6.5
)

primary_extended_model <- models_by_year[[as.character(primary_year)]][[
  "extended"
]]


# Here we predict the net survival for some specific reference individuals in the population.
profile_predictions <- bind_rows(
  predict_profile_curve(
    model = primary_extended_model,
    poptab = poptab,
    anchor_year = primary_year,
    age_years = 12,
    age_band = "5-17",
    sex = "female",
    school_member = "Yes",
    worker = "No",
    vaccination_day = Inf,
    profile_label = "School-aged child, never vaccinated"
  ),
  predict_profile_curve(
    model = primary_extended_model,
    poptab = poptab,
    anchor_year = primary_year,
    age_years = 50,
    age_band = "40-64",
    sex = "female",
    school_member = "No",
    worker = "Yes",
    vaccination_day = Inf,
    profile_label = "Working-age adult, never vaccinated"
  ),
  predict_profile_curve(
    model = primary_extended_model,
    poptab = poptab,
    anchor_year = primary_year,
    age_years = 75,
    age_band = "65-85",
    sex = "female",
    school_member = "No",
    worker = "No",
    vaccination_day = Inf,
    profile_label = "Older adult, never vaccinated"
  ),
  predict_profile_curve(
    model = primary_extended_model,
    poptab = poptab,
    anchor_year = primary_year,
    age_years = 75,
    age_band = "65-85",
    sex = "female",
    school_member = "No",
    worker = "No",
    vaccination_day = 45,
    profile_label = "Older adult, vaccinated at day 45"
  )
)

profile_prediction_summary <- profile_predictions |>
  filter(time == 90) |>
  select(profile, cumulative_excess_mortality, net_survival, overall_survival)

write_csv(profile_predictions, "profile_predictions_2020.csv")
write_csv(profile_prediction_summary, "profile_prediction_summary_2020.csv")

forest_plot_data <- model_coefficients |>
  filter(
    anchor_year == primary_year,
    component == "covariate",
    model_name %in% c("minimal", "extended")
  ) |>
  mutate(
    model_name = recode(model_name, minimal = "Minimal", extended = "Extended"),
    term_label = factor(term_label, levels = rev(unique(term_label)))
  )

forest_plot <- ggplot(
  forest_plot_data,
  aes(x = hazard_ratio, y = term_label, color = model_name)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(
    aes(xmin = hazard_ratio_low, xmax = hazard_ratio_high),
    orientation = "y",
    width = 0.18,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2) +
  scale_x_log10() +
  labs(
    title = "Excess hazard ratios from the 2020 relative survival models",
    x = "Excess hazard ratio",
    y = NULL,
    color = NULL
  )

profile_plot <- ggplot(
  profile_predictions,
  aes(x = time, y = cumulative_excess_mortality, color = profile)
) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Model-based cumulative excess mortality for representative profiles",
    x = "Days of follow-up",
    y = "Cumulative excess mortality",
    color = NULL
  )

save_plot("model_forest_plot.pdf", forest_plot, width = 10, height = 6.5)
save_plot("model_profile_predictions.pdf", profile_plot, width = 10, height = 6)

# =====================================================================
# 4. SENSITIVITY ANALYSIS ----
# =====================================================================

sensitivity_table <- model_coefficients |>
  filter(
    component == "covariate",
    model_name %in% c("minimal", "extended"),
    term %in% c("vaccinated_td", "age_band65.85", "sexmale")
  ) |>
  mutate(
    model_name = recode(model_name, minimal = "Minimal", extended = "Extended")
  ) |>
  select(
    anchor_year,
    model_name,
    term_label,
    hazard_ratio,
    hazard_ratio_low,
    hazard_ratio_high
  )

write_csv(sensitivity_table, "year_sensitivity_summary.csv")

sensitivity_plot <- sensitivity_table |>
  filter(term_label == "Vaccinated (time-dependent)") |>
  ggplot(aes(
    x = factor(anchor_year),
    y = hazard_ratio,
    color = model_name,
    group = model_name
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(
    aes(ymin = hazard_ratio_low, ymax = hazard_ratio_high),
    width = 0.12,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(size = 2.4, position = position_dodge(width = 0.4)) +
  geom_line(position = position_dodge(width = 0.4)) +
  scale_y_log10() +
  labs(
    title = "Vaccination effect across ratetable anchor years",
    x = "Anchor year",
    y = "Excess hazard ratio",
    color = NULL
  )

save_plot(
  "sensitivity_vaccination_effect.pdf",
  sensitivity_plot,
  width = 8.5,
  height = 5.5
)

message("Excess mortality analysis outputs written to ", here())
