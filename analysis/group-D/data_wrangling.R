library(dplyr)
library(tidyr)
library(survival)
library(purrr)
require(ggplot2)
require(tidyverse)

# Load data
demo <- readRDS("data/demographic_data.rds")
vax  <- readRDS("data/vaccination_data.rds")
inc  <- readRDS("data/incidence_data_reduced.rds")

# Clean up snapshot time column (not needed for analysis)
demo <- select(demo, -time)
vax  <- select(vax,  -time)
inc  <- select(inc,  -time)

# Merge: start with full population, left join vaccination and infection times
df <- demo |>
  left_join(vax, by = "ind_id") |>
  left_join(inc, by = "ind_id") |>
  mutate(
    vaccinated = !is.na(time_of_vaccination),
    infected   = !is.na(time_of_reporting)
  )

glimpse(df)

# ----- Create synthetic cohort by random sampling -----
# Set n_sample to desired cohort size; set seed for reproducibility
# To switch to stratified sampling later, replace slice_sample() with
# a group_by() + slice_sample() call (e.g. group_by(vaccinated))

n_sample <- 10000
set.seed(42)

cohort <- df |>
  slice_sample(n = n_sample)

cat("\nCohort size:", nrow(cohort), "\n")
cat("Vaccinated:", sum(cohort$vaccinated), "\n")
cat("Infected:  ", sum(cohort$infected), "\n")

# ----- Weighted daily incidence by age group and vaccination status -----
# Computed on the full population (df), not just the cohort sample

study_end  <- 93L  # max observed reporting time
age_breaks <- c(0, 5, 18, 50, 65, Inf)
age_labels  <- c("0-4", "5-17", "18-49", "50-64", "65+")

df_aged <- df |>
  mutate(age_group = cut(age, breaks = age_breaks, labels = age_labels, right = FALSE))

# Numerator: new cases on each day by (age group, vaccination status at reporting)
cases_daily <- df_aged |>
  filter(infected) |>
  mutate(
    day        = as.integer(time_of_reporting),
    vax_status = as.integer(!is.na(time_of_vaccination) & time_of_vaccination < time_of_reporting)
  ) |>
  count(day, age_group, vax_status, name = "n_cases")

# Denominator: at-risk population on each day by (age group, vaccination status)
# A person is at risk on day t if they have not yet been infected (time_of_reporting >= t or NA)
# Their vaccination status on day t: vaccinated if time_of_vaccination < t
atrisk_daily <- map_dfr(1:study_end, \(t) {
  df_aged |>
    filter(is.na(time_of_reporting) | time_of_reporting >= t) |>
    mutate(vax_status = as.integer(!is.na(time_of_vaccination) & time_of_vaccination < t)) |>
    count(age_group, vax_status, name = "n_at_risk") |>
    mutate(day = t)
})

# Join and compute incidence rate and population weight
incidence_weighted <- atrisk_daily |>
  left_join(cases_daily, by = c("day", "age_group", "vax_status")) |>
  replace_na(list(n_cases = 0)) |>
  mutate(
    incidence_rate = n_cases / n_at_risk,
    weight         = n_at_risk
  ) |>
  arrange(day, age_group, vax_status)

glimpse(incidence_weighted)
cat("\nRows:", nrow(incidence_weighted), "\n")
cat("Unique (day × age group × vax status) combinations:", nrow(incidence_weighted), "\n")

# ----- Time-to-event with time-varying vaccination covariate -----
# Counting process format: (tstart, tstop, event)
# Vaccinated individuals are split into two intervals:
#   [0, time_of_vaccination)     vax_status = 0
#   [time_of_vaccination, tstop] vax_status = 1

study_end <- 90

cohort_base <- cohort |>
  mutate(
    tstop = if_else(infected, time_of_reporting, as.numeric(study_end)),
    event = as.integer(infected)
  )

# Step 1: initialise one interval per person (0, tstop)
cohort_tv <- tmerge(
  cohort_base |> select(ind_id, age, sex, tstop, event),
  cohort_base |> select(ind_id, tstop, event),
  id    = ind_id,
  event = event(tstop, event)
)

# Step 2: split at vaccination time — vax_status switches 0 → 1
cohort_tv <- tmerge(
  cohort_tv,
  cohort_base |> filter(vaccinated) |> select(ind_id, time_of_vaccination),
  id         = ind_id,
  vax_status = tdc(time_of_vaccination)
)

glimpse(cohort_tv)
cat("\nIntervals per person:\n")
cohort_tv |> count(ind_id) |> count(n, name = "n_individuals")

# ----- Join background incidence onto individual linelist -----
# incidence_weighted gives daily incidence_rate by (day, age_group, vax_status)
# cohort_tv has intervals (tstart, tstop); we expand to one row per person-day
# so the daily rate can be attached as a time-varying covariate

cohort_tv <- cohort_tv |>
  mutate(age_group = cut(age, breaks = age_breaks, labels = age_labels, right = FALSE))

# Expand each interval (tstart, tstop] to individual days, then join incidence.
# The event flag is 1 only on the last day of an interval that has event = 1.
cohort_tv_daily <- cohort_tv |>
  mutate(day = map2(as.integer(tstart) + 1L, as.integer(tstop), seq)) |>
  unnest(day) |>
  mutate(
    event  = if_else(day == as.integer(tstop) & event == 1L, 1L, 0L),
    tstart = day - 1L,
    tstop  = as.numeric(day)
  ) |>
  left_join(
    incidence_weighted |> select(day, age_group, vax_status, incidence_rate),
    by = c("day", "age_group", "vax_status")
  ) |>
  select(ind_id, age, sex, age_group, day, tstart, tstop, event, vax_status, incidence_rate) %>%
  left_join(., cohort_base %>% select(ind_id, time_of_vaccination), by = "ind_id")

glimpse(cohort_tv_daily)
cat("\nRows:", nrow(cohort_tv_daily), "\n")
cat("Missing incidence_rate:", sum(is.na(cohort_tv_daily$incidence_rate)), "\n")


# Calculating the prevalence (number of infections per day by vax and age) with initialy assumptions regarding
# the duration of infection and reporting delay -- to be updated once these characteristics are estimated
duration.infection <- 7
reporting.delay <- 2

prevalence.data <- cases_daily %>%
    rowwise() %>%
    mutate(inf.days = list(enframe(seq(day-reporting.delay, day + (duration.infection - reporting.delay))))) %>%
    select(-day) %>%
    unnest(inf.days) %>%
    group_by(age_group, vax_status, value) %>%
    summarise(prevalence.daily = sum(n_cases)) %>%
    rename(day = value)

prevalence.data.overall <- cases_daily %>%
  rowwise() %>%
  mutate(inf.days = list(enframe(seq(day-reporting.delay, day + (duration.infection - reporting.delay))))) %>%
  select(-day) %>%
  unnest(inf.days) %>%
  group_by(value) %>%
  summarise(prevalence.daily = sum(n_cases)) %>%
  rename(day = value)

cohort_tv_daily <- cohort_tv_daily |>
  left_join(
    prevalence.data.overall,
    by = c("tstart" = "day")
  ) |>
  rename(I = prevalence.daily) |>
  left_join(
    prevalence.data |> pivot_wider(names_from = age_group, values_from = prevalence.daily),
    by = c("vax_status", "tstart" = "day")
  )
names(cohort_tv_daily)[12:16] <- paste("Ia", 1:5, sep="")
