library(dplyr)
library(tidyr)
library(purrr)

# socialmixr provides POLYMOD contact matrices
if (!requireNamespace("socialmixr", quietly = TRUE)) install.packages("socialmixr")
library(socialmixr)

# ---- Run data wrangling pipeline -----------------------------------------
source("analysis/group-D/data_wrangling.R")
# Objects now available:
#   cohort_tv_daily  - person-day survival data (cohort of 10 000)
#   prevalence.data  - daily prevalence by (age_group, vax_status, day)
#   df_aged          - full population with age_group

# ---- Parameters (assumptions) --------------------------------------------
theta        <- 0.5   # relative reduction in infectiousness for vaccinated
study_end    <- 90L
age_breaks   <- c(0, 5, 18, 50, 65, Inf)
age_labels   <- c("0-4", "5-17", "18-49", "50-64", "65+")
K            <- length(age_labels)

# ---- POLYMOD contact matrix for Netherlands --------------------------------
# Aggregate POLYMOD (5-year bands) to the 5 age groups used in this study
polymod_raw <- socialmixr::contact_matrix(
  socialmixr::polymod,
  countries    = "Netherlands",
  age.limits   = c(0, 5, 18, 50, 65),
  symmetric    = TRUE,
  estimated.contact.age = "mean",
  missing.contact.age   = "remove"
)
# contact_matrix() returns a list; $matrix rows = participant, cols = contact
C_mat <- polymod_raw$matrix          # K × K, rows = from, cols = to
rownames(C_mat) <- age_labels
colnames(C_mat) <- age_labels

# ---- Population size per age group (from full simulated population) -------
N_age <- df_aged |>
  count(age_group, name = "N_a") |>
  arrange(age_group)

N_vec <- setNames(N_age$N_a, N_age$age_group)   # named vector, length K

# ---- Compute Ī(t_c, a) for each (day × age_group) -----------------------
# Ī(t_c, a) = Σ_{a'} [I_u^{a'}(t_c-1) + (1-theta)*I_v^{a'}(t_c-1)] * C_{a',a} / N_a
#
# prevalence.data has columns: age_group, vax_status, day, prevalence.daily

prev_wide <- prevalence.data |>
  ungroup() |>
  pivot_wider(names_from = vax_status, values_from = prevalence.daily,
              names_prefix = "I_", values_fill = 0) |>
  rename(I_u = I_0, I_v = I_1) |>
  mutate(I_eff = I_u + (1 - theta) * I_v)   # effective infectious per (a', day)

# For each calendar day and recipient age group, sum weighted contacts
foi_long <- map_dfr(seq_len(study_end), function(tc) {
  # Use prevalence from the previous day (lag = 1); day 1 gets zeros
  prev_tc <- prev_wide |> filter(day == tc - 1)

  map_dfr(age_labels, function(a) {
    Ibar_val <- if (nrow(prev_tc) == 0) {
      0
    } else {
      sum(map_dbl(age_labels, function(ap) {
        I_ap <- prev_tc |> filter(age_group == ap) |> pull(I_eff)
        I_ap <- if (length(I_ap) == 0) 0 else I_ap
        C_mat[ap, a] / N_vec[a] * I_ap
      }))
    }
    tibble(day = tc, age_group = a, Ibar = pmax(Ibar_val, 1e-10))
  })
})

# ---- Aggregate person-days for Stan --------------------------------------
# Daily bins for time since vaccination.
# Bin 1 = unvaccinated; bins 2..(study_end+1) = day 1..study_end since vaccination
T_vax <- study_end + 1L   # bin 1 reserved for unvaccinated

stan_input_raw <- cohort_tv_daily |>
  mutate(
    time_since_vax = if_else(
      vax_status == 1 & !is.na(time_of_vaccination),
      as.integer(floor(tstop - time_of_vaccination)),
      NA_integer_
    ),
    # bin 1 = unvaccinated; bin t+1 = day t since vaccination (capped at study_end)
    tvax_day = if_else(
      vax_status == 0,
      1L,
      pmin(time_since_vax, study_end) + 1L
    )
  ) |>
  left_join(foi_long, by = c("tstop" = "day", "age_group" = "age_group")) %>%
  mutate(day = tstop)

# Aggregate to cells
stan_cells <- stan_input_raw |>
  group_by(age_group, vax_status, day, tvax_day) |>
  summarise(
    Y          = sum(event),
    n_persdays = n(),
    Ibar       = first(Ibar),
    .groups    = "drop"
  ) |>
  filter(!is.na(Ibar))

N_cells <- nrow(stan_cells)
cat("Stan data cells:", N_cells, "\n")
cat("Total events:   ", sum(stan_cells$Y), "\n")

stan_data <- list(
  N           = N_cells,
  K           = K,
  T_vax       = T_vax,
  Y           = as.integer(stan_cells$Y),
  log_Ibar    = log(stan_cells$Ibar),
  log_persdays = log(stan_cells$n_persdays),
  vax         = as.integer(stan_cells$vax_status),
  tvax_bin    = as.integer(stan_cells$tvax_day)
)

# ---- No-age-structure FoI ------------------------------------------------
# Simple alternative: overall effective prevalence / total population,
# the same for every age group on a given day.
N_total <- nrow(df_aged)

foi_simple <- prev_wide |>
  group_by(day) |>
  summarise(I_eff_total = sum(I_eff), .groups = "drop") |>
  mutate(
    day_lag  = day + 1L,   # attach lagged prevalence to the following calendar day
    Ibar_simple = pmax(I_eff_total / N_total, 1e-10)
  ) |>
  select(day = day_lag, Ibar_simple)

stan_cells_simple <- stan_input_raw |>
  left_join(foi_simple, by = c("day" = "day")) |>
  group_by(age_group, vax_status, day, tvax_day) |>
  summarise(
    Y            = sum(event),
    n_persdays   = n(),
    Ibar_simple  = first(Ibar_simple),
    .groups      = "drop"
  ) |>
  filter(!is.na(Ibar_simple))

stan_data_simple <- list(
  N            = nrow(stan_cells_simple),
  K            = K,
  T_vax        = T_vax,
  Y            = as.integer(stan_cells_simple$Y),
  log_Ibar     = log(stan_cells_simple$Ibar_simple),
  log_persdays = log(stan_cells_simple$n_persdays),
  vax          = as.integer(stan_cells_simple$vax_status),
  tvax_bin     = as.integer(stan_cells_simple$tvax_day)
)

# ---- Save for use in run_model.R -----------------------------------------
saveRDS(stan_data,          "analysis/group-D/stan_data.rds")
saveRDS(stan_cells,         "analysis/group-D/stan_cells.rds")
saveRDS(stan_data_simple,   "analysis/group-D/stan_data_simple.rds")
saveRDS(stan_cells_simple,  "analysis/group-D/stan_cells_simple.rds")
saveRDS(C_mat,              "analysis/group-D/contact_matrix_NL.rds")

cat("Saved stan_data.rds (age-structured) and stan_data_simple.rds (no age structure)\n")
cat("T_vax bins:", T_vax, "(bin 1 = unvaccinated, bins 2-", T_vax, "= days 1-", study_end, "since vaccination)\n")
