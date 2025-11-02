devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")

num_threads <- as.integer(Sys.getenv("APM_MAX_THREADS"))

library(tidyverse)
library(arrow)

FIRST_YEAR <- 1998
LAST_YEAR <- 2001
MAX_RANK <- 2
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- 50

clustered_panels <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.parquet"))) |>
    .to_arrow_dataset() |>
    mutate(log_avg_weekly_earnings = log10(avg_weekly_earnings)) |>
    collect()

chosen_k <- 14

full_panel_df <- clustered_panels |> 
    .to_arrow_dataset() |> 
    filter(k == chosen_k) |> 
    collect()

panel <- UnbalancedPanel$new(full_panel_df, 
    "worker", "outcome", "log_avg_weekly_earnings", 
    model_rank = MAX_RANK, min_cohort_size = MIN_COHORT_SIZE,
    subset_to_largest_super_cohort = TRUE)

est_specs <- list(
    "twfe" = list(
        factor_model_estimator = "twfe",
        include_outcome_fes = TRUE,
        r = 1
    ),
    "pc_1_no_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = FALSE,
        r = 1
    ),
    "pc_1_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = TRUE,
        r = 1
    ),
    "pc_2_no_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = FALSE,
        r = 2
    ),
    "pc_2_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = TRUE,
        r = 2
    )
)

wb <- get_weighted_bootstrap_draws(panel$get_num_units(), 100, type = "bayesian", seed = 123L)

timing <- system.time({
    cohort_specific_params <- est_cohort_specific_params(
        panel, est_specs, 
        num_threads = num_threads
    )
})
print(timing)

timing <- system.time({
    target_param_components <- est_target_param_components(panel, 
        est_specs = est_specs, 
        bootstrap = wb,
        num_threads = num_threads, 
        est_outcome_means_via_imputation = TRUE,
        imputation_options = list(
            tol = 1e-8
        )
    )
})
print(timing)