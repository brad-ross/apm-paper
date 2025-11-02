devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")

num_threads <- as.integer(Sys.getenv("APM_MAX_THREADS"))
if (is.na(num_threads)) {
    num_threads <- NULL
}

library(tidyverse)
library(arrow)

FIRST_YEAR <- as.integer(Sys.getenv("FIRST_YEAR"))
LAST_YEAR <- as.integer(Sys.getenv("LAST_YEAR"))
MAX_RANK <- as.integer(Sys.getenv("MAX_RANK"))
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- as.integer(Sys.getenv("MIN_COHORT_SIZE"))

clustered_panels <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.parquet"))) |>
    .to_arrow_dataset() |>
    mutate(log_avg_weekly_earnings = log10(avg_weekly_earnings)) |>
    collect()

CHOSEN_K <- as.integer(Sys.getenv("CHOSEN_K"))

print(str_glue("Running leave-cohort-outcome-out evals for FIRST_YEAR={FIRST_YEAR}, LAST_YEAR={LAST_YEAR}, MAX_RANK={MAX_RANK}, YEAR_CLUSTER_SIZE={YEAR_CLUSTER_SIZE}, MIN_COHORT_SIZE={MIN_COHORT_SIZE}, CHOSEN_K={CHOSEN_K}"))

full_panel_df <- clustered_panels |>
    .to_arrow_dataset() |>
    filter(k == CHOSEN_K) |>
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
    # "pc_1_no_fes" = list(
    #     factor_model_estimator = "principal_components",
    #     include_outcome_fes = FALSE,
    #     r = 1
    # ),
    "pc_1_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = TRUE,
        r = 1
    ),
    # "pc_2_no_fes" = list(
    #     factor_model_estimator = "principal_components",
    #     include_outcome_fes = FALSE,
    #     r = 2
    # ),
    "pc_2_fes" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = TRUE,
        r = 2
    )
)

wb <- get_weighted_bootstrap_draws(panel$get_num_units(), 500, type = "bayesian", seed = 9001L)

MIN_MASKED_COHORT_SIZE <- 1000

observed_outcome_indices <- panel$get_observed_outcome_indices()
cohort_sizes <- panel$get_cohort_sizes()
cohort_order <- order(cohort_sizes, decreasing = TRUE)

masked_metrics <- NULL

for (cohort_pos in seq_along(cohort_order)) {
    cohort_id <- cohort_order[[cohort_pos]]
    cohort_outcomes <- observed_outcome_indices[[cohort_id]]

    if (length(cohort_outcomes) == 0L) {
        next
    }

    cohort_size <- cohort_sizes[[cohort_id]]

    if (cohort_size < MIN_MASKED_COHORT_SIZE) {
        next
    }

    for (outcome_id in cohort_outcomes) {
        print(str_glue("Checking if factors are still identified when masking cohort {cohort_id} ({cohort_size} workers), outcome {outcome_id}"))

        mask <- setNames(list(as.integer(outcome_id)), as.character(cohort_id))
        masked_ooi <- get_masked_observed_outcome_indices(observed_outcome_indices, mask)

        if (!aligned_factors_identified(masked_ooi, MAX_RANK)) {
            next
        }

        print(str_glue("Estimating error metrics for masked cohort {cohort_id}, outcome {outcome_id} mean"))

        timing <- system.time({
            comps <- tryCatch(
                est_target_param_components(panel,
                    est_specs = est_specs,
                    bootstrap = wb,
                    num_threads = num_threads,
                    cohort_outcomes_to_mask = mask,
                    est_outcome_means_via_imputation = TRUE
                ),
                error = function(err) {
                    warning(str_glue("est_target_param_components failed for masked cohort {cohort_id}, outcome {outcome_id}: {err$message}"))
                    NULL
                }
            )
        })

        if (is.null(comps)) {
            next
        }

        metrics <- tryCatch(
            est_masked_outcome_mean_err_metrics(comps),
            error = function(err) {
                warning(str_glue("est_masked_outcome_mean_err_metrics failed for masked cohort {cohort_id}, outcome {outcome_id}: {err$message}"))
                NULL
            }
        )

        if (is.null(metrics) || nrow(metrics) == 0L) {
            next
        }

        metrics <- metrics |> mutate(
            total_time = timing["elapsed"],
            cohort_size = cohort_size
        )
        masked_metrics <- masked_metrics |> bind_rows(metrics)

        # Checkpoint after every (cohort, outcome) pair
        write_parquet(masked_metrics,
            file.path(RESULTS_PATH, str_glue("masked_metrics_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.parquet")))
    }
}