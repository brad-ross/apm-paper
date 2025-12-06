# Shared estimation helper functions for replication scripts
#
# This file provides common functions for loading data and setting up estimations.
# Requires: tidyverse, arrow, and the apm package to be loaded.
# Requires: env_config.R and get_io_paths.R to be sourced first.

library(arrow)

# Load clustered panels with log earnings transformation
load_clustered_panels <- function(first_year = FIRST_YEAR, 
                                   last_year = LAST_YEAR,
                                   year_cluster_size = YEAR_CLUSTER_SIZE,
                                   max_rank = 2L) {
    read_parquet(file.path(PROCESSED_DATA_PATH, 
        str_glue("panels_with_clustered_outcomes_start_year={first_year}_end_year={last_year}_year_cluster_size={year_cluster_size}_max_rank={max_rank}.parquet"))) |>
        .to_arrow_dataset() |>
        mutate(log_avg_weekly_earnings = log10(avg_weekly_earnings)) |>
        collect()
}

# Filter clustered panels to a specific k value
filter_panel_by_k <- function(clustered_panels, chosen_k = CHOSEN_K) {
    clustered_panels |>
        .to_arrow_dataset() |>
        filter(k == chosen_k) |>
        collect()
}

# Create an UnbalancedPanel for estimation
create_estimation_panel <- function(panel_df,
                                     unit_col = "worker",
                                     outcome_col = "outcome", 
                                     value_col = "log_avg_weekly_earnings",
                                     model_rank = MAX_RANK,
                                     min_cohort_size = MIN_COHORT_SIZE,
                                     subset_to_largest_super_cohort = TRUE) {
    UnbalancedPanel$new(panel_df,
        unit_col, outcome_col, value_col,
        model_rank = model_rank, 
        min_cohort_size = min_cohort_size,
        subset_to_largest_super_cohort = subset_to_largest_super_cohort)
}

# Load firm clusters data
load_firm_clusters <- function(first_year = FIRST_YEAR,
                                last_year = LAST_YEAR,
                                year_cluster_size = YEAR_CLUSTER_SIZE) {
    read_parquet(file.path(PROCESSED_DATA_PATH, 
        str_glue("firm_clusters_start_year={first_year}_end_year={last_year}_year_cluster_size={year_cluster_size}.parquet")))
}

# Get default estimation specifications (TWFE and PC)
get_default_est_specs <- function(include_extended = FALSE) {
    specs <- list(
        "twfe" = list(
            factor_model_estimator = "twfe",
            include_outcome_fes = TRUE,
            r = 1
        ),
        "pc_1_no_fes_equal_weights" = list(
            factor_model_estimator = "principal_components",
            include_outcome_fes = FALSE,
            r = 1,
            cohort_weighting = "equal"
        )
    )
    
    if (include_extended) {
        specs <- c(specs, list(
            "pc_1_fes_equal_weights" = list(
                factor_model_estimator = "principal_components",
                include_outcome_fes = TRUE,
                r = 1,
                cohort_weighting = "equal"
            ),
            "pc_2_no_fes_equal_weights" = list(
                factor_model_estimator = "principal_components",
                include_outcome_fes = FALSE,
                r = 2,
                cohort_weighting = "equal"
            ),
            "pc_2_fes_equal_weights" = list(
                factor_model_estimator = "principal_components",
                include_outcome_fes = TRUE,
                r = 2,
                cohort_weighting = "equal"
            )
        ))
    }
    
    specs
}

# Get default bootstrap draws
get_default_bootstrap <- function(panel, n_draws = 1000, seed = 9001L) {
    get_weighted_bootstrap_draws(panel$get_num_units(), n_draws, type = "bayesian", seed = seed)
}

# Compute number of firms by outcome for weighting
compute_num_firms_by_outcome <- function(firm_clusters, panel_df, panel, chosen_k = CHOSEN_K) {
    firm_clusters |>
        filter(k == chosen_k) |>
        group_by(province, cluster) |>
        summarize(
            firms_for_outcome = n(),
            .groups = "drop"
        ) |>
        inner_join(panel_df |> select(province, cluster, outcome) |> distinct(), 
                   by = c("province", "cluster")) |>
        mutate(outcome_idx = panel$get_outcome_to_index()[outcome]) |>
        select(outcome_idx, firms_for_outcome) |>
        arrange(outcome_idx) |>
        pull(firms_for_outcome)
}