devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")
source("code/estimation_helpers.R")

library(tidyverse)
library(arrow)

print(str_glue("Running leave-cohort-outcome-out evals for FIRST_YEAR={FIRST_YEAR}, LAST_YEAR={LAST_YEAR}, MAX_RANK={MAX_RANK}, YEAR_CLUSTER_SIZE={YEAR_CLUSTER_SIZE}, MIN_COHORT_SIZE={MIN_COHORT_SIZE}, CHOSEN_K={CHOSEN_K}"))

# Load data using shared helpers
clustered_panels <- load_clustered_panels()
full_panel_df <- filter_panel_by_k(clustered_panels)
panel <- create_estimation_panel(full_panel_df)

# Get extended estimation specifications (includes additional PC specs)
est_specs <- get_default_est_specs(include_extended = TRUE)

# Compute cohort-specific parameters and check eigenvalue gaps
cohort_specific_params <- est_cohort_specific_params(panel, est_specs, num_threads=num_threads)
for (spec in names(est_specs)) {
    fmes <- lapply(cohort_specific_params$cohort_specific_factor_ests[[spec]], \(fme) fme$G())
    cohort_weights <- cohort_specific_params$cohort_weights[[spec]]$cohort_weights()[, 1]
    agg_fme <- compute_aggregated_projection_matrix(fmes, panel$get_observed_outcome_indices(), cohort_weights)
    r <- dim(fmes[[1]])[2]
    lambdas <- sort(eigen(agg_fme, symmetric = TRUE, only.values = TRUE)$values)
    lambda_r_plus_1 <- lambdas[r + 1L]
    lambda_r <- lambdas[r]
    print(str_glue("apm for spec {spec} (r = {r}) has:\n lambda_(1:r+2) = [{paste(lambdas[1:(r+2L)], collapse = ', ')}],\n lambda_(r+1) - lambda_r = {lambda_r_plus_1 - lambda_r},\n lambda_(r+1) / lambda_r = {lambda_r_plus_1/lambda_r}"))
}

# Get bootstrap draws
wb <- get_default_bootstrap(panel)

MIN_MASKED_COHORT_SIZE <- 1000

observed_outcome_indices <- panel$get_observed_outcome_indices()
cohort_sizes <- panel$get_cohort_sizes()
cohort_order <- order(cohort_sizes, decreasing = TRUE)

masked_metrics_path <- file.path(RESULTS_PATH, str_glue("masked_metrics_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.parquet"))
masked_metrics <- data.frame(cohort = integer(), outcome = integer(), spec = character())
if (file.exists(masked_metrics_path)) {
    masked_metrics <- read_parquet(masked_metrics_path)
}

for (cohort_pos in seq_along(cohort_order)) {
    cohort_id <- cohort_order[[cohort_pos]]
    cohort_outcomes <- observed_outcome_indices[[cohort_id]]

    if (length(cohort_outcomes) == 0L) {
        next
    }

    cohort_size <- cohort_sizes[[cohort_id]]

    if (cohort_size < MIN_MASKED_COHORT_SIZE) {
        break
    }

    for (outcome_id in cohort_outcomes) {
        print(str_glue("Checking if factors are still identified when masking cohort {cohort_id} ({cohort_size} workers), outcome {outcome_id}"))

        mask <- setNames(list(as.integer(outcome_id)), as.character(cohort_id))
        masked_ooi <- get_masked_observed_outcome_indices(observed_outcome_indices, mask)

        existing_est_specs_for_cohort_outcome <- masked_metrics |>
            filter(cohort == cohort_id, outcome == outcome_id) |>
            pull(spec)
        non_existing_est_specs_for_cohort_outcome <- setdiff(names(est_specs), existing_est_specs_for_cohort_outcome)

        if (!aligned_factors_identified(masked_ooi, MAX_RANK) || length(non_existing_est_specs_for_cohort_outcome) == 0) {
            next
        }

        timing <- system.time({
            comps <- tryCatch(
                est_target_param_components(panel,
                    est_specs = est_specs[setdiff(names(est_specs), existing_est_specs_for_cohort_outcome)],
                    bootstrap = wb,
                    num_threads = num_threads,
                    cohort_outcomes_to_mask = mask,
                    est_outcome_means_via_imputation = TRUE,
                    imputation_options = list(solver = "lsmr", lsmr_atol = 1e-8, lsmr_btol = 1e-8, lsmr_max_iters = 1000, lsmr_conlim = Inf)
                ),
                error = function(err) {
                    warning(str_glue("est_target_param_components failed for masked cohort {cohort_id}, outcome {outcome_id}: {err$message}"))
                    NULL
                }
            )
        })

        print(str_glue("Estimated masked cohort {cohort_id}, outcome {outcome_id} mean using est_specs {paste(non_existing_est_specs_for_cohort_outcome, collapse = ', ')} in {timing['elapsed']} seconds"))

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
        write_parquet(masked_metrics, masked_metrics_path)
    }
}

print(str_glue("Leave-cohort-outcome-out evaluation results saved to {masked_metrics_path}"))