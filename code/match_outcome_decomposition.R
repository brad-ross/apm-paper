devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")
source("code/estimation_helpers.R")

library(tidyverse)
library(arrow)

print(str_glue("Running match outcome decomposition for FIRST_YEAR={FIRST_YEAR}, LAST_YEAR={LAST_YEAR}, MAX_RANK={MAX_RANK}, YEAR_CLUSTER_SIZE={YEAR_CLUSTER_SIZE}, MIN_COHORT_SIZE={MIN_COHORT_SIZE}, CHOSEN_K={CHOSEN_K}"))

# Load data using shared helpers
clustered_panels <- load_clustered_panels()
full_panel_df <- filter_panel_by_k(clustered_panels)
panel <- create_estimation_panel(full_panel_df)

# Load firm clusters and compute outcome weights
firm_clusters <- load_firm_clusters()
num_firms_by_outcome <- compute_num_firms_by_outcome(firm_clusters, full_panel_df, panel)

# Get estimation specifications
est_specs <- get_default_est_specs(include_extended = FALSE)

# Get bootstrap draws
wb <- get_default_bootstrap(panel)

# Compute outcome indices by province for the target year
outcome_indices_by_province <- {
    outcome_ids <- panel$get_outcome_ids()
    components <- strsplit(outcome_ids, ":", fixed = TRUE)
    years <- vapply(components, `[`, character(1), 1)
    provinces <- vapply(components, `[`, character(1), 2)
    target_year <- as.character(LAST_YEAR - YEAR_CLUSTER_SIZE + 1L)
    matching_indices <- which(years == target_year)
    split(matching_indices, provinces[matching_indices])
}

print("Estimating target parameter components...")

timing_target_param_components <- system.time({
    target_param_components <- tryCatch(
        est_target_param_components(panel,
            est_specs = est_specs,
            bootstrap = wb,
            num_threads = num_threads,
            est_outcome_means_via_imputation = TRUE,
            imputation_options = list(solver = "lsmr", lsmr_atol = 1e-8, lsmr_btol = 1e-8, lsmr_max_iters = 1000, lsmr_conlim = Inf)
        ),
        error = function(err) {
            warning(str_glue("est_target_param_components failed: {err$message}"))
            NULL
        }
    )
})

print(str_glue("Estimated target parameters in {timing_target_param_components['elapsed']} seconds"))

print("Estimating match outcome decomposition...")

match_outcome_decomposition <- est_avg_fgw_bipartite_match_outcome_diff_params(
    target_param_components$outcome_means,
    panel$get_observed_outcome_indices(),
    suff_stats = target_param_components$cohort_outcome_mean_ests,
    outcome_groupings = outcome_indices_by_province,
    outcome_weights = num_firms_by_outcome,
    num_threads = num_threads
)

target_param_diff <- get_target_param_diff_ests(
    match_outcome_decomposition$pc_1_no_fes_equal_weights,
    match_outcome_decomposition$twfe
)

combined_target_param_ests <- combine_target_param_ests(list(
    match_outcome_decomposition$pc_1_no_fes_equal_weights$subset(1),
    match_outcome_decomposition$twfe$subset(1),
    target_param_diff$subset(1)
))

print("Performing bootstrap inference...")

bootstrap_inference <- target_param_inference(combined_target_param_ests, panel, sig_level = 0.05)

target_param_table_data <- bootstrap_inference$as_data_frame() |>
    mutate(spec = c("APM", "TWFE", "Difference"), parameter = "Province") |>
    select(spec, parameter, everything())

write_parquet(target_param_table_data, file.path(RESULTS_PATH, "match_outcome_decomp_results.parquet"))

print("Match outcome decomposition results saved to parquet.")