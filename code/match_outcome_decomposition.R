devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")

library(tidyverse)
library(arrow)
library(kableExtra)

num_threads <- as.integer(Sys.getenv("APM_MAX_THREADS"))
if (is.na(num_threads)) {
    num_threads <- NULL
}

FIRST_YEAR <- as.integer(Sys.getenv("FIRST_YEAR"))
LAST_YEAR <- as.integer(Sys.getenv("LAST_YEAR"))
MAX_RANK <- as.integer(Sys.getenv("MAX_RANK"))
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- as.integer(Sys.getenv("MIN_COHORT_SIZE"))
CHOSEN_K <- as.integer(Sys.getenv("CHOSEN_K"))

clustered_panels <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank=2.parquet"))) |>
    .to_arrow_dataset() |>
    mutate(log_avg_weekly_earnings = log10(avg_weekly_earnings)) |>
    collect()

print(str_glue("Running match outcome decomposition for FIRST_YEAR={FIRST_YEAR}, LAST_YEAR={LAST_YEAR}, MAX_RANK={MAX_RANK}, YEAR_CLUSTER_SIZE={YEAR_CLUSTER_SIZE}, MIN_COHORT_SIZE={MIN_COHORT_SIZE}, CHOSEN_K={CHOSEN_K}"))

full_panel_df <- clustered_panels |>
    .to_arrow_dataset() |>
    filter(k == CHOSEN_K) |>
    collect()

panel <- UnbalancedPanel$new(full_panel_df,
    "worker", "outcome", "log_avg_weekly_earnings",
    model_rank = MAX_RANK, min_cohort_size = MIN_COHORT_SIZE,
    subset_to_largest_super_cohort = TRUE)

firm_clusters <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet")))
num_firms_by_outcome <- firm_clusters |>
    filter(k == CHOSEN_K) |>
    group_by(province, cluster) |>
    summarize(
        firms_for_outcome = n()
    ) |>
    inner_join(full_panel_df |> select(province, cluster, outcome) |> distinct(), by = c("province", "cluster")) |>
    ungroup() |>
    mutate(outcome_idx = panel$get_outcome_to_index()[outcome]) |>
    select(outcome_idx, firms_for_outcome) |>
    arrange(outcome_idx) |>
    pull(firms_for_outcome)

est_specs <- list(
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

wb <- get_weighted_bootstrap_draws(panel$get_num_units(), 1000, type = "bayesian", seed = 9001L)

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

target_param_table_data <- read_parquet(file.path(RESULTS_PATH, "match_outcome_decomp_results.parquet"))

ci_label <- "95\\% CI"
ci_sym <- rlang::sym(ci_label)

formatted_target_params <- target_param_table_data |>
    mutate(spec = factor(spec, levels = c("APM", "TWFE", "Difference"))) |>
    group_by(parameter, spec) |>
    summarize(
        estimate = first(estimate),
        std_error = first(std_error),
        p_value = first(p_value),
        ci_lb = first(ci_lb),
        ci_ub = first(ci_ub),
        .groups = "drop"
    ) |>
    mutate(
        `Point Est.` = formatC(estimate, format = "f", digits = 3),
        `Std. Err.` = formatC(std_error, format = "f", digits = 3),
        `$p$-value` = formatC(p_value, format = "f", digits = 3),
        !!ci_sym := str_glue("[{formatC(ci_lb, format = 'f', digits = 3)}, {formatC(ci_ub, format = 'f', digits = 3)}]")
    )

spec_order <- levels(formatted_target_params$spec)
province_stat_order <- c("Point Est.", "Std. Err.", "$p$-value", ci_label)

province_stats <- formatted_target_params |>
    filter(parameter == "Province")

extract_spec_values <- function(df, spec_name, stat_order) {
    spec_row <- df[df$spec == spec_name, , drop = FALSE]
    if (nrow(spec_row) == 0) {
        rep(NA_character_, length(stat_order))
    } else {
        as.character(unlist(spec_row[1, stat_order], use.names = FALSE))
    }
}

target_param_table <- tibble(
    Parameter = rep("Province", length(province_stat_order)),
    Statistic = province_stat_order
)

for (spec_name in spec_order) {
    target_param_table[[spec_name]] <- extract_spec_values(province_stats, spec_name, province_stat_order)
}

kbl_data <- target_param_table |>
    mutate(
        Statistic = if_else(Statistic == "Point Est.", "Province Share", Statistic)
    ) |>
    select(Statistic, all_of(spec_order))

target_param_table_kable <- kbl_data |>
    kbl(
        format = "latex",
        booktabs = TRUE,
        col.names = c("Statistic", spec_order),
        align = c("l", rep("c", length(spec_order))),
        escape = FALSE
    )

target_param_table_latex <- target_param_table_kable |>
    as.character() |>
    str_replace_all("\\{\\}\\[", "[")

writeLines(target_param_table_latex, file.path(TABLES_PATH, "match_outcome_decomp_results.tex"))