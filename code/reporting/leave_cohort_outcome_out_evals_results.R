# Generate tables for leave-cohort-outcome-out evaluation results
#
# This script reads the parquet output from leave_cohort_outcome_out_evals.R
# and generates comparison tables for the paper.

source("code/lib/vwh_data_helpers.R")
source("code/lib/get_io_paths.R")
source("code/lib/env_config.R")
source("code/lib/estimation_helpers.R")
source("code/lib/text_formatting_helpers.R")

library(tidyverse)
library(arrow)
library(kableExtra)

clustered_panels <- load_clustered_panels()
full_panel_df <- filter_panel_by_k(clustered_panels)

firm_clusters <- load_firm_clusters()
firms_per_cluster <- compute_firms_per_cluster(firm_clusters) |>
    filter(k == CHOSEN_K)
firm_cluster_sizes <- firms_per_cluster |>
    inner_join(full_panel_df |> select(province, cluster, outcome) |> distinct(), 
               by = c("province", "cluster")) |>
    rename(firms_for_outcome = num_firms) |>
    ungroup()

# We need panel to get outcome_to_index mapping - create a minimal version
# This requires the apm package to be loaded
library(apm)

panel <- create_estimation_panel(full_panel_df)

firm_cluster_sizes <- firm_cluster_sizes |>
    mutate(outcome_idx = panel$get_outcome_to_index()[outcome]) |>
    select(outcome_idx, firms_for_outcome)

est_specs <- get_default_est_specs(include_extended = TRUE)

masked_metrics_path <- file.path(RESULTS_PATH, str_glue("masked_metrics_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.parquet"))
masked_metrics <- read_parquet(masked_metrics_path) |>
    inner_join(firm_cluster_sizes, by = c("outcome" = "outcome_idx")) |>
    filter(spec %in% names(est_specs)) |>
    mutate(cohort_outcome_weight = cohort_size * firms_for_outcome)

long_masked_metrics <- masked_metrics |>
    pivot_longer(c(bias, se, rmse), names_to = "err_stat")

twfe_comparisons <- long_masked_metrics |>
    filter(spec != "twfe") |>
    inner_join(long_masked_metrics |> filter(spec == "twfe") |> select(cohort, outcome, err_stat, value), 
        by = c("cohort", "outcome", "err_stat"), suffix = c("_other", "_twfe")) |>
    mutate(value_diff = abs(value_other) - abs(value_twfe),
           value_ratio = abs(value_other) / abs(value_twfe))

est_specs_to_table <- tibble(
    spec = names(est_specs),
    r = sapply(est_specs, \(spec) spec$r),
    include_outcome_fes = sapply(est_specs, \(spec) spec$include_outcome_fes)
)

eval_table_data <- twfe_comparisons |>
    group_by(spec, err_stat) |>
    summarize(
        geom_mean_value_ratio = exp(sum(log(value_ratio)*cohort_outcome_weight)/sum(cohort_outcome_weight)),
        share_better = sum((value_ratio <= 1)*cohort_outcome_weight)/sum(cohort_outcome_weight),
        .groups = "drop"
    ) |>
    mutate(avg_pct_better = 1 - geom_mean_value_ratio) |>
    pivot_longer(-c(spec, err_stat), names_to = "metric") |>
    inner_join(est_specs_to_table, by = "spec")

spec_meta <- eval_table_data |>
    distinct(spec, r, include_outcome_fes) |>
    arrange(desc(include_outcome_fes), r, spec)

spec_levels <- spec_meta$spec

eval_table_wide <- eval_table_data |>
    filter(metric %in% c("geom_mean_value_ratio", "share_better")) |>
    mutate(
        err_stat = factor(err_stat, levels = c("bias", "se", "rmse")),
        metric = factor(metric, levels = c("geom_mean_value_ratio", "share_better"))
    ) |>
    arrange(err_stat, metric) |>
    select(err_stat, metric, spec, value) |>
    pivot_wider(names_from = spec, values_from = value) |>
    select(err_stat, metric, all_of(spec_levels))

table_for_display <- eval_table_wide |>
    mutate(
        err_stat_label = recode(err_stat,
            "bias" = "Bias",
            "se" = "Std. Err.",
            "rmse" = "RMSE"
        ),
        metric_label = recode(metric,
            "geom_mean_value_ratio" = "Geometric mean of statistic ratio",
            "share_better" = "Share with statistic ratio $\\leq 1$"
        )
    ) |>
    arrange(metric_label, err_stat_label)

metric_group_info <- table_for_display |>
    mutate(row_id = row_number()) |>
    group_by(metric_label) |>
    summarize(row_start = min(row_id), row_end = max(row_id), .groups = "drop")

row_metrics <- table_for_display$metric_label

kbl_data <- table_for_display |>
    select(err_stat_label, all_of(spec_levels))

value_matrix <- kbl_data |>
    select(-err_stat_label)

kbl_data <- kbl_data |>
    mutate(
        across(all_of(spec_levels), \(x) sprintf("%.2f", round(x, 2)))
    )

spec_cols <- spec_levels

geom_rows <- which(row_metrics == "Geometric mean of statistic ratio")
if (length(geom_rows) > 0) {
    for (row in geom_rows) {
        row_vals <- as.numeric(value_matrix[row, spec_cols])
        if (all(is.na(row_vals))) {
            next
        }
        min_val <- min(row_vals, na.rm = TRUE)
        min_idx <- which(row_vals == min_val)
        for (col_idx in min_idx) {
            col_name <- spec_cols[[col_idx]]
            kbl_data[row, col_name] <- str_glue("\\textbf{{{kbl_data[row, col_name]}}}")
        }
    }
}

share_rows <- which(row_metrics == "Share with statistic ratio $\\leq 1$")
if (length(share_rows) > 0) {
    for (row in share_rows) {
        row_vals <- as.numeric(value_matrix[row, spec_cols])
        if (all(is.na(row_vals))) {
            next
        }
        max_val <- max(row_vals, na.rm = TRUE)
        max_idx <- which(row_vals == max_val)
        for (col_idx in max_idx) {
            col_name <- spec_cols[[col_idx]]
            kbl_data[row, col_name] <- str_glue("\\textbf{{{kbl_data[row, col_name]}}}")
        }
    }
}

colnames(kbl_data) <- c("Stat.", str_glue("$r = {spec_meta$r}$"))

fes_header_rle <- rle(spec_meta$include_outcome_fes)
header_fes <- c(
    " " = 1,
    setNames(
        fes_header_rle$lengths,
        ifelse(fes_header_rle$values, "FEs", "No FEs")
    )
)

eval_table_latex <- kableExtra::kbl(
        kbl_data,
        format = "latex",
        booktabs = TRUE,
        linesep = "",
        digits = 2,
        escape = FALSE,
        align = "llrrr"
    ) |>
    kableExtra::add_header_above(header_fes, escape = FALSE)

for (i in seq_len(nrow(metric_group_info))) {
    eval_table_latex <- eval_table_latex |>
        kableExtra::group_rows(
            metric_group_info$metric_label[[i]],
            metric_group_info$row_start[[i]],
            metric_group_info$row_end[[i]],
            escape = FALSE
        )
}

eval_table_latex <- eval_table_latex |>
    as.character()

writeLines(eval_table_latex, file.path(TABLES_PATH, "leave_cohort_outcome_out_evals_results.tex"))

print(str_glue("Table written to {file.path(TABLES_PATH, 'leave_cohort_outcome_out_evals_results.tex')}"))

# ============================================================================
# Output single-statistic result snippets for evaluation metrics
# ============================================================================

# Output each cell value from the evaluation table
# Format: {spec}_{err_stat}_{metric}.txt
for (i in seq_len(nrow(eval_table_data))) {
    row <- eval_table_data[i, ]
    filename <- str_glue("{row$spec}_{row$err_stat}_{row$metric}.txt")
    write_result_snippet(format_pct(row$value), filename)
}

print(str_glue("Evaluation result snippets saved to {SNIPPETS_PATH}"))
