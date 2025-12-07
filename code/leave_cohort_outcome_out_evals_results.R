# Generate tables for leave-cohort-outcome-out evaluation results
#
# This script reads the parquet output from leave_cohort_outcome_out_evals.R
# and generates comparison tables for the paper.

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")
source("code/estimation_helpers.R")
source("code/text_formatting_helpers.R")

library(tidyverse)
library(arrow)
library(kableExtra)

# Weighted quantile function
# Computes weighted quantiles with optional interpolation
wtd_quantile <- function(x, weights = NULL, probs = c(0.25, 0.5, 0.75),
                         na.rm = FALSE, interpolate = FALSE, names = TRUE) {
    if (is.null(weights)) weights <- rep(1, length(x))
    if (!is.numeric(x)) stop("`x` must be numeric.")
    if (length(weights) != length(x)) stop("`weights` must have same length as `x`.")
    if (any(!is.finite(weights))) stop("`weights` must be finite.")
    if (any(weights < 0)) stop("`weights` must be non-negative.")
    if (!is.numeric(probs) || anyNA(probs) || any(probs < 0 | probs > 1))
        stop("`probs` must be numeric in [0, 1].")

    if (na.rm) {
        keep <- is.finite(x) & is.finite(weights) & !is.na(x) & !is.na(weights)
        x <- x[keep]; weights <- weights[keep]
    } else {
        if (anyNA(x) || anyNA(weights) || any(!is.finite(x)))
            stop("Missing/non-finite values present; set `na.rm = TRUE` to drop them.")
    }
    if (!length(x)) return(rep(NA_real_, length(probs)))

    # Drop zero-weight observations
    keep <- weights > 0
    x <- x[keep]; weights <- weights[keep]
    if (!length(x)) stop("All weights are zero after filtering.")

    # Order by x
    ord <- order(x)
    x <- x[ord]
    w <- weights[ord]

    # Weighted CDF
    cdf <- cumsum(w) / sum(w)

    # For each p, j = number of cdf values <= p
    j <- findInterval(probs, cdf, rightmost.closed = TRUE)
    n <- length(x)

    if (!interpolate) {
        # Step-function inverse CDF: smallest x with CDF >= p
        idx <- pmin(j + 1L, n)
        q <- x[idx]
    } else {
        # Piecewise linear interpolation between adjacent x's
        left_cdf  <- ifelse(j == 0, 0, cdf[j])
        right_cdf <- cdf[pmin(j + 1L, n)]
        left_x    <- ifelse(j == 0, x[1], x[j])
        right_x   <- x[pmin(j + 1L, n)]
        t <- (probs - left_cdf) / pmax(right_cdf - left_cdf, .Machine$double.eps)
        q <- left_x + t * (right_x - left_x)
    }

    if (names) names(q) <- paste0(format(100 * probs, trim = TRUE), "%")
    q
}

# Load panel data to get outcome indices
clustered_panels <- load_clustered_panels()
full_panel_df <- filter_panel_by_k(clustered_panels)

# Load firm clusters and compute sizes
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
devtools::load_all("../apm-package/r")

panel <- create_estimation_panel(full_panel_df)

firm_cluster_sizes <- firm_cluster_sizes |>
    mutate(outcome_idx = panel$get_outcome_to_index()[outcome]) |>
    select(outcome_idx, firms_for_outcome)

# Get estimation specifications
est_specs <- get_default_est_specs(include_extended = TRUE)

# Read masked metrics
masked_metrics_path <- file.path(RESULTS_PATH, str_glue("masked_metrics_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.parquet"))
masked_metrics <- read_parquet(masked_metrics_path) |>
    inner_join(firm_cluster_sizes, by = c("outcome" = "outcome_idx")) |>
    filter(spec %in% names(est_specs)) |>
    mutate(cohort_outcome_weight = cohort_size * firms_for_outcome)

# Get 15 largest cohorts for whom we have masked metrics
largest_cohorts <- masked_metrics |>
    select(cohort, cohort_size) |>
    distinct() |>
    arrange(desc(cohort_size)) |>
    pull(cohort) |>
    head(15)

masked_metrics <- masked_metrics |>
    filter(cohort %in% largest_cohorts)

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
    pivot_longer(-c(spec, err_stat), names_to = "metric") |>
    inner_join(est_specs_to_table, by = "spec")

spec_meta <- eval_table_data |>
    distinct(spec, r, include_outcome_fes) |>
    arrange(desc(include_outcome_fes), r, spec)

spec_levels <- spec_meta$spec

eval_table_wide <- eval_table_data |>
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

colnames(kbl_data) <- c("Statistic", str_glue("$r = {spec_meta$r}$"))

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
    write_result_snippet(format_decimal(row$value), filename)
}

print(str_glue("Evaluation result snippets saved to {SNIPPETS_PATH}"))