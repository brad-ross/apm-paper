devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")

library(kableExtra)

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
CHOSEN_K <- as.integer(Sys.getenv("CHOSEN_K"))

clustered_panels <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank=2.parquet"))) |>
    .to_arrow_dataset() |>
    mutate(log_avg_weekly_earnings = log10(avg_weekly_earnings)) |>
    collect()

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
    "pc_1_no_fes_equal_weights" = list(
        factor_model_estimator = "principal_components",
        include_outcome_fes = FALSE,
        r = 1,
        cohort_weighting = "equal"
    ),
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
)

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

wb <- get_weighted_bootstrap_draws(panel$get_num_units(), 1000, type = "bayesian", seed = 9001L)

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

        # break

        # Checkpoint after every (cohort, outcome) pair
        write_parquet(masked_metrics, masked_metrics_path)
    }

    # if (!is.null(masked_metrics) && nrow(masked_metrics) > 0) {
    #     break
    # }
}

firm_clusters <- read_parquet(file.path(PROCESSED_DATA_PATH, str_glue("firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet")))
firm_cluster_sizes <- firm_clusters |>
    filter(k == CHOSEN_K) |>
    group_by(province, cluster) |>
    summarize(
        firms_for_outcome = n()
    ) |>
    inner_join(full_panel_df |> select(province, cluster, outcome) |> distinct(), by = c("province", "cluster")) |>
    ungroup() |>
    mutate(outcome_idx = panel$get_outcome_to_index()[outcome]) |>
    select(outcome_idx, firms_for_outcome)

masked_metrics <- read_parquet(masked_metrics_path) |>
    inner_join(firm_cluster_sizes, by = c("outcome" = "outcome_idx")) |>
    filter(spec %in% names(est_specs)) |>
    mutate(cohort_outcome_weight = cohort_size * firms_for_outcome)

# get 15 largest cohorts for whom we have masked metrics
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
        # ratio_q10 = wtd_quantile(value_ratio, 0.1, weights = cohort_outcome_weight),
        # ratio_q25 = wtd_quantile(value_ratio, 0.25, weights = cohort_outcome_weight),
        # ratio_q50 = wtd_quantile(value_ratio, 0.5, weights = cohort_outcome_weight),
        # ratio_q75 = wtd_quantile(value_ratio, 0.75, weights = cohort_outcome_weight),
        # ratio_q90 = wtd_quantile(value_ratio, 0.9, weights = cohort_outcome_weight),
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
    kableExtra::kable_styling(full_width = FALSE, latex_options = c("hold_position"))

print(eval_table_latex)