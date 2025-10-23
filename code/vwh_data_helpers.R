library(tidyverse)
library(arrow)

# Internal helper: coerce inputs to an Arrow Dataset for consistent lazy execution
.to_arrow_dataset <- function(x) {
  # Leave Arrow objects as-is
  if (inherits(x, c("Dataset", "ArrowDplyrQuery", "Table", "RecordBatchReader", "ArrowTabular"))) {
    return(x)
  }
  # Wrap tibbles/data.frames as in-memory Arrow Table (lazy dplyr supported)
  if (inherits(x, c("data.frame", "tbl_df"))) {
    return(arrow::arrow_table(x))
  }
  x
}

# Paths and constants
DATA_PATH <- Sys.getenv("DATA_PATH")
RAW_DATA_PATH <- file.path(DATA_PATH, "raw_data")
CLEAN_DATA_PATH <- file.path(DATA_PATH, "clean_data")
VENETO_PROVINCES <- c("BL", "PD", "RO", "TV", "VE", "VR", "VI")

# Read raw earnings as Arrow-backed query or tibble
read_raw_earnings_data <- function(collect = FALSE) {
  ds <- arrow::open_dataset(file.path(RAW_DATA_PATH, "earnings.parquet")) |>
    transmute(
      worker     = arrow::cast(cod_pgr, arrow::int64()),
      firm       = arrow::cast(matr_az, arrow::int64()),
      year       = as.integer(anno),
      paid_weeks = as.integer(sett_r),
      paid_days  = as.integer(gior_r),
      earnings   = as.numeric(retrib03),
      province   = if_else(prov_l == "", NA_character_, prov_l)
    )
  if (collect) {
    old_opt <- options(arrow.int64_downcast = TRUE)
    on.exit(options(old_opt), add = TRUE)
    ds <- collect(ds)
  }
  ds
}

# Read raw firm data as Arrow-backed query or tibble
read_raw_firm_data <- function(collect = FALSE) {
  ds <- arrow::open_dataset(file.path(RAW_DATA_PATH, "firms.parquet")) |>
    transmute(
      firm          = arrow::cast(matr_az, arrow::int64()),
      postal_code   = arrow::cast(cap, arrow::int64()),
      city          = comune,
      province      = prov,
      taxation_code = arrow::cast(csc, arrow::int64())
    )
  if (collect) {
    old_opt <- options(arrow.int64_downcast = TRUE)
    on.exit(options(old_opt), add = TRUE)
    ds <- collect(ds)
  }
  ds
}

# Join earnings and firms, filter by year/province and positive weeks
filter_and_join_match_data <- function(raw_earnings, raw_firms, min_year, max_year, 
    provinces = NULL, collect = FALSE) {
  raw_earnings <- .to_arrow_dataset(raw_earnings)
  raw_firms    <- .to_arrow_dataset(raw_firms)

  out <- raw_earnings |>
    select(-province) |>
    left_join(raw_firms, by = "firm") |>
    filter(
      year >= min_year, year <= max_year,
      (is.null(provinces) | (!is.na(province) & province %in% provinces)),
      paid_weeks > 0
    )
  if (collect) out <- collect(out)
  out
}

# Compute start year bins for a year column in a dataset/table
group_years <- function(df, years_col, bin_width, collect = FALSE) {
  stopifnot(is.character(years_col), length(years_col) == 1)

  # Ensure Arrow-backed dataset for lazy execution
  df <- .to_arrow_dataset(df)

  col_sym <- rlang::sym(years_col)

  # Compute dataset-wide minimum year in Arrow, pull scalar
  min_year <- df |>
    summarise(min_year = min(!!col_sym, na.rm = TRUE)) |>
    collect() |>
    pull(min_year)

  out <- df |>
    mutate(!!col_sym := as.integer(min_year + floor((!!col_sym - min_year) / bin_width) * bin_width))

  if (collect) {
    old_opt <- options(arrow.int64_downcast = TRUE)
    on.exit(options(old_opt), add = TRUE)
    out <- collect(out)
  }

  out
}

# Keep only workers and firms present in all years
filter_to_always_present_workers_and_firms <- function(df, collect = FALSE) {
  df <- .to_arrow_dataset(df)
  stopifnot("year" %in% names(df))

  # Number of distinct years that define "always present"
  total_years <- df |>
    distinct(year) |>
    tally(name = "total_years") |>
    collect() |>
    pull(total_years)

  prev_n_workers <- Inf
  prev_n_firms   <- Inf
  iter <- 0L
  max_iter <- 100L

  repeat {
    iter <- iter + 1L
    print(str_glue("Iteration {iter} with {prev_n_workers} workers and {prev_n_firms} firms"))
    # Firms present in all years given current subset
    always_present_firms <- df |>
      distinct(firm, province, year) |>
      count(firm, province, name = "firm_num_years") |>
      filter(firm_num_years == total_years) |>
      select(firm, province)

    # Restrict to those firms
    df_new <- df |>
      inner_join(always_present_firms, by = c("firm", "province"))

    # Workers present in all years given current subset
    always_present_workers <- df_new |>
      distinct(worker, year) |>
      count(worker, name = "worker_num_years") |>
      filter(worker_num_years == total_years) |>
      select(worker)

    # Restrict to those workers as well
    df_new <- df_new |>
      inner_join(always_present_workers, by = "worker") |>
      collect()

    # Convergence check: sizes of kept sets stop changing
    n_firms <- always_present_firms |>
      summarise(n = dplyr::n()) |>
      collect() |>
      dplyr::pull(n)
    n_workers <- always_present_workers |>
      summarise(n = dplyr::n()) |>
      collect() |>
      dplyr::pull(n)

    df <- .to_arrow_dataset(df_new)
    if ((n_firms == prev_n_firms && n_workers == prev_n_workers) || iter >= max_iter) break
    prev_n_firms   <- n_firms
    prev_n_workers <- n_workers
  }

  out <- df |>
    filter(!is.na(earnings))

  if (collect) out <- collect(out)

  gc()

  out
}

# Average weekly earnings by worker and additional grouping variables
constr_avg_weekly_wages_by_group <- function(df, group_vars, collect = FALSE) {
  # Ensure character vector for tidy-eval
  group_vars <- as.character(group_vars)

  # Ensure Arrow-backed dataset for lazy execution
  df <- .to_arrow_dataset(df)

  out <- df |>
    # Use tidy-eval to build grouping compatible with Arrow translation
    group_by(!!!rlang::syms(c("worker", group_vars))) |>
    summarise(
      earnings = sum(earnings, na.rm = TRUE),
      paid_weeks    = sum(paid_weeks, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(avg_weekly_earnings = earnings / paid_weeks) |>
    filter(!is.na(avg_weekly_earnings))

  if (collect) {
    old_opt <- options(arrow.int64_downcast = TRUE)
    on.exit(options(old_opt), add = TRUE)
    out <- collect(out)
  }

  out
}