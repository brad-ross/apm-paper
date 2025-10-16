library(tidyverse)
library(arrow)

# Silence R CMD check NOTES for non-standard evaluation column names used with dplyr
utils::globalVariables(c(
  "cod_pgr", "matr_az", "anno", "sett_r", "gior_r", "retrib03", "prov_l",
  "cap", "comune", "prov", "csc", "province", "year", "paid_weeks",
  "paid_days", "earnings", "worker", "firm", "firm_num_years",
  "worker_num_years", "avg_weekly_earnings", "earnings_sum", "weeks_sum"
))

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
RAW_DATA_PATH <- file.path(Sys.getenv("DATA_PATH"), "raw_data")
CLEAN_DATA_PATH <- file.path(Sys.getenv("DATA_PATH"), "clean_data")
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

  total_years <- df |>
    distinct(year) |>
    tally(name = "total_years") |>
    collect() |>
    pull(total_years)

  firm_counts <- df |>
    distinct(firm, year) |>
    count(firm, name = "firm_num_years")

  worker_counts <- df |>
    distinct(worker, year) |>
    count(worker, name = "worker_num_years")

  out <- df |>
    left_join(firm_counts, by = "firm") |>
    left_join(worker_counts, by = "worker") |>
    filter(
      firm_num_years == total_years,
      worker_num_years == total_years
    ) |>
    select(-firm_num_years, -worker_num_years) |>
    filter(!is.na(avg_weekly_earnings))

  if (collect) out <- collect(out)
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
      earnings_sum = sum(earnings, na.rm = TRUE),
      weeks_sum    = sum(paid_weeks, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(avg_weekly_earnings = earnings_sum / weeks_sum) |>
    select(-earnings_sum, -weeks_sum) |>
    filter(!is.na(avg_weekly_earnings))

  if (collect) {
    old_opt <- options(arrow.int64_downcast = TRUE)
    on.exit(options(old_opt), add = TRUE)
    out <- collect(out)
  }

  out
}