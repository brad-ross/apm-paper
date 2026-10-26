# Private, explicitly requested preparation step; public rendering never sources it.
# Recovered from eed1a22:code_old/visualize_panels.R:279-320.
source("code/get_io_paths.R")
source("code/driver_panel_helpers.R")

input_path <- file.path(RAW_DATA_PATH, "riders.csv")
if (!file.exists(input_path)) {
  stop("Private Figure 1a preparation requires riders.csv in RAW_DATA_PATH.",
       call. = FALSE)
}
required_columns <- c("rider_id", "first_week_in_panel_dateweek_id", "began_active_week_id")
# Read identifiers as strings, and never write or print them.
drivers <- read.csv(input_path, colClasses = "character", check.names = FALSE)
if (!all(required_columns %in% names(drivers))) {
  stop("Private riders.csv is missing required Figure 1a columns.", call. = FALSE)
}
drivers <- drivers[required_columns]
drivers <- drivers[complete.cases(drivers), , drop = FALSE]
if (!nrow(drivers)) stop("No complete driver records for Figure 1a.", call. = FALSE)
if (any(!nzchar(trimws(drivers$rider_id)))) {
  stop("Private riders.csv contains an empty driver identifier.", call. = FALSE)
}
first_week <- suppressWarnings(as.numeric(drivers$first_week_in_panel_dateweek_id))
treatment_week <- suppressWarnings(as.numeric(drivers$began_active_week_id))
if (any(!is.finite(first_week) | !is.finite(treatment_week)) ||
    any(first_week != round(first_week) | treatment_week != round(treatment_week)) ||
    any(first_week > treatment_week)) {
  stop("Driver weeks must be finite integers with entry no later than treatment.",
       call. = FALSE)
}
drivers <- unique(data.frame(rider_id = drivers$rider_id,
                             first_week = first_week, treatment_week = treatment_week))
if (anyDuplicated(drivers$rider_id)) {
  stop("A driver has conflicting entry/treatment cohorts; resolve private input.",
       call. = FALSE)
}
# Historical str_pad(width = 2) plus arrange(cohort) sets the vertical order.
pad_week <- function(week) {
  value <- as.character(week)
  ifelse(nchar(value) < 2, paste0("0", value), value)
}
drivers$cohort <- paste(pad_week(drivers$first_week),
                        pad_week(drivers$treatment_week), sep = " - ")
cohorts <- aggregate(list(cohort_size = rep(1, nrow(drivers))),
                     drivers[c("cohort", "first_week", "treatment_week")], sum)
cohorts <- cohorts[cohorts$cohort_size >= 100, , drop = FALSE]
cohorts <- cohorts[order(cohorts$cohort), , drop = FALSE]
if (!nrow(cohorts)) stop("No driver cohort meets the minimum size of 100.", call. = FALSE)
# geom_tile's default width is the smallest spacing of observed week centers.
# Inclusive multiweek intervals guarantee spacing 1; singleton-only panels may not.
week_centers <- sort(unique(cohorts$first_week))
tile_width <- if (any(cohorts$treatment_week > cohorts$first_week) ||
                  length(week_centers) < 2) 1 else min(diff(week_centers))
last_row <- cumsum(cohorts$cohort_size)
grid <- data.frame(xmin = cohorts$first_week - tile_width / 2,
                   xmax = cohorts$treatment_week + tile_width / 2,
                   ymin = c(0, head(last_row, -1)) + 0.5,
                   ymax = last_row + 0.5)
validate_driver_panel_grid(grid)
grid_path <- file.path(CLEAN_DATA_PATH, "driver_panel_grid.csv")
write.csv(grid, grid_path, row.names = FALSE)
message("Saved Figure 1a cohort rectangles to ", grid_path)
