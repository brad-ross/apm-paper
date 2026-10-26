# Figure 1a's public input contains rectangle bounds only, never driver records.
driver_grid_columns <- c("xmin", "xmax", "ymin", "ymax")

validate_driver_panel_grid <- function(grid) {
  invalid <- function(reason) {
    stop("Invalid Figure 1a plotting grid: ", reason,
         ". Restore the supplied data/clean_data/driver_panel_grid.csv; ",
         "public rendering does not regenerate it from private records.",
         call. = FALSE)
  }
  if (!identical(names(grid), driver_grid_columns)) {
    invalid("expected exactly the columns xmin, xmax, ymin, ymax")
  }
  if (!nrow(grid)) invalid("no cohort rectangles")
  if (!all(vapply(grid, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(grid)))) {
    invalid("all bounds must be finite numbers")
  }
  if (any(as.matrix(grid) * 2 != round(as.matrix(grid) * 2))) {
    invalid("bounds must lie on integer or half-integer coordinates")
  }
  if (any(grid$xmax <= grid$xmin) ||
      any((grid$xmax - grid$xmin) != round(grid$xmax - grid$xmin))) {
    invalid("week intervals must have positive whole-number widths")
  }
  sizes <- grid$ymax - grid$ymin
  if (any(sizes < 100 | sizes != round(sizes))) {
    invalid("each cohort must contain at least 100 whole display rows")
  }
  if (grid$ymin[1] != 0.5 ||
      any(grid$ymin != c(0.5, head(grid$ymax, -1)))) {
    invalid("anonymous display rows must be contiguous and start at row 1")
  }
  invisible(grid)
}
