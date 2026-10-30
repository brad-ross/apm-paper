# Public Figure 1a rendering: its only data input is the supplied plotting grid.
source("code/lib/get_io_paths.R")
source("code/lib/driver_panel_helpers.R")

grid_path <- file.path(PROCESSED_DATA_PATH, "driver_panel_grid.csv")
if (!file.exists(grid_path)) {
  stop("Incomplete replication package: missing Figure 1a plotting grid at ",
       grid_path, ". Restore the supplied data/processed_data/driver_panel_grid.csv; ",
       "public rendering does not regenerate it from private records.", call. = FALSE)
}
grid <- tryCatch(read.csv(grid_path, check.names = FALSE), error = function(e) {
  stop("Invalid Figure 1a plotting grid: cannot read CSV. Restore the supplied ",
       "data/processed_data/driver_panel_grid.csv.", call. = FALSE)
})
validate_driver_panel_grid(grid)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Figure 1a rendering requires ggplot2; run the public figure setup.", call. = FALSE)
}
panel_plot <- ggplot2::ggplot(grid) +
  ggplot2::theme_test() +
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     fill = "#7a0a0a") +
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#f5e7e7"),
                 panel.border = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3)) +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0), labels = function(x) {
    parse(text = ifelse(x == 0, "0",
      sub("e[+]?", " %*% 10^", format(x, scientific = TRUE, trim = TRUE))))
  }) +
  ggplot2::labs(x = "Week", y = "Driver")
ggplot2::ggsave(file.path(FIGURES_PATH, "panel_plot_drivers_times.pdf"),
                panel_plot, width = 3.0625, height = 3.25, units = "in")
