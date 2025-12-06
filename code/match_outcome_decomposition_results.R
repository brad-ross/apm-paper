# Generate tables for match outcome decomposition results
#
# This script reads the parquet output from match_outcome_decomposition.R
# and generates LaTeX tables for the paper.

source("code/get_io_paths.R")

library(tidyverse)
library(arrow)
library(kableExtra)

# Read the results from parquet
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

print(str_glue("Table written to {file.path(TABLES_PATH, 'match_outcome_decomp_results.tex')}"))