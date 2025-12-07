# Generate figures for clustered panel analysis
#
# This script reads the parquet output from construct_clustered_panels.R
# and generates figures for the paper.

devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")
source("code/estimation_helpers.R")
source("code/text_formatting_helpers.R")

library(tidyverse)
library(arrow)
library(kableExtra)

# Read clustering data
firm_clusters <- load_firm_clusters()

panels_with_clustered_outcomes <- load_clustered_panels()

# ============================================================================
# Code to plot distributions of cohort sizes and number of cohorts across clusterings
# ============================================================================

cohorts_by_k <- panels_with_clustered_outcomes |>
    group_by(k) |>
    group_modify(\(panel_with_k_clustered_outcomes, keys) {
        k <- keys$k[1]
        print(paste("# clusters:", k))

        clustered_panel_cohorts <- construct_cohorts_from_panel(
            panel_with_k_clustered_outcomes,
            "worker", "outcome", "avg_weekly_earnings",
            model_rank = 1,
            min_cohort_size = MIN_COHORT_SIZE,
            subset_to_largest_super_cohort = FALSE
        )

        cohort_sizes <- clustered_panel_cohorts$cohort_sizes

        data.frame(
            cohort_size = cohort_sizes
        )
    })

cohort_counts_by_k <- cohorts_by_k |>
    ungroup() |>
    count(k, name = "num_cohorts")

cohort_counts_named <- setNames(cohort_counts_by_k$num_cohorts, as.character(cohort_counts_by_k$k))

cohort_size_dist_and_num_cohorts_by_k_plot <- (ggplot(cohorts_by_k |> mutate(k = as.factor(k)))
    + theme_bw()
    + geom_freqpoly(aes(x = cohort_size, color = k), size = 1.2, bins = 30) 
    + scale_color_viridis_d(end = 0.9, labels = function(x) paste0(x, " (", cohort_counts_named[as.character(x)], ")"))
    + scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(MIN_COHORT_SIZE, NA))
    + labs(
        x = "# Workers in Cohort",
        y = "# Cohorts",
        color = "# Clusters per Province\n(Total # Cohorts)",
    )
    + guides(color = guide_legend(ncol = 3))
)

ggsave(file.path(FIGURES_PATH, str_glue("cohort_size_dist_and_num_cohorts_by_k_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.pdf")), cohort_size_dist_and_num_cohorts_by_k_plot, width = PAPER_FIG_WIDTH, height = PAPER_FIG_HEIGHT)

# ============================================================================
# Code to plot the largest super cohort size across clusterings for a given rank
# ============================================================================

firms_per_firm_cluster <- compute_firms_per_cluster(firm_clusters)

outcomes_per_outcome_cluster_df <- build_outcomes_with_firm_counts_df(
    panels_with_clustered_outcomes, 
    firms_per_firm_cluster
)

comp_largest_super_cohort_sizes_across_clusterings <- function(panels_with_clustered_outcomes, rank) {
    panels_with_clustered_outcomes |>
    group_by(k) |>
    group_modify(\(panel_with_k_clustered_outcomes, keys) {
        curr_k <- keys$k[1]
        num_units <- length(unique(panel_with_k_clustered_outcomes$worker))

        print(paste("# clusters:", curr_k, "; num units:", num_units))

        clustered_panel_cohorts <- construct_cohorts_from_panel(
            panel_with_k_clustered_outcomes,
            "worker", "outcome", "avg_weekly_earnings",
            model_rank = rank,
            min_cohort_size = MIN_COHORT_SIZE,
            subset_to_largest_super_cohort = FALSE
        )

        ooi <- clustered_panel_cohorts$observed_outcome_indices
        cohort_sizes <- clustered_panel_cohorts$cohort_sizes
        num_cohorts <- length(cohort_sizes)

        outcomes_per_outcome_cluster <- get_outcome_weights_for_k(
            outcomes_per_outcome_cluster_df, curr_k, clustered_panel_cohorts$outcome_ids)

        num_outcomes <- length(clustered_panel_cohorts$outcome_ids)
        total_outcome_weight <- sum(outcomes_per_outcome_cluster)

        id_summary <- summarize_identification(ooi, cohort_sizes, rank, outcome_weights = outcomes_per_outcome_cluster)

        largest_super_cohort_df <- NULL
        for (iter in 1:id_summary$num_o3_iterations) {
            iter_id_summary <- summarize_identification(ooi, cohort_sizes, rank, 
                iter = iter, outcome_weights = outcomes_per_outcome_cluster)

            print(str_glue("O³ iteration: {iter}; largest super cohort share workers: {iter_id_summary$largest_super_cohort_size/num_units}"))
            print(str_glue("O³ iteration: {iter}; largest super cohort size firms: {iter_id_summary$total_outcome_weight_in_largest_super_cohort}"))
            print(str_glue("O³ iteration: {iter}; largest super cohort share firm factors identified: {iter_id_summary$total_outcome_weight_in_largest_super_cohort/total_outcome_weight}"))

            new_largest_super_cohort_df <- data.frame(
                o3_iter = as.character(iter - 1),
                num_units = num_units,
                num_cohorts = num_cohorts,
                num_outcomes = num_outcomes,
                num_firms = total_outcome_weight,
                num_units_largest_super_cohort = iter_id_summary$largest_super_cohort_size,
                largest_super_cohort_share = iter_id_summary$largest_super_cohort_size/num_units,
                min_cohort_size_in_largest_super = iter_id_summary$min_cohort_size_in_largest_super,
                num_outcomes_in_largest_super_cohort = iter_id_summary$num_outcomes_in_largest_super_cohort,
                total_outcome_weight_in_largest_super_cohort = iter_id_summary$total_outcome_weight_in_largest_super_cohort,
                share_outcomes_in_largest_super_cohort = iter_id_summary$share_outcomes_in_largest_super_cohort,
                share_outcome_weight_in_largest_super_cohort = iter_id_summary$share_outcome_weight_in_largest_super_cohort
            )

            largest_super_cohort_df <- bind_rows(largest_super_cohort_df, new_largest_super_cohort_df)
        }

        largest_super_cohort_df
    })
}

comp_block_missingness_identification_across_clusterings <- function(panels_with_clustered_outcomes, rank) {
    panels_with_clustered_outcomes |>
        group_by(k) |>
        group_modify(\(panel_with_k_clustered_outcomes, keys) {
            curr_k <- keys$k[1]
            num_units <- length(unique(panel_with_k_clustered_outcomes$worker))

            print(paste("# clusters:", curr_k, "; num units:", num_units))

            clustered_panel_cohorts <- construct_cohorts_from_panel(
                panel_with_k_clustered_outcomes,
                "worker", "outcome", "avg_weekly_earnings",
                model_rank = rank,
                min_cohort_size = MIN_COHORT_SIZE,
                subset_to_largest_super_cohort = FALSE
            )

            cohort_sizes <- clustered_panel_cohorts$cohort_sizes

            outcomes_per_outcome_cluster <- get_outcome_weights_for_k(
                outcomes_per_outcome_cluster_df, curr_k, clustered_panel_cohorts$outcome_ids)

            num_outcomes <- length(clustered_panel_cohorts$outcome_ids)
            total_outcome_weight <- sum(outcomes_per_outcome_cluster)

            outcomes_with_rank_overlap <- count_outcomes_with_rank_overlap_per_cohort(
                clustered_panel_cohorts$observed_outcome_indices, rank)
            total_outcome_weight_with_rank_overlap <- count_outcomes_with_rank_overlap_per_cohort(
                clustered_panel_cohorts$observed_outcome_indices, rank,
                outcome_weights = outcomes_per_outcome_cluster)

            data.frame(
                num_units = num_units,
                num_cohorts = length(cohort_sizes),
                num_outcomes = num_outcomes,
                num_firms = total_outcome_weight,
                cohort_size = cohort_sizes,
                outcomes_with_rank_overlap = outcomes_with_rank_overlap,
                total_outcome_weight_with_rank_overlap = total_outcome_weight_with_rank_overlap
            )
        })
}

create_largest_super_cohort_plot <- function(super_df,
    connected_df,
    y_var,
    y_label,
    linetype_scale_labels,
    linetype_scale_values,
    show_legend = TRUE) {

    y_var <- rlang::ensym(y_var)

    final_iteration <- super_df |>
        group_by(k) |>
        summarize(
            y_value = (!!y_var)[which.max(as.numeric(o3_iter))],
            .groups = "drop"
        )

    initial_iteration <- super_df |>
        filter(as.numeric(o3_iter) == 0)

    plot <- (ggplot() +
        theme_bw() +
        geom_line(data = final_iteration, mapping = aes(x = k, y = y_value, linetype = linetype_scale_labels[1])) +
        geom_line(data = connected_df, mapping = aes(x = k, y = !!y_var, linetype = linetype_scale_labels[2])) +
        geom_line(data = initial_iteration, mapping = aes(x = k, y = !!y_var, linetype = linetype_scale_labels[3])) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_linetype_manual(values = linetype_scale_values) +
        labs(
            x = "Number of Firm Clusters per Province",
            y = y_label,
            linetype = NULL
        )) 
        
    if (show_legend) {
        plot <- plot + guides(
            linetype = guide_legend(
                order = 1,
                keyheight = grid::unit(1.5, "lines"),
                label.theme = element_text(margin = margin(t = 4, b = 4))
            )
        ) +
        theme(legend.spacing.y = grid::unit(0.1, "lines"))
    } else {
        plot <- plot + theme(legend.position = "none")
    }

    plot
}

create_block_missingness_ecdf_plot <- function(df, x_expr, x_label) {
    x_expr <- rlang::enquo(x_expr)

    ggplot(df, aes(x = !!x_expr, color = as.factor(k))) +
        theme_bw() +
        stat_ecdf(aes(weight = cohort_size)) +
        scale_color_viridis_d(end = 0.9) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(
            x = x_label,
            y = "Share of Workers",
            color = "Number of Clusters\nper Province"
        ) +
        guides(color = guide_legend(ncol = 3))
}

largest_connected_component_sizes_across_clusterings <- comp_largest_super_cohort_sizes_across_clusterings(panels_with_clustered_outcomes, 1) |>
    filter(o3_iter == 1) |>
    mutate(o3_iter = "r=1 final iteration")

for (rank in 2:MAX_RANK) {
    print(paste("Rank:", rank))

    largest_super_cohort_sizes_across_clusterings <- comp_largest_super_cohort_sizes_across_clusterings(panels_with_clustered_outcomes, rank)

    linetype_scale_values <- c("solid", "dashed", "dotted")
    linetype_scale_labels <- c(str_glue("Final Largest\nSuper Cohort\n(r = {rank})"), "Largest\nConnected\nComponent\n(r = 1)", "Initial Largest\nCohort")
    linetype_scale_values <- setNames(linetype_scale_values, linetype_scale_labels)

    largest_super_cohort_share_workers_plot <- create_largest_super_cohort_plot(
        super_df = largest_super_cohort_sizes_across_clusterings,
        connected_df = largest_connected_component_sizes_across_clusterings,
        y_var = largest_super_cohort_share,
        y_label = "Share of Workers",
        linetype_scale_labels = linetype_scale_labels,
        linetype_scale_values = linetype_scale_values,
        show_legend = FALSE
    )
    ggsave(
        file.path(
            FIGURES_PATH,
            str_glue("largest_super_cohort_share_workers_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")
        ),
        largest_super_cohort_share_workers_plot,
        width = PAPER_FIG_HALF_WIDTH,
        height = PAPER_FIG_HEIGHT
    )

    largest_super_cohort_share_firms_plot <- create_largest_super_cohort_plot(
        super_df = largest_super_cohort_sizes_across_clusterings,
        connected_df = largest_connected_component_sizes_across_clusterings,
        y_var = share_outcome_weight_in_largest_super_cohort,
        y_label = "Share of Firms",
        linetype_scale_labels = linetype_scale_labels,
        linetype_scale_values = linetype_scale_values,
        show_legend = TRUE
    )
    ggsave(
        file.path(
            FIGURES_PATH,
            str_glue("largest_super_cohort_share_firms_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")
        ),
        largest_super_cohort_share_firms_plot,
        width = PAPER_FIG_HALF_WIDTH,
        height = PAPER_FIG_HEIGHT
    )

    num_o3_iterations_needed_plot <- (ggplot(largest_super_cohort_sizes_across_clusterings |>
        group_by(k) |>
        summarize(num_o3_iterations = n()))
        + theme_bw()
        + geom_line(aes(x = k, y = num_o3_iterations))
        + xlim(c(min(largest_super_cohort_sizes_across_clusterings$k), max(largest_super_cohort_sizes_across_clusterings$k)))
        + labs(
            x = "# of Firm Clusters per Province",
            y = "# of O³ Iterations until Convergence",
        ))
    ggsave(file.path(FIGURES_PATH, str_glue("num_o3_iterations_needed_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), num_o3_iterations_needed_plot, width = PAPER_FIG_WIDTH, height = PAPER_FIG_HEIGHT)

    # Plot of distribution of share of cohort outcome means identified by block missingness across cohorts

    block_missingness_identification_across_clusterings <- comp_block_missingness_identification_across_clusterings(
        panels_with_clustered_outcomes,
        rank
    )

    block_missingness_id_plot <- create_block_missingness_ecdf_plot(
        df = block_missingness_identification_across_clusterings,
        x_expr = outcomes_with_rank_overlap / num_outcomes,
        x_label = "Share of Cohort's Outcome Means Identified by Block Missingness"
    )

    ggsave(file.path(FIGURES_PATH, str_glue("block_missingness_id_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), block_missingness_id_plot, width = PAPER_FIG_WIDTH, height = PAPER_FIG_HEIGHT)

    firm_weighted_block_missingness_id_plot <- create_block_missingness_ecdf_plot(
        df = block_missingness_identification_across_clusterings,
        x_expr = total_outcome_weight_with_rank_overlap / num_firms,
        x_label = "Share of Cohort's Avg. Wages Across Firms Identified by Block Missingness"
    )

    ggsave(file.path(FIGURES_PATH, str_glue("block_missingness_id_firm_weighted_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), firm_weighted_block_missingness_id_plot, width = PAPER_FIG_WIDTH, height = PAPER_FIG_HEIGHT)
}

# ============================================================================
# Generate summary statistics table comparing panel subsets
# ============================================================================

# Load raw data for full panel statistics
raw_earnings <- read_raw_earnings_data()
raw_firms <- read_raw_firm_data()

# Column 1: Full panel (all data from FIRST_YEAR to LAST_YEAR)
# Collapse multiple spells to one observation per worker-firm-year cluster
full_panel <- filter_and_join_match_data(raw_earnings, raw_firms, FIRST_YEAR, LAST_YEAR, VENETO_PROVINCES) |>
    group_years("year", YEAR_CLUSTER_SIZE) |>
    constr_avg_weekly_wages_by_group(c("year", "province", "firm"), collect = TRUE) |>
    filter(avg_weekly_earnings > 0)

full_panel_stats <- list(
    num_observations = nrow(full_panel),
    num_workers = n_distinct(full_panel$worker),
    num_firms = n_distinct(full_panel$firm)
)

# Column 2: Always present workers and firms
# Collapse multiple spells to one observation per worker-firm-year cluster, then filter
always_present_panel <- filter_and_join_match_data(raw_earnings, raw_firms, FIRST_YEAR, LAST_YEAR, VENETO_PROVINCES) |>
    group_years("year", YEAR_CLUSTER_SIZE) |>
    constr_avg_weekly_wages_by_group(c("year", "province", "firm")) |>
    filter(avg_weekly_earnings > 0) |>
    filter_to_always_present_workers_and_firms(collect = TRUE)

always_present_stats <- list(
    num_observations = nrow(always_present_panel),
    num_workers = n_distinct(always_present_panel$worker),
    num_firms = n_distinct(always_present_panel$firm)
)

# Column 3: Largest super cohort for k = CHOSEN_K
panel_with_chosen_k <- panels_with_clustered_outcomes |>
    filter(k == CHOSEN_K)

cohorts_for_chosen_k <- construct_cohorts_from_panel(
    panel_with_chosen_k,
    "worker", "outcome", "avg_weekly_earnings",
    model_rank = MAX_RANK,
    min_cohort_size = MIN_COHORT_SIZE,
    subset_to_largest_super_cohort = TRUE
)

# Get the workers in the largest super cohort
workers_in_largest_super_cohort <- unique(cohorts_for_chosen_k$unit_cohorts$unit_id)

# Get the firms in the largest super cohort by matching outcomes back to firm clusters
outcomes_in_largest_super_cohort <- cohorts_for_chosen_k$outcome_ids
outcome_parts <- strsplit(outcomes_in_largest_super_cohort, ":")
outcome_df <- data.frame(
    year = as.integer(sapply(outcome_parts, `[[`, 1)),
    province = sapply(outcome_parts, `[[`, 2),
    cluster = as.integer(sapply(outcome_parts, `[[`, 3))
)

firms_in_largest_super_cohort <- firm_clusters |>
    filter(k == CHOSEN_K) |>
    inner_join(outcome_df, by = c("province", "cluster")) |>
    distinct(firm) |>
    pull(firm)

largest_super_cohort_panel <- panel_with_chosen_k |>
    filter(worker %in% workers_in_largest_super_cohort)

largest_super_cohort_stats <- list(
    num_observations = nrow(largest_super_cohort_panel),
    num_workers = length(workers_in_largest_super_cohort),
    num_firms = length(firms_in_largest_super_cohort)
)

# Compute share of workers with k observations
compute_obs_shares <- function(panel, worker_col = "worker") {
    obs_per_worker <- panel |>
        group_by(across(all_of(worker_col))) |>
        summarize(n_obs = n(), .groups = "drop")
    
    total_workers <- nrow(obs_per_worker)
    
    list(
        share_2_obs = sum(obs_per_worker$n_obs == 2) / total_workers,
        share_3_obs = sum(obs_per_worker$n_obs == 3) / total_workers,
        share_4_obs = sum(obs_per_worker$n_obs == 4) / total_workers,
        share_5plus_obs = sum(obs_per_worker$n_obs >= 5) / total_workers
    )
}

# Compute observation shares for each panel
full_panel_obs_shares <- compute_obs_shares(full_panel)
always_present_obs_shares <- compute_obs_shares(always_present_panel)
largest_super_cohort_obs_shares <- compute_obs_shares(largest_super_cohort_panel)

# Create the summary statistics table
summary_stats_df <- tibble(
    Statistic = c(
        "Num. Employment Spells",
        "Num. Workers",
        "Num. Firms",
        "\\% Workers w/ 2 Spells",
        "\\% Workers w/ 3 Spells",
        "\\% Workers w/ 4 Spells",
        "\\% Workers w/ 5+ Spells"
    ),
    `Full Panel` = c(
        format_num(full_panel_stats$num_observations),
        format_num(full_panel_stats$num_workers),
        format_num(full_panel_stats$num_firms),
        format_pct(full_panel_obs_shares$share_2_obs),
        format_pct(full_panel_obs_shares$share_3_obs),
        format_pct(full_panel_obs_shares$share_4_obs),
        format_pct(full_panel_obs_shares$share_5plus_obs)
    ),
    `Always-Present` = c(
        format_num(always_present_stats$num_observations),
        format_num(always_present_stats$num_workers),
        format_num(always_present_stats$num_firms),
        format_pct(always_present_obs_shares$share_2_obs),
        format_pct(always_present_obs_shares$share_3_obs),
        format_pct(always_present_obs_shares$share_4_obs),
        format_pct(always_present_obs_shares$share_5plus_obs)
    ),
    `Largest Super Cohort` = c(
        format_num(largest_super_cohort_stats$num_observations),
        format_num(largest_super_cohort_stats$num_workers),
        format_num(largest_super_cohort_stats$num_firms),
        format_pct(largest_super_cohort_obs_shares$share_2_obs),
        format_pct(largest_super_cohort_obs_shares$share_3_obs),
        format_pct(largest_super_cohort_obs_shares$share_4_obs),
        format_pct(largest_super_cohort_obs_shares$share_5plus_obs)
    )
)

# Generate LaTeX table
summary_stats_table_latex <- summary_stats_df |>
    kbl(
        format = "latex",
        booktabs = TRUE,
        align = c("l", "r", "r", "r"),
        escape = FALSE,
        linesep = ""
    ) |>
    as.character()

writeLines(summary_stats_table_latex, file.path(TABLES_PATH, 
    str_glue("panel_summary_stats_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.tex")))

print(str_glue("Summary statistics table saved to {TABLES_PATH}"))

# ============================================================================
# Dot plot of number of firms by province and cluster for k = CHOSEN_K
# ============================================================================

firms_by_province_cluster <- firms_per_firm_cluster |>
    filter(k == CHOSEN_K) |>
    mutate(
        province = factor(province, levels = sort(unique(province))),
        cluster = as.character(cluster),
        is_total = FALSE
    )

# Add province totals
province_totals <- firms_by_province_cluster |>
    group_by(province) |>
    summarize(num_firms = sum(num_firms), .groups = "drop") |>
    mutate(cluster = "Total", is_total = TRUE)

firms_by_province_cluster <- bind_rows(firms_by_province_cluster, province_totals) |>
    mutate(cluster = factor(cluster, levels = c(as.character(sort(unique(as.integer(
        firms_by_province_cluster$cluster[firms_by_province_cluster$cluster != "Total"])))), "Total")))

firms_by_province_cluster_dot_plot <- ggplot(firms_by_province_cluster, 
    aes(x = cluster, y = num_firms, color = province, shape = is_total)) +
    theme_bw() +
    geom_point(size = 3, alpha = 0.8) +
    facet_wrap(~ province, nrow = 1, scales = "free_x") +
    scale_color_viridis_d(end = 0.9) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    labs(
        x = "Firm Cluster",
        y = "Number of Firms",
        color = "Province"
    ) +
    theme(
        legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(size = 8)
    )

ggsave(
    file.path(FIGURES_PATH, 
        str_glue("firms_by_province_cluster_dot_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_k={CHOSEN_K}.pdf")),
    firms_by_province_cluster_dot_plot,
    width = PAPER_FIG_WIDTH,
    height = PAPER_FIG_HEIGHT
)

print(str_glue("Firms by province and cluster dot plot saved to {FIGURES_PATH}"))

# ============================================================================
# Output single-statistic result snippets to text files
# ============================================================================

# Number of employment spells for largest super cohort panel
write_result_snippet(format_num(largest_super_cohort_stats$num_observations), "num_obs.txt")

# Number of workers for largest super cohort panel
write_result_snippet(format_num(largest_super_cohort_stats$num_workers), "num_workers.txt")

# Number of firms for largest super cohort panel
write_result_snippet(format_num(largest_super_cohort_stats$num_firms), "num_firms.txt")

# Number of outcomes in the CHOSEN_K largest super cohort panel
write_result_snippet(format_num(length(cohorts_for_chosen_k$outcome_ids)), "num_outcomes.txt")

# Number of cohorts in the CHOSEN_K largest super cohort panel
write_result_snippet(format_num(length(cohorts_for_chosen_k$cohort_sizes)), "num_cohorts.txt")

# CHOSEN_K
write_result_snippet(CHOSEN_K, "chosen_k.txt")

# MIN_COHORT_SIZE
write_result_snippet(format_num(MIN_COHORT_SIZE), "min_cohort_size.txt")

# FIRST_YEAR
write_result_snippet(FIRST_YEAR, "first_year.txt")

# LAST_YEAR
write_result_snippet(LAST_YEAR, "last_year.txt")

# YEAR_CLUSTER_SIZE
write_result_snippet(YEAR_CLUSTER_SIZE, "year_cluster_size.txt")

print(str_glue("Result snippets saved to {SNIPPETS_PATH}"))

print("All figures and tables saved successfully.")