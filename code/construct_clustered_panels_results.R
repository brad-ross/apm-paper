# Generate figures for clustered panel analysis
#
# This script reads the parquet output from construct_clustered_panels.R
# and generates figures for the paper.

devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")
source("code/estimation_helpers.R")

library(tidyverse)
library(arrow)

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

ggsave(file.path(FIGURES_PATH, str_glue("cohort_size_dist_and_num_cohorts_by_k_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.pdf")), cohort_size_dist_and_num_cohorts_by_k_plot, width = 8, height = 3.25)

# ============================================================================
# Code to plot the largest super cohort size across clusterings for a given rank
# ============================================================================

firms_per_firm_cluster <- firm_clusters |>
    group_by(k, province, cluster) |>
    summarize(num_firms = n_distinct(firm), .groups = "drop")

outcomes_per_outcome_cluster_df <- panels_with_clustered_outcomes |>
    .to_arrow_dataset() |>
    select(k, year, province, cluster, outcome) |>
    distinct() |>
    left_join(firms_per_firm_cluster, by = c("k", "province", "cluster")) |>
    collect()

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

        outcomes_per_outcome_cluster <- outcomes_per_outcome_cluster_df |>
            filter(k == curr_k) |>
            mutate(outcome_idx = match(outcome, clustered_panel_cohorts$outcome_ids)) |>
            arrange(outcome_idx) |>
            pull(num_firms)

        print(length(outcomes_per_outcome_cluster))

        num_outcomes <- length(clustered_panel_cohorts$outcome_ids)
        print(num_outcomes)
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

            outcomes_per_outcome_cluster <- outcomes_per_outcome_cluster_df |>
                filter(k == curr_k) |>
                mutate(outcome_idx = match(outcome, clustered_panel_cohorts$outcome_ids)) |>
                arrange(outcome_idx) |>
                pull(num_firms)

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

largest_connected_component_sizes_across_clusterings <- comp_largest_super_cohort_sizes_across_clusterings(panels_with_clustered_outcomes, 1) |>
    filter(o3_iter == 1) |>
    mutate(o3_iter = "r=1 final iteration")

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

    ggplot(df) +
        theme_bw() +
        stat_ecdf(aes(
            x = !!x_expr,
            weight = cohort_size,
            color = as.factor(k)
        )) +
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
        width = 4,
        height = 3.25
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
        width = 4,
        height = 3.25
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
    ggsave(file.path(FIGURES_PATH, str_glue("num_o3_iterations_needed_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), num_o3_iterations_needed_plot, width = 8, height = 3.25)

    # Plot of distribution of share of cohort outcome means identified across cohorts

    block_missingness_identification_across_clusterings <- comp_block_missingness_identification_across_clusterings(
        panels_with_clustered_outcomes,
        rank
    )

    block_missingness_id_plot <- create_block_missingness_ecdf_plot(
        df = block_missingness_identification_across_clusterings,
        x_expr = outcomes_with_rank_overlap / num_outcomes,
        x_label = "Share of Cohort's Outcome Means Identified by Block Missingness"
    )

    ggsave(file.path(FIGURES_PATH, str_glue("block_missingness_id_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), block_missingness_id_plot, width = 8, height = 3.25)

    firm_weighted_block_missingness_id_plot <- create_block_missingness_ecdf_plot(
        df = block_missingness_identification_across_clusterings,
        x_expr = total_outcome_weight_with_rank_overlap / num_firms,
        x_label = "Share of Cohort's Avg. Wages Across Firms Identified by Block Missingness"
    )

    ggsave(file.path(FIGURES_PATH, str_glue("block_missingness_id_firm_weighted_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), firm_weighted_block_missingness_id_plot, width = 8, height = 3.25)
}

print("All figures saved successfully.")