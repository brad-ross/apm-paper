devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")

library(tidyverse)

earnings <- read_raw_earnings_data()
firms <- read_raw_firm_data()

FIRST_YEAR <- as.integer(Sys.getenv("FIRST_YEAR"))
LAST_YEAR <- as.integer(Sys.getenv("LAST_YEAR"))
MIN_K <- as.integer(Sys.getenv("MIN_K"))
MAX_K <- as.integer(Sys.getenv("MAX_K"))
MAX_RANK <- as.integer(Sys.getenv("MAX_RANK"))
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- as.integer(Sys.getenv("MIN_COHORT_SIZE"))

avg_weekly_earn_by_firm_year <- filter_and_join_match_data(earnings, firms, FIRST_YEAR, LAST_YEAR, VENETO_PROVINCES) |>
    group_years("year", YEAR_CLUSTER_SIZE) |>
    constr_avg_weekly_wages_by_group(c("year", "province", "firm")) |>
    filter(avg_weekly_earnings > 0) |>
    filter_to_always_present_workers_and_firms(collect = TRUE)

firm_clusters = NULL
for (prov in unique(avg_weekly_earn_by_firm_year$province)) {
    data_construction_time <- system.time({
        avg_weekly_earn_by_firm_prov_first_year <- avg_weekly_earn_by_firm_year |>
            filter(province == prov, year == FIRST_YEAR)
        prov_panel <- apm::UnbalancedPanel$new(avg_weekly_earn_by_firm_prov_first_year, "worker", "firm", "avg_weekly_earnings")
    })
    print(paste("Panel creation time for", prov, ":", data_construction_time["elapsed"], "seconds, with", length(prov_panel$get_outcome_ids()), "workers and", length(prov_panel$get_unit_ids()), "firms"))

    clustering_time <- system.time({
        raw_firm_clusters_for_prov <- comp_outcome_clusterings(prov_panel, 100, MIN_K, MAX_K, n_inits = 10, num_threads = Sys.getenv("APM_MAX_THREADS"))
    })
    print(paste("Clustering time for", prov, ":", clustering_time["elapsed"], "seconds"))
    
    firm_clusters_for_prov = NULL
    for (k in MIN_K:MAX_K) {
        k_idx <- k - MIN_K + 1
        firm_clusters_for_prov = bind_rows(firm_clusters_for_prov, data.frame(
            province = prov,
            firm = prov_panel$get_outcome_ids(),
            k = k,
            cluster = setNames(raw_firm_clusters_for_prov[, k_idx], NULL)
        ))
    }

    firm_clusters = bind_rows(firm_clusters, firm_clusters_for_prov)
}

panels_with_clustered_outcomes <- avg_weekly_earn_by_firm_year |>
    .to_arrow_dataset() |>
    inner_join(firm_clusters, by = c("province", "firm")) |>
    constr_avg_weekly_wages_by_group(c("year", "province", "k", "cluster")) |>
    mutate(outcome = paste(year, province, cluster, sep = ":")) |>
    collect()

workers_in_panels_with_T_c_gte_max_rank <- panels_with_clustered_outcomes |>
    .to_arrow_dataset() |>
    select(worker, k, outcome) |>
    distinct() |>
    group_by(worker, k) |>
    summarize(T_c = n()) |>
    filter(T_c >= MAX_RANK) |>
    select(worker, k)

panels_with_clustered_outcomes <- panels_with_clustered_outcomes |>
    .to_arrow_dataset() |>
    inner_join(workers_in_panels_with_T_c_gte_max_rank, by = c("worker", "k")) |>
    collect()

write_parquet(panels_with_clustered_outcomes, file.path(PROCESSED_DATA_PATH, 
    str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.parquet")))

panels_with_clustered_outcomes <- read_parquet(file.path(PROCESSED_DATA_PATH, 
    str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}.parquet"))) |>
    collect()

# Code to plot distributions of cohort sizes and number of cohorts across clusterings:

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
    # + scale_color_viridis_d(end = 0.9)
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

# Code to plot the largest super cohort size across clusterings for a given rank:

comp_largest_super_cohort_sizes_across_clusterings <- function(panels_with_clustered_outcomes, rank) {
    panels_with_clustered_outcomes |>
    group_by(k) |>
    group_modify(\(panel_with_k_clustered_outcomes, keys) {
        print(paste("# clusters:", keys$k[1]))

        num_units <- length(unique(panel_with_k_clustered_outcomes$worker))

        clustered_panel_cohorts <- construct_cohorts_from_panel(
            panel_with_k_clustered_outcomes,
            "worker", "outcome", "avg_weekly_earnings",
            model_rank = rank,
            min_cohort_size = MIN_COHORT_SIZE,
        )

        ooi <- clustered_panel_cohorts$observed_outcome_indices
        cohort_sizes <- clustered_panel_cohorts$cohort_sizes

        id_summary <- summarize_identification(ooi, cohort_sizes, rank)

        largest_super_cohort_df <- NULL
        for (iter in 1:id_summary$num_o3_iterations) {
            iter_id_summary <- summarize_identification(ooi, cohort_sizes, rank, iter = iter)
            
            largest_super_cohort_df <- bind_rows(largest_super_cohort_df, data.frame(
                o3_iter = iter - 1,
                num_units = num_units,
                num_units_largest_super_cohort = iter_id_summary$largest_super_cohort_size,
                largest_super_cohort_share = iter_id_summary$largest_super_cohort_size/num_units
            ))
        }

        largest_super_cohort_df
    })
}

largest_connected_component_sizes_across_clusterings <- comp_largest_super_cohort_sizes_across_clusterings(panels_with_clustered_outcomes, 1) |>
    filter(o3_iter == 1) |>
    mutate(o3_iter = "r=1 final iteration")

for (rank in 2:MAX_RANK) {
    print(paste("Rank:", rank))

    largest_super_cohort_sizes_across_clusterings <- comp_largest_super_cohort_sizes_across_clusterings(panels_with_clustered_outcomes, rank)

    largest_super_cohort_sizes_across_clusterings_for_plot <- largest_super_cohort_sizes_across_clusterings |> 
        mutate(o3_iter = as.character(o3_iter))

    firm_clustering_plot <- (ggplot()
        + theme_bw()
        + geom_line(data = largest_super_cohort_sizes_across_clusterings |>
            group_by(k) |>
            summarize(largest_super_cohort_share = largest_super_cohort_share[which.max(o3_iter)]), 
            aes(x = k, y = largest_super_cohort_share, linetype="Final Largest\nSuper Cohort"))
        + geom_point(data = largest_super_cohort_sizes_across_clusterings |> 
            filter(o3_iter > 0) |>
            mutate(o3_iter = as.factor(o3_iter)),
            aes(x = k, y = largest_super_cohort_share, color = o3_iter))
        + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        + scale_color_viridis_d(end = 0.9)
        + geom_line(data = largest_connected_component_sizes_across_clusterings, 
            aes(x = k, y = largest_super_cohort_share, linetype="r=1 Largest\nConnected\nComponent"))
        + geom_line(data = largest_super_cohort_sizes_across_clusterings |> filter(o3_iter == 0), 
            aes(x = k, y = largest_super_cohort_share, linetype="Initial Largest\nCohort"))
        + scale_linetype_manual(values = c("r=1 Largest\nConnected\nComponent" = "dashed", "Initial Largest\nCohort" = "dotted", "Final Largest\nSuper Cohort" = "solid"))
        + labs(
            x = "Number of Firm Clusters per Province", 
            y = "% of Workers in Largest Super Cohort", 
            color = "O³ Algorithm\nIteration:",
            linetype = NULL)
        + guides(
            color = guide_legend(order = 2),
            linetype = guide_legend(
                order = 1,
                keyheight = grid::unit(1.5, "lines"),
                label.theme = element_text(margin = margin(t = 4, b = 4))
            ))
        + theme(legend.spacing.y = grid::unit(0.1, "lines"))
    )
    ggsave(file.path(FIGURES_PATH, str_glue("firm_cluster_largest_super_cohort_share_plot_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_rank={rank}.pdf")), firm_clustering_plot, width = 8, height = 3.25)
}