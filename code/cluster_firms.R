devtools::load_all("../../../code/apm/r")

source("code/vwh_data_helpers.R")

library(tidyverse)

DATA_PATH <- Sys.getenv("DATA_PATH")
PROCESSED_DATA_PATH = file.path(DATA_PATH, "processed_data")
if (!dir.exists(PROCESSED_DATA_PATH)) {
  dir.create(PROCESSED_DATA_PATH, recursive = TRUE)
}

FIGURES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "figures")
if (!dir.exists(FIGURES_PATH)) {
  dir.create(FIGURES_PATH, recursive = TRUE)
}

earnings <- read_raw_earnings_data()
firms <- read_raw_firm_data()

FIRST_YEAR <- 1998
LAST_YEAR <- 2001
MIN_K <- 2
MAX_K <- 10
MAX_RANK <- 2
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- 50

avg_weekly_earn_by_firm_year <- filter_and_join_match_data(earnings, firms, FIRST_YEAR, LAST_YEAR, VENETO_PROVINCES) |>
    group_years("year", YEAR_CLUSTER_SIZE) |>
    constr_avg_weekly_wages_by_group(c("year", "province", "firm")) |>
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
        raw_firm_clusters_for_prov <- comp_outcome_clusterings(prov_panel, 100, MIN_K, MAX_K, n_inits = 2, num_threads = 10)
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

write_parquet(panels_with_clustered_outcomes, file.path(PROCESSED_DATA_PATH, 
    str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_max_rank={MAX_RANK}.parquet")))

largest_super_cohort_sizes_across_clusterings <- panel_with_clustered_outcomes |>
    group_by(k) |>
    group_modify(\(panel_with_k_clustered_outcomes, keys) {
        print(paste("# clusters:", keys$k[1]))

        num_units <- length(unique(panel_with_k_clustered_outcomes$worker))

        clustered_panel_cohorts <- construct_cohorts_from_panel(
            panel_with_k_clustered_outcomes,
            "worker", "outcome", "avg_weekly_earnings",
            model_rank = MAX_RANK,
            min_cohort_size = MIN_COHORT_SIZE,
        )

        ooi <- clustered_panel_cohorts$observed_outcome_indices
        cohort_sizes <- clustered_panel_cohorts$cohort_sizes

        id_summary <- summarize_identification(ooi, cohort_sizes, MAX_RANK)

        largest_super_cohort_df <- NULL
        for (iter in 1:id_summary$num_o3_iterations) {
            iter_id_summary <- summarize_identification(ooi, cohort_sizes, MAX_RANK, iter = iter)
            
            largest_super_cohort_df <- bind_rows(largest_super_cohort_df, data.frame(
                o3_iter = iter - 1,
                num_units = num_units,
                num_units_largest_super_cohort = iter_id_summary$largest_super_cohort_size,
                largest_super_cohort_share = iter_id_summary$largest_super_cohort_size/num_units
            ))
        }

        largest_super_cohort_df
    })

# Create full grid of all unique k and o3_iter
k_vals <- sort(unique(largest_super_cohort_sizes_across_clusterings$k))
o3_iters <- sort(unique(largest_super_cohort_sizes_across_clusterings$o3_iter))
all_combos <- expand.grid(k = k_vals, o3_iter = o3_iters, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# Join to ensure all combos exist
largest_super_cohort_sizes_across_clusterings <- all_combos |>
    left_join(largest_super_cohort_sizes_across_clusterings, by = c("k", "o3_iter")) |>
    group_by(k) |>
    # For missing combos, fill with entry for largest o3_iter with that k
    mutate(
        largest_super_cohort_share = if_else(
            is.na(largest_super_cohort_share),
            largest_super_cohort_share[which.max(o3_iter[!is.na(largest_super_cohort_share)])],
            largest_super_cohort_share
        ),
        num_units_largest_super_cohort = if_else(
            is.na(num_units_largest_super_cohort),
            num_units_largest_super_cohort[which.max(o3_iter[!is.na(num_units_largest_super_cohort)])],
            num_units_largest_super_cohort
        ),
        num_units = ifelse(
            is.na(num_units),
            num_units[which.max(o3_iter[!is.na(num_units)])],
            num_units
        )
    ) |>
    ungroup()

firm_clustering_plot <- (ggplot(
    largest_super_cohort_sizes_across_clusterings |> 
        mutate(o3_iter = as.factor(o3_iter)))
    + theme_bw()
    + geom_line(aes(x = k, y = largest_super_cohort_share, group = o3_iter, color = o3_iter))
    + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    + scale_color_viridis_d(limits = as.factor(sort(unique(largest_super_cohort_sizes_across_clusterings$o3_iter))), end = 0.9)
    + labs(x = "Number of Firm Clusters per Province", y = "% of Units in Largest Super Cohort", color = "O³\nAlgorithm\nIteration:")
)
ggsave(file.path(FIGURES_PATH, "firm_clustering_plot.pdf"), firm_clustering_plot, width = 8, height = 3.25)