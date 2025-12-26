devtools::load_all("../apm-package/r")

source("code/vwh_data_helpers.R")
source("code/get_io_paths.R")
source("code/env_config.R")

library(tidyverse)

earnings <- read_raw_earnings_data()
firms <- read_raw_firm_data()

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

write_parquet(firm_clusters, file.path(PROCESSED_DATA_PATH, 
    str_glue("firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet")))

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
    str_glue("panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank={MAX_RANK}.parquet")))

print(str_glue("Clustered panels saved to {PROCESSED_DATA_PATH}"))