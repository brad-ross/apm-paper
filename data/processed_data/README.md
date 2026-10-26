# Processed VWH records — restricted

`code/construct_clustered_panels.R` writes these files under `DATA_PATH`:

- `firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet`
- `panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank={MAX_RANK}.parquet`

They feed panel reports, evaluation and decomposition, remain restricted, and
are ignored by Git and excluded from the public archive. No derived-data
disclosure permission is inferred. Only this README is included. The public
Figure 1a grid belongs in `clean_data`. See `docs/replication/data-sources.md`.
