# Processed data

## Restricted VWH datasets

`code/data/construct_clustered_panels.R` writes these files in `DATA_PATH/processed_data/`:

- `firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet`
- `panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank={MAX_RANK}.parquet`

The panel reporting, estimator evaluation, and wage decomposition programs read these VWH datasets. The datasets remain restricted and are excluded from the public archive.

Only `driver_panel_grid.csv` is approved for distribution. This permission does not extend to the VWH datasets in this directory. The public archive therefore includes only this README and the specified CSV file from this directory.

## Public Figure 1a grid

`driver_panel_grid.csv` contains the rectangle coordinates used to plot driver observation patterns in Figure 1a. The underlying data come from Ater, Ross, Shany and Vasserman (2026), *Road Pricing Under Heterogeneity*. See that paper for details about the driver data.

Each row represents a cohort defined by an entry week and a treatment week, containing at least 100 distinct drivers. The grid contains no driver identifiers or individual records. Run `make figure_1a` from the package root to render the supplied grid; the figure is saved to `output/result_summaries/figures/panel_plot_drivers_times.pdf` with the default settings.

## Grid codebook

The CSV has exactly four numerical columns and one row per retained cohort:

| Column | Meaning and units |
| --- | --- |
| `xmin` | Left edge of the cohort's inclusive week interval, in week coordinates. |
| `xmax` | Right edge of the cohort's inclusive week interval, in week coordinates. |
| `ymin` | Lower edge of the cohort's anonymous consecutive driver display rows. |
| `ymax` | Upper edge of the cohort's anonymous consecutive driver display rows. |

Each rectangle covers the cohort’s drivers and the weeks between entry and treatment, including both endpoints. Week centers are integer week IDs. The width assigned to a week is the smallest spacing between retained week centers, or 1 when there is a single center. The width is 1 whenever a retained cohort spans multiple weeks. Thus `xmin = min(first_week, treatment_week) - width/2` and `xmax = max(first_week, treatment_week) + width/2`, including singleton-week cohorts. Driver row centers are the integers starting at 1; their tiles have height 1. Consequently `ymin = first_display_row - 0.5`, `ymax = last_display_row + 0.5`, and `ymax - ymin` is the number of drivers in the cohort (at least 100). Adjacent rectangles share vertical boundaries; the first begins at `ymin = 0.5`.

Rows are ordered by cohort labels: entry and treatment week IDs are separately padded on the left with zeros to at least two characters, then joined by ` - ` and sorted lexicographically. No driver identifier or individual record is exported. Anonymous display-row ordering within a cohort has no effect on the occupied cells or the figure.

## Grid construction

[construct_driver_panel_grid.R](../../code/data/construct_driver_panel_grid.R) constructs the grid from private `RAW_DATA_PATH/riders.csv`, using `rider_id`, `first_week_in_panel_dateweek_id`, and `began_active_week_id`. Records with missing identifiers or week fields are omitted, and duplicate driver/cohort records are counted once. Blank week fields are treated as missing. A driver must belong to only one cohort, and week intervals must be valid. Only cohorts with at least 100 distinct drivers are retained.

Intervals include both endpoints regardless of the ordering of the entry and treatment weeks. Distinct cohorts can therefore have the same horizontal bounds on different display rows. The supplied grid is sufficient to reproduce Figure 1a; its construction from individual records requires access to the private driver data.
