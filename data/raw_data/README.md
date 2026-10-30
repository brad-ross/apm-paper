# Raw data

This directory contains the source records used by the replication programs and the converted VWH files. Its location is `DATA_PATH/raw_data/`, which defaults to `data/raw_data/` within the package. None of these data files is included in the public archive; only this README is distributed.

## Veneto Workers Histories (VWH)

Obtain the three source files through the application process described in the [main README](../../README.md#veneto-workers-histories-vwh) and place them in this directory. All three are required for the VWH analysis.

| File | Contents/source | Availability |
|---|---|---|
| `anagr.dta` | VWH worker registry | Obtain from provider |
| `azien.dta` | VWH firm registry | Obtain from provider |
| `contr.dta` | VWH employment/contribution records | Obtain from provider |

The `conversion` stage runs [raw_vwh_dta_to_parquet.R](../../code/data/raw_vwh_dta_to_parquet.R) and writes `workers.parquet`, `firms.parquet`, and `earnings.parquet`, respectively, in this directory. The remaining VWH stages read these converted files. Running `make vwh-all` from the package root includes conversion automatically; Stata is not required.

### Character encoding

The conversion program reads `azien.dta` and `contr.dta` using Latin-1 encoding and `anagr.dta` using Haven's default encoding.

## Driver data and Figure 1a

The private source file, `riders.csv`, comes from *Road Pricing Under Heterogeneity* by Ater, Ross, Shany and Vasserman (2026). See that paper for details about the driver data. These records cannot be shared, including with the replication team, and are used only to construct the plotting grid for Figure 1a.

The optional `figure_1a_data_prep` stage runs [construct_driver_panel_grid.R](../../code/data/construct_driver_panel_grid.R). It reads `riders.csv` from this directory, using `rider_id`, `first_week_in_panel_dateweek_id`, and `began_active_week_id`, and writes `driver_panel_grid.csv` to `DATA_PATH/processed_data/`.

Public replication uses the supplied grid, documented in the [processed-data README](../processed_data/README.md). Running `make figure_1a` reads that grid without requiring individual driver records.

## Distribution restrictions

The source and converted records in this directory are excluded from the public archive. VWH files must not be redistributed without permission, and the individual driver records cannot be shared. The only data file approved for distribution with this package is `driver_panel_grid.csv` in the processed-data directory.
