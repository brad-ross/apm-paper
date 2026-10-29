# Restricted raw and converted data

Users with authorized access to VWH should place `anagr.dta` (workers), `azien.dta` (firms) and
`contr.dta` (employment/contributions) here under `DATA_PATH`. All three files are
required. The conversion program writes `workers.parquet`, `firms.parquet` and
`earnings.parquet` here. FRDB access instructions are in the root `README.md`.

The conversion program reads the legacy `azien.dta` and `contr.dta` with explicit Latin-1
decoding: they contain bytes undefined in Haven's default Windows-1252
encoding (including firm-name/address fields in `azien.dta`). Decoding these bytes as Unicode characters preserves the text without
discarding or replacing characters.
The workers file, `anagr.dta`, uses Haven's default.

The separate driver-data preparation program reads `riders.csv`, using `rider_id`,
`first_week_in_panel_dateweek_id` and `began_active_week_id`. See Ater, Ross, Shany and Vasserman (2026), *Road Pricing Under Heterogeneity*,
for details about the driver data. These private records cannot be shared.
The public plotting program reads neither this file nor any other individual
driver records, and it does not reconstruct a missing grid.

Neither the source files nor their converted versions are approved for
distribution. The public archive includes only this README from this directory.
Anonymizing identifiers or aggregating records does not itself confer
permission to distribute the resulting data.
