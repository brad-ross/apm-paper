# Raw and converted records — restricted

Authorized VWH users place `anagr.dta` (workers), `azien.dta` (firms) and
`contr.dta` (employment/contributions) here under `DATA_PATH`. All three are
required. Conversion writes `workers.parquet`, `firms.parquet` and
`earnings.parquet` here. FRDB access instructions and unconfirmed access costs,
delays and provenance are in `docs/replication/data-sources.md`.

Conversion reads the legacy `azien.dta` and `contr.dta` with explicit Latin-1
decoding: they contain bytes undefined in Haven's default Windows-1252
encoding (including firm-name/address fields in `azien.dta`). This preserves
those bytes as Unicode characters; it does not discard or replace text.
The workers file, `anagr.dta`, uses Haven's default.

Private driver preparation separately reads `riders.csv`, using `rider_id`,
`first_week_in_panel_dateweek_id` and `began_active_week_id`. Driver access and
source arrangements remain unconfirmed. Public rendering never reads this
file or automatically regenerates a missing grid.

Neither source nor converted files are approved for distribution. Only this
README is included from this directory. Anonymous identifiers or aggregation
do not establish disclosure rights.
