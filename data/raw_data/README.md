# Raw and converted records — restricted

Authorized VWH users place `anagr.dta` (workers), `azien.dta` (firms) and
`contr.dta` (employment/contributions) here under `DATA_PATH`. All three are
required. Conversion writes `workers.parquet`, `firms.parquet` and
`earnings.parquet` here. FRDB access instructions and unconfirmed access costs,
delays and provenance are in `docs/replication/data-sources.md`.

Private driver preparation separately reads `riders.csv`, using `rider_id`,
`first_week_in_panel_dateweek_id` and `began_active_week_id`. Driver access and
source arrangements remain unconfirmed. Public rendering never reads this
file or automatically regenerates a missing grid.

Neither source nor converted files are approved for distribution. Only this
README is included from this directory. Anonymous identifiers or aggregation
do not establish disclosure rights.
