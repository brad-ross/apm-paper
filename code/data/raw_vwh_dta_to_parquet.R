library(tidyverse)
library(haven)
library(arrow)

RAW_DATA_PATH = file.path(Sys.getenv("DATA_PATH"), "raw_data")

workers_data = read_dta(file.path(RAW_DATA_PATH, "anagr.dta"))
write_parquet(workers_data, file.path(RAW_DATA_PATH, "workers.parquet"))

# These legacy firms/earnings files contain bytes undefined in Haven's
# Windows-1252 default. Latin-1 preserves them without dropping/replacing text.
firms_data = read_dta(file.path(RAW_DATA_PATH, "azien.dta"), encoding = "latin1")
write_parquet(firms_data, file.path(RAW_DATA_PATH, "firms.parquet"))

earnings_data = read_dta(file.path(RAW_DATA_PATH, "contr.dta"), encoding = "latin1")
write_parquet(earnings_data, file.path(RAW_DATA_PATH, "earnings.parquet"))
