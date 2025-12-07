DATA_PATH <- Sys.getenv("DATA_PATH")

RAW_DATA_PATH <- file.path(DATA_PATH, "raw_data")
if (!dir.exists(RAW_DATA_PATH)) {
  dir.create(RAW_DATA_PATH, recursive = TRUE)
}

CLEAN_DATA_PATH <- file.path(DATA_PATH, "clean_data")
if (!dir.exists(CLEAN_DATA_PATH)) {
  dir.create(CLEAN_DATA_PATH, recursive = TRUE)
}

PROCESSED_DATA_PATH = file.path(DATA_PATH, "processed_data")
if (!dir.exists(PROCESSED_DATA_PATH)) {
  dir.create(PROCESSED_DATA_PATH, recursive = TRUE)
}

RESULTS_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "simulation_results")
if (!dir.exists(RESULTS_PATH)) {
  dir.create(RESULTS_PATH, recursive = TRUE)
}

FIGURES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "figures")
if (!dir.exists(FIGURES_PATH)) {
  dir.create(FIGURES_PATH, recursive = TRUE)
}

TABLES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "tables")
if (!dir.exists(TABLES_PATH)) {
  dir.create(TABLES_PATH, recursive = TRUE)
}

RESULT_SNIPPETS_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "result_snippets")
if (!dir.exists(RESULT_SNIPPETS_PATH)) {
  dir.create(RESULT_SNIPPETS_PATH, recursive = TRUE)
}