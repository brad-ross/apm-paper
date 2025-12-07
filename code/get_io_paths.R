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

RESULTS_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "empirical_results")
if (!dir.exists(RESULTS_PATH)) {
  dir.create(RESULTS_PATH, recursive = TRUE)
}

SUMMARY_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "result_summaries")
if (!dir.exists(SUMMARY_PATH)) {
  dir.create(SUMMARY_PATH, recursive = TRUE)
}

FIGURES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "result_summaries", "figures")
if (!dir.exists(FIGURES_PATH)) {
  dir.create(FIGURES_PATH, recursive = TRUE)
}

TABLES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "result_summaries", "tables")
if (!dir.exists(TABLES_PATH)) {
  dir.create(TABLES_PATH, recursive = TRUE)
}

SNIPPETS_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "result_summaries", "snippets")
if (!dir.exists(SNIPPETS_PATH)) {
  dir.create(SNIPPETS_PATH, recursive = TRUE)
}