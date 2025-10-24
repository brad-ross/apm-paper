DATA_PATH <- Sys.getenv("DATA_PATH")
PROCESSED_DATA_PATH = file.path(DATA_PATH, "processed_data")
if (!dir.exists(PROCESSED_DATA_PATH)) {
  dir.create(PROCESSED_DATA_PATH, recursive = TRUE)
}

FIGURES_PATH = file.path(Sys.getenv("OUTPUT_PATH"), "figures")
if (!dir.exists(FIGURES_PATH)) {
  dir.create(FIGURES_PATH, recursive = TRUE)
}

