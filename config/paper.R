# Defaults for the paper. Environment overrides are recorded in run manifests.
paper_defaults <- c(FIRST_YEAR="1998", LAST_YEAR="2001", MAX_RANK="2",
                    MIN_COHORT_SIZE="75", MIN_K="1", MAX_K="25", CHOSEN_K="3",
                    APM_SEED="9001", APM_MAX_THREADS="1")
for (name in names(paper_defaults)) {
  if (!nzchar(Sys.getenv(name))) do.call(Sys.setenv, setNames(list(paper_defaults[[name]]), name))
}
if (!nzchar(Sys.getenv("DATA_PATH"))) Sys.setenv(DATA_PATH=file.path(getwd(), "data"))
if (!nzchar(Sys.getenv("OUTPUT_PATH"))) Sys.setenv(OUTPUT_PATH=file.path(getwd(), "output"))
if (!nzchar(Sys.getenv("REPLICATION_LIBRARY"))) {
  Sys.setenv(REPLICATION_LIBRARY=file.path(getwd(), ".replication", "library"))
}
