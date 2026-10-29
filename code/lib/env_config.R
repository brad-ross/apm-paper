source("config/paper.R")

# Shared environment configuration for replication scripts
#
# This file centralizes environment variable parsing used across multiple scripts.
# Source this file to get access to common configuration variables.

# Shared random seed for clustering and bootstrap draws
APM_SEED <- suppressWarnings(as.numeric(Sys.getenv("APM_SEED", unset = "9001")))
if (is.na(APM_SEED) || !is.finite(APM_SEED) ||
    APM_SEED < 0 || APM_SEED > .Machine$integer.max ||
    APM_SEED != floor(APM_SEED)) {
    stop("APM_SEED must be a nonnegative integer no greater than ", .Machine$integer.max)
}
APM_SEED <- as.integer(APM_SEED)

# Year range configuration
FIRST_YEAR <- as.integer(Sys.getenv("FIRST_YEAR"))
LAST_YEAR <- as.integer(Sys.getenv("LAST_YEAR"))

# Model configuration
MAX_RANK <- as.integer(Sys.getenv("MAX_RANK"))
YEAR_CLUSTER_SIZE <- as.integer(ceiling((LAST_YEAR - FIRST_YEAR) / MAX_RANK))
MIN_COHORT_SIZE <- as.integer(Sys.getenv("MIN_COHORT_SIZE"))

# Clustering configuration
MIN_K <- as.integer(Sys.getenv("MIN_K"))
MAX_K <- as.integer(Sys.getenv("MAX_K"))
CHOSEN_K <- as.integer(Sys.getenv("CHOSEN_K"))

# Threading configuration
num_threads <- as.integer(Sys.getenv("APM_MAX_THREADS"))
if (is.na(num_threads)) {
    num_threads <- NULL
}

# Figure dimensions for paper
PAPER_FIG_HEIGHT <- 3.5
PAPER_FIG_WIDTH <- 8
PAPER_FIG_HALF_WIDTH <- PAPER_FIG_WIDTH / 2

for (name in c("FIRST_YEAR", "LAST_YEAR", "MAX_RANK", "MIN_COHORT_SIZE", "MIN_K", "MAX_K", "CHOSEN_K", "num_threads")) {
  value <- get(name)
  if (length(value) != 1L || is.na(value) || value < 1L) stop("Invalid positive integer setting: ", name)
}
if (LAST_YEAR < FIRST_YEAR || MIN_K > MAX_K || CHOSEN_K < MIN_K || CHOSEN_K > MAX_K) stop("Inconsistent paper configuration")
