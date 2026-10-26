# Run from the repository root: Rscript tests/test_seed_config.R
# This dry run requires no packages or empirical data.
local({
    original_seed <- Sys.getenv("APM_SEED", unset = NA_character_)
    on.exit(if (is.na(original_seed)) Sys.unsetenv("APM_SEED") else
        Sys.setenv(APM_SEED = original_seed))
    original_threads <- Sys.getenv("APM_MAX_THREADS", unset = NA_character_)
    on.exit(if (is.na(original_threads)) Sys.unsetenv("APM_MAX_THREADS") else
        Sys.setenv(APM_MAX_THREADS = original_threads), add = TRUE)

    read_config <- function(seed = NULL) {
        if (is.null(seed)) Sys.unsetenv("APM_SEED") else
            Sys.setenv(APM_SEED = seed)
        config <- new.env()
        sys.source("code/env_config.R", envir = config)
        config
    }

    stopifnot(identical(read_config()$APM_SEED, 9001L))
    for (seed in c("0", "1234", "2147483647")) {
        stopifnot(identical(read_config(seed)$APM_SEED, as.integer(seed)))
    }
    for (seed in c("", "abc", "NA", "NaN", "Inf", "-1", "1.5", "2147483648")) {
        error <- tryCatch(read_config(seed), error = identity)
        stopifnot(inherits(error, "error"), grepl("APM_SEED", conditionMessage(error)))
    }

    # Stub the module command; inspect an exported seed in a child process.
    for (setting in c("unset APM_SEED", "export APM_SEED=1234")) {
        command <- paste(setting,
            "ml() { :; }; source code/init_server.sh; printenv APM_SEED",
            sep = "; ")
        actual <- system2("bash", c("-c", shQuote(command)), stdout = TRUE)
        expected <- if (setting == "unset APM_SEED") "9001" else "1234"
        stopifnot(identical(actual, expected))
    }

    Sys.unsetenv("APM_MAX_THREADS")
    config <- read_config("1234")
    stopifnot(is.null(config$num_threads))
    # Load only the helper definition, avoiding its unrelated arrow dependency.
    helper <- Filter(function(expr) is.call(expr) &&
        identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("get_default_bootstrap")),
        as.list(parse("code/estimation_helpers.R")))
    stopifnot(length(helper) == 1L)
    eval(helper[[1]], envir = config)
    config$get_weighted_bootstrap_draws <- function(n_units, n_draws, type, seed) {
        list(n_units = n_units, n_draws = n_draws, type = type, seed = seed)
    }
    panel <- list(get_num_units = function() 7L)
    stopifnot(identical(config$get_default_bootstrap(panel),
        list(n_units = 7L, n_draws = 1000, type = "bayesian", seed = 1234L)))
    stopifnot(identical(config$get_default_bootstrap(panel, seed = 42L)$seed, 42L))

    # Exercise the clustering call without running the empirical script.
    find_clustering_calls <- function(expr) {
        if (!is.call(expr) && !is.expression(expr)) return(list())
        if (is.call(expr) && identical(expr[[1]], as.name("comp_outcome_clusterings"))) {
            return(list(expr))
        }
        unlist(lapply(as.list(expr), find_clustering_calls), recursive = FALSE)
    }
    calls <- find_clustering_calls(parse("code/construct_clustered_panels.R"))
    stopifnot(length(calls) == 1L)
    config$prov_panel <- panel
    config$MIN_K <- 1L
    config$MAX_K <- 3L
    config$comp_outcome_clusterings <- function(..., seed, num_threads) {
        list(seed = seed, num_threads = num_threads)
    }
    stopifnot(identical(eval(calls[[1]], envir = config),
        list(seed = 1234L, num_threads = NULL)))

    cat("Seed configuration dry run passed.\n")
})
