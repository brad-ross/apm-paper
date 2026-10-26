#!/usr/bin/env Rscript
build_archive <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) stop('Run from the package root: Rscript tools/build_replication_archive.R (no arguments)')
  root <- normalizePath(getwd(), mustWork = TRUE)
  lib <- Sys.getenv('REPLICATION_LIBRARY', file.path(root, '.replication', 'library'))
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  for (package in c('digest', 'jsonlite')) {
    if (!requireNamespace(package, quietly = TRUE)) stop('Missing ', package, '; run Rscript setup.R --scope public first')
  }
  sha <- function(path) digest::digest(file = path, algo = 'sha256')
  # Check every component before normalization: normalizePath alone hides symlinks.
  regular <- function(path) {
    if (!nzchar(path) || grepl('[\\\\[:space:]:]', path) || startsWith(path, '/') ||
        any(strsplit(path, '/', fixed = TRUE)[[1]] %in% c('', '.', '..'))) stop('Unsafe path: ', path)
    parts <- strsplit(path, '/', fixed = TRUE)[[1]]
    current <- root
    for (part in parts) {
      current <- file.path(current, part)
      link <- Sys.readlink(current)
      if (!is.na(link) && nzchar(link)) stop('Symlink is forbidden: ', path)
    }
    if (!file.exists(current) || !file_test('-f', current) || isTRUE(file.info(current)$isdir)) stop('Missing required regular file: ', path)
    canonical <- normalizePath(current, mustWork = TRUE)
    if (!startsWith(canonical, paste0(root, '/'))) stop('Path escapes project: ', path)
    canonical
  }
  manifest_path <- regular('release/public-files.txt')
  files <- trimws(readLines(manifest_path, warn = FALSE))
  files <- files[nzchar(files) & !startsWith(files, '#')]
  if (!length(files) || anyDuplicated(files)) stop('Allowlist must contain unique explicit file paths')
  required <- c('README.md', 'README.pdf', 'LICENSE', 'Makefile', 'setup.R', 'run.R',
                'renv.lock', 'config/public.lock', 'config/paper.R',
                'software/apm-source.json', 'software/README.md',
                'vendor/apm_0.1.0.tar.gz', 'data/clean_data/driver_panel_grid.csv',
                'data/clean_data/README.md', 'docs/replication/verification.md',
                'docs/replication/data-sources.md', 'docs/replication/exhibits.csv',
                'release/public-files.txt', 'tools/build_replication_archive.R')
  missing <- setdiff(required, files)
  if (length(missing)) stop('Required entries missing from allowlist: ', paste(missing, collapse = ', '))
  data_docs <- c('data/README.md', 'data/raw_data/README.md',
                 'data/processed_data/README.md', 'data/clean_data/README.md')
  for (path in files) {
    regular(path)
    if (any(startsWith(strsplit(path, '/', fixed = TRUE)[[1]], '.'))) stop('Hidden/private path forbidden: ', path)
    if (startsWith(path, 'data/') && !path %in% c(data_docs, 'data/clean_data/driver_panel_grid.csv')) stop('Unapproved data path: ', path)
    if (grepl('\\.csv$', path, ignore.case = TRUE) && startsWith(path, 'data/') && path != 'data/clean_data/driver_panel_grid.csv') stop('Unapproved data CSV: ', path)
    allowed <- path %in% c('README.md', 'README.pdf', 'LICENSE', 'Makefile', 'setup.R', 'run.R', 'renv.lock') ||
      path %in% c(data_docs, 'data/clean_data/driver_panel_grid.csv', 'vendor/apm_0.1.0.tar.gz') ||
      grepl('^(code/[^/]+\\.(R|sh|slurm)|code/jobs/[^/]+\\.sh|config/(paper\\.R|public\\.lock)|software/(README\\.md|apm-source\\.json)|docs/replication/[^/]+\\.(md|csv)|release/public-files\\.txt|tools/build_replication_archive\\.R)$', path)
    if (!allowed) stop('Path outside approved package namespaces: ', path)
  }
  meta <- jsonlite::read_json(regular('software/apm-source.json'), simplifyVector = TRUE)
  if (!identical(meta$source_tarball, 'vendor/apm_0.1.0.tar.gz') ||
      !identical(sha(regular(meta$source_tarball)), meta$source_sha256)) stop('Frozen APM source asset identity/checksum mismatch')
  files <- sort(files, method = 'radix')
  hashes <- vapply(files, function(path) sha(regular(path)), character(1))
  dist <- file.path(root, 'dist')
  link <- Sys.readlink(dist)
  if (!is.na(link) && nzchar(link)) stop('dist must not be a symlink')
  if (!dir.exists(dist) && !dir.create(dist)) stop('Cannot create dist directory')
  if (!startsWith(normalizePath(dist, mustWork = TRUE), paste0(root, '/'))) stop('dist escapes project')
  build <- tempfile('build-', tmpdir = dist)
  if (!dir.create(build)) stop('Cannot create staging directory')
  package <- file.path(build, 'replication')
  dir.create(package)
  for (i in seq_along(files)) {
    source <- regular(files[[i]])
    destination <- file.path(package, files[[i]])
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, destination, copy.mode = TRUE) || !identical(sha(destination), hashes[[i]])) stop('Staging copy/hash failed: ', files[[i]])
  }
  writeLines(paste(hashes, files, sep = '  '), file.path(package, 'SHA256SUMS'))
  writeLines(c(files, 'SHA256SUMS'), file.path(build, 'files.txt'))
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(build)
  if (!nzchar(Sys.which('zip'))) stop('Archive assembly requires the zip utility; install it and retry. Replication itself does not need zip.')
  zip_archive <- file.path(build, 'replication.zip')
  status <- utils::zip(zip_archive, files = paste0('replication/', c(files, 'SHA256SUMS')), flags = '-q')
  if (status != 0L || !file.exists(zip_archive)) stop('ZIP archive creation failed')
  zip_listing <- utils::unzip(zip_archive, list = TRUE)$Name
  if (!setequal(zip_listing, paste0('replication/', c(files, 'SHA256SUMS'))) || anyDuplicated(zip_listing)) stop('ZIP listing differs from staged allowlist')
  writeLines(paste(sha(zip_archive), 'replication.zip', sep = '  '), file.path(build, 'zip.sha256'))
  archive <- file.path(build, 'replication.tar.gz')
  utils::tar(archive, files = 'replication', compression = 'gzip', tar = 'internal')
  listing <- utils::untar(archive, list = TRUE, tar = 'internal')
  expected <- paste0('replication/', c(files, 'SHA256SUMS'))
  archived_files <- listing[!endsWith(listing, '/')]
  if (!setequal(archived_files, expected) || anyDuplicated(archived_files)) stop('Archive listing differs from staged allowlist')
  writeLines(listing, file.path(build, 'archive-listing.txt'))
  writeLines(paste(sha(archive), 'replication.tar.gz', sep = '  '), file.path(build, 'archive.sha256'))
  file.copy(file.path(package, 'SHA256SUMS'), file.path(build, 'SHA256SUMS'))
  cat('Local ZIP archive:', zip_archive, '\nLocal tar archive:', archive, '\nListing:', file.path(build, 'archive-listing.txt'),
      '\nHashes:', file.path(build, 'SHA256SUMS'), '\nArchive checksum:',
      file.path(build, 'archive.sha256'), '\nNo files were published.\n')
}
build_archive()
