#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) stop('Run Rscript setup.R without arguments')
source('config/paper.R')
lib <- Sys.getenv('REPLICATION_LIBRARY')
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
lib <- normalizePath(lib, mustWork=TRUE)
Sys.setenv(REPLICATION_LIBRARY=lib, RENV_PATHS_ROOT=file.path(getwd(),'.replication','renv'))
lock_dir <- paste0(lib,'.setup-lock')
if (!dir.create(lock_dir, showWarnings=FALSE)) stop('Another setup owns ',lock_dir,'. Retry after it finishes; remove only a confirmed stale lock.')
# Top-level scripts have no on.exit scope; use a function for reliable cleanup.
setup <- function() {
  on.exit(unlink(lock_dir,recursive=TRUE),add=TRUE)
  previous <- .libPaths()
  .libPaths(c(lib, previous))
  if (!dir.exists(file.path(lib,'renv'))) {
    existing <- find.package('renv',quiet=TRUE)
    if (length(existing) && as.character(packageVersion('renv')) == '1.3.1') {
      if (!file.copy(existing,lib,recursive=TRUE)) stop('Cannot copy renv into project library')
    } else {
      url <- 'https://cran.r-project.org/src/contrib/Archive/renv/renv_1.3.1.tar.gz'
      archive <- tempfile(fileext='.tar.gz')
      # Current CRAN location first; the archive remains valid after a newer release.
      ok <- tryCatch({download.file(sub('/Archive/renv','',url),archive,quiet=TRUE);TRUE},error=function(e) FALSE)
      if (!ok) download.file(url,archive,quiet=TRUE)
      install.packages(archive,lib=lib,repos=NULL,type='source')
    }
  }
  if (!requireNamespace('renv',quietly=TRUE)) stop('renv bootstrap failed')
  lockfile <- 'renv.lock'
  lock <- renv::lockfile_read(lockfile)
  current_r <- paste(R.version$major,R.version$minor,sep='.')
  if (current_r != lock$R$Version) stop('Use R ',lock$R$Version,' to match this replication environment (found ',current_r,').')
  matches <- function(name, record) {
    path <- file.path(lib,name,'DESCRIPTION')
    file.exists(path) && identical(unname(read.dcf(path,fields='Version')[1,1]),record$Version)
  }
  needed <- names(Filter(function(x) !matches(x$Package,x),lock$Packages))
  if(!requireNamespace('jsonlite',quietly=TRUE)) {
    renv::restore(project=getwd(),library=lib,lockfile=lockfile,packages='jsonlite',prompt=FALSE)
  }
  if(!file.exists('vendor/apm_source.json')) stop('Missing bundled APM source manifest')
  source_meta <- jsonlite::read_json('vendor/apm_source.json',simplifyVector=TRUE)
  asset <- source_meta$source_tarball
  if(is.null(asset)) stop('APM manifest must contain source_tarball')
  install_apm <- !file.exists(file.path(lib,'apm','replication-source.json'))
  if(!install_apm) {
    installed <- jsonlite::read_json(file.path(lib,'apm','replication-source.json'),simplifyVector=TRUE)
    install_apm <- !identical(installed$asset_sha256,source_meta$source_sha256)
  }
  for(name in needed) {
    for(base in previous) {
      src <- file.path(base,name)
      desc <- file.path(src,'DESCRIPTION')
      if(file.exists(desc) && identical(unname(read.dcf(desc,fields='Version')[1,1]),lock$Packages[[name]]$Version)) {
        unlink(file.path(lib,name),recursive=TRUE)
        if(!file.copy(src,lib,recursive=TRUE)) stop('Could not copy dependency ',name)
        break
      }
    }
  }
  needed <- names(Filter(function(x) !matches(x$Package,x),lock$Packages))
  if(length(needed)) renv::restore(project=getwd(),library=lib,lockfile=lockfile,packages=needed,prompt=FALSE)
  if(!all(vapply(names(lock$Packages),function(n) matches(n,lock$Packages[[n]]),logical(1)))) stop('Dependency restoration incomplete')
  source('code/lib/apm_source.R')
  asset <- ensure_apm_source(source_meta)
  if(install_apm) {
    if(!nzchar(Sys.which('make'))) stop('Installing APM requires Make and a C++17 compiler; see README')
    status <- system2(file.path(R.home('bin'),'R'), c('CMD','INSTALL',paste0('--library=',shQuote(lib)),shQuote(asset)))
    if(status != 0L) stop('APM installation failed; previous installation was not certified')
    if(!file.exists(file.path(lib,'apm','DESCRIPTION'))) stop('APM installation failed')
    jsonlite::write_json(list(asset_sha256=source_meta$source_sha256,commit=source_meta$upstream_commit),file.path(lib,'apm','replication-source.json'),auto_unbox=TRUE,pretty=TRUE)
  }
  if(!requireNamespace('apm',quietly=TRUE,lib.loc=lib)) stop('Installed APM does not load')
  cat('Setup ready: environment in',lib,'\n')
}
setup()
