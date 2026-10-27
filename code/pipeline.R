pipeline_stages <- c(conversion='raw_vwh_dta_to_parquet.R', panels='construct_clustered_panels.R',
  `panel-reports`='construct_clustered_panels_results.R', evals='leave_cohort_outcome_out_evals.R',
  `eval-reports`='leave_cohort_outcome_out_evals_results.R', decomposition='match_outcome_decomposition.R',
  `decomp-reports`='match_outcome_decomposition_results.R', `prepare-figure`='construct_driver_panel_grid.R',
  `figure-1a`='driver_panel_results.R')
pipeline_threads <- c(conversion=4L,panels=40L,`panel-reports`=16L,evals=250L,`eval-reports`=4L,decomposition=250L,`decomp-reports`=4L,`prepare-figure`=4L,`figure-1a`=1L)
pipeline_control_root <- file.path(getwd(),'.replication')
pipeline_mode <- 'local'
pipeline_hash_cache <- new.env(parent=emptyenv())
pipeline_package_cache <- new.env(parent=emptyenv())
pipeline_deps <- list(conversion=character(), panels='conversion', `panel-reports`='panels',
  evals='panels', `eval-reports`='evals', decomposition='panels', `decomp-reports`='decomposition',
  `prepare-figure`=character(), `figure-1a`=character())
pipeline_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, version=2)
  unname(tools::md5sum(f))
}
pipeline_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop('Missing required input/output: ', paste(missing,collapse=', '),
    '. Restore supplied public grid or obtain documented authorized inputs; see README.')
  present <- paths[file.exists(paths)]
  if (!length(present)) return(setNames(character(),character()))
  if (any(file.info(present)$isdir)) stop('Expected regular files: ',paste(present,collapse=', '))
  hashes <- vapply(present,function(path) {
    info <- file.info(path)
    key <- paste(path,info$size,as.numeric(info$mtime),sep='|')
    if (!exists(key,pipeline_hash_cache,inherits=FALSE)) assign(key,unname(tools::md5sum(path)),pipeline_hash_cache)
    get(key,pipeline_hash_cache,inherits=FALSE)
  },character(1))
  setNames(unname(hashes),present)
}
pipeline_atomic <- function(object, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  f <- tempfile(tmpdir=dirname(path)); saveRDS(object,f)
  if (!file.rename(f,path)) stop('Cannot atomically save ',path)
}
pipeline_lock <- function(path, wait=0) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  deadline <- Sys.time() + wait
  repeat {
    if (dir.create(path,showWarnings=FALSE)) break
    if (Sys.time() >= deadline) stop('Writer lock exists: ',path,
    '. Check its owner; remove only after confirming the process/job has ended.')
    Sys.sleep(0.1)
  }
  writeLines(c(paste('host',Sys.info()[['nodename']]),paste('pid',Sys.getpid()),
    paste('job',Sys.getenv('SLURM_JOB_ID')),as.character(Sys.time())),file.path(path,'owner'))
  path
}
pipeline_outputs <- function(stage) {
  core <- sprintf('start_year=%s_end_year=%s_year_cluster_size=%s',FIRST_YEAR,LAST_YEAR,YEAR_CLUSTER_SIZE)
  study <- sprintf('%s_min_cohort_size=%s',core,MIN_COHORT_SIZE)
  chosen <- sprintf('%s_k=%s',study,CHOSEN_K)
  raw <- function(x) file.path(RAW_DATA_PATH,x)
  fig <- function(x) file.path(FIGURES_PATH,x)
  snippets <- function(x) file.path(SNIPPETS_PATH,x)
  switch(stage,
    conversion=raw(paste0(c('workers','firms','earnings'),'.parquet')),
    panels=file.path(PROCESSED_DATA_PATH,c(paste0('firm_clusters_',core,'.parquet'),
      sprintf('panels_with_clustered_outcomes_%s_max_rank=%s.parquet',core,MAX_RANK))),
    evals=file.path(RESULTS_PATH,paste0('masked_metrics_',chosen,'.parquet')),
    decomposition=file.path(RESULTS_PATH,'match_outcome_decomp_results.parquet'),
    `prepare-figure`=file.path(CLEAN_DATA_PATH,'driver_panel_grid.csv'),
    `figure-1a`=fig('panel_plot_drivers_times.pdf'),
    `panel-reports`=c(fig(paste0(c('cohort_size_dist_plot_','match_outcome_panel_structure_plot_',
      'block_missingness_id_plot_','block_missingness_id_firm_weighted_plot_'),chosen,'.pdf')),
      fig(unlist(lapply(2:MAX_RANK,function(r) paste0(c('largest_super_cohort_share_workers_plot_',
        'largest_super_cohort_share_firms_plot_','num_o3_iterations_needed_plot_'),study,'_rank=',r,'.pdf')))),
      fig(sprintf('firms_by_province_cluster_dot_plot_%s_k=%s.pdf',core,CHOSEN_K)),
      file.path(TABLES_PATH,paste0('panel_summary_stats_',chosen,'.tex')),
      snippets(paste0(c('num_obs','num_workers','num_firms','num_outcomes','num_cohorts','chosen_k',
        'min_cohort_size','first_year','last_year','year_cluster_size'),'.txt'))),
    `eval-reports`=c(file.path(TABLES_PATH,'leave_cohort_outcome_out_evals_results.tex'),
      snippets(paste0(apply(expand.grid(c('pc_1_no_fes_equal_weights','pc_1_fes_equal_weights',
        'pc_2_no_fes_equal_weights','pc_2_fes_equal_weights'),c('bias','se','rmse'),
        c('geom_mean_value_ratio','share_better','avg_pct_better'),stringsAsFactors=FALSE),1,paste,collapse='_'),'.txt'))),
    `decomp-reports`=c(file.path(TABLES_PATH,'match_outcome_decomp_results.tex'),
      snippets(paste0(apply(expand.grid('fgw_decomp_province',c('apm','twfe','difference'),
        c('estimate','std_error','p_value','ci_lb','ci_ub'),stringsAsFactors=FALSE),1,paste,collapse='_'),'.txt'))),
    stop('Unknown stage: ',stage))
}
pipeline_inputs <- function(stage) {
  switch(stage, conversion=file.path(RAW_DATA_PATH,c('anagr.dta','azien.dta','contr.dta')),
    `prepare-figure`=file.path(RAW_DATA_PATH,'riders.csv'),
    `figure-1a`=file.path(CLEAN_DATA_PATH,'driver_panel_grid.csv'),
    unique(unlist(lapply(pipeline_deps[[stage]],pipeline_outputs))))
}
pipeline_destination_key <- function(stage) {
  destination <- normalizePath(dirname(pipeline_outputs(stage)[1]),mustWork=FALSE)
  pipeline_hash(c(stage,destination))
}
pipeline_record <- function(stage,kind='success') {
  if (kind == 'submission') {
    return(file.path(pipeline_control_root,'submissions',paste0(pipeline_destination_key(stage),'.rds')))
  }
  file.path(pipeline_state,paste0(stage,'-',kind,'.rds'))
}
pipeline_registry_lock <- function() pipeline_lock(file.path(pipeline_control_root,'submissions','registry.lock'),wait=60)
pipeline_read <- function(path) if (file.exists(path)) tryCatch(readRDS(path),error=function(e) NULL) else NULL
pipeline_intent <- function(stage) {
  deps <- pipeline_deps[[stage]]
  private <- !stage %in% c('figure-1a','prepare-figure')
  helpers <- if (private) c('code/vwh_data_helpers.R','code/estimation_helpers.R','code/text_formatting_helpers.R') else 'code/driver_panel_helpers.R'
  software <- if (private) c('renv.lock','software/apm-source.json') else 'config/public.lock'
  settings <- as.list(Sys.getenv(c('FIRST_YEAR','LAST_YEAR','MAX_RANK','MIN_COHORT_SIZE','MIN_K','MAX_K',
      'CHOSEN_K','APM_SEED','APM_MAX_THREADS','DATA_PATH','OUTPUT_PATH')))
  if (pipeline_mode %in% c('submit','stage')) settings$APM_MAX_THREADS <- as.character(pipeline_threads[[stage]])
  lockfile <- if (private) 'renv.lock' else 'config/public.lock'
  package_key <- paste(lockfile,pipeline_files(lockfile),Sys.getenv('REPLICATION_LIBRARY'))
  if (!exists(package_key,pipeline_package_cache,inherits=FALSE)) {
    installed <- installed.packages(lib.loc=Sys.getenv('REPLICATION_LIBRARY'))[,c('Package','Version','Built'),drop=FALSE]
    package_names <- names(renv::lockfile_read(lockfile)$Packages)
    if (private) package_names <- c(package_names,'apm')
    assign(package_key,installed[installed[,'Package'] %in% package_names,,drop=FALSE],pipeline_package_cache)
  }
  installed <- get(package_key,pipeline_package_cache,inherits=FALSE)
  if (private) software <- c(software,file.path(Sys.getenv('REPLICATION_LIBRARY'),'apm','replication-source.json'))
  inputs <- if (!length(deps)) pipeline_files(pipeline_inputs(stage)) else lapply(deps,pipeline_intent)
  pipeline_hash(list(schema=1,stage=stage,inputs=inputs,
    code=pipeline_files(c('code/pipeline.R','config/paper.R','code/env_config.R','code/get_io_paths.R',
      file.path('code',pipeline_stages[[stage]]),helpers)),
    config=settings,
    software=pipeline_files(software), R=R.version.string,
    installed=installed))
}
pipeline_identity <- function(stage) list(intent=pipeline_intent(stage),inputs=pipeline_files(pipeline_inputs(stage)),
  upstream=lapply(pipeline_deps[[stage]],function(s) pipeline_read(pipeline_record(s))$generation))
pipeline_valid <- function(stage) {
  if (any(!vapply(pipeline_deps[[stage]],pipeline_valid,logical(1)))) return(FALSE)
  record <- pipeline_read(pipeline_record(stage))
  if (is.null(record)) return(FALSE)
  tryCatch(identical(record$identity,pipeline_identity(stage)) &&
    identical(record$outputs,pipeline_files(pipeline_outputs(stage))) &&
    all(file.info(pipeline_outputs(stage))$size > 0),error=function(e) FALSE)
}
pipeline_scheduler <- function(id) {
  sq <- Sys.which('squeue'); sa <- Sys.which('sacct')
  if (!nzchar(sq) || !nzchar(sa)) stop('Cannot establish scheduler status: squeue and sacct are required. No work submitted.')
  query <- function(bin,args,allow_absent=FALSE) {
    value <- suppressWarnings(system2(bin,args,stdout=TRUE,stderr=TRUE))
    if (!is.null(attr(value,'status'))) {
      if (allow_absent && identical(as.integer(attr(value,'status')),1L) &&
          any(grepl('Invalid job id specified',value,fixed=TRUE))) return(character())
      stop('Scheduler query failed for job ',id,': ',paste(value,collapse=' '))
    }
    trimws(value[nzchar(trimws(value))])
  }
  state <- query(sq,c('-h','-j',id,'-o','%T'),allow_absent=TRUE)
  if (!length(state)) state <- query(sa,c('-n','-X','-j',id,'--format=State','--parsable2'))
  if (!length(state)) stop('Scheduler has no status for recorded job ',id,'. Resolve registry before retrying.')
  state <- sub('[| +].*$','',state)
  active <- c('PENDING','RUNNING','CONFIGURING','COMPLETING','SUSPENDED','REQUEUED','RESIZING','REQUEUE_FED','REQUEUE_HOLD','SIGNALING','STAGE_OUT')
  terminal <- c('COMPLETED','FAILED','CANCELLED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
  if (any(!state %in% c(active,terminal))) stop('Unrecognized scheduler state for ',id,': ',paste(state,collapse=', '))
  any(state %in% active)
}
pipeline_active <- function(stage) {
  path <- pipeline_record(stage,'submission')
  if (!file.exists(path)) return(NULL)
  entry <- tryCatch(readRDS(path),error=function(e) stop('Unreadable submission record: ',path,'. Resolve it against scheduler records before retrying.'))
  scalar <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x)
  if (!is.list(entry) || !scalar(entry$id) || !grepl('^[0-9]+$',entry$id) ||
      !scalar(entry$intent) || !is.character(entry$dependencies) || anyNA(entry$dependencies))
    stop('Invalid submission record: ',path,'. Resolve it against scheduler records before retrying.')
  if (identical(as.character(entry$id),Sys.getenv('SLURM_JOB_ID'))) return(NULL)
  if (pipeline_scheduler(entry$id)) entry else NULL
}
pipeline_backup <- function(stage,attempt) {
  files <- pipeline_outputs(stage); files <- files[file.exists(files)]
  if (!length(files)) return(invisible(NULL))
  backup <- file.path(pipeline_state,'history',stage,attempt)
  dir.create(backup,recursive=TRUE,showWarnings=FALSE)
  for (i in seq_along(files)) {
    target <- file.path(backup,paste0(i,'-',basename(files[i])))
    if (!file.rename(files[i],target)) {
      if (!file.copy(files[i],target,copy.mode=TRUE,copy.date=TRUE)) stop('Cannot preserve baseline ',files[i])
      if (unlink(files[i]) != 0L) stop('Cannot clear preserved output ',files[i])
    }
  }
  pipeline_atomic(setNames(files,paste0(seq_along(files),'-',basename(files))),file.path(backup,'original-paths.rds'))
}
pipeline_execute <- function(stage,force=FALSE,intent=NULL) {
  decision <- pipeline_registry_lock()
  lock <- NULL
  tryCatch({
    active <- pipeline_active(stage)
    if (!is.null(active)) stop('Stage ',stage,' has active job ',active$id,'. Wait or resolve it before local/forced execution.')
    lock <- pipeline_lock(file.path(pipeline_control_root,'writer-locks',pipeline_destination_key(stage)))
  },finally=unlink(decision,recursive=TRUE))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  if (!is.null(intent) && !identical(intent,pipeline_intent(stage))) stop('Inputs, code, software or configuration changed since submission: ',stage)
  if (any(!vapply(pipeline_deps[[stage]],pipeline_valid,logical(1)))) stop('Prerequisite stage incomplete or stale for ',stage)
  if (!force && pipeline_valid(stage)) {message('Reuse ',stage,': validated identity and complete output hashes');return(invisible(NULL))}
  identity <- pipeline_identity(stage)
  old_attempt <- pipeline_read(pipeline_record(stage,'attempt'))
  resume <- !force && is.null(pipeline_read(pipeline_record(stage))) && !is.null(old_attempt) && identical(old_attempt$identity,identity)
  attempt <- if (resume) old_attempt$id else paste0(format(Sys.time(),'%Y%m%dT%H%M%S'),'-',Sys.getpid(),'-',basename(tempfile()))
  if (!resume) pipeline_backup(stage,attempt)
  unlink(pipeline_record(stage))
  pipeline_atomic(list(id=attempt,identity=identity,started=Sys.time()),pipeline_record(stage,'attempt'))
  message(if (resume) 'Resume ' else 'Run ',stage,' (attempt ',attempt,')')
  status <- system2('bash',c('code/jobs/run_stage.sh',shQuote(getwd()),shQuote(file.path('code',pipeline_stages[[stage]]))))
  if (status != 0L) stop('Stage failed: ',stage,'; attempt retained for compatible resume')
  outputs <- pipeline_files(pipeline_outputs(stage))
  if (any(file.info(names(outputs))$size <= 0)) stop('Empty output from ',stage)
  if (!identical(identity,pipeline_identity(stage))) stop('Stage inputs changed during execution; completion withheld: ',stage)
  pipeline_atomic(list(identity=identity,outputs=outputs,generation=attempt,finished=Sys.time()),pipeline_record(stage))
  message('Completed ',stage,': ',length(outputs),' outputs validated')
}
pipeline_main <- function(args) {
  pipeline_hash_cache <<- new.env(parent=emptyenv())
  pipeline_package_cache <<- new.env(parent=emptyenv())
  option <- function(name,default=NULL) {
    i <- match(name,args); if (is.na(i)) return(default)
    if (i == length(args) || startsWith(args[i+1],'--')) stop('Missing value for ',name)
    args[i+1]
  }
  unknown <- args[startsWith(args,'--') & !args %in% c('--target','--mode','--force','--stage','--intent')]
  if (length(unknown)) stop('Unknown option: ',paste(unknown,collapse=', '))
  target <- option('--target','figure-1a'); mode <- option('--mode','local')
  pipeline_mode <<- mode
  if (!mode %in% c('local','submit','stage')) stop('Mode must be local, submit or stage')
  if (!target %in% c(names(pipeline_stages),'vwh-all','private-figure')) stop('Unknown target: ',target)
  pipeline_deps[['figure-1a']] <<- character()
  if (target == 'private-figure') pipeline_deps[['figure-1a']] <<- 'prepare-figure'
  closure <- function(s) unique(c(unlist(lapply(pipeline_deps[[s]],closure)),s))
  requested <- if (target == 'vwh-all') names(pipeline_stages)[1:7] else if (target == 'private-figure') c('prepare-figure','figure-1a') else closure(target)
  force <- strsplit(option('--force',''),',',fixed=TRUE)[[1]]; force <- force[nzchar(force)]
  if (any(!force %in% c('all',names(pipeline_stages)))) stop('Unknown force stage: ',paste(setdiff(force,c('all',names(pipeline_stages))),collapse=', '))
  if (length(setdiff(force,c('all',requested)))) stop('Forced stage outside requested graph: ',paste(setdiff(force,c('all',requested)),collapse=', '))
  rerun <- setNames(rep(FALSE,length(requested)),requested)
  for (s in requested) rerun[s] <- 'all' %in% force || s %in% force || any(rerun[pipeline_deps[[s]]])
  source('config/paper.R'); source('code/env_config.R'); source('code/get_io_paths.R')
  scope <- if (all(requested %in% c('prepare-figure','figure-1a'))) 'public' else 'vwh'
  if (mode != 'stage') {
    status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','setup.R','--scope',scope))
    if (status != 0L) stop('Setup failed; no stages executed or submitted')
  }
  lib <- Sys.getenv('REPLICATION_LIBRARY'); Sys.setenv(R_LIBS_USER=lib); .libPaths(c(lib,.libPaths()))
  users <- paste0(lib,'.users'); dir.create(users,recursive=TRUE,showWarnings=FALSE)
  marker <- file.path(users,paste0(Sys.info()[['nodename']],'-',Sys.getpid()))
  library_lock <- pipeline_lock(paste0(lib,'.setup-lock'),wait=60)
  tryCatch(writeLines(as.character(Sys.time()),marker),finally=unlink(library_lock,recursive=TRUE))
  on.exit(unlink(marker),add=TRUE)
  pipeline_state <<- file.path(Sys.getenv('OUTPUT_PATH'),'.replication')
  dir.create(pipeline_state,recursive=TRUE,showWarnings=FALSE)
  for (s in requested) if (!length(pipeline_deps[[s]])) pipeline_files(pipeline_inputs(s))
  if (mode == 'stage') {
    stage <- option('--stage'); if (is.null(stage) || !stage %in% requested) stop('Invalid --stage')
    allocation <- suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')))
    if (!nzchar(Sys.getenv('SLURM_JOB_ID')) || is.na(allocation) || allocation != pipeline_threads[[stage]] ||
        num_threads != allocation) stop('Stage CPU allocation differs from the recorded profile; submit through the matching Slurm wrapper')
    pipeline_execute(stage,rerun[[stage]],option('--intent')); return(invisible(NULL))
  }
  if (mode == 'local') {
    lock <- pipeline_lock(file.path(pipeline_state,'local-pipeline.lock'))
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    for (s in requested) pipeline_execute(s,rerun[[s]])
    return(invisible(NULL))
  }
  registry <- pipeline_registry_lock()
  on.exit(unlink(registry,recursive=TRUE),add=TRUE)
  dir.create(file.path(pipeline_state,'logs'),recursive=TRUE,showWarnings=FALSE)
  ids <- setNames(rep('',length(requested)),requested)
  schedule <- setNames(rep(FALSE,length(requested)),requested)
  intents <- setNames(lapply(requested,pipeline_intent),requested)
  active <- setNames(lapply(requested,pipeline_active),requested)
  for (s in requested) {
    schedule[s] <- rerun[s] || any(schedule[pipeline_deps[[s]]]) || !pipeline_valid(s)
    if (!is.null(active[[s]])) {
      if (rerun[s] || !identical(active[[s]]$intent,intents[[s]])) stop('Conflicting active job ',active[[s]]$id,' for ',s,'. Resolve before retrying; force does not cancel jobs.')
      queued_deps <- pipeline_deps[[s]][schedule[pipeline_deps[[s]]]]
      if (length(queued_deps) && (any(!nzchar(ids[queued_deps])) ||
          !identical(unname(active[[s]]$dependencies[queued_deps]),unname(ids[queued_deps]))))
        stop('Active downstream job ',active[[s]]$id,' has incompatible upstream jobs for ',s,'. Resolve before retrying.')
      ids[s] <- active[[s]]$id; schedule[s] <- TRUE
    }
    if (dir.exists(file.path(pipeline_control_root,'writer-locks',pipeline_destination_key(s))) && is.null(active[[s]])) stop('Active local writer for ',s,'; wait before submitting')
  }
  for (s in requested) {
    if (nzchar(ids[s])) {message('Reuse active ',s,': job ',ids[s]);next}
    if (!schedule[s]) {message('Skip ',s,': validated identity and output inventory');next}
    deps <- ids[pipeline_deps[[s]]]; deps <- deps[nzchar(deps)]
    argv <- c('--parsable',paste0('--chdir=',shQuote(getwd())),
      paste0('--output=',shQuote(file.path(pipeline_state,'logs',paste0(s,'-%j.out')))),
      paste0('--error=',shQuote(file.path(pipeline_state,'logs',paste0(s,'-%j.err')))),
      if(length(deps)) paste0('--dependency=afterok:',paste(deps,collapse=':')),
      shQuote(file.path('code',paste0(s,'.slurm'))),shQuote(getwd()),'--target',target,'--stage',s,
      '--intent',intents[[s]],if(rerun[s]) c('--force',paste(force,collapse=',')))
    reply <- suppressWarnings(system2('sbatch',argv,stdout=TRUE,stderr=TRUE))
    if (!is.null(attr(reply,'status')) || length(reply)!=1 || !grepl('^[0-9]+(;[^[:space:]]+)?$',reply)) stop('Submission failed for ',s,': ',paste(reply,collapse=' '),'. Earlier successful submissions remain recorded.')
    ids[s] <- sub(';.*$','',reply)
    pipeline_atomic(list(id=ids[[s]],intent=intents[[s]],dependencies=ids[pipeline_deps[[s]]],submitted=Sys.time()),pipeline_record(s,'submission'))
    message('Submitted ',s,': job ',ids[s],if(length(deps)) paste0(' afterok ',paste(deps,collapse=':')) else '')
  }
  invisible(ids)
}
