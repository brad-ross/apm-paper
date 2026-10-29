# Replication package

Lihua Lei and Brad Ross, [“Estimating Counterfactual Matrix Means with Short Panel Data”](https://arxiv.org/abs/2312.07520).

This package supplies empirical code and approved numerical plotting inputs for Figure 1a. Public reproduction renders Figure 1a without restricted records. Full Veneto Workers Histories (VWH) reproduction requires separately authorized access. No data-availability exemption has been granted for this package; no journal reproducibility check or public deposit is claimed.

## Requirements and setup

Run commands from the extracted package root on a Unix-like system with Bash, Make and **R 4.3.0**. `renv.lock` pins the VWH dependencies; `config/public.lock` pins the smaller rendering environment. Setup restores packages into `.replication/library`, leaving matching installations in place. Internet access is needed for unavailable dependencies. Set `REPLICATION_LIBRARY` to an absolute path to use another library. Compute nodes must see the same library and files as the submission host. Setup serializes installation. Do not run setup or install packages into a library while jobs use it; wait for those jobs or select a different library.

Execution and submission targets automatically run scoped setup before work. Standalone setup is available:

```sh
make setup
Rscript setup.R --scope public
Rscript setup.R --scope vwh
```

Public setup requires only rendering dependencies. VWH setup also downloads, verifies and installs the exact `vendor/apm_0.1.0.tar.gz` when it is not already supplied. Its SHA-256, full upstream commit and packaging/install evidence are in `vendor/apm-source.json`; version `0.1.0` alone does not identify the build. No sibling checkout is required. The source asset is pinned to [apm v0.1.0](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication). See `vendor/README.md` for standalone installation and the source asset's MIT notices.

Source installation needs R development headers, Make, a C++17 compiler, BLAS/LAPACK and Fortran linking support. RcppParallel supplies APM's threading interface. Full VWH installation also needs system support for the locked Arrow package and other source-built dependencies. Lockfiles do not install operating-system libraries or compilers. Resolve installation errors against the relevant package's system requirements before retrying. Figure 1a rendering takes roughly one second on the tested Linux/R 4.3.0 environment. The full evaluation stage took roughly nine hours with the documented 250-CPU allocation; timings depend on hardware and scheduling.

`config/paper.R` defaults `DATA_PATH` to the project's `data` directory and `OUTPUT_PATH` to `output`. Use absolute environment overrides for external storage. Scientific defaults are years 1998–2001, two-year bins, rank 2, minimum cohort size 75, k=1–25 and chosen k=3. `APM_SEED` defaults to `9001`, with integer overrides from 0 through 2147483647. Clustering and the 1,000-draw bootstrap receive this seed; individual bootstrap calls can override it. `APM_MAX_THREADS` controls local thread allocation. Historical output seeds/environments are unconfirmed; the current default is not evidence of their seed.

We cluster all firms in the cleaned data within each province, then restrict
the simulation panel to the largest super-cohort identified by the O³ algorithm.
The simulation exercise evaluates the 15 largest cohorts with at least three
observed outcomes.

## Code organization

- `code/data/`: Stata conversion, clustered panels, and private driver-grid preparation.
- `code/analysis/`: evaluations and outcome decomposition.
- `code/reporting/`: figures, tables, and manuscript snippets.
- `code/lib/`: shared R configuration, data access, and analysis helpers.
- `code/jobs/`: shell stage runner, cluster environment, and Slurm submission.
- `code/slurm/`: resource requests and wrappers for individual jobs.

Run commands from the repository root. The Makefile owns the dependency graph;
R scripts contain the scientific work, and shell scripts handle execution.

## Public Figure 1a

The archive supplies `data/processed_data/driver_panel_grid.csv`; `data/processed_data/README.md` defines its rectangle geometry. Run:

```sh
make figure-1a
make figure-1a FORCE=all
```

The output is `output/result_summaries/figures/panel_plot_drivers_times.pdf` under default paths. Both commands read the supplied grid. Missing or invalid grids are incomplete-package errors: restore the supplied input. Rendering never opens `riders.csv`, generates the grid, or requires a private preparation record, including when forced or when preparation code is newer.

Separately authorized private preparation uses `make private-figure`, or the one-off cluster wrapper `sbatch code/slurm/driver_panel_figure.slurm` from the project root. `code/data/construct_driver_panel_grid.R` reads private `DATA_PATH/raw_data/riders.csv`, saves the grid, then `code/reporting/driver_panel_results.R` reads it in a separate process. Preparation failure stops rendering. This job is independent of VWH and unnecessary for public users. `make submit-figure` provides the explicit cluster submission interface.

Geometry and style were recovered from `eed1a22:code_old/visualize_panels.R:279-320`, correcting a dropped week-coordinate column. The public input contains cohort rectangles, not driver records. Public rendering verifies plotting from supplied inputs, not their upstream construction. Driver provider, collection period, acquisition version, citation and access arrangements require author confirmation. Grid release permission does not extend to other data.

## VWH access and execution

Apply through [Fondazione Rodolfo Debenedetti (FRDB)](https://www.frdb.org/en/dati/dati-inps-carriere-lavorative-in-veneto/): complete and sign its linked form and email it to `info@frdb.org`. Eligibility, approval, fees, delays, delivery and journal verification access require provider confirmation. The authors' acquisition date, extract version and agreement are unconfirmed. FRDB requests credit: “Dataset developed by the Economics Department in Università Ca’ Foscari Venezia under the supervision of Giuseppe Tattara”. See `docs/replication/data-sources.md` for coverage and provenance caveats.

Place `anagr.dta`, `azien.dta` and `contr.dta` in `DATA_PATH/raw_data/`. All three are required. Conversion writes `workers.parquet`, `firms.parquet` and `earnings.parquet` there. Firm clusters and clustered panels go under `processed_data`. Raw, converted and derived VWH records are restricted and excluded from the archive.

Run synchronously on the machine intended for analysis:

```sh
APM_MAX_THREADS=8 make vwh-all
make vwh-all FORCE=conversion
make vwh-all FORCE=all
```

The graph is conversion → panels, followed by panel reports and two branches: evaluation → evaluation reports, and decomposition → decomposition reports. Named targets are `conversion`, `panels`, `panel-reports`, `evals`, `eval-reports`, `decomposition` and `decomp-reports`; individual targets include prerequisites. Results go to `output/empirical_results` and `output/result_summaries/{figures,tables,snippets}`. `docs/replication/exhibits.csv` maps exact filenames, scripts and inputs.

On a configured Slurm cluster, submit from a login/submission host:

```sh
make submit-vwh
make submit-vwh FORCE=conversion
make submit-vwh FORCE=all
```

Submission loads the configured R module/toolchain before setup, prepares a shared project library and queues `afterok` dependencies. `REPLICATION_R_MODULE` defaults to `R/4.3.0`; adapt module configuration and resource profiles to your installation. Compute wrappers explicitly pass allocated CPUs to initialization. Local targets do not submit jobs. Submission returns job IDs and does not establish completion; setup failure must yield zero submissions.

Requested cluster resources (CPUs / memory / time) are conversion 4 / 16 GB / 4 hours; panels 40 / 50 GB / 24 hours; panel reports 16 / 50 GB / 12 hours; evaluation and decomposition each 250 / 50 GB / 2 days; their reports each 4 / 16 GB / 2 hours; private grid preparation 4 / 16 GB / 1 hour; public rendering 1 / 4 GB / 30 minutes. These requests are informed by existing cluster allocations where available, not measured successful runtimes. The scheduler must admit the requested per-job resources, including 250 CPUs. Local execution runs branches sequentially within the configured allocation, including under parallel Make. Do not start overlapping local runs or submissions against the same data/output directories.

The Makefile defines the dependency graph for both local and cluster execution. Successful stages write `OUTPUT_PATH/.replication/stamps/<stage>.done`; Make skips them until an upstream completion marker is newer or a stage is explicitly forced. Submission uses this same Make plan and queues downstream jobs with `afterok` dependencies. A queued prerequisite causes its requested dependents to be queued even if their old completion markers exist.

`FORCE=conversion` reruns conversion and its requested dependents; `FORCE=all` reruns the requested graph. Comma-separated stage names are supported. Interrupted evaluations can resume their saved batches when the panel has not rerun; forcing evaluations or rebuilding panels starts fresh. The evaluation script still checks that all expected results were produced before succeeding. A supplied driver grid skips private preparation unless `prepare-figure` is explicitly forced; public rendering never regenerates it.

Changes to code, seeds, settings, package versions, and output contents are **not** detected automatically. Force the affected stage when making such changes. If an output was removed, force its stage (or remove its `.done` marker). Markers belong to one data/output configuration: use a separate output directory or force a rebuild when switching inputs. Old RDS success records are ignored; the first invocation under this design runs stages without the new markers.

There is no automatic backup or active-job registry. Archive outputs yourself before forcing a rerun. Do not repeat a submission while its jobs are pending/running; use the saved job list to cancel blocked jobs before resubmitting. Each submission writes `OUTPUT_PATH/.replication/jobs-<timestamp>-<pid>.tsv` and job logs. A submission error may leave already-submitted jobs in the queue; their IDs remain in that file. Keep these private execution records out of the public deposit.

Preview the cluster plan without setup or submission:

```sh
bash code/jobs/submit_pipeline.sh --target vwh-all --dry-run
bash code/jobs/submit_pipeline.sh --target eval-reports --force evals --dry-run
```

`make submit-figure` is the explicit private preparation/rendering path; `make submit-public-figure` queues rendering from the supplied grid only.

## Verification, coverage and rights

Two full runs of the frozen scientific pipeline reproduced the same numerical results and rendered figures. Some results differ from the historical archive. The simplified Make/shell orchestration is checked with manual dry runs and synthetic stage/submission commands; the full scientific pipeline was not rerun for this orchestration change.

Authors must reconcile final manuscript/appendix coverage: Figure 1a is identified, while other exact references, external simulations and illustrations remain open in the exhibit map. The author confirmed on 2026-10-26 that the three legacy figure filenames are not referenced by the paper sources; they are excluded from required reproduction. Overleaf access is not required to run the package. Authors separately upload approved outputs using existing relative paths and compile the final manuscript to check imports; local rendering does not establish that check.

The [Econometric Society policy](https://www.econometricsociety.org/publications/es-data-editor-website/data-and-code-availability-policy) and [checklist](https://www.econometricsociety.org/uploads/reports%20Editorial/ES_Data_Editor_Website/Checklist.pdf) govern exemptions, transformation programs, verification and deposit. Authors must resolve the restricted-data arrangement with the Editorial Office/Data Editor before claiming submission readiness. The [Kline–Saggio–Sølvsten archive](https://eml.berkeley.edu/~pkline/papers/KSS2020_COMP.zip) is a documentation precedent for separating access, construction and analysis; its exemption does not transfer to this paper. No synthetic substitute or hosted driver-data access has been approved here.

The replication scripts are licensed under MIT; see `LICENSE`. This grant excludes data, numerical plotting inputs, documentation and third-party software. The APM asset retains its own MIT notices. The author approved the Figure 1a numerical plotting input only. Dataset citations, other disclosure permissions, final manuscript mapping and journal arrangements remain author release decisions.

## Local archive

After setup, with the approved grid present and Pandoc, pdfLaTeX, and `zip` installed on the packaging machine, run:

```sh
Rscript tools/build_replication_archive.R
```

The R builder uses `tools/public-files.txt`, rejects missing/unsafe/symlink paths, downloads the pinned APM asset if absent, verifies it and stages a fresh package under `dist`. It creates a journal-format ZIP archive and a tarball, per-file SHA-256 hashes, a listing and archive checksums. The builder generates the journal-required top-level `README.pdf` from `README.md` using Pandoc and pdfLaTeX; the PDF is not tracked in Git. Archive assembly uses R plus these document tools and the system `zip` utility. None of those packaging tools is needed to run replication. The [journal checklist](https://www.econometricsociety.org/uploads/reports%20Editorial/ES_Data_Editor_Website/Checklist.pdf) explicitly requests a PDF README; Markdown alone does not satisfy that checklist. Only allowlisted code, configuration/locks, documentation, the public CSV and frozen APM source are included. Raw/restricted records, outputs, Git history, development plans/verification logs, and local state are excluded. Building never publishes. Authors must review the exact listing and verification limitations before release; the intended journal deposit follows reproducibility checks and ES Journals Zenodo instructions. No DOI or deposit is claimed here.
