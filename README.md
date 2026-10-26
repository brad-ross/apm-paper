# Replication package

Lihua Lei and Brad Ross, [“Estimating Counterfactual Matrix Means with Short Panel Data”](https://arxiv.org/abs/2312.07520).

This package supplies empirical code and approved numerical plotting inputs for Figure 1a. Public reproduction renders Figure 1a without restricted records. Full Veneto Workers Histories (VWH) reproduction requires separately authorized access. No data-availability exemption has been granted for this package; no journal reproducibility check or public deposit is claimed.

## Requirements and setup

Run commands from the extracted package root on a Unix-like system with Bash, Make and **R 4.3.0**. `renv.lock` pins the VWH dependencies; `config/public.lock` pins the smaller rendering environment. Setup restores packages into `.replication/library`, leaving matching installations in place. Internet access is needed for unavailable dependencies. Set `REPLICATION_LIBRARY` to an absolute path to use another library. Compute nodes must see the same library and files as the submission host. Setup serializes installation; do not replace a library used by active jobs. Wait or select a different library if setup reports a conflict.

Execution and submission targets automatically run scoped setup before work. Standalone setup is available:

```sh
make setup
Rscript setup.R --scope public
Rscript setup.R --scope vwh
```

Public setup requires only rendering dependencies. VWH setup also verifies and installs the exact `vendor/apm_0.1.0.tar.gz`. Its SHA-256, full upstream commit and packaging/install evidence are in `software/apm-source.json`; version `0.1.0` alone does not identify the build. No sibling checkout is required. A public release tag/URL remains unassigned. See `software/README.md` for standalone installation and the source asset's MIT notices.

Source installation needs R development headers, Make, a C++17 compiler, BLAS/LAPACK and Fortran linking support. RcppParallel supplies APM's threading interface. Full VWH installation also needs system support for the locked Arrow package and other source-built dependencies. Lockfiles do not install operating-system libraries or compilers. Resolve installation errors against the relevant package's system requirements before retrying. Recorded hardware, platform, session information and measured Figure 1a runtime are in `docs/replication/verification.md`. VWH runtimes are unmeasured here.

`config/paper.R` defaults `DATA_PATH` to the project's `data` directory and `OUTPUT_PATH` to `output`. Use absolute environment overrides for external storage. Scientific defaults are years 1998–2001, two-year bins, rank 2, minimum cohort size 75, k=1–25 and chosen k=3. `APM_SEED` defaults to `9001`, with integer overrides from 0 through 2147483647. Clustering and the 1,000-draw bootstrap receive this seed; individual bootstrap calls can override it. `APM_MAX_THREADS` controls local thread allocation. Historical output seeds/environments are unconfirmed; the current default is not evidence of their seed.

## Public Figure 1a

The archive supplies `data/clean_data/driver_panel_grid.csv`; `data/clean_data/README.md` defines its rectangle geometry. Run:

```sh
make figure-1a
make figure-1a FORCE=all
```

The output is `output/result_summaries/figures/panel_plot_drivers_times.pdf` under default paths. Both commands read the supplied grid. Missing or invalid grids are incomplete-package errors: restore the supplied input. Rendering never opens `riders.csv`, generates the grid, or requires a private preparation record, including when forced or when preparation code is newer.

Separately authorized private preparation uses `make private-figure`, or the one-off cluster wrapper `sbatch code/driver_panel_figure.slurm` from the project root. `code/construct_driver_panel_grid.R` reads private `DATA_PATH/raw_data/riders.csv`, saves the grid, then `code/driver_panel_results.R` reads it in a separate process. Preparation failure stops rendering. This job is independent of VWH and unnecessary for public users. `make submit-figure` provides the explicit cluster submission interface.

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

Requested cluster resources (CPUs / memory / time) are conversion 4 / 16 GB / 4 hours; panels 40 / 50 GB / 24 hours; panel reports 16 / 50 GB / 12 hours; evaluation and decomposition each 250 / 50 GB / 2 days; their reports each 4 / 16 GB / 2 hours; private grid preparation 4 / 16 GB / 1 hour; public rendering 1 / 4 GB / 30 minutes. These requests are informed by existing cluster allocations where available, not measured successful runtimes. The scheduler must admit the requested per-job resources, including 250 CPUs. Local execution runs branches sequentially within the configured allocation, including under parallel Make; conflicting local invocations are blocked.

A stage is reusable only when success records, input/code/software/configuration identities and expected output hashes validate. Partial or outdated stages rerun; upstream reruns invalidate requested dependents. `FORCE=conversion` forces conversion and its dependents; `FORCE=all` forces the requested graph. A comma-separated list of recognized stage names is also supported. Force creates a fresh attempt instead of resuming old checkpoints and preserves historical attempts. Matching queued/running jobs are reused through submission records. Conflicting active writers block even forced requests; force never cancels jobs. Unknown scheduler state causes failure rather than duplicate submission. Keep project `.replication` state, `OUTPUT_PATH/.replication` run records and private logs out of a public deposit.

## Verification, coverage and rights

`docs/replication/verification.md` records dated checks and actual outcomes. Real-data verification is limited to Figure 1a preparation/rendering; other checks concern installation and temporary non-empirical dry runs. **The revised VWH conversion, clustering, estimation and reporting pipeline was not rerun or submitted. Numerical equivalence to historical paper outputs remains unverified.** No final full VWH run is implied. `docs/replication/baseline.md` distinguishes filename inventory from validated results and describes future comparison methods.

Authors must reconcile final manuscript/appendix coverage: Figure 1a is identified, while other exact references, three historical figures, external simulations and illustrations remain open in the exhibit map. Overleaf access is not required to run the package. Authors separately upload approved outputs using existing relative paths and compile the final manuscript to check imports; local rendering does not establish that check.

The [Econometric Society policy](https://www.econometricsociety.org/publications/es-data-editor-website/data-and-code-availability-policy) and [checklist](https://www.econometricsociety.org/uploads/reports%20Editorial/ES_Data_Editor_Website/Checklist.pdf) govern exemptions, transformation programs, verification and deposit. Authors must resolve the restricted-data arrangement with the Editorial Office/Data Editor before claiming submission readiness. The [Kline–Saggio–Sølvsten archive](https://eml.berkeley.edu/~pkline/papers/KSS2020_COMP.zip) is a documentation precedent for separating access, construction and analysis; its exemption does not transfer to this paper. No synthetic substitute or hosted driver-data access has been approved here.

The replication scripts are licensed under MIT; see `LICENSE`. This grant excludes data, numerical plotting inputs, documentation and third-party software. The APM asset retains its own MIT notices. The author approved the Figure 1a numerical plotting input only. Dataset citations, other disclosure permissions, final software release metadata, final manuscript mapping and journal arrangements remain author release decisions.

## Local archive

After setup, with `README.pdf` and the approved grid present, run:

```sh
Rscript tools/build_replication_archive.R
```

The R builder uses `release/public-files.txt`, rejects missing/unsafe/symlink paths, verifies the frozen APM asset and stages a fresh package under `dist`. It creates a journal-format ZIP archive and a tarball, per-file SHA-256 hashes, a listing and archive checksums. Archive assembly uses R plus the system `zip` utility; neither Python nor `zip` is needed to run the replication pipeline. Only allowlisted code, configuration/locks, documentation, the public CSV and frozen APM source are included. Raw/restricted records, outputs, Git history and local state are excluded. Building never publishes. Authors must review the exact listing and verification limitations before release; the intended journal deposit follows reproducibility checks and ES Journals Zenodo instructions. No DOI or deposit is claimed here.
