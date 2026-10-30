## Overview

This package reproduces the empirical results from [“Estimating Counterfactual Matrix Means with Short Panel Data”](https://arxiv.org/abs/2312.07520) by Lihua Lei and Brad Ross using the Veneto Workers Histories (VWH) employment dataset and driver behavior data from *Road Pricing Under Heterogeneity* by Ater, Ross, Shany and Vasserman. The VWH analysis groups similar firms, evaluates predictions of wages withheld from estimation, and estimates how much of the average wage differences between provinces can be attributed to province effects. The programs produce the figures, tables, and numerical results used in Figure 1b, Section 5, and Appendix G.

Replication has two parts. The VWH analysis requires data obtained separately from the provider and runs with `make vwh-all`. Figure 1a uses the supplied public plotting data and runs with `make figure_1a`; the underlying driver microdata are not needed. Both commands install the required R packages before computation. The sections below describe the data inputs, computing requirements, and execution steps. With the resource allocations reported below, the full analysis takes approximately 10 hours, excluding software installation and queue time.

## Data Availability and Provenance Statements

### Rights and availability

The public package excludes raw, converted, and derived VWH records. Access to the source data requires permission from Fondazione Rodolfo Debenedetti (FRDB), as described below. Permission is typically granted to economics researchers; see for example [Kline, Saggio, Solvsten (2020)](https://doi.org/10.3982/ECTA16410), also published in Econometrica.

The individual driver records underlying Figure 1a cannot be shared. We instead distribute the numerical plotting grid needed to reproduce that figure. These driver data are not used elsewhere in the analysis.

### Veneto Workers Histories (VWH)

VWH is derived from Italian social security (INPS) records covering people who ever lived or worked in any of Veneto's seven provinces during 1975–2001. Records identify worker–firm–year employment, weeks worked and total pay. The dataset was developed by the Economics Department at Università Ca’ Foscari Venezia under the supervision of Giuseppe Tattara.

The underlying VWH records cannot be distributed with this package. Researchers wishing to replicate the analysis should request them through the [Fondazione Rodolfo Debenedetti (FRDB) VWH access page](https://www.frdb.org/en/dati/dati-inps-carriere-lavorative-in-veneto/) and follow the provider's current application instructions and access conditions.

Obtain `contr.dta` (job spells), `anagr.dta` (worker demographics) and `azien.dta` (employer information), and place them in `data/raw_data/`. The conversion stage described below reads these files directly; Stata is not needed.

### Driver data and Figure 1a

The private source file is `riders.csv`. It cannot be shared, including with the replication team. Only the approved aggregate plotting input, `driver_panel_grid.csv`, is supplied. The raw records are used solely to construct Figure 1a and are not inputs to the VWH analysis.

The driver data come from *Road Pricing Under Heterogeneity* by Ater, Ross, Shany and Vasserman. See that paper for details about the driver data. The underlying private records are unavailable through this package.

The grid contains anonymous cohort rectangles with at least 100 drivers per cohort. Its columns and construction are documented in [data/processed_data/README.md](data/processed_data/README.md). The plotting program reads this grid directly. A separate preparation program constructs it from the private records.

## Dataset List

Paths below are relative to `DATA_PATH`, which defaults to `data/`.

| File | Contents/source | Availability |
|---|---|---|
| `raw_data/anagr.dta` | VWH worker registry | Obtain from provider |
| `raw_data/azien.dta` | VWH firm registry | Obtain from provider |
| `raw_data/contr.dta` | VWH employment/contribution records | Obtain from provider |
| `raw_data/riders.csv` | Individual driver records | Not provided |
| `processed_data/driver_panel_grid.csv` | Figure 1a rectangle coordinates | Provided |

The converted and derived VWH files must not be redistributed without permission. The driver-grid codebook is included. Request documentation for the original VWH files from the provider. See [data/raw_data/README.md](data/raw_data/README.md) for the conversion encoding and input-field notes.

## Description of Code

### Stage dependencies

The diagram shows the job dependencies. An arrow means that the preceding stage must finish successfully before the next stage runs. Setup installs the required R packages before any computation.
```text
setup
  |
  +--> conversion --> panels --+--> panel_reports
  |                           |
  |                           +--> evals --> eval_reports
  |                           |
  |                           +--> decomposition --> decomp_reports
  |
  +--> figure_1a  <-- supplied driver_panel_grid.csv

Optional private preparation (requires riders.csv):
setup --> figure_1a_data_prep --> driver_panel_grid.csv --> figure_1a
```
`make vwh-all` runs the VWH branch, and `make figure_1a` renders the driver figure from the supplied grid. Make automatically runs the required stages in dependency order; individual scripts do not need to be run manually. It skips completed stages according to the completion and rerun rules below. Local execution is sequential. The Slurm submission commands use the same graph and queue jobs with dependencies, allowing independent branches to run concurrently.

The private preparation branch is used only by `make private-figure` or `make figure_1a_data_prep`. An existing grid skips preparation unless explicitly forced. Public Figure 1a replication does not require private preparation or access to `riders.csv`.

### Stage descriptions and code

Stage names below are Make targets. Each stage runs through a shell script in `code/jobs/`; the corresponding `.slurm` script in `code/slurm/` specifies its cluster resources. The resource table above lists those script pairs. Stage names match script filenames, for example `panel_reports` uses `panel_reports.sh`.

**Setup.** `setup.R` installs the dependencies in `renv.lock` and the specified APM source package into `.replication/library/` (or `REPLICATION_LIBRARY` if set). Downloaded APM source is saved as `vendor/apm_0.1.0.tar.gz`. Make invokes setup automatically before running computation.

**Conversion (`conversion`).** `code/data/raw_vwh_dta_to_parquet.R` reads the three raw VWH Stata files and writes `workers.parquet`, `firms.parquet`, and `earnings.parquet` in `data/raw_data/`. The remaining VWH stages read these converted files rather than the Stata originals.

**Panel construction (`panels`).** `code/data/construct_clustered_panels.R` reads the converted employment and firm records and writes two files in `data/processed_data/`: `firm_clusters_start_year=1998_end_year=2001_year_cluster_size=2.parquet` and `panels_with_clustered_outcomes_start_year=1998_end_year=2001_year_cluster_size=2_max_rank=2.parquet` with the default settings.

As described in Section 5.1 and Appendix G.1, we use observations from 1998–2001 on workers employed by Veneto firms. We aggregate observations into two-year periods and iteratively retain firms and workers observed in both periods. Within each province (BL, PD, RO, TV, VE, VI and VR), we cluster firms into three types using their 1998–1999 empirical wage distributions, evaluated at 100 quantiles. We retain the best solution from 10 clustering initializations. An outcome is the log of a worker’s average weekly wage for a two-year period (1998–1999 or 2000–2001), province, and firm type. A cohort consists of workers with the same set of observed outcomes. We drop cohorts with fewer than 75 workers and restrict estimation to the largest group of cohorts connected by the O³ Algorithm (Algorithm 2 in the paper) after clustering firms. This group, called a super cohort, determines which counterfactual means can be identified under the model.

**Panel reporting (`panel_reports`).** `code/reporting/construct_clustered_panels_results.R` reads the constructed panels, firm assignments, and converted VWH records. It produces Figure 1b, Figures G.1–G.5, Table G.1, and sample-description snippets, including worker and cohort counts. Figures are saved in `output/result_summaries/figures/`, the table in `output/result_summaries/tables/`, and numerical snippets in `output/result_summaries/snippets/`. The plots describe observation patterns, firm cluster sizes, cohort sizes, and identification across different cluster counts.

**Estimator evaluation (`evals`).** `code/analysis/leave_cohort_outcome_out_evals.R` reads the panels and firm assignments and saves the evaluation results to `output/empirical_results/masked_metrics_start_year=1998_end_year=2001_year_cluster_size=2_min_cohort_size=75_k=3.parquet` with the default settings. This file is also updated as evaluation batches complete.

In the semi-synthetic simulation study of Section 5.2, we mask one observed outcome at a time in the 15 largest eligible cohorts with at least three observed outcomes. We compare the Two-Way Fixed Effects (TWFE) estimator with four factor-model specifications: ranks 1 and 2, with and without outcome fixed effects. We evaluate bias, standard error, and root mean squared error (RMSE) across 1,000 Bayesian bootstrap draws, treating the observed but masked cohort-outcome mean as the target.

**Evaluation reporting (`eval_reports`).** `code/reporting/leave_cohort_outcome_out_evals_results.R` summarizes all saved evaluation results. It produces Table 1a and in-text snippets comparing the factor-model estimators with TWFE using geometric mean error ratios and the shares of targets for which each estimator performs better. The table is saved to `output/result_summaries/tables/leave_cohort_outcome_out_evals_results.tex`, and the snippets to `output/result_summaries/snippets/`.

**Wage decomposition (`decomposition`).** `code/analysis/match_outcome_decomposition.R` reads the panels and firm assignments and estimates the contributions of province effects and worker sorting to average wage differences between province pairs (Section 5.3). It saves the estimates and bootstrap results to `output/empirical_results/match_outcome_decomp_results.parquet`.

**Decomposition reporting (`decomp_reports`).** `code/reporting/match_outcome_decomposition_results.R` reads the saved decomposition results and produces Table 1b and in-text snippets for province-effect attribution estimates and uncertainty. The table is saved to `output/result_summaries/tables/match_outcome_decomp_results.tex`, and the snippets to `output/result_summaries/snippets/`.

**Private driver-grid preparation (`figure_1a_data_prep`).** `code/data/construct_driver_panel_grid.R` reads private `data/raw_data/riders.csv` and writes `data/processed_data/driver_panel_grid.csv`. The grid records the rectangles needed to display driver observation patterns without exporting individual driver records. This optional stage is not needed when the supplied grid is present.

**Driver figure (`figure_1a`).** `code/reporting/driver_panel_results.R` reads `data/processed_data/driver_panel_grid.csv` and writes Figure 1a to `output/result_summaries/figures/panel_plot_drivers_times.pdf`.

The reporting stages write figures, tables, and in-text snippets under `output/result_summaries/`. The exhibit mapping below gives the exact output filenames and distinguishes outputs reproducible from the supplied data from those requiring VWH access. Paths above are relative to the package root and use the default settings. Setting `DATA_PATH` or `OUTPUT_PATH` replaces the `data/` or `output/` prefix, respectively. Filenames containing analysis settings change when those settings change.

### Shared configuration and supporting code

`Makefile` defines the stage dependencies and rerun rules. `config/paper.R` specifies paths, analysis settings, and the random seed. Shared R functions are in `code/lib/`. The shell scripts use `code/jobs/common.sh` to configure execution and record stage completion; `code/jobs/submit_pipeline.sh` submits the selected stages to Slurm, and `code/jobs/init_server.sh` initializes the cluster R environment. `tools/build_replication_archive.R` assembles the distribution archive and is not part of the analysis graph.

Replication scripts are licensed under MIT; see [LICENSE](LICENSE). APM retains its own MIT notices. The code license does not grant rights to data, numerical plotting inputs, documentation or other third-party material.

## Computational Requirements

### Software

The analysis requires **R 4.3.0**, Bash, and GNU Make on a Unix-like system. All stages were run in a shared Ubuntu Linux environment on x86_64 servers. The recorded operating-system version was Ubuntu 22.04.5 LTS. Stata is not required because the R package Haven reads the original `.dta` files.

Install the system tools needed to compile R packages if they are not already present before starting replication: R development headers, a C++17 compiler (we used GNU g++ 11.4.0), Fortran linking support, and BLAS/LAPACK (we used R's bundled BLAS and LAPACK 3.10.0). Arrow and other R packages may require additional system libraries.

The replication commands then run `setup.R` to install the needed R package dependencies at the versions listed in `renv.lock` automatically. Setup also installs the [apm v0.1.0 release](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication), which implements the paper's Aggregated Projection Matrix (APM) estimator. Its source archive and checksum are specified in [vendor/apm_source.json](vendor/apm_source.json). The replication package already includes the APM source; setup downloads it from the [associated GitHub release](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication) if absent and checks its checksum before installation. As such, a separate copy of the APM repository is not required. `setup.R` installs packages into `.replication/library` by default and requires internet access for dependencies that are not already available. Instructions for changing the library location are given under “Paths, reruns and recovery.”

The distributed README includes an [R package-version appendix](#appendix-r-package-versions) generated from `renv.lock` during archive assembly; the lockfile remains the authoritative dependency record.

Slurm is optional. It is used only for the cluster execution commands below; local execution requires no scheduler. Both execution modes use the same analysis programs and R packages.

### Hardware, runtime and storage

The complete pipeline ran on the same shared computing system, using the software environment described above. The table reports the CPU and memory allocations, time limits, and observed runtimes for each stage. This is a known working configuration, not an established minimum hardware requirement. CPU counts refer to scheduler allocations and need not correspond to physical cores.

Script names below are relative to `code/jobs/` for `.sh` files and `code/slurm/` for `.slurm` files.

| Stage | Shell / Slurm scripts | Requested CPUs | Requested memory | Time limit | Observed elapsed time |
|---|---|---:|---:|---:|---|
| Conversion | `conversion{.sh/.slurm}` | 4 | 16 GB | 4 hours | 18 min 40 sec |
| Clustering / panels | `panels{.sh/.slurm}` | 40 | 50 GB | 24 hours | 4 min 32 sec |
| Panel reports | `panel_reports{.sh/.slurm}` | 16 | 50 GB | 12 hours | 4 min 38 sec |
| Evaluations | `evals{.sh/.slurm}` | 250 | 50 GB | 2 days | 9 hr 10 min 26 sec |
| Evaluation reports | `eval_reports{.sh/.slurm}` | 4 | 16 GB | 2 hours | 27 sec |
| Decomposition | `decomposition{.sh/.slurm}` | 250 | 50 GB | 2 days | 5 min 37 sec |
| Decomposition reports | `decomp_reports{.sh/.slurm}` | 4 | 16 GB | 2 hours | 11 sec |
| Figure 1a | `figure_1a{.sh/.slurm}` | 1 | 4 GB | 30 min | Approximately 1–2 sec |
| Private grid preparation | `figure_1a_data_prep{.sh/.slurm}` | 4 | 16 GB | 1 hour | Less than 1 sec in the measured R session |

The resource requests specify the limits configured for these jobs; we have not measured minimum memory requirements. Local execution is sequential. Under Slurm, independent stages can run concurrently after clustering. Total elapsed time also depends on the job queue, software installation, and hardware.

A working data directory with the raw and generated data occupied approximately 3.7 GB, and the installed R library approximately 535 MB. We recommend allowing at least 10 GB of working disk space, with additional space for package compilation, caches, and archived runs. This recommendation does not reflect a measurement of peak temporary storage.

## Instructions to Replicators

Run all commands from the package root, the directory containing `Makefile` and `setup.R`. Choose either local execution or Slurm execution. Both follow the same sequence: convert the raw VWH files, construct the panels, run the analyses, and generate the reports. Figure 1a is independent of these stages.

### Local execution

1. Obtain the VWH inputs and put the three `.dta` files in `data/raw_data/`. Leave the supplied grid in `data/processed_data/`.
2. Run the VWH pipeline with `APM_MAX_THREADS` set to the number of CPUs available to your run. For example, if you allocate eight CPUs:

   ```sh
   APM_MAX_THREADS=8 make vwh-all
   ```

3. Render Figure 1a:

   ```sh
   make figure_1a
   ```

The runtime table reports timings with the listed per-stage cluster resource allocations. Make automatically runs setup but standalone environment installation is available with `make setup` or `Rscript setup.R`.

After panel construction, panel reporting, estimator evaluation, and wage decomposition can proceed independently. Each analysis is followed by its reporting stage. Make includes the prerequisites of each requested stage.

The programs write outputs to `output/empirical_results/` and `output/result_summaries/{figures,tables,snippets}/`. Upload the latter directory to the manuscript project with its relative paths intact. Manuscript compilation takes place separately; this package does not require Overleaf access.

### Cluster execution

Place the input data in the same locations described for local execution. The submission host and compute nodes must share the code, data, output, and R library paths. Use a shell where the `module` command is initialized, set `REPLICATION_R_MODULE` to your installation's R module name (default `R/4.3.0`), and adjust resource requests in `code/slurm/` for your scheduler. Then run:

```sh
make submit-vwh
make submit-public-figure
```

The submission program loads `REPLICATION_R_MODULE` (default `R/4.3.0`), installs the required software, and submits stages with Slurm `afterok` dependencies. The Slurm scripts pass the allocated CPU count to the initialization script. Submission returns job identifiers before computation finishes.

Preview the stages without installing packages or submitting jobs:

```sh
bash code/jobs/submit_pipeline.sh --target vwh-all --dry-run
```

### Paths, reruns and recovery

The project directory and configured data, output, and library paths must not contain whitespace.

To use different storage locations, set `DATA_PATH`, `OUTPUT_PATH`, or `REPLICATION_LIBRARY` to absolute paths before running Make. Do not install packages into a library that running jobs are using; wait for those jobs to finish or select a separate library. The analysis settings described below can also be changed through environment variables.

Each successfully completed stage creates `OUTPUT_PATH/.replication/stamps/<stage>.done`. Make uses this record to omit completed stages unless a prerequisite has been rerun or a stage is explicitly forced using `FORCE`:

```sh
make vwh-all FORCE=conversion
make vwh-all FORCE=all
make eval_reports FORCE=evals
make figure_1a FORCE=all
```

Comma-separated stage names are accepted. Interrupted evaluations can resume from saved batches if panel construction has not been rerun. Forcing evaluation or reconstructing the panels starts evaluation again from the beginning. Completion requires all expected evaluations.

Changes to code, settings, package versions, data contents or missing outputs are not detected automatically. Force the affected stages when those change. Use a separate output directory or force a rebuild when switching datasets.

Archive outputs before a forced rerun if you need to retain them. Do not overlap runs or submissions against the same paths. There is no active-job registry or automatic backup. Submission writes job IDs to `OUTPUT_PATH/.replication/jobs-<timestamp>-<pid>.tsv` and logs under `OUTPUT_PATH/.replication/logs/`. A submission failure can leave earlier jobs queued; inspect the saved job list and cancel blocked jobs before resubmitting.

### Random seeds and analysis settings

`config/paper.R` sets the default seed used everywhere in the replication package to `APM_SEED=9001`; `code/lib/env_config.R` parses and validates it. Firm clustering and the 1,000-draw bootstrap receive that seed. Keep the seed, package versions and thread settings fixed when comparing runs.

The default settings are years 1998–2001, two-year periods, model rank 2, a minimum cohort size of 75, and firm cluster counts from 1 to 25. The main analysis uses three firm clusters per province. We cluster firms within provinces before restricting the estimation panel to the largest super cohort of the O³ Algorithm. The evaluation targets the 15 largest cohorts with at least three observed outcomes.

### Private driver-grid preparation

This step is unnecessary with the supplied grid and requires access to private `data/raw_data/riders.csv`:

```sh
make private-figure
# Cluster alternative:
make submit-figure
```

Preparation saves the grid before a separate rendering process reads it. Preparation is omitted when the grid already exists, unless `FORCE=figure_1a_data_prep` or `FORCE=all` is supplied to the private target. The public `figure_1a` target reads only the supplied grid.

## List of Outputs and Programs

Only Figure 1a is reproducible from the public data supplied here; the VWH outputs require data access. This table lists only exhibits used in the manuscript. Figure 2 is a TikZ illustration maintained in the manuscript source and requires no replication-package output.

Names are relative to `output/result_summaries/`. `{P}` denotes `start_year=1998_end_year=2001_year_cluster_size=2`, `{S}` denotes `{P}_min_cohort_size=75`, and `{K}` denotes `{S}_k=3`. These are abbreviations in this table only.

| Manuscript reference | Output | Reporting program in `code/reporting/` |
|---|---|---|
| Figure 1a | `figures/panel_plot_drivers_times.pdf` | `driver_panel_results.R` |
| Figure 1b | `figures/match_outcome_panel_structure_plot_{K}.pdf` | `construct_clustered_panels_results.R` |
| Table 1a | `tables/leave_cohort_outcome_out_evals_results.tex` | `leave_cohort_outcome_out_evals_results.R` |
| Table 1b | `tables/match_outcome_decomp_results.tex` | `match_outcome_decomposition_results.R` |
| Figure G.1 | `figures/firms_by_province_cluster_dot_plot_{P}_k=3.pdf` | `construct_clustered_panels_results.R` |
| Figure G.2 | `figures/cohort_size_dist_plot_{K}.pdf` | `construct_clustered_panels_results.R` |
| Table G.1 | `tables/panel_summary_stats_{K}.tex` | `construct_clustered_panels_results.R` |
| Figure G.3a | `figures/largest_super_cohort_share_workers_plot_{S}_rank=2.pdf` | `construct_clustered_panels_results.R` |
| Figure G.3b | `figures/largest_super_cohort_share_firms_plot_{S}_rank=2.pdf` | `construct_clustered_panels_results.R` |
| Figure G.4 | `figures/num_o3_iterations_needed_plot_{S}_rank=2.pdf` | `construct_clustered_panels_results.R` |
| Figure G.5 | `figures/block_missingness_id_firm_weighted_plot_{K}.pdf` | `construct_clustered_panels_results.R` |

Upload `result_summaries/` with its structure intact to include them in a new compilation of the paper.

### In-text results

The following inventory locates the numerical results imported into the manuscript text. Locating phrases identify the surrounding prose, with ellipses in place of computed values. Filenames are relative to `output/result_summaries/snippets/`. Repeated uses are listed separately so that each occurrence can be located. All results in this inventory require VWH access.

#### Section 5.1: Setting and Sample Construction

Generated by `code/reporting/construct_clustered_panels_results.R`.

| Locating phrase and statistic | Snippet file |
|---|---|
| “yielding T = … outcomes” — number of outcomes | `num_outcomes.txt` |
| “Dropping cohorts with fewer than … workers” — minimum cohort size | `min_cohort_size.txt` |
| “leaves N = … workers” — estimation sample size | `num_workers.txt` |
| “in C = … cohorts” — number of cohorts | `num_cohorts.txt` |
| “with … observed outcomes” — number of observed worker-outcome records | `num_obs.txt` |
| “an average of … observed outcomes per worker” | `avg_observed_outcomes_per_worker.txt` |
| “converges after … iterations” — O³ convergence for ranks 1 and 2 | `num_o3_iterations_rank=1.txt`, `num_o3_iterations_rank=2.txt` |
| “comprising all cohorts with at least … workers” — repeated cohort-size threshold | `min_cohort_size.txt` |

#### Section 5.2: Semi-Synthetic Simulation Study of Estimator Performance

Generated by `code/reporting/leave_cohort_outcome_out_evals_results.R`.

| Locating phrase and statistic | Snippet file |
|---|---|
| “has lower bias than the TWFE estimator for …” — share of targets with lower bias | `pc_1_no_fes_equal_weights_bias_share_better.txt` |
| “and lower RMSE for …” — share of targets with lower RMSE | `pc_1_no_fes_equal_weights_rmse_share_better.txt` |
| “with average bias and RMSE reductions of …” — bias reduction | `pc_1_no_fes_equal_weights_bias_avg_pct_better.txt` |
| “with average bias and RMSE reductions of …” — RMSE reduction | `pc_1_no_fes_equal_weights_rmse_avg_pct_better.txt` |

#### Section 5.3: Match Outcome Decomposition

Generated by `code/reporting/match_outcome_decomposition_results.R`.

| Locating phrase and statistic | Snippet file |
|---|---|
| “provinces causally contribute … percentage points more” — difference between APM and TWFE attribution estimates | `fgw_decomp_province_difference_estimate.txt` |

#### Section 6.2: From Short to Sparse Panels

Generated by `code/reporting/construct_clustered_panels_results.R`.

| Locating phrase and statistic | Snippet file |
|---|---|
| “has T = … outcomes” — repeated outcome count | `num_outcomes.txt` |
| “but C = … cohorts” — repeated cohort count | `num_cohorts.txt` |
| “of workers in the largest cohort” — largest cohort’s worker share | `largest_cohort_worker_share.txt` |
| “cohorts with fewer than 100 workers” — number of small retained cohorts | `num_cohorts_below_100_workers.txt` |

#### Appendix G.1: Dataset Construction Details

Generated by `code/reporting/construct_clustered_panels_results.R`.

| Locating phrase and statistic | Snippet file |
|---|---|
| “cluster the firms located in each province into k = … types” | `chosen_k.txt` |

The worker share and small-cohort count in Section 6.2 describe cohorts retained after the minimum-size restriction and before selection of the largest super cohort, matching Figure G.2. The share uses retained workers as its denominator. Figure G.3a instead divides by all workers in the clustered panel before the minimum-cohort-size restriction.

Other empirical statements summarize results displayed in tables or figures rather than importing a single number:

| Location | Locating phrase | Supporting output |
|---|---|---|
| Section 5.1 | “a single super cohort” and “the entire C × T matrix … is identified” | Figures G.3a–G.3b and G.4 |
| Section 5.1 | “no embedded block missingness patterns” | Figure G.5 |
| Section 5.2 | “best-performing variant,” “despite higher variance,” and “nearly double the bias” | Table 1a |
| Section 5.3 | “similarly imprecise estimates” and “the 95% bootstrap confidence interval” | Table 1b |
| Section 6.2 | “the largest cohort shrinks rapidly” and decreasing shares in the largest terminal super cohort | Figures G.3a–G.3b |

The exhibit mapping above identifies the programs and filenames for these outputs.

Intermediate estimates in `output/empirical_results/` are `masked_metrics_start_year=1998_end_year=2001_year_cluster_size=2_min_cohort_size=75_k=3.parquet` and `match_outcome_decomp_results.parquet`, generated by the evaluation and decomposition analysis scripts respectively. The reporting programs read these estimates to construct the manuscript outputs; the manuscript does not import them directly.

## References

- Università Ca’ Foscari Venezia, Economics Department, under the supervision of Giuseppe Tattara. n.d. *Veneto Workers Histories (VWH)* [dataset], based on INPS administrative records; distributed through Fondazione Rodolfo Debenedetti. [Data access page](https://www.frdb.org/en/dati/dati-inps-carriere-lavorative-in-veneto/).
- Ater, Itai, Brad Ross, Adi Shany, and Shoshana Vasserman. 2026. *Road Pricing Under Heterogeneity*. Working paper.
- Kline, Patrick, Raffaele Saggio, and Mikkel Sølvsten. 2020. “Leave-Out Estimation of Variance Components.” *Econometrica*, 88(5), 1859–1898. [Replication package](https://eml.berkeley.edu/~pkline/papers/KSS2020_COMP.zip).
- Lei, Lihua, and Brad Ross. *Estimating Counterfactual Matrix Means with Short Panel Data*. [arXiv:2312.07520](https://arxiv.org/abs/2312.07520).
- Ross, Brad, and Lihua Lei. *APM: The Aggregated Projection Matrix Method*, R package 0.1.0. [Source release](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication).

## Archive assembly for authors

With the approved figure 1a grid data present as `data/processed_data/driver_panel_grid.csv`, run:

```sh
make archive
```

This target runs setup and then `tools/build_replication_archive.R` to build the replication package `.zip` file.

Packaging additionally requires Pandoc, pdfLaTeX and the system `zip` utility; these are not required for analysis. The builder generates a top-level `README.pdf` from this Markdown file as well. The archive program reads the file list in `tools/public_files.txt`, verifies the APM checksum, and writes ZIP/tar archives and file checksums under `dist/`.
