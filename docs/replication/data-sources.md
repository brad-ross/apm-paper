# Data sources, access, and distribution status

Inventory date: 2026-10-26. This document separates obtaining source files, constructing analysis data, and running analysis, following the structure of the Kline–Saggio–Sølvsten (KSS) replication documentation linked below. KSS's exemption does not apply to this paper. No exemption is currently approved for this package. Actual access agreements, access dates, source versions, costs, delays, and journal verification arrangements require author/provider confirmation.

## VWH: obtain source files

Source: INPS administrative archives, Veneto Workers Histories (VWH), distributed through Fondazione Rodolfo Debenedetti (FRDB). The [FRDB source and access page](https://www.frdb.org/en/dati/dati-inps-carriere-lavorative-in-veneto/) was checked on 2026-10-26; the page is dated 2015-03-04. Follow its attached application form: complete it, sign it, and send it to `info@frdb.org`. This instruction is an access application, not a guarantee of access. Eligibility, approval criteria, fees, expected waiting time, delivery method, journal access, and permitted use/disclosure must be confirmed with FRDB. The page does not establish these facts for this package.

The three source files are:

| File under configured `DATA_PATH/raw_data/` | Source-file role | Conversion output in the same directory |
|---|---|---|
| `anagr.dta` | Worker registry | `workers.parquet` |
| `azien.dta` | Firm registry | `firms.parquet` |
| `contr.dta` | Contribution/employment records | `earnings.parquet` |

FRDB describes Treviso/Vicenza firms and workers followed beyond those provinces through 2001. The code permits BL, PD, RO, TV, VE, VR, and VI. Authors must confirm the acquired extract, its coverage and version before treating those descriptions as equivalent. Use FRDB's requested credit: “Dataset developed by the Economics Department in Università Ca’ Foscari Venezia under the supervision of Giuseppe Tattara”. [Source](https://www.frdb.org/en/dati/dati-inps-carriere-lavorative-in-veneto/)

The authors' acquisition date, agreement, extract identifier, file versions, and exact data citation/identifier are **unconfirmed**. The web-page inspection date above is not the date the authors obtained data. Do not replace missing provenance with the page date or code's analysis years.

## VWH: construct analysis data

`code/raw_vwh_dta_to_parquet.R` reads all three Stata files and writes the three Parquet files listed above. `workers.parquet` is produced even though no downstream worker-registry read was found in the inspected reporting/estimation scripts. The full workflow still requires all three source files.

`code/construct_clustered_panels.R`, using `code/vwh_data_helpers.R`, reads earnings and firm Parquet files; filters years/provinces and positive paid weeks; bins years; computes weekly earnings; keeps positive earnings and iteratively retains always-present workers/firms; clusters firms by province using the first year bin; then constructs panels with clustered outcomes. It writes under `DATA_PATH/processed_data/`:

- `firm_clusters_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}.parquet`
- `panels_with_clustered_outcomes_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_max_rank={MAX_RANK}.parquet`

Treat raw, converted, and derived VWH records as restricted. They are excluded from a public deposit. Anonymous identifiers, clustering, or aggregation do not establish redistribution rights. Existing figures, tables, snippets, logs, and fingerprints also require a disclosure decision; this inventory is not approval to release their contents.

## VWH: run analysis and reports

Panel reports consume the two processed files; their summary-statistics section additionally reads raw earnings/firm Parquet files. Evaluation consumes clustered panels and writes `OUTPUT_PATH/empirical_results/masked_metrics_start_year={FIRST_YEAR}_end_year={LAST_YEAR}_year_cluster_size={YEAR_CLUSTER_SIZE}_min_cohort_size={MIN_COHORT_SIZE}_k={CHOSEN_K}.parquet`. Its report also requires firm clusters. Decomposition consumes both processed files and writes `OUTPUT_PATH/empirical_results/match_outcome_decomp_results.parquet`. Each result script produces the exact exhibit files in [exhibits.csv](exhibits.csv).

No VWH stage was run for this inventory. Existing output filenames are evidence of candidate files only. They do not establish successful execution, correct content, matching data/software versions, or manuscript use.

## Driver example: private source and approved public plotting input

Private source: `DATA_PATH/raw_data/riders.csv`. Original source/provider, collection period, acquisition date/version, citation, license, and access/contact arrangements are **unconfirmed**; no hosted access route is promised. The author permits release of the numerical plotting inputs for Figure 1a, not the driver records. Raw-to-grid transformation code must accompany the package; approved plotting inputs do not by themselves resolve the journal's raw-data requirements.

Historical source is `eed1a22:code_old/visualize_panels.R:279-320`, deleted in `de93961`. It reads `rider_id`, `first_week_in_panel_dateweek_id`, and `began_active_week_id`; removes missing values; expands inclusive entry-to-treatment weeks; retains cohorts with at least 100 drivers; and orders anonymous display rows by the historical cohort labels. The historical script drops `week_id` before plotting; the recovered code repairs this while preserving the intended geometry.

The package supplies `code/construct_driver_panel_grid.R` for explicit private preparation and `code/driver_panel_results.R` for rendering. The supplied `DATA_PATH/clean_data/driver_panel_grid.csv` contains rectangle geometry, with its codebook in `data/clean_data/README.md`. Public reproduction starts at that file and writes `OUTPUT_PATH/result_summaries/figures/panel_plot_drivers_times.pdf`; it must not open private input or automatically regenerate a missing grid. Private preparation and public plotting are separate claims: successful plotting does not verify the construction of the supplied grid. See the scoped run evidence in `verification.md`.

## Editorial and coverage decisions still required

The [Econometric Society policy](https://www.econometricsociety.org/publications/es-data-editor-website/data-and-code-availability-policy), reviewed in the plan on 2026-10-22, requires an approved restricted-data arrangement; this inventory has no evidence of one. Direct policy retrieval on 2026-10-26 was blocked by the site's robots policy, so this document does not claim a fresh policy review. Authors must resolve exemptions/omissions and the journal's verification route before submission. No synthetic-data substitute has been approved or provided.

Authors must supply manuscript and appendix references, identify simulations/illustrations outside the current code, and approve exact distributable files. The author confirmed on 2026-10-26 that the three legacy figure filenames listed in `baseline.md` are not referenced by the paper sources. The [KSS archive](https://eml.berkeley.edu/~pkline/papers/KSS2020_COMP.zip), inspected for the plan, is a documentation precedent only; its historical URL/access and exemption are not present-day authorization. The [Social Science Data Editors README template](https://social-science-data-editors.github.io/template_README/) provides the intended final documentation structure; this inventory does not certify completed journal compliance.
