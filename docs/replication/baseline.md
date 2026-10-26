# Reference-output baseline and coverage inventory

Inspection date: 2026-10-26. Scope: tracked source code, historical driver plotting source, configuration, and local output filenames. No restricted file content, credentials, empirical values, PDFs, table contents, or snippet contents were opened. No VWH conversion, clustering, estimation, reporting, or Slurm job was run. Figure 1a implementation and run verification are recorded separately in `verification.md`.

## Source identities and scientific settings

The inspected replication revision is `9c8cac6dbd7aad675fff62edc005359915c1775a`; the sibling APM revision is `264907b0a2db57bf4c49285313b6eea4681f2b95`. APM DESCRIPTION records version `0.1.0` and the merged seed/install fixes. These identify inspected source, not the build that generated historical outputs. The bundled APM source identity, checksum and installation evidence are recorded in `software/apm-source.json`. Historical output-generating replication/APM revisions, build options, platform, BLAS/toolchain, seeds, threads, data fingerprints, timestamps of successful runs, and numerical agreement remain unconfirmed.

At the inspected replication revision, `code/init_server.sh` specifies first/last year 1998/2001, maximum rank 2, minimum cohort size 75, k from 1 to 25, and chosen k=3. `code/env_config.R` derives year-bin width `ceiling((LAST_YEAR-FIRST_YEAR)/MAX_RANK)=2`. The shared `APM_SEED` defaults to 9001, supports an explicit environment override, and is passed to firm clustering and the 1,000-draw Bayesian bootstrap. This default must not be assigned retrospectively to old output. Clustering uses 100 as its third positional argument and `n_inits=10`; preserve the source call rather than guessing its historical interpretation. Thread count is configured; deterministic clustering fixes remove the former single-thread requirement, but no empirical thread-equivalence comparison was performed here.

Panel/estimation helpers use log10 weekly earnings. The evaluation source masks cohorts with at least 1,000 workers when identification permits, writes checkpoints after cohort/outcome pairs, and the report uses the 15 largest cohorts represented in the metrics. All five specifications are evaluated: TWFE and rank-1/rank-2 principal components with/without outcome fixed effects, using equal cohort weights for PC. Reports compare the four PC specifications against TWFE. Decomposition reports Province shares for APM, TWFE and their difference with 95% intervals. Confirm these specifications against the manuscript before freezing a release.

## Exhibit map and baseline availability

[exhibits.csv](exhibits.csv) records exact baseline paths, dynamic filename patterns, parameter coverage, producing stages/scripts, direct input dependencies, upstream stages, validation methods, access scope, and unresolved manuscript references. Paths use the settled author-confirmed upload-root mapping `output/result_summaries/{figures,tables,snippets}/` (relative to `OUTPUT_PATH` when configured elsewhere). Preserve the singular `result_summaries` and existing filenames.

The map contains:

- 9 current/recovered figure outputs: 8 VWH figures and Figure 1a. Rank-dependent families produce one file per rank in `2:MAX_RANK` (one rank at the baseline); within each such figure k spans 1–25. The two block-missingness plots each contain rank-1 and rank-2 curves. Province facets/cluster totals stay within the firms dot plot.
- 3 LaTeX tables.
- 61 in-text snippet files: 10 panel/configuration snippets, 36 evaluation snippets (4 specifications × 3 error statistics × 3 summary metrics), and 15 decomposition snippets (1 parameter × 3 specifications × 5 statistics). Rows expand every baseline combination; filename-pattern columns document the dynamic families. Fixed snippet/table filenames can be overwritten by another study configuration.
- 3 unmapped historical figure filenames and 3 explicit gaps covering final manuscript reconciliation, external simulations, and other computational illustrations/appendices.

The initial filename-only inspection found 75 files beneath `output/result_summaries`: 11 figures, 3 tables, 61 snippets. All 75 have rows. Figure 1a was absent at that inventory point and is included as the recovered target; subsequent generation is not a historical baseline. Presence means only a filename was seen. Historical files remain private and unchanged; none is certified as the authoritative manuscript reference or approved for redistribution.

The three extra historical filenames have no generating write site in the current tracked code:

- `figures/block_missingness_id_plot_start_year=1998_end_year=2001_year_cluster_size=2_min_cohort_size=75_rank=2.pdf`
- `figures/block_missingness_id_firm_weighted_plot_start_year=1998_end_year=2001_year_cluster_size=2_min_cohort_size=75_rank=2.pdf`
- `figures/cohort_size_dist_and_num_cohorts_by_k_plot_start_year=1998_end_year=2001_year_cluster_size=2_min_cohort_size=75.pdf`

The first two coexist with current `_k=3.pdf` names. No equivalence or replacement is inferred. Author confirmation or producing historical source is required if any are manuscript imports. Directory names under `output/old_sim_results/` include `bea_sims` and several `leave_one_out_sims*` variants. Names alone do not establish whether they are simulations, empirical resampling, manuscript exhibits, or obsolete experiments. No standalone simulation generator was found among the current tracked computational sources. This is a coverage gap, not a claim that the manuscript has no simulations.

## Evidence inspected

- All 12 tracked R files under `code/`, including every `ggsave`, `writeLines`, `write_result_snippet`, and Parquet write site, plus `code/init_server.sh` and `README.md`.
- Historical `eed1a22:code_old/visualize_panels.R:279-320` for the driver section; the deletion/recovery provenance and author decisions in the unchanged plan.
- Filename-only listing of `output/result_summaries/` and directory-only identification of historical simulation-output locations. No source data were opened.
- The author-approved implementation scope, output paths and rights/coverage decisions.
- Current FRDB access page on 2026-10-26; see [data-sources.md](data-sources.md) for source links and provenance limits.

No Overleaf source was requested or accessed. Except Figure 1a, exact manuscript numbers and panel/in-text references remain unconfirmed. The source inventory is complete for the inspected write sites; manuscript-wide coverage remains open.

## Author reference handoff and comparison protocol

Authors must identify the exact final manuscript revision and map each panel, table, appendix exhibit, and nontrivial in-text number to these rows, adding missing external code/outputs. For each authoritative reference, retain a private immutable copy outside the public artifact with its filename, checksum, generating revision/configuration/data version, known runtime/environment evidence, and provenance uncertainty. No private checksum list or logs were generated for this inventory. Do not replace old files during future comparisons.

When separately authorized, compare future VWH runs against those references using unrounded underlying values, declared tolerances, exact integer/count agreement where appropriate, stable cohort/outcome ordering, and formatted table/snippet agreement. Visual checks should compare geometry, scales, labels and dimensions; PDF byte equality is not required. Record seeded clustering changes and require an author decision for scientific differences. Do not invent a historical seed or claim the newer software reproduces historical estimates without evidence.

Runtime/resource evidence for conversion, clustering, panel reports, evaluation, evaluation reports, decomposition and decomposition reports is **unmeasured/unverified here**. Timing-print statements and scheduler allocations are not measured successful runtimes. Final submission readiness also requires confirmed data citations/rights, editorial exemption and verification arrangements, source coverage for any external material, and an approved archive manifest. Implementation verification remains limited to Figure 1a and non-empirical setup/manual checks; this protocol does not authorize rerunning VWH.
