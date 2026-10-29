# Verification update — 2026-10-28

Two independent full runs of the frozen scientific pipeline reproduced exactly
the same evaluation estimates, decomposition, tables/snippets and rendered
figures. Generated panels and assignments matched after row alignment across
all 25 cluster-count settings. This establishes repeatability under that fixed
software/configuration, not equality to the historical paper outputs.

The current Make/shell workflow replaces the former hash-based R controller.
Temporary manual checks cover stage ordering, successful-marker reuse, force
propagation, failed stages, evaluation checkpoint restart/resume, private/public
grid handling, Slurm dependency IDs, side-effect-free planning, and setup failure
before submission. No permanent tests were added and no empirical stage was
rerun for this orchestration change. The code folders and archive allowlist now
use data/analysis/reporting/lib/jobs/slurm subdirectories. The reorganized source
passes R parsing and shell syntax checks. Temporary ZIP/tar archives built and
extracted successfully with the revised allowlist. Public Figure 1a rendered
from the extracted package, a repeat invocation reused its completion marker,
and explicit force rendered it again. No restricted inputs were used in those
archive checks. The driver grid now lives at
`data/processed_data/driver_panel_grid.csv`; its bytes match the archived grid.
The revised archive rendered successfully after extraction without creating a
`clean_data` folder. A temporary attempt to allowlist another processed-data
file was rejected as an unapproved data path.

The dated record below describes the earlier implementation and its checks;
its hash/registry/lock behavior is not part of the current workflow.

# Verification record — 2026-10-26

This record distinguishes executed checks from intended functionality. No VWH
conversion, firm clustering, estimation, or reporting was run or submitted. No
historical VWH output was regenerated. Numerical equivalence to the paper remains
unverified. No journal check, exemption, release, or deposit is claimed.

## Figure 1a

Private preparation and public rendering were executed with R 4.3.0 on Linux
x86_64 (Ubuntu 22.04.5 LTS; AMD EPYC 7763 processor), using ggplot2 4.0.1. The private source file's SHA-256 was compared before
and after preparation and was unchanged; the fingerprint itself is not published.
Preparation took 0.236 seconds and rendering 0.958 seconds in the measured R
session. A separate command-line render took 1.44 seconds including R startup,
with peak resident memory 137,036 KB (about 134 MiB). These are observations on
this host, not hardware minima or Slurm resource recommendations.

The saved grid was compared to an independent in-memory reconstruction of the
historical `seq(first_week_id, treatment_week_id)` driver/week tiles, correcting
only the historical dropped `week_id` column. All retained cohort rectangle
bounds matched exactly. Blank week fields are omitted, and reverse intervals
retain the historical inclusive sequence behavior. The PDF was visually
inspected for axes, colors and cohort geometry. Its dimensions are 3.0625 by
3.25 inches. No manuscript reference PDF was available for visual comparison.

Temporary synthetic checks outside the repository covered the 99/100 cohort
threshold, missing values, repeated records, ordering, singleton intervals,
saved-grid readback, invalid/missing public inputs and preparation failure.
Renderer checks did not use private data. Only the four geometry columns are
exported; the grid remains ignored by Git and is included explicitly in the
local archive.

## Software and execution

The installed R dependencies are pinned in the two lockfiles. Setup checks use
project-local libraries, including an automatic public-only dependency scope.
APM source and installation evidence is recorded in `vendor/apm-source.json`;
its upstream commit is not evidence of historical numerical equivalence. The
standalone installation used GNU g++ 11.4.0, the R 4.3.0 bundled BLAS and system
LAPACK 3.10.0. These describe the checked host, not portable minimum versions.

Runner verification uses temporary synthetic input/output and checkpoint records,
no-op stage substitutes, and fake scheduler commands. These checks establish
control flow and state handling, not successful empirical estimation or cluster
execution. Real Slurm submission and compute-node access have not been tested.
The manual fixture checks passed for complete-output reuse, missing/partial
outputs, changed seed/input/code/software/specification, matching checkpoint
resume, incompatible/forced fresh attempts, descendant invalidation, selective
force, scheduler errors, active-job reuse/conflicts, and writer-lock exclusion.
Additional synthetic checks passed for matching/mismatched compute-job allocation,
stale intent, partial submission recovery, accounting fallback for jobs absent
from the live queue, scheduler outages and corrupted submission records. Invalid
cases produced no stage execution, new submission or success record. Extracted
evaluation-loop fixtures checked successful completion and estimator, metrics,
empty-result and partial-specification failures: successful checkpoints survived,
and incomplete eligible results could not be certified complete.

The explicit private Make target, including FORCE=all, prepared the saved grid
and then rendered it in a separate R process. The public Make target rendered
from that grid and then reused its verified output.

A fresh archive extraction in a separate directory restored a new public project
library and rendered Figure 1a without raw data, private preparation records,
APM installation or a sibling checkout. Exact locked dependencies already
installed on this host were available to copy; this was not an offline or fresh
operating-system test. Repeated execution reused the figure, FORCE=all rerendered
from the grid, and removing the grid caused an error rather than preparation.
The supplied grid was unchanged by public execution. The extracted and working-copy
PDFs rendered to identical PNG images. The journal-format ZIP was also extracted,
all listed SHA-256 checksums passed, and its public Make target rendered successfully.

## Archive and coverage

Manual archive-builder checks rejected path traversal, unapproved raw data,
symlink files, symlink parent directories, and a missing README PDF. Unlisted
private-file and credential sentinels were excluded. Checks used an isolated
synthetic fixture outside the repository, including a project path with spaces.
No test scripts or fixtures are shipped.

`exhibits.csv` inventories source write sites and historical filenames; it does
not establish final manuscript numbering or certify historical results. See
`baseline.md` and `data-sources.md` for the remaining author confirmations.
