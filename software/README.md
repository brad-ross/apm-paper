# APM source distribution

The replication archive uses `vendor/apm_0.1.0.tar.gz`, a standalone R
source package built from the separate [APM repository](https://github.com/brad-ross/apm).
`apm-source.json` records the full upstream Git commit, the packaging tool and
license hashes, and the exact tarball SHA-256. The package version alone does
not identify this build. No software release tag or release URL has been
published for this asset; those manifest fields remain `null`.

The baseline is `264907b0a2db57bf4c49285313b6eea4681f2b95`, including the
merged seed and install fixes and match attribution. Its numerical equivalence
to historical paper outputs has not been verified. Building or loading this
package does not rerun empirical analysis.

## Build in the APM project

Run its `tools/build_source_package.R` with R and the build-only packages
`digest` and `jsonlite`. The tool stages `r/` and places `core/src` under the
staged package's `src/core/`. It rewrites binding/header include paths and both
platforms' Makevars paths, and lists every staged C/C++ translation unit in
`OBJECTS`, including match attribution. It excludes tests, compiled objects,
shared libraries, and saved R objects. It includes the MIT notices naming Brad
Ross and Lihua Lei. Scientific algorithms are unchanged.

From the APM project directory, supply absolute destinations:

```sh
Rscript tools/build_source_package.R \
  --output-dir /absolute/path/to/replication/vendor \
  --manifest /absolute/path/to/replication/software/apm-source.json
```

The source tarball includes `inst/APM-SOURCE.json`. After installation,
`system.file("APM-SOURCE.json", package = "apm")` locates its provenance.
This includes SHA-256 hashes of every staged source file, the build tool and
licenses, and the tracked scientific-source patch relative to the recorded
upstream commit. The outer manifest separately hashes the finished tarball.
Local packaging files may be uncommitted: their hashes, rather than a claimed
software release commit, identify the exact packaging changes. Build output
must be frozen and reused; repeated `R CMD build` runs can change tarball bytes.

## Install and verify independently

Use R 4.3.0 for the captured environment, a C++17 compiler, R's development
headers, and BLAS/LAPACK/Fortran linking support. Package dependencies are
Rcpp, RcppArmadillo, RcppParallel, BH, data.table, R6, and xxhashlite (plus
R's bundled parallel package). RcppParallel supplies the threading interface;
the Unix Makevars retains the upstream system-oneTBB fallback and serial
fallback. Windows requires the upstream RcppParallel configuration and a
matching R build toolchain. The source package does not bundle third-party
libraries; their separately installed distributions supply their own licenses.

Verify the SHA-256 against `apm-source.json` before installation. Copy the
tarball to an empty temporary directory, extract there, inspect the source
paths and install into an empty temporary library. Run installation from that
directory so no build step can resolve a developer-relative sibling path:

```sh
R CMD INSTALL --library=/absolute/path/to/temporary-library apm
Rscript -e '.libPaths(c("/absolute/path/to/temporary-library", .libPaths())); library(apm); stopifnot(nzchar(system.file("APM-SOURCE.json", package="apm")))'
```

Record install/load results in `apm-source.json` only after those checks pass.
These are package checks, with no data loading or numerical analysis. Release
publication remains a separate author action; the bundled asset and a future
release asset must be byte-identical.
