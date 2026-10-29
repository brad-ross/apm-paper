# Pinned APM source package

`apm_0.1.0.tar.gz` is the standalone R source package used by this replication.
`apm-source.json` records its download URL, tarball path, SHA-256, and upstream commit.
Setup downloads this file if absent, verifies its checksum and installs it; no
sibling checkout is needed.
The package contains its MIT license notices and an internal source manifest.

The pinned APM commit is `a94bfba2f8865cce6d711fa2b0b9c9dc7ede39f8` in
[brad-ross/apm](https://github.com/brad-ross/apm). Version `0.1.0` alone does not
identify the source. Installation requires the R dependencies in `renv.lock`,
Make, a C++17 compiler, and BLAS/LAPACK/Fortran support described in the root README.

The exact tested asset is published at
[apm v0.1.0](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication).
The manifest pins its download URL and checksum. The tarball is an ignored local
cache, not a tracked Git file. Journal archives bundle it so APM source need not
be downloaded again; package dependencies still require setup. Archive assembly
also downloads the asset if it is absent. A mismatched cached or downloaded file
stops setup/packaging rather than silently replacing the pinned source.

To produce a replacement from an explicitly selected APM checkout, run in that
repository:

```sh
Rscript r/tools/build_source_package.R \
  --output-dir /absolute/path/to/replication/vendor \
  --manifest /absolute/path/to/replication/vendor/apm-source.json
```

The builder stages the R package and C++ sources as an installable source archive,
preserving license notices. Review the changed source revision and checksum before
replacing the pinned asset. Rebuilding the same source can change tarball bytes;
reuse the verified artifact rather than rebuilding it during ordinary setup.
