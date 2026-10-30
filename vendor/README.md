# APM source package

`apm_0.1.0.tar.gz` contains the R source package used for this replication.
`apm_source.json` records its download URL, tarball path, SHA-256 checksum, and
upstream commit. The setup program downloads the file if it is absent, verifies
its checksum, and installs it. Installation does not require a separate copy
of the APM repository.
The package contains its MIT license notices and an internal source manifest.

The pinned APM commit is `a94bfba2f8865cce6d711fa2b0b9c9dc7ede39f8` in
[brad-ross/apm](https://github.com/brad-ross/apm). Version `0.1.0` alone does not
identify the source. Installation requires the R dependencies in `renv.lock`,
Make, a C++17 compiler, and BLAS/LAPACK/Fortran support described in the root README.

The source archive used for the replication is available from
[apm v0.1.0](https://github.com/brad-ross/apm/releases/tag/r-v0.1.0-replication).
The manifest fixes the download URL and checksum. Git excludes the locally
cached tarball, but the journal archive includes it so that installation does
not require another download of the APM source. Package dependencies still
require setup. The archive program also downloads the source if it is absent.
Both setup and archive assembly stop if the checksum of a cached or downloaded
file differs from the recorded checksum.

To construct a replacement from a selected revision of the APM repository,
run the following command in that repository:

```sh
Rscript r/tools/build_source_package.R \
  --output-dir /absolute/path/to/replication/vendor \
  --manifest /absolute/path/to/replication/vendor/apm_source.json
```

The build program combines the R package and C++ sources into an installable
source archive and preserves the license notices. Review the source revision
and checksum before replacing the recorded archive. Rebuilding the same source
can change the bytes in the tarball, so ordinary setup uses the verified archive
without rebuilding it.
