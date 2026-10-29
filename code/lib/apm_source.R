# Fetch the pinned standalone source asset; shared by setup and archive assembly.
ensure_apm_source <- function(meta) {
  asset <- meta$source_tarball
  if (!identical(asset, 'vendor/apm_0.1.0.tar.gz')) stop('Unexpected APM asset path')
  matches <- function(path) identical(digest::digest(file = path, algo = 'sha256'), meta$source_sha256)
  if (file.exists(asset)) {
    if (!matches(asset)) stop('APM source checksum mismatch: ', asset)
    return(asset)
  }
  url <- meta$source_url
  if (is.null(url) || length(url) != 1L || !startsWith(url, 'https://')) {
    stop('Missing pinned HTTPS download URL in vendor/apm_source.json')
  }
  dir.create(dirname(asset), recursive = TRUE, showWarnings = FALSE)
  download <- tempfile('apm-download-', tmpdir = dirname(asset))
  on.exit(unlink(download), add = TRUE)
  download.file(url, download, mode = 'wb', quiet = TRUE)
  if (!matches(download)) stop('Downloaded APM source checksum mismatch; asset not installed')
  if (!file.rename(download, asset)) stop('Cannot save verified APM source: ', asset)
  asset
}
