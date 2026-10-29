# Data directories

`DATA_PATH` defaults to this directory. `raw_data` holds private source records
and converted VWH files; `processed_data` holds restricted clusters and panels.
`processed_data/driver_panel_grid.csv` is the sole approved public data file, with
its codebook in `processed_data/README.md`.

Data files remain ignored by Git. The archive explicitly includes the approved
grid and directory READMEs, without authorizing any other data. See
the root `README.md` for access, citations, rights and provenance.
Never copy the data directory wholesale into a public archive.
