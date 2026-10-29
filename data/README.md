# Data directories

`DATA_PATH` defaults to this directory. The `raw_data` directory contains private
source records and converted VWH files, while `processed_data` contains the
derived firm clusters and worker panels. Of these data files, only
`processed_data/driver_panel_grid.csv` is approved for public distribution. Its
codebook is provided in `processed_data/README.md`.

Git excludes the data files. The public archive includes the approved plotting
grid and the directory READMEs; permission to distribute the grid does not
extend to the other data files. The data directory should therefore not be
copied in full into a public archive. See the root `README.md` for data access
instructions, citations, rights, and provenance.
