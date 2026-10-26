# Clean data and Figure 1a plotting grid

`driver_panel_grid.csv` is the approved numerical plotting input for Figure 1a.
It is supplied with the public replication archive, but remains ignored by Git.
Other clean data in this directory may be restricted and are not approved for
public distribution by this Figure 1a permission.

Public reproduction runs `code/driver_panel_results.R` from the repository root.
The renderer reads this CSV under configured `CLEAN_DATA_PATH` and writes
`panel_plot_drivers_times.pdf` under configured `FIGURES_PATH`. Missing or invalid
grids are package errors: restore the supplied file. Rendering never constructs
the grid, reads driver records, or requires a private preparation success record.

## Grid codebook

The CSV has exactly four numerical columns and one row per retained cohort:

| Column | Meaning and units |
| --- | --- |
| `xmin` | Left edge of the cohort's inclusive week interval, in week coordinates. |
| `xmax` | Right edge of the cohort's inclusive week interval, in week coordinates. |
| `ymin` | Lower edge of the cohort's anonymous consecutive driver display rows. |
| `ymax` | Upper edge of the cohort's anonymous consecutive driver display rows. |

Each rectangle is the union of its historical driver/week tiles. Week centers
are integer week IDs. Tile width follows historical `geom_tile()` defaults:
the smallest spacing between retained week centers, or 1 for a single center.
Any multiweek cohort makes this width 1. Thus `xmin = first_week - width/2` and
`xmax = treatment_week + width/2`, including singleton-week cohorts. Driver row
centers are the integers starting at 1; their tiles have height 1. Consequently
`ymin = first_display_row - 0.5`, `ymax = last_display_row + 0.5`, and
`ymax - ymin` is the number of drivers in the cohort (at least 100). Adjacent
rectangles share vertical boundaries; the first begins at `ymin = 0.5`.

Rows follow ascending historical cohort labels: entry and treatment week IDs
are separately padded on the left with zeros to at least two characters, then
joined by ` - ` and sorted lexicographically. No driver identifier or individual
record is exported. Anonymous display-row ordering within a cohort has no effect
on the occupied cells or the figure.

## Private preparation and provenance

Authorized preparation runs `code/construct_driver_panel_grid.R` separately.
It reads private `RAW_DATA_PATH/riders.csv`, using only `rider_id`,
`first_week_in_panel_dateweek_id`, and `began_active_week_id`. Records missing
any of these fields are omitted. Identical driver/cohort records count once;
conflicting cohorts for one driver and invalid week intervals fail explicitly.
Cohorts are entry/treatment-week pairs with at least 100 distinct drivers.
Preparation writes the geometry before a separate rendering process reads it.

The transformation and style come from
`eed1a22:code_old/visualize_panels.R:279-320`. Recovery retains the inclusive
`seq(first_week_id, treatment_week_id)` coverage that the historical code
accidentally removed with `select()` before plotting `week_id`. Compact
rectangles preserve that intended tile coverage without expanding or exporting
individual driver/week records. The figure retains the historical colors,
`theme_test()`, Week/Driver labels, zero axis expansion, and dimensions
3.0625 by 3.25 inches.

Public rendering verifies the figure from supplied plotting inputs. Verification
of the upstream construction requires authorized access to private driver data;
the public grid does not supply that access.
