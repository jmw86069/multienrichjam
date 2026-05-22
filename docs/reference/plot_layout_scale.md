# Plot layout scale by percentage of coordinate range

Plot layout scale by percentage of coordinate range

## Usage

``` r
plot_layout_scale(
  layout,
  grid_colors = c("grey80", "grey70"),
  grid_lty = c("dotted", "solid"),
  big_tick = 10,
  small_tick = 2.5,
  ...
)
```

## Arguments

- layout:

  `matrix` or `data.frame` with at least two columns, only the first two
  columns are used for the grid.

- grid_colors:

  `character` colors used for the small and big grid lines,
  respectively.

- grid_lty:

  `integer` or `character` line type used for the small and big grid
  lines, respectively.

- big_tick, small_tick:

  `numeric` values in percent, the step size between big grid lines, and
  small grid lines, respectively.

- ...:

  additional arguments are ignored.

## Details

The `layout` argument supplied coordinates, and largest numeric range of
any column is used to define 100 percent scale. A grey grid is drawn on
a base R plot to indicate the big and small steps across the range,
using `big_tick` and `small_tick`, respectively.

## See also

Other jam plot functions:
[`adjust_polygon_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md)
