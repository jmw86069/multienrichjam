# MultiEnrichMap color legend

MultiEnrichMap color legend

## Usage

``` r
mem_legend(
  mem,
  x = "bottomleft",
  y = NULL,
  bg = "#FFFFFF99",
  box.col = "transparent",
  title = "Color Key",
  cex = 0.8,
  ncol = 1,
  pch = 21,
  pt.cex = 2,
  pt.lwd = 2,
  inset = 0,
  do_directional = NULL,
  directional_column = c("same", "added-bottom", "added-top"),
  directional_colors = c(`up-regulated` = "firebrick3", `no change` = "grey80",
    `down-regulated` = "blue"),
  ...
)
```

## Arguments

- mem:

  `list` object output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
  specifically expected to contain element `"colorV"`.

- x, y, bg, box.col, title, cex, ncol, pch, pt.cex, pt.lwd, inset:

  arguments are passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html). Note
  `pt.lwd` is mostly relevant with `do_direction=TRUE`, which adds open
  circles to the legend, whose line width has default `pt.lwd=2`.

- do_directional:

  `logical` indicating whether to include directional colors defined in
  `directional_colors`, indicated only as the border color of nodes.

- directional_column:

  `character` indicating how to add the directional colors to columns of
  legend, with two options:

  - "same": appends `directional_colors` to legend colors using the
    defined `ncol` number of columns.

  - "added-bottom": appends `directional_colors` as a new column so the
    resulting legend with have `ncol+1` columns. In this case,
    intervening empty rows are filled with blank space, and the
    `directional_colors` are shown in the bottom-most rows in the far
    right column of the legend.

  - "added-top": appends `directional_colors` as a new column so the
    resulting legend with have `ncol+1` columns. In this case,
    intervening empty rows are filled with blank space, and the
    `directional_colors` are shown in the top-most rows in the far right
    column of the legend.

- directional_colors:

  `character` vector of R colors, named by the label to be shown in the
  legend, displayed in order (top to bottom, left to right) they appear
  in this vector. To remove the entry `"no change"`, supply a new
  vector:
  `directional_colors=c(up.regulated="firebrick3", down.regulated="blue")`

- ...:

  additional arguments are passed to
  [`legend()`](https://rdrr.io/r/graphics/legend.html).

## Details

This function is a simple wrapper around
[`legend()`](https://rdrr.io/r/graphics/legend.html) to add a color key
to a plot, typically for `igraph` plots.

## See also

Other custom plot functions:
[`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md),
[`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
[`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
