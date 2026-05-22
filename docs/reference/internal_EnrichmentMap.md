# EnrichmentMap internal function

EnrichmentMap internal function

## Usage

``` r
internal_EnrichmentMap(
  x,
  do_plot,
  legend_x = "bottomleft",
  legend_y = NULL,
  params = list(repulse = 3.5, width = 30, group = "default", mark.expand = 4, do_legend
    = TRUE),
  ...
)
```

## Arguments

- x:

  `Mem` or `MemPlotFolio` object.

- do_plot:

  `logical` whether to render the `igraph`, which is done using
  `jam_igraph(Emap, ...)`.

- legend_x, legend_y:

  passed to
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  to render a color legend for the plot.

- params:

  `list` of parameters specific to EnrichmentMap:

  - `'repulse'`: `numeric` default 3.5 to define the initial network
    layout.

  - `'width'`: `integer` character width to apply word-wrap to node
    labels, and to community or nodegroup labels.

  - `'group'`: `character` string indicating the type of node grouping
    to use:

    - `'default'`: will use, in order of preference: if supplied a
      `MemPlotFolio` object, it will use cluster labels if available,
      then cluster titles; if supplied a `Mem` object it will use
      `igraph` community detection and the corresponding top keywords.

    - `'clusters'`: requires `MemPlotFolio` input and uses cluster
      titles.

    - `'cluster_labels'`: requires `MemPlotFolio` input and cluster
      labels if available, otherwise cluster titles.

    - `'community'`: will use `igraph` community detection.

    - `'none'`: will not use node grouping.

  - `'mark.expand'`: `numeric` passed to
    [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
    to define node group expansion for the shaded hull around node
    groups, default 4. Units are percentage of plot dimensions.

  - `'do_legend'`: `logical` default TRUE, whether to render the color
    legend.

  - Other parameters will be added to the function arguments with
    default values when implemented.

- ...:

  additional arguments are passed through to
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
  and
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md).

## Value

`igraph` Enrichment Map network is returned. When `do_plot=TRUE` the
network is drawn using `jam_igraph(Emap, ...)`.
