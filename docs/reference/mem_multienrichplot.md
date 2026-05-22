# MultiEnrichMap plot

MultiEnrichMap plot, deprecated, use
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

## Usage

``` r
mem_multienrichplot(
  mem,
  overlap = 0.1,
  overlap_count = 2,
  do_plot = TRUE,
  do_legend = TRUE,
  remove_blanks = TRUE,
  remove_singlets = TRUE,
  spread_labels = TRUE,
  y_bias = 1,
  label_edges = c("overlap_count", "count", "overlap", "label", "none"),
  edge_cex = 1,
  node_cex = 0.8,
  node_size = 5,
  edge_color = "#55555588",
  frame_color = "#55555500",
  shape = "pie",
  repulse = 3.7,
  sets = NULL,
  rescale = TRUE,
  edge_bundling = "connections",
  main = "MultiEnrichMap\noverlap >= {overlap}, overlap_count >= {overlap_count}",
  ...
)
```

## Arguments

- mem:

  legacy `list` mem from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
  specifically containing elements `"multiEnrichMap","multiEnrichMap2"`
  which are expected to be `igraph` objects.

- overlap:

  numeric value between 0 and 1, indicating the Jaccard overlap filter
  to use for edge nodes. The value is used to delete edges whose values
  `E(g)$overlap` are below this threshold.

- overlap_count:

  numeric value indicating the minimum overlap count below which edges
  are removed. The value `E(g)$overlap_count` is used for this filter.

- do_plot:

  logical indicating whether to plot the final result.

- do_legend:

  logical indicating whether to print a color legend, valid when
  `do_plot=TRUE`. Arguments `...` are also passed to
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md),
  for example `x="bottomleft"` can be overriden.

- remove_blanks:

  logical indicating whether to call
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
  which removes blank/empty colors in igraph nodes.

- remove_singlets:

  logical indicating whether to remove singlet nodes, which are nodes
  that have no connections to other nodes. By default, singlets are
  removed, in order to help visualize the node connections that remain
  after filtering by `overlap`.

- spread_labels:

  logical indicating whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  which places node labels at an angle offset from the node, in order to
  improve default label positions.

- y_bias:

  numeric value passed to
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  when `spread_labels=TRUE`.

- repulse:

  numeric value used for network layout when
  [`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md)
  is used.

- sets:

  optional character vector of enriched sets to include, all other sets
  will be excluded. These values are matched with `V(g)$name`.

- rescale:

  logical indicating whether the `igraph` layout coordinates are scaled
  to range `c(-1, 1)` before plotting. In practice, when `rescale=FALSE`
  the function
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  is called because it does much better at properly handling other
  settings during the change. The effect is mainly to keep layout aspect
  intact, in cases where the x- and y-axis ranges are not roughly the
  same size, for example a short-wide layout.

- main:

  character string used as the title when `do_plot=TRUE`. This character
  string is passed to
  [`glue::glue()`](https://glue.tidyverse.org/reference/glue.html) in
  order to include certain argument values such as `overlap`.

- ...:

  additional arguments are passed to
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
  [`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
  and
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  as needed.

## Value

invisibly returns the `igraph` object used for plotting, a by-product of
this function when `do_plot=TRUE` is that the igraph object is also
visualized. All custom plot elements are updated in the `igraph` object,
so in principle a simple call to `plot(...)` should suffice.

## Details

This function is replaced by
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
and is deprecated.

This function takes output from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
and produces customized "multiple enrichMap" plots using an igraph
network. It differs from the data provided from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
mostly by enabling different overlap filters, and by automating several
steps that help with network layout, and node label placement.

For the most flexible exploration of data, run
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
using a lenient `overlapThreshold`, for example `overlapThreshold=0.1`.
Then call this function with increasing `overlap` until the
visualization has insightful structure.

## See also

Other jam deprecated functions:
[`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md),
[`enrichMapJam()`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md),
[`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md),
[`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
