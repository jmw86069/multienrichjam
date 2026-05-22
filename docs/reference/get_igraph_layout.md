# Obtain or create layout for igraph object

Obtain or create layout for igraph object

Set the node layout for an igraph object

## Usage

``` r
get_igraph_layout(
  g,
  layout = NULL,
  default_layout = igraph::layout_nicely,
  verbose = FALSE,
  ...
)

set_igraph_layout(
  g,
  layout = NULL,
  default_layout = igraph::layout_nicely,
  verbose = FALSE,
  prefer = c("graph_attr"),
  spread_labels = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object

- layout:

  is always applied when not NULL, even when layout exists in `g`. Input
  should be one of:

  - `numeric` matrix of layout coordinates, with `nrow(layout)` equal to
    the number of nodes `igraph::vcount(g)`.

  - `function` that takes `igraph` input, and returns `numeric` matrix
    of layout coordinates.

  - `NULL`: default, uses the graph attribute 'layout' if it exists,
    `igraph::graph_attr(g, "layout")` if it exists. If it does not
    exist, it follows `default_layout`.

- default_layout:

  is only applied when `layout` is NULL, and no layout is defined in
  graph `g`. Input should be `function` or NULL:

  - Default
    [`igraph::layout_nicely()`](https://r.igraph.org/reference/layout_nicely.html)
    is used for consistency with igraph conventions. This algorithm
    should be safe for very large graphs, which may be extremely
    inefficient with
    [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
    for example.

  - When NULL, it will return NULL unless layout is defined in the input
    `g` graph.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to any layout function called.

- prefer:

  `character` vector with preferred method of storage:

  - 'graph_attr': store the `matrix` in graph attribute 'layout'.

  - 'vertex_attr': store coordinates in vertex attributes 'x', 'y', and
    optionally 'z' when defined.

- spread_labels:

  `logical` default FALSE, whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  to re-position node labels radially away from incoming edges.

  - Note that
    [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
    also by default `do_reorder=TRUE` which will re-order nodes by
    color, border, label, and name.

  - It is used to use TRUE when node labels were previously spread, so
    that the angle of offset is updated per the new layout coordinates.

## Value

`get_igraph_layout()` returns a `numeric` matrix when: the input graph
`g` contains layout as a graph attribute as either numeric matrix or
function, or coordinates as x,y,z vertex attributes. The layout will
contain colnames that begin 'x', 'y', 'z', for consistency with
downstream use.

- When there is no layout defined in `g` and `default_layout` is NULL,
  it returns NULL. This logic is intended when it is preferable to avoid
  returning a layout if it does not already exist.

- When `matrix` is returned, the number of rows matches the input graph
  `g` using `igraph::vcount(g)`, and in that order. All
  [`rownames()`](https://rdrr.io/r/base/colnames.html) are defined to
  match vertex name when it exists, using `igraph::V(g)$name`.

`set_igraph_layout()` returns `igraph` object with layout defined per
function arguments.

## Details

This function is a simple helper function intended to retrieve the node
layout for an `igraph` object.

The layout is defined with the following priority:

1.  When `layout` is supplied as an argument, it is used. When it is a
    `function` it is applied to `g` to produce numeric `matrix`;
    otherwise it should be a numeric `matrix` and is used as-is.

2.  When graph attribute 'layout' is defined, it is used as described
    for argument `layout` above, accepting either `function` or `matrix`
    values.

3.  When vertex attributes 'x' and 'y' are defined, optionally 'z',
    their values are used to produce a numeric `matrix`.

4.  When `default_layout` is supplied as a `function` it is applied to
    graph `g` to produce a numeric `matrix`.

5.  Finally, when `default_layout` is NULL, this function returns NULL.
    This fallback is intended only when it is desirable not to apply a
    new layout function, useful for large graphs.

### Additional rules

- When `layout` is defined as a matrix with rownames, the rownames are
  matched to vertex attribute 'name' if it exists, using
  `igraph::V(g)$name`. This step is intended to help ensure nodes the
  layout can be supplied in any order without regard to the order
  defined in graph `g`.

  - When the `layout` rownames do not match vertex names, this function
    will [`stop()`](https://rdrr.io/r/base/stop.html).

- When `layout` is defined as a `function`, or when any layout function
  is applied as relevant, it is expected to return a `numeric` `matrix`,
  or data which can be coerced to a `matrix` with
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) or
  `as(x, "matrix")`.

  - The matrix rownames are matched to vertex names as described above.

  - Note that `data.frame` rownames are only retained at this step when
    they were already explicitly defined before coersion to `matrix`.

### Suggested Usage

- Define `layout` as a function in order to force the use of that
  function to produce layout coordinates. This step would always ignore
  pre-existing layout coordinates in graph `g`.

- Define `layout` as NULL, and `default_layout` as a function, to use an
  existing layout stored in graph `g`, then to apply the default layout
  function only if no layout already existed in graph `g`.

- Define `layout` as NULL, and `default_layout` as NULL, to return an
  existing layout stored in graph `g`, otherwise to return NULL without
  applying any layout. This option would avoid computationally expensive
  layout for large graphs, for example.

This function is a simple wrapper to `get_igraph_layout()` which also
defines the resulting layout in the graph `g`.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md),
[`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md),
[`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md),
[`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md),
[`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
[`find_enrich_colnames()`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md),
[`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md),
[`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md),
[`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md),
[`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md),
[`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md),
[`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md),
[`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md),
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
g <- make_cnet_test(2, c(12, 5))
gl <- get_igraph_layout(g, verbose=TRUE)
#> ##  (11:02:40) 22May2026:   get_igraph_layout(): using graph_attr layout. 
jam_igraph(g)


# apply repulse=4
gl2 <- get_igraph_layout(g, layout=layout_with_qfrf(repulse=4), verbose=TRUE)
#> ##  (11:02:40) 22May2026:   get_igraph_layout(): creating layout using layout() 
jam_igraph(g, layout=gl2)


igraph::graph_attr(g, "layout") <- layout_with_qfrf(repulse=4)
gl3 <- get_igraph_layout(g, verbose=TRUE)
#> ##  (11:02:40) 22May2026:   get_igraph_layout(): using graph_attr layout. 
#> ##  (11:02:40) 22May2026:   get_igraph_layout(): applying graph_attr layout() 
identical(gl2, gl3)
#> [1] TRUE
#> TRUE

g2 <- set_igraph_layout(g, spread_labels=TRUE)
jam_igraph(g2)

```
