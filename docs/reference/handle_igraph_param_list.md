# Handle igraph attribute parameter list

Handle igraph attribute parameter list, internal function for
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

## Usage

``` r
handle_igraph_param_list(
  x,
  attr,
  factor_l,
  i_values = NULL,
  attr_type = c("node", "vertex", "edge"),
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `igraph` object

- attr:

  `character` name of the attribute to update in `x`.

- factor_l:

  `list` or `numeric` vector or `function`:

  - `list` of `numeric` vectors where `names(factor_l)` correspond to
    attribute names, and the names of numeric vectors are attribute
    values The attribute names and attribute values are used to match
    relevant entities of type `attr_type`. For matching entities,
    attribute values are used as defined by attribute name `attr`, and
    are multiplied by the matching numeric value in `factor_l`.

  - `numeric` vector which is directly multiplied by `i_values` to
    produce an adjusted output vector `i_values`.

  - `function` which is used to modify `i_values` by calling
    `factor_l(i_values)` to produce adjusted output `i_values`.

- i_values:

  `vector` of attribute values that represent the current attribute
  values in `x` for the attribute `attr`.

- attr_type:

  `character` string indicating the type of entity being adjusted in
  `x`:

  - `"node"` or `"vertex"` refers to
    [`igraph::vertex_attr()`](https://r.igraph.org/reference/vertex_attr.html)

  - `"edge"` refers to
    [`igraph::edge_attr()`](https://r.igraph.org/reference/edge_attr.html)

- ...:

  additional arguments are ignored.

## Value

`vector` of attribute values representing `attr`.

## Details

This mechanism is intended to help update `igraph` attributes in bulk
operations by the attribute values associated with nodes or edges. Most
commonly, the argument `factor_l` is multiplied by numeric attributes to
scale attribute values, for example label font size, or node size.

For example:

    handle_igraph_param_list(x,
       attr="size",
       factor_l=list(nodeType=c(Gene=1, Set=2)),
       i_values=rep(1, igraph::vcount(x)),
       attr_type="node")

This function call will match node attribute `nodeType`, the size of
nodes with attribute value `nodeType="Set"` are multiplied `size * 2`,
`nodeType="Gene"` are multiplied `size * 1`.

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
[`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md),
[`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md),
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
