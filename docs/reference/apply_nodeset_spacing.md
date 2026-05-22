# Apply minimum node spacing for each Cnet node set

Apply minimum node spacing for each Cnet node set

## Usage

``` r
apply_nodeset_spacing(
  cnet,
  percent_spacing = 3,
  apply_negative = FALSE,
  metric = c("median", "minimum"),
  nodesets = NULL,
  tolerance = 1.05,
  cnet_nodesets = NULL,
  debug = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- cnet:

  `igraph` object that contains attribute `"nodeType"` and values
  `nodeType=c("Gene", "Set")`, consistent with use of
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  and
  [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md).
  The node layout is expected to be stored in
  `igraph::graph_attr(cnet, "layout")`.

- percent_spacing:

  `numeric` indicating the target spacing in relation to the maximum
  x-axis and y-axis span for the layout coordinates. For example
  `percent_spacing=4` will require that the median spacing between nodes
  in each nodeset is at least 4 percent of the larger of the x-axis and
  y-axis ranges for the layout of input `cnet`.

- apply_negative:

  `logical` indicating whether to shrink node spacing when the spacing
  is larger than the target `percent spacing`. By default
  `apply_negative=FALSE` which allows nodes to exceed the spacing, and
  only adjusts node spacing when a nodeset is below this threshold.

- metric:

  `character` string indicating the metric used to summarize nodeset
  spacing:

  - `"median"` - takes the median spacing between nodes in the nodeset,
    so the end result has at least half nodes with `percent_spacing`
    distance, "typical distance between nodes". This option is useful
    when there may be one node pair close together, but most nodes have
    consistent distance from their nearest node neighbor.

  - `"minimum"` - takes the minimum spacing between nodes in the
    nodeset, so the end result requires all nodes to have at least
    `percent_spacing` distance.

  - `function` - a `function` can be provided which takes a `numeric`
    vector of node distances, and produces a single `numeric` summary
    value, which will be scaled relative to the maximum xy axis range.
    This option would allow using
    [`mean()`](https://rdrr.io/r/base/mean.html), or a quantile function
    such as `function(x){quantile(x, 0.3)}`.

- nodesets:

  `character` with optional subset of nodesets to adjust, or the default
  `nodesets=NULL` will adjust all nodesets. Values can be specified
  either in form or nodeset names, or any one or more members of a
  nodeset.

- tolerance:

  `numeric` value used to apply a small tolerance to the distance
  threshold, for example the default `tolerance=1.05` requires nodes to
  be more than 5% lower than the actual threshold before adjusting
  nodes. This additional activation threshold is experimental and may
  not be necessary. It is mainly intended to reduce calculations for
  nodesets that are only a small rounding error from being acceptable,
  which can happen when the overall plot range is adjusted during rounds
  of small node coordinate adjustments like with this function,
  [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md),
  `apply_nodeset_spacing()`, or
  [`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md).

- cnet_nodesets:

  `list` optional to provide the known nodesets upfront, in order to
  avoid re-determining the nodesets.

- debug:

  `logical` indicating whether to print debug information without
  adjusting the nodeset spacing. This option is useful to see summary
  information about nodeset spacing before choosing the
  `percent_spacing` threshold.

- verbose:

  `logical` indicating whether to print verbose output. The output
  includes text printed with `debug=TRUE` except that nodes are adjusted
  when `debug=FALSE`.

- ...:

  additional arguments are ignored.

## Value

`igraph` Cnet object, whose nodes are adjusted according to the
arguments given.

## Details

This function automates the process of spacing nodesets in a Cnet plot
by calling
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
on all nodesets whose median nearest node distance is above or below a
threshold.

By default nodes must be at least `percent_spacing` distance apart,
relative to the maximum x-axis and y-axis layout coordinate range.
Roughly speaking, nodes are clearly visible with distance between nodes
of at least 3 percent of the coordinate range.

Nodes are not compressed down to this distance unless using the argument
`apply_negative=TRUE`.

Note that this function does not detect when two clusters overlap. For
other specific layout adjustments use
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md).
In future, this function may attempt to re-position nodesets relative to
each other, but that process tends to be iterative and is not guaranteed
to have a solution.

When a nodeset contains only one node, there is no calculated distance
so it is reported as `NA` when `verbose=TRUE` or `debug=TRUE`. However
for internal purposes the nodeset is treated as if its nodeset spacing
is exactly `percent_spacing` so that no subsequent adjustments are
performed on the nodeset.

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)
