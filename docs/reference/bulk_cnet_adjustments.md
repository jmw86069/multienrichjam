# Bulk Cnet plot adjustments

Bulk Cnet plot adjustments

## Usage

``` r
bulk_cnet_adjustments(
  g,
  nodeset_adj,
  node_adj,
  apply_negative = TRUE,
  nodesets = NULL,
  ...
)
```

## Arguments

- g:

  `igraph` object expected to be a Cnet object with vertex attribute
  'nodeType' containing values 'Gene' and 'Set'.

- nodeset_adj:

  `data.frame` recognizing five colnames:

  - nodeset: character name of each nodeset

  - x,y: numeric x,y coordinate adjustments to add to existing
    coordinates

  - rotate_degrees: numeric angle in degrees to rotate existing nodes

  - percent_spacing: absolute spacing to enforce for nodes in the
    nodeset.

- node_adj:

  `data.frame` recognizing five colnames:

  - node: `character` name of each node

  - x,y: `numeric` x,y coordinate adjustments to add to existing
    coordinates

  - label.dist: `numeric` value to add to existing label distance.

  - label.degree: `numeric` angle in degrees to add to the existing
    angle.

- apply_negative:

  `logical` whether to allow shrinking nodeset spacing.

- nodesets:

  `list` as output from `get_cnet_nodesets()`, or any list named by
  nodeset, containing `character` vector of node names.

- ...:

  additional arguments are ignored.

## Value

`igraph` object after applying adjustments.

## Details

This function applies a series of adjustments to nodesets and nodes in a
Cnet `igraph` object. Nodesets are adjusted first, then nodes are
adjusted afterward. This order is important because the x,y coordinates
may be altered by nodeset first, before the node adjustment is applied.

All x,y coordinate adjustments are performed using scaled units, using
the max x- or y-range for the plot layout. All Cnet `igraph` plots are
assumed to use `asp=1` to enforce aspect ratio 1, in order to maintain
the intended visual distances used by the layout algorithms.

Therefore it is possible to expand a plot layout such that the overall
scale is wider than before, thereby affecting the scaled unit space.

In any event, this function is intended to replicated the effect done by
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md).

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)
