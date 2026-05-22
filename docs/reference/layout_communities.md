# Layout igraph communities

Layout igraph communities, applies layout to each node community,
experimental.

## Usage

``` r
layout_communities(g, layout = layout_with_qfrf(repulse = 3), ...)
```

## Arguments

- g:

  `igraph` with graph attributes that include either 'mark.groups' or
  'nodegroups'.

- layout:

  `function` to apply to each igraph community.

- ...:

  additional arguments are ignored.

## Value

`igraph`

## Details

This function mimics
[`igraph::layout_components()`](https://r.igraph.org/reference/merge_coords.html),
in that it splits the input `igraph` into a `list`, applies layout to
each graph, then merges the layouts together.

This layout is experimental, as it currently uses methods in `igraph`
which assume each component is not connected, therefore the placement of
disjoint nodesets is not at all based upon any relationship to other
nodesets. In this case, communities may be connected to each other, and
so placement could usually be improved. This layout is an area for
future development.

## See also

Other jam igraph layouts:
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
