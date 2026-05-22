# Convert communities object to nodegroups list format

Convert communities object to nodegroups list format

## Usage

``` r
communities2nodegroups(wc, sep = ",", ...)
```

## Arguments

- wc:

  `communities` object as returned by `igraph` functions such as
  `cluster_optimal()`, `cluster_walktrap()`, or
  `cluster_leading_eigen()`.

  - Alternatively, a `list` object with 'membership' and 'names'.

  - Alternatively, a `list` object which mimics the output of this
    function, in which case the input `wc` is returned as-is.

- sep:

  `character` string, default ',' (comma), used when cluster_names are
  defined as a `list` with vector of names for each cluster. Passed to
  [`jamba::cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html).

- ...:

  additional arguments are passed to
  [`jamba::cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html)
  when relevant.

## Value

`list` of `character` vectors, where each vector contains names of
`igraph` nodes. When `algorithm` is defined in the input object, it is
included as an attribute of the output `list`, accessible with
`attr(out, "algorithm")`.

When optional value `"cluster_names"` is present in the `communities`
object, they are used to define the output `list` names.

## Details

Note that this function is "lossy", in that the output `list` does not
contain all the information necessary to reconstitute the input
`communities` object in detail. However, the output `list` can be
converted to a `communities` object that will be accepted by most
`igraph` related functions that require that object type as an input
value.

Alternatively, this function can be used to confirm the input is already
appropriate nodegroups `list` format, by supplying `list` format
upfront.

## See also

Other jam igraph functions:
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
[`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md),
[`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md),
[`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
