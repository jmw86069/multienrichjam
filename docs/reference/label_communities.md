# Assign labels to igraph communities

Assign labels to igraph communities

## Usage

``` r
label_communities(
  wc,
  labels = NULL,
  add_catchwords = NULL,
  num_keep_terms = 3,
  keep_terms_sep = ",\n",
  do_fixSetLabels = TRUE,
  ...
)
```

## Arguments

- wc:

  `communities` object, or `list` in form of nodegroups, which is a
  `list` of `character` vectors that contain `igraph` node names.

- labels:

  `character` vector of optional labels to assign directly to community
  clusters. When not defined, the auto-detection method is used.

- add_catchwords:

  `character` of optional words to include as catchwords, to be excluded
  from use in the final label.

- num_keep_terms:

  `integer` maximum number of terms to be included in the final output
  label, when auto-detection is used.

- keep_terms_sep:

  `character` string used as a delimited to separate each term when
  multiple terms are concatenated together to form the cluster label.

- ...:

  additional arguments are ignored.

- do_fix_terms:

  `logical` default TRUE, whether to apply
  [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  on the resulting words.

## Value

`communities` or `list` format matching the input `wc` format.

- When `communities` is input, additional value `cluster_names` will
  contain a `character` vector of names corresponding to each integer
  index in `wc$membership`.

- When `nodegroups` is input, the `list` names will be a `character`
  vector of cluster labels.

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
[`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md),
[`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
