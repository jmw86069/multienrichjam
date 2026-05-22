# Flip direction of igraph edges

Flip direction of igraph edges

## Usage

``` r
flip_edges(g, edge_idx, verbose = FALSE, ...)
```

## Arguments

- g:

  `igraph` object

- edge_idx:

  `integer` index of edges in the order they are stored in
  `igraph::E(g)`, or what igraph calls an "edge sequence" which is a
  character name for each node, defined as "node1\|node2". For example
  "D\|A" would define an edge from node name "D" to node name "A". When
  `verbose=TRUE` a summary table is printed out to show which edges were
  flipped.

- verbose:

  `logical` indicating whether to print verbose output. When
  `verbose=TRUE` a summary table is printed with these columns:

  - `edge_seq`: the input edge sequence, for example when `edge_idx` is
    provided as a `character` vector, the input vector is printed here.

  - `edge_seq_matched`: edge sequence that matched the `g` object. For
    example, when `edge_idx` input is a `character` vector, only the
    edges that match the `g` input are included here.

  - `edge_idx`: the integer index values of edges flipped. An `NA` value
    indicates the edge was not flipped, which should only happen when
    input `edge_idx` is provided as a `character` vector and some edges
    do not match the `g` input.

- ...:

  additional arguments are ignored.

## Details

This function simply flips the direction of igraph edges, keeping all
other node and edge attributes.

Note that this function will flip the order of nodes for each edge
defined by `edge_idx`, regardless whether the `igraph` itself is a
directed graph.

When `edge_idx` is provided as a `character` vector edge sequence, any
entries that do not match edges in `g` are ignored. A summary table is
printed when `verbose=TRUE`.

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
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

## Examples

``` r
am <- matrix(ncol=5, nrow=5, byrow=TRUE,
   data=c(0,0,0,0,0,
      1,0,0,0,0,
      1,0,0,0,0,
      1,0,0,0,0,
      1,0,0,0,0),
   dimnames=list(head(LETTERS, 5),
      head(LETTERS, 5)))
am;
#>   A B C D E
#> A 0 0 0 0 0
#> B 1 0 0 0 0
#> C 1 0 0 0 0
#> D 1 0 0 0 0
#> E 1 0 0 0 0
g1 <- igraph::graph_from_adjacency_matrix(am)
plot(g1);

g2 <- flip_edges(g1, 3:4);
plot(g2);

```
