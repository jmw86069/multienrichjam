# Relayout each nodegroup in a bipartite (Cnet) graph, experimental

Relayout each nodegroup in a bipartite (Cnet) graph, experimental

## Usage

``` r
relayout_nodegroups(
  cnet,
  nodegroups = NULL,
  repulse = 3.5,
  fix_set_nodes = TRUE,
  spread_labels = TRUE,
  add_edges = TRUE,
  edge_factor = 2,
  do_final_relayout = NULL,
  final_repulse = 3.5,
  apply_by_size = TRUE,
  byCols = c("-num_terms", "-num_nodes", "nodegroup"),
  verbose = FALSE,
  ...
)
```

## Arguments

- cnet:

  `igraph` object with node layout already defined

- nodegroups:

  `list` of node names, or `communities` object, passed to
  [`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md).
  When NULL, it checks for supporting data in this order:

  1.  If graph attribute 'mark.groups' is defined, it is used.

  2.  If vertex attribute 'nodeType' exists, it calls
      [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md).

  3.  Finally, it calls
      [`igraph::cluster_optimal()`](https://r.igraph.org/reference/cluster_optimal.html),
      hoping this method will be appropriate for the graph size.

- repulse:

  `numeric` default 3.5, passed to
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)

- fix_set_nodes:

  `logical` default TRUE, whether to fix all nodes with nodeType=='Set'
  to prevent them from moving.

- spread_labels:

  `logical` default TRUE, whether to apply
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  after the relayout iterations are complete.

- add_edges:

  `logical` default TRUE, whether to add edges within nodes of each
  nodegroup. The same is accomplished by setting `edge_factor=0`.

- edge_factor:

  `numeric` default 2, used as the numerator in new edge weights with
  equation `edge_factor/sqrt(n)` where 'n' is the number of nodes in the
  nodegroup.

  - Set `edge_factor=0` or `add_edges=FALSE` to skip this step.

- do_final_relayout:

  `logical` default NULL, whether to apply one more
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)
  after each nodegroup is adjusted. It uses `final_repulse`.

- final_repulse:

  `numeric` used when `do_final_relayout` is TRUE, used as the 'repulse'
  argument in
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md).

- apply_by_size:

  `logical` default TRUE, whether to apply the relayout to nodegroups
  ordered by size, using `byCols` to sort.

  - The default applies layout such that nodegroups with the most terms
    in `names(nodegroups)` are applied first, then largest to smallest
    nodegroups.

  - For Cnet plots, Gene nodes in nodegroups connected to the most Set
    nodes are adjusted first, then largest nodegroups, then sorted by
    nodegroup name.

  - It is unclear if the order is useful, future iterations of this
    approach may "move" other nodegroups aside first, then re-introduce
    each nodegroup into the layout one by one. Otherwise nodes could
    become "tangled" in the center, with no ideal method to optimize
    separation by nodegroup.

- byCols:

  `character` vector used when `apply_by_size` is TRUE.

  - '-num_terms': reverse-order by the number of terms in each nodegroup
    name, assuming comma-delimited terms.

  - '-num_nodes': reverse-order by the number of nodes in each
    nodegroup.

  - 'nodegroup': alphnumeric sort of the `names(nodegroups)`.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are passed to
  [`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md)
  for argument 'sep', or to
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  for arguments regarding node ordering, etc.

## Value

`igraph` object with updated layout.

## Details

This function iteratively re-applies a layout function to each nodegroup
in a graph, constraining the position of all other nodes for each
iteration. Currently the layout uses
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
in future it may use any layout function.

The purpose is to "encourage" nodes in a nodegroup to become bundled
together.

Strategy:

- Each nodegroup is isolated, and
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)
  is called on nodes in each nodegroup, while constraining all other
  nodes so they cannot move.

- In theory, using the same `repulse`, the nodes would not move at all.
  Changing the repulse force could encourage nodes to stay together.

- By default `add_edges=TRUE` which adds phantom edges to connect all
  nodes in a node group.

  - The edge weight is scaled down by the number of nodes.

  - These phantom edges are intended to help 'encourage' the nodegroup
    nodes to group together.

  - Otherwise, most layout algorithms are only focused on specific edge
    forces, and not secondary forces which are common in Cnet plots.

  - For example, nodes in a nodegroup all share the same network
    connections, however they are not otherwise attracted to each other
    in a network layout. In absence of any repulsive force, they would
    all be co-located. But with some repulsive force, they are repelled
    from each other, and sometimes end up radially positioned around the
    plot, and not grouped together.

  - The phantom edges add some minimal force for nodes to be grouped
    closer together, and are removed once the layout is complete.

- As a final polishing step, `do_final_relayout=TRUE` enables a final
  round of global node layout, with `final_repulse`.

  - We observed that sometimes the nodegroups are too clumped, in a big
    circular "ball" due to the phantom edge process above. The final
    relayout is helpful to allow nodes to space out somewhat.

  - It may be helpful to pass `niter` to control the number of layout
    iterations in this final step. The default is 500.

### Todo

Bonus points:

- Determine if each nodegroup is "split" around another nodegroup. If
  so, iterate that nodegroup using varying 'repulse' or other strategies
  to attempt to minimize the issue.

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md)

## Examples

``` r
cnet <- make_cnet_test();
ns <- get_cnet_nodeset(cnet)
# mark.groups: highlight just one nodegroup
# nodegroups: enables the edge_bundling to work properly
jam_igraph(cnet, mark.groups=ns["SetA,SetB"], nodegroups=ns,
   main="One nodegroup is split")


cnet2 <- relayout_nodegroups(cnet, nodegroups=ns["SetA,SetB"], do_final_layout=TRUE, verbose=TRUE)
#> ##  (11:03:42) 22May2026:   relayout_nodegroups(): Using nodegroups as supplied. 
#> ##  (11:03:42) 22May2026:   relayout_nodegroups(): Applying to nodegroup 'SetA,SetB' with members: BH,BT,BU 
jam_igraph(cnet2, mark.groups=ns["SetA,SetB"], nodegroups=ns,
   main="This nodegroup is 'encouraged' to re-group")


# by default, it is applied to all nodes
cnet3 <- relayout_nodegroups(cnet)
jam_igraph(cnet3, mark.groups=ns, nodegroups=ns,
   main="Re-layout across all nodegroups")

```
