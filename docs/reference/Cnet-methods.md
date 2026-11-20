# Cnet-methods for MemPlotFolio

Methods to create Concept Network (Cnet) plots from a `MemPlotFolio`
object.

## Arguments

- x:

  `MemPlotFolio` object

- type:

  `character` used for
  [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- num:

  `integer` number of exemplar pathways, used with
  [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- cluster:

  `character` name of pathway cluster which matches an entry in
  `Clusters(x)`.

- do_plot:

  `logical` whether to render the igraph using
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

- legend_x, legend_y:

  passed to
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md).
  The default `legend_x='bottomleft'`. Use `legend_x=FALSE` to hide the
  legend.

- main:

  `character` with plot title.

  - For
    [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
    if `ClusterLabels(x)` is defined, it will use the corresponding
    label as plot title, otherwise if it uses `names(Clusters(x))`. To
    hide the title, use main=FALSE.

- ...:

  additional arguments may be recognized.

## Value

`igraph` Cnet object as defined for each method above.

## Details

The following methods are defined for `MemPlotFolio`:

- `CnetCollapsed(x, type, do_plot, legend_x, legend_y)`:

  Returns `igraph` Cnet object with pathways collapsed using
  `Clusters()`. Argument `type` controls the pathway cluster label:
  'title' uses the cluster title; 'set' uses the abbreviated set names;
  'cluster' uses
  [`ClusterLabels()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  if defined. Append **`'2'`** suffix to hide gene labels, for example
  `type='title2'` will display pathway cluster title, and hide gene
  labels.

- `CnetExemplar(x, num, do_plot, legend_x, legend_y, main)`:

  Returns `igraph` Cnet object with exemplar pathways for each cluster
  in `Clusters()`. Argument 'num' defines the number of exemplar
  pathways used for each cluster. It can be a single value, or vector of
  values applied in order of `Clusters()`.

- `CnetCluster(x, cluster, do_plot, legend_x, legend_y, main)`:

  Returns `igraph` Cnet object foran individual cluster defined in
  `Clusters()`. Argument 'cluster' defines the cluster, and must match
  an entry in `Clusters()` or be an `integer` index to the cluster.

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)
