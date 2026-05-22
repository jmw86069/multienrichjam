# Collapse Multienrichment clusters

Collapse Multienrichment clusters

## Usage

``` r
collapse_mem_clusters(
  mem,
  clusters,
  choose = NULL,
  per_cluster = Inf,
  byCols = c("composite_rank", "minp_rank", "gene_count_rank"),
  return_type = c("cnet", "Mem", "mem"),
  max_labels = 4,
  max_nchar_labels = 25,
  include_cluster_title = TRUE,
  cluster_color_min_fraction = 0.5,
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `Mem` or legacy `list` mem object

- clusters:

  `list` named by cluster name, containing `character` vectors all of
  whose values should be present in `sets(mem)`.

- choose:

  `character` default NULL, passed to
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
  to choose a subset of clusters to return. Values must be found in
  `names(clusters)`.

- per_cluster:

  `integer` default 'Inf' with the maximum sets to include in each
  cluster. Use an integer like `per_cluster=1` to choose one "exemplar"
  per cluster, ranked according to `byCols` which is described in
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md).

- byCols:

  `character` vector of summary values to use when sorting sets within
  each cluster, passed to and described in
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md).

- return_type:

  `character` default "cnet" with the output format:

  - "cnet" returns an `igraph` Cnet plot.

  - "Mem" returns a `Mem` object after collapsing by cluster.

  - "mem" returns legacy `list` mem output.

- max_labels:

  `integer` default 4, the number of pathway labels to include in each
  cluster, ordered according to 'byCols'.

- max_nchar_labels:

  `integer` default 25, max character length for each pathway label.

- include_cluster_title:

  `logical` default TRUE, whether to include the cluster title as a
  prefix. For example if the cluster name is 'A' the new title would be:
  'A: pathway 1; pathway 2; pathway 3; pathway 4'.

- cluster_color_min_fraction:

  `numeric` value default 0.5, the fraction of pathways in a cluster
  which has "non blank" assigned color in
  [`enrichIMcolors()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  in order for the color to be assigned to the resulting node, in the
  case of Cnet output with `return_type="cnet"`.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are passed to `cnet2mem()` which includes
  thresholds used by
  [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
  such as `direction_col_fn` (color function applied to
  enrichIMdirection), `gene_direction_col_fn` (color function applied to
  geneIMdirection), and others. To avoid applying directional color to
  Set nodes (nodeType='Set') use `direction_cutoff=100` and
  `direction_max=101`, which will consider any enrichIMdirection value
  at or below 100 to be non-directional.

## Value

By default `return_type="cnet"` and this function returns an augmented
Cnet plot. The labels of each cluster are defined by the input
`names(clusters)`, however an `igraph` attribute `"set_labels"` includes
an abbreviated label of the top ranked sets for each cluster. This label
is probably the closest thing to summarizing the composition of each
cluster.

## Details

This function is similar to
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
in that it starts with `mem` results from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
and a list of `clusters` of pathways/sets. Instead of ranking and
choosing exemplars from each clusters, it simply collapses each cluster
into one super-set that contains union of all genes.

It does also run
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
in the event one would want to collapse only the top `per_cluster`
entries for each cluster, but the default is to include them all.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
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
