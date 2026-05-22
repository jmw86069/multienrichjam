# Subset mem multienrichment object

Subset mem multienrichment object

## Usage

``` r
subset_mem(
  mem,
  includeSets = NULL,
  includeGenes = NULL,
  min_gene_ct = 1,
  min_set_ct = 1,
  p_cutoff = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `list` object as returned by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- includeSets:

  `character` vector with specific sets to retain, all other sets will
  be dropped.

- includeGenes:

  `character` vector with specific genes to retain, all other genes will
  be dropped.

- min_gene_ct:

  `numeric` filter applied to genes representing the minimum number of
  occurrences across sets in the `mem$memIM` incidence matrix. The
  default value `min_gene_ct=1` effectively requires a gene to be
  present in at least one set, which is useful after filtering by
  `includeSets`.

- min_set_ct:

  `numeric` filter applied to sets representing the minimum number of
  occurrences of genes in the `mem$memIM` incidence matrix. The default
  value `min_set_ct=1` effectively requires a set to contain at least
  one gene.

- p_cutoff:

  `numeric` optional enrichment P-value filter to apply to
  `mem$enrichIM` enrichment P-values. It is intended to apply optionally
  higher stringency by using a lower `p_cutoff` than used by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` object of `mem` data

## Details

This function is intended to subset the incidence matrix data contained
in a `mem` object by heuristics. It does not update other data in the
`mem` object such as `enrichList` and `multiEnrichDF`, nor any `igraph`
objects. It is intended mainly to subset by sets (pathways), or genes,
then also subset other corresponding incidence matrix data consistently.

## See also

Other jam utility functions:
[`ashape`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)`()`,
[`avg_angles`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md)`()`,
[`avg_colors_by_list`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)`()`,
[`bulk_cnet_adjustments`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md)`()`,
[`call_fn_ellipsis_deprecated`](https://jmw86069.github.io/multienrichjam/reference/call_fn_ellipsis_deprecated.md)`()`,
[`cell_fun_bivariate`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md)`()`,
[`collapse_mem_clusters`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)`()`,
[`colorRamp2D`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)`()`,
[`deconcat_df2`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md)`()`,
[`display_colorRamp2D`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md)`()`,
[`enrichList2geneHitList`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md)`()`,
[`filter_mem_genes`](https://jmw86069.github.io/multienrichjam/reference/filter_mem_genes.md)`()`,
[`filter_mem_sets`](https://jmw86069.github.io/multienrichjam/reference/filter_mem_sets.md)`()`,
[`find_colname`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)`()`,
[`find_enrich_colnames`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md)`()`,
[`get_hull_data`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md)`()`,
[`get_igraph_layout`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)`()`,
[`gsubs_remove`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md)`()`,
[`handle_igraph_param_list`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md)`()`,
[`isColorBlank`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md)`()`,
[`make_legend_bivariate`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md)`()`,
[`make_point_hull`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)`()`,
[`mem_find_overlap`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)`()`,
[`order_colors`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md)`()`,
[`rank_mem_clusters`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)`()`,
[`rotate_coordinates`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md)`()`,
[`subgraph_jam`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)`()`,
[`summarize_node_spacing`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md)`()`,
[`xyAngle`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)`()`
