# Rank Multienrichment clusters

Rank Multienrichment clusters

## Usage

``` r
rank_mem_clusters(
  mem,
  clusters,
  choose = NULL,
  per_cluster = Inf,
  byCols = c("minp_rank", "composite_rank", "gene_count_rank"),
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `Mem` or legacy `list` mem output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- clusters:

  `list` containing set names, that must match `colnames(mem$memIM)` and
  `rownames(mem$enrichIM)`.

- choose:

  `character` vector with optional subset of clusters to return,
  matching `names(clusters)`; or `integer` vector referring to
  `clusters` by index position. When `choose` is `NULL`, all clusters
  are returned.

- per_cluster:

  `integer` vector with the number of entries to return for each
  cluster. Values will be recycled to the length of the clusters to be
  returned, defined by `choose` or by `length(clusters)` when `choose`
  is `NULL`.

- byCols:

  `character` vector used to sort the resulting `data.frame` within each
  cluster. This argument is passed directly to
  [`jamba::mixedSortDF()`](https://jmw86069.github.io/jamba/reference/mixedSortDF.html).
  Default 'minp_rank'. Recognized columns:

  - '**minp_rank**': the lowest P-value for all enrichments.

  - '**composite_rank**': a composite score which sorts by the order of
    magnitude of the lowest P-value, then by number of genes.

  - '**gene_count_rank**': the set with the most genes.

  Note that any column order can be reversed using "-" prefix, for
  example '-gene_count_rank' would return the set with the lowest gene
  count.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`data.frame` sorted by the criteria defined by `byCols`, with colname
`"set"` to indicate the pathway/set name. It includes additional columns
which may be useful in filtering or sorting.

- 'gene_count': the total number of genes involved in enrichment

- 'minp': the lowest enrichment P-value

- 'gene_count_rank': the rank by descending gene count per cluster.

- 'minp_rank': the rank by ascending minimum P-value.

- 'composite_rank': the rank of the composite score, ascending. The
  composite rank uses the order of magnitude of minimum P-value, then
  descending number of genes.

## Details

This function takes `list` output from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
and a `list` of clusters, and returns a `data.frame` that contains
several rank order metrics. It is intended to be used with column
clusters following
[`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
see examples.

The argument `per_cluster` is intended to make it convenient to pick the
top exemplar pathways, especially when argument `byCols` is defined so
that it sorts by the rank columns.

The argument `choose` is intended to make it easy to retrieve pathways
from specific clusters.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
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
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
## Start with mem
# mem <- multiEnrichMap(...);
# gp_hm <- mem_gene_path_heatmap(mem, column_split=4);
## Retrieve clusters from the Heatmap output, there should be 4 clusters
# clusters <- heatmap_column_order(gp_hm)
# clusters_df <- rank_mem_clusters(mem, clusters)
```
