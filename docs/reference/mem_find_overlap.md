# Find recommended overlap threshold for EnrichMap, experimental

Find recommended overlap threshold for EnrichMap, experimental

## Usage

``` r
mem_find_overlap(
  mem,
  overlap_range = c(0.1, 0.99),
  max_cutoff = 0.4,
  adjust = -0.01,
  debug = FALSE,
  ...
)
```

## Arguments

- mem:

  `list` output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- overlap_range:

  `numeric` range of Jaccard overlap values, default `0.1, 0.99` using
  step `0.01`.

- max_cutoff:

  `numeric` value between 0 and 1, to define the maximum fraction of
  nodes in the largest connected component, compared to the total number
  of non-singlet nodes.

- adjust:

  `numeric` used to adjust the final overlap, default `-0.01` will use
  the overlap one step before the max O score.

- debug:

  `logical` indicating whether to return full debug data, which is used
  internally to determine the best overlap cutoff to use.

- ...:

  additional arguments are passed to
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

## Value

`numeric` value with recommended Jaccard overlap coefficient.

## Details

It implements a straightforward approach to determine a reasonable
Jaccard overlap threshold for Enrichment Map data, and is still very
much open to improvement after more experience using it on varied
datasets.

The premise is that two pathways that have Jaccard overlap above a
threshold are connected by a network "edge".

- With extremely low threshold, most pathways would be connected, even
  if they have only one gene in common.

- With an extremely high threshold, pathways would only be connected if
  nearly all genes were in common.

- A moderate threshold is intended to balance the two extremes.

- The aesthetic and biological interesting threshold appears to be
  dependent upon the type and number of pathways returned from
  enrichment analysis. For example, immunology pathways may favor a
  different threshold than metabolic pathways. (Purely hypothetical.)

- As a result, this function is intended to find a middle ground based
  upon the pathway data used for analysis at the time, where some but
  not all pathways are connected.

The method finds the overlap threshold at which the first connected
component is no more than `max_cutoff` fraction of the whole network.
This fraction is defined by the number of nodes in the largest connected
component, divided by the total number of non-singlet nodes.

We found that `max_cutoff=0.4`, the point at which the largest connected
component contains no more than 40% of all nodes, seems to be a
reasonably good threshold.

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
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)
