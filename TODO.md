# TODO

This document describes plans for enhancements to the
multienrichjam R package.

## Apply min_count in `multiEnrichMap()`

* Currently `multiEnrichMap()` does not filter by number of
genes involved in enrichment, it only filters by enrichment P-value.
New argument `min_count` is applied only when `topEnrichN` is used,
but nothing else downstream is aware of filtering by `min_count`.
The corresponding argument in `mem_plot_folio()` is
`min_set_ct_each`, which requires a set (pathway) to contain
at least this many entries in at least one enrichment result
which also meets `p_cutoff` criteria for enrichment P-value.


## Optional highlight genes

* `mem_plot_folio()` and subsequent plots, optional argument
`highlight_genes` which would effectively hide all gene labels
except `highlight_genes` -- to help especially crowded plots.
Would propagate to new function `plot_cnet_heatmaps()` which
displays expression heatmaps around central cnet plot.

## add `plot_cnet_heatmaps()`

* `plot_cnet_heatmaps()` is in development, and arranges expression
heatmaps around a central Cnet plot, using genes in each cnet
cluster. It is intended with `collapse_mem_clusters()` to be
used with Cnet clusters.

## Improve reorderIgraphNotes() -- DONE needs more testing, verification

* `reorderIgraphNodes()` when it encounters attributes with multiple
colors per node, such as `"pie.colors"` and `"coloredrect.colors"`,
calls `avg_colors_by_list()` to generate one blended color per node,
then it sorts those colors using `sort_colors()`. However, the
color order ultimately does not match the order in the color legend,
and other plots such as heatmaps.
* Future idea is to convert node attributes with multiple colors per
node into a `data.frame` where each color is in a separate column.
Then convert each column to a factor whose levels match the color
order from `colorV` (argument to `multiEnrichMap()`). The end result
should sort nodes consistent with the order of colors.
* Note this process assumes that all node colors are from a limited
set, which ideally should match `colorV`. Therefore if any color
gradient is applied to nodes with `nodeType="gene"`, this process
will not work.
