# TODO

This document describes plans for enhancements to the
multienrichjam R package.

## Improve reorderIgraphNotes()

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
