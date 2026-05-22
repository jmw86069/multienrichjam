# MultiEnrichment Heatmap of enrichment P-values

MultiEnrichment Heatmap of enrichment P-values

## Usage

``` r
mem_enrichment_heatmap(
  mem,
  style = c("dotplot_inverted", "dotplot", "heatmap"),
  apply_direction = FALSE,
  p_cutoff = NULL,
  min_count = 1,
  p_floor = 1e-10,
  point_size_factor = 1,
  point_size_max = 8,
  point_size_min = 2,
  row_method = "euclidean",
  column_method = "euclidean",
  name = "-log10P",
  row_dend_reorder = TRUE,
  row_dend_width = grid::unit(30, "mm"),
  row_fontsize = NULL,
  row_cex = 1,
  row_split = NULL,
  row_gap = grid::unit(2, "mm"),
  cluster_rows = TRUE,
  column_fontsize = NULL,
  column_cex = 1,
  cluster_columns = FALSE,
  sets = NULL,
  color_by_column = FALSE,
  cex.axis = 1,
  lens = 3,
  cexCellnote = 1,
  column_title = NULL,
  row_names_max_width = grid::unit(300, "mm"),
  column_names_max_height = grid::unit(300, "mm"),
  heatmap_legend_param = NULL,
  hm_cell_size = NULL,
  legend_height = grid::unit(6, "cm"),
  legend_cex = 1,
  direction_cutoff = 0,
  gene_count_max = NULL,
  top_annotation = NULL,
  outline = TRUE,
  show_enrich = NULL,
  use_raster = FALSE,
  do_plot = TRUE,
  ...
)
```

## Arguments

- mem:

  `Mem` or legacy `list` mem created by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

- style:

  `character` string indicating the style of heatmap: `"heatmap"`
  produces a regular heatmap, shaded by `log10(Pvalue)`; `"dotplot"`
  produces a dotplot, where the dot size is proportional to the number
  of genes. See function description for details on how to include the
  point size legend beside the heatmap. The main benefit of using
  "dotplot" style is that it also indicates the relative number of genes
  involved in the enrichment.

- apply_direction:

  `logical`, default FALSE, whether to define a bivariate color scheme
  which uses `mem$enrichIMdirection` when defined. The color scheme is
  intended to indicate both the directional strength (usually with some
  type of z-score) and the statistical enrichment (usually with the
  enrichment P-value or adjusted P-value).

- p_cutoff:

  `numeric` value of the enrichment P-value cutoff, by default this
  value is obtained from `mem$p_cutoff` to be consistent with the
  original
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  analysis. P-values above `p_cutoff` are not colored, and are therefore
  white. This behavior is intended to indicate pathways with P-value
  above this threshold did not meet the threshold, instead of pathways
  with similar P-values displaying with similar color.

- min_count:

  `numeric` number of genes required for a pathway to be considered
  dysregulated.

- p_floor:

  `numeric` minimum P-value used for the color gradient. P-values below
  this floor are colored with the maximum color gradient. This value is
  intended to be used in cases where one enrichment P-value is very low
  (e.g. 1e-36) to prevent all other P-values from being colored pale
  red-white and not be noticeable.

- point_size_factor:

  `numeric` used to adjust the legend point size, since the heatmap
  point size is dependent upon the number of rows, the legend may
  require some manual adjustment to make sure the legend matches the
  heatmap.

- point_size_min, point_size_max:

  `numeric` values which define the minimum and maximum point sizes,
  effectively defining the range permitted when `style="dotplot"`.

- row_method:

  `character` string of the distance method to use for row and column
  clustering. The clustering is performed by
  [`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html).

- name:

  `character` value passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  used as a label above the heatmap color legend.

- row_dend_reorder:

  `logical` indicating whether to reorder the row dendrogram using the
  method described in
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
  The end result is minor reshuffling of leaf nodes on the dendrogram
  based upon mean signal in each row, which can sometimes look cleaner.

- row_fontsize, column_fontsize:

  optional `numeric` arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  to size row and column labels.

- cluster_columns:

  `logical` indicating whether to cluster heatmap columns, by default
  columns are not clustered.

- sets:

  `character` vector of sets (pathways) to include in the heatmap, all
  other sets will be excluded.

- color_by_column:

  `logical` indicating whether to colorize the heatmap using
  `mem$colorV` as defined for each comparison. This option is currently
  experimental, and produces a base R heatmap using
  [`jamba::imageByColors()`](https://jmw86069.github.io/jamba/reference/imageByColors.html).

- cex.axis:

  `numeric` adjustment for axis labels, passed to
  [`jamba::imageByColors()`](https://jmw86069.github.io/jamba/reference/imageByColors.html)
  only when `color_by_column=TRUE`.

- lens:

  `numeric` value used in color gradients, defining the extent the color
  gradient is enhanced in the mid-ranges (positive `lens`), or
  diminished in the mid-ranges (negative `lens`). There is no
  quantitative standard measure for color gradient changes, so this
  option is intended to help adjust and improve the visual perception of
  the current data.

- cexCellnote:

  `numeric` character expansion value used only when
  `color_by_column=TRUE`, used to adjust the P-value label size inside
  each heatmap cell.

- column_title:

  optional `character` string with title to display above the heatmap.

- row_names_max_width, column_names_max_height, heatmap_legend_param:

  arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  and provided here for convenience.

- hm_cell_size:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html) or `numeric`, default
  `NULL`, to define an optional fixed heatmap cell size, useful to
  define absolute square heatmap cells. When `numeric` it is interpreted
  as "mm" units. Note that the heatmap total height is determined by the
  number of cells, and the total row gaps defined by the number of row
  gaps with `row_split` multiplied by `row_gap`.

- legend_height:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html), default 6 cm (60
  mm), to define the absolute height of the color gradient in the color
  key. This value is only used when `heatmap_legend_param` is not
  defined.

- legend_cex:

  `numeric` default 1, used to scale the legend fontsize relative to the
  default fontsize `10`. This value is only used when
  `heatmap_legend_param` is not defined.

- top_annotation:

  `HeatmapAnnotation` as produced by
  [`ComplexHeatmap::HeatmapAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
  or `NULL`, used to display customized annotation at the top of the
  heatmap. The order of columns must match the order of columns in the
  data displayed in the heatmap.

- outline:

  `logical` default TRUE, whether to draw an outline for each heatmap
  cell. Note: The outline is not drawn for `style="dotplot"` which
  already adds lines through the middle of each cell, not the border of
  each cell.

- show_enrich:

  `numeric` default NULL, indicating which of the enrichment metrics to
  show as a label in each cell. When only one type is shown, there is no
  prefix, but for multiple types, a prefix is shown for each. The
  metrics in order include:

  1.  "-log10P"

  2.  "z-score"

  3.  "number of genes"

- use_raster:

  `logical` passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  indicating whether to rasterize the heatmap output, used when
  `style="heatmap"`. Rasterization is not relevant to dotplot output
  since dotplot is drawn using an individual circle in each heatmap
  cell.

- do_plot:

  `logical` indicating whether to display the plot with
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
  or
  [`jamba::imageByColors()`](https://jmw86069.github.io/jamba/reference/imageByColors.html)
  as relevant. The underlying data is returned invisibly.

- ...:

  additional arguments are passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  for customization.

## Details

Note: It is recommended to call
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
with `do_which=1` in order to utilize the gene-pathway content during
clustering, which is more effective at clustering similar pathways by
gene content. Otherwise pathways are clustered using only the
`-log10(p)` enrichment P-value.

This function is a lightweight wrapper to
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
intended to visualize the enrichment P-values from multiple enrichment
results. The P-value threshold is used to colorize every cell whose
P-value meets the threshold, while all other cells are therefore white.

### Drawing Dotplot with Point Legend

The `style` argument controls whether a heatmap or dotplot is created.

- `style="dotplot"`: each heatmap cell is not filled, and the color is
  drawn as a circle with size proportional to the number of genes
  involved in enrichment. A separate point legend is returned as an
  attribute of the heatmap object.

- `style="dotplot_inverted"`: each heatmap cell is filled, and a circle
  is drawn with size proportional to the number of genes involved in
  enrichment. A separate point legend is returned as an attribute of the
  heatmap object.

To draw the dotplot heatmap including the point legend, use this
command:

    ComplexHeatmap::draw(hm,
       annotation_legend_list=attr(hm, "annotation_legend_list"))

Generally, the clustering using the gene-pathway incidence matrix is
more effective at representing biologically-driven pathway clusters.

## See also

Other custom plot functions:
[`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
[`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md),
[`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
