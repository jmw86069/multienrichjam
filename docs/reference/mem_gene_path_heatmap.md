# MultiEnrichment Heatmap of Genes and Pathways

MultiEnrichment Heatmap of Genes and Pathways

## Usage

``` r
mem_gene_path_heatmap(
  mem,
  genes = NULL,
  sets = NULL,
  min_gene_ct = 1,
  min_set_ct = 1,
  min_set_ct_each = 4,
  column_fontsize = NULL,
  column_cex = 1,
  row_fontsize = NULL,
  row_cex = 1,
  row_method = "binary",
  column_method = "binary",
  enrich_im_weight = 0.3,
  gene_im_weight = 0.5,
  gene_annotations = c("im", "direction", "default"),
  annotation_suffix = c(im = "hit", direction = "dir"),
  simple_anno_size = grid::unit(5, "mm"),
  cluster_columns = NULL,
  cluster_rows = NULL,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  name = NULL,
  p_cutoff = NULL,
  p_floor = 1e-10,
  row_split = NULL,
  column_split = NULL,
  auto_split = TRUE,
  column_title = LETTERS,
  row_title = letters,
  row_title_rot = 0,
  colorize_by_gene = TRUE,
  na_col = "white",
  rotate_heatmap = FALSE,
  colramp = "Reds",
  column_names_max_height = grid::unit(8, "cm"),
  column_names_rot = 90,
  show_gene_legend = FALSE,
  show_pathway_legend = TRUE,
  show_heatmap_legend = 8,
  use_raster = FALSE,
  seed = 123,
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `Mem` or `list` object created by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

- genes:

  `character` vector of genes to include in the heatmap, all other genes
  will be excluded.

- sets:

  `character` vector of sets (pathways) to include in the heatmap, all
  other sets will be excluded.

- min_gene_ct:

  `integer` minimum number of occurrences of each gene across the
  pathways, all other genes are excluded.

- min_set_ct:

  `integer` minimum number of genes required for each set, all other
  sets are excluded.

- min_set_ct_each:

  `integer` minimum number of genes required for each set, required for
  at least one enrichment test.

- column_fontsize, row_fontsize:

  `numeric` passed as `fontsize` to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  to define a specific fontsize for column and row labels. When `NULL`
  the nrow/ncol of the heatmap are used to infer a reasonable starting
  point fontsize, which can be adjusted with `column_cex` and `row_cex`.

- row_method, column_method:

  `character` string with distance method, default is 'binary' which
  considers all non-zero values to be 1. It is used for row and column
  hierarchical clustering by
  [`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html). It
  offers more methods than
  [`hclust()`](https://rdrr.io/r/stats/hclust.html) hence its use here.

- enrich_im_weight:

  `numeric` value between 0 and 1 (default 0.3), the relative weight of
  enrichment `-log10 P-value` and overall gene-pathway incidence matrix
  when clustering pathways.

  - When `enrich_im_weight=0` then only the gene-pathway incidence
    matrix is used for pathway clustering.

  - When `enrich_im_weight=1` then only the pathway significance
    (`-log10 P-value`) is used for pathway clustering.

  - The default `enrich_im_weight=0.3` balances the combination of the
    enrichment P-value matrix, with the gene-pathway incidence matrix.

- gene_im_weight:

  `numeric` value between 0 and 1 (default 0.5), the relative weight of
  the `mem$geneIM` gene incidence matrix, and overall gene-pathway
  incidence matrix when clustering genes.

  - When `gene_im_weight=0` then only the gene-pathway incidence matrix
    is used for gene clustering.

  - When `gene_im_weight=1` then only the gene incidence matrix
    (`mem$geneIM`) is used for gene clustering.

  - The default `_im_weight=0.5` balances the gene incidence matrix with
    the gene-pathway incidence matrix, giving each matrix equal weight
    (since values are typically all `(0, 1)`.

- gene_annotations:

  `character` string indicating which annotation(s) to display alongside
  the gene axis of the heatmap.

  - Default is `"im", "direction"`, and `"default"` which will hide the
    `"direction"` if all non-zero values are positive.

  - When `"im"` is present, the colored incidence matrix is displayed.

  - When `"direction"` is present, the directional matrix is displayed
    using colors defined by `colorjam::col_div_xf(1.2)`.

  - When both `"im"` and `"direction"` are present, they are displayed
    in the order defined.

  - When no values are given, the gene annotation is not displayed.

- annotation_suffix:

  `character` vector named by values permitted by `gene_annotations`,
  with optional suffix to add to the annotation labels. For example, as
  by default, it may be helpful to add "hit" or "dir" to distinguish the
  enrichment labels.

- name:

  `character` value passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  used as a label above the heatmap color legend. Default `NULL` uses
  "Gene Hits by Enrichment".

- p_cutoff:

  `numeric` value of the enrichment P-value cutoff, above which P-values
  are not colored, and are therefore white. The enrichment P-values are
  displayed as an annotated heatmap at the top of the main heatmap. Any
  cell that has a color meets at least the minimum P-value threshold.
  The default is taken from input `mem` for consistency with the input
  multienrichment analysis.

- column_split, row_split:

  row and column split, default NULL, detects an appropriate number of
  clusters. Passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

  - To turn off split use `1`.

  - To specify a fixed number of clusters, use an `integer` value. The
    value may be changed if the underlying data does not support that
    number of clusters.

  - To specify fixed clusters by name, use an atomic vector named by the
    rownames `genes(Mem)` or colnames `sets(Mem)` with values which will
    become the cluster names. When supplying a `factor` the factor level
    order will be maintained.

  - Alternatively, supply a `data.frame` whose rownames match the
    `genes(Mem)` rows or `sets(Mem)` columns, respectively. In this
    case, column values are combined using
    [`jamba::pasteByRowOrdered()`](https://jmw86069.github.io/jamba/reference/pasteByRowOrdered.html)
    which also maintains factor level order.

  - When supplying either an atomic vector or `data.frame` the actual
    names will be used to subset the resulting `Mem` data if there are
    fewer names provided than exist in `Mem`. Note that the `Mem` data
    are not subset again by other filter criteria, for example
    `min_set_ct`, `min_gene_ct`, etc.

  - The `column_title` argument is used when

- auto_split:

  `logical` whether to determine clusters when column_split or row_split
  is NULL, default TRUE.

- column_title:

  optional `character` string or vector to display above the column
  splits. Default uses `LETTERS` with length to match the number of
  clusters. When there is only one column cluster, it is not named
  unless `column_title` also has length 1.

- row_title:

  optional `character` string or vector to display to the left side of
  the heatmap row splits. The default uses `letters` to match the number
  of row clusters. When there is only one row cluster, it is not named
  unless `row_title` also has length 1.

- row_title_rot:

  `numeric` value, default 0, the rotation of `row_title` text, where
  `0` is not rotated, and `90` is rotated 90 degrees.

- colorize_by_gene:

  `logical` default TRUE, whether to color the main heatmap body using
  the colors from `geneIM` to indicate each enrichment in which a given
  gene is involved. Colors are blended using
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html),
  using colors from `colorV` in the `geneIMcolors(Mem)`.

- na_col:

  `character` string indicating the color to use for NA or missing
  values. Typically this argument is only used when
  `colorize_by_gene=TRUE`, where entries with no color are recognized as
  `NA` by
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- rotate_heatmap:

  `logical` indicating whether the entire heatmap should be rotated so
  that pathway names are displayed as rows, and genes as columns. When
  enabled, arguments referring to columns and rows are flipped, so
  "column" arguments will continue to affect pathways/sets, and "row"
  arguments will continue to affect genes. This includes `column_method`
  and `row_method` as of 0.0.90.900.

  - Exceptions:

  - `row_title_rot` is only applied to rows, due to its purpose.

  - `column_names_rot` is only applied to columns, also due to its
    purpose.

- colramp:

  `character`, default "Reds", with name of color, color gradient, or a
  vector of colors, anything that can be converted to a color gradient
  by
  [`jamba::getColorRamp()`](https://jmw86069.github.io/jamba/reference/getColorRamp.html).

- column_names_max_height:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html) passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
  When supplied as `numeric` it is converted to units in "mm". Default
  180 mm.

- column_names_rot:

  `numeric` passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- show_gene_legend, show_pathway_legend:

  `logical`, whether to show the gene IM and pathway IM legends,
  respectively.

  - The gene IM legend is FALSE by default, since it only describes the
    color used for each column, and is somewhat redundant with the
    pathway IM legend.

  - The pathway IM default is TRUE, it displays the color scale
    including the range of enrichment P-values colorized.

- show_heatmap_legend:

  `numeric` or `logical`, (default 8), the maximum number of labels to
  use for the heatmap color legend. When 'colorize_by_gene' is TRUE, the
  heatmap legend would include all possible blended colors using gene IM
  data, which frankly can become too much.

  - When `logical`, `TRUE` is converted to `8` by default.

  - When there are more legend items than than `show_heatmap_legend` the
    color legend will only display singlet colors, which means only one
    color per individual set defined in colorV.

  - This legend can be created and extracted from the output `Heatmap`
    object to be displayed independently for publishable figures if
    necessary.

- use_raster:

  `logical` passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  default FALSE, whether to rasterize the heatmap body. This option is
  recommended FALSE when 'colorize_by_gene' is TRUE, due to the way the
  rasterization is handled at matrix level. For very large heatmaps you
  may try 'colorize_by_gene=FALSE' and 'use_raster=TRUE' to reduce the
  figure size with PDF and SVG output, and will improve visual fidelity
  in some cases.

- seed:

  `numeric` value passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) to define a
  reproducible random seed during row and column clustering.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  for customization. However, if `...` causes an error, the same
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function is called without `...`, which is intended to allow
  overloading `...` for different functions.

## Value

`Heatmap` object defined in
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
with custom attributes with the method caption:

- `"caption"`: `character` string with method details.

- `"caption_legendlist"`:
  [`ComplexHeatmap::Legends`](https://rdrr.io/pkg/ComplexHeatmap/man/Legends.html)
  object suitable to be included with Heatmap legends by
  `draw(hm, annotation_legend_list=caption_legendlist)`, or or drawn
  directly `grid::grid.draw(caption_legendlist)`.

- `"draw_caption"` - `function` that will draw the caption in the
  bottom-right corner of the device by default, to be called with
  `attr(hm, "draw_caption")()` or `draw_caption()`.

The `Heatmap` row and column order can be retrieved:

1.  [`jamba::heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html) -
    returns a `list` of vectors of rownames in the order they appear in
    the heatmap, with list names defined by row split.

2.  [`jamba::heatmap_column_order()`](https://jmw86069.github.io/jamba/reference/heatmap_column_order.html) -
    returns a `list` of vectors of colnames in the order they appear in
    the heatmap, with list names defined by row split.

## Details

This function takes the `Mem` output from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
and creates a gene-by-pathway incidence matrix heatmap, using
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
The major output of this function is to define pathway clusters which
influences other figures produced by
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md):
the enrichment heatmap; Cnet cluster plots; Cnet exemplars.

It uses three basic sources of data to annotate the heatmap:

1.  `memIM` the gene-set incidence matrix

2.  `geneIM` the gene incidence matrix by dataset

3.  `enrichIM` the pathway enrichment P-value matrix by dataset

### Pathway and gene clusters

When `column_split` (pathways) or `row_split` (genes) are not defined, a
reasonable number is assigned to split columns and rows, respectively.
The specific number is controlled by defining an integer number. The
splits are named using `column_title` and `row_title`, which defaults
`LETTERS` and `letters`, respectively.

Custom pathway clusters (columns) can be defined by defining
`column_split` as a vector named by values in `sets(Mem)`, with values
used to define the name for each cluster. When provided as a `factor` it
will honor the factor level order.

When a subset of `sets(Mem)` are provided, it will also subset the `Mem`
object to show only the matching sets, however it does not re-apply row
and column filtering described above. It is recommended to use argument
'sets' to subset the pathway gene sets upfront, then use `column_split`
with names that match 'sets'.

Alternatively, `column_split` can be supplied as a `data.frame` with
rownames that match `sets(Mem)`, in which case column values are
combined using
[`jamba::pasteByRowOrdered()`](https://jmw86069.github.io/jamba/reference/pasteByRowOrdered.html)
which also keeps factor level order.

Custom gene clusters (rows) can be defined with `row_split` using the
same mechanism as with `column_split` described above, with names that
match `genes(Mem)` as appropriate.

#### Column pathway clustering

Columns (pathways) are clustered using a combination of the gene-pathway
incidence matrix, and the -log10 P-values shown along the pathway axis
(usually columns). The purpose is to allow the enrichment P-value to
influence the clustering together with the gene content. The balance is
adjusted using `enrich_im_weight`, default 0.3, and where 0.0 will
ignore the enrichment P-values during clustering. Higher values will
tend to create clusters that represent shared/unique *significant*
pathways, and less determined based solely on the gene content.

Note that the enrichment P-value matrix used during clustering is
adjusted by using `p_cutoff` and `p_floor`. Values above the `p_cutoff`
are converted to 1 so they do not influence clustering. Values below
`p_floor` are converted to `p_floor` so they influence clustering only
at the level of other values at `p_floor`.

The default `cluster_columns=TRUE` will employ
[`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html) as a
convenient and efficient one-step distance-hierarchical clustering
approach, using `cluster_method` as the distance method. A custom
function can be supplied with `cluster_columns` as long as the output is
'hclust' or can be coerced to 'hclust'.

The clustering distance method `column_method` uses default 'binary',
which treats any non-zero value as 1. In this case, the relative weight
is not important since all non-zero values are equivalent. However,
`column_method='euclidean'`, current default in
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
may improve the output when adjusting the relative weight of enrichment
P-values with incidence matrix.

#### Gene clustering

Rows (genes) are clustered using a combination of the gene-pathway
incidence matrix, and the 'geneIM' or 'geneIMdirection' matrix data
shown, as defined with `gene_annotations`. The relative weight of these
matrices is controlled with `gene_im_weight` with default 0.5, which
gives equal weight to each matrix. By default, when directional gene
values are shown the directional matrix is used with clustering.

The default `cluster_rows=TRUE` will employ
[`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html) as
described for Pathway clustering, or `cluster_rows` can be a custom
function that produces 'hclust' or can be coerced to 'hclust'.

When `row_method='binary'` the relative weight of gene incidence matrix
and the pathway-gene matrix is not important, since all non-zero values
are treated as 1 during clustering. The
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
default uses 'euclidean' to improve the effect of the relative weight of
these two matrices.

Gene clusters are not often used for downstream analysis, for example
they do not form clusters in the Cnet plots, and are not (yet) used for
other analysis in multienrichjam. However, gene clusters are quite
useful when interpreting pathway-gene data.

It is helpful in practice to refer to a gene cluster by name: "The genes
in cluster 'd' all appear to be cytokines." (Also maybe there should be
better names, but that's for another time.)

Gene clusters may form what we call "hot spots", where most of a gene
cluster is colorized and associated with one or more pathway clusters. A
hot spot indicates a set of genes shared across multiple pathways, a
"core" set of genes which may have serve an important functional basis
in several pathways.

An example might be mitogen-activated protein kinase (MAPK) genes, which
typically involve a multi-step kinase call signaling cascade. Pathways
which involve one MAPK would very often involve each MAPK at subsequent
steps in the cascade. In fact, MAPK genes are often the core signaling
mechanism of numerous apparently unrelated pathways - we refer to it as
the "internal wiring" of the signaling in a cell, as an analogy to a
building which may use wiring to pass along any number of messages.
Different cell types may employ the MAPK cascade to send a message, and
this message may have different meaning across cell types, and indeed
may have meaning based upon the cell state.

### Filtering

A subset of pathways can be defined with argument `sets`, which may be
useful to prepare this plot for a set of "exemplar pathways" which
represent key pathways of interest for a study.

Similar for genes with argument `genes`, however this option is less
commonly used.

Pathways and genes can be subset by number of occurrences of each, for
example:

- `min_gene_ct`: minimum occurrences of a gene across pathways, which
  also means the number of pathways in which a gene occurs.

- `min_set_ct`: minimum occurrences of a set (pathway) across genes,
  which means the number of genes in the set across all enrichments.

- `min_set_ct_each`: minimum occurrences of a set (pathway) across
  genes. It requires at least `min_set_ct_each` genes in at least one
  enrichment, which is consistent with `min_count` in
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

When pathways are filtered by `min_gene_ct`, `min_set_ct`, and
`min_set_ct_each`, the order of operations is as follows:

1.  `min_set_ct_each`, `min_set_ct` - these filters are applied before
    filtering genes, in order to ensure all genes are present from the
    start.

2.  `min_gene_ct` - genes are filtered after pathway filtering, in order
    to remove pathways which were not deemed "significant" based upon
    the required number of genes. Only after those pathways are removed
    can the number of occurrences of each gene be judged appropriately.

## See also

Other custom plot functions:
[`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md),
[`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md),
[`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
