# Multienrichment folio of summary plots

Multienrichment folio of summary plots

## Usage

``` r
mem_plot_folio(
  mem,
  do_which = NULL,
  mpf = NULL,
  p_cutoff = NULL,
  p_floor = 1e-10,
  main = "",
  use_raster = FALSE,
  min_gene_ct = 1,
  min_set_ct = 1,
  min_set_ct_each = NULL,
  column_method = "euclidean",
  cluster_columns = NULL,
  row_method = "euclidean",
  cluster_rows = NULL,
  exemplar_range = c(1, 2, 3),
  pathway_column_split = NULL,
  pathway_column_title = LETTERS,
  gene_row_split = NULL,
  gene_row_title = letters,
  edge_color = NULL,
  cex.main = 2,
  cex.sub = 1.5,
  row_cex = 1,
  column_cex = 1,
  max_labels = 4,
  max_nchar_labels = 25,
  include_cluster_title = TRUE,
  repulse = 4,
  use_shadowText = FALSE,
  color_by_column = FALSE,
  style = "dotplot_inverted",
  enrich_im_weight = 0.3,
  gene_im_weight = 0.5,
  colorize_by_gene = TRUE,
  cluster_color_min_fraction = 0.4,
  byCols = c("composite_rank", "minp_rank", "gene_count_rank"),
  edge_bundling = "connections",
  apply_direction = NULL,
  rotate_heatmap = FALSE,
  row_anno_padding = grid::unit(3, "mm"),
  column_anno_padding = grid::unit(3, "mm"),
  returnType = c("MemPlotFolio", "list"),
  do_plot = TRUE,
  do_rmd_tabs = FALSE,
  rmd_tab_level = 2,
  rmd_tab_title = "Mem Plot Folio",
  verbose = FALSE,
  ...
)

prepare_folio(
  mem,
  do_which = NULL,
  mpf = NULL,
  p_cutoff = NULL,
  p_floor = 1e-10,
  main = "",
  use_raster = FALSE,
  min_gene_ct = 1,
  min_set_ct = 1,
  min_set_ct_each = NULL,
  column_method = "euclidean",
  cluster_columns = NULL,
  row_method = "euclidean",
  cluster_rows = NULL,
  exemplar_range = c(1, 2, 3),
  pathway_column_split = NULL,
  pathway_column_title = LETTERS,
  gene_row_split = NULL,
  gene_row_title = letters,
  edge_color = NULL,
  cex.main = 2,
  cex.sub = 1.5,
  row_cex = 1,
  column_cex = 1,
  max_labels = 4,
  max_nchar_labels = 25,
  include_cluster_title = TRUE,
  repulse = 4,
  use_shadowText = FALSE,
  color_by_column = FALSE,
  style = "dotplot_inverted",
  enrich_im_weight = 0.3,
  gene_im_weight = 0.5,
  colorize_by_gene = TRUE,
  cluster_color_min_fraction = 0.4,
  byCols = c("composite_rank", "minp_rank", "gene_count_rank"),
  edge_bundling = "connections",
  apply_direction = NULL,
  rotate_heatmap = FALSE,
  row_anno_padding = grid::unit(3, "mm"),
  column_anno_padding = grid::unit(3, "mm"),
  returnType = c("MemPlotFolio", "list"),
  do_rmd_tabs = FALSE,
  rmd_tab_level = 2,
  do_plot = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `Mem` or `list` object from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

- do_which:

  `integer` vector of plots to produce. When `do_which` is `NULL`, then
  all plots are produced. This argument is intended to produce only a
  subset of plots.

- mpf:

  `MemPlotFolio`, default NULL, used only to re-apply the same settings
  as another `MemPlotFolio`.

  - Note: `pathway_column_split` and `gene_row_split` are also taken
    from `thresholds(mpf)`, unless defined in function arguments.

  - Note: All other values in `thresholds(mpf)` are used, and not taken
    from corresponding function arguments.

- p_cutoff:

  `numeric` value, default NULL is taken from `mem`, indicating the
  enrichment P-value threshold.

- p_floor:

  `numeric` with the lowest enrichment P-value used in the color
  gradient on the Enrichment Heatmap. The purpose is to prevent very low
  P-values from shifting the color gradient too far from the `p_cutoff`
  causing those colors to be pale and nearly white.

- main:

  `character` string used as a title on Cnet plots.

- use_raster:

  `logical` default FALSE, deprecated, whether to use raster heatmaps,
  passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

  - Note that `use_raster=TRUE` may produce visual artifacts especially
    with argument `colorize_by_gene=TRUE` in
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md).
    Changing this argument is no longer supported

- min_gene_ct, min_set_ct:

  `integer` values passed to
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md).
  The `min_gene_ct` requires each set to contain `min_gene_ct` genes,
  and `min_set_ct` requires each gene to be present in at least
  `min_set_ct` sets.

- min_set_ct_each:

  `integer`, default NULL, minimum genes per set in at least one
  enrichment test.

  - Default NULL uses `thresholds(mem)$min_count` to use the same
    criteria.

  - The distinction from `min_set_ct` is that this threshold requires
    this number of genes in one enrichment, while `min_set_ct` applies
    the threshold to the combined multi-enrichment data.

- column_method, row_method:

  `character` arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  which indicate the distance method used to cluster columns and rows,
  respectively.

- cluster_columns, cluster_rows:

  `logical`, default NULL, whether to cluster columns (Sets) and rows
  (Genes), respectively. When NULL it uses default clustering with
  [`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html) and
  applies `column_method` or `row_method`, respectively.

- exemplar_range:

  `integer` vector (or `NULL`) used to create Cnet exemplar plots, using
  this many exemplars per cluster.

- pathway_column_split, gene_row_split:

  `integer` value passed as `column_split` and `row_split`,
  respectively, to
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
  indicating the number of pathway clusters, and gene clusters, to
  create in the gene-pathway heatmap. When either value is `NULL` then
  auto-split logic is used.

  - `integer` input will define this number of clusters

  - `list` input should be named by pathway cluster, in order the
    clusters will be split; and have `character` vectors that match
    pathways (sets) in the 'mem' input data.

  - `character` input should be named by pathway (set) and have values
    with the name of the pathway cluster. It can be `factor`.

- pathway_column_title, gene_row_title:

  `character` vectors passed to
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  as `column_title` and `row_title`, respectively. When one value is
  supplied, it is displayed and centered across all the respective
  splits. When multiple values are supplied, values are used to the
  number of splits, and recycled as needed. In that case, repeated
  values are made unique by
  [`jamba::makeNames()`](https://jmw86069.github.io/jamba/reference/makeNames.html).

- cex.main, cex.sub:

  `numeric` values passed to
  [`title()`](https://rdrr.io/r/graphics/title.html) which size the
  default title and sub-title in Cnet plots.

- row_cex, column_cex:

  `numeric` character expansion factor, used to adjust the relative size
  of row and column labels, respectively. A value of `1.1` will make row
  font size 10% larger.

- color_by_column:

  `logical` indicating whether to colorize the enrichment heatmap
  columns using `colorV` in the input `mem`. This argument is only
  relevant when `do_which` include `1`.

- enrich_im_weight, gene_im_weight:

  `numeric` value between 0 and 1, passed to
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
  used to apply relative weight to clustering columns and rows,
  respectively, when combining the gene-pathway incidence matrix with
  either column enrichment P-values, or row gene incidence matrix data.

- colorize_by_gene:

  `logical` passed to
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  indicating whether the heatmap body for the gene-pathway heatmap will
  be colorized using the enrichment colors for each gene.

- cluster_color_min_fraction:

  `numeric` value passed to
  [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  used to determine which enrichment colors to associate with each Cnet
  cluster.

- byCols:

  `character` vector describing how to sort the pathways within Cnet
  clusters. This argument is passed to
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md).

- edge_bundling:

  `character` string passed to
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  to control edge bundling. The default `edge_bundling="connections"`
  will bundle Cnet plot edges for genes that share the same pathway
  connections.

- apply_direction:

  `logical` or `NULL` indicating whether to indicate directionality in
  the
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  which is the first plot in the series. The default
  `apply_direction=NULL` will auto-detect whether there is
  directionality present in the data, and will set
  `apply_direction=TRUE` only when there are non-NA values that differ
  from zero.

- rotate_heatmap:

  `logical` passed to `mem_gene_path_heatmap()()` and only this
  function, default `FALSE`. It indicates whether to rotate the heatmap
  to have gene columns and pathway rows. If you find most people tilt
  their head to read the pathways, it might be preferable.

- row_anno_padding, column_anno_padding:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html) or `numeric` which
  will be converted to "mm" units. These values control the space
  between the heatmap body and row/column annotations, respectively,
  only relevant for `mem_gene_path_heatmap()()`. The value is only
  applied during `draw()` and cannot be defined in the `Heatmap` object
  itself, which is why it is included here and not
  `mem_gene_path_heatmap()()`.

- do_plot:

  `logical` indicating whether to render each plot. Default `TRUE` for
  `mem_plot_folio()`, and `FALSE` for `prepare_folio()`. In either case,
  plot data are created and returned, but `do_plot=TRUE` will draw each
  plot on a unique page, suitable for use with PDF output and
  `onefile=TRUE` for example.

- do_rmd_tabs:

  `logical` default FALSE, whether to print Rmarkdown tabs to STDOUT
  using [`cat()`](https://rdrr.io/r/base/cat.html), and only used when
  `do_plot=TRUE`. When active, before each plot is drawn, an Rmd tab
  header is printed in the form '## Plot Type .tabset', where the '##'
  heading is defined by 'rmd_tab_level'.

- rmd_tab_level:

  `integer` default 2, used when `do_rmd_tabs=TRUE` to define the
  heading level.

- rmd_tab_title:

  `character` default 'Mem Plot Folio' which is printed using
  [`cat()`](https://rdrr.io/r/base/cat.html) when `do_rmd_tabs=TRUE`. To
  skip the title, use `rmd_tab_title=NULL`, or set to ”. The title is
  printed as '## Title .tabset' and all plots are printed as
  sub-headings. The main reason not to use 'rmd_tab_title' is to to
  create the title in the Rmd (or Qmd) file itself, however when doing
  so, be sure to include ' .tabset' at the end of you heading, for
  example '# My Mem Plots .tabset'. The default here is to include a
  title to ensure all plot sub-headings also appear as tabs.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to downstream functions. Some useful
  examples:

  - `sets` is passed to
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
    which allows one to define a specific subset of sets to use in the
    gene-pathway heatmap.

  - `cell_size` is passed to
    [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
    with the option to define square cell size in the heatmap dotplot.
    However, the resulting heatmap will be at least `ncol * cell_height`
    width, and `nrow * cell_size[2]` height, in addition to the heights
    of the title and column labels, and widths of the color key and
    dendrogram.

## Value

`MemPlotFolio` object using
[`invisible()`](https://rdrr.io/r/base/invisible.html), containing each
plot object enabled by the argument `do_which`. The `MemPlotFolio-class`
data are accessible using common functions:

- [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- `Clusters()`

- `GeneClusters()`

- [`thresholds()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.html)

- [`metadata()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`Caption()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- [`CaptionLegendList()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

## Details

`prepare_folio()` and `mem_plot_folio()` both prepare data
visualizations from `Mem` data input. However, `prepare_folio()` does
not render figures, while `mem_plot_folio()` does render each resulting
figure.

Multiple figures can be added to a single PDF file using
`pdf(file, onefile=TRUE)` or `cairo_pdf(filename, onefile=TRUE)`.

The data are returned as `MemPlotFolio-class` which can be used to
create figures.

Pathways are clustered using the gene-pathway incidence matrix, then
used to define pathway clusters. This step can be customized by
supplying `pathway_column_split` as a `list` of `character` vectors
containing pathway (set) names.

### Pathway Clustering

Pathways are hierarchical clustered by default using
[`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html), with
`column_method="euclidean"`. The resulting hclust/dendrogram is split
using `pathway_column_split` with an `integer` number of sub-clusters.
The default (NULL) will determine a default based upon the total Sets in
the analysis, and is intended to be customized by the analyst.

The data to be clustered include the
[`memIM()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
matrix of Genes (rows) and Sets (columns), itself an incidence matrix.
It could therefore be clustered using `column_method="binary"`. However,
the default behavior is to append `-log10()` values from the `enrichIM`
matrix, in order to cluster both the incidence matrix, and the
enrichment P-values together. The matrix data are weighted relative to
one another using `enrich_im_weight=0.3`. Note that the enrichment
P-value matrix applies the `p_cutoff` such that values that do not meet
the threshold are considered '1' and become '0' with `-log10()`
transformation.

In order to cluster the gene-pathway incidence matrix with no influence
of the enrichment P-value matrix, use `enrich_im_weight=0`. In this case
it often works well to use `column_method="binary"`.

Caveat: Clustering is itself "imperfect": The results are not
definitive, and there is generally no one "true" answer. However, the
output is mainly intended to help organize information already present
in the data, and not to declare ground truth. As such, the results are
considered stable for the methods and parameters used, and the
interpretation is performed in that context.

The clustering methods and parameters are included in the legend
alongside the
[`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
and the
[`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
for clarity.

Our experience is that clusters do not need to be perfect to be useful,
informative, and valid.

#### Custom Cluster Function

It is also possible to supply a custom function to perform clustering,
by using argument `cluster_columns` which is passed via '...' to the
underlying function
[`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md).
This function should take a `matrix` of `numeric` values as input, and
return either hclust or dendrogram output, suitable for use with
[`cutree()`](https://rdrr.io/r/stats/cutree.html).

An interesting custom function is `cluster::daisy(x, method='gower')`
whose output can be converted to `hclust`. The 'gower' method accepts
mixed input types such as ordinal, numeric, signed, or other values.

When supplying a custom function via `cluster_columns`, we recommend
defining a label using argument `column_method`, for example
`custom_method='gower'`. This label will appear in the legend of the
resulting visualizations.

The default
[`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html) was
chosen because it is fast, includes both the
[`dist()`](https://rdrr.io/r/stats/dist.html) and
[`hclust()`](https://rdrr.io/r/stats/hclust.html) steps in one function,
and includes novel distance metrics such as `method='correlation'` which
have proven useful in other contexts.

#### Pathway Grouping

Pathways can be grouped manually by defining `pathway_column_split` as a
`list` of `character` vectors, where those vectors match pathway names
returned by `sets(Mem)` - which are also colnames of the
[`memIM()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
incidence matrix. This step may also use a subset of sets, for example
including only pathways of relevance to the downstream analysis.

### Gene Clustering

Genes in the Gene-Pathway Heatmap are clustered using similar logic as
described for Pathways. Note that gene clusters are not used directly in
other analysis steps (yet).

However,
[`score_gene_path_clusters()`](https://jmw86069.github.io/multienrichjam/reference/score_gene_path_clusters.md)
offers some metric to identify "hot spots" where a majority of genes in
a gene cluster are represented in the majority of pathways in a pathway
cluster. These "hot spots" are subject to further study, as potential
hubs for interpreting core components across similar biological
pathways.

Custom gene groups can be supplied with `gene_row_split` as a `list` of
`character` vectors whose values match `genes(Mem)`.

### Enrichment Heatmap

The Gene-Pathway Heatmap is used to define the pathway clusters (or
pathway groups), and therefore the order of pathways shown in the
Gene-Pathway Heatmap.

The same pathway order is used with the Enrichment Heatmap, for
consistency.

An important point is that the Enrichment Heatmap does not cluster nor
order pathways using enrichment P-values, and instead re-uses the
dendrogram (if relevant) from the Gene-Pathway Heatmap.

Enrichment P-values are not specific to the genes contained in each
pathway, and therefore clustering by enrichment P-value can be
misleading by grouping pathways together which are not related by the
underlying biological data.

### Cnet Collapsed Plot

The pathways in each pathway cluster (or pathway group) are collapsed
into one virtual pathway, then used to create a Concept network (Cnet)
plot. This Cnet plot is intended to balance the motivation to show
complete pathway enrichment data, with the motivation to reduce
redundant information.

The pathway clusters (or pathway groups) are the key components of the
Cnet Collapsed output. Therefore adjusting the clustering options
`enrich_im_weight` or `column_method` are the most common ways to
influence and optimize the outcome.

### Cnet Exemplar Plot

There are two general paradigms for Concept networks (Cnets).

1.  Pathways in clusters:
    [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

2.  Exemplar Pathways per cluster:
    [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

The Cnet Exemplar plot includes one pathway for each pathway cluster
`num=1` but may include two `num=2` or three `num=3` pathways per
cluster if relevant. The exemplar pathway is intended as a
representative of the cluster, which can help simplify the number of
genes displayed, while also focusing the meaning of those genes to the
pathway shown.

This option may be preferred when pathway clusters are not visibly
distinctive, for example when one pathway cluster does not appear to
have a common set of genes shared with other pathways in the same
cluster.

The alternative occurs when pathways in each cluster do share many of
the same genes, and in this case the
[`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
plot may be more effective.

Finally, a custom Cnet exemplar plot may be warranted when displaying a
few pathways with particular relevance to a research study or
experiment. The steps are described in the README.Rmd and README html
page for the multienrichjam package:
<https://jmw86069.github.io/multienrichjam>

### Cnet Cluster Plot

The last set of Concept network (Cnet) plots include all pathways in
each individual pathway cluster. These plots are useful for pathway
clusters which may have different functional sub-components, or when
pathways in a cluster do not appear to be very cohesive.

The Cnet Cluster plot helps display the relationship of pathways within
one pathway cluster, which may help provide the basis for interpreting
the results.

### Order of Plots

When using `mem_plot_folio()` plots are produced in a specific order.
The argument `do_which` may also be useful to focus the output only on a
particular plot, thereby skipping the preparatory steps used for other
plot types.

1.  **Enrichment Heatmap** (Plot \#1),
    [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
    Detailed arguments are described in
    [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md).

2.  **Gene-Pathway Heatmap** (Plot \#2)
    [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
    Detailed arguments are described in
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md).

3.  **Cnet Cluster Plot** (Plots \#3,#4,#5)
    [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
    Detailed arguments are described in
    [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
    and
    [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).

    - For `mem_plot_folio()` three styles are produced, with different
      node labeling strategies.

    - Plot \#3 uses pathway cluster titles.

    - Plot \#4 uses pathway names.

    - Plot \#5 uses pathway names, and hides gene labels.

4.  **Cnet Exemplar Plots** (Plots \#6,#7,#8)
    [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
    The number of exemplars uses `exemplar_range=c(1, 2, 3)`.

5.  **Cnet Cluster Plots** (Plots \#9,#10,#11,etc.)
    [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
    Detailed arguments are described in
    [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).

### Final Points

To define your own pathway cluster labels, define `pathway_column_title`
as a vector with length equal to `pathway_column_split`. These labels
become network node labels in subsequent plots, and in the resulting
`igraph` object.

Metadata are stored in `metadata(Mpf)` for MemPlotFolio output, and
includes:

- 'colorV': the color vector used

- 'hasDirection': whether the input data contained directionality,
  defined by any negative value in 'geneIMdirection',
  'enrichIMdirection', or any non-blank column for headers(Mem)
  'directionColname'.

## See also

Other multienrichjam core functions:
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

## Examples

``` r
data(Memtest)
Mpf <- prepare_folio(Memtest)
GenePathHeatmap(Mpf)

# take a subset of clusters and re-Mpf
MpfSub <- prepare_folio(Memtest,
   pathway_column_split=Clusters(Mpf)[c("A", "C")])
GenePathHeatmap(MpfSub)

# adjust set names
Memtest2 <- fixSetLabels(Memtest)
Mpf2 <- prepare_folio(Memtest2,
   row_names_rot=0,
   column_names_rot=60)
GenePathHeatmap(Mpf2)

# Now add a temporary buffer for the heatmap dimnames
with_ht_opts(list(DIMNAME_PADDING=grid::unit(c(2), "mm")), {
   GenePathHeatmap(Mpf2)
})

```
