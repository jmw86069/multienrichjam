# Plot function for MemPlotFolio objects

Plot function for MemPlotFolio objects

## Usage

``` r
plot_mpf(
  Mpf,
  plot_which = c("all", "EnrichmentHeatmap", "eh", "GenePathHeatmap", "gphm",
    "CnetCollapsed", "cc", "CnetExemplar", "ce", "CnetCluster", "c", "EnrichmentMap",
    "em"),
  do_md_tabs = FALSE,
  md_tab_open = NULL,
  md_tab_suffix = NULL,
  md_tab_close = NULL,
  md_tab_level = 2,
  md_title = "Multi-Enrichment Folio",
  cc_type = NULL,
  ce_num = NULL,
  c_cluster = NULL,
  em_group = NULL,
  Mem = NULL,
  params = list(em = list(repulse = 3.5, width = 30)),
  do_newpage = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- Mpf:

  `MemPlotFolio` object as returned by
  [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  or
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md).

- plot_which:

  `character` vector with one or more plots, default 'all' renders all
  plots defined in 'Mpf'. Plots will be created in the order provided,
  except for 'all' which follows the default order. The following values
  are recognized:

  - 'all': All available plot types are included.

  - 'EnrichmentHeatmap' or 'eh':
    [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - 'GenePathHeatmap' or 'gp':
    [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - 'CnetCollapsed' or 'cc':
    [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - 'CnetExemplar' or 'ce':
    [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - 'CnetCluster' or 'c':
    [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - 'EnrichmentMap' or 'em':
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).
    Note this plot requires either 'Mem' argument, or 'Mem' being
    included in the `Mpf@metadata$Mem`.

- do_md_tabs:

  `logical` default FALSE, whether to print markdown-compatible tab
  headers before rendering each plot, done using
  [`cat()`](https://rdrr.io/r/base/cat.html) to STDOUT. To enable
  specific Rmarkdown or Quarto (Qmd) style, provide a `character`
  string, which will override all other checks, and will define
  'md_tab_open' and 'md_tab_close'.

  - 'Rmd' or 'rmd' to enable Rmarkdown '.tabset' style

  - 'Qmd' or 'qmd' or 'Quarto' or 'quarto' to enable Quarto ':::
    .panel-tabset' style, with footer ':::'.

- md_tab_open:

  `character` string to define the tab open style. When NULL, it will
  auto-detect an appropriate style.

  - Rmarkdown with HTML output: `''`

  - Quarto with HTML output: `'\\n\\n::: {.panel-tabset}\\n\\n'`.

  - All others: `''`.

- md_tab_suffix:

  `character` string to define the tab suffix style. When NULL, it will
  auto-detect an appropriate style.

  - Rmarkdown with HTML output: `' {.tabset}\\n\\n'`

  - All others: `''`.

- md_tab_close:

  `character` string to define the tab close style. When NULL, it will
  auto-detect an appropriate style.

  - Rmarkdown with HTML output: `''`

  - Quarto with HTML output: `'\\n\\n:::\\n\\n'`.

  - All others: `''`

  In Quarto, to keep the tabset open, define `md_tab_close=''`.

- md_tab_level:

  `integer` heading level to begin, default 2 uses heading '##'.

- md_title:

  `character` default 'Multi-Enrichment Folio' will print a title
  heading using `md_tab_level`, then create the tabset underneath at a
  lower heading.

  - Use md_title=FALSE or md_title=NULL to suppress this header.

- cc_type, ce_num, c_cluster:

  passed to corresponding functions, and supports multiple values.

  - cc_type: passed to CnetCollapsed(Mpf, type=cc_type). When NULL it
    will use 'cluster' if `ClusterLabels(Mpf)` is available, otherwise
    'set'. If more than 500 nodes, it appends '2' to the end, which
    hides gene node labels.

  - ce_num: passed to CnetExemplar(Mpf, num=ce_num). When NULL it will
    use ce_num=1.

  - c_cluster: passed to CnetCluster(Mpf, cluster=c_cluster). When NULL
    it will iterate each cluster in Mpf.

- em_group:

  `character` default NULL uses community detection in
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

  - 'community': uses community detection in
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

  - 'clusters': uses `Clusters(Mpf)`.

  - 'cluster_labels': uses `ClusterLabels(Mpf)`.

- Mem:

  `Mem` default NULL, only included for backwards support of Mpf objects
  which do not have 'Mem' included in metadata. Experimental.

- params:

  `list` named by shorthand plot type ('eh','gphm','em', 'cc', 'ce',
  'c') each containing a list of optional parameters relevant to each
  plot. Limited functionality currently, however check here for
  recognized arguments as they become available.

  - 'em' EnrichmentMap:

    - repulse: apply via `layout_with_qfr(g, repulse=repulse)`, default
      3.5.

    - width: `integer` word-wrap character width, default 30. Word-wrap
      is applied using
      [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
      however it uses `do_abbreviations=FALSE` and `removeGrep=NULL` so
      the primary effect is to apply word-wrap, and capitalization.

- do_newpage:

  `logical` default NULL sets TRUE when knitr is running, whether to
  call
  [`grid::grid.newpage()`](https://rdrr.io/r/grid/grid.newpage.html)
  after certain plots to encourage Rmd and Qmd to recognize end of a
  plot for the purpose of placing figures inside the correct tabset
  regions. Experimental. In some cases, knitr appears to embed a figure
  only when it detects that figure is "complete", and sometimes it does
  so after the next markdown tabset has been printed to STDOUT via
  [`cat()`](https://rdrr.io/r/base/cat.html). Using
  [`grid::grid.newpage()`](https://rdrr.io/r/grid/grid.newpage.html)
  appears to force knitr to recognize a new plot has started, therefore
  embeds the previous figure image. Also, knitr skips blank plots, so
  calling newpage, then having another plot call a newpage is not
  problematic.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` plot data, named by type, returned invisibly.

## Details

Main purpose is to permit independent plotting for MemPlotFolio objects,
which should allow custom settings for things like Rmd or Qmd tabs, or
for multi-page PDF output.

Note argument `plot_which` takes `character` vector, to be distinct from
`plot_which` which takes `integer` as used historically by
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md).

Todo:

- Future work will allow customization of each plot.

### Markdown tab output

When `do_md_tabs=TRUE` the default is to autodetect appropriate values
for `md_tab_open` and `md_tab_close`. These values are defined if
running in knitr context, with HTML output, using recognized input file
type 'Rmd', 'rmd', 'Qmd', 'qmd'. In all other cases, it sets tabs to
empty character ”.

You can set `do_md_tabs='Rmd'` to force Rmarkdown tabset style, or
`do_md_tabs='Qmd'` to force Quarto tabset style. In either case, the
tabs are defined using defaults, where: Rmd uses '.tabset' and no
footer; and Qmd uses '::: .panel-tabset' and ':::' footer.

Setting `do_md_tabs=TRUE` then all other values `''` blank will print
markdown-style headings for each plot type, without adding tabset tags.
This output may be preferred for a "long form" report where each plot
has its own distinct header, which may also appear in a table of
contents.

Secondary goals:

- Decouple all the sub-plot customizations so they can be applied
  without re-creating the object.

- Permit subsetting the Mpf data as relevant, perhaps by `Clusters()` as
  a drill-down technique.

Active development:

- Each plot will support more customization options over time,
  particularly to adjust Cnet network layouts, and other custom labeling
  options.

- Notably, support for
  [`ClusterLabels()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  is propagating through various plots, and
  [`ClusterData()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  may provide interactive capability, for example detailed mouse-over or
  event-driven data for each pathway or gene cluster.

## See also

Other custom plot functions:
[`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md),
[`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
[`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
