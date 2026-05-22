# Create a cnetplot igraph object, deprecated

Create a cnetplot igraph object, deprecated, use
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).

## Usage

``` r
cnetplotJam(
  x,
  showCategory = 5,
  categorySize = "geneNum",
  nodeLabel = c("Name", "Description", "ID"),
  foldChange = NULL,
  fixed = TRUE,
  doPlot = FALSE,
  categoryColor = "#E5C494",
  geneColor = "#B3B3B3",
  normalizeGeneSize = TRUE,
  labelCex = 0.45,
  colorSub = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  enrichResults object

- showCategory:

  integer number of categories to include in the resulting Cnet plot.

- categorySize:

  character value indicating how to size the pathway nodes, where
  `"geneNum"` sizes nodes by the number of genes in that pathway, and
  `"pvalue"` sizes nodes by the enrichment P-value using the format
  `-log10(pvalue)`.

- nodeLabel:

  character value indicating which colname in `as.data.frame(x)` to use
  as a node label. Depending upon the source of data, there may be
  alternative colnames that are more suitable as node labels.

- foldChange:

  numeric vector named by gene, or NULL. When supplied, the vector names
  must use the same nomenclature as the `x` input object, which can be
  inspected with `print(head(x@gene))`.

- fixed:

  optional argument passed to `netplot`.

- doPlot:

  logical indicating whether to plot the result.

- categoryColor, geneColor:

  character color, used to colorize category (pathway) nodes, or gene
  nodes, respectively.

- normalizeGeneSize:

  logical indicating whether to re-scale the gene node sizes so the mean
  gene node size is no larger than the median category node size. This
  option is intended to help gene and category node sizes to be in
  relative proportion.

- labelCex:

  numeric value to re-scale label text in the igraph object, and is
  applied directly to the igraph object.

- colorSub:

  character vector of valid R colors, whose names are compared to node
  names, for example `V(g)$name %in% names(colorSub)`.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

The purpose of this function is to mimic the steps in
`DOSE:::cnetplot()` except not plot the output, and provide some
customizations.

This function calls `cnetplot_internalJam()`, which among other things
adds a node attribute to the resulting `igraph`, `"nodeType"`, where
`nodeType="Gene"` identifies gene nodes, and `nodeType="Set"` identifies
pathway/gene set nodes.

## See also

Other jam deprecated functions:
[`enrichMapJam()`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md),
[`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md),
[`mem_multienrichplot()`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md),
[`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
