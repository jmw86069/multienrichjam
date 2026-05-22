# Multi-enrichment analysis with multienrichjam

Methods to analyze multiple gene set enrichment results, visualize
gene-pathway relationships, cluster pathways, create concept network
plots, and render igraph networks with edge bundling for bipartite
graphs. Extends EnrichMap (Merico, Bader) and clusterProfiler
(Guangchuang Yu) functions specifically to explore and compare multiple
enrichment results.

## Details

This package aims to enable analysis for one or multiple pathway or
functional enrichment datasets, with flexible options for publication of
data visualizations as figures.

The core workflow follows these steps:

1.  Import data for use.

2.  Call
    [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
    using relevant parameters.

3.  Visualize with
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    for review.

4.  Iterate then polish final figures.

The analysis produces:

- `Mem-class`: multi-enrichment data used for visualization.

- Data visualizations from
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  in `list` format.

Concept network (Cnet) plots can be customized using ShinyCat:

- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md) -
  interactive network adjustments

## See also

Useful links:

- <https://jmw86069.github.io/multienrichjam/>

- Report bugs at <http://github.com/jmw86069/multienrichjam/issues>

## Author

**Maintainer**: James M. Ward <jmw86069@gmail.com>
([ORCID](https://orcid.org/0000-0002-9510-2848)) \[copyright holder\]
