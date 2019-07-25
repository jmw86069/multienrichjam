
<!-- README.md is generated from README.Rmd. Please edit that file -->
multienrichjam
==============

The goal of multienrichjam is to implement MultiEnrichMap, an extension to EnrichMap and Enrichment Map, originally envisioned and implemented as a Cytoscape plugin by the lab of Dr. Gary Bader

-   Bader lab [EnrichmentMap](https://www.baderlab.org/Software/EnrichmentMap), and
-   reference: [Merico,et al, PLoS One, 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0013984)),

and re-implemented into a larger analysis workflow in R by Dr. Guangchuang Yu in the R packages:

-   [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), and
-   [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html).

The MultiEnrichMap extension is a relatively straightforward set of functions to help enable comparison of enrichment results from multiple gene set enrichment tests.

Package Reference
-----------------

A full online function reference is available via the pkgdown documentation:

[Full multienrichjam command reference](https://jmw86069.github.io/multienrichjam)

How to install
--------------

Install using the R package `devtools` and this command:

    devtools::install_github("jmw86069/multienrichjam",
       dependencies=TRUE);

If you do not have the `devtools` package, you can install it from CRAN:

    install.packages("devtools");
