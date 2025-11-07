
# multienrichjam-data objects

#' Reese 2019 cohort multienrichment data
#' 
#' Data used for multi-enrichment analysis as described in
#' [Reese et al 2019](https://doi.org/10.1016/j.jaci.2018.11.043)
#' doi:10.1016/j.jaci.2018.11.043,
#' **Epigenome-wide meta-analysis of DNA methylation and childhood asthma**.
#' It is provided for convenient testing of multienrichjam.
#' 
#' Memtest is a `Mem` object with enrichment results
#' for 'Newborns' and 'OlderChildren'.
#' 
#' A minor change was made to published data:
#' * Gene symbols were updated to current Entrez gene symbols
#' as of November 2025. This process changed 'INADL' to 'PATJ'
#' in the enrichment results.
#' * Two additional genes were updated in 'geneHitIM' for genes
#' not present in the IPA pathway results.
#' * The `Reese_genes` data were also updated for consistency,
#' and to ensure that the `names(Reese_genes)` using Entrez ID
#' would be consistent with the gene symbols shown.
#' Finally, the genes in `Reese_genes` should be usable for
#' other pathway enrichment tests.
#' 
#' @docType data
#' @family Mem
#' @format `Memtest` is a `Mem` multi-enrichment object including two
#'    enrichment tests 'Newborns' and 'OlderChildren'
#' 
#' @source [Reese et al 2019](https://doi.org/10.1016/j.jaci.2018.11.043)
#' doi:10.1016/j.jaci.2018.11.043,
#' **Epigenome-wide meta-analysis of DNA methylation and childhood asthma**
#' 
#' @examples
#' # Mem data for Reese 2019 cohorts
#' data(Memtest)
#' # enrichment names
#' names(Memtest)
#' # genes represented
#' genes(Memtest)
#' # gene sets, pathways represented
#' sets(Memtest)
#' 
#' @name Mem-data
#' @aliases Memtest
"Memtest"

# describeIn Mem-data `Memtest` contains `Mem` data from Reese et al 2019.

#' @format `Reese_genes` is a `list` with names 'Newborns' and 'OlderChildren',
#'    each with a `character` vector of gene symbols, named by ENTREZID,
#'    representing genes associated with each cohort.
#'    For example, the `names(Reese_genes$Newborns)` may be useful for
#'    pathway enrichment methods that expect Entrez gene ID as input.
#' 
#' @examples
#' # Reese 2019 genes by cohort
#' data(Reese_genes)
#' # Number of genes per cohort
#' lengths(Reese_genes)
#' 
#' # first 6 gene symbols for Newborns
#' head(Reese_genes$Newborns)
#' 
#' # first 6 Entrez gene IDs for OlderChildren
#' head(names(Reese_genes$OlderChildren))
#' 
#' @rdname Mem-data
"Reese_genes"

