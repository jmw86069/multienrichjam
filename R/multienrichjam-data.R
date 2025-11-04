
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
#' for 'Newborns' and 'OlderChildren' as described in
#' 
#' @docType data
#' @family Mem
#' @format Memtest is a `Mem` multi-enrichment object including two
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
#' @name Memtest
"Memtest"


#' @format `Reese_genes` is a `list` with names 'Newborns' and 'OlderChildren',
#'    with a `character` vectors with genes associated with each cohort.
#' 
#' @examples
#' # Reese 2019 genes by cohort
#' data(Reese_genes)
#' # Number of genes per cohort
#' lengths(Reese_genes)
#' 
#' @rdname Memtest
"Reese_genes"

