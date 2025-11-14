
#' Score Gene-Path Clusters in MemPlotFolio
#' 
#' Score Gene-Path Clusters in MemPlotFolio, still in development
#' 
#' This function provides a quick summary of the "cohesiveness"
#' of gene-pathway clustering, using a basic metric to summarize
#' the fraction of genes represented in each pathway cluster.
#' 
#' The summary may be useful to identify "hot spots" in a gene-pathway
#' heatmap, which may itself be useful to identify functional
#' hubs which may be broadly shared across numerous pathways.
#' 
#' The simple metric appears effective at identifying gene clusters
#' which are populated by the majority of genes, and therefore could
#' be the basis for understanding the nature of common genes associated
#' with multiple pathways.
#' It also appears interesting to observe pathway clusters,
#' some of which have one or more gene clusters with high gene-path cluster
#' scores, and some which have moderate gene-path cluster scores across
#' numerous gene clusters. The latter case may be an artifact of
#' small pathway clusters, in which case a better metric should
#' be evaluated.
#' 
#' @family MemPlotFolio
#' 
#' @returns `matrix` of numeric values ranging from 0 to 1, reflecting
#'    the fraction of genes with at least `set_fraction` representation
#'    in each pathway cluster. The columns are pathway clusters, and
#'    rows are gene clusters.
#' 
#' @param Mpf `MemPlotFolio` with `Clusters()` and `GeneClusters()` used
#'    to define pathway (set) clusters, and gene clusters, respectively.
#'    The data are taken from the `GenePathHeatmap()` heatmap data,
#'    which are equivalent to using the `memIM(Mem)` gene-pathway
#'    incidence matrix data.
#' @param set_fraction `numeric` value, default 0.5, indicating the fraction
#'    of pathways in a pathway cluster at or above which a gene is required
#'    to be present for that gene to contribute to the gene-pathway cluster
#'    score.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' Mpf <- prepare_folio(fixSetLabels(Memtest), do_plot=FALSE, do_which=2, column_cex=0.5)
#' GenePathHeatmap(Mpf)
#' gpscores <- score_gene_path_clusters(Mpf)
#' 
#' 
#' # gpscores could even be used to create a network
#' gpnet <- mem2cnet((gpscores > 0.2) * 1)
#' jam_igraph(gpnet)
#' 
#' @export
score_gene_path_clusters <- function
(Mpf,
 set_fraction=0.5,
 ...)
{
   #
   gp_hm <- GenePathHeatmap(Mpf, do_plot=FALSE)
   gpm <- gp_hm@matrix;
   gro <- GeneClusters(Mpf)
   pco <- Clusters(Mpf)
   
   gpscores <- jamba::rbindList(lapply(gro, function(i){
      pset <- sapply(pco, function(j){
         # jamba::printDebug("head(i): ", head(i), ", head(j): ", head(j));
         gpm1 <- gpm[i, j, drop=FALSE]
         gpm1[] <- 1 * (!gpm1 %in% c(0, NA))
         ggood <- rowSums(gpm1) >= (ncol(gpm1) * set_fraction);
         sum(ggood) / length(ggood)
      })
   }))
   gpscores
}
