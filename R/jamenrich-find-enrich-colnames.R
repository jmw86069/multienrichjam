
#' Find enrichment colnames
#'
#' Find enrichment colnames
#'
#' @returns `character` of recognized colnames named by the type
#'    of column, using `NA` for any column types not found.
#'
#' @family jam utility functions
#'
#' @param x `data.frame`, `enrichList`, `Mem`, or `list` of `data.frame`
#'    objects.
#' @param keyColname `character` default 'ID' indicating the primary
#'    identifier for each set. This column may be a numeric identifier.
#' @param nameColname `character` default 'Name', with the set name,
#'    typically a short name for each set.
#' @param descriptionColname `character` default 'Description' with the
#'    longer set description. It will use `nameColname` or `keyColname`
#'    as needed.
#' @param geneColname `character` default 'geneID' containing delimited
#'    genes associated with each enrichment result.
#' @param countColname `character` default 'Count' with the number of
#'    genes in the `geneColname` column. It will be calculated as needed.
#' @param geneRatioColname `character` default 'GeneRatio' with the
#'    numeric (decimal) ratio of test genes to pathway genes, or
#'    a character indication in the form '6/24' with 'tested/pathway'
#'    gene counts.
#' @param pvalueColname `character` default 'padjust' with the best available
#'    column to use for statistical significance of enrichment.
#' @param directionColname `character` default 'zscore' with a directional
#'    score, typically a z-score or some other reasonably scaled
#'    numeric value where the sign indicates directionality, with '+'
#'    meaning activated and '-' meaning suppressed.
#' @param pathGenes `character` default 'setSize' indicating the number
#'    of genes in each set as tested for enrichment. This number is not
#'    always reported, however it is not used by 'multienrichjam', but
#'    is required by some `clusterProfiler` functions.
#' @param geneHits `character` default 'Count' indicating the number of
#'    genes in the `geneColname` column. It will be calculated as needed.
#' @param verbose `logical` default FALSE, whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @returns `list` named by each column name argument, with one `character`
#'    value or `NULL` in each entry.
#'
#' @examples
#' newborn_txt <- system.file("extdata",
#'    "Newborns-IPA.txt",
#'    package="multienrichjam");
#' ipa_dfs <- importIPAenrichment(newborn_txt);
#' find_enrich_colnames(ipa_dfs[[1]])
#'
#' er <- enrichDF2enrichResult(ipa_dfs[[1]])
#' find_enrich_colnames(er)
#'
#' @export
find_enrich_colnames <- function
(x,
 keyColname=c("ID",
    "Name",
    "pathway",
    "itemsetID",
    "Description"),
 nameColname=c("Name",
    "pathway",
    "Description",
    "itemsetID",
    "ID"),
 descriptionColname=c("Description",
    "Name",
    "Pathway",
    "ID"),
 geneColname=c("geneID",
    "geneNames",
    "Genes"),
 countColname=c("gene_count",
    "count",
    "geneHits"),
 geneRatioColname=c("GeneRatio",
    "^Ratio"),
 pvalueColname=c("padjust",
    "p.adjust",
    "adjp",
    "padj",
    "qvalue",
    "qval",
    "q.value",
    "pvalue",
    "p.value",
    "pval",
    "FDR"),
 directionColname=c("activation.z.{0,1}score",
    "z.{0,1}score"),
 pathGenes=c("setSize",
    "pathGenes",
    "Count"),
 geneHits=c("Count",
    "geneHits",
    "gene_count"),
 verbose=FALSE,
 ...)
{
   #
   if (inherits(x, "Mem")) {
      x <- Mem@enrichList;
   }
   if (!inherits(x, c("data.frame", "list", "enrichResult"))) {
      stop("Input must be one of: 'Mem', 'data.frame', 'enrichResult', 'list'")
   }

   # find colnames
   keyColname <- find_colname(keyColname, x);
   nameColname <- find_colname(nameColname, x);
   descriptionColname <- find_colname(descriptionColname, x);
   geneColname <- find_colname(geneColname, x);
   countColname <- find_colname(countColname, x);
   geneRatioColname <- find_colname(geneRatioColname, x);
   pvalueColname <- find_colname(pvalueColname, x);
   directionColname <- find_colname(directionColname, x);
   pathGenes <- find_colname(pathGenes, x);
   geneHits <- find_colname(geneHits, x);

   # optional verbose output
   if (verbose) {
      jamba::printDebug("find_enrich_colnames(): ",
         "keyColname:", keyColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "nameColname:", nameColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "descriptionColname:", descriptionColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "geneColname:", geneColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "countColname:", countColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "geneRatioColname:", geneRatioColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "pvalueColname:", pvalueColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "directionColname:", directionColname);
      jamba::printDebug("find_enrich_colnames(): ",
         "pathGenes:", pathGenes);
      jamba::printDebug("find_enrich_colnames(): ",
         "geneHits:", geneHits);
   }

   # assemble a list
   colname_list <- list(
      keyColname=keyColname,
      nameColname=nameColname,
      descriptionColname=descriptionColname,
      geneColname=geneColname,
      countColname=countColname,
      geneRatioColname=geneRatioColname,
      pvalueColname=pvalueColname,
      directionColname=directionColname,
      pathGenes=pathGenes,
      geneHits=geneHits)
   return(colname_list)
}
