
#' Add directionality to pathway enrichment
#'
#' Add directionality to pathway enrichment
#'
#' This function uses pathway enrichment data provided as an
#' `enrichResult` object. Data is added as a new column (default "z-score")
#' with `numeric` values indicating the direction of change.
#'
#' It requires `gene_hits` with the direction of change for each gene,
#' which is used to determine the up- and down-regulated genes
#' in each gene set.
#'
#' The default equation uses the IPA z-score of activation, since the
#' z-score has beneficial properties in downstream visualization functions.
#' It is symmetric around zero, and implies reasonable thresholds.
#'
#' Note that when no genes involved in enrichment are present in the
#' `gene_hits` data, the resulting directional value is expected
#' to be `NA`, which will be treated as "non-directional" in downstream
#' visualization functions.
#'
#' ## IPA z-score activation calculation
#'
#' The IPA Activation `"z-score"` uses the following formula:
#'
#' ```R
#' z <- (N_genes_up - N_genes_down) / sqrt(N_genes_up + N_genes_down)
#' ```
#'
#' Citation for IPA formula: [Kramer 2014](https://doi.org/10.1093/bioinformatics/btt703)
#' https://doi.org/10.1093/bioinformatics/btt703,
#' as referenced: [IPA FAQ - Statistical Calculations](https://qiagen.my.salesforce-sites.com/KnowledgeBase/KnowledgeNavigatorPage?id=kA41i000000L5nQCAS&categoryName=BioX)
#'
#' ## Direction calculation (deprecated)
#'
#' The deprecated calculation for `"direction"` is shown below:
#'
#' ```R
#' direction_1 <- (N_genes_up - N_genes_down) / (N_genes);
#' # then log2(1.2 + x) transformed
#' direction <- log2(1.2 + abs(direction_1)) * sign(direction_1)
#' ```
#'
#' @family jam enrichment functions
#'
#' @returns `enrichResult` object, with directional values populated
#'    into the result column name defined by `dir_colname`.
#'    When no genes in a gene set are present in `gene_hits` the
#'    expected directional value is `NA`.
#'
#' @param x `enrichResult` object as returned by `clusterProfiler::enricher()`
#' @param gene_hits `numeric` vector, with values `c(-1, 1)` to
#'    indicate up or down, respectively, and whose `names(gene_hits)`
#'    should match the gene symbols (or values) in `x@result$geneID`.
#'    If no `names(gene_hits)` match values in the `geneID` column,
#'    the direction score will ultimately be populated with `NA` values.
#' @param delim `character` string with delimiter used in the
#'    `enrichResult` object, default: `"/"`
#' @param dir_colname `character` string indicating the column name to
#'    create, default `"z-score"`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
add_pathway_direction <- function
(x,
 gene_hits,
 delim="/",
 dir_colname="z-score",
 verbose=FALSE,
 ...)
{
   #
   result_genes <- strsplit(x@result$geneID, delim);
   result_genes_im <- venndir::list2im_opt(result_genes,
      do_sparse=FALSE)
   if (is.list(gene_hits)) {
      if (length(gene_hits) > 1) {
         stop("gene_hits may be a list of length one, or a single vector.")
      }
      gene_hits <- gene_hits[[1]];
   }
   if (!is.numeric(gene_hits)) {
      gene_hits <- as.numeric(gene_hits);
   }
   if (any(is.na(gene_hits))) {
      stop("gene_hits must be numeric, coercible to numeric, and must not contain NA values.")
   }
   gene_hits_im <- venndir::list2im_value(list(hits=gene_hits),
      do_sparse=FALSE)

   if (verbose > 1) {
      jamba::printDebug("add_pathway_direction(): ",
         "dim(result_genes_im):");
      print(dim(result_genes_im));
      print(head(result_genes_im[,head(colnames(result_genes_im), 10)], 3));
      jamba::printDebug("add_pathway_direction(): ",
         "dim(gene_hits_im):");
      print(dim(gene_hits_im));
      print(head(gene_hits_im, 3));
   }

   shared_genes <- intersect(rownames(result_genes_im),
      rownames(gene_hits_im));
   result_genes_im <- result_genes_im[shared_genes, , drop=FALSE];
   gene_hits_im <- gene_hits_im[shared_genes, , drop=FALSE];

   if (verbose > 1) {
      jamba::printDebug("add_pathway_direction(): ",
         "head(shared_genes, 20):");
      print(head(shared_genes, 20));
      jamba::printDebug("add_pathway_direction(): ",
         "dim(result_genes_im):");
      print(dim(result_genes_im));
      print(head(result_genes_im[,head(colnames(result_genes_im), 10)], 3));
      jamba::printDebug("add_pathway_direction(): ",
         "dim(gene_hits_im):");
      print(dim(gene_hits_im));
      print(head(gene_hits_im, 3));
   }


   results_hits_im <- (as.vector(gene_hits_im) * result_genes_im);
   if (verbose) {
      jamba::printDebug("add_pathway_direction(): ",
         "head(result_genes_im):");
      print(results_hits_im[head(seq_len(nrow(results_hits_im))),
         head(seq_len(ncol(results_hits_im))), drop=FALSE]);
   }

   # calculate "direction" - DEPRECATED
   # x@result$direction <- (colSums(results_hits_im) /
   #    colSums(abs(results_hits_im)))
   # x@result$direction <- jamba::log2signed(x@result$direction,
   #    offset=1.2);

   # calculate activation z-score defined by IPA
   x@result[[dir_colname]] <- (colSums(results_hits_im, na.rm=TRUE) /
         sqrt(colSums(abs(results_hits_im))));

   return(x)
}
