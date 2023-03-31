
#' Subset mem multienrichment object
#'
#' Subset mem multienrichment object
#'
#' This function is intended to subset the incidence matrix data contained
#' in a `mem` object by heuristics. It does not update other data in
#' the `mem` object such as `enrichList` and `multiEnrichDF`, nor
#' any `igraph` objects. It is intended mainly to subset by sets (pathways),
#' or genes, then also subset other corresponding incidence matrix
#' data consistently.
#'
#' @return `list` object of `mem` data
#'
#' @family jam utility functions
#'
#' @param mem `list` object as returned by `multiEnrichMap()`
#' @param includeSets `character` vector with specific sets to retain,
#'    all other sets will be dropped.
#' @param includeGenes `character` vector with specific genes to retain,
#'    all other genes will be dropped.
#' @param min_gene_ct `numeric` filter applied to genes representing the
#'    minimum number of occurrences across sets in the `mem$memIM`
#'    incidence matrix. The default value `min_gene_ct=1` effectively
#'    requires a gene to be present in at least one set, which is useful
#'    after filtering by `includeSets`.
#' @param min_set_ct `numeric` filter applied to sets representing the
#'    minimum number of occurrences of genes in the `mem$memIM`
#'    incidence matrix. The default value `min_set_ct=1` effectively
#'    requires a set to contain at least one gene.
#' @param p_cutoff `numeric` optional enrichment P-value filter to apply
#'    to `mem$enrichIM` enrichment P-values.
#'    It is intended to apply optionally higher stringency by using a
#'    lower `p_cutoff` than used by `multiEnrichMap()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
subset_mem <- function
(mem,
 includeSets=NULL,
 includeGenes=NULL,
 min_gene_ct=1,
 min_set_ct=1,
 p_cutoff=NULL,
 verbose=FALSE,
 ...)
{
   # filter for specific sets
   if (length(includeSets) > 0) {
      use_sets <- intersect(colnames(mem$memIM), includeSets)
      if (length(use_sets) == 0) {
         use_sets <- colnames(mem$memIM)[
            tolower(colnames(mem$memIM)) %in% tolower(includeSets)]
      }
      # if (length(use_sets) == 0) {
      #    stop("includeSets did not match any entries in memIM")
      # }
      # update all relevant set entries
      if (length(use_sets) < ncol(mem$memIM)) {
         if (verbose) {
            jamba::printDebug("subset_mem(): ",
               "subset by includeSets to ",
               jamba::formatInt(length(use_sets)),
               " sets.")
         }
         mem <- filter_mem_sets(mem,
            includeSets=use_sets)
      }
   }

   # filter sets by enrichment P-value
   if (length(p_cutoff) == 1 && p_cutoff < 1) {
      keep_set_rows <- rowSums((mem$enrichIM <= p_cutoff) * 1) > 0;
      if (any(!keep_set_rows)) {
         keep_sets <- rownames(mem$enrichIM)[keep_set_rows];
         # update all relevant set entries
         if (verbose) {
            jamba::printDebug("subset_mem(): ",
               "subset by p_cutoff<=",
               p_cutoff,
               " to ",
               jamba::formatInt(length(keep_sets)),
               " sets.")
         }
         mem <- filter_mem_sets(mem,
            includeSets=keep_sets)
      }
   }

   # filter for specific genes
   if (length(includeGenes) > 0) {
      use_genes <- intersect(rownames(mem$memIM), includeGenes)
      if (length(use_genes) == 0) {
         use_genes <- rownames(mem$memIM)[
            tolower(rownames(mem$memIM)) %in% tolower(includeGenes)]
      }
      # if (length(use_genes) == 0) {
      #    stop("includeGenes did not match any entries in memIM")
      # }
      # update all relevant set entries
      if (length(use_genes) < nrow(mem$memIM)) {
         if (verbose) {
            jamba::printDebug("subset_mem(): ",
               "subset by includeGenes to ",
               jamba::formatInt(length(use_genes)),
               " genes.")
         }
         mem <- filter_mem_genes(mem,
            includeGenes=use_genes)
      }
   }

   # filter by min_set_ct
   if (length(min_set_ct) == 1 && min_set_ct > 0) {
      set_ct <- colSums((mem$memIM > 0) * 1);
      if (any(set_ct < min_set_ct)) {
         keep_sets <- colnames(mem$memIM)[set_ct >= min_set_ct]
         # update all relevant set entries
         if (verbose) {
            jamba::printDebug("subset_mem(): ",
               "subset by min_set_ct>=",
               min_set_ct,
               " to ",
               jamba::formatInt(length(keep_sets)),
               " sets.")
         }
         mem <- filter_mem_sets(mem,
            includeSets=keep_sets)
      }
   }

   # filter by min_gene_ct
   if (length(min_gene_ct) == 1 && min_gene_ct > 0) {
      gene_ct <- rowSums((mem$memIM > 0) * 1);
      if (any(gene_ct < min_gene_ct)) {
         keep_genes <- rownames(mem$memIM)[gene_ct >= min_gene_ct]
         # update all relevant set entries
         if (verbose) {
            jamba::printDebug("subset_mem(): ",
               "subset by min_gene_ct>=",
               min_gene_ct,
               " to ",
               jamba::formatInt(length(keep_genes)),
               " genes.")
            printDebug("keep_genes: ", keep_genes);
         }
         mem <- filter_mem_genes(mem,
            includeGenes=keep_genes)
      }
   }

   return(mem)
}

#' Filter mem multienrichment object by Set names
#'
#' Filter mem multienrichment object by Set names
#'
#' This is intended to be an internal function. It simply takes a
#' `character` vector of set names, and subsets incidence matrix data
#' in the `mem` object. It performs no other filtering.
#' This function is called by `subset_mem()`, the recommended method
#' to subset a `mem` object.
#'
#' @return `list` object of `mem` data content
#'
#' @family jam utility functions
#'
#' @param mem `list` object returned by `multiEnrichMap()`
#' @param includeSets `character` vector of sets, matched with
#'    `colnames(mem$memIM)`
#' @param ... additional arguments are ignored
#'
filter_mem_sets <- function
(mem,
 includeSets=NULL,
 ...)
{
   #
   includeSets <- intersect(includeSets,
      colnames(mem$memIM));
   # update all relevant column set entries
   col_entries <- vigrep("^memIM", names(mem))
   for (i in col_entries) {
      col_match <- match(includeSets, colnames(mem[[i]]))
      mem[[i]] <- mem[[i]][, col_match, drop=FALSE]
   }
   # update all relevant row set entries
   row_entries <- vigrep("^enrichIM", names(mem))
   for (i in row_entries) {
      row_match <- match(includeSets, rownames(mem[[i]]))
      mem[[i]] <- mem[[i]][row_match, , drop=FALSE]
   }
   return(mem)
}

#' Filter mem multienrichment object by Gene names
#'
#' Filter mem multienrichment object by Gene names
#'
#' This is intended to be an internal function. It simply takes a
#' `character` vector of gene names, and subsets incidence matrix data
#' in the `mem` object. It performs no other filtering.
#' This function is called by `subset_mem()`, the recommended method
#' to subset a `mem` object.
#'
#' @return `list` object of `mem` data content
#'
#' @family jam utility functions
#'
#' @param mem `list` object returned by `multiEnrichMap()`
#' @param includeGenes `character` vector of genes, matched with
#'    `rownames(mem$memIM)`
#' @param ... additional arguments are ignored
#'
filter_mem_genes <- function
(mem,
 includeGenes=NULL,
 ...)
{
   #
   includeGenes <- intersect(includeGenes,
      rownames(mem$memIM));
   # no column entries are relevant for genes
   #
   # update all relevant column set entries
   row_entries <- vigrep("^memIM|^geneIM", names(mem))
   for (i in row_entries) {
      row_match <- match(includeGenes, rownames(mem[[i]]))
      mem[[i]] <- mem[[i]][row_match, , drop=FALSE]
   }
   return(mem)
}
