# methods-Mem.R
#
# Todo:
# as.list()
# dimnames()<-

# Future expansion:
# genesPerSet()
# setsWithGene()
# geneInCategory(): mimic generic function in `DOSE::geneInCategory()`
#
# cbind(Mem1, Mem2): combine Mem by adding new enrichment data
#    for example adding enrichments into the same Mem
# rbind(Mem1, Mem2): combine Mem by adding new pathways
#    for example, same enrichments but new pathway source

#' Extract some elements of potentially long vector for display
#'
#' Extract some elements of potentially long vector for display,
#' using ellipses when there are more than 'maxToShow' entries.
#' 
#' @keywords internal
#' @noRd
some_vector <- function
(x,
 maxToShow=5,
 ellipsis="...",
 ellipsisPos=c("middle", "start", "end"),
 quote=FALSE,
 sep=", ",
 ...)
{
   #
   ellipsisPos <- match.arg(ellipsisPos)
   if (length(x) == 0 || !is.atomic(x)) {
      return(x)
   }
   if (is.character(x) && TRUE %in% quote) {
      x <- sQuote(x)
   }
   len <- length(x)
   if (length(maxToShow) == 0 || head(maxToShow, 1) < 3) {
      maxToShow <- 3
   }
   if (len > maxToShow) {
      maxToShow <- maxToShow - 1;
      if ("end" %in% ellipsisPos) {
         y <- c(head(x, maxToShow - 1), ellipsis)
      } else if ("start" %in% ellipsisPos) {
         y <- c(ellipsis, head(x, maxToShow - 1))
      } else {
         bot <- ceiling(maxToShow / 2)
         top <- len - (maxToShow - bot - 1)
         y <- c(as.character(x[seq_len(bot)]),
            ellipsis,
            as.character(x[seq(from=top, to=len)]))
      }
   } else {
      y <- x
   }
   jamba::cPaste(y, sep=sep)
}

#' Subset a Mem object
#'
#' @importFrom methods new setClass setGeneric setMethod
#' @importFrom methods setReplaceMethod setValidity prototype
#' 
#' @param x `Mem` object
#' @param i `character` or `integer` with gene names
#' @param j `character` or `integer` with pathway names
#' @param ... k for '[' `character` or `integer` with enrichment names
#' @param drop `logical` always set to FALSE for Mem objects.
#' @docType methods
#' @describeIn Mem-class Subset a Mem object, `Mem[i, j, k]`,
#'    where `i` are genes, `j` are pathways, and `k` are enrichments.
#' @export
setMethod("[",
   signature=c(x="Mem"),
   definition=function(x, i, j, k, drop=FALSE) {
      if (missing(i)) i <- NULL;
      if (missing(j)) j <- NULL;
      if (missing(k)) k <- NULL;
      if (missing(drop)) drop <- FALSE;
      
      ivals <- genes(x);
      jvals <- sets(x);
      kvals <- enrichments(x);
      
      # row subset by gene name
      if (length(i) > 0) {
         if (is.numeric(i)) {
            if (any(abs(i) > length(ivals))) {
               stop("i is out of bounds, does not match set indices.");
            }
            i <- ivals[i];
         } else if (inherits(i, c("character", "factor"))) {
            i <- as.character(i);
            if (!all(i %in% ivals)) {
               stop("i is out of bounds, does not match character gene names.")
            }
         }
         # update if provided values of different length, or
         # with any change in values
         i_update <- length(i) > 0 && (
            !all(ivals %in% i) ||
               (length(i) == length(ivals) && !all(ivals == i)) );
         if (i_update) {
            islot_cols <- NULL
            islot_rows <- c("memIM",
               "geneIM",
               "geneIMdirection",
               "geneIMcolors")
            # subset by colnames
            for (islot in islot_cols) {
               imatch <- match(i, colnames(slot(x, islot)))
               slot(x, islot) <- slot(x, islot)[, imatch, drop=FALSE];
            }
            # subset by rownames
            for (islot in islot_rows) {
               imatch <- match(i, rownames(slot(x, islot)))
               slot(x, islot) <- slot(x, islot)[imatch, , drop=FALSE];
            }
            # optionally subset enrichList,multiEnrichDF,multiEnrichResult
            # for matching set names
         }
      }
      # row subset by set name
      if (length(j) > 0) {
         if (is.numeric(j)) {
            if (any(abs(j) > length(jvals))) {
               stop("j is out of bounds, does not match set indices.");
            }
            j <- jvals[j];
         } else if (inherits(j, c("character", "factor"))) {
            j <- as.character(j);
            if (!all(j %in% jvals)) {
               stop("j is out of bounds, does not match character set names.")
            }
         }
         j_update <- length(j) > 0 && (
            !all(jvals %in% j) ||
               (length(j) == length(jvals) && !all(jvals == j)) );
         if (j_update) {
            jslot_cols <- c("memIM")
            jslot_rows <- c("enrichIM",
               "enrichIMcolors",
               "enrichIMdirection",
               "enrichIMgeneCount")
            # subset by colnames
            for (jslot in jslot_cols) {
               jmatch <- match(j, colnames(slot(x, jslot)))
               slot(x, jslot) <- slot(x, jslot)[, jmatch, drop=FALSE];
            }
            # subset by rownames
            for (jslot in jslot_rows) {
               jmatch <- match(j, rownames(slot(x, jslot)))
               slot(x, jslot) <- slot(x, jslot)[jmatch, , drop=FALSE];
            }
            # optionally subset enrichList,multiEnrichDF,multiEnrichResult
            # for matching set names
         }
      }
      # column subset by enrichment name
      if (length(k) > 0) {
         if (is.numeric(k)) {
            if (any(abs(k) > length(kvals))) {
               stop("k is out of bounds, does not match set indices.");
            }
            k <- kvals[k];
         } else if (inherits(k, c("character", "factor"))) {
            k <- as.character(k);
            if (!all(k %in% kvals)) {
               stop("k is out of bounds, does not match character enrichment names.")
            }
         }
         k_update <- length(k) > 0 && (
            !all(kvals %in% k) ||
               (length(k) == length(kvals) && !all(kvals == k)) );
         if (k_update) {
            kslot_names <- c("enrichList",
               "enrichLabels",
               "colorV",
               "geneHitList")
            kslot_cols <- c(
               "enrichIM",
               "enrichIMcolors",
               "enrichIMdirection",
               "enrichIMgeneCount",
               "geneHitIM",
               "geneIM",
               "geneIMcolors",
               "geneIMdirection")
            kslot_rows <- NULL;
            # subset by names
            for (kslot in kslot_names) {
               kmatch <- match(k, names(slot(x, kslot)))
               slot(x, kslot) <- slot(x, kslot)[kmatch];
            }
            # subset by colnames
            for (kslot in kslot_cols) {
               kmatch <- match(k, colnames(slot(x, kslot)))
               slot(x, kslot) <- slot(x, kslot)[, kmatch, drop=FALSE];
            }
            # subset by rownames
            for (kslot in kslot_rows) {
               kmatch <- match(k, rownames(slot(x, kslot)))
               slot(x, kslot) <- slot(x, kslot)[kmatch, , drop=FALSE];
            }
            # Update memIM with new contents in enrichIM, geneIM
            geneColname <- x@headers[["geneColname"]];
            nameColname <- x@headers[["nameColname"]];
            geneDelim <- "[,/ ]+";
            # This step seems problematic, to use multiEnrichDF
            # which is not updated effectively by sets()<- and genes()<-.
            # Long-term it might be useful, for now use simpler approach
            if (FALSE) {
               new_memIM <- list2im(
                  keepCounts=TRUE,
                  strsplit(
                     jamba::nameVector(
                        as.character(x@multiEnrichDF[[geneColname]]),
                        x@multiEnrichDF[[nameColname]]),
                     geneDelim));
            } else {
               # strategy: multiply rowSums of the current geneIM
               # which will properly update values to reflect the number
               # of enrichments in which the gene was found.
               geneIMrowsum <- rowSums(
                  jamba::rmNA(naValue=0, abs(geneIM(x))),
                  na.rm=TRUE);
               new_memIM <- sign(memIM(x)) * geneIMrowsum;
            }
            # empty the existing memIM since we cannot yet subset the columns
            # instead populate with columns that remain,
            # then consider subsetting columns afterward
            x@memIM <- x@memIM * 0;
            kmatchrow <- match(rownames(new_memIM), rownames(x@memIM));
            kmatchcol <- match(colnames(new_memIM), colnames(x@memIM));
            x@memIM[kmatchrow, kmatchcol] <- new_memIM;
            #
            # Todo:
            # (Optional?) Subset memIM columns to remove empty pathways
         }
      }
      x;
   }
)

#' @describeIn Mem-class Show summary of a Mem object, dimensions defined by
#'    genes, sets, enrichments, showing up to 5 entries of each.
#'    Also includes a list of analysis parameters: topN, and thresholds.
setMethod("show", "Mem",
   function(object) {
      nenrich <- jamba::formatInt(length(object@enrichList))
      nsets <- jamba::formatInt(nrow(object@enrichIM))
      ngenes <- jamba::formatInt(nrow(object@geneIM))
      # ngenehits <- nrow(object@geneHitIM);
      if ("p_cutoff" %in% names(object@thresholds)) {
         p_cutoff <- object@thresholds$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(object@thresholds)) {
         p_cutoff <- object@thresholds$cutoffRowMinP;
      } else {
         p_cutoff <- 1;
      }
      if ("topEnrichN" %in% names(object@thresholds)) {
         top_enrich_n <- object@thresholds$topEnrichN;
      } else {
         top_enrich_n <- object@thresholds$top_enrich_n;
      }
      min_count <- object@thresholds$min_count;
      summary_v <- c(
         paste0("class: Mem"),
         paste0("dim: ", nsets, " ", nenrich)
      )
      
      summary_v <- c(
         paste0("class: Mem"),
         paste0("dim: ",
            nenrich, " enrichments, ",
            nsets, " sets, ",
            ngenes, " genes"),
         paste0("- enrichments (", nenrich, "): ",
            some_vector(names(object@enrichList), sep=", ")),
         paste0("- sets (", nsets, "): ",
            some_vector(colnames(object@memIM), sep=", ")),
         paste0("- genes (", ngenes, "): ",
            some_vector(rownames(object@memIM), sep=", ")),
         paste0("Analysis parameters:"),
         paste0("- top N per enrichment: ", top_enrich_n),
         paste0("- significance threshold: ", p_cutoff,
            " (colname: ",
            jamba::cPaste(object@headers$pvalueColname, sep=", "),
            ")"),
         paste0("- min gene count: ", min_count)
      )
      if ("directionColname" %in% names(object@headers) &&
      		length(object@headers[["directionColname"]]) > 0) {
         summary_v <- c(summary_v,
            paste0("- direction colname: ",
               object@headers$directionColname))
      }
      summary_text <- jamba::cPaste(summary_v, sep="\n")
      cat(summary_text, sep="");
      # jamba::sdim(object)[, c("enrichIM", "geneIM", "memIM"), drop=FALSE]
   }
)



#' @describeIn Mem-class Coerce S4 `Mem` to legacy `list` mem format
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#'
#' @returns `Mem_to_list()` returns a `list` suitable for other mem functions.
#'
#' @examples
#' # examples for Mem_to_list and as(Mem, "list")
#' data(Memtest)
#' mem1 <- Mem_to_list(Memtest)
#' jamba::sdim(mem1)
#'
#' @export
Mem_to_list <- function
(x,
 ...)
{
   #
   if (!inherits(x, "Mem")) {
      stop("Input must be 'Mem'")
   }
   
   use_slots=c(
      enrichList="list",
      enrichLabels="character",
      colorV="character",
      geneHitList="list",
      geneHitIM="matrix",
      memIM="matrix",
      geneIM="matrix",
      enrichIM="matrix",
      multiEnrichDF="data.frame",
      multiEnrichResult="enrichResult",
      thresholds="list",
      headers="list",
      enrichIMcolors="matrix",
      enrichIMdirection="matrix",
      enrichIMgeneCount="matrix",
      geneIMcolors="matrix",
      geneIMdirection="matrix"
   )
   use_slotnames <- names(use_slots);
   names(use_slotnames) <- use_slotnames;
   new_list <- lapply(use_slotnames, function(slotname){
      slot(x, slotname)
   })
   # rename 'headers' to 'colnames'
   hmatch <- match("headers", names(new_list));
   if (length(hmatch) == 1) {
      names(new_list)[hmatch] <- "colnames";
   }
   # et voila
   new_list
}

#' @importFrom S4Vectors as.list
#' @param x `Mem` object
#' @docType methods
#' @describeIn Mem-class Coerce S4 `Mem` to legacy `list` mem format
#' @returns `as.list(Mem)` returns a `list` suitable for other mem functions.
#' @export
setMethod("as.list", "Mem", function(x) Mem_to_list(x))

#' The enrichment names from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Names for the enrichment tests in a Mem object
#' @returns `names()` returns the `character` vector of enrichment names
#' @export
setMethod("names", "Mem", function(x) names(x@enrichList))

#' The enrichment names from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Names for the enrichment tests in a Mem object
#' @returns `enrichments()` returns the `character` vector of enrichment names
#' @export
setMethod("enrichments", "Mem", function(x) names(x@enrichList))

#' Assign enrichment names from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign names for the enrichment tests in a Mem object
#' @export
setMethod("enrichments<-", "Mem", function(x, value) {
   kslot_names <- c("enrichList",
      "enrichLabels",
      "colorV",
      "geneHitList")
   kslot_cols <- c(
      "enrichIM",
      "enrichIMcolors",
      "enrichIMdirection",
      "enrichIMgeneCount",
      "geneHitIM",
      "geneIM",
      "geneIMcolors",
      "geneIMdirection")
   for (use_k in kslot_names) {
      xobj <- slot(x, use_k);
      if (length(value) != length(xobj)) {
         stop(paste0("length is not equal to length(x@", use_k, ")"));
      }
      names(xobj) <- value;
      slot(x, use_k, check=FALSE) <- xobj;
   }
   for (use_k in kslot_cols) {
      xobj <- slot(x, use_k);
      if (length(value) != ncol(xobj)) {
         stop(paste0("length is not equal to ncol(x@", use_k, ")"));
      }
      colnames(xobj) <- value;
      slot(x, use_k, check=FALSE) <- xobj;
   }
   x
})

#' List pathway gene set names represented in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List pathway gene set names in a Mem object
#' @returns `sets()` returns the `character` vector of pathway gene sets
#' @export
setMethod("sets", "Mem",
   function(x, ...) rownames(x@enrichIM)
)

#' Assign pathway gene set names represented in a Mem object
#'
#' @param x `Mem` object
#' @param value `character` vector of names to assign to pathway gene sets
#' @docType methods
#' @describeIn Mem-class Assign pathway gene set names to a Mem object
#' @export
setMethod("sets<-", "Mem",
   function(x, value) {
      jslot_cols <- c("memIM")
      jslot_rows <- c("enrichIM",
         "enrichIMcolors",
         "enrichIMdirection",
         "enrichIMgeneCount")
      for (use_j in jslot_cols) {
         if (length(value) != ncol(slot(x, use_j))) {
            stop(paste0("length is not equal to ncol(x@", use_j, ")"));
         }
         xmatrix <- slot(x, use_j);
         colnames(xmatrix) <- value;
         slot(x, use_j, check=FALSE) <- xmatrix;
      }
      for (use_j in jslot_rows) {
         if (length(value) != nrow(slot(x, use_j))) {
            stop(paste0("length is not equal to nrow(x@", use_j, ")"));
         }
         xmatrix <- slot(x, use_j);
         rownames(xmatrix) <- value;
         slot(x, use_j, check=FALSE) <- xmatrix;
      }
      x
   }
)

#' List genes represented in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List genes represented
#' @export
setMethod("genes", "Mem",
   function(x, ...) rownames(x@geneIM)
)

#' Assign gene names in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign gene names in a Mem object
#' @export
setMethod("genes<-", "Mem",
   function(x, value) {
      islot_rows <- c("memIM",
         "geneIM",
         "geneIMdirection",
         "geneIMcolors")
      
      # attempt to update geneHitIM for rows that match genes(x)
      genes_input <- genes(x);
      match_genehitim <- match(genes_input, rownames(x@geneHitIM))
      match_genehitim_nonna <- !is.na(match_genehitim);
      if (any(match_genehitim_nonna)) {
         # as long as some values match, do the updates
         k <- match_genehitim[match_genehitim_nonna]
         rownames(x@geneHitIM)[k] <- value[match_genehitim_nonna];
      }
      if (!all(match_genehitim_nonna)) {
         # consider warning that geneHitIM is not aligned with geneIM
         warn_msg <- paste0("genes<- is not properly aligned with ",
            "geneHitIM, ", 
            "because genes(x) are not represented in rownames(geneHitIM).")
         warning(warn_msg);
      }
      
      for (use_i in islot_rows) {
         xobj <- slot(x, use_i);
         if (length(value) != nrow(xobj)) {
            stop(paste0("length is not equal to nrow(x@", use_i, ")"));
         }
         rownames(xobj) <- value;
         slot(x, use_i, check=FALSE) <- xobj;
      }
      x
   }
)

#' The matrix of genes tested versus enrichment tests from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The gene-enrichment matrix of genes represented in
#'    enrichment tests.
#' @export
setMethod("geneIM", "Mem",
   function(x, ...) x@geneIM
)

#' The matrix of genes tested versus enrichment tests from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The matrix of genes tested versus enrichment tests, with directionality.
#' @export
setMethod("geneIMdirection", "Mem",
   function(x, ...) x@geneIMdirection
)

#' Assign directional matrix of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign directional matrix to the gene-enrichment matrix
#' @export
setMethod("geneIMdirection<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(genes(x) %in% rownames(value))) {
         stop_msg <- paste0("genes(x) must all be represented in ",
            "rownames(value).")
         stop(stop_msg);
      }
      # enrichments(x) must be present in colnames(value)
      if (!all(enrichments(x) %in% colnames(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "colnames(value).")
         stop(stop_msg);
      }
      imatch <- match(genes(x), rownames(value));
      jmatch <- match(enrichments(x), colnames(value));
      x@geneIMdirection <- value[imatch, jmatch, drop=FALSE];
      x
   }
)

#' The matrix of genes tested versus enrichment tests from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The matrix of colors indicating genes in each enrichment test
#' @export
setMethod("geneIMcolors", "Mem",
   function(x, ...) x@geneIMcolors
)

#' Assign a color matrix of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign a color matrix to the gene-enrichment matrix
#' @export
setMethod("geneIMcolors<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(genes(x) %in% rownames(value))) {
         stop_msg <- paste0("genes(x) must all be represented in ",
            "rownames(value).")
         stop(stop_msg);
      }
      # enrichments(x) must be present in colnames(value)
      if (!all(enrichments(x) %in% colnames(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "colnames(value).")
         stop(stop_msg);
      }
      imatch <- match(genes(x), rownames(value));
      jmatch <- match(enrichments(x), colnames(value));
      x@geneIMcolors <- value[imatch, jmatch, drop=FALSE];
      x
   }
)

#' The list of enrichResult data in an Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The list of enrichResult data in an Mem object
#' @export
setMethod("enrichList", "Mem",
   function(x, ...) x@enrichList
)

#' The pathway/P-value matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The pathway/P-value matrix
#' @export
setMethod("enrichIM", "Mem",
   function(x, ...) x@enrichIM
)

#' The pathway-enrichment color matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The pathway-enrichment directional color matrix
#' @export
setMethod("enrichIMcolors", "Mem",
   function(x, ...) x@enrichIMcolors
)

#' Assign the pathway-enrichment directional color matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign the pathway-enrichment directional color matrix
#' @export
setMethod("enrichIMcolors<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(sets(x) %in% rownames(value))) {
         stop_msg <- paste0("sets(x) must all be represented in ",
            "rownames(value).")
         stop(stop_msg);
      }
      # enrichments(x) must be present in colnames(value)
      if (!all(enrichments(x) %in% colnames(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "colnames(value).")
         stop(stop_msg);
      }
      imatch <- match(sets(x), rownames(value));
      jmatch <- match(enrichments(x), colnames(value));
      x@enrichIMcolors <- value[imatch, jmatch, drop=FALSE];
      x
   }
)

#' The pathway-enrichment directional score matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The pathway-enrichment directional score matrix
#' @export
setMethod("enrichIMdirection", "Mem",
   function(x, ...) x@enrichIMdirection
)

#' Assign the pathway-enrichment directional score matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign the pathway-enrichment directional score matrix
#' @export
setMethod("enrichIMdirection<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(sets(x) %in% rownames(value))) {
         stop_msg <- paste0("sets(x) must all be represented in ",
            "rownames(value).")
         stop(stop_msg);
      }
      # enrichments(x) must be present in colnames(value)
      if (!all(enrichments(x) %in% colnames(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "colnames(value).")
         stop(stop_msg);
      }
      imatch <- match(sets(x), rownames(value));
      jmatch <- match(enrichments(x), colnames(value));
      x@enrichIMdirection <- value[imatch, jmatch, drop=FALSE];
      x
   }
)

#' The pathway-enrichment gene count matrix with gene counts from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Pathway-enrichment gene count matrix, genes involved
#'    in enrichment of each pathway, for each enrichment test.
#' @export
setMethod("enrichIMgeneCount", "Mem",
   function(x, ...) x@enrichIMgeneCount
)

#' The gene-pathway incidence matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The gene-pathway incidence matrix
#' @export
setMethod("memIM", "Mem",
   function(x, ...) x@memIM
)

#' The matrix of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The matrix of genes tested for enrichment, including
#'    genes not associated with enrichment results.
#' @export
setMethod("geneHitIM", "Mem",
   function(x, ...) x@geneHitIM
)


#' Assign matrix of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign the matrix of genes tested in a Mem object
#' @export
setMethod("geneHitIM<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(genes(x) %in% rownames(value))) {
         stop_msg <- paste0("genes(x) must all be represented in ",
            "rownames(value).")
         stop(stop_msg);
      }
      # enrichments(x) must be present in colnames(value)
      if (!all(enrichments(x) %in% colnames(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "colnames(value).")
         stop(stop_msg);
      }
      x@geneHitIM <- value;
      x
   }
)

#' The list of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The list of genes tested for enrichment, including
#'    genes not associated with enrichment results.
#' @export
setMethod("geneHitList", "Mem",
   function(x, ...) x@geneHitList
)

#' Assign the list of genes tested for enrichment in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign the list of genes tested for enrichment,
#'    including genes not associated with enrichment results.
#' @export
setMethod("geneHitList<-", "Mem",
   function(x, value) {
      # enforce genes(x) must be present in rownames(value)
      if (!all(enrichments(x) %in% names(value))) {
         stop_msg <- paste0("enrichments(x) must all be represented in ",
            "names(value).")
         stop(stop_msg);
      }
      x@geneHitList <- value;
      x
   }
)

#' The column headers mapped in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Enrichment data column headers associated to
#'    enrichResult data in a Mem object
#' @export
setMethod("headers", "Mem", function(x) x@headers)

#' List colors assigned to enrichment tests in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List colors assigned to  represented
#' @export
setMethod("colorV", "Mem",
   function(x, ...) x@colorV
)

#' Assign colors to enrichment tests in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Assign colors to enrichments in a Mem object. If
#'    supplied 'value' matches the number of 'enrichments(x)', colors are
#'    assigned in the order supplied.
#'    Otherwise, 'names(value)' must be present in 'enrichments(x)',
#'    and values are assigned in matching order.
#' @export
setMethod("colorV<-", "Mem",
   function(x, value) {
      if (length(colorV) == length(x@colorV)) {
         x@colorV[] <- colorV;
         names(x@colorV) <- enrichments(x)
      } else if (length(names(colorV)) > 0 &&
            all(names(colorV) %in% enrichments(x))) {
         x@colorV[names(colorV)] <- colorV;
      } else {
         stopmsg <- paste0("colorV must be correct length, or ",
            "have names(colorV) that match enrichments(x)")
         stop(stopmsg)
      }
      # confirm names(colorV) match enrichments(x)
      names(x@colorV) <- enrichments(x)
      x
   }
)

#' List thresholds defined in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List thresholds defined
#' @export
setMethod("thresholds", "Mem",
   function(x, ...) x@thresholds
)

#' Assign thresholds defined in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class thresholds defined in a Mem object.
#'    Note that thresholds do not trigger any other updates to the
#'    Mem object.
#' @export
setMethod("thresholds<-", "Mem",
   function(x, value) {
      if (!is.list(value)) {
         if (is.atomic(value)) {
            if (length(names(value)) == 0) {
               stop("value must have names.");
            }
            value <- as.list(value)
         }
      }
      # Add some checks for required threshold names?
      x@thresholds <- value;
      x
   }
)

#' Describe dimensions of a Mem object
#'
#' @param x `Mem` object
#' @docType methods
#' @describeIn Mem-class dimensions in order of genes, sets, and enrichments
#' @export
setMethod("dim", "Mem",
   function(x) c(nrow(x@memIM), ncol(x@memIM), ncol(x@enrichIM))
)

#' Dimension names for a Mem object
#'
#' @param x `Mem` object
#' @docType methods
#' @describeIn Mem-class dimension names for each genes, sets, and enrichments
#' @export
setMethod("dimnames", "Mem",
   function(x) 
      list(genes=rownames(x@memIM),
         sets=colnames(x@memIM),
         enrichments=colnames(x@enrichIM))
)

#' @name as
#' @rdname Mem-class
#' @aliases coerce
setAs("list", "Mem",
   function(from) {
      list_to_Mem(from)
   })

#' @name as
#' @rdname Mem-class
#' @aliases coerce
setAs("Mem", "list",
   function(from) {
      Mem_to_list(from)
   })

#' @returns `list_to_Mem()` returns a `Mem` S4 object, from 'list' or
#'    'Mem' input.
#'
#' @describeIn Mem-class Convert legacy `list` mem format to S4 `Mem`
#' @param mem `list` output from `multiEnrichMap()`
#' 
#' @examples
#' # list_to_Mem examples
#' data(Memtest)
#' mem <- as(Memtest, "list")
#' Mem <- list_to_Mem(mem)
#'
#' @export
list_to_Mem <- function
(mem)
{
   #
   use_slots=c(
      enrichList="list",
      enrichLabels="character",
      colorV="character",
      geneHitList="list",
      geneHitIM="matrix",
      memIM="matrix",
      geneIM="matrix",
      enrichIM="matrix",
      multiEnrichDF="data.frame",
      multiEnrichResult="enrichResult",
      thresholds="list",
      headers="list",
      enrichIMcolors="matrix",
      enrichIMdirection="matrix",
      enrichIMgeneCount="matrix",
      geneIMcolors="matrix",
      geneIMdirection="matrix"
   )
   # return Mem unchanged if it is already 'Mem'
   if (inherits(mem, "Mem")) return(mem);
   # rename 'colnames' to 'headers'
   hmatch <- match("colnames", names(mem));
   if (length(hmatch) == 1) {
      names(mem)[hmatch] <- "headers";
   }
   
   # Some slots can be added with defaults
   if (!"geneIMdirection" %in% names(mem)) {
      mem$geneIMdirection <- (abs(mem$geneIM) > 0) * 1;
   }
   
   if (!all(names(use_slots) %in% names(mem))) {
      missing_slots <- setdiff(names(use_slots), names(mem));
      stop(paste0(
         "Not all required slot names were provided:\n",
         jamba::cPaste(missing_slots, sep=", ")))
   }
   # print(jamba::sdim(x[names(use_slots)]));# debug
   Mem <- new("Mem",
      enrichList=mem$enrichList,
      enrichLabels=mem$enrichLabels,
      colorV=mem$colorV,
      geneHitList=mem$geneHitList,
      geneHitIM=mem$geneHitIM,
      memIM=mem$memIM,
      geneIM=mem$geneIM,
      enrichIM=mem$enrichIM,
      multiEnrichDF=mem$multiEnrichDF,
      multiEnrichResult=mem$multiEnrichResult,
      thresholds=mem$thresholds,
      headers=mem$headers,
      enrichIMcolors=mem$enrichIMcolors,
      enrichIMdirection=mem$enrichIMdirection,
      enrichIMgeneCount=mem$enrichIMgeneCount,
      geneIMcolors=mem$geneIMcolors,
      geneIMdirection=mem$geneIMdirection
   )
   return(Mem);
}


#' Mem updateObject
#' 
#' @returns `updateObject()` returns a `Mem` S4 object with the current
#'    class version, and current valid class definition. It should be
#'    able to update any serialized 'Mem' object from previous
#'    class versions.
#'
#' @describeIn Mem-class Confirm Mem object meets current class definition.
#' @export
setMethod("updateObject",
   methods::signature(object="Mem"),
   definition=function(object, ..., verbose=FALSE) {
      # be open about what's happening
      if (verbose) message("updateObject(object = 'Mem')");
      
      object <- methods::callNextMethod()
      
      # if current, return as-is
      if (Biobase::isCurrent(object)["Mem"]) return(object)
      
      ## Create an updated instance.
      if (!Biobase::isVersioned(object)) {
         if (verbose) message("Re-creating Mem to include Versioned.")
         ## Radical surgery – create a new, up-to-date instance
         object <- new("Mem",
            enrichList=object@enrichList,
            enrichLabels=object@enrichLabels,
            colorV=object@colorV,
            geneHitList=object@geneHitList,
            geneHitIM=object@geneHitIM,
            memIM=object@memIM,
            geneIM=object@geneIM,
            geneIMcolors=object@geneIMcolors,
            geneIMdirection=object@geneIMdirection,
            enrichIM=object@enrichIM,
            enrichIMcolors=object@enrichIMcolors,
            enrichIMdirection=object@enrichIMdirection,
            enrichIMgeneCount=object@enrichIMgeneCount,
            headers=object@headers,
            thresholds=object@thresholds,
            multiEnrichDF=object@multiEnrichDF,
            multiEnrichResult=object@multiEnrichResult
            # annotation = updateObject(annotation(object), ..., verbose=verbose)
         )
         use_version <- Biobase::classVersion(object="Mem")["Mem"];
         Biobase::classVersion(object) <- use_version;
      } else if (Biobase::classVersion(object)["Mem"] < "1.0.0") {
         # Do something for specific versions as necessary
         # For now, do nothing.
      } else {
         ## Make minor changes, and update version by consulting class definition
         if (verbose) message ("Updating Mem class version to current.")
         use_version <- Biobase::classVersion(object="Mem")["Mem"];
         Biobase::classVersion(object)["Mem"] <- use_version;
      }
      object
   })


#' Genes in each gene set (category)
#'
#' @param x `Mem` object
#' @docType methods
#' @describeIn Mem-class Genes in each pathway set (category)
#' @returns `geneInCategory()` returns a `list` named by pathway, containing
#'    `character` vectors with genes in each pathway.
#' @export
setMethod("geneInCategory", "Mem", function(x) im2list(memIM(x)))


#' Pathway gene sets associated with each gene
#'
#' @param x `Mem` object
#' @docType methods
#' @describeIn Mem-class Pathways gene sets associated with each gene
#' @returns `geneInCategory()` returns a `list` named by gene,
#'    containing `character` vectors with associated gene sets.
#' @export
setMethod("setsByGene", "Mem", function(x) im2list(t(memIM(x))))
