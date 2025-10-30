
# Mem object


#' Check Mem object
#'
#' Check whether a Mem object is valid.
#'
#' It requires:
#'
#' * `slotName(object)`: `"memIM", "enrichIM", "geneIM", "enrichList"`
#' * `object@memIM` to inherit `"matrix"`
#' * `object@enrichIM` to inherit `"matrix"`
#' * `object@geneIM` to inherit `"matrix"`
#' * `colnames(object@memIM)` to equal `rownames(object@enrichIM)`
#' * `rownames(object@memIM)` to equal `colnames(object@geneIM)`
#' * `names(enrichList)` to equal `colnames(object@enrichIM)`
#'
#' General guidance for Mem objects:
#'
#' * TBD
#'
#' @importClassesFrom DOSE enrichResult
#' 
#' @param object `Mem` object
#'
#' @family Mem
#'
#' @export
check_Mem <- function
(object)
{
   # check for valid content
   required_slots <- c(
      "enrichList",
      "enrichLabels",
      "colorV",
      "geneHitList",
      "geneHitIM",
      "memIM",
      "geneIM",
      "geneIMcolors",      # optional?
      "geneIMdirection",   # optional?
      "enrichIM",
      "enrichIMcolors",    # optional?
      "enrichIMdirection", # optional?
      "enrichIMgeneCount", # optional?
      "multiEnrichDF",
      "multiEnrichResult",
      "thresholds",
      "headers")
   test0 <- all(required_slots %in% slotNames(object))
   if (!test0) {
      print(paste0("Required slots were not present: ",
         jamba::cPaste(setdiff(required_slots, slotNames(object)), sep=", ")))
      return(test0)
   }
   # confirm rownames and colnames match across objects
   criteria_set <- c(
      `memIM, enrichIM`=all(
         colnames(object@memIM) %in% rownames(object@enrichIM)),
      `memIM, geneIM`=all(
         rownames(object@memIM) %in% rownames(object@geneIM)),
      `enrichList, enrichIM`=all(
         names(object@enrichList) %in% colnames(object@enrichIM)),
      `enrichList, geneIM`=all(
         names(object@enrichList) %in% colnames(object@geneIM)),
      `enrichList, colorV names`=all(
         names(object@enrichList) %in% names(object@colorV)),
      `enrichList, geneHitList names`=all(
         names(object@enrichList) %in% names(object@geneHitList)),
      `enrichList, geneHitIM`=all(
         names(object@enrichList) %in% colnames(object@geneHitIM))
   )
   if (!all(criteria_set)) {
      jamba::printDebug(
         "Mem Invalid:\n- ",
         jamba::cPaste(names(criteria_set)[!criteria_set],
            sep=",\n- "));
      return(FALSE);
   }
   return(TRUE)

   # Todo:
   # - require colnames be present in all enrichList entries
   # - require colnames be present in multiEnrichDF
   # - require colnames be present in multiEnrichResult@result
}


#' Mem class
#'
#' Mem class containing results from `multiEnrichMap()`.
#'
#' ## Slots
#'
#' * **enrichList**: `list` of `enrichResults` named by enrichment test.
#' * **enrichLabels**: `character` vector with custom enrichment labels,
#' named by enrichment test.
#' * **colorV**: `character` vector of R colors, named by enrichment test.
#' * **geneHitList**: `list` of vectors named by enrichment test.
#' Vectors are either:<br>
#' `character` vector of genes tested as hits for enrichment, or<br>
#' `integer` vector named by genes with values indicating
#' the directionality, using `1` or `-1`.
#' * **geneHitIM**: `numeric` matrix with enrichment tests as colnames, and
#' genes as rownames.<br>
#' Values are `1` or `0` to indicate gene hits tested,<br>
#' optionally `1` or `-1` to indicate directionality.
#' * **memIM**: `numeric` matrix with gene sets / pathways as colnames, and
#' genes as rownames.<br>
#' Values are `1` or `0` indicating genes present in each pathway.
#' * **geneIM**: `numeric` matrix with enrichment test as colnames, and genes
#' as rownames.<br>
#' This matrix differs from 'geneHitIM' in that it only includes
#' genes also in 'memIM', and in that order.
#' * **enrichIM**: `numeric` matrix with enrichment test as colnames, and
#' gene sets / pathways as rownames.<br>
#' Matrix values are enrichment P-values,
#' using the appropriate adjusted P-value or FDR to represent significance.
#' * **multiEnrichDF**: `data.frame` with single combined enrichment table.<br>
#' Data include union of genes involved in each pathway, and
#' lowest P-value observed across all enrichment tests.
#' * **multiEnrichResult**: `enrichResult` equivalent of 'multiEnrichDF'.<br>
#' It is intended to facilitate use in `enrichplot` functions.
#' * **thresholds**: `list` of thresholds used when calling `multiEnrichMap()`.
#' * **headers**: `list` to associate actual colnames (headers) in 'enrichList'
#' with the conceptual names used in multienrichjam.
#'
#' * **geneIMcolors**: `character` matrix with enrichment test colnames, and
#' gene rownames.<br>
#' Values are colors assigned using 'colorV'.
#' * **geneIMdirection**: `numeric` matrix equivalent to 'geneIM' with values
#' `1` and `-1` indicating directionality.<br>
#' When no directionality is provided, all values will be `1` or `0`.
#'
#' * **enrichIMcolors**: `character` matrix with enrichment test colnames, and
#' gene set / pathway rownames.<br>
#' Values are color gradients using the
#' appropriate P-value threshold used with `multiEnrichMap()`.
#' * **enrichIMdirection**: `numeric` matrix equivalent to 'enrichIM' with values
#' indicating directional score.<br>
#' Directional scores may include a IPA z-score or other directional
#' score where available.<br>
#' When not available, all values are `1`.
#' * **enrichIMgeneCount**: `integer` matrix equivalent to 'enrichIM'.<br>
#' Values are the number of genes involved in the enrichment of each pathway.
#'
#' @family Mem
#' @rdname Mem-class
#' @aliases Mem
#' @examples
#' # See vignettes for a full walk-through
#'
setClass("Mem",
   slots=c(
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
);
setValidity("Mem", method=check_Mem)

# setGeneric("[", function(x, i, j, ..., drop = TRUE) standardGeneric("["))

#' Subset a Mem object
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
      
      ivals <- rownames(x@memIM);
      jvals <- rownames(x@enrichIM);
      kvals <- colnames(x@enrichIM);
      
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
         if (length(i) > 0 && !all(i == ivals)) {
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
         if (length(j) > 0 && !all(j == jvals)) {
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
            jamba::printDebug("character k:");print(k);# debug
         } else if (inherits(k, c("character", "factor"))) {
            k <- as.character(k);
            if (!all(k %in% kvals)) {
               stop("k is out of bounds, does not match character enrichment names.")
            }
         }
         if (length(k) > 0 && !all(k == kvals)) {
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
            new_memIM <- list2im(
               keepCounts=TRUE,
               strsplit(
                  jamba::nameVector(
                     as.character(x@multiEnrichDF[[geneColname]]),
                     x@multiEnrichDF[[nameColname]]),
                  geneDelim));
            # empty the existing memIM since we cannot yet subset the columns
            # instead populate with columns that remain,
            # then consider subsetting columns afterward
            x@memIM <- x@memIM * 0;
            kmatchrow <- match(rownames(new_memIM), rownames(x@memIM));
            kmatchcol <- match(colnames(new_memIM), colnames(x@memIM));
            x@memIM[kmatchrow, kmatchcol] <- new_memIM;
            #
            # Todo:
            # Subset memIM columns to remove empty pathways
         }
      }
      x;
   }
)


#' @describeIn Mem-class Convert legacy `list` mem format to S4 `Mem`
#'
#' @returns `list_to_Mem()` returns a `Mem` S4 object, from 'list' or 'Mem' input.
#'
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


if (!isGeneric("summary")) {
   setGeneric("summary",
      function(object, ...) standardGeneric("summary"))
}

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
      if ("directionColname" %in% names(object@headers)) {
         summary_v <- c(summary_v,
            paste0("- direction colname: ",
               object@headers$directionColname))
      }
      summary_text <- jamba::cPaste(summary_v, sep="\n")
      cat(summary_text, sep="");
      # jamba::sdim(object)[, c("enrichIM", "geneIM", "memIM"), drop=FALSE]
   }
)

# setAs() to convert Mem back to mem using list format
# as(Mem, "list")

#' @name as
#' @rdname Mem-class
#' @aliases coerce
setAs("Mem", "list",
   function(from) {
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
         slot(from, slotname)
      })
      # rename 'headers' to 'colnames'
      hmatch <- match("headers", names(new_list));
      if (length(hmatch) == 1) {
         names(new_list)[hmatch] <- "colnames";
      }
      # et voila
      new_list
   })

#' @describeIn Mem-class Coerce S4 `Mem` to legacy `list` mem format
#' @aliases coerce
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#'
#' @returns `Mem_to_list()` returns a `list` suitable for ther mem functions.
#'
#' @examples
#' # examples for Mem_to_list and as(Mem, "list")
#' data(Memtest)
#' mem1 <- Mem_to_list(Memtest)
#' mem2 <- as(Memtest, "list")
#' identical(mem1, mem2)
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
   as(x, "list")
}



setGeneric("enrichments", signature="x",
   function(x, ...) standardGeneric("enrichments")
)

#' The enrichment names from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Names for the enrichment tests in a Mem object
#' @export
setMethod("enrichments", "Mem", function(x) names(x@enrichList))

#' The enrichment names from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Names for the enrichment tests in a Mem object
#' @export
setMethod("names", "Mem", function(x) names(x@enrichList))

setGeneric("enrichments<-", signature="x",
   function(x, value) standardGeneric("enrichments<-")
)

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

setGeneric("sets", signature="x",
   function(x, ...) standardGeneric("sets")
)

#' List pathway gene set names represented in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List pathway gene set names in a Mem object
#' @export
setMethod("sets", "Mem",
   function(x, ...) rownames(x@enrichIM)
)

setGeneric("sets<-", signature="x",
   function(x, value) standardGeneric("sets<-")
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

setGeneric("genes", signature="x",
   function(x, ...) standardGeneric("genes")
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

setGeneric("genes<-", signature="x",
   function(x, value) standardGeneric("genes<-")
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

setGeneric("geneIM", signature="x",
   function(x, ...) standardGeneric("geneIM")
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

setGeneric("geneIMdirection", signature="x",
   function(x, ...) standardGeneric("geneIMdirection")
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

setGeneric("geneIMcolors", signature="x",
   function(x, ...) standardGeneric("geneIMcolors")
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

setGeneric("enrichList", signature="x",
   function(x, ...) standardGeneric("enrichList")
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

setGeneric("enrichIM", signature="x",
   function(x, ...) standardGeneric("enrichIM")
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

setGeneric("enrichIMcolors", signature="x",
   function(x, ...) standardGeneric("enrichIMcolors")
)

#' The pathway-enrichment color matrix from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The pathway/directional enrichment color matrix
#' @export
setMethod("enrichIMcolors", "Mem",
   function(x, ...) x@enrichIMcolors
)

setGeneric("enrichIMdirection", signature="x",
   function(x, ...) standardGeneric("enrichIMdirection")
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

setGeneric("enrichIMgeneCount", signature="x",
   function(x, ...) standardGeneric("enrichIMgeneCount")
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

setGeneric("memIM", signature="x",
   function(x, ...) standardGeneric("memIM")
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

setGeneric("geneHitIM", signature="x",
   function(x, ...) standardGeneric("geneHitIM")
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

setGeneric("headers", signature="x",
   function(x, ...) standardGeneric("headers")
)

#' The column headers mapped in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class Enrichment data column headers associated to
#'    enrichResult data in a Mem object
#' @export
setMethod("headers", "Mem", function(x) names(x@headers))


setGeneric("colorV", signature="x",
   function(x, ...) standardGeneric("colorV")
)

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

setGeneric("colorV<-", signature="x",
   function(x, value) standardGeneric("colorV<-")
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


setGeneric("thresholds", signature="x",
   function(x, ...) standardGeneric("thresholds")
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


# future expansion:
#
# - cbind: combine two Mem objects, concatenating the enrichments,
#      and making new geneIM, enrichIM, memIM

# Consider: cbind() to combine two Mem objects, different enrichments
# Consider: rbind() to combine two Mem objects, same enrichments, new pathways
# - probably need to require identical colnames and thresholds

# setMethod("cbind", "Mem",
#    function(..., deparse.level=1) {
#       #
#       args <- unname(list(...))
#       .cbind.SummarizedExperiment(args)
#    })

# Accessors to consider:
# - names
# - geneIM
# - enrichIM
# - memIM
