
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
      "colnames")
   test0 <- all(required_slots %in% slotNames(object))
   if (!test0) {
      print(paste0("Required slots were not present: ",
         jamba::cPaste(setdiff(required_slots, slotNames(object)), sep=", ")))
      return(test0)
   }
   criteria_set <- c(
      `memIM, enrichIM`=all(colnames(object@memIM) %in% rownames(object@enrichIM)),
      `memIM, geneIM`=all(rownames(object@memIM) %in% rownames(object@geneIM)),
      `enrichList, enrichIM`=all(names(object@enrichList) %in% colnames(object@enrichIM)),
      `enrichList, geneIM`=all(names(object@enrichList) %in% colnames(object@geneIM)),
      `enrichList, colorV`=all(names(object@enrichList) %in% names(object@colorV)),
      `enrichList, geneHitList`=all(names(object@enrichList) %in% names(object@geneHitList)),
      `enrichList, geneHitIM`=all(names(object@enrichList) %in% colnames(object@geneHitIM))
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
#' * **colnames**: `list` to associate actual colnames in 'enrichList'
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
      colnames="list",
      enrichIMcolors="matrix",
      enrichIMdirection="matrix",
      enrichIMgeneCount="matrix",
      geneIMcolors="matrix",
      geneIMdirection="matrix"
   )
);
setValidity("Mem", method=check_Mem)

#' Subset a Mem object
#'
#' @param x `Mem` object
#' @param i `character` or `integer` with gene set / pathway names
#' @param j `character` or `integer` with enrichment names
#' @param ... additional arguments are ignored for '['.
#' @docType methods
#' @describeIn Mem-class Subset a Mem object
#' @export
setMethod("[",
   signature=c(x="Mem",
      i="ANY",
      j="ANY"),
   definition=function(x, i=NULL, j=NULL, ...) {
      if (missing(i)) {
         i <- NULL;
      }
      if (missing(j)) {
         j <- NULL;
      }
      ivals <- rownames(x@enrichIM);
      jvals <- colnames(x@enrichIM);
      # row subset by set name
      if (length(i) > 0) {
         if (is.numeric(i)) {
            if (any(abs(i) > length(ivals))) {
               stop("i is out of bounds, does not match set indices.");
            }
            i <- ivals[i];
            # x@polygons <- x@polygons[i, , drop=FALSE]
         } else if (inherits(i, c("character", "factor"))) {
            i <- as.character(i);
            if (!all(i %in% ivals)) {
               stop("i is out of bounds, does not match character set names.")
            }
            # imatch <- match(as.character(i), ivals);
            # x@polygons <- x@polygons[imatch, , drop=FALSE]
         }
         if (length(i) > 0 && !all(i == ivals)) {
            islot_cols <- c("memIM")
            islot_rows <- c("enrichIM",
               "enrichIMcolors",
               "enrichIMdirection",
               "enrichIMgeneCount",
               "geneHitIM")
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
      # column subset by enrichment name
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
            jslot_names <- c("enrichList",
               "enrichLabels",
               "colorV",
               "geneHitList")
            jslot_cols <- c(
               "enrichIM",
               "enrichIMcolors",
               "enrichIMdirection",
               "enrichIMgeneCount",
               "geneHitIM",
               "geneIM",
               "geneIMcolors",
               "geneIMdirection")
            jslot_rows <- NULL;
            # subset by names
            for (jslot in jslot_names) {
               jmatch <- match(j, names(slot(x, jslot)))
               slot(x, jslot) <- slot(x, jslot)[jmatch];
            }
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
            # Update memIM with new contents in enrichIM, geneIM
            geneColname <- x@colnames[["geneColname"]];
            nameColname <- x@colnames[["nameColname"]];
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
            jmatchrow <- match(rownames(new_memIM), rownames(x@memIM));
            jmatchcol <- match(colnames(new_memIM), colnames(x@memIM));
            x@memIM[jmatchrow, jmatchcol] <- new_memIM;
            #
            # Todo:
            # Subset memIM columns to remove empty pathways
         }
      }
      x;
   }
)

#' Memtest S4 Mem object for testing in multienrichjam
#'
#' Mem object for testing in multienrichjam, with 'Newborns' and
#' 'OlderChildren' described in
#' Reese et al 2019 [https://doi.org/10.1016/j.jaci.2018.11.043].
#'
#' @family Mem
#' @format `Mem` object containing two enrichment tests
#' @rdname Memtest
"Memtest"

#' @describeIn Mem-class Convert list to Mem S4 object
#'
#' @returns `list_to_Mem()` returns a `Mem` S4 object.
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
      colnames="list",
      enrichIMcolors="matrix",
      enrichIMdirection="matrix",
      enrichIMgeneCount="matrix",
      geneIMcolors="matrix",
      geneIMdirection="matrix"
   )

   # Some slots can be added with defaults
   if (!"geneIMdirection" %in% names(mem)) {
      x$geneIMdirection <- (abs(x$geneIM) > 0) * 1;
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
      colnames=mem$colnames,
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
      if ("cutoffRowMinP" %in% names(object@thresholds)) {
         p_cutoff <- object@thresholds$cutoffRowMinP;
      } else {
         p_cutoff <- object@thresholds$p_cutoff;
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
            jamba::cPaste(object@colnames$pvalueColname, sep=", "),
            ")"),
         paste0("- min gene count: ", min_count)
      )
      if ("directionColname" %in% names(object@colnames)) {
         summary_v <- c(summary_v,
            paste0("- direction colname: ",
               object@colnames$directionColname))
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
         colnames="list",
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
      new_list
   })

#' @describeIn Mem-class Coerce `Mem` object to legacy 'mem' `list` format
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


setGeneric("geneIM", signature="x",
   function(x, ...) standardGeneric("geneIM")
)

#' The matrix of genes tested versus enrichment tests from a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class The matrix of genes tested versus enrichment tests
#' @export
setMethod("geneIM", "Mem",
   function(x, ...) x@geneIM
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

setGeneric("sets", signature="x",
   function(x, ...) standardGeneric("sets")
)

#' List gene sets represented in a Mem object
#'
#' @param x `Mem` object
#' @param ... additional arguments are ignored
#' @docType methods
#' @describeIn Mem-class List gene sets represented
#' @export
setMethod("sets", "Mem",
   function(x, ...) rownames(x@enrichIM)
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
