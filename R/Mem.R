
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
#' * enrichList
#' * enrichLabels
#' * colorV
#' * geneHitList
#' * geneHitIM
#' * memIM
#' * geneIM
#' * enrichIM
#' * multiEnrichDF
#' * multiEnrichResult
#' * thresholds
#' * colnames
#'
#' * geneIMcolors
#' * geneIMdirection
#'
#' * enrichIMcolors
#' * enrichIMdirection
#' * enrichIMgeneCount
#'
#' @family Mem
#'
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

#' Subset Mem object
#'
#' @docType methods
#' @rdname Mem-methods
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

#' Coerce legacy multiEnrichMap list format to Mem object
#'
#' Coerce legacy multiEnrichMap list format to `Mem` S4 object
#'
#' @family Mem
#'
#' @param x `list` output from `multiEnrichMap()`
#' @param ... additional arguments are ignored
#'
#' @examples
#' jamsession::load_object("test_mem")
#' test_mem$thresholds <- list(p_cutoff=0.05, min_count=3)
#' Mem <- list_to_Mem(test_mem)
#'
#' @export
list_to_Mem <- function
(x,
 ...)
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
   if (!"geneIMdirection" %in% names(x)) {
      x$geneIMdirection <- (abs(x$geneIM) > 0) * 1;
   }

   if (!all(names(use_slots) %in% names(x))) {
      missing_slots <- setdiff(names(use_slots), names(x));
      stop(paste0(
         "Not all required slot names were provided:\n",
         jamba::cPaste(missing_slots, sep=", ")))
   }
   # print(jamba::sdim(x[names(use_slots)]));# debug
   Mem <- new("Mem",
      enrichList=x$enrichList,
      enrichLabels=x$enrichLabels,
      colorV=x$colorV,
      geneHitList=x$geneHitList,
      geneHitIM=x$geneHitIM,
      memIM=x$memIM,
      geneIM=x$geneIM,
      enrichIM=x$enrichIM,
      multiEnrichDF=x$multiEnrichDF,
      multiEnrichResult=x$multiEnrichResult,
      thresholds=x$thresholds,
      colnames=x$colnames,
      enrichIMcolors=x$enrichIMcolors,
      enrichIMdirection=x$enrichIMdirection,
      enrichIMgeneCount=x$enrichIMgeneCount,
      geneIMcolors=x$geneIMcolors,
      geneIMdirection=x$geneIMdirection
   )
   return(Mem);
}

if (!isGeneric("summary")) {
   setGeneric("summary",
      function(object, ...) standardGeneric("summary"))
}

#' Extract some elements of potentially long vector for display
#'
#' Extract some elements of potentially long vector for display
#'
#' @docType methods
#' @rdname Mem-methods
#' @export
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

setMethod("names", "Mem", function(x) names(x@enrichList))


setGeneric("geneIM", signature="x",
   function(x, ...) standardGeneric("geneIM")
)
setMethod("geneIM", "Mem",
   function(x, ...) x@geneIM
)

setGeneric("enrichIM", signature="x",
   function(x, ...) standardGeneric("enrichIM")
)
setMethod("enrichIM", "Mem",
   function(x, ...) x@enrichIM
)

setGeneric("memIM", signature="x",
   function(x, ...) standardGeneric("memIM")
)
setMethod("memIM", "Mem",
   function(x, ...) x@memIM
)
