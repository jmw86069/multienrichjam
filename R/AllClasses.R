# AllClasses.R


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
#' * The aim is to support all user-facing tasks using S4 methods,
#' to prevent the need to access slots directly.
#'
#' @importClassesFrom DOSE enrichResult
#' @importClassesFrom Biobase Versioned
#' @importFrom BiocGenerics updateObject
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
         names(object@enrichList) %in% colnames(object@geneHitIM)),
      `geneIM, geneHitIM rownames`=all(
         rownames(object@geneIM) %in% rownames(object@geneHitIM))
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


#' Mem S4 class, accessors, getters, and setters
#'
#' Mem class containing results from `multiEnrichMap()`,
#' current class version "1.0.0".
#'
#' @family Mem
#' @family jam Mem utilities
#' @rdname Mem-class
#' @aliases Mem
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
   ),
   contains="Versioned",
   prototype=prototype(
      new("Versioned",
         versions=c(Mem="1.0.0")))
);
setValidity("Mem", method=check_Mem)



#' Mem description of slots
#'
#' Mem description of slots, for reference, slots may change at any time
#' 
#' Note that the Mem S4 object is intended to be used via the accessor
#' functions as described in `Mem-class`.
#' 
#' However, for reference, and in case *really necessary* the slots are
#' described in more detail below. These slots may change at any time,
#' as the 'Mem' S4 object evolves after more usage, so relying upon
#' specific slots to be stable may not be ideal.
#' 
#' The descriptions below should be edited whenever the 'Mem' S4 slot
#' structure changes, but be prepared there may sometimes be some delay
#' even under the best of intentions.
#' 
#' # Enrichment data input
#' 
#' * **enrichList**: `list` of `enrichResults` named by enrichment test.
#'    This list represents the input enrichment data.
#' * **enrichLabels**: `character` vector with custom enrichment labels,
#'    named by enrichment test.
#' * **geneHitList**: `list` of vectors named by enrichment test.
#' Vectors are either:<br>
#' `character` vector of genes tested as hits for enrichment, or<br>
#' `integer` vector named by genes with values indicating
#' the directionality, using `1` or `-1`.
#' * **geneHitIM**: `numeric` matrix with enrichment tests as colnames, and
#' genes as rownames.<br>
#' Values are `1` or `0` to indicate gene hits tested,<br>
#' optionally `1` or `-1` to indicate directionality.
#' 
#' # Incidence matrix data
#' 
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
#' 
#' # Combined enrichment data
#' 
#' Enrichment data are combined into singular objects, mainly intended
#' to be used with `clusterProfiler` related functions.
#' Also the `data.frame` could be used as a supplemental table.
#' 
#' * **multiEnrichDF**: `data.frame` with single combined enrichment table.<br>
#' Data include union of genes involved in each pathway, and
#' lowest P-value observed across all enrichment tests.
#' * **multiEnrichResult**: `enrichResult` equivalent of 'multiEnrichDF'.<br>
#' It is intended to facilitate use in `enrichplot` functions.
#' 
#' ## Analysis details
#' 
#' * **colorV**: `character` vector of R colors, named by enrichment test.
#' * **thresholds**: `list` of thresholds used when calling `multiEnrichMap()`.
#' 
#'    * **p_cutoff**: `numeric` P-value threshold for significant enrichment.
#'    * **min_count**: `integer` number of genes involved in a pathway
#'    for the pathway to be considered significant.
#' 
#' * **headers**: `list` to associate actual colnames (headers) in 'enrichList'
#' with the conceptual names used in multienrichjam.
#' 
#' ## Additional incidence matrix data
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
#' 
#' @name Mem-slots
NULL

