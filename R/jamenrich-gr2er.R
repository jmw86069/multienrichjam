
#' Convert gseaResult to enrichResult equivalent result content
#' 
#' Convert gseaResult to enrichResult equivalent result content,
#' focused mainly on providing minimal result columns for downstream
#' work.
#' 
#' This function is intended for internal use, and is not a full
#' conversion of `gseaResult` to a proper `enrichResult`.
#' It simply makes minor modifications to colnames for convenience
#' within `multiEnrichMap()` and other supporting functions.
#' 
#' Specific changes:
#' 
#' * Renames `'core_enrichment'` to `'geneID'` if `'geneID'`
#' does not exist.
#' * Uses `'geneID'` to create `'Count'` column with number of
#' genes that contribute to enrichment.
#' * Adds new attribute `'genes_signed'` which contains a subset
#' of the gene rank scores stored in attribute `'geneList'`,
#' only including those genes that contribute to enrichment,
#' which is defined by the presence of a gene in the column
#' `'core_enrichment'`.
#' 
#' @returns `gseaResult`
#' 
#' @param er `gseaResult`
#' @param ... additional arguments are ignored.
#' 
#' @keywords internal
gr2er <- function(
   er,
   ...
) {
	#
	if (inherits(er, "enrichResult")) {
		return(er)
	}
	if (!inherits(er, "gseaResult")) {
		stop("Input must be gseaResult or enrichResult.");
	}
	if (!"geneID" %in% colnames(er@result)) {
		if (!"core_enrichment" %in% colnames(er@result)) {
			cli::cli_abort(c(
            "Input {.var er@result} must contain ",
            "{.val core_enrichment} or {.val geneID}"
         ))
		}
		er@result <- jamba::renameColumn(er@result,
			from=c("core_enrichment"),
			to="geneID")
	}
   # Split delimited genes into vector of values
   er_geneID <- strsplit(er@result$geneID, "[/]");
	if (!"Count" %in% colnames(er@result)) {
		er@result$Count <- lengths(er_geneID);
	}

   # Check for 'geneList' with rank scores
   if ("geneList" %in% names(attributes(er))) {
      used_geneList <- attr(er, "geneList");
      kept_geneID <- intersect(
         names(used_geneList),
         unique(unlist(er_geneID))
      )
      used_geneList <- used_geneList[kept_geneID]
      attr(er, "genes_signed") <- used_geneList;
   }

   # Add direction as a column if not present
	if (!"direction" %in% colnames(er@result)) {
		if ("NES" %in% colnames(er@result)) {
			er@result$direction <- er@result$NES;
		} else {
         cli::cli_abort(c(
            "Input {.var er@result} must contain ",
            "{.val NES} or {.val directioin}"
         ))
      }
	}
	return(er)
}
