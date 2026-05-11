
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
#' * Renames 'core_enrichment' to 'geneID' if 'geneID' does not exist.
#' * Uses 'geneID' to create 'Count' column with number of genes.
#' 
#' @returns `gseaResult`
#' 
#' @param er `gseaResult`
#' @param ... additional arguments are ignored.
#' 
#' @keywords internal
gr2er <- function
(er,
 ...)
{
	#
	if (inherits(er, "enrichResult")) {
		return(er)
	}
	if (!inherits(er, "gseaResult")) {
		stop("Input must be gseaResult or enrichResult.");
	}
	if (!"geneID" %in% colnames(er@result)) {
		if (!"core_enrichment" %in% colnames(er@result)) {
			stop("Input er@result must contain 'core_enrichment' or 'geneID'");
		}
		er@result <- jamba::renameColumn(er@result,
			from=c("core_enrichment"),
			to="geneID")
	}
	if (!"Count" %in% colnames(er@result)) {
		gc <- lengths(strsplit(er@result$geneID, "[/]"));
		er@result$Count <- gc;
	}
	if (!"direction" %in% colnames(er@result)) {
		if ("NES" %in% colnames(er@result)) {
			er@result$direction <- er@result$NES;
		}
	}
	return(er)
}
