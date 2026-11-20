
# Design idea:
# Include pathway rank per enrichment,
# as a way of showing the top N for each enrichment

#' Return the Mem enrichment matrix using column ranks
#' 
#' Return the Mem enrichment matrix using column ranks,
#' applying thresholds 'topEnrichN', 'min_count', and 'p_cutoff'.
#' 
#' Ties are broken first by using `enrichIMgeneCount` then by
#' using the current row and column order.
#' 
#' @family jam Mem utilities
#' 
#' @param x `Mem` object
#' @param topEnrichN,p_cutoff,min_count `numeric` values obtained
#'    from `thresholds(x)` by default, however they can be customized
#'    to determine updated ranks using the custom thresholds.
#' @param ... additional arguments are ignored.
#' 
#' @returns `numeric` matrix of `numeric` ranks by column.
#'    There should be no ties, all values should be `integer`.
#' 
#' @examples
#' enrichIMrank(Memtest, min_count=2)
#' 
#' enrichIMrank(Memtest, topEnrichN=20)
#' @export
enrichIMrank <- function(
   x,
   topEnrichN=NULL,
   p_cutoff=NULL,
   min_count=NULL,
   ...
) {
   # Determine thresholds using x by default
   x_thresholds <- thresholds(x);
   if (length(p_cutoff) == 0) {
      p_cutoff <- x_thresholds$p_cutoff;
      if (length(p_cutoff) == 0) {
         p_cutoff <- x_thresholds$cutoffRowMinP;
      }
   }
   if (length(topEnrichN) == 0) {
      topEnrichN <- x_thresholds$topEnrichN
   }
   if (length(min_count) == 0) {
      min_count <- x_thresholds$min_count;
   }

   # calculate column rank after applying thresholds,
   # and breaking ties using enrichIMgeneCount then current order
   use_enrichIM <- jamba::noiseFloor(
      enrichIM(x),
      ceiling=p_cutoff,
      newCeiling=NA);
   use_enrichIMgeneCount <- (jamba::noiseFloor(
      enrichIMgeneCount(x),
      minimum=min_count,
      newValue=NA
   ) > 0) * 1;
   enrichIMtoRank <- (
      -log10(use_enrichIM) * (use_enrichIMgeneCount) +
         enrichIMgeneCount(x) * 1e-10 +
         seq_len(prod(dim(enrichIM(x)))) * 1e-15)
   enrichIMrank <- apply(enrichIMtoRank, 2, function(i){
      irank <- rank(-i);
      irank[irank > topEnrichN | is.na(i)] <- NA;
      irank
   })
   enrichIMrank
}

