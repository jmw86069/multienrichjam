
#' Get Mem input data
#' 
#' Get Mem data for use in multienrichjam functions, intended for
#' compatibility when transitioning from 'list' mem to S4 'Mem' object
#' as the primary data format.
#' 
#' @keywords internal
#' @noRd
get_mem <- function
(Mem,
 output=c("Mem", "mem"),
 ...)
{
   #
   output <- match.arg(output);
   if ("Mem" %in% output) {
      if (inherits(Mem, "Mem")) {
         return(Mem)
      }
      return(list_to_Mem(Mem))
   }
   if ("mem" %in% output) {
      if (is.list(Mem)) {
         return(Mem);
         # optionally validate the data contents
         # - consider converting to Mem and back?
         if (all(slotNames(getClass("Mem"))) %in% names(Mem)) {
            return(Mem)
         }
      } else if (inherits(Mem, "Mem")) {
         return(Mem_to_list(Mem))
      } else {
         stop("Unrecognized input.")
      }
   }
   return(Mem)
}
