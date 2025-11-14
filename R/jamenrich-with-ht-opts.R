
# withr helper functions

#' Withr mimic of with_options() for ComplexHeatmap options
#' 
#' Withr mimic of with_options() for ComplexHeatmap options
#' 
#' This function uses the same mechanism as used by `withr::with_options()`
#' except focuses on options defined by `ComplexHeatmap::ht_opt()`.
#' 
#' @param new `list` named by heatmap option, with corresponding values
#'    as described in `ComplexHeatmap::ht_opt()`.
#' @param code `any` valid R code to execute in the temporary environment.
#' 
#' @returns `any` The results of the evaluation of the `code` argument.
#' 
#' @family jam utility functions
#' @examples
#' ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#' #> 1mm
#' with_ht_opts(list(COLUMN_ANNO_PADDING=grid::unit(3, "mm")),
#'    ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"))
#' #> 3mm
#' ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#' #> 1mm
#' 
#' test_local <- function() {
#'    local_ht_opts(list(COLUMN_ANNO_PADDING=grid::unit(3, "mm")))
#'    print(ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"))
#' }
#' test_local()
#' #> 3mm
#' ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#' #> 1mm
#' @export
with_ht_opts <- function
(new,
 code)
{
   # define existing options
   old_ht_opt <- lapply(jamba::nameVectorN(new), function(iname){
      ComplexHeatmap::ht_opt(iname)
   })
   # define the reset
   reset_ht_opt <- function(old) {
      ComplexHeatmap::ht_opt(old)
   }
   on.exit(reset_ht_opt(old_ht_opt))
   
   # set new options
   ComplexHeatmap::ht_opt(new)
   
   # run the code
   force(code)
}

#' @rdname with_ht_opts
#' @export
local_ht_opts <- function
(.new=list(),
 ...,
 .local_envir=parent.frame())
{
   # combine .new and '...'
   lhs <- list(...)
   for (nme in names(lhs)) {
      .new[nme] <- lhs[nme];
   }
   
   # define existing options
   old_ht_opt <- lapply(jamba::nameVectorN(.new), function(iname){
      ComplexHeatmap::ht_opt(iname)
   })
   
   # define the reset
   reset_ht_opt <- function(old) {
      ComplexHeatmap::ht_opt(old)
   }
   # defer the reset
   withr::defer(reset_ht_opt(old_ht_opt), envir=.local_envir)
   
   # set new options
   ComplexHeatmap::ht_opt(.new)
}
