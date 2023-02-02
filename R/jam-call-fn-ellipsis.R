
#' Call function using safe ellipsis arguments
#'
#' Call function using safe ellipsis arguments
#'
#' This function is deprecated, instead please use `jamba::call_fn_ellipsis()`.
#'
#' This function is a wrapper function intended to help
#' pass ellipsis arguments `...` from a parent function
#' to an external function in a safe way. It will only
#' include arguments from `...` that are recognized by
#' the external function.
#'
#' When the external function `FUN` arguments `formals()` includes
#' ellipsis `...`, then the `...` will be passed as-is without
#' change.
#'
#' When the external function `FUN` arguments `formals()` does not
#' include ellipsis `...`, then only named arguments in `...` that
#' are recognized by `FUN` will be passed, as defined by
#' `names(formals(FUN))`.
#'
#' Note that arguments must be named.
#'
#' @family jam utility functions
#'
#' @return output from `FUN()` when called with relevant named arguments
#'    from ellipsis `...`
#'
#' @param FUN `function` that should be called with arguments in `...`
#' @param ... arguments are passed to `FUN()` in safe manner.
#'
#' @examples
#' new_mean <- function(x, trim=0, na.rm=FALSE) {
#'    mean(x, trim=trim, na.rm=na.rm)
#' }
#' x <- c(1, 3, 5, NA);
#' new_mean(x, na.rm=TRUE);
#' # throws an error as expected (below)
#' # new_mean(x, na.rm=TRUE, color="red");
#'
#' call_fn_ellipsis_deprecated(new_mean, x=x, na.rm=TRUE, color="red")
#' # throws an error as expected (below)
#' # call_fn_ellipsis_deprecated(new_mean, x=x, color="red")
#'
#' @export
call_fn_ellipsis_deprecated <- function
(FUN,
 ...)
{
   FUN_argnames <- names(formals(FUN));
   if ("..." %in% FUN_argnames) {
      FUN(...)
   } else {
      arglist <- list(...)
      argkeep <- which(names(arglist) %in% FUN_argnames);
      arguse <- arglist[argkeep]
      do.call(FUN, arguse)
   }
}
