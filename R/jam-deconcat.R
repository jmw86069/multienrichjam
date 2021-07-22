
#' Deconcatenate delimited column values in a data.frame
#'
#' Deconcatenate delimited column values in a data.frame
#'
#' This function deconcatenates delimited values in a column
#' of a `data.frame` by calling `strsplit()` on column values,
#' and repeating values in all other columns to match `lengths()`
#' following `strsplit()`.
#'
#' This function includes a correction for cases where `strsplit()`
#' would otherwise return zero-length entries, and which would
#' otherwise be dropped from the output. From this function,
#' zero-length entries are replaced with `blank=""` so these
#' rows are not dropped from the output.
#'
#' @family jam utility functions
#'
#' @param x `data.frame` or compatible object
#' @param column `character` vector with one or more `colnames(x)`
#'    that should be de-concatenated.
#' @param split `character` pattern used by `strsplit()` to split
#'    multiple values in each column.
#' @param blank `character` string used to replace entries that
#'    would otherwise be zero-length as returned by `strsplit()`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' df <- data.frame(one=c("AB", "BC", "AC"),
#'    two=c("a,b", "b,c", "a,c"));
#' deconcat_df2(df, column="two")
#'
#' @export
deconcat_df2 <- function
(x,
   column,
   split="[,; |/]+",
   blank="",
   ...)
{
   ## Purpose is to expand multiple values by delimiter in a
   ## data.frame column by repeating rows for each value
   if (length(column) == 0) {
      stop("column must be supplied");
   }
   if (!all(column %in% colnames(x))) {
      stop("column must be present in colnames(x).");
   }
   if (length(column) > 1) {
      for (i in column) {
         x <- deconcat_df2(x,
            column=i,
            split=split,
            blank=blank,
            ...)
      }
      return(x)
   }
   if (!jamba::igrepHas(split, as.character(x[[column]]))) {
      return(x);
   }
   xv <- strsplit(as.character(x[[column]]),
      split=split);
   xvlen <- lengths(xv);
   if (any(xvlen == 0) && length(blank) >= 1) {
      xv[xvlen == 0] <- as.list(rep(blank,
         length.out=sum(xvlen == 0)));
   }
   xn <- rep(seq_len(nrow(x)), lengths(xv));
   xnew <- x[xn,,drop=FALSE];
   xnew[[column]] <- unname(unlist(xv));
   if (length(rownames(x)) > 0) {
      xrownames <- jamba::makeNames(rep(rownames(x), lengths(xv)));
      rownames(xnew) <- xrownames;
   }
   return(xnew);
}

