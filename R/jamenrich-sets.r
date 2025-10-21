# jamenrich-sets.r
# R functions handling sets and list

#' convert list to incidence matrix
#' 
#' convert list to incidence matrix
#' 
#' This function converts a list of vectors into an incidence matrix, where
#' the rows are the vector items and the columns are the list names.
#' 
#' @family jam list functions
#' 
#' @returns `numeric` matrix with values 0 or 1 indicating absense or
#'    presence of each item (row) in each set (column).
#'    The sets are derived from `names(x)`.
#' 
#' @param x `list` of item vectors, where vector values are item names
#' @param empty `numeric` used for missing or empty values, default 0.
#' @param do_sparse `logical` default FALSE, whether to return sparse Matrix.
#'    When `do_sparse=TRUE`, it sets `keepCounts=FALSE` since a sparse Matrix
#'    does not represent counts.
#' @param keepCounts `logical` default FALSE, whether to include the count
#'    of any duplicated items.
#' @param sort_rows `logical` default FALSE, whether to sort rows using
#'    `jamba::mixedOrder()`, or `function` can be supplied to use
#'    another sort algorithm.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' L1 <- list(A=c("C","A","B","A"),
#'    D=c("D","E","F","D"),
#'    A123=c(1:8,3,5),
#'    B=c("A1", "A2", "A10"),
#'    T=LETTERS[7:9]);
#' # Default behavior is to make items unique
#' list2im(L1);
#'
#' # Option to report the counts
#' list2im(L1, keepCounts=TRUE);
#' 
#' # Option to sort item rows with jamba::mixedOrder()
#' list2im(L1, sort_rows=TRUE);
#' 
#' # Option to sort item rows with vanilla sort()
#' list2im(L1, sort_rows=sort);
#' 
#' @export
list2im <- function
(x,
 empty=0,
 do_sparse=FALSE,
 keepCounts=FALSE,
 sort_rows=FALSE,
 ...)
{
   # optionally count items present multiple times
   if (TRUE %in% do_sparse) {
      keepCounts <- FALSE;
   }
   if (TRUE %in% keepCounts) {
      xCt <- jamba::rmNULL(lapply(x, jamba::tcount, minCount=2));
      if (length(xCt) == 0) {
         keepCounts <- FALSE;
      }
   }
   
   # generate unique items
   setnamesunion <- Reduce("union", x)
   if (length(empty) == 0) {
      empty <- NA
   } else {
      empty <- head(empty, 1)
   }
   setlistim <- do.call(cbind, lapply(x, function(i) {
      i_match <- match(i, setnamesunion)
      j <- rep(empty, length(setnamesunion))
      j[i_match] <- 1
      j
   }))
   rownames(setlistim) <- setnamesunion
   if (TRUE %in% do_sparse && requireNamespace("Matrix", quietly=TRUE)) {
      setlistim <- as(as(as(setlistim, "nMatrix"), "generalMatrix"), 
         "CsparseMatrix")
   }
   if (TRUE %in% keepCounts) {
      for (i in names(xCt)) {
         setlistim[names(xCt[[i]]), i] <- xCt[[i]];
      }
   }
   
   # optionally sort rows
   if (nrow(setlistim) > 1) {
      if (inherits(sort_rows, "function")) {
         use_rows <- sort_rows(rownames(setlistim))
         setlistim <- setlistim[match(use_rows, rownames(setlistim)), ];
      } else if (TRUE %in% sort_rows) {
         setlistim <- setlistim[jamba::mixedOrder(rownames(setlistim)), ];
      }
   }
   
   return(setlistim)
}


#' convert list to directional incidence matrix
#'
#' convert list to directional incidence matrix
#'
#' This function extends `list2im()` in that it stores the
#' value associated with each element in the list. As such, the input
#' format is a named vector, where the names of the vector are the items,
#' and the numeric values are the values to be stored in the
#' incidence matrix.
#'
#' A common scenario is to generate a vector of genes, with values
#' `c(-1, 0, 1)` indicating the direction of gene expression changes,
#' named by the gene symbol. Each vector in the list represents one
#' statistical test. Here, `list2imSigned()` will convert this list
#' into a directional matrix representing the gene changes across the
#' comparisons.
#'
#' Note that this function currently does not combine multiple values,
#' instead only the last occurring value is stored in the resulting
#' matrix. This decision is partly due to efficiency, and partly because
#' there are multiple possible methods to combine multiple values.
#' For example, taking the `mean(x)` for a given gene, which has
#' a value `1` and `-1` would result in `0` and might suggest the
#' gene is not a statistical hit. Instead, when multiple values
#' are anticipated per named vector entry, use functions in a
#' package like `data.table` or `dplyr` to apply a function to
#' combine values.
#'
#' @family jam list functions
#'
#' @returns `numeric` matrix with rownames defined by vector names from
#'    each vector in the input list. The colnames are defined by
#'    names of the input list if they exist.
#'    list. Values in the matrix are values from each vector.
#'
#' @param x `list` of named vectors, where the names are used to
#'    identify each element, and become rownames in the output
#'    incidence matrix. The vector values become values in the
#'    incidence matrix.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `list2im()`.
#'
#' @examples
#' L1 <- list(A=c("C","A","B","A"),
#'    D=c("D","E","F","D"),
#'    A123=c(1:8,3,5),
#'    T=LETTERS[7:9]);
#' L1;
#'
#' # Convert each vector to a signed vector
#' set.seed(123);
#' L2 <- lapply(L1, function(i){
#'    i <- unique(i);
#'    jamba::nameVector(sample(c(-1,1), size=length(i), replace=TRUE), i);
#' });
#' L2;
#'
#' # Convert to signed incidence matrix
#' list2imSigned(L2);
#'
#' @export
list2imSigned <- function
(x,
 emptyValue=NA,
 verbose=FALSE,
 ...)
{
   ## Purpose is to extend list2im() except maintain the directionality
   ## in the form of the sign for each entry.
   ##
   ## Input is expected to be a list of named numeric vectors, whose
   ## names are the entities to compare across sets.
   ##
   if (!jamba::igrepHas("list", class(x))) {
      if (jamba::igrepHas("array", class(x))) {
         x <- as.list(x);
      } else {
         stop("Input is expected to be a list class.");
      }
   }
   xClass <- sapply(x, function(i){
      jamba::igrepHas("numeric|integer|float|long|array", class(i)) &
         !is.null(names(i))
   });
   if (!all(TRUE %in% xClass)) {
      stop("Input is expected to be a list of named numeric vectors.");
   }
   #if (length(names(x)) == 0) {
   #   names(x) <- paste0("set", (seq_along(x)));
   #}
   ## For this step, only use unique elements, since we overwrite the value with the sign anyway
   imx <- list2im(lapply(x, names),
      emptyValue=NA,
      keepCounts=FALSE,
      ...);
   imx[] <- imx[] * 0;
   ## TODO: handle multiple values somehow... but not yet.
   ## It means we need to decide how to combine multiple signs,
   ## do we add them, average them, comma-delimit?
   for (i in seq_along(x)) {
      xi <- x[[i]];
      xi <- xi[!xi %in% c(NA)];
      imx[names(xi),i] <- xi;
   }
   #for (iName in names(x)) {
   #   imx[names(x[[iName]]),iName] <- x[[iName]];
   #}
   if (length(emptyValue) == 1 && !is.na(emptyValue) && any(is.na(imx))) {
      imx[is.na(imx)] <- emptyValue;
   }
   return(imx);
}

#' Convert list to concordance matrix
#'
#' Convert list to concordance matrix
#'
#' This function calculates pairwise concordance using
#' Kruskal concordance coefficient (ref) using the following equation:
#'
#' * (number_agree - number_disagree) / (total_shared)
#'
#' The equation is applied to each pair of named vectors in the input
#' list `x`, and reflects the degree of agreement in direction (+ or -)
#' between shared named elements, with +1 being perfect concordance
#' (agreement), and -1 being perfect discordance (disagreement.) Values
#' of zero indicate equal agreement and disagreement, and therefore
#' reflect no concordance nor discordance. Values
#' of `NA` occur when no named entries are shared.
#'
#' This function calls `list2imSigned()` to produce a signed
#' incidence matrix, which is then used with `base::crossprod()`
#' to calculate the full matrix of values.
#'
#' @family jam list functions
#'
#' @param x `list` of named numerical vectors, where the sign (positive
#'    or negative sign) indicates directionality, and is used to calculate
#'    concordance, which is a measure of the agreement of the overall set
#'    of directions shared between each pair of vectors.
#' @param naValue value passed to `jamba::rmNA()` used to replace any
#'    `NaN` values in the output matrix. The `NaN` values result when a
#'    pair of vectors has no shared non-zero named entry.
#' @param makeSigned logical indicating whether to force the vectors in
#'    the input list `x` to contain only values `c(-1,0,1)`, by calling
#'    `base::sign()`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `list2imSigned()`.
#'
#' @examples
#' set.seed(123);
#' l123 <- lapply(jamba::nameVector(1:3), function(i){
#'    jamba::nameVector(
#'       sample(c(-1,-1,0,1), replace=TRUE, size=15),
#'       letters[1:15]
#'    )
#' });
#' list2concordance(l123);
#'
#' # observe the signed incidence matrix
#' list2imSigned(l123);
#'
#' @export
list2concordance <- function
(x,
 naValue=NA,
 makeSigned=TRUE,
 verbose=FALSE,
 ...)
{
   # check if any values are non-sign
   if (makeSigned && !all(unique(unlist(x)) %in% c(-1,0,1))) {
      x <- lapply(x, sign);
   }

   # convert to signed incidence matrix
   imSigned <- list2imSigned(x,
      verbose=verbose,
      ...);

   # matrix math method
   # sum of the product is equal to (agreement - disagreement)
   # sum of the product of absolute values is equal to (number non-zero)
   concordM <- crossprod(imSigned) / crossprod(abs(imSigned));
   if (!identical(naValue, NaN)) {
      concordM <- jamba::rmNA(concordM,
         naValue=naValue);
   }
   return(concordM);
}

#' convert incidence matrix to list
#'
#' convert incidence matrix to list
#'
#' This function converts an incidence `matrix`, or equivalent
#' `data.frame`, to a list. The `matrix` should contain either
#' numeric values such as `c(0, 1)`, or logical values such
#' as `c(TRUE,FALSE)`, otherwise values are considered either
#' zero == `FALSE`, or non-zero == `TRUE`.
#'
#' The resulting list will be named by `colnames(x)` of the input,
#' and will contain members named by `rownames(x)` which are
#' either non-zero, or contain `TRUE`.
#'
#' Values of `NA` are converted to zero `0` and therefore ignored.
#'
#' @family jam list functions
#'
#' @param x `matrix` or equivalent object with `colnames(x)` indicating
#'    list set names, and `rownames(x)` indicating list contents.
#' @param empty `character` vector of incidence matrix values that
#'    should be considered "empty" and therefore do not indicate
#'    the row in `x` is present for the given column in `x`.
#'    All other items are considered to be present.
#' @param ... additional arguments are ignored.
#'
#' @return `list` of `character vectors`, where list names
#'    are defined by `colnames(x)`, and list elements are vectors
#'    that contain values from `rownames(x)`.
#'
#' @examples
#' im <- matrix(c(0,1,-1,1,1,NA,-1,0,1),
#'    ncol=3,
#'    nrow=3,
#'    dimnames=list(letters[1:3], LETTERS[1:3]))
#' print(im);
#' # matrix input
#' im2list(im);
#'
#' # data.frame
#' imdf <- data.frame(im);
#' print(imdf);
#' im2list(im);
#'
#' # logical input
#' imtf <- (!im == 0);
#' print(imtf);
#' im2list(imtf);
#'
#' @export
im2list <- function
(x,
 empty=c(NA, "", 0, FALSE),
 ...)
{
   # the reciprocal of list2im()
   x_rows <- rownames(x);
   x_cols <- colnames(x);

   # vicious bug when options("warn"=2) forcing warnings into errors
   # For now, force to max warn=1.
   # Who would do such a thing.
   if (getOption("warn", -1) > 1) {
      options("warn", 1)
   }

   l <- lapply(jamba::nameVector(x_cols), function(i){
      # i_empty <- as(empty, class(x[, i]));
      # has_value <- (!x[, i] %in% i_empty);
      has_value <- (!x[, i] %in% empty);
      x_rows[has_value];
   });
   return(l);
}


#' convert signed incidence matrix to list
#'
#' convert signed incidence matrix to list
#'
#' This function converts an signed incidence `matrix`
#' that contains positive and negative values, or equivalent
#' `data.frame`, to a list of named vectors containing values
#' `c(-1, 1)` to indicate signed direction.
#' The input `matrix` should contain numeric values where
#' positive and negative values indicate directionality.
#' When the input contains only logical values `c(TRUE,FALSE)`
#' the direction is assumed to be `+1` positive.
#'
#' Values of `NA` are converted to zero `0` and therefore ignored.
#'
#' Values that are `logical` with `TRUE` and `FALSE` are converted
#' to `numeric` before output.
#'
#' @family jam list functions
#'
#' @return `list` of named numeric vectors, where list names
#'    are defined by `colnames(x)`, and vector names are derived
#'    from `rownames(x)`. Values in each vector indicate the
#'    signed direction, `c(-1,1)`.
#'
#' @param x `matrix` or equivalent object with `colnames(x)` indicating
#'    list set names, and `rownames(x)` indicating list contents.
#' @param empty `character` vector of incidence matrix values that
#'    should be considered "empty" and therefore do not indicate
#'    the row in `x` is present for the given column in `x`.
#'    All other items are considered to be present, and are assigned
#'    direction based upon the value in that cell of `x`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' im <- matrix(c(0,1,-1,1,1,NA,-1,0,1),
#'    ncol=3,
#'    nrow=3,
#'    dimnames=list(letters[1:3], LETTERS[1:3]))
#' print(im);
#' # matrix input
#' im2list(im);
#' imSigned2list(im);
#' imSigned2list(im != 0);
#'
#' @export
imSigned2list <- function
(x,
 empty=c(NA, "", 0, FALSE),
 ...)
{
   # the reciprocal of list2im_value()
   x_rows <- rownames(x);
   x_cols <- colnames(x);
   l <- lapply(jamba::nameVector(x_cols), function(i){
      has_value <- (!x[,i] %in% empty);
      if (is.logical(x[,i])) {
         jamba::nameVector(
            as.numeric(x[has_value,i]),
            x_rows[has_value],
            makeNamesFunc=c);
      } else {
         jamba::nameVector(x[has_value,i],
            x_rows[has_value],
            makeNamesFunc=c);
      }
   });
   return(l);
}
