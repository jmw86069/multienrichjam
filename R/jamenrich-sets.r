# jamenrich-sets.r
# R functions handling sets and list

#' convert list to incidence matrix
#'
#' convert list to incidence matrix
#'
#' This function converts a list of vectors into an incidence matrix, where
#' the rows are the vector items and the columns are the list names.
#' It uses an object from the `arules` package called
#' `arules::transactions` which offers highly efficient methods
#' for interconverting from list to matrix. The
#' \code{\link[arules]{transactions}} class is itself an enhanced data matrix,
#' which stores data using sparse matrix object type from the
#' \code{\link{Matrix}} package, but also associates a data.frame to both
#' the rows and columns of the matrix to offer additional row and column
#' annotation, as needed.
#'
#' Performance benchmarks showed high speed of converting a list to a matrix,
#' but also that the resulting matrix was substantially smaller (5-20 times)
#' then comparable methods producing a data matrix.
#'
#' @family jam list functions
#'
#' @param x list of vectors
#' @param makeUnique logical indicating whether to enforce uniqueness on each
#'    vector in the list. For an extremely long list, the uniqueness step can
#'    become the rate-limiting step, and yet the function used to coerce list
#'    to matrix requires each vector have a unique set of items. When
#'    \code{makeUnique=FALSE]} the value returned depends upon
#'    \code{keepCounts} which optionally reports the frequency of each
#'    item.
#' @param keepCounts boolean indicating whether to return values indicating
#'    the number of occurrences of each item, only useful when
#'    \code{makeUnique=FALSE}.
#' @param verbose boolean indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @return numeric matrix whose rownames were vector items of the input list,
#'    and whole colnames were list names.
#'
#' @examples
#' L1 <- list(A=c("C","A","B","A"),
#'    D=c("D","E","F","D"),
#'    A123=c(1:8,3,5),
#'    T=LETTERS[7:9]);
#' # Default behavior is to make items unique
#' list2im(L1);
#'
#' # When not using unique entries, by default it reports the counts
#' list2im(L1, makeUnique=FALSE);
#'
#' # The highest efficiency is seen when makeUnique=FALSE, keepCounts=FALSE,
#' # however the difference is seen typically only with lists longer
#' # than 10,000.
#' list2im(L1, makeUnique=FALSE, keepCounts=FALSE);
#'
#' @export
list2im <- function
(x,
 makeUnique=TRUE,
 keepCounts=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert a list of vectors into an incident matrix
   ## using the arules package
   if (!suppressPackageStartupMessages(require(arules))) {
      stop("list2im() requires the arules package.");
   }
   makeListUnique <- function
   (x,
    verbose=FALSE,
    ...) {
      if (suppressPackageStartupMessages(require(S4Vectors))) {
         if (verbose) {
            jamba::printDebug("makeListUnique():",
               "Calling as.list(unique(S4Vectors::SimpleList(x)))");
         }
         x <- as.list(unique(SimpleList(x)));
      } else {
         if (verbose) {
            jamba::printDebug("makeListUnique():",
               "Calling lapply(x, unique)");
         }
         x <- lapply(x, unique);
      }
      return(x);
   }
   if (makeUnique) {
      if (keepCounts) {
         xCt <- rmNULL(lapply(x, jamba::tcount, minCount=2));
         if (length(xCt) == 0) {
            if (verbose) {
               jamba::printDebug("list2im():",
                  "No duplicate values observed.");
            }
            keepCounts <- FALSE;
         }
      } else {
         x <- makeListUnique(x,
            verbose=verbose);
      }
   } else {
      keepCounts <- FALSE;
   }

   ## Convert to transactions
   xT <- as(x, "transactions");

   ## Extract the matrix
   xM <- as.matrix(xT@data*1);
   if (ncol(xT@itemsetInfo) > 0) {
      colnames(xM) <- xT@itemsetInfo[,1];
   }
   if (ncol(xT@itemInfo) > 0) {
      rownames(xM) <- xT@itemInfo[,1];
   }
   if (keepCounts) {
      if (verbose) {
         printDebug("list2im(): ",
            "Applying item counts to the incidence matrix ",
            format(big.mark=",", length(xCt)),
            " items.");
      }
      for (i in seq_along(xCt)) {
         xM[names(xCt[[i]]),i] <- xCt[[i]];
      }
   }
   return(xM);
}

#' convert list to signeddirectional incidence matrix
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
#' @return numeric matrix with rownames defined by vector names from
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
#'    nameVector(sample(c(-1,1), size=length(i), replace=TRUE), i);
#' });
#' L2;
#'
#' # Convert to signed incidence matrix
#' list2imSigned(L2);
#'
#' @export
list2imSigned <- function
(x,
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
   if (!all(xClass)) {
      stop("Input is expected to be a list of named numeric vectors.");
   }
   #if (length(names(x)) == 0) {
   #   names(x) <- paste0("set", (seq_along(x)));
   #}
   ## For this step, only use unique elements, since we overwrite the value with the sign anyway
   imx <- list2im(lapply(x, names),
      makeUnique=TRUE,
      keepCounts=FALSE,
      verbose=verbose,
      ...);
   imx[] <- 0;
   ## TODO: handle multiple values somehow... but not yet.
   ## It means we need to decide how to combine multiple signs,
   ## do we add them, average them, comma-delimit?
   for (i in seq_along(x)) {
      xi <- x[[i]];
      xi <- xi[xi != 0];
      imx[names(xi),i] <- xi;
   }
   #for (iName in names(x)) {
   #   imx[names(x[[iName]]),iName] <- x[[iName]];
   #}
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
   # check method
   method <- match.arg(method);
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
