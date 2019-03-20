# jamenrich-sets.r
# R functions handling sets and list

#' convert list to incidence matrix
#'
#' convert list to incidence matrix
#'
#' This function converts a list of vectors into an incidence matrix, where
#' the rows are the vector items and the columns are the list names.
#' It uses an object from the \code{\line{arules}} package called
#' \code{\link[arules]{transactions}} which offers highly efficient methods
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
#' @param x list of vectors
#' @param makeUnique boolean indicating whether to enforce uniqueness on each
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
   makeListUnique <- function(x, ...) {
      if (suppressPackageStartupMessages(require(S4Vectors))) {
         x <- as.list(unique(SimpleList(x)));
      } else {
         x <- lapply(x, unique);
      }
      return(x);
   }
   if (makeUnique) {
      x <- makeListUnique(x);
      keepCounts <- FALSE;
   } else {
      if (keepCounts) {
         xCt <- rmNULL(lapply(x, tcount, minCount=2));
         if (length(xCt) == 0) {
            keepCounts <- FALSE;
         }
         x <- makeListUnique(x);
      }
   }

   ## Convert to transactions
   xT <- as(x, "transactions");

   ## Extract the matrix
   xM <- as.matrix(xT@data*1);
   colnames(xM) <- xT@itemsetInfo[,1];
   rownames(xM) <- xT@itemInfo[,1];
   if (keepCounts) {
      if (verbose) {
         printDebug("list2im(): ",
            "Putting item counts into the incidence matrix for ",
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
#' This function extends \code{\link{list2im}} in that it stores the
#' value associated with each element in the list. As such, the input
#' format is a named vector, where the names of the vector are the items,
#' and the numeric values are the values to be stored in the
#' incidence matrix.
#'
#' A common scenario is to generate a vector of genes, with values
#' \code{c(-1, 0, 1)} indicating the direction of gene expression changes,
#' named by the gene symbol. Each vector in the list represents one
#' statistical test. Here, \code{list2imSigned()} will convert this list
#' into a directional matrix representing the gene changes across the
#' comparisons.
#'
#' @return numeric matrix with rownames defined by names from each vector
#'    in the input list, and colnames defined by the names of the input
#'    list. Values in the matrix are values from each vector.
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
 makeUnique=TRUE,
 doTimer=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to extend list2im() except maintain the directionality
   ## in the form of the sign for each entry.
   ##
   ## Input is expected to be a list of named numeric vectors, whose
   ## names are the entities to compare across sets.
   ##
   if (!igrepHas("list", class(x))) {
      if (igrepHas("array", class(x))) {
         x <- as.list(x);
      } else {
         stop("Input is expected to be a list class.");
      }
   }
   xClass <- sapply(x, function(i){
      igrepHas("numeric|integer|float|long|array", class(i)) &
         !is.null(names(i))
   });
   if (!all(xClass)) {
      stop("Input is expected to be a list of named numeric vectors.");
   }
   if (length(names(x)) == 0) {
      names(x) <- paste0("set", (seq_along(x)));
   }
   ## For this step, only use unique elements, since we overwrite the value with the sign anyway
   imx <- list2im(lapply(x, names),
      makeUnique=TRUE,
      verbose=verbose,
      ...);
   ## TODO: handle multiple values somehow... which means we need to decide
   ## how to combine multiple signs, do we add them, average them, what?
   for (iName in names(x)) {
      imx[names(x[[iName]]),iName] <- x[[iName]];
   }
   return(imx);
}
