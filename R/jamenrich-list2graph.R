#' list2df conversion
#' @noRd
list2df_ggt <- function(inputList) {
   ldf <- lapply(seq_len(length(inputList)), function(i) {
      data.frame(
         categoryID=rep(names(inputList[i]),
         length(inputList[[i]])),
         Gene=inputList[[i]])
   })
   
   do.call('rbind', ldf)
}

#' list2graph conversion
#' @noRd
list2graph_ggt <- function
(inputList,
 directed=FALSE,
 add_size=FALSE,
 ...)
{
   x <- list2df(inputList)
   g <- igraph::graph_from_data_frame(x,
      directed=directed)
   
   # igraph::V(g)$.isCategory <- igraph::V(g)$name %in% names(inputList)
   
   if (TRUE %in% add_size) {
      size <- vapply(inputList, length, FUN.VALUE=numeric(1))
      igraph::V(g)$size <- igraph::degree(g)
   }
   return(g)
}
