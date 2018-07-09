# jamenrich-base.r

.onLoad <- function
(libname,
 pkgname)
{
   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("coloredrectangle",
      clip=igraph::shape_noclip,
      plot=shape.coloredrectangle.plot)

}

#' plot function for igraph vertex shape coloredrectangle
#'
#' plot function for igraph vertex shape coloredrectangle
#'
#' This function defines the plotting function for custom igraph vertex
#' shape coloredrectangle. The coloredrectangle shape is described as:
#'
#' \itemize{
#'    \item a vertex drawn as a rectangle, filled with squares which are
#'   each individually colored.
#'   \item the squares are arrayed into a number of columns and rows,
#'   which are defined for each vertex using \code{coloredrect.ncol} and
#'   \code{coloredrect.nrow}, respectively.
#'   \item the vector of colors is arrayed as values of a matrix, therefore
#'   \code{coloredrect.byrow} is boolean to indicate whether colors should fill
#'   the matrix by row, similar to how \code{byrow=TRUE} is used.
#'   \item the colors for each vertex are defined by
#'   \code{coloredrect.color}. When this value does not exist, this function
#'   will attempt to use values from \code{pie.color} if they exist.
#'   \item the size of the rectangle is defined by \code{size2}, which is
#'   typically the height of a \code{rectangle} node. The intent was to use
#'   a size parameter which is already convenient, but which does not
#'   conflict with the size used for vertex shapes such as \code{'circle'}.
#'   \item when a vertex has 3 columns, and 2 rows, and only 5 colors,
#'   the colors are not cycled to fill the complete node. Instead, the
#'   last color is \code{"transparent"} to render no color in that square.
#'   \item the frame around each node is described with \code{frame.color} and
#'   \code{frame.lwd} which draws one rectangular border around the full
#'   vertex. The border around all coloredrect vertices is drawn before the
#'   squares are drawn for all coloredrect vertices, so that the border will
#'   not overlap the squares. It has the visible effect of showing vertex
#'   borders behind the vertex colors, while allowing for vectorized drawing,
#'   which is substantially more efficient than per-vertex rendering.
#'   The \code{\link{symbols}} function is used to render rectangles,
#'   it accepts only one value for \code{lwd} and \code{lty}, so rectangles
#'   are split by \code{lwd} and \code{lty} and rendered in order. As a
#'   result, if each vertex frame had a different lwd,lty combination, the
#'   rendering would effectively be non-vectorized. For large numbers of
#'   vertices and distinct lwd values, it would be preferred to reduce the
#'   lwd values to limited significant digits (e.g. with
#'   \code{signif(...,digits=2)}). That said, line width is not an optimal
#'   way to convey a quantitative measurement.
#'   \item a colored border can be drawn around each square in the vertex,
#'   defined by \code{coloredrect.border} and \code{coloredrect.lwd},
#'   however it is recommended to use
#'   this color only as a visible break, for example the default is
#'   \code{"grey30"} which draws a simple border. Note that squares are
#'   rendered after the vertex frame color, so each square color border will
#'   be drawn over (on top of) the frame border.
#'   \item when \code{coloredrect.ncol} does not exist, it will attempt
#'   to use two rows of colors if \code{coloredrect.color} contains
#'   two or more values.
#'   \item when \code{coloredrect.nrow} does not exist, it will use a value
#'   based upon \code{coloredrect.ncol} and the length of
#'   \code{coloredrect.color}.
#' }
#'
#' Currently, the plotting of coloredrectangle vertices does not define
#' clipping, which means edges are drawn to the center of each vertex,
#' and the coloredrectangle vertex shapes are drawn on top of the edges.
#' Transparent nodes will therefore show edges beneath them.
#'
#' @return invisible \code{list} of \code{data.frame} objects which were
#'    used to draw the rectangle objects.
#'    However the purpose of this function is the by-product that it
#'    draws rectangles onto an igraph graph.
#'
#' @seealso \code{\link[igraph]{shapes}}
#'
#' @param coords two-column numeric matrix with x- and y-coordinates,
#'    respectively.
#' @param v optional ids of the vertices to plot. It should match
#'    the number of rows in the \code{coords} argument.
#' @param params a function object that can be called to query
#'    vertex/edge/plot graphical parameters. The first argument of
#'    the function is \code{'vertex'}, \code{'edge'} or \code{'plot'} to decide
#'    the type of the parameter, the second is a character string
#'    giving the name of the parameter.
#'
#' @export
shape.coloredrectangle.plot <- function
(coords,
 v=NULL,
 params)
{
   ## Purpose is to extend igraph:::.igraph.shape.rectangle.plot
   ## to handle the custom colored rectangle node type.
   ## Shape "coloredrectangle" has these features for each vertex:
   ## * the vertex is a rectangle filled with squares, where each
   ##   square has its own assigned color.
   ## * the squares are arrayed into a number of columns defined
   ##   by "coloredrect.ncol", and a number of rows defined by
   ##   "coloredrect.nrow". Each vertex can have its own number
   ##   of columns and rows, and its own number of colors.
   ## * the size of the overall rectangle is defined by "size2",
   ##   which is normally the height of a rectangle shape. The
   ##   intention was to allow using "size" for shapes such as
   ##   "circle" but "size2" for
   ##
   ## This function combines some logic from
   ## igraph:::.igraph.shape.pie.plot() as well.

   getparam <- function(pname) {
      p <- params("vertex", pname)
      if (length(p) != 1 && !is.null(v)) {
         p <- p[v]
      }
      p;
   }
   vertex.color <- getparam("color");

   vertex.frame.color <- getparam("frame.color");
   if (length(vertex.frame.color) == 0) {
      vertex.frame.color <- "grey30";
   }
   vertex.frame.lwd <- getparam("frame.lwd");
   if (length(vertex.frame.lwd) == 0) {
      vertex.frame.lwd <- 1;
   }

   vertex.coloredrect.color <- getparam("coloredrect.color");
   if (length(vertex.coloredrect.color) == 0) {
      vertex.coloredrect.color <- getparam("pie.color");
   }

   vertex.coloredrect.border <- getparam("coloredrect.border");
   if (length(vertex.coloredrect.border) == 0) {
      vertex.coloredrect.border <- "grey30";
      #vertex.coloredrect.border <- vertex.frame.color;
   }

   vertex.coloredrect.byrow <- getparam("coloredrect.byrow");
   if (length(vertex.coloredrect.byrow) == 0) {
      vertex.coloredrect.byrow <- TRUE;
   }

   vertex.coloredrect.ncol <- getparam("coloredrect.ncol");
   if (length(vertex.coloredrect.ncol) == 0) {
      vertex.coloredrect.ncol <- ceiling(length(vertex.coloredrect)/2);
   }
   vertex.coloredrect.nrow <- getparam("coloredrect.nrow");
   if (length(vertex.coloredrect.nrow) == 0) {
      vertex.coloredrect.nrow <- ceiling(length(vertex.coloredrect)/vertex.coloredrect.ncol);
   }
   vertex.coloredrect.lty <- getparam("coloredrect.lty");
   if (length(vertex.coloredrect.lty) == 0) {
      vertex.coloredrect.lty <- 1;
   }
   vertex.coloredrect.lwd <- getparam("coloredrect.lwd");
   if (length(vertex.coloredrect.lwd) == 0) {
      vertex.coloredrect.lwd <- 1;
   }

   vertex.size1 <- getparam("size2");
   vertex.size2 <- getparam("size");
   #printDebug("vertex.size1:", head(vertex.size1, 10),
   #   ", vertex.size2:", head(vertex.size2, 10));
   vertex.size1 <- rep(1/200 * vertex.size1, length.out=nrow(coords));
   vertex.size2 <- rep(1/200 * vertex.size2, length.out=nrow(coords));

   ## Use size1 (height) to define the size of each square, then apply that
   ## to calculate size2 (width)
   vertex.size1 <- vertex.size1 * 5 / vertex.coloredrect.nrow;
   vertex.size2 <- (vertex.size1 *
         (vertex.coloredrect.nrow / vertex.coloredrect.ncol));

   #printDebug("vertex.size1:", head(vertex.size1, 10),
   #   ", vertex.size2:", head(vertex.size2, 10));
   vertex.size <- cbind(vertex.size1, vertex.size2);

   ## Define custom function to help vectorize drawing, by creating a
   ## data.frame of coordinates for each square and rectangle
   mycoloredrectangle <- function
   (x,
    y,
    size1,
    size2,
    col,
    ncol,
    nrow,
    border,
    lty,
    lwd=1,
    frame.lwd=1,
    frame.color="grey30",
    byrow=TRUE,
    ...)
   {
      ## Purpose is to draw a rectangle filles with multi-color squares
      nrow <- rep(nrow, length.out=length(x));
      ncol <- rep(ncol, length.out=length(x));
      border <- rep(border, length.out=length(x));
      lty <- rep(lty, length.out=length(x));
      lwd <- rep(lwd, length.out=length(x));
      frame.lwd <- rep(frame.lwd, length.out=length(x));
      frame.color <- rep(frame.color, length.out=length(x));
      byrow <- rep(byrow, length.out=length(x));

      ## Iterate each vertex, create a data.frame describing
      ## frame and square colors, then combine into one large
      ## data.frame for vectorized drawing.
      rectDF <- rbindList(lapply(seq_along(x), function(k){
         xk <- x[[k]];
         yk <- y[[k]];
         size1k <- size1[[k]];
         size2k <- size2[[k]];
         nrowk <- nrow[[k]];
         ncolk <- ncol[[k]];
         byrowk <- byrow[[k]];
         ## Individual rectangles
         x01 <- seq(from=xk-(size1k/2),
            to=xk+(size1k/2),
            length.out=ncolk+1);
         y01 <- seq(from=yk-(size2k/2),
            to=yk+(size2k/2),
            length.out=nrowk+1);

         numk <- ncolk * nrowk;
         size1v <- rep(min(diff(x01)), numk);
         size2v <- rep(min(diff(y01)), numk);
         if (byrowk) {
            x01v <- rep(head(x01, -1), nrowk)+(size1v/2);
            y01v <- rep(rev(tail(y01, -1)), each=ncolk)-(size2v/2);
         } else {
            x01v <- rep(head(x01, -1), each=nrowk)+(size1v/2);
            y01v <- rep(rev(tail(y01, -1)), ncolk)-(size2v/2);
         }
         ## Make one large rectangle border
         x01all <- mean(x01v);
         y01all <- mean(y01v);
         ##
         colk <- rep(col[[k]], length.out=length(x01v));
         borderk <- rep(border[[k]], length.out=length(colk));
         ltyk <- rep(lty[[k]], length.out=length(colk));
         lwdk <- rep(lwd[[k]], length.out=length(colk));
         if (length(getOption("debug"))>0) {
            printDebug("k:", k,
               ", numk:", numk,
               ", byrowk:", byrowk,
               ", ncolk:", ncolk,
               ", nrowk:", nrowk,
               ", xk:", xk,
               ", x01v:", signif(digits=3, x01v),
               ", x01:", signif(digits=3, x01),
               ", yk:", yk,
               ", y01v:", signif(digits=3, y01v),
               ", y01:", signif(digits=3, y01),
               ", colk:", colk,
               ", size1v:", signif(digits=3, size1v),
               ", size2v:", signif(digits=3, size2v)
            );
         }
         kDF <- data.frame(x=c(xk, x01v),
            y=c(yk, y01v),
            bg=c("transparent", colk),
            fg=c(head(frame.color[[k]],1),
               borderk),
            rectx=c(head(size1v,1)*ncolk, size1v),
            recty=c(head(size2v,1)*nrowk, size2v),
            lty=c(head(ltyk,1), ltyk),
            lwd=c(head(frame.lwd[[k]],1)*2.5,
               lwdk/4),
            rect_type=rep(factor(c("frame","square")),
               c(1,length(borderk)))
         );
         kDF;
      }));

      ## Split into a list of data.frames, because symbols()
      ## can only use one value for lwd and lty.
      rectDFL <- split(rectDF,
         pasteByRowOrdered(rectDF[,c("rect_type","lwd","lty")]));

      for (rectDFi in rev(rectDFL)) {
         symbols(x=rectDFi$x,
            y=rectDFi$y,
            bg=rectDFi$bg,
            fg=rectDFi$fg,
            rectangles=as.matrix(rectDFi[,c("rectx","recty")]),
            add=TRUE,
            inches=FALSE,
            lty=rectDFi$lty,
            lwd=rectDFi$lwd);
      }
      return(rectDFL);
   }

   ## For reference, the code below draws a single rectangle
   if (1 == 2) {
      symbols(x=coords[, 1],
         y=coords[, 2],
         bg=vertex.color,
         fg=vertex.frame.color,
         rectangles=2*vertex.size,
         add=TRUE,
         inches=FALSE);
   }

   mcr <- mycoloredrectangle(x=coords[,1],
      y=coords[,2],
      size1=vertex.size1,
      size2=vertex.size2,
      col=vertex.coloredrect.color,
      ncol=vertex.coloredrect.ncol,
      nrow=vertex.coloredrect.nrow,
      border=vertex.coloredrect.border,
      lty=vertex.coloredrect.lty,
      lwd=vertex.coloredrect.lwd,
      frame.lwd=vertex.frame.lwd,
      frame.color=vertex.frame.color,
      byrow=vertex.coloredrect.byrow);

}
