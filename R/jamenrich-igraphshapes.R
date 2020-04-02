
#' plot function for igraph vertex shape coloredrectangle
#'
#' plot function for igraph vertex shape coloredrectangle
#'
#' This function defines the plotting function for custom igraph vertex
#' shape coloredrectangle. The coloredrectangle shape is described as:
#'
#' \itemize{
#'    \item{a vertex drawn as a rectangle, filled with squares which are
#'   each individually colored.}
#'   \item{the squares are arrayed into a number of columns and rows,
#'   which are defined for each vertex using `coloredrect.ncol` and
#'   `coloredrect.nrow`, respectively.}
#'   \item{the vector of colors is arrayed as values of a matrix, therefore
#'   `coloredrect.byrow` is logical to indicate whether colors should fill
#'   the matrix by row, similar to how `byrow=TRUE` is used.}
#'   \item{the colors for each vertex are defined by
#'   `coloredrect.color`. When this value does not exist, this function
#'   will attempt to use values from `pie.color` if they exist.}
#'   \item{the size of the rectangle is defined by `size2`, which is
#'   typically the height of a `rectangle` node. The intent was to use
#'   a size parameter which is already convenient, but which does not
#'   conflict with the size used for vertex shapes such as `"circle"`.}
#'   \item{when a vertex has 3 columns, and 2 rows, and only 5 colors,
#'   the colors are not cycled to fill the complete node. Instead, the
#'   last color is `"transparent"` to render no color in that square.}
#'   \item{the frame around each node is described with `frame.color` and
#'   `frame.lwd` which draws one rectangular border around the full
#'   vertex. The border around all coloredrect vertices is drawn before the
#'   squares are drawn for all coloredrect vertices, so that the border will
#'   not overlap the squares. It has the visible effect of showing vertex
#'   borders behind the vertex colors, while allowing for vectorized drawing,
#'   which is substantially more efficient than per-vertex rendering.
#'   The `graphics::symbols()` function is used to render rectangles,
#'   it accepts only one value for `lwd` and `lty`, so rectangles
#'   are split by `lwd` and `lty` and rendered in order. As a
#'   result, if each vertex frame had a different lwd,lty combination, the
#'   rendering would effectively be non-vectorized. For large numbers of
#'   vertices and distinct lwd values, it would be preferred to reduce the
#'   lwd values to limited significant digits (e.g. with
#'   `signif(...,digits=2)`). That said, line width is not an optimal
#'   way to convey a quantitative measurement.}
#'   \item{a colored border can be drawn around each square in the vertex,
#'   defined by `coloredrect.border` and `coloredrect.lwd`,
#'   however it is recommended to use
#'   this color only as a visible break, for example the default is
#'   `"grey30"` which draws a simple border. Note that squares are
#'   rendered after the vertex frame color, so each square color border will
#'   be drawn over (on top of) the frame border.}
#'   \item{when `coloredrect.ncol` does not exist, it will attempt
#'   to use two rows of colors if `coloredrect.color` contains
#'   two or more values.}
#'   \item{when `coloredrect.nrow` does not exist, it will use a value
#'   based upon `coloredrect.ncol` and the length of
#'   `coloredrect.color`.}
#' }
#'
#' Currently, the plotting of coloredrectangle vertices does not define
#' clipping, which means edges are drawn to the center of each vertex,
#' and the coloredrectangle vertex shapes are drawn on top of the edges.
#' Transparent nodes will therefore show edges beneath them.
#'
#' @return invisible `list` of `data.frame` objects which were
#'    used to draw the rectangle objects.
#'    However the purpose of this function is the by-product that it
#'    draws rectangles onto an igraph graph.
#'
#' @seealso `igraph::shapes()`
#'
#' @family jam igraph shapes
#'
#' @param coords two-column numeric matrix with x- and y-coordinates,
#'    respectively.
#' @param v optional ids of the vertices to plot. It should match
#'    the number of rows in the `coords` argument.
#' @param params a function object that can be called to query
#'    vertex/edge/plot graphical parameters. The first argument of
#'    the function is `'vertex'`, `'edge'` or `'plot'` to decide
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

   verbose <- getOption("verbose");
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
   if (verbose) {
      jamba::printDebug("shape.coloredrectangle.plot(): ",
         "vertex.coloredrect.color:");
      print(vertex.coloredrect.color);
   }

   vertex.coloredrect.byrow <- getparam("coloredrect.byrow") > 0;
   if (length(vertex.coloredrect.byrow) == 0) {
      vertex.coloredrect.byrow <- TRUE;
   }

   vertex.coloredrect.ncol <- getparam("coloredrect.ncol");
   if (length(vertex.coloredrect.ncol) == 0) {
      vertex.coloredrect.ncol <- ceiling(length(vertex.coloredrect)/2);
   }
   vertex.coloredrect.nrow <- getparam("coloredrect.nrow");
   if (length(vertex.coloredrect.nrow) == 0) {
      vertex.coloredrect.nrow <- ceiling(
         length(vertex.coloredrect) /
            vertex.coloredrect.ncol);
   }
   vertex.coloredrect.lty <- getparam("coloredrect.lty");
   if (length(vertex.coloredrect.lty) == 0) {
      vertex.coloredrect.lty <- 1;
   }
   vertex.coloredrect.lwd <- getparam("coloredrect.lwd");
   if (length(vertex.coloredrect.lwd) == 0) {
      vertex.coloredrect.lwd <- 1;
   }

   vertex.size2 <- getparam("size");
   vertex.size1 <- getparam("size2");
   if (any(is.na(vertex.size1))) {
      vertex.size1[is.na(vertex.size1)] <- vertex.size2[is.na(vertex.size1)] / 2;
   }
   vertex.size1 <- rep(1/200 * vertex.size1, length.out=nrow(coords));
   vertex.size2 <- rep(1/200 * vertex.size2, length.out=nrow(coords));

   ## Use size1 (height) to define the size of each square, then apply that
   ## to calculate size2 (width)
   #vertex.size1 <- vertex.size1 * 5 / vertex.coloredrect.nrow;
   vertex.size1 <- vertex.size1 * 5 * pmax(vertex.coloredrect.nrow,
      vertex.coloredrect.ncol) / vertex.coloredrect.nrow;
   vertex.size2 <- (vertex.size1 *
         (vertex.coloredrect.nrow) / vertex.coloredrect.ncol);

   if (verbose) {
      jamba::printDebug("shape.coloredrectangle.plot(): ",
         "vertex.size1:", vertex.size1,
         ", vertex.coloredrect.nrow:", vertex.coloredrect.nrow);
      jamba::printDebug("shape.coloredrectangle.plot(): ",
         "vertex.size2:", vertex.size2,
         ", vertex.coloredrect.ncol:", vertex.coloredrect.ncol);
   }

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
      frame.lwd=0.5,
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
      rectDF <- jamba::rbindList(lapply(seq_along(x), function(k){
         xk <- x[[k]];
         yk <- y[[k]];
         size1k <- size1[[k]];
         size2k <- size2[[k]];
         nrowk <- jamba::rmNULL(nullValue=1, nrow[[k]]);
         ncolk <- jamba::rmNULL(nullValue=1, ncol[[k]]);
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
            jamba::printDebug("k:", k,
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
         jamba::pasteByRowOrdered(rectDF[,c("rect_type","lwd","lty")]));
      if (verbose) {
         jamba::printDebug("shape.coloredrectangle.plot(): ",
            "names(rectDFL):", names(rectDFL));
      }

      rect_order <- jamba::provigrep(c("frame", "square", "."), names(rectDFL));
      for (rectDFi in rectDFL[rect_order]) {
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

#' plot function for igraph vertex shape ellipse
#'
#' plot function for igraph vertex shape ellipse
#'
#' This function defines the plotting function for custom igraph vertex
#' shape ellipse.
#'
#' @family jam igraph shapes
#'
#' @export
shape.ellipse.plot <- function
(coords,
   v=NULL,
   params)
{
   vertex.color <- params("vertex", "color");
   if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v];
   }
   vertex.frame.color <- params("vertex", "frame.color");
   if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v];
   }
   vertex.frame.width <- params("vertex", "frame.width");
   if (length(vertex.frame.width) == 0) {
      vertex.frame.width <- 1;
   }
   if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v];
   }
   vertex.size <- 1/200 * params("vertex", "size");
   if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v];
   }
   vertex.ellipse.ratio <- params("vertex", "ellipse.ratio");
   if (length(vertex.ellipse.ratio) == 0) {
      vertex.ellipse.ratio <- 2;
   }
   if (length(vertex.ellipse.ratio) != 1 && !is.null(v)) {
      vertex.ellipse.ratio <- vertex.ellipse.ratio[v];
   }

   drawEllipse(x=coords[,1],
      y=coords[,2],
      a=vertex.size,
      b=vertex.size/vertex.ellipse.ratio,
      col=vertex.color,
      border=vertex.frame.color,
      lwd=vertex.frame.width,
      draw=TRUE);
}

#' plot function for igraph vertex shape jampie
#'
#' plot function for igraph vertex shape jampie
#'
#' This function is a vectorized replacement for plotting
#' vertex shape `"pie"` in much more efficient manner.
#'
#' @family jam igraph shapes
#'
#' @export
shape.jampie.plot <- function
(coords,
 v = NULL,
 params)
{
   getparam <- function(pname) {
      p <- params("vertex", pname)
      if (length(p) != 1 && !is.null(v)) {
         p <- p[v]
      }
      p
   }
   vertex.color <- getparam("color")
   vertex.frame.color <- getparam("frame.color")
   vertex.size <- rep(1/200 * getparam("size"), length = nrow(coords))
   vertex.pie <- getparam("pie")
   vertex.pie.color <- getparam("pie.color")
   vertex.pie.angle <- getparam("pie.angle")
   vertex.pie.density <- getparam("pie.density")
   vertex.pie.lty <- getparam("pie.lty")
   ## Convert for loop to lapply that obtains polygon coordinates
   if (1 == 1) {
      #jamba::printDebug("Calculating pie node polygons using ",
      #   "jam_mypie()");
      poly_df <- jamba::rbindList(lapply(seq_len(nrow(coords)), function(i){
         pie <- if (length(vertex.pie) == 1) {
            vertex.pie[[1]]
         } else {
            vertex.pie[[i]]
         }
         col <- if (length(vertex.pie.color) == 1) {
            vertex.pie.color[[1]]
         } else {
            vertex.pie.color[[i]]
         }
         jam_mypie(x = coords[i, 1],
            y = coords[i, 2],
            pie,
            radius = vertex.size[i],
            edges = 200,
            col = col,
            angle = na.omit(vertex.pie.angle[c(i, 1)])[1],
            density = na.omit(vertex.pie.density[c(i, 1)])[1],
            border = na.omit(vertex.frame.color[c(i, 1)])[1],
            lty = na.omit(vertex.pie.lty[c(i, 1)])[1])
      }));
      poly_xv <- jamba::cPaste(as(poly_df$x, "list"));
      poly_yv <- jamba::cPaste(as(poly_df$y, "list"));
      poly_x <- as.numeric(unlist(strsplit(paste(poly_xv, collapse=",NA"), ",")));
      poly_y <- as.numeric(unlist(strsplit(paste(poly_yv, collapse=",NA"), ",")));
      #jamba::printDebug("Plotting pie node polygons.");
      polygon(x=poly_x,
         y=poly_y,
         density = poly_df$density,
         angle = poly_df$angle,
         border = poly_df$border,
         col = poly_df$col,
         lty = poly_df$lty)
   } else {
      ## Legacy for loop below
      for (i in seq_len(nrow(coords))) {
         pie <- if (length(vertex.pie) == 1) {
            vertex.pie[[1]]
         } else {
            vertex.pie[[i]]
         }
         col <- if (length(vertex.pie.color) == 1) {
            vertex.pie.color[[1]]
         } else {
            vertex.pie.color[[i]]
         }
         igraph:::mypie(x = coords[i, 1],
            y = coords[i, 2],
            pie,
            radius = vertex.size[i],
            edges = 200,
            col = col,
            angle = na.omit(vertex.pie.angle[c(i, 1)])[1],
            density = na.omit(vertex.pie.density[c(i, 1)])[1],
            border = na.omit(vertex.frame.color[c(i, 1)])[1],
            lty = na.omit(vertex.pie.lty[c(i, 1)])[1])
      }
   }
}

#' Vectorized mypie() function for igraph vertex pie polygons
#'
#' Vectorized mypie() function for igraph vertex pie polygons
#'
#' This function is a light rewrite of `igraph:::mypie()`, except
#' that this function determines polygon coordinates without
#' drawing them, instead returns the polygon coordinates to
#' the calling function `shape.jampie.plot()` which in turn
#' draws all polygons once using the vectorized approach
#' described for `graphics::polygon()`.
#'
#' One small additional change, pie shapes with only one large
#' 100% wedge no longer display the small line from origin.
#'
#' @family jam igraph shapes
#'
#' @export
jam_mypie <- function
(x,
 y,
 values,
 radius,
 edges = 200,
 col = NULL,
 angle = 45,
 density = NULL,
 border = NULL,
 lty = NULL,
 init.angle = 90,
 ...)
{
   values <- c(0, cumsum(values)/sum(values))
   dx <- diff(values)
   nx <- length(dx)
   twopi <- 2 * pi
   if (is.null(col)) {
      if (is.null(density)) {
         col <- c("white",
            "lightblue",
            "mistyrose",
            "lightcyan",
            "lavender",
            "cornsilk")
      } else {
         col <- par("fg")
      }
   }
   col <- rep(col, length.out = nx)
   border <- rep(border, length.out = nx)
   lty <- rep(lty, length.out = nx)
   angle <- rep(angle, length.out = nx)
   density <- rep(density, length.out = nx)
   t2xy <- function(t) {
      t2p <- 2 * pi * t + init.angle * pi/180;
      list(x = radius * cos(t2p),
         y = radius * sin(t2p))
   }
   ## convert this loop to return matrix of polygon coordinates
   ## then stitch multiple polygons with empty row between each
   ## so polygon() can be called once on the whole set of polygons
   ##
   ## bonus points: for one-section pie graph, draw a circle without segment
   poly_df <- jamba::rbindList(lapply(seq_len(nx), function(i){
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(seq.int(values[i],
         values[i + 1],
         length.out = n))
      if (nx == 1) {
         xvals <- (x + c(P$x));
         yvals <- (y + c(P$y));
      } else {
         xvals <- (x + c(P$x, 0));
         yvals <- (y + c(P$y, 0));
      }
      data.frame(x=I(list(xvals)),
         y=I(list(yvals)),
         density = density[i],
         angle = angle[i],
         border = border[i],
         col = col[i],
         lty = lty[i])
   }));
   return(poly_df);
   ## Legacy for loop below
   for (i in 1:nx) {
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(seq.int(values[i], values[i + 1], length.out = n))
      polygon(x + c(P$x, 0),
         y + c(P$y, 0),
         density = density[i],
         angle = angle[i],
         border = border[i],
         col = col[i],
         lty = lty[i],
         ...)
   }
}
