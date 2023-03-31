
#' custom igraph vertex shape coloredrectangle
#'
#' custom igraph vertex shape coloredrectangle
#'
#' This function defines the plotting function for custom igraph vertex
#' shape coloredrectangle. The coloredrectangle shape is described as:
#'
#'   * Each vertex is drawn as a rectangle, filled with squares which are
#'   each individually colored.
#'   * The squares are arrayed into a number of columns and rows,
#'   which are defined for each vertex using `coloredrect.ncol` and
#'   `coloredrect.nrow`, respectively.
#'   * The vector of colors is arrayed as values of a matrix, therefore
#'   `coloredrect.byrow` is logical to indicate whether colors should fill
#'   the matrix by row, similar to how `byrow=TRUE` is used.
#'   * The colors for each vertex are defined by `coloredrect.color`,
#'   which is a `list` of `character` color vectors.
#'   When `coloredrect.color` does not exist, values from `pie.color`
#'   will be used if they exist.
#'   * Any missing colors are displayed as `NA` which applies no color.
#'   For example, when a vertex has 3 columns, 2 rows, and only 5 colors,
#'   the colors are not recycled. Instead, the last color is `NA` to
#'   render no color in that position.
#'   * Each square may also have an optional border, defined by
#'   `coloredrect.border`, `coloredrect.lwd`, and `coloredrect.lty`.
#'   These borders are drawn as inner borders so they do not overlap
#'   optional frame border.
#'   * Each node may have a frame border, defined by `frame.color`,
#'   `frame.lwd`, and `frame.lty`. The frame border is drawn as an outer
#'   border so it does not overlap inner borders, adding some height and width
#'   to the final node. The frame border is drawn before the
#'   node squares are drawn.
#'   * The size of each square inside the rectangle is defined by `size2`,
#'   such that the rectangle width is `coloredrect.ncol * size2` and
#'   the rectangle height is `coloredrect.nrow * size2`.
#'   * Node sizes may be adjusted by enabling `equalize_sizes`
#'   in one of two ways:
#'
#'      1. `options(coloredrectangle.equalize_sizes=TRUE)`
#'      (priority); or
#'      2. `vertex.coloredrect.equalize_sizes=TRUE`,
#'      which is equivalent to `igraph::V(g)$coloredrect.equalize_sizes=TRUE`.
#'      in the latter case, only the first value is recognized.
#'
#'   * The behavior of `equalize_sizes` is described below:
#'
#'      * `equalize_sizes=FALSE` (default), nodes are exactly `size2`
#'      multiples of `coloredrect.ncol` and `coloredrect.nrow`.
#'      * `equalize_sizes=TRUE` or `1`, the shortest side of each node
#'      is scaled to `size2`. This option is useful to ensure nodes are never
#'      smaller than `size2` width or height, but can be larger.
#'      * `equalize_sizes=2`, the longest side of each node is scaled to
#'      `size2`. This option is useful to ensure nodes fit insides a square
#'      with `size2` sides.
#'
#'   * Nodes are drawn using vectorized processes where possible.
#'   However, the primary function `graphics::symbols()` only permits
#'   vectorized plotting for one `lwd`/`lty` combination at a time,
#'   so rendering is split into unique combinations of `lwd`/`lty`.
#'   Vectorized rendering is substantially faster than iterative rendering,
#'   and the majority of circumstances uses only one `lwd`/`lty` combination
#'   for all nodes. Line width is not an optimal way to convey a quantitative
#'   measurement, however it can be useful to highlight particular nodes
#'   of interest.
#'   * When `coloredrect.ncol` does not exist, the `ncol` will be defined
#'   by allowing up to two rows by default, enough to accommodate
#'   the number of colors in `coloredrect.color`.
#'   * When `coloredrect.ncol * coloredrect.nrow` will not fit all
#'   colors in `coloredrect.color`, the `coloredrect.ncol` will be
#'   extended to create enough positions to display all colors.
#'   * When `coloredrect.nrow` does not exist, it will use a value
#'   based upon `coloredrect.ncol` and the length of `coloredrect.color`.
#'
#' The values `coloredrect.color` and other variables described above refer to
#' `igraph` vertex attributes, and can be accessed for a given `igraph`
#' object `g` as follows:
#'
#' * `igraph::V(g)$coloredrect.color`
#' * `igraph::vertex_attr(g, "coloredrect.color")`
#' * or during plotting the value can be defined using the syntax
#' `igraph::plot(g, vertex.coloredrect.color=list(...))`
#'
#' Note that blank positions inside coloredrectangle nodes can be removed
#' via `removeIgraphBlanks()`, which also has the effect of modifying the
#' `coloredrect.ncol` and `coloredrect.nrow`, by applying the appropriate
#' logic.
#'
#' @return the plot function returns an invisible `list` of
#'    `data.frame` objects which were used to draw the rectangle objects.
#'    However the purpose of this function is the by-product that it
#'    draws rectangles onto an igraph graph.
#'    The clip function returns coordinates corresponding to the outer
#'    node edge based upon argument `end`:
#'    * `end="both"` returns 4 columns, the x,y 'edge from' and
#'    x,y 'edge to' coordinates, respectively
#'    * `end="from"` returns 2 columns, the x,y 'edge from' coordinates
#'    * `end="to"` returns 2 columns, the x,y 'edge to' coordinates
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
#' @examples
#' # prepare example igraph object
#' am <- matrix(ncol=5, nrow=5,
#'    data=0,
#'    dimnames=list(LETTERS[1:5], LETTERS[1:5]))
#' am[2:5, 1] <- 1;
#' g1 <- igraph::graph_from_adjacency_matrix(am)
#' igraph::graph_attr(g1, "layout") <- cbind(x=c(0, 1, 1, -1, -1),
#'    y=c(0, 1, -0.5, 0.5, -1))
#' colorset <- c("firebrick3", "gold", "deepskyblue",
#'    "mediumpurple3", "orchid1")
#' colorset <- c("firebrick3", "dodgerblue3");
#' vseq <- seq_len(igraph::vcount(g1));
#' vsizes <- c(3, 2, 2, 2, 1);
#' set.seed(1);
#' igraph::V(g1)$coloredrect.border <- lapply(vseq, function(i){
#'    sample(colorset,
#'       replace=TRUE,
#'       size=vsizes[i])
#' })
#' igraph::V(g1)$coloredrect.color <- lapply(vseq, function(i){
#'    jamba::alpha2col(alpha=0.5,
#'       igraph::V(g1)$coloredrect.border[[i]])
#' })
#' igraph::V(g1)$coloredrect.ncol <- c(2, 2, 2, 1, 1);
#' igraph::V(g1)$coloredrect.nrow <- c(2, 1, 1, 2, 1);
#' igraph::V(g1)$coloredrect.lwd <- rep(3, igraph::vcount(g1))
#' igraph::V(g1)$frame.lwd <- c(2, 1, 1, 1, 1);
#' igraph::V(g1)$frame.color <- "black"
#' igraph::V(g1)$size2 <- 10;
#' igraph::V(g1)$shape <- "coloredrectangle";
#'
#' plot(g1, vertex.label="")
#' title(font.main=1, line=1.5, main=paste0(
#'    "Each square is consistent size by vertex.size2\n"))
#' title(font.main=1, cex.main=1, line=0.5, main=paste0(
#'    "'vertex.coloredrect.equalize_sizes=FALSE' or\n",
#'    "'options(coloredrectangle.equalize_sizes=FALSE'"))
#'
#' # equalize shortest side to size2
#' plot(g1, vertex.label="", vertex.coloredrect.equalize_sizes=TRUE)
#' title(font.main=1, line=1.5, main=paste0(
#'    "The shortest side is fixed by vertex.size2\n"))
#' title(font.main=1, cex.main=1, line=0.5, main=paste0(
#'    "'vertex.coloredrect.equalize_sizes=TRUE' or\n",
#'    "'options(coloredrectangle.equalize_sizes=TRUE'"))
#'
#' # equalize longest side to size2
#' plot(g1, vertex.label="", vertex.coloredrect.equalize_sizes=2)
#' title(font.main=1, line=1.5, main=paste0(
#'    "The longest side is fixed by vertex.size2\n"))
#' title(font.main=1, cex.main=1, line=0.5, main=paste0(
#'    "'vertex.coloredrect.equalize_sizes=2' or\n",
#'    "'options(coloredrectangle.equalize_sizes=2'"))
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

   vertex.coloredrect.equalize_sizes <- getparam("coloredrect.equalize_sizes");
   if (length(vertex.coloredrect.equalize_sizes) == 0) {
      vertex.coloredrect.equalize_sizes <- FALSE;
   }
   vertex.coloredrect.equalize_sizes <- head(vertex.coloredrect.equalize_sizes,
      1);

   # get relevant global options
   verbose <- getOption("verbose", FALSE);
   equalize_sizes <- getOption("coloredrectangle.equalize_sizes",
      vertex.coloredrect.equalize_sizes)

   vertex.color <- getparam("color");

   vertex.frame.color <- getparam("frame.color");
   if (length(vertex.frame.color) == 0) {
      vertex.frame.color <- "grey30";
   }
   vertex.frame.lwd <- getparam("frame.lwd");
   if (length(vertex.frame.lwd) == 0) {
      vertex.frame.lwd <- rep(1,
         length.out=length(vertex.frame.color));
   }
   vertex.frame.lwd <- ifelse(
      is.na(vertex.frame.color) | jamba::col2alpha(vertex.frame.color) < 0.01,
      0,
      vertex.frame.lwd);
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
   } else if (is.list(vertex.coloredrect.border)) {
      if (all(lengths(vertex.coloredrect.border) == 0)) {
         vertex.coloredrect.border <- rep(
            list("grey30"),
            length.out=length(vertex.coloredrect.border));
      }
   }
   if (!is.list(vertex.coloredrect.border)) {
      vertex.coloredrect.border <- rep(
         list(vertex.coloredrect.border),
         length.out=length(vertex.coloredrect.color))
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
      vertex.coloredrect.ncol <- ceiling(lengths(vertex.coloredrect.color)/2);
   }
   vertex.coloredrect.nrow <- getparam("coloredrect.nrow");
   if (length(vertex.coloredrect.nrow) == 0) {
      vertex.coloredrect.nrow <- ceiling(
         lengths(vertex.coloredrect.color) /
            vertex.coloredrect.ncol);
   }
   # validate the number of colors will fit inside ncol,nrow
   vertex.coloredrect.ncells <- vertex.coloredrect.ncol * vertex.coloredrect.nrow;
   if (any(lengths(vertex.coloredrect.color) > vertex.coloredrect.ncells)) {
      # expand ncol to accommodate all colors
      vertex.coloredrect.ncol <- pmax(vertex.coloredrect.ncol,
         ceiling(lengths(vertex.coloredrect.color) / vertex.coloredrect.nrow))
   }

   vertex.coloredrect.lty <- getparam("coloredrect.lty");
   if (length(vertex.coloredrect.lty) == 0) {
      vertex.coloredrect.lty <- 1;
   }
   vertex.coloredrect.lwd <- getparam("coloredrect.lwd");
   if (length(vertex.coloredrect.lwd) == 0) {
      vertex.coloredrect.lwd <- 1;
   }
   vertex.coloredrect.lwd <- rep(vertex.coloredrect.lwd,
      length.out=length(vertex.coloredrect.color))
   vertex.coloredrect.lwd <- lapply(seq_along(vertex.coloredrect.border), function(i){
      iborder <- vertex.coloredrect.border[[i]];
      ilwd <- vertex.coloredrect.lwd[[i]];
      ifelse(
         is.na(iborder) | jamba::col2alpha(iborder) < 0.01,
         0,
         ilwd)
   })
   vertex.coloredrect.lwd.max <- sapply(vertex.coloredrect.lwd, max, na.rm=TRUE)

   # version 0.0.68.900: refactor size1, size2 calculations
   vertex.size1 <- getparam("size");
   vertex.size2 <- getparam("size2");
   if (length(vertex.size2) == 0) {
      vertex.size2 <- vertex.size1;
   }
   if (any(is.na(vertex.size2))) {
      vertex.size2 <- ifelse(is.na(vertex.size2), vertex.size1, vertex.size2);
   }
   # replace NA size1 with size2 / 2 (why divide by 2?)
   size1_na <- is.na(vertex.size1)
   if (any(size1_na)) {
      vertex.size1[size1_na] <- vertex.size2[size1_na] / 2;
   }
   # convert size to graph coordinates
   vertex.size1 <- rep(1/200 * vertex.size1, length.out=nrow(coords));
   vertex.size2 <- rep(1/200 * vertex.size2, length.out=nrow(coords));
   # Use size2 to define the size of each square
   # with consistent square size for all nodes
   vertex.size1 <- vertex.size2 * 5 * vertex.coloredrect.ncol;
   vertex.size2 <- vertex.size2 * 5 * vertex.coloredrect.nrow;

   # optionally adjust nodes for consistent max height/width
   if (equalize_sizes) {
      if (2 %in% equalize_sizes) {
         # method 2: the longest side is fixed at size2
         vertex.size1 <- vertex.size1 /
            pmax(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
         vertex.size2 <- vertex.size2 /
            pmax(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
      } else {
         # method 3: the shortest side is fixed at size2
         vertex.size1 <- vertex.size1 /
            pmin(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
         vertex.size2 <- vertex.size2 /
            pmin(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
      }
   }

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
      ncells <- nrow * ncol;
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
         ncellk <- ncells[k];
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
         size1v <- rep(min(diff(unique(x01))), numk);
         size2v <- rep(min(diff(unique(y01))), numk);
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
         # colk <- rep(col[[k]], length.out=length(x01v));
         colk <- rep(NA, length.out=length(x01v));
         colk[seq_along(col[[k]])] <- col[[k]];
         borderk <- rep(border[[k]], length.out=length(colk));
         # borderk <- rep(NA, length.out=length(x01v));
         # borderk[seq_along(border[[k]])] <- border[[k]];

         ltyk <- rep(lty[[k]], length.out=length(colk));
         lwdk <- rep(lwd[[k]], length.out=length(colk));
         if (TRUE %in% getOption("debug", FALSE)) {
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
         kDF <- data.frame(
            stringsAsFactors=FALSE,
            x=c(xk, x01v),
            y=c(yk, y01v),
            bg=c("transparent", colk),
            fg=c(head(frame.color[[k]],1),
               borderk),
            rectx=c(head(size1v,1)*ncolk, size1v),
            recty=c(head(size2v,1)*nrowk, size2v),
            lty=c(head(ltyk,1), ltyk),
            ## lwd=c(head(frame.lwd[[k]], 1)*2.5, wdk/4),
            # lwd=c(head(frame.lwd[[k]], 1) * 1, lwdk),
            lwd=c(frame.lwd[[k]], lwdk),
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
         rectDFi$lwd <- ifelse(is.na(rectDFi$fg), 1, rectDFi$lwd);
         rectDFi$lwd <- ifelse(is.na(rectDFi$lwd) | rectDFi$lwd == 0,
            1, rectDFi$lwd);
         # jamba::printDebug("rectDFi$lwd:");print(rectDFi$lwd);
         # jamba::printDebug("head(rectDFi):");print(head(rectDFi));
         if ("square" %in% rectDFi$rect_type[1]) {
            # jamba::printDebug("inner lwd: ", rectDFi$lwd);
            inner_coords <- adjust_rect_border(
               x=rectDFi$x,
               y=rectDFi$y,
               rectangles=as.matrix(rectDFi[,c("rectx","recty")]),
               type="inner",
               lwd=rectDFi$lwd);
            rectDFi$rectx <- inner_coords$rectangles[,1];
            rectDFi$recty <- inner_coords$rectangles[,2];
         } else {
            # jamba::printDebug("outer lwd: ", rectDFi$lwd);
            outer_coords <- adjust_rect_border(
               x=rectDFi$x,
               y=rectDFi$y,
               rectangles=as.matrix(rectDFi[,c("rectx","recty")]),
               type="outer",
               lwd=rectDFi$lwd);
            rectDFi$rectx <- outer_coords$rectangles[,1];
            rectDFi$recty <- outer_coords$rectangles[,2];
         }
         graphics::symbols(x=rectDFi$x,
            y=rectDFi$y,
            bg=rectDFi$bg,
            fg=rectDFi$fg,
            rectangles=as.matrix(rectDFi[,c("rectx","recty")]),
            add=TRUE,
            inches=FALSE,
            lty=rectDFi$lty,
            lwd=rectDFi$lwd,
            xpd=TRUE);
      }
      return(rectDFL);
   }

   ## For reference, the code below draws a single rectangle
   if (FALSE) {
      graphics::symbols(x=coords[, 1],
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


#' @details
#' Todo: The clip function should adjust the node boundary to account
#' for `frame.color` and `frame.lwd` when present, since the frame is
#' drawn as an outer border and slightly increases the size of the node.
#'
#' @family jam igraph shapes
#'
#' @rdname shape.coloredrectangle.plot
#'
#' @export
shape.coloredrectangle.clip <- function
(coords,
 el,
 params,
 end=c("both", "from", "to"))
{
   end <- match.arg(end)
   if (length(coords) == 0) {
      return(coords)
   }


   getparam <- function(pname) {
      p <- params("vertex", pname)
      p;
   }

   vertex.coloredrect.equalize_sizes <- getparam("coloredrect.equalize_sizes");
   if (length(vertex.coloredrect.equalize_sizes) == 0) {
      vertex.coloredrect.equalize_sizes <- FALSE;
   }
   vertex.coloredrect.equalize_sizes <- head(vertex.coloredrect.equalize_sizes,
      1);

   # get relevant global options
   verbose <- getOption("verbose", FALSE);
   equalize_sizes <- getOption("coloredrectangle.equalize_sizes", vertex.coloredrect.equalize_sizes)

   vertex.coloredrect.color <- getparam("coloredrect.color");
   if (length(vertex.coloredrect.color) == 0) {
      vertex.coloredrect.color <- getparam("pie.color");
   }
   vertex.coloredrect.ncol <- getparam("coloredrect.ncol");
   if (length(vertex.coloredrect.ncol) == 0) {
      vertex.coloredrect.ncol <- ceiling(lengths(vertex.coloredrect.color)/2);
   }
   vertex.coloredrect.nrow <- getparam("coloredrect.nrow");
   if (length(vertex.coloredrect.nrow) == 0) {
      vertex.coloredrect.nrow <- ceiling(
         lengths(vertex.coloredrect.color) /
            vertex.coloredrect.ncol);
   }

   # validate the number of colors will fit inside ncol,nrow
   vertex.coloredrect.ncells <- vertex.coloredrect.ncol * vertex.coloredrect.nrow;
   if (any(lengths(vertex.coloredrect.color) > vertex.coloredrect.ncells)) {
      # expand ncol to accommodate all colors
      vertex.coloredrect.ncol <- pmax(vertex.coloredrect.ncol,
         ceiling(lengths(vertex.coloredrect.color) / vertex.coloredrect.nrow))
   }

   # version 0.0.68.900: refactor size1, size2 calculations
   vertex.size1 <- getparam("size");
   vertex.size2 <- getparam("size2");
   if (length(vertex.size2) == 0) {
      vertex.size2 <- vertex.size1;
   }
   if (any(is.na(vertex.size2))) {
      vertex.size2 <- ifelse(is.na(vertex.size2), vertex.size1, vertex.size2);
   }
   # replace NA size1 with size2 / 2 (why divide by 2?)
   size1_na <- is.na(vertex.size1)
   if (any(size1_na)) {
      vertex.size1[size1_na] <- vertex.size2[size1_na] / 2;
   }
   # convert size to graph coordinates
   vertex.size1 <- rep(1/200 * vertex.size1, length.out=nrow(coords));
   vertex.size2 <- rep(1/200 * vertex.size2, length.out=nrow(coords));
   # Use size2 to define the size of each square
   # with consistent square size for all nodes
   vertex.size1 <- vertex.size2 * 2.5 * vertex.coloredrect.ncol;
   vertex.size2 <- vertex.size2 * 2.5 * vertex.coloredrect.nrow;

   # optionally adjust nodes for consistent max height/width
   if (equalize_sizes) {
      if (2 %in% equalize_sizes) {
         # method 2: the longest side is fixed at size2
         vertex.size1 <- vertex.size1 /
            pmax(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
         vertex.size2 <- vertex.size2 /
            pmax(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
      } else {
         # method 1: the shortest side is fixed at size2
         vertex.size1 <- vertex.size1 /
            pmin(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
         vertex.size2 <- vertex.size2 /
            pmin(vertex.coloredrect.ncol, vertex.coloredrect.nrow);
      }
   }

   vertex.size <- vertex.size1;
   if (verbose) {
      jamba::printDebug("vertex.size1: ", vertex.size1);
      jamba::printDebug("vertex.size2: ", vertex.size2);
   }

   rec.shift <- function(x0, y0, x1, y1, vsize, vsize2) {
      m <- (y0 - y1)/(x0 - x1)
      l <- cbind(
         x1 - vsize/m,
         y1 - vsize2,
         x1 - vsize,
         y1 - vsize * m,
         x1 + vsize2/m,
         y1 + vsize2,
         x1 + vsize,
         y1 + vsize * m)
      v <- cbind(
         x1 - vsize <= l[, 1] &
            l[, 1] <= x1 + vsize &
            y1 - vsize2 <= l[, 2] &
            l[, 2] <= y1 + vsize2,
         x1 - vsize <= l[, 3] &
            l[, 3] <= x1 + vsize &
            y1 - vsize2 <= l[, 4] &
            l[, 4] <= y1 + vsize2,
         x1 - vsize <= l[, 5] &
            l[, 5] <= x1 + vsize &
            y1 - vsize2 <= l[, 6] &
            l[, 6] <= y1 + vsize2,
         x1 - vsize <= l[, 7] &
            l[, 7] <= x1 + vsize &
            y1 - vsize2 <= l[, 8] &
            l[, 8] <= y1 + vsize2)
      d <- cbind(
         (l[, 1] - x0)^2 + (l[, 2] - y0)^2,
         (l[, 3] - x0)^2 + (l[, 4] - y0)^2,
         (l[, 5] - x0)^2 + (l[, 6] - y0)^2,
         (l[, 7] - x0)^2 + (l[, 8] - y0)^2)
      t(sapply(seq_len(nrow(l)), function(x) {
         d[x, ][!v[x, ]] <- Inf
         m <- which.min(d[x, ])
         l[x, c(m * 2 - 1, m * 2)]
      }))
   }
   if (end %in% c("from", "both")) {
      vsize <- if (length(vertex.size) == 1) {
         vertex.size
      } else {
         vertex.size[el[, 1]]
      }
      if (length(vertex.size2) == 1) {
         vsize2 <- vertex.size2
      } else {
         vsize2 <- vertex.size2[el[, 1]]
      }
      res <- res1 <- rec.shift(
         coords[, 3],
         coords[, 4],
         coords[, 1],
         coords[, 2],
         vsize,
         vsize2)
   }
   if (end %in% c("to", "both")) {
      if (length(vertex.size) == 1) {
         vsize <- vertex.size
      } else {
         vsize <- vertex.size[el[, 2]]
      }
      if (length(vertex.size2) == 1) {
         vsize2 <- vertex.size2
      } else {
         vsize2 <- vertex.size2[el[, 2]]
      }
      res <- res2 <- rec.shift(
         coords[, 1],
         coords[, 2],
         coords[, 3],
         coords[, 4],
         vsize,
         vsize2)
   }
   if (end == "both") {
      res <- cbind(res1, res2)
   }
   res
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

   # angle of rotation of the ellipse, clockwise in degrees
   vertex.ellipse.angle <- params("vertex", "ellipse.angle");
   if (length(vertex.ellipse.angle) == 0) {
      vertex.ellipse.angle <- 0;
   }
   if (length(vertex.ellipse.angle) != 1 && !is.null(v)) {
      vertex.ellipse.angle <- vertex.ellipse.angle[v];
   }

   drawEllipse(x=coords[,1],
      y=coords[,2],
      a=vertex.size,
      b=vertex.size/vertex.ellipse.ratio,
      col=vertex.color,
      border=vertex.frame.color,
      lwd=vertex.frame.width,
      angle=vertex.ellipse.angle,
      draw=TRUE);
}

#' clip function for igraph vertex shape ellipse
#'
#' clip function for igraph vertex shape ellipse
#'
#' This function defines the clipping function for custom igraph vertex
#' shape ellipse.
#'
#' @family jam igraph shapes
#'
#' @export
shape.ellipse.clip <- function
(coords,
 el,
 params,
 end=c("both", "from", "to"))
{
   end <- match.arg(end)
   if (length(coords) == 0) {
      return(coords)
   }
   # vertex size
   vertex.size <- 1/200 * params("vertex", "size");
   # vertex ellipse ratio, height:width
   vertex.ellipse.ratio <- params("vertex", "ellipse.ratio");
   if (length(vertex.ellipse.ratio) == 0) {
      vertex.ellipse.ratio <- 2;
   }
   vertex.ellipse.angle <- params("vertex", "ellipse.angle");
   if (length(vertex.ellipse.angle) == 0) {
      vertex.ellipse.angle <- 0;
   }
   # angle of rotation of the ellipse, clockwise in degrees
   vertex.ellipse.angle <- params("vertex", "ellipse.angle");
   if (length(vertex.ellipse.angle) == 0) {
      vertex.ellipse.angle <- 0;
   }
   vertex_count <- max(c(
      length(vertex.size),
      length(vertex.ellipse.ratio),
      length(vertex.ellipse.angle)))
   vertex.size <- rep(vertex.size, length.out=vertex_count);
   vertex.ellipse.ratio <- rep(vertex.ellipse.ratio, length.out=vertex_count);
   vertex.ellipse.angle <- rep(vertex.ellipse.angle, length.out=vertex_count);

   # ellipse width (a) and height (b) for un-rotated ellipses
   # a <- vertex.size;
   # b <- vertex.size / vertex.ellipse.ratio;
   if (end %in% c("from", "both")) {
      if (length(vertex.size) == 1) {
         vsize.from <- vertex.size
      } else {
         vsize.from <- vertex.size[el[, 1]]
      }
      if (length(vertex.ellipse.ratio) == 1) {
         vratio.from <- vertex.ellipse.ratio
      } else {
         vratio.from <- vertex.ellipse.ratio[el[, 1]]
      }
      if (length(vertex.ellipse.angle) == 1) {
         vangle.from <- vertex.ellipse.angle;
      } else {
         vangle.from <- vertex.ellipse.angle[el[, 1]]
      }
      a <- vsize.from;
      b <- vsize.from / vratio.from;

      x <- coords[, 1];
      y <- coords[, 2];
      x1 <- coords[, 3];
      y1 <- coords[, 4];
      angle <- rep(-vangle.from,
         length.out=nrow(coords));
      # rotation matrix adjustments
      co <- cos(-jamba::deg2rad(angle));
      si <- sin(-jamba::deg2rad(angle));
      # points rotated along with ellipse, around ellipse center
      x2 <- co * (x1 - x) - si * (y1 - y) + x;
      y2 <- si * (x1 - x) + co * (y1 - y) + y;
      # points relative to ellipse center
      x20 <- x2 - x;
      y20 <- y2 - y;
      # intercept with ellipse
      x_int <- (a * b) / sqrt(a^2 * y20^2 + b^2 * x20^2) * x20 + x
      y_int <- (a * b) / sqrt(a^2 * y20^2 + b^2 * x20^2) * y20 + y
      # un-rotate intersecting point
      co2 <- cos(-jamba::deg2rad(-angle));
      si2 <- sin(-jamba::deg2rad(-angle));
      x_int_rot <- co2 * (x_int - x) - si2 * (y_int - y) + x;
      y_int_rot <- si2 * (x_int - x) + co2 * (y_int - y) + y;
      x_from <- x_int_rot;
      y_from <- y_int_rot
   } else {
      x_from <- coords[, 1];
      y_from <- coords[, 2];
   }

   # end="to"
   if (end %in% c("to", "both")) {
      if (length(vertex.size) == 1) {
         vsize.to <- vertex.size
      } else {
         vsize.to <- vertex.size[el[, 2]]
      }
      if (length(vertex.ellipse.ratio) == 1) {
         vratio.to <- vertex.ellipse.ratio
      } else {
         vratio.to <- vertex.ellipse.ratio[el[, 2]]
      }
      if (length(vertex.ellipse.angle) == 1) {
         vangle.to <- vertex.ellipse.angle
      } else {
         vangle.to <- vertex.ellipse.angle[el[, 2]]
      }
      a <- vsize.to;
      b <- vsize.to / vratio.to;

      x <- coords[, 3];
      y <- coords[, 4];
      x1 <- coords[, 1];
      y1 <- coords[, 2];
      angle <- rep(-vangle.to,
         length.out=nrow(coords));
      # rotation matrix adjustments
      co <- cos(-jamba::deg2rad(angle));
      si <- sin(-jamba::deg2rad(angle));
      # points rotated along with ellipse, around ellipse center
      x2 <- co * (x1 - x) - si * (y1 - y) + x;
      y2 <- si * (x1 - x) + co * (y1 - y) + y;
      # points relative to ellipse center
      x20 <- x2 - x;
      y20 <- y2 - y;
      # intercept with ellipse
      x_int <- (a * b) / sqrt(a^2 * y20^2 + b^2 * x20^2) * x20 + x
      y_int <- (a * b) / sqrt(a^2 * y20^2 + b^2 * x20^2) * y20 + y
      # un-rotate intersecting point
      co2 <- cos(-jamba::deg2rad(-angle));
      si2 <- sin(-jamba::deg2rad(-angle));
      x_int_rot <- co2 * (x_int - x) - si2 * (y_int - y) + x;
      y_int_rot <- si2 * (x_int - x) + co2 * (y_int - y) + y;
      x_to <- x_int_rot;
      y_to <- y_int_rot
   } else {
      x_to <- coords[, 1];
      y_to <- coords[, 2];
   }
   if (end == "from") {
      res <- cbind(
         x_from,
         y_from);
   } else if (end == "to") {
      res <- cbind(
         x_to,
         y_to);
   } else if (end == "both") {
      res <- cbind(
         x_from,
         y_from,
         x_to,
         y_to);
   }
   res
}

#' custom igraph vertex shape jampie
#'
#' custom igraph vertex shape jampie
#'
#' This function is a vectorized replacement for plotting
#' vertex shape `"pie"` in much more efficient manner.
#'
#' It is substantially faster to use `shape.jampie.plot()` than
#' default igraph plotting, even for only 20 pie nodes, the speed
#' becomes even more dramatically faster for larger networks with
#' 200+ nodes. Minutes reduced to 1-2 seconds rendering time.
#'
#' Pie nodes with only one large
#' 100% wedge no longer display the small line from origin,
#' which is a change and improvement from default `igraph` rendering.
#'
#' Attribute `vertex.pie.border` can be used to draw a border around
#' each pie wedge, for each node. It should be a `list` with
#' `lengths(vertex.pie.border)` equal to `lengths(vertex.pie)`.
#' To disable, use `pie.border=NA` on the entire attribute, or individual
#' nodes.
#'
#' Attribute `vertex.frame.color` can be used to draw a single circular
#' border around the entire pie node. The `length(vertex.frame.color)`
#' should equal the number of nodes in the graph, for example
#' determined with `igraph::vcount(g)`.
#' Note that `frame.color` is drawn for each node after the pie
#' wedges, on top of `pie.border` if defined, so it is
#' recommended to use only one form of border for each node.
#'
#' Each pie node is drawn completely, in order: pie wedges including optional
#' `pie.border` outline for each pie wedge, then `frame.color`
#' around the entire node circle; then the next pie node is drawn.
#' This ordering ensures each entire pie node will overlap, or be
#' overlapped by other nodes, without artifacts of the `frame.color`
#' being shown on top of pie nodes that are otherwise beneath
#' visibility.
#'
#' To disable `pie.border` set to `NA` with `vertex.pie.border=NA`
#' or `V(g)[[2]]$pie.border <- NA`.
#'
#' To disable `frame.color` set to `NA` with `vertex.frame.color=NA`
#' or `V(g)[2]$frame.color <- NA`.
#'
#' @family jam igraph shapes
#'
#' @examples
#' # prepare example igraph object
#' am <- matrix(ncol=5, nrow=5,
#'    data=0,
#'    dimnames=list(LETTERS[1:5], LETTERS[1:5]))
#' am[2:5, 1] <- 1;
#' g1 <- igraph::graph_from_adjacency_matrix(am)
#' igraph::graph_attr(g1, "layout") <- cbind(x=c(0, 1, 1, -1, -1),
#'    y=c(0, 1, -0.5, 0.5, -1))
#' colorset <- c("firebrick3", "dodgerblue3");
#' vseq <- seq_len(igraph::vcount(g1));
#' vsizes <- c(3, 2, 2, 2, 1);
#' igraph::V(g1)$pie <- lapply(vseq, function(i){
#'    rep(1, vsizes[i])
#' })
#' set.seed(1);
#' igraph::V(g1)$pie.border <- lapply(vseq, function(i){
#'    sample(colorset,
#'       replace=TRUE,
#'       size=vsizes[i])
#' })
#' igraph::V(g1)$pie.color <- lapply(vseq, function(i){
#'    jamba::alpha2col(alpha=0.5,
#'       igraph::V(g1)$pie.border[[i]])
#' })
#' igraph::V(g1)$pie.lwd <- rep(5, igraph::vcount(g1))
#' igraph::V(g1)$frame.lwd <- c(2, 1, 1, 1, 1)*2;
#' igraph::V(g1)$frame.color <- NA
#' igraph::V(g1)$size <- c(45, 25, 25, 25, 25);
#' igraph::V(g1)$shape <- "jampie";
#'
#' par("mar"=c(2, 2, 3, 2))
#' igraph::plot.igraph(g1, vertex.label="",
#'    vertex.shape="pie", vertex.frame.color="grey45",
#'    main="shape='pie'\nigraph::plot.igraph()")
#'
#' jam_igraph(g1, vertex.label="",
#'    main="shape='jampie', frame.color=NA\njam_igraph()")
#'
#'
#' jam_igraph(g1, vertex.label="",
#'    vertex.frame.lwd=2,
#'    main="shape='jampie', frame.color='black'\njam_igraph()",
#'    vertex.shape="jampie", vertex.frame.color="black")
#'
#' jam_igraph(g1, vertex.label="", vertex.frame.color="black")
#'
#' # print pie.color data
#' print(igraph::vertex_attr(g1)$pie.color);
#'
#' @export
shape.jampie.plot <- function
(coords,
 v=NULL,
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
   vertex.frame.lwd <- getparam("frame.lwd")
   vertex.size <- rep(1/200 * getparam("size"),
      length.out=nrow(coords))
   vertex.pie <- getparam("pie")

   vertex.pie.color <- getparam("pie.color")

   # fill empty vertex.pie with uniform value=1
   if (length(vertex.pie) == 0 &&
         length(vertex.pie.color) > 0) {
      vertex.pie <- lapply(vertex.pie.color, function(i){
         rep(1, length(i))
      })
   }

   vertex.pie.border <- getparam("pie.border")
   if (!is.list(vertex.pie.border)) {
      if (length(vertex.pie.border) == 1) {
         vertex.pie.border <- rep(vertex.pie.border,
            length.out=length(vertex.pie))
      }
      vertex.pie.border <- rep(vertex.pie.border,
         lengths(vertex.pie));
      vertex.pie.border <- split(vertex.pie.border,
         rep(seq_along(vertex.pie),
            lengths(vertex.pie)))
   } else {
      vertex.pie.border <- rep(vertex.pie.border,
         length.out=length(vertex.pie))
   }

   vertex.pie.angle <- getparam("pie.angle")
   vertex.pie.density <- getparam("pie.density")

   vertex.pie.lty <- getparam("pie.lty")
   if (length(vertex.pie.lty) == 0) {
      vertex.pie.lty <- default_igraph_values()$vertex$pie.lty;
   }
   if (!is.list(vertex.pie.lty)) {
      if (length(vertex.pie.lty) == 1) {
         vertex.pie.lty <- rep(vertex.pie.lty,
            length.out=length(vertex.pie))
      }
      vertex.pie.lty <- rep(vertex.pie.lty,
         lengths(vertex.pie))
      vertex.pie.lty <- split(vertex.pie.lty,
         rep(seq_along(vertex.pie),
            lengths(vertex.pie)))
   } else {
      vertex.pie.lty <- rep(vertex.pie.lty,
         length.out=length(vertex.pie))
   }

   vertex.pie.lwd <- getparam("pie.lwd")
   if (length(vertex.pie.lwd) == 0) {
      vertex.pie.lwd <- default_igraph_values()$vertex$pie.lwd;
   }
   if (!is.list(vertex.pie.lwd)) {
      if (length(vertex.pie.lwd) == 1) {
         vertex.pie.lwd <- rep(vertex.pie.lwd,
            length.out=length(vertex.pie))
      }
      vertex.pie.lwd <- rep(vertex.pie.lwd,
         lengths(vertex.pie))
      vertex.pie.lwd <- split(vertex.pie.lwd,
         rep(seq_along(vertex.pie),
            lengths(vertex.pie)))
   } else {
      vertex.pie.lwd <- rep(vertex.pie.lwd,
         length.out=length(vertex.pie))
   }
   vertex.pie.lwd.max <- sapply(vertex.pie.lwd, max, na.rm=TRUE);

   # calculate one lwd to rule them all
   # because polygon() does not accept multiple lwd values in one call
   overall_lwd <- rep(
      pmax(
         vertex.frame.lwd,
         sapply(vertex.pie.lwd, max, na.rm=TRUE),
         na.rm=TRUE),
      length.out=nrow(coords))

   ## Convert for loop to lapply that obtains polygon coordinates
   if (TRUE) {
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
         } else if (length(vertex.pie.color) == 0) {
            NA
         } else {
            vertex.pie.color[[i]]
         }
         border <- if (length(vertex.pie.border) == 1) {
            vertex.pie.border[[1]]
         } else if (length(vertex.pie.border) == 0) {
            NA
         } else {
            vertex.pie.border[[i]]
         }
         frame.color <- if (length(vertex.frame.color) == 1) {
            vertex.frame.color[[1]]
         } else if (length(vertex.frame.color) == 0) {
            NA
         } else {
            vertex.frame.color[[i]]
         }
         frame.lwd <- if (length(vertex.frame.lwd) == 1) {
            vertex.frame.lwd[[1]]
         } else if (length(vertex.frame.lwd) == 0) {
            1
         } else {
            vertex.frame.lwd[[i]]
         }
         jam_mypie(
            x=coords[i, 1],
            y=coords[i, 2],
            pie,
            radius=vertex.size[i],
            edges=200,
            col=col,
            angle=na.omit(vertex.pie.angle[c(i, 1)])[1],
            density=na.omit(vertex.pie.density[c(i, 1)])[1],
            border=border,
            frame.color=frame.color,
            frame.lwd=frame.lwd,
            # frame.lwd=overall_lwd[i],
            lty=vertex.pie.lty[[i]],
            lwd=vertex.pie.lwd[[i]])
            # lty=na.omit(vertex.pie.lty[c(i, 1)])[1],
            # lwd=vertex.pie.lwd[c(i, 1)])
            # lwd=overall_lwd[i])
            # lwd = na.omit(vertex.pie.lwd[c(i, 1)])[1])
      }));

      # enable basic vectorized lwd
      poly_dfs <- split(poly_df, poly_df$lwd);
      for (poly_df in poly_dfs) {
         poly_x <- unlist(lapply(poly_df$x, function(ix1){
            c(ix1, NA)
         }))
         poly_y <- unlist(lapply(poly_df$y, function(ix1){
            c(ix1, NA)
         }))
         polygon(x=poly_x,
            y=poly_y,
            density=poly_df$density,
            angle=poly_df$angle,
            border=poly_df$border,
            col=poly_df$col,
            lty=poly_df$lty,
            lwd=ifelse(poly_df$lwd == 0, 1, poly_df$lwd),
            ljoin="mitre",
            lend="square")
      }
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


#' clip function for igraph vertex shape jampie
#'
#' clip function for igraph vertex shape jampie
#'
#' This function defines the clipping function for custom igraph vertex
#' shape jampie.
#'
#' @rdname shape.jampie.plot
#'
#' @family jam igraph shapes
#'
#' @export
shape.jampie.clip <- function
(coords,
 el,
 params,
 end=c("both", "from", "to"))
{
   end <- match.arg(end)
   if (length(coords) == 0) {
      return(coords)
   }
   vertex.pie <- params("vertex", "pie")
   vertex.size <- 1/200 * params("vertex", "size")
   # vmax <- max(c(el[,1], el[,2]));
   vmax <- length(vertex.pie);
   vertex.size <- rep(vertex.size, length.out=vmax);

   vertex.frame.lwd <- params("vertex", "frame.lwd")
   vertex.frame.color <- params("vertex", "frame.color")
   vertex.pie.lwd <- params("vertex", "pie.lwd")
   if (length(vertex.pie.lwd) == 0) {
      vertex.pie.lwd <- 1;
   }
   if (!is.list(vertex.pie.lwd)) {
      if (length(vertex.pie.lwd) == 1) {
         vertex.pie.lwd <- rep(vertex.pie.lwd,
            length.out=length(vertex.pie))
      }
      vertex.pie.lwd <- rep(vertex.pie.lwd,
         lengths(vertex.pie))
   }
   vertex.pie.lwd.max <- sapply(vertex.pie.lwd, max)

   if (length(vertex.frame.lwd) == 0) {
      vertex.frame.lwd <- 1
   }
   vertex.frame.lwd <- rep(vertex.frame.lwd,
      length.out=vmax);
   if (length(vertex.frame.color) == 0) {
      vertex.frame.color <- NA
   }
   vertex.frame.color <- rep(vertex.frame.color,
      length.out=vmax);
   vertex.frame.lwd <- ifelse(
      is.na(vertex.frame.color) | jamba::col2alpha(vertex.frame.color) < 0.01,
      0,
      vertex.frame.lwd);
   new.frame.lwd <- ifelse(vertex.frame.lwd > vertex.pie.lwd.max & vertex.pie.lwd.max > 0,
      vertex.pie.lwd.max,
      vertex.frame.lwd)
   vertex.frame.lwd <- new.frame.lwd;

   overall_lwd <- rep(
      pmax(
         vertex.frame.lwd,
         sapply(vertex.pie.lwd, max, na.rm=TRUE),
         na.rm=TRUE),
      length.out=length(vertex.size))

   # note that jam_mypie() adds radius to account for frame width for consistency
   inner_cex <- 1.01;
   fig_adj <- diff(par("usr")[1:2]) / par("pin")[1];
   inner_adjustment <- (vertex.pie.lwd.max / 96) * fig_adj;
   # inner_adjustment <- (overall_lwd / 96) * fig_adj;
   outer_adjustment <- (vertex.frame.lwd / 96) * fig_adj;
   radius_add <- inner_adjustment*1.01 - outer_adjustment;

   vertex.size <- vertex.size + radius_add;

   if (end %in% c("both", "from")) {
      phi <- atan2(coords[, 4] - coords[, 2],
         coords[, 3] - coords[, 1])
      if (length(vertex.size) == 1) {
         vsize.from <- vertex.size
      } else {
         vsize.from <- vertex.size[el[, 1]]
      }
      x_from <- coords[, 1] + vsize.from * cos(phi);
      y_from <- coords[, 2] + vsize.from * sin(phi);
   } else {
      x_from <- coords[, 1];
      y_from <- coords[, 2];
   }

   if (end %in% c("both", "to")) {
      # phi <- atan2(coords[, 4] - coords[, 2],
      #    coords[, 3] - coords[, 1])
      phi <- atan2(coords[, 4] - coords[, 2],
         coords[, 3] - coords[, 1]) + pi;
      # phi <- atan2(coords[, 2] - coords[, 4],
      #    coords[, 1] - coords[, 3])
      # r <- sqrt(
      #    (coords[, 3] - coords[, 1])^2 +
      #    (coords[, 4] - coords[, 2])^2)
      if (length(vertex.size) == 1) {
         vsize.to <- vertex.size
      } else {
         vsize.to <- vertex.size[el[, 2]]
      }
      # x_to <- coords[, 1] + (r - vsize.to*2) * cos(phi);
      # y_to <- coords[, 2] + (r - vsize.to*2) * sin(phi);
      x_to <- coords[, 3] + vsize.to * cos(phi);
      y_to <- coords[, 4] + vsize.to * sin(phi);
   } else {
      x_to <- coords[, 3];
      y_to <- coords[, 4];
   }

   if (end == "from") {
      res <- cbind(
         x_from,
         y_from);
   } else if (end == "to") {
      res <- cbind(
         x_to,
         y_to);
   } else if (end == "both") {
      res <- cbind(
         x_from,
         y_from,
         x_to,
         y_to);
   }
   res
}

#' Vectorized mypie() function for igraph vertex pie polygons
#'
#' Vectorized mypie() function for igraph vertex pie polygons
#'
#' This is a function called internally by `shape.jampie.plot()`,
#' not intended for direct use.
#'
#' See `reorder_igraph_nodes()` for visual examples.
#'
#' This function is a light rewrite of `igraph:::mypie()`, except
#' that this function determines polygon coordinates without
#' drawing them, instead returns the polygon coordinates to
#' the calling function `shape.jampie.plot()` which in turn
#' draws all polygons once using the vectorized approach
#' described for `graphics::polygon()`.
#'
#' See `shape.jampie.plot()` for more detail.
#'
#' To disable `pie.border` set to `NA` with `vertex.pie.border=NA`
#' or `V(g)[[2]]$pie.border <- NA`.
#'
#' To disable `frame.color` set to `NA` with `vertex.frame.color=NA`
#' or `V(g)[2]$frame.color <- NA`.
#'
#' @return `data.frame` with columns suitable for use in `polygon()`
#'    after each column is expanded to vector form. Columns `x` and `y`
#'    are stored `AsIs()` in `list` format, and should be combined
#'    with one `NA` value between each `numeric` vector. The `NA` values
#'    cause `polygon()` to draw a series of separate polygons in
#'    one vectorized step, much faster than drawing a series of
#'    polygons in an iterative loop.
#'
#' @family jam igraph shapes
#'
#' @param x,y `numeric` coordinate for the center of each `igraph` node.
#' @param values `numeric` vector of relative pie wedge sizes.
#' @param radius `numeric` radius of pie
#' @param edges `integer` number of edges to make a circle
#' @param col `character` vector of R colors to fill each pie wedge, in order
#' @param angle,density used to draw lines to fill each pie node, passed
#'    to `polygon()`
#' @param border `character` vector of R colors for each pie wedge border
#' @param frame.color `character` R color used around the entire pie circle.
#' @param lty `numeric` or `character` line type
#' @param init.angle `numeric` angle in degrees (0 to 360) where `0` is the
#'    top of the circle, proceeding clockwide.
#' @param inner_pie_border `logical` whether to apply `pie.border` colors
#'    only along the inside of each pie wedge polygon, so that adjacent
#'    colors will be seen beside each other without overlapping adjacent
#'    borders. This method is currently in development.
#' @param overclose_polygons `logical` indicating whether to close polygon
#'    coordinates with an extra couple points, to ensure the line join
#'    `"ljoin"` is properly called when drawing each polygon.
#' @param ... additional arguments are passed to `polygon()`
#'
#' @examples
#' # See reorder_igraph_nodes() for visual examples.
#'
#' @export
jam_mypie <- function
(x,
 y,
 values,
 radius,
 edges=200,
 col=NULL,
 angle=45,
 density=NULL,
 border=NULL,
 frame.color=NULL,
 frame.lwd=par("lwd"),
 lty=NULL,
 lwd=par("lwd"),
 init.angle=90,
 inner_pie_border=getOption("inner_pie_border", TRUE),
 inner_cex=1.01,
 xy_max=NULL,
 overclose_polygons=FALSE,
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
   col <- rep(jamba::rmNULL(nullValue=NA, col), length.out=nx)
   border <- rep(jamba::rmNULL(nullValue=NA, border), length.out=nx)
   frame.color <- rep(jamba::rmNULL(nullValue=NA, frame.color), length.out=nx)
   frame.lwd <- rep(jamba::rmNULL(nullValue=1, frame.lwd), length.out=nx)
   lty <- rep(jamba::rmNULL(nullValue=NA, lty), length.out=nx)
   lwd <- rep(jamba::rmNULL(nullValue=1, lwd), length.out=nx)
   lwd.max <- max(lwd, na.rm=TRUE);

   angle <- rep(jamba::rmNULL(nullValue=NA, angle), length.out=nx)
   density <- rep(jamba::rmNULL(nullValue=NA, density), length.out=nx)
   t2xy <- function(t, init.angle=90, radius=1) {
      t2p <- 2 * pi * t + init.angle * pi/180;
      list(
         x=radius * cos(t2p),
         y=radius * sin(t2p))
   }

   # extra adjustment for frame.color
   # when frame.lwd=0, the frame.color is set to NA to prevent
   # rendering any color
   frame.color <- ifelse(frame.lwd <= 0 | is.na(frame.lwd),
      NA,
      frame.color)

   # special adjustment for frame.lwd
   frame.lwd <- ifelse(
      is.na(frame.color) | jamba::col2alpha(frame.color) < 0.01,
      0,
      frame.lwd);
   new.frame.lwd <- ifelse(frame.lwd > lwd.max & lwd.max > 0,
      lwd.max,
      frame.lwd)
   frame.lwd <- new.frame.lwd;
   # frame border is adjusted based upon frame.lwd
   # however for frame border to be drawn with proper vectorization
   # order alongside pie wedges, it should share the same wedge lwd
   # since the wedge is drawn over the frame border, hiding all
   # but the effective frame.lwd.
   use.frame.lwd <- ifelse(lwd.max > 0 & frame.lwd < lwd.max,
      lwd.max, frame.lwd);

   if (TRUE %in% inner_pie_border) {
      # assume line width is "point"
      # each line width is 1/96 inches, so: 96 points / inch
      # 96 * (device_width inches) / x_width = points per x-coordinate
      # 1inch/96pt * x_width/dev_inch = x_width/pt
      # inner_adjustment <- (par("lwd")/96) * diff(par("usr")[1:2]) / par("pin")[1]
      inner_adjustment <- (lwd / 96) * diff(par("usr")[1:2]) / par("pin")[1]
      outer_adjustment <- (frame.lwd / 96) * diff(par("usr")[1:2]) / par("pin")[1]

      # radii <- radius + (inner_adjustment) - max(outer_adjustment);
      radii <- radius + (inner_adjustment)*(max(inner_adjustment)/inner_adjustment) - max(outer_adjustment);
      radius <- radius + max(inner_adjustment) - max(outer_adjustment);
      if (length(frame.color) == 0 || all(is.na(frame.color))) {
         # no frame is displayed, therefore make the node slightly larger
         # so the output is consistent for all nodes
         # radius <- radius + inner_adjustment;
      }
   }

   ## convert this loop to return matrix of polygon coordinates
   ## then stitch multiple polygons with empty row between each
   ## so polygon() can be called once on the whole set of polygons
   ##
   ## bonus points: for one-section pie graph, draw a circle without segment
   poly_df <- jamba::rbindList(lapply(seq_len(nx), function(i){
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(
         t=seq.int(values[i],
            values[i + 1],
            length.out = n),
         init.angle=init.angle,
         # radius=radius)
         radius=radii[i])
      if (nx == 1) {
         xvals <- (x + c(P$x));
         yvals <- (y + c(P$y));
      } else {
         xvals <- (x + c(P$x, 0));
         yvals <- (y + c(P$y, 0));
      }

      # original polygon
      polydf <- data.frame(
         stringsAsFactors=FALSE,
         x=I(list(xvals)),
         y=I(list(yvals)),
         density=density[i],
         angle=angle[i],
         border=border[i],
         col=col[i],
         lty=lty[i],
         lwd=lwd[i])

      # optional inner border polygon
      if (TRUE %in% inner_pie_border &&
            length(border[i]) == 1 &&
            !is.na(border[i])) {
         # slightly shrink each polygon
         polym <- cbind(x=xvals, y=yvals);
         polym <- rbind(polym, head(polym, 1));
         psf <- sf::st_polygon(list(polym));
         # jamba::printDebug("inner_adjustment:", inner_adjustment);
         use_adjustment <- -inner_adjustment[i]/(2*inner_cex);
         # use_adjustment <- -inner_adjustment[i]/(2*inner_cex) *
         #    (max(inner_adjustment) / inner_adjustment[i]);
         psf1 <- sf::st_buffer(x=psf,
            dist=use_adjustment);
         polym2 <- as(psf1, "matrix");
         if (TRUE %in% overclose_polygons) {
            xvals <- c(polym2[,1], head(polym2[,1], 2))
            yvals <- c(polym2[,2], head(polym2[,2], 2))
         } else {
            xvals <- head(polym2[,1], -1);
            yvals <- head(polym2[,2], -1);
         }
         # version 0.0.88.910: col=NA changed to col="#FFFFFF01"
         polydf2 <- data.frame(
            stringsAsFactors=FALSE,
            x=I(list(xvals)),
            y=I(list(yvals)),
            density=NA,
            angle=angle[i],
            border=border[i],
            # col=NA,
            col="#FFFFFF01",
            lty=lty[i],
            lwd=lwd[i])
         polydf[1, "border"] <- NA;
         polydf <- rbind(polydf, polydf2);
      }

      polydf
   }));

   if (length(frame.color) >= 1 && !all(is.na(frame.color))) {
      # new method re-calculates the full circle
      n <- edges
      P <- t2xy(
         t=seq.int(0,
            1,
            length.out = n),
         init.angle=init.angle,
         radius=radius)
      border_x <- (x + c(P$x));
      border_y <- (y + c(P$y));

      if (TRUE %in% inner_pie_border) {
         # slightly enlarge each polygon to draw the border outside
         polym <- cbind(x=border_x, y=border_y);
         polym <- rbind(polym, head(polym, 1));
         psf <- sf::st_polygon(list(polym));
         psf1 <- sf::st_buffer(x=psf,
            dist=outer_adjustment/(2*inner_cex));
         polym2 <- as(psf1, "matrix");
         border_x <- polym2[,1];
         border_y <- polym2[,2];
      }

      # previous method
      # border_x <- unlist(lapply(poly_df[[1]], function(i){head(i, -2)}));
      # border_y <- unlist(lapply(poly_df[[2]], function(i){head(i, -2)}));
      border_df <- data.frame(
         stringsAsFactors=FALSE,
         x=I(list(border_x)),
         y=I(list(border_y)),
         density=-1,
         angle=45,
         border=head(frame.color, 1),
         # border=frame.color[i],
         # col=NA,
         col="#FFFFFF01",
         lty=head(lty, 1),
         lwd=head(use.frame.lwd, 1))
         # lwd=head(frame.lwd, 1))
      return_df <- jamba::rbindList(list(
         border_df,
         poly_df
         ));
      return(return_df);
   }
   # more bonus points: optionally draw frame.color around each circle
   return(poly_df);
}
