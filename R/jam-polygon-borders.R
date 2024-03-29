
#' Adjust polygon border to inner or outer edge
#'
#' Adjust polygon border to inner or outer edge
#'
#' This function intends to allow displaying a polygon border
#' where the line itself only covers the inner or outer edge
#' of the border itself. As a result, each polygon border for
#' two adjacent polygons can be displayed without overlapping the other.
#'
#' @return `list` with elements `x` and `y` representing adjusted
#'    polygon coordinates. In the event that input `x` and `y`
#'    coordinates represent multiple polygons, the output `x` and `y`
#'    will also have `NA` values positioned between each polygon.
#'
#' @family jam plot functions
#'
#' @param x `numeric` polygon coordinates along the x-axis. Multiple
#'    polygons can be supplied by separating each polygon with `NA`,
#'    similar to how `polygon()` breaks coordinates into separate
#'    polygons at the position of missing values.
#' @param y `numeric` polygon coordinates along the x-axis.
#' @param type `character` string indicating where the border will
#'    be displayed:
#'    * `"inner"`: the polygon is shrunk by half the border line width
#'    * `"outer"`: the polygon is expanded by half the border line width
#' @param lwd_inch `numeric` conversion of line width `lwd` units per inch,
#'    by default `lwd=1` is defined as 1/96 inch.
#' @param inner_cex `numeric` minor adjustment to calculations
#' @param lwd `numeric` to define the intended line width for the resulting
#'    polygon border. By default it uses `par("lwd")`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' a <- head(seq(from=0, to=360, length.out=6), 3) + 90;
#' xl <- round(10 * c(cos(a/180 * pi) * 5, 0)) / 10;
#' y <- round(10 * sin(a/180 * pi) * 5 / 10);
#' y <- y - min(y);
#' yl <- c(y, y[3]);
#' xr <- xl * -1;
#' yr <- yl;
#' lwdl <- 6;
#' lwdr <- 6;
#' lwdo <- 5;
#' par("ljoin"="mitre",
#'    "lend"="round",
#'    "mar"=c(3, 3, 1, 1))
#' plot(NULL, bty="L",
#'    xlim=c(-5, 30), ylim=c(-15, 10),
#'    #xlim=c(20, 30), ylim=c(0, 10),
#'    type="n", asp=1, xlab="", ylab="")
#' abline(h=c(0, -15), lwd=0.5, col="black", xpd=FALSE)
#'
#' polygon(x=xl, y=yl,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xr, y=yr,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' jamba::drawLabels(x=0, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: Two adjacent polygon\nborders overlap.")
#'
#' xy2l <- adjust_polygon_border(x=xl, y=yl, type="inner", lwd=lwdl)
#' xy2r <- adjust_polygon_border(x=xr, y=yr, type="inner", lwd=lwdr)
#'
#' # draw properly resized polygon with inner border
#' polygon(x=xy2l$x + 25, y=xy2l$y,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 25, y=xy2r$y,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' jamba::drawLabels(x=25, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Adjust polygons so adjacent\ninner borders are both visible.")
#'
#' # draw border around the pentagon
#' # except notice the border also overlaps the original polygon borders
#' polygon(x=xy2l$x, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' xy <- list(x=c(head(xl, 3), rev(head(xr, 3))),
#'    y=c(head(yl, 3), rev(head(yr, 3))) - 15)
#' polygon(x=xy$x, y=xy$y, col=NA, border="gold", lwd=lwdo)
#' jamba::drawLabels(x=0, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: The overall polygon border\noverlaps the inner polygon borders.")
#'
#' # adjust for outer border
#' xy2 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo)
#' # draw properly resized polygon with inner border
#' polygon(x=xy2l$x + 25, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 25, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' # now draw outer border
#' polygon(x=xy2$x + 25, y=xy2$y,
#'    col=NA, border="gold", lwd=lwdo)
#' jamba::drawLabels(x=25, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Both the inner and outer\nborders are visible.")
#'
#'
#' # another example using different border line widths
#' lwdl <- 10;
#' lwdr <- 5;
#' lwdo <- 6;
#' plot(NULL, bty="L",
#'    xlim=c(-5, 30), ylim=c(-15, 10),
#'    type="n", asp=1, xlab="", ylab="")
#' abline(h=c(0, -15), lwd=0.5, col="black", xpd=FALSE)
#'
#' polygon(x=xl, y=yl,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xr, y=yr,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' jamba::drawLabels(x=0, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: Two adjacent borders\noverlap, and the line widths\nare not aligned.")
#'
#' xy2l <- adjust_polygon_border(x=xl, y=yl, type="inner", lwd=lwdl)
#' xy2r <- adjust_polygon_border(x=xr, y=yr, type="inner", lwd=lwdr)
#'
#' # draw properly resized polygon with inner border
#' polygon(x=xy2l$x + 25, y=xy2l$y,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 25, y=xy2r$y,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' jamba::drawLabels(x=25, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Inner borders are both visible\nand well-aligned.")
#'
#' # draw border around the pentagon
#' # except notice the border also overlaps the original polygon borders
#' polygon(x=xy2l$x, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' xy <- list(x=c(head(xl, 3), rev(head(xr, 3))),
#'    y=c(head(yl, 3), rev(head(yr, 3))) - 15)
#' polygon(x=xy$x, y=xy$y, col=NA, border="gold", lwd=lwdo)
#' jamba::drawLabels(x=0, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: The overall border still\noverlaps the inner polygon borders.")
#'
#' # adjust for outer border
#' xy2 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo)
#' # draw properly resized polygon with inner border
#' polygon(x=xy2l$x + 25, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 25, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' # now draw outer border
#' polygon(x=xy2$x + 25, y=xy2$y,
#'    col=NA, border="gold", lwd=lwdo)
#' jamba::drawLabels(x=25, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: The inner and outer\nborders are all visible.")
#'
#'
#' # final example showing how to draw an additional border
#' lwdl <- 10;
#' lwdr <- 5;
#' lwdo <- 6;
#' plot(NULL, bty="L",
#'    xlim=c(-5, 30), ylim=c(-15, 10),
#'    type="n", asp=1, xlab="", ylab="")
#' abline(h=c(0, -15), lwd=0.5, col="black", xpd=FALSE)
#'
#' # vanilla polygon border
#' polygon(x=xy$x, y=xy$y + 15,
#'    col="indianred3", border="darkorchid4", lwd=lwdo)
#' jamba::drawLabels(x=0, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: A vanilla polygon is drawn,\nwe want to add a border.")
#' # draw the same polygon
#' polygon(x=xy$x + 25, y=xy$y + 15,
#'    col="indianred3", border="darkorchid4", lwd=lwdo)
#' # add another outer border
#' xy3 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo,
#'    lwd_buffer=lwdo / 2)
#' polygon(x=xy3$x + 25, y=xy3$y + 15,
#'    col=NA, border="skyblue", lwd=lwdo)
#' xy4 <- adjust_polygon_border(x=xy$x, y=xy$y, type="inner", lwd=lwdo,
#'    lwd_buffer=lwdo / 2)
#' polygon(x=xy4$x + 25, y=xy4$y + 15,
#'    col=NA, border="yellow", lwd=lwdo)
#' jamba::drawLabels(x=25, y=10, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Use 'lwd_buffer' to adjust\nfor the existing border line width.")
#'
#' # more complicated example with inner and outer borders
#' # draw properly resized polygon with inner border
#' polygon(x=xy2l$x, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' # now draw outer border
#' xy5 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo)
#' polygon(x=xy5$x, y=xy5$y,
#'    col=NA, border="skyblue", lwd=lwdo)
#' jamba::drawLabels(x=0, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: A polygon with inner and\nouter borders is drawn,\nwe want to add more borders.")
#'
#' # draw same polygons
#' polygon(x=xy2l$x + 25, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 25, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' polygon(x=xy5$x + 25, y=xy5$y,
#'    col=NA, border="skyblue", lwd=lwdo)
#' # add inner borders
#' xy6 <- adjust_polygon_border(x=xl, y=yl, type="inner", lwd=lwdl,
#'    lwd_buffer=lwdl)
#' # note col="#FFFFFF01" is required for the final sharp polygon corner
#' polygon(x=c(xy6$x, xy6$x[c(0)]) + 25,
#'    y=c(xy6$y, xy6$y[c(0)]) - 15,
#'    col="#FFFFFF01", border="yellow", lwd=lwdl)
#' points(x=c(xy6$x, xy6$x[c(0)]) + 25,lend="butt",
#'    y=c(xy6$y, xy6$y[c(0)]) - 15,
#'    pch=as.character(seq_along(c(xy6$x, xy6$x[c(0)]))),
#'    col="black", cex=0.5)
#' xy7 <- adjust_polygon_border(x=xr, y=yr, type="inner", lwd=lwdr,
#'    lwd_buffer=lwdr)
#' polygon(x=c(xy7$x, xy7$x[c(0)]) + 25,
#'    y=c(xy7$y, xy7$y[c(0)]) - 15,
#'    col="#FFFFFF01", border="yellow", lwd=lwdr)
#' # new outer border
#' xy8 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo,
#'    lwd_buffer=lwdo)
#' polygon(x=xy8$x + 25, y=xy8$y,
#'    col="#FFFFFF01", border="red", lwd=lwdo)
#' jamba::drawLabels(x=25, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Use 'lwd_buffer' with\nproper line widths.\n(red and yellow)")
#'
#' # show the effect of col=NA
#' plot(NULL, bty="L",
#'    xlim=c(-10, 25), ylim=c(-15, -5),
#'    type="n", asp=1, xlab="", ylab="")
#' abline(h=c(0, -15), lwd=0.5, col="black", xpd=FALSE)
#' # note that adjusted borders must be recalculated for new xlim,ylim
#' xy2l <- adjust_polygon_border(x=xl, y=yl, type="inner", lwd=lwdl)
#' xy2r <- adjust_polygon_border(x=xr, y=yr, type="inner", lwd=lwdr)
#' xy5 <- adjust_polygon_border(x=xy$x, y=xy$y, type="outer", lwd=lwdo)
#' # draw polygons
#' polygon(x=xy2l$x, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' polygon(x=xy5$x, y=xy5$y,
#'    col=NA, border="skyblue", lwd=lwdo)
#' polygon(x=xy2l$x + 15, y=xy2l$y - 15,
#'    col="indianred3", border="navy", lwd=lwdl)
#' polygon(x=xy2r$x + 15, y=xy2r$y - 15,
#'    col="dodgerblue", border="darkorchid3", lwd=lwdr)
#' polygon(x=xy5$x + 15, y=xy5$y,
#'    col=NA, border="skyblue", lwd=lwdo)
#' # add inner borders
#' xy6 <- adjust_polygon_border(x=xl, y=yl, type="inner", lwd=lwdl*2,
#'    lwd_buffer=lwdl)
#' # note col="#FFFFFF01" is required for the final sharp polygon corner
#' polygon(x=c(xy6$x, xy6$x[c(0)]) + 0,
#'    y=c(xy6$y, xy6$y[c(0)]) - 15,
#'    col=NA, border="yellow", lwd=lwdl*2)
#' polygon(x=c(xy6$x, xy6$x[c(0)]) + 15,
#'    y=c(xy6$y, xy6$y[c(0)]) - 15,
#'    col="#FFFFFF01", border="yellow", lwd=lwdl*2)
#' jamba::drawLabels(x=0, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Problem: Using col=NA causes the polygon corner\nto be 'ended', not 'joined'.")
#' jamba::drawLabels(x=15, y=-5, adjPreset="top",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt="Solution: Use col='#FFFFFF01' to force\nproper sharp polygon corners.")
#'
#'
#'
#' # example with multiple polygons vectorized
#' # generate an n-pointed star
#' n <- 3;
#' stars <- lapply(seq_len(n), function(i){
#'    angle <- seq(from=(i-1)*pi*2/n, to=i*pi*2/n, length=3);
#'    angle <- c(angle, angle[c(1, 1)], NA) + (n-2) * pi/(n * 2);
#'    xp <- cos(angle) * c(1, 2, 1, 0, 1, 1)
#'    yp <- sin(angle) * c(1, 2, 1, 0, 1, 1)
#'    cbind(x=xp, y=yp)
#' })
#' star_colors <- colorjam::rainbowJam(n)
#' plot(NULL, bty="L",
#'    type="n", asp=1, xlim=c(-2, 7), ylim=c(-6, 2.2),
#'    xlab="", ylab="")
#' polygon(jamba::rbindList(stars),
#'    col=jamba::alpha2col(star_colors, alpha=0.5),
#'    border=star_colors,
#'    lwd=10)
#' points(jamba::rbindList(stars), pch=20)
#' jamba::drawLabels(x=0, y=-1.5, adjPreset="bottom",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt=paste0("Problem: Base R polygons\n",
#'       "have overlapping borders.\n",
#'       "Black points show each corner."))
#'
#' stars_xy <- jamba::rbindList(stars)
#' stars_xy[,1] <- stars_xy[,1] + 5;
#' stars_xy_inner <- adjust_polygon_border(x=stars_xy, lwd=10)
#' polygon(stars_xy_inner,
#'    col=jamba::alpha2col(star_colors, alpha=0.5),
#'    border=star_colors,
#'    lwd=10)
#' jamba::drawLabels(x=5, y=-1.5, adjPreset="bottom",
#'    drawBox=FALSE,
#'    labelCex=0.8,
#'    txt=paste0("Solution: Polygon inner borders\n",
#'       "display each border color.\n",
#'       "And the triangle is properly sized!)"))
#'
#' stars_xy_inner[[2]] <- stars_xy_inner[[2]] - 5;
#' polygon(x=jamba::rbindList(stars)[,1] + 5,
#'    y=jamba::rbindList(stars)[,2] - 5,
#'    col="white",
#'    border="black",
#'    lwd=2)
#' polygon(stars_xy_inner,
#'    col=jamba::alpha2col(star_colors, alpha=0.5),
#'    border=star_colors,
#'    lwd=10)
#' # the inner border adjustment can be vectorized
#' stars_xy_inner2 <- adjust_polygon_border(x=stars_xy_inner[[1]],
#'    y=stars_xy_inner[[2]], lwd=5, lwd_buffer=5)
#' polygon(stars_xy_inner2,
#'    col=NA,
#'    border="#FFFFEE88",
#'    lwd=5)
#'
#' @export
adjust_polygon_border <- function
(x,
 y=NULL,
 type="inner",
 lwd_inch=1/96,
 inner_cex=1.001,
 lwd=par("lwd"),
 lwd_buffer=0,
 verbose=FALSE,
 ...)
{
   #
   # check for open device
   if (!requireNamespace("sf")) {
      stop("The sf package is required.")
   }

   # confirm a device is already open
   if (length(dev.cur()) == 0) {
      stop("This function requires an open device.")
   }

   # validate arguments
   if (length(type) == 0) {
      type <- "inner"
   }
   if (!all(type %in% c("inner", "outer"))) {
      stop("type must contain either 'inner' or 'outer'")
   }
   if (length(y) == 0 && length(dim(x)) >= 2) {
      y <- x[,2];
      x <- x[,1];
   }
   if (length(lwd_inch) == 0) {
      lwd_inch <- 1/96;
   }
   if (length(inner_cex) == 0) {
      inner_cex <- 1.01;
   }
   if (length(lwd) == 0) {
      lwd <- 1;
   }

   # check coordinates for NA breaks between multiple polygons
   if (any(is.na(x))) {
      x_cuts <- sort(unique(c(0,
         which(is.na(x)) + 0.5,
         which(is.na(x)) - 0.5,
         length(x)+0.5)))
      x_list <- split(x, as.numeric(cut(seq_along(x), breaks=x_cuts)));
      y_list <- split(y, as.numeric(cut(seq_along(x), breaks=x_cuts)));
      x_keep <- sapply(x_list, function(i){!all(is.na(i))});
      x_set <- x_list[x_keep];
      y_set <- y_list[x_keep];
      lwd <- rep(lwd, length.out=length(x_set));
      lwd_buffer <- rep(lwd_buffer, length.out=length(x_set));
      type <- rep(type, length.out=length(x_set));
      xy_list <- lapply(seq_along(x_set), function(i){
         i_xy <- adjust_polygon_border(x=x_set[[i]],
            y=y_set[[i]],
            lwd=lwd[[i]],
            lwd_buffer=lwd_buffer[[i]],
            type=type[[i]],
            verbose=verbose,
            ...);
         rbind(do.call(cbind, i_xy),
            matrix(ncol=2, c(NA, NA)))
      });
      new_xy <- jamba::rbindList(xy_list);
      retvals <- list(
         x=new_xy[,1],
         y=new_xy[,2]);
      return(retvals);
   }

   # inner adjustment converts device size and coordinate range,
   # and assumes one line width (lwd=1) is 1/96 inch.
   inner_adjustment <- ((lwd + lwd_buffer * 2) * lwd_inch) *
      diff(par("usr")[1:2]) / par("pin")[1];
   if (verbose) {
      jamba::printDebug("adjust_polygon_border(): ",
         "adjustment (inner):",
         inner_adjustment)
   }

   adjust_sign <- -1;
   if ("outer" %in% type) {
      adjust_sign <- 1;
   }
   if (verbose) {
      jamba::printDebug("adjust_polygon_border(): ",
         "adjust_sign:",
         adjust_sign)
   }

   # slightly shrink each polygon
   # so that adding the outer border makes the width the correct total size
   polym <- cbind(x=x, y=y);
   polym <- rbind(polym,
      head(polym, 1));

   # create sf polygon
   psf <- sf::st_polygon(list(polym));

   # use sf to create buffer using adjusted border line width
   psf1 <- sf::st_buffer(x=psf,
      dist=adjust_sign * inner_adjustment/(2*inner_cex));

   # convert back to coordinate matrix
   polym2 <- as(psf1, "matrix");
   xvals <- head(polym2[,1], -1);
   yvals <- head(polym2[,2], -1)
   retvals <- list(
      x=xvals,
      y=yvals);
   return(retvals);
}

#' @rdname adjust_polygon_border
#'
#' `adjust_rect_border()` is a convenience wrapper to `adjust_polygon_border()`
#' intended for use with `graphics::symbols()` in the form
#' `graphics::symbols(x, y, rectangles)`.
#'
#' @export
adjust_rect_border <- function
(x,
 y,
 rectangles=NULL,
 type="inner",
 lwd=par("lwd"),
 lwd_buffer=0,
 lwd_inch=1/96,
 verbose=FALSE,
 ...)
{
   #
   # confirm a device is already open
   if (length(dev.cur()) == 0) {
      stop("This function requires an open device.")
   }

   if (length(lwd) == 0) {
      lwd <- 1;
   }
   if (length(lwd_buffer) == 0) {
      lwd_buffer <- 0;
   }
   if (length(type) == 0) {
      type <- "inner";
   }
   lwd <- rep(lwd, length.out=length(x));
   lwd_buffer <- rep(lwd_buffer, length.out=length(x));
   type <- rep(type, length.out=length(x));

   # attempt direct adjustment
   # inner adjustment converts device size and coordinate range,
   # and assumes one line width (lwd=1) is 1/96 inch.
   inner_adjustment <- ((lwd + lwd_buffer * 2) * lwd_inch) *
      diff(par("usr")[1:2]) / par("pin")[1];

   adjust_sign <- ifelse(type %in% "inner", -1, 1)
   rectangles[,1] <- rectangles[,1] + adjust_sign * inner_adjustment;
   rectangles[,2] <- rectangles[,2] + adjust_sign * inner_adjustment;
   retvals <- list(
      x=x,
      y=y,
      rectangles=rectangles,
      lwd=lwd)
   return(retvals);
}
