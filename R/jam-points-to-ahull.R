#' Make alpha hull from points
#'
#' Make alpha hull from points
#'
#' This function makes an alpha hull around points, calling
#' `alphahull::ashape()` then piecing together the somewhat
#' random set of outer edges into a coherent polygon.
#'
#' @return `numeric` matrix with polygon coordinates, where
#'    each polygon is separated by one row that contains `NA`
#'    values. This output is sufficient for vectorized plotting
#'    in base R graphics using `graphics::polygon()`.
#'
#' @family jam utility functions
#'
#' @param x `numeric` matrix with 2 columns that contains the
#'    coordinate of each point.
#' @param expand `numeric` value indicating the buffer width around
#'    each point, scaled based upon the total range of coordinates,
#'    used only when `buffer` is not supplied.
#' @param buffer `numeric` value indicating the absolute buffer width
#'    around each point. This value is used if provided, otherwise
#'    `expand` is used to derive a value for `buffer`.
#' @param alpha `numeric` value passed to `alphahull::ashape()` when
#'    hull_method is `"alphahull"`. This value determines the level of
#'    detail of the resulting hull.
#' @param seed `numeric` seed used with `set.seed()` to define reproducible
#'    behavior.
#' @param color,border `character` colors used when `do_plot=TRUE` to draw
#'    the resulting hull polygon.
#' @param lwd,lty line width and line type parameters, respectively.
#' @param max_iterations `integer` number of attempts to call
#'    `alphahull::ashape()` with varying values of `alpha`. Each iteration
#'    checks to confirm the resulting polygon includes all input points.
#' @param do_plot `logical` indicating whether to plot the polygon
#'    output.
#' @param add `logical` used when `do_plot=TRUE` to indicate whether
#'    the hull should be drawn onto an existing plot device, or whether
#'    to open a new plot prior to drawing the hull.
#' @param hull_method `character` string indicating the hull method to use:
#'    * `"default"` - will use `"alphahull"` if the `alphahull` R package
#'    is available.
#'    * `"alphahull"` - use `alphahull::ashape()` which is the preferred
#'    method, in fact the only available option that will allow a concave
#'    shape in the output.
#'    * `"igraph"` - calls hidden function `igraph:::convex_hull()` as
#'    used when drawing `mark.groups` around grouped nodes.
#'    * `"sf"` - calls `sf::st_convex_hull()`, with same effective output
#'    as `"igraph"`.
#'    * `"chull"` - calls `grDevices::chull()`, again with same effective
#'    output as `"igraph"`, but with benefit of not incurring additional
#'    R package dependencies.
#' @param smooth `logical` indicating whether to smooth the final polygon
#'    using `graphics::xspline()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' set.seed(12)
#' n <- 22
#' xy <- cbind(x=sample(seq_len(n), size=n, replace=TRUE),
#'    y=sample(seq_len(n), size=n, replace=TRUE));
#' xy <- rbind(xy, xy[1,,drop=FALSE])
#' x4 <- sf::st_multipoint(xy)
#'
#' if (jamba::check_pkg_installed("alphahull")) {
#'    plot(x4, col="red", pch=20, cex=3,
#'       main="hull_method='alphahull'")
#'    phxy <- make_point_hull(x=xy, expand=0.05, do_plot=TRUE,
#'       hull_method="alphahull",
#'       add=TRUE, xpd=TRUE)
#' }
#'
#' plot(x4, col="red", pch=20, cex=3,
#'    main="hull_method='chull'")
#' phxy2 <- make_point_hull(x=xy, expand=0.05, do_plot=TRUE,
#'    add=TRUE, verbose=TRUE, xpd=TRUE, hull_method="chull")
#'
#' plot(x4, col="red", pch=20, cex=3,
#'    main="hull_method='igraph'")
#' phxy2 <- make_point_hull(x=xy, expand=0.05, do_plot=TRUE,
#'    add=TRUE, verbose=TRUE, xpd=TRUE, hull_method="igraph")
#'
#' plot(x4, col="red", pch=20, cex=3,
#'    main="hull_method='sf'")
#' phxy2 <- make_point_hull(x=xy, expand=0.05, do_plot=TRUE,
#'    add=TRUE, verbose=TRUE, xpd=TRUE, hull_method="sf")
#'
#' @export
make_point_hull <- function
(x,
 expand=0.1,
 buffer=NULL,
 alpha=NULL,
 seed=123,
 col="#FF000033",
 border="#FF0000FF",
 lwd=2,
 lty=1,
 max_iterations=10,
 do_plot=FALSE,
 add=FALSE,
 hull_method=c("default",
    "alphahull",
    "igraph",
    "sf",
    "chull"),
 smooth=TRUE,
 shape=1/2,
 verbose=FALSE,
 ...)
{
   # validate hull_method
   hull_method <- match.arg(hull_method);
   if ("default" %in% hull_method) {
      if (jamba::check_pkg_installed("alphahull")) {
         hull_method <- "alphahull"
      } else {
         hull_method <- "chull"
      }
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "hull_method: ", hull_method)
      }
   }

   # ensure at least 3 points
   set.seed(seed);
   x <- unique(x);
   if (nrow(unique(x)) < 3) {
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "Expanding input points to ", nrow(x)*3, " rows.");
      }
      x <- jamba::rbindList(list(
         x,
         x + rnorm(prod(dim(x))),
         x + rnorm(prod(dim(x)))
      ));
      x <- head(x, 3);
   }

   set.seed(seed);
   # default size
   xy_max <- max(apply(apply(x, 2, range, na.rm=TRUE), 2, diff, na.rm=TRUE));
   if (length(alpha) == 0 || any(is.na(alpha))) {
      alpha <- xy_max * 0.5;
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "calculated alpha: ", format(alpha, digits=3));
      }
   }

   # custom internal function
   get_hull_data <- function
   (x,
    verbose=FALSE,
    hull_method="alphahull")
   {
      hiA <- list(alpha=NULL);
      if ("alphahull" %in% hull_method) {
         hiA <- alphahull::ashape(x, alpha=alpha)
         if (verbose) {
            jamba::printDebug("make_point_hull(): ",
               "sdim(hiA):");
            print(jamba::sdim(hiA));
         }
         hiAedges <- data.frame(hiA$edges);
         if (verbose) {
            jamba::printDebug("make_point_hull(): ",
               "alpha: ", hiA$alpha);
         }

         #hiAedges$length <- sqrt((hiAedges$x1 - hiAedges$x2)^2 + (hiAedges$y1 - hiAedges$y2)^2);
         hiAedges$coord <- paste0(round(hiAedges$x1*10),
            "_",
            round(hiAedges$y1*10));
         hiAedges$coord2 <- paste0(round(hiAedges$x2*10),
            "_",
            round(hiAedges$y2*10));
         hiAedges2 <- jamba::rbindList(lapply(seq_len(nrow(hiAedges)), function(irow){
            ix <- jamba::mixedSortDF(data.frame(x=c(hiAedges$x1[irow], hiAedges$x2[irow]),
               y=c(hiAedges$y1[irow], hiAedges$y2[irow])))
            data.frame(x1=ix$x[1], x2=ix$x[2],
               y1=ix$y[1], y2=ix$y[2])
         }))
         coordu <- unique(c(hiAedges$coord, hiAedges$coord2));
         #match(hiAedges$coord, hiAedges$coord2)
         #subset(hiAedges, coord %in%names(tcount(hiAedges$coord, 2)))
         hiAedges2 <- jamba::mixedSortDF(hiAedges, byCols=c("x1", "y1", "x2", "y2"));
         if (1 == 2) {
            plot(hiAedges2[,c("x1", "y1")], pch=c(0:9, LETTERS))
            plot(
               x=hiAedges2$x1+c(-3, 3),
               y=hiAedges2$y1+c(-3, 3), pch="1")
            points(
               x=hiAedges2$x2+c(-3, 3),
               y=hiAedges2$y2+c(3, -3), pch="2", col="red")
            segments(
               x0=hiAedges2$x1+c(-3, 3),
               y0=hiAedges2$y1+c(-3, 3),
               x1=hiAedges2$x2+c(-3, 3),
               y1=hiAedges2$y2+c(3, -3),
               lwd=c(2, 4),
               col=colorjam::rainbowJam(8))
         }

         # make a properly connected polygon
         icurr <- 1;
         ixy <- hiAedges[icurr, c("x1", "y1")]
         i2start <- hiAedges$coord[icurr]
         i2end <- hiAedges$coord2[icurr]
         for (i1 in seq_len(nrow(hiAedges))) {
            i3 <- head(
               which((hiAedges$coord %in% i2end | hiAedges$coord2 %in% i2end) &
                     seq_along(hiAedges$coord) != icurr),
               1);
            #printDebug("icurr:", icurr, ", i3:", i3);
            if (length(i3) == 0) {
               next;
            }
            if (hiAedges$coord[i3] == i2end) {
               ixy <- rbind(ixy, hiAedges[i3, c("x1", "y1")])
               i2start <- hiAedges$coord[i3]
               i2end <- hiAedges$coord2[i3]
            } else {
               ixy <- rbind(ixy,
                  jamba::renameColumn(hiAedges[i3, c("x2", "y2")],
                     from=c("x2", "y2"), to=c("x1", "y1")))
               i2start <- hiAedges$coord2[i3]
               i2end <- hiAedges$coord[i3]
            }
            icurr <- i3;
         }
      } else if ("igraph" %in% hull_method) {
         ixy1 <- igraph:::convex_hull(x)
         ixy <- rbind(ixy1$rescoords,
            ixy1$rescoords[1, , drop=FALSE]);
      } else if ("sf" %in% hull_method) {
         ixy <- as(sf::st_convex_hull(sf::st_multipoint(x)), "matrix");
      } else {
         # chull
         xy_index <- grDevices::chull(x)
         ixy <- x[c(xy_index, head(xy_index, 1)), , drop=FALSE];
      }
      # # convert to sf polygon
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "head(ixy):");
         print(head(ixy));
         jamba::printDebug("get_hull_data(): ",
            "class(ixy):",
            class(ixy));
      }
      hull_sf <- sf::st_polygon(
         list(as.matrix(ixy)))
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "class(hull_sf):",
            class(hull_sf));
      }

      # expand polygon using sp buffer
      xy_max <- max(apply(apply(x, 2, range, na.rm=TRUE), 2, diff, na.rm=TRUE));
      if (length(buffer) == 0) {
         buffer <- xy_max * expand;
      }
      # hull_sp_exp <- venndir::get_sp_buffer(hull_sp,
      #    sp_buffer=buffer,
      #    relative_size=FALSE)
      # sf buffer
      hull_sf_exp <- sf::st_buffer(x=hull_sf,
         dist=buffer)
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "class(hull_sf_exp):",
            class(hull_sf_exp));
      }
      # extract matrix data
      mxys <- as(hull_sf_exp, "matrix");
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "class(mxys):",
            class(mxys));
      }

      # check if all points are inside the polygon(s)
      sfpt <- sf::st_multipoint(x);
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "class(sfpt):",
            class(sfpt));
      }
      pts_inpoly <- sf::st_contains(
         x=sf::st_polygon(list(mxys)),
         y=sfpt,
         sparse=FALSE);
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "class(pts_inpoly):",
            class(pts_inpoly));
      }

      # if points are not all inside the polygon re-run with higher alpha
      if (verbose) {
         jamba::printDebug("pts_inpoly:",
            pts_inpoly, fgText=c("darkorange2", "dodgerblue"));
      }
      # if (!all(apply(pts_inpoly, 1, max, na.rm=TRUE) > 0)) {
      if (!all(pts_inpoly)) {
         hull_failed <- TRUE;
         if (verbose) {
            jamba::printDebug("Points were not inside the polygon");
         }
         mxys <- NULL;
      }
      if (length(mxys) > 0) {
         attr(mxys, "alpha") <- hiA$alpha;
         # attr(mxys, "sf") <- hull_sf_exp;
      }
      return(mxys);
   }
   ##################### end internal function


   #### iterate attempts with increasing alpha
   hull_failed <- FALSE;
   for (iteration in seq_len(max_iterations)) {
      mxys <- tryCatch({
         get_hull_data(x,
            verbose=(verbose > 1),
            hull_method=hull_method)
      }, error=function(e){
         if (verbose) {
            jamba::printDebug("Error:");print(e);
         }
         return(NULL);
      });
      if (length(mxys) > 0 && nrow(mxys) >= 3) {
         break;
      }
      alpha <- alpha * 1.5;
   }
   if (length(mxys) > 0) {
      attr(mxys, "iteration") <- iteration;
   }

   # optionally smooth the final polygon points
   if (TRUE %in% smooth) {
      mxys <- tryCatch({
         mxys1 <- do.call(cbind, graphics::xspline(
            x=mxys,
            shape=shape,
            open=FALSE,
            draw=FALSE));
         rbind(mxys1, head(mxys1, 1));
      }, error=function(e){
         print(e)
         mxys;
      });
   }

   if (do_plot && length(mxys) > 0) {
      hull_sf_exp <- sf::st_polygon(list(mxys));
      plot(hull_sf_exp,
         add=add,
         col=col,
         border=border,
         lwd=lwd,
         lty=lty,
         ...)
   }

   return(invisible(mxys))
}

