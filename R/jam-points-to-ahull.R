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
#' @param label_preset `character` (default `NULL`) indicating the side
#'    to place a label, when `label` is provided.
#'    Recognized values: `"bottom", "top", "left", "right"`.
#'    When `NULL` it detects the offset from the plot center.
#' @param label_adj_preset `character` (default label_preset) indicating
#'    the label adjustment relative to the position of the label. In
#'    most cases it should equal `label_preset`.
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
#'       label="alphahull",
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
 max_iterations=100,
 do_plot=FALSE,
 add=FALSE,
 hull_method=c("default",
    "alphahull",
    "igraph",
    "sf",
    "chull"),
 smooth=TRUE,
 shape=1/2,
 label=NULL,
 label.cex=1,
 label.x.nudge=0,
 label.y.nudge=0,
 label_preset=NULL,
 label_adj_preset=label_preset,
 verbose=FALSE,
 ...)
{
   # validate hull_method
   if (length(label_preset) > 0) {
      label_preset <- head(intersect(label_preset,
         c("bottom", "top", "left", "right")), 1)
   }
   if (length(label_adj_preset) > 0) {
      label_adj_preset <- head(intersect(label_adj_preset,
         c("bottom", "top", "left", "right")), 1)
   }
   if (TRUE %in% verbose) {
      jamba::printDebug("make_point_hull(): ",
         "label.y.nudge:",
         label.y.nudge);
   }
   hull_method <- match.arg(hull_method);
   if ("default" %in% hull_method) {
      if (jamba::check_pkg_installed("alphahull")) {
         hull_method <- "alphahull"
      } else {
         hull_method <- "chull"
      }
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "hull_method: ",
            hull_method)
      }
   }

   # ensure at least 3 points
   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }
   x <- unique(x);

   npoints <- nrow(x);
   if (npoints < 3) {
      xy_max <- max(apply(apply(x, 2, range, na.rm=TRUE), 2, diff, na.rm=TRUE));
      x <- jamba::rbindList(list(
         x,
         x + rnorm(prod(dim(x))) * xy_max / 100,
         x + rnorm(prod(dim(x))) * xy_max / 100
      ));
      ## Ensure polygon is "closed"?
      # x <- rbind(x, x[1, , drop=FALSE])
      ## Take only the first N rows?
      # x <- head(x, 4);
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "Expanding input points to ", nrow(x), " rows.");
         print(x);
      }
   }

   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }
   # default size
   xy_max <- max(apply(apply(x, 2, range, na.rm=TRUE), 2, diff, na.rm=TRUE));
   if (length(alpha) == 0 || any(is.na(alpha))) {
      alpha <- xy_max * 0.5;
      if (verbose) {
         jamba::printDebug("make_point_hull(): ",
            "calculated alpha: ", format(alpha, digits=3));
      }
   }


   #### iterate attempts with increasing alpha
   hull_failed <- FALSE;
   if (verbose) {
      jamba::printDebug("make_point_hull(): ",
         "Iterating alpha values up to ", max_iterations, " times.");
   }
   for (iteration in seq_len(max_iterations)) {
      mxys <- tryCatch({
         get_hull_data(x,
            verbose=(verbose > 1),
            buffer=buffer,
            alpha=alpha,
            expand=expand,
            npoints=npoints,
            hull_method=hull_method)
      }, error=function(e){
         if (verbose) {
            jamba::printDebug("Error:");print(e);
         }
         return(NULL);
      });
      if (length(mxys) > 0 && nrow(mxys) >= 3) {
         attr(mxys, "alphahull_alpha") <- alpha;
         break;
      }
      alpha <- alpha * 1.2;
   }
   if (length(mxys) == 0 || nrow(mxys) < 3) {
      return(NULL)
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
         print("mxys:");print(mxys);
         print("shape:");print(shape);
         print("x:");print(x);
         mxys;
      });
   }
   # Replace alpha attribute
   attr(mxys, "alphahull_alpha") <- alpha;

   if (do_plot && length(mxys) > 0) {
      hull_sf_exp <- sf::st_polygon(list(mxys));
      if (any(lwd == 0)) {
         border <- NA;
         lwd <- 1;
      }
      plot(hull_sf_exp,
         add=add,
         col=col,
         border=border,
         lwd=lwd,
         lty=lty,
         ...)
      if (length(label) > 0) {
         if (length(label) > 1) {
            label <- jamba::cPaste(label, sep="\n")
         }
         # determine the hull placement from plot center
         parusr <- par("usr");
         ## Todo: consider using point farthest from plot center
         hull_center <- c(x=mean(range(mxys[,1], na.rm=TRUE)),
            y=mean(range(mxys[,2], na.rm=TRUE)));
         plot_center <- c(x=mean(parusr[1:2], na.rm=TRUE),
            y=mean(parusr[3:4], na.rm=TRUE));
         hull_degree <- 360 - (round(
            jamba::rad2deg(atan2(
               y=(hull_center[2] - plot_center[2]),
               x=(hull_center[1] - plot_center[1]))) + 270) %% 360);
         if (length(label_preset) == 0) {
            if (hull_degree >= 315 || hull_degree <= 45) {
               label_preset <- "top";
            } else if (hull_degree >= 135 && hull_degree <= 225) {
               label_preset <- "bottom";
            } else if (hull_degree <= 135) {
               label_preset <- "right"
            } else {
               label_preset <- "left"
            }
         }
         if ("top" %in% label_preset) {
            use_mxys <- which.max(mxys[,2]);
         } else if ("bottom" %in% label_preset) {
            use_mxys <- which.min(mxys[,2]);
         } else if ("right" %in% label_preset) {
            use_mxys <- which.max(mxys[,1]);
         } else {
            use_mxys <- which.min(mxys[,1]);
         }
         if (length(label_adj_preset) == 0) {
            label_adj_preset <- label_preset;
         }
         # calculate preset position
         label_xy <- jamba::coordPresets(
            preset=label_preset,
            adjPreset=label_adj_preset);
         # now scale relative to polygon and not plot coordinates
         if (FALSE) {
            label_x <- jamba::normScale(x=label_xy[,1],
               low=parusr[1], high=parusr[2],
               from=min(mxys[,1], na.rm=TRUE),
               to=max(mxys[,1], na.rm=TRUE))
            label_y <- jamba::normScale(x=label_xy[,2],
               low=parusr[3], high=parusr[4],
               from=min(mxys[,2], na.rm=TRUE),
               to=max(mxys[,2], na.rm=TRUE))
         } else {
            label_x <- mxys[use_mxys, 1];
            label_y <- mxys[use_mxys, 2];
         }
         if (length(label.x.nudge) == 1 && !label.x.nudge == 0) {
            label_x <- label_x + label.x.nudge;
         }
         if (length(label.y.nudge) == 1 && !label.y.nudge == 0) {
            label_y <- label_y + label.y.nudge;
         }
         jamba::drawLabels(txt=label,
            labelCex=label.cex,
            boxCexAdjust=c(2.2, 2.2),
            x=label_x, y=label_y,
            adjX=label_xy$adjX,
            adjY=label_xy$adjY,
            drawSegments=FALSE,
            drawBox=FALSE)
      }
   }

   return(invisible(mxys))
}

#' Get data for alpha hull (internal)
#'
#' Get data for alpha hull (internal)
#'
#' This function is intended for internal use by `make_point_hull()`
#' and is exported for convenient re-use by other functions as
#' relevant.
#'
#' @param x `numeric` matrix with x,y coordinates
#' @param verbose `logical` indicating whether to print verbose output
#' @param hull_method `character` string with the preferred hull method.
#' @param alpha `numeric` passed to `alphahull:ashape()`
#' @param expand `numeric` used to define `alpha` based upon coordinate range.
#' @param ... additional arguments are ignored.
#'
#' @family jam utility functions
#'
#' @export
get_hull_data <- function
(x,
 verbose=FALSE,
 hull_method="alphahull",
 buffer=NULL,
 alpha=NULL,
 expand=0.1,
 ...)
{
   hiA <- list(alpha=NULL);
   if ("alphahull" %in% hull_method) {
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "hull_method: ", hull_method);
      }
      hiA <- alphahull::ashape(x,
         alpha=alpha)
      if (nrow(hiA$edges) <= 2) {
         # 2 or 0 edges means the hull is not a true polygon hull
         return(NULL)
      }
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "sdim(hiA):");
         print(jamba::sdim(hiA));
      }
      hiAedges <- data.frame(hiA$edges);
      if (verbose) {
         jamba::printDebug("get_hull_data(): ",
            "alpha: ", hiA$alpha);
      }

      #hiAedges$length <- sqrt((hiAedges$x1 - hiAedges$x2)^2 + (hiAedges$y1 - hiAedges$y2)^2);
      hiAedges$coord <- paste0(round(hiAedges$x1*10),
         "_",
         round(hiAedges$y1*10));
      hiAedges$coord2 <- paste0(round(hiAedges$x2*10),
         "_",
         round(hiAedges$y2*10));
      # assemble edges into polygon coordinates
      # - Note: There could be a better method, in hindsight not sure
      #   how this approach even works, haha.
      # - Edges are supplied as ind1, ind2, in apparently random order.
      #   Each unique index value can be assigned its corresponding coordinate
      #   Then the ind1,ind2 can be traversed in order.
      #   Example data:
      #   ind1 ind2
      #   1    2
      #   3    1
      #   3    2
      #   Define 1-to-2, then 2-to-3, then 3-to-1.
      #   Note each point must appear exactly twice.
      #
      #   This situation defines a "star" shape, one central point with spokes.
      #   (Should probably fail a validation check.)
      #   ind1 ind2
      #   4    2
      #   4    1
      #   3    4
      ind1 <- hiAedges$ind1;
      ind2 <- hiAedges$ind2;
      if (!all(table(c(ind1, ind2)) == 2)) {
         # points are not duplicated meaning hull is incorrect
         return(NULL)
      }
      # Insert validation that hull is one full piece?
      #
      if (verbose > 1) {
         jamba::printDebug("get_hull_data(): ",
            "hiAedges:");print(hiAedges);# debug
      }
      pointlist <- list();
      for (irow in seq_len(nrow(hiAedges))) {
         ind1i <- hiAedges$ind1[[irow]];
         ind2i <- hiAedges$ind2[[irow]];
         pointlist[[as.character(ind1i)]] <- c(hiAedges$x1[irow], hiAedges$y1[irow])
         pointlist[[as.character(ind2i)]] <- c(hiAedges$x2[irow], hiAedges$y2[irow])
      }
      # iterate each edge, convert to list of coordinates
      coordlist <- list()
      used_rows <- NULL;
      for (irow in seq_len(nrow(hiAedges))) {
         if (irow == 1) {
            ind1i <- hiAedges$ind1[[irow]];
            ind2i <- hiAedges$ind2[[irow]];
            # jamba::printDebug("row: ", 1, ", ind1i:", ind1i, ", ind2i:", ind2i);
            coordlist <- c(pointlist[as.character(ind1i)], pointlist[as.character(ind2i)])
            nextrow <- ind2i;
            used_rows <- c(used_rows, 1);
         } else {
            match1 <- setdiff(which(hiAedges$ind1 %in% nextrow), used_rows)
            match2 <- setdiff(which(hiAedges$ind2 %in% nextrow), used_rows)
            if (length(match1) == 1) {
               used_rows <- c(used_rows, match1);
               ind1i <- hiAedges$ind1[[match1]];
               ind2i <- hiAedges$ind2[[match1]];
               # jamba::printDebug("match1: ", match1, ", ind1i:", ind1i, ", ind2i:", ind2i);
               coordlist <- c(coordlist, pointlist[as.character(ind2i)])
               nextrow <- ind2i
            } else if (length(match2) == 1) {
               used_rows <- c(used_rows, match2);
               ind1i <- hiAedges$ind2[[match2]];
               ind2i <- hiAedges$ind1[[match2]];
               # jamba::printDebug("match2: ", match2, ", ind1i:", ind1i, ", ind2i:", ind2i);
               coordlist <- c(coordlist, pointlist[as.character(ind2i)])
               nextrow <- ind2i
            } else {
               if (verbose) {
                  jamba::printDebug("get_hull_data(): ",
                     "Hull had disconnected polygons.");
               }
               break;
            }
         }
      }
      ixy <- do.call(rbind, coordlist)
      if (FALSE) {
         hiAedges2 <- jamba::rbindList(lapply(seq_len(nrow(hiAedges)), function(irow){
            ix <- jamba::mixedSortDF(data.frame(
               x=c(hiAedges$x1[irow], hiAedges$x2[irow]),
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
         "head(ixy, 20):");
      print(head(ixy, 20));
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
