
# various igraph internal functions that are required,
# and which CRAN does not permit to be called directly.

#' Parse igraph plot params
#'
#' This function mimics the internal function `igraph:::i.parse.plot.params()`.
#'
#' @family jam igraph internal functions
#'
parse_igraph_plot_params <- function
(graph,
 params)
{
   # initialize empty list
   p <- list(vertex=list(),
      edge=list(),
      plot=list())

   # iterate each element of the list
   for (n in names(params)) {
      if (substr(n, 1, 7) == "vertex.") {
         nn <- substring(n, 8)
         p[["vertex"]][[nn]] <- params[[n]]
      } else if (substr(n, 1, 5) == "edge.") {
         nn <- substring(n, 6)
         p[["edge"]][[nn]] <- params[[n]]
      } else {
         p[["plot"]][[n]] <- params[[n]]
      }
   }
   # create wrapper function to return each relevant value
   func <- function
   (type,
    name,
    range=NULL,
    dontcall=FALSE)
   {
      if (!type %in% names(p)) {
         stop("Invalid plot option type")
      }
      ret <- function() {
         v <- p[[type]][[name]]
         if (is.function(v) && !dontcall) {
            v <- v(graph)
         }
         if (length(range) == 0) {
            return(v)
         } else {
            if (length(v) == 1) {
               return(rep(v,
                  length(range)))
            } else {
               return(rep(v,
                  length=max(range) + 1)[[range + 1]])
            }
         }
      }
      if (name %in% names(p[[type]])) {
         return(ret())
      } else {
         if (type == "vertex" && name %in% igraph::vertex_attr_names(graph)) {
            p[[type]][[name]] <- igraph::vertex_attr(graph, name)
            return(ret())
         } else if (type == "edge" && name %in% igraph::edge_attr_names(graph)) {
            p[[type]][[name]] <- igraph::edge_attr(graph, name)
            return(ret())
         } else if (type == "plot" && name %in% igraph::graph_attr_names(graph)) {
            p[[type]][[name]] <- igraph::graph_attr(graph, name)
            return(ret())
         } else {
            n <- paste0(type, ".", name);
            v <- igraph::igraph_opt(n)
            if (!is.null(v)) {
               p[[type]][[name]] <- v
               return(ret())
            }
            p[[type]][[name]] <- default_igraph_values()[[type]][[name]]
            return(ret())
         }
      }
   }
   return(func)
}

#' Default igraph parameter values
#'
#' Default igraph parameter values
#'
#' @family jam igraph internal functions
#'
#' @return `list` of default igraph plotting and data parameters,
#'    including `"plot"`, `"vertex"`, and `"edge"`.
#'
#' @export
default_igraph_values <- function
()
{
   #
   paramnames <- ls(igraph:::i.default.values);
   names(paramnames) <- paramnames;
   paramvalues <- lapply(paramnames, function(paramname){
      get(paramname, envir=igraph:::i.default.values)
   })

   # plot
   default_plot_params <- list(
      palette=c("#E69F00",
         "#56B4E9",
         "#009E73",
         "#F0E442",
         "#0072B2",
         "#D55E00",
         "#CC79A7",
         "#999999"),
      layout=function(graph, dim=2) {
         if ("layout" %in% igraph::graph_attr_names(graph)) {
            lay <- igraph::graph_attr(graph, "layout")
            if (is.function(lay)) {
               lay <- lay(graph,
                  ...)
            } else {
               lay
            }
         } else if (all(c("x", "y") %in% igraph::vertex_attr_names(graph))) {
            if ("z" %in% igraph::vertex_attr_names(graph)) {
               lay <- cbind(igraph::V(graph)$x,
                  igraph::V(graph)$y,
                  igraph::V(graph)$z)
            } else {
               lay <- cbind(igraph::V(graph)$x,
                  igraph::V(graph)$y)
            }
         } else if (igraph::vcount(graph) < 1000) {
            lay <- igraph::layout_with_fr(graph,
               dim=dim)
         } else {
            lay <- igraph::layout_with_drl(graph,
               dim=dim)
         }
         # new in multienrichjam: rownames use V(graph)$name
         if ("name" %in% igraph::vertex_attr_names(graph)) {
            rownames(lay) <- igraph::vertex_attr(graph, "name")
         }
         lay
      },
      margin=c(0, 0, 0, 0),
      rescale=FALSE,
      asp=1,
      frame=FALSE,
      main=function(graph) {
         if (igraph::igraph_opt("annotate.plot")) {
            n <- graph$name[1]
            n
         } else {
            ""
         }
      },
      sub="",
      xlab=function(graph) {
         if (igraph::igraph_opt("annotate.plot")) {
            paste(igraph::vcount(graph), "vertices,",
               igraph::ecount(graph), "edges")
         } else {
            ""
         }
      },
      ylab= ""
   );

   # vertex
   default_vertex_params <- list(
      color=1,
      size=15,
      size2=15,
      label=function(graph, labels=NULL) {
         if (is.null(labels)) {
            if ("name" %in% igraph::vertex_attr_names(graph)) {
               labels <- igraph::vertex_attr(graph, "name")
            } else {
               labels <- seq_len(igraph::vcount(graph))
            }
         }
         labels
      },
      label.degree=-0.785,
      label.color="darkblue",
      label.dist=0,
      # label.family="serif",
      label.family="sans",
      label.font=1,
      label.cex=1,
      frame.color="black",
      frame.lwd=1,
      shape="circle",
      pie=1,
      pie.color=list("white",
         "lightblue",
         "mistyrose",
         "lightcyan",
         "lavender",
         "cornsilk"),
      pie.border=list("grey30"),
      pie.angle=45,
      pie.density=-1,
      pie.lty=1,
      pie.lwd=1,
      coloredrect.lwd=1,
      coloredrect.border="grey30"
   );

   # edge
   default_edge_params <- list(
      color="darkgrey",
      label=function(graph, edge.labels=NULL) {
         if (length(edge.labels) == 0) {
            edge.labels <- rep(NA,
               igraph::ecount(graph))
         }
         edge.labels
      },
      lty=1,
      width=1,
      loop.angle=0,
      loopangle2=0,
      label.family="serif",
      label.font=1,
      label.cex=1,
      label.color="darkblue",
      label.x=NULL,
      label.y=NULL,
      arrow.size=1,
      arrow.mode=function(graph, arrow.mode=NULL){
         if (is.character(arrow.mode) &&
               length(arrow.mode) == 1 &&
               substr(arrow.mode, 1, 2) == "a:") {
            arrow.mode <- igraph::vertex_attr(graph,
               substring(arrow.mode, 3))
         }
         if (is.character(arrow.mode)) {
            tmp <- numeric(length(arrow.mode))
            tmp[arrow.mode %in% c("<", "<-")] <- 1
            tmp[arrow.mode %in% c(">", "->")] <- 2
            tmp[arrow.mode %in% c("<>", "<->")] <- 3
            arrow.mode <- tmp
         }
         if (length(arrow.mode) == 0) {
            if (igraph::is_directed(graph)) {
               arrow.mode <- 2
            } else {
               arrow.mode <- 0
            }
         }
         arrow.mode
      },
      curved=function(graph, start=0.5){
         el <- apply(igraph::as_edgelist(graph, names=FALSE), 1, paste, collapse=":")
         ave(rep(NA, length(el)), el, FUN=function(x) {
            if (length(x) == 1) {
               return(0)
            } else {
               return(seq(from=-start,
                  to=start,
                  length=length(x)))
            }
         })
      },
      arrow.width=1
   );
   default_igraph_params <- list(
      plot=default_plot_params,
      vertex=default_vertex_params,
      edge=default_edge_params);
}

#' Render igraph arrows
#'
#' Render igraph arrows
#'
#' This function is a mimic of the internal `igraph:::igraph.Arrows()`
#' which is not permitted to be called directly for CRAN-approved
#' R packages.
#'
#' @family jam igraph internal functions
#'
#' @param x1,y1,x2,y2 `numeric` coordinates for initial and final x and y
#'    coordinates.
#' @param code `integer` indicating the position of arrow:
#'    * `code=1` arrow is positioned on the line end
#'    * `code=2` arrow is positioned on the line start
#'    * `code=3` arrow is positioned on both ends of the line
#' @param size `numeric` scaled size of the arrow head, which is applied to
#'    both the length and width of the arrow head.
#' @param width `numeric` scalar for the arrow head width, which is only
#'    applied to the relative arrow width.
#' @param open `logical` indicating whether the arrow head should be a filled
#'    polygon, otherwise only the outer "V" lines are drawn.
#' @param sh.adj `numeric` adjustment for segment length, where:
#'    * `sh.adj=0`  will extend the edge line (underneath the arrow head)
#'    to the end of the line
#'    * `sh.adj=1` will extend the edge line only to the base of the arrow head
#'    * `sh.adj=1.1` will leave a gap approximately 10% the arrow head length,
#'     between the edge line and the start of the arrow head.
#' @param sh.lwd `numeric` line width of main edge line
#' @param sh.col `character` color of main edge line
#' @param sh.lty `numeric` line type of main edge line
#' @param h.col,h.col.bo `character` arrow head color and arrow head border
#'    color, respectively.
#' @param h.lwd `numeric` arrow head line width
#' @param h.lty `numeric` arrow head line type
#' @param arrows_only `logical` indicating whether to draw only arrows,
#'    or `arrows_only=FALSE` to draw arrows and edge lines.
#' @param curved `logical` indicating whether to draw curved edges
#' @param verbose `logical` indicating whether to print verbose output.
#'
#' @examples
#' plot(NULL, xlim=c(-3, 3), ylim=c(-4, 4), type="n", xlab="", ylab="", bty="n")
#' jam_igraph_arrows(-2, 3, 2, 3, code=1, open=FALSE, sh.col="blue", sh.lwd=2)
#' jam_igraph_arrows(-2, 2, 2, 2, code=2, open=FALSE, sh.col="red", sh.lwd=2)
#' jam_igraph_arrows(-2, 1, 2, 1, code=3, open=FALSE, sh.col="gold", sh.lwd=2)
#' jam_igraph_arrows(-2, 0, 2, 0, code=3, arrows_only=TRUE, open=FALSE, sh.col="purple4", sh.lwd=2)
#'
#' jam_igraph_arrows(-2, -1, 2, -1, code=1, open=FALSE, sh.col="blue", h.col="#FF000055", sh.lwd=2, size=2, sh.adj=0.1)
#' jam_igraph_arrows(-2, -2, 2, -2, code=1, open=FALSE, sh.col="blue", h.col="#FF000055", sh.lwd=2, size=2, sh.adj=1.1)
#' jam_igraph_arrows(-2, -3, 2, -3, code=2, open=FALSE, sh.col="blue", h.col="#FF000055", sh.lwd=2, size=2, sh.adj=1.1)
#' jam_igraph_arrows(-2, -4, 2, -4, code=3, open=FALSE, sh.col="blue", h.col="#FF000055", sh.lwd=2, size=2, sh.adj=1.1)
#' text(x=rep(0, 8), y=seq(from=3, to=-4)+0.2,
#'    labels=c("code=1",
#'       "code=2",
#'       "code=3",
#'       "code=3, arrows_only=TRUE",
#'       "code=1, size=2, sh.adj=0.1",
#'       "code=1, size=2, sh.adj=1.1",
#'       "code=2, size=2, sh.adj=1.1",
#'       "code=3, size=2, sh.adj=1.1"))
#'
#' @export
jam_igraph_arrows <- function
(x1,
 y1,
 x2,
 y2,
 code=2,
 size=1,
 width=1,
 open=TRUE,
 sh.adj=0.1,
 sh.lwd=1,
 sh.col=if (is.R()) par("fg") else 1,
 sh.lty=1,
 h.col=sh.col,
 h.col.bo=sh.col,
 h.lwd=sh.lwd,
 h.lty=sh.lty,
 arrows_only=FALSE,
 curved=FALSE,
 verbose=FALSE)
{
   cin <- head(size * par("cin")[2], 1);
   width <- head(width * (1.2/4/cin), 1);

   if (verbose) {
      jamba::printDebug("jam_igraph_arrows()");
   }

   if (is.R()) {
      uin <- 1/xyinch()
   } else {
      uin <- par("uin")
   }
   x <- sqrt(seq(0, cin^2, length=floor(35 * cin) + 2))

   delta <- sqrt(h.lwd) * par("cin")[2] * 0.005
   x.arr <- c(-rev(x), -x)
   wx2 <- width * x^2
   y.arr <- c(-rev(wx2 + delta), wx2 + delta)
   deg.arr <- c(atan2(y.arr, x.arr), NA)
   r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
   bx1 <- x1
   bx2 <- x2
   by1 <- y1
   by2 <- y2
   lx <- length(x1)
   r.seg <- rep(cin * sh.adj, lx)
   theta1 <- atan2(
      (y1 - y2) * uin[2],
      (x1 - x2) * uin[1])
   th.seg1 <- theta1 + rep(atan2(0, -cin), lx)
   theta2 <- atan2(
      (y2 - y1) * uin[2],
      (x2 - x1) * uin[1])
   th.seg2 <- theta2 + rep(atan2(0, -cin), lx)
   x1d <- y1d <- x2d <- y2d <- 0
   # if (code %in% c(1, 3)) {
   if (code %in% c(2, 3)) {
      x2d <- r.seg * cos(th.seg2)/uin[1]
      y2d <- r.seg * sin(th.seg2)/uin[2]
   }
   # if (code %in% c(2, 3)) {
   if (code %in% c(1, 3)) {
      x1d <- r.seg * cos(th.seg1)/uin[1]
      y1d <- r.seg * sin(th.seg1)/uin[2]
   }

   # edge line drawn between arrow heads
   if (is.logical(curved) && all(!curved) ||
         is.numeric(curved) && all(!curved)) {
      # straight line between arrow heads
      if (!arrows_only) {
         segments(
            x1 + x1d,
            y1 + y1d,
            x2 + x2d,
            y2 + y2d,
            lwd=sh.lwd,
            col=sh.col,
            lty=sh.lty)
      }
      phi <- atan2(
         y1 - y2,
         x1 - x2)
      r <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
      lc.x <- x2 + 2/3 * r * cos(phi)
      lc.y <- y2 + 2/3 * r * sin(phi)
   } else {
      if (is.numeric(curved)) {
         lambda <- curved
      } else {
         lambda <- as.logical(curved) * 0.5
      }
      lambda <- rep(lambda, length.out = length(x1))
      c.x1 <- x1 + x1d
      c.y1 <- y1 + y1d
      c.x2 <- x2 + x2d
      c.y2 <- y2 + y2d
      midx <- (x1 + x2)/2
      midy <- (y1 + y2)/2
      spx <- midx - lambda * 1/2 * (c.y2 - c.y1)
      spy <- midy + lambda * 1/2 * (c.x2 - c.x1)
      sh.col <- rep(sh.col, length.out=length(c.x1))
      sh.lty <- rep(sh.lty, length.out=length(c.x1))
      sh.lwd <- rep(sh.lwd, length.out=length(c.x1))
      lc.x <- lc.y <- numeric(length(c.x1))
      for (i in seq_len(length(c.x1))) {
         if (lambda[i] == 0) {
            if (!arrows_only) {
               segments(
                  c.x1[i],
                  c.y1[i],
                  c.x2[i],
                  c.y2[i],
                  lwd=sh.lwd[i],
                  col=sh.col[i],
                  lty=sh.lty[i])
            }
            phi <- atan2(
               y1[i] - y2[i],
               x1[i] - x2[i])
            r <- sqrt((x1[i] - x2[i])^2 + (y1[i] - y2[i])^2)
            lc.x[i] <- x2[i] + 2/3 * r * cos(phi)
            lc.y[i] <- y2[i] + 2/3 * r * sin(phi)
         } else {
            spl <- xspline(
               x=c(c.x1[i], spx[i], c.x2[i]),
               y=c(c.y1[i], spy[i], c.y2[i]),
               shape=1,
               draw=FALSE)
            if (!arrows_only) {
               lines(spl,
                  lwd=sh.lwd[i],
                  col=sh.col[i],
                  lty=sh.lty[i])
            }
            if (code %in% c(2, 3)) {
               x1[i] <- spl$x[3 * length(spl$x)/4]
               y1[i] <- spl$y[3 * length(spl$y)/4]
            }
            if (code %in% c(1, 3)) {
               x2[i] <- spl$x[length(spl$x)/4]
               y2[i] <- spl$y[length(spl$y)/4]
            }
            lc.x[i] <- spl$x[2/3 * length(spl$x)]
            lc.y[i] <- spl$y[2/3 * length(spl$y)]
         }
      }
   }
   if (code %in% c(2, 3)) {
      # head border outline
      theta <- atan2(
         (by2 - y1) * uin[2],
         (bx2 - x1) * uin[1])
      Rep <- rep(length(deg.arr), lx)
      p.x2 <- rep(bx2, Rep)
      p.y2 <- rep(by2, Rep)
      ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
      r.arr <- rep(r.arr, lx)
      if (open) {
         # head arrow with no color fill
         lines(
            (p.x2 + r.arr * cos(ttheta)/uin[1]),
            (p.y2 + r.arr * sin(ttheta)/uin[2]),
            lwd=h.lwd,
            col=h.col.bo,
            lty=h.lty)
      } else {
         # head arrow with color fill
         polygon(p.x2 + r.arr * cos(ttheta)/uin[1],
            p.y2 + r.arr * sin(ttheta)/uin[2],
            col=h.col,
            lwd=h.lwd,
            border=h.col.bo,
            lty=h.lty)
      }
   }
   if (code %in% c(1, 3)) {
      x1 <- bx1
      y1 <- by1
      tmp <- x1
      x1 <- x2
      x2 <- tmp
      tmp <- y1
      y1 <- y2
      y2 <- tmp
      theta <- atan2(
         (y2 - y1) * uin[2],
         (x2 - x1) * uin[1])
      lx <- length(x1)
      Rep <- rep(length(deg.arr), lx)
      p.x2 <- rep(x2, Rep)
      p.y2 <- rep(y2, Rep)
      ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
      r.arr <- rep(r.arr, lx)
      if (open) {
         lines(
            (p.x2 + r.arr * cos(ttheta)/uin[1]),
            (p.y2 + r.arr * sin(ttheta)/uin[2]),
            lwd=h.lwd,
            col=h.col.bo,
            lty=h.lty)
      } else {
         polygon(
            p.x2 + r.arr * cos(ttheta)/uin[1],
            p.y2 + r.arr * sin(ttheta)/uin[2],
            col=h.col,
            lwd=h.lwd,
            border=h.col.bo,
            lty=h.lty)
      }
   }
   list(lab.x=lc.x,
      lab.y=lc.y)
}

#' Get igraph arrow mode
#'
#' Get igraph arrow mode
#'
#' This function mimics the internal function `igraph:::i.get.arrow.mode()`.
#'
#' @family jam igraph internal functions
#'
get_igraph_arrow_mode <- function
(graph,
 arrow.mode=NULL)
{
   if (is.character(arrow.mode) &&
         length(arrow.mode) == 1 &&
         substr(arrow.mode, 1, 2) == "a:") {
      arrow.mode <- igraph::vertex_attr(graph,
         substring(arrow.mode, 3))
   }
   if (is.character(arrow.mode)) {
      tmp <- numeric(length(arrow.mode))
      tmp[arrow.mode %in% c("<", "<-")] <- 1
      tmp[arrow.mode %in% c(">", "->")] <- 2
      tmp[arrow.mode %in% c("<>", "<->")] <- 3
      arrow.mode <- tmp
   }
   if (length(arrow.mode) == 0) {
      if (igraph::is_directed(graph)) {
         arrow.mode <- 2
      } else {
         arrow.mode <- 0
      }
   }
   arrow.mode
}
