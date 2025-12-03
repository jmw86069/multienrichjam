
# alphahull functions


#' Alpha shape calculation
#' 
#' Alpha shape calculation for a set of points, and alpha threshold
#' 
#' This function is primarily intended to be called by `make_point_hull()`,
#' since that function also iterates `alpha` values until it finds a
#' suitable, and successful, threshold.
#' 
#' @family jam utility functions
#' 
#' @returns `ashape` object, which is a `list` containing:
#' * edges: x,y coordinates of Delauney triangulation of the alpha-shape.
#' * length: length of the alpha-shape.
#' * alpha: value of alpha used.
#' * alpha.extremes: `integer` index of points which were alpha-extremes.
#' * delvor.obj: `delvor` object with Delauney/Voronoi supporting data.
#' * x: x,y coordinates of input data
#' 
#' @param x,y `numeric` vector with coordinate points.
#' @param alpha `numeric` with alpha threshold to use.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' n <- 300
#' theta <- runif(n, 0, 2*pi)
#' r <- sqrt(runif(n, 0.25^2, 0.5^2))
#' x <- cbind(0.5+r*cos(theta), 0.5+r*sin(theta))
#' alpha <- 0.1
#' ashape.obj <- ashape(x, alpha=alpha)
#' 
#' plot(ashape.obj$x, asp=1)
#' segments(x0=ashape.obj$edges[, "x1"], x1=ashape.obj$edges[, "x2"],
#'    y0=ashape.obj$edges[, "y1"], y1=ashape.obj$edges[, "y2"])
#' 
#' @export
ashape <- function
(x,
 y=NULL,
 alpha,
 ...)
{
   if (alpha < 0) {
      stop("Parameter alpha must be greater or equal to zero")
   }
   # run Delauney triangulation Voronoi diagram as needed
   if (!inherits(x, "delvor")) {
      dd.obj <- delvor(x, y)
   } else {
      dd.obj <- x
   }
   
   # extract mesh
   xy.data <- dd.obj$x;
   mesh <- dd.obj$mesh;
   
   dm1 <- sqrt(
      (mesh[, "x1"] - mesh[, "mx1"])^2 +
      (mesh[, "y1"] - mesh[, "my1"])^2);
   dm2 <- sqrt(
      (mesh[, "x1"] - mesh[, "mx2"])^2 +
      (mesh[, "y1"] - mesh[, "my2"])^2);
   dm1[mesh[, "bp1"] == 1] <- Inf;
   dm2[mesh[, "bp2"] == 1] <- Inf;
   n <- nrow(xy.data);
   ind <- seq_len(n);
   ind.on <- grDevices::chull(xy.data)
   ind.in <- ind[-ind.on]
   n.on <- length(ind.on)
   n.in <- length(ind.in)
   
   if (n.in > 0) {
      aux <- rbind(
         cbind(mesh[, c("ind1", "ind2"), drop=FALSE], dm1),
         cbind(mesh[, c("ind1", "ind2"), drop=FALSE], dm2),
         cbind(mesh[, c("ind2", "ind1"), drop=FALSE], dm1),
         cbind(mesh[, c("ind2", "ind1"), drop=FALSE], dm2))
      # calculate by factor
      # - useful to keep the original order of values?
      fc <- factor(aux[, 1])
      aux <- tapply(aux[, 3], fc, max)
      alpha.max <- cbind(aux, as.numeric(levels(fc)));
      use_match <- match(
         alpha.max[alpha < alpha.max[, 1], 2],
         ind.in)
      alpha.ext <- c(ind.on,
         ind.in[na.omit(use_match)]);
   } else {
      alpha.ext <- ind.on
   }
   
   n.edges <- nrow(mesh);
   is.edge <- numeric(0);
   ind.is <- 0;
   i1 <- match(mesh[, 1], alpha.ext);
   i2 <- match(mesh[, 2], alpha.ext);
   is.edge <- which(i1 & i2);
   aux <- subset(mesh, (i1 & i2) %in% TRUE);
   # aux <- mesh[is.edge, , drop=FALSE];
   
   # fiddle with triangle vertices ordering
   n.pos <- nrow(aux);
   pm.x <- (aux[, "x1"] + aux[, "x2"]) * 0.5;
   pm.y <- (aux[, "y1"] + aux[, "y2"]) * 0.5;
   dm <- sqrt(
      (aux[, "x1"] - aux[, "x2"])^2 +
      (aux[, "y1"] - aux[, "y2"])^2) * 0.5;
   betw <- rep(NA, n.pos)
   for (i in 1:n.pos) {
      if (aux[i, "mx1"] == aux[i, "mx2"]) {
         test_rank_y <- rank(c(
            aux[i, "my1"],
            aux[i, "my2"],
            pm.y[i]));
         if (test_rank_y[3] == 2) {
            betw[i] <- 1
         }
      } else {
         test_rank_x <- rank(c(
            aux[i, "mx1"],
            aux[i, "mx2"],
            pm.x[i]));
         if (test_rank_x[3] == 2) {
            betw[i] <- 1
         }
      }
   }
   
   # rowMin and rowMax
   l_matrix <- cbind(dm1[is.edge], dm2[is.edge], dm * betw);
   l.min <- apply(l_matrix, 1, min, na.rm=TRUE)
   l.max <- apply(l_matrix, 1, max, na.rm=TRUE)
   
   # attempt to simplify (?)
   in.ashape <- (l.min <= alpha & alpha <= l.max);
   edges <- subset(aux, in.ashape);
   # edges <- matrix(t(aux[in.ashape, ]),
   #    byrow TRUE,
   #    ncol=12)
   # colnames(edges) <- colnames(aux)
   
   # return object
   ashape.obj <- list(edges=edges,
      length=sum(2 * dm[in.ashape]),
      alpha=alpha,
      alpha.extremes=alpha.ext,
      delvor.obj=dd.obj,
      x=xy.data);
   class(ashape.obj) <- "ashape";
   invisible(ashape.obj)
}


#' alphahull delvor port Delauney triangulation and Voronoi diagram
#' 
#' @keywords internal
#' @noRd
delvor <- function
(x,
 y=NULL,
 ...) 
{
   X <- xy.coords(x, y)
   x <- cbind(X$x, X$y)
   if (dim(x)[1] <= 2) {
      stop("At least three non-collinear points are required")
   }
   # convert to Delauney triangulation
   # - consider jitter=TRUE
   tri.obj <- interp::tri.mesh(X);
   # extract triangle vertices
   tri <- interp::triangles(tri.obj);
   
   nt <- nrow(tri)
   circenter <- matrix(nrow=nt, ncol=2) 
   colnames(circenter) <- c("circumx", "circumy")
   for (i in 1:nt){
      aux <- interp::circum(
         x=c(x[tri[i, 1], 1],
            x[tri[i, 2], 1],
            x[tri[i, 3], 1]),
         y=c(x[tri[i, 1], 2],
            x[tri[i, 2], 2],
            x[tri[i, 3], 2]))
      circenter[i, ] <- c(aux$x, aux$y)
   }
   tri.info <- cbind(tri, circenter)
   
   n.tri <- nrow(tri.info)
   n.arc <- max(tri.info[, 7:9])
   
   # simplified previous logic
   aux1 <- cbind(
      tri.info[, c("arc1", "node2", "node3"), drop=FALSE],
      seq_len(n.tri),
      tri.info[, "tr1"])
   aux2 <- cbind(
      tri.info[, c("arc2", "node1", "node3"), drop=FALSE],
      seq_len(n.tri),
      tri.info[, "tr2"])
   aux3 <- cbind(
      tri.info[, c("arc3", "node1", "node2"), drop=FALSE],
      seq_len(n.tri),
      tri.info[, "tr3"])
   # combine into tall form
   aux <- rbind(aux1, aux2, aux3)
   # jamba::printDebug("aux:");print(aux);# debug
   
   # remove duplicates
   aux <- subset(aux, !duplicated(aux[, 1]))
   # repeated <- duplicated(aux[, 1])
   # aux <- aux[!repeated, ]
   
   # assign colnames
   colnames(aux) <- c("arc", "ind1", "ind2", "indm1", "indm2")
   bp1 <- (aux[, "indm1"] == 0)
   bp2 <- (aux[, "indm2"] == 0)
   is.dummy <- which(bp2)
   n.dummy <- length(is.dummy)
   circumcentres <- tri.info[, c("circumx", "circumy"), drop=FALSE]
   
   away <- max(diff(range(x[, 1])),
      diff(range(x[, 2])))
   for (i in is.dummy) {
      n.tri <- n.tri + 1;
      dum <- dummycoor(tri.obj,
         x[aux[i, "ind1"], ],
         x[aux[i, "ind2"], ],
         tri.info[aux[i, "indm1"], c("circumx", "circumy")],
         away)
      circumcentres <- rbind(circumcentres, dum);
      aux[i, "indm2"] <- n.tri;
   }
   
   # create mesh
   mesh <- cbind(
      aux[, c("ind1", "ind2")],
      x[aux[, "ind1"], ],
      x[aux[, "ind2"], ],
      circumcentres[aux[, "indm1"], ],
      circumcentres[aux[, "indm2"], ],
      bp1,
      bp2)
   colnames(mesh) <- c("ind1", "ind2", "x1", "y1", "x2", "y2", 
      "mx1", "my1", "mx2", "my2", "bp1", "bp2")
   delvor.obj <- list(mesh=mesh,
      x=x,
      tri.obj=tri.obj)
   class(delvor.obj) <- "delvor";
   invisible(delvor.obj)
}


#' alphahull dummycoor port
#' 
#' @keywords internal
#' @noRd
dummycoor <- function
(tri.obj,
 l1,
 l2,
 m,
 away)
{
   v <- l2 - l1;
   v <- c(v[2], -v[1]);
   norm <- sum(v^2);
   if (norm > 0) {
      v <- v/norm;
   }
   mp <- (l1 + l2)/2;
   eps <- 1e-05;
   test <- mp + eps * v;
   inconv <- interp::in.convex.hull(tri.obj=tri.obj,
      x=test[1],
      y=test[2]);
   
   if (inconv) {
      dum <- m - away * v;
   } else {
      dum <- m + away * v;
   }
   return(dum)
}
