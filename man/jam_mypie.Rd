% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamenrich-igraphshapes.R
\name{jam_mypie}
\alias{jam_mypie}
\title{Vectorized mypie() function for igraph vertex pie polygons}
\usage{
jam_mypie(
  x,
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
  ...
)
}
\description{
Vectorized mypie() function for igraph vertex pie polygons
}
\details{
This function is a light rewrite of \code{igraph:::mypie()}, except
that this function determines polygon coordinates without
drawing them, instead returns the polygon coordinates to
the calling function \code{shape.jampie.plot()} which in turn
draws all polygons once using the vectorized approach
described for \code{graphics::polygon()}.

One small additional change, pie shapes with only one large
100\% wedge no longer display the small line from origin.
}
\seealso{
Other jam igraph shapes: 
\code{\link{shape.coloredrectangle.plot}()},
\code{\link{shape.ellipse.plot}()},
\code{\link{shape.jampie.plot}()}
}
\concept{jam igraph shapes}