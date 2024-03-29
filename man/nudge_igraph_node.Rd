% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamenrich-cnet.R
\name{nudge_igraph_node}
\alias{nudge_igraph_node}
\title{Nudge igraph layout by node}
\usage{
nudge_igraph_node(
  g,
  nodes = NULL,
  x = 0,
  y = 0,
  nodes_xy = NULL,
  use_grep = FALSE,
  aspect = 1,
  debug = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{g}{\code{igraph} object that contains layout coordinates stored as
graph attribute \code{"layout"}, for example \code{graph_attr(g, "layout")}.}

\item{nodes}{\code{character} vector indicating which nodes in \code{g} should
be nudged.}

\item{x, y}{\code{numeric} values indicating the amount to move each node
defined by \code{nodes}. These values are relative to the x- and y-axis
ranges of the layout coordinates, and based upon \code{aspect} below.}

\item{nodes_xy}{\code{list} alternative to using arguments \verb{nodes,x,y}.
This argument assumes a \code{list} of \code{numeric} adjustments to \code{x} and \code{y},
where the \code{names(nodes_xy)} are node names. For example:
\itemize{
\item \verb{nodes_xy=list(APOE=c(x=0.05, y=-0.01), GAPDH=c(x=0.1, y=0.0)}
however the \code{numeric} vector does not need to contain names.
}}

\item{use_grep}{\code{logical} indicating whether to match values in \code{nodes}
to \code{V(g)$name} and \code{V(g)$label} using \code{jamba::provigrep()}, which
follows case-insensitive \code{grep()}. When \code{use_grep=FALSE} the
values in either \code{V(g)$name} or \code{V(g)$label} must be identical to
\code{nodes}.}

\item{aspect}{\code{numeric} indicating the aspect ratio of the output
plot. Any value other than \code{aspect=1} will use the observed x-axis
and y-axis range of the layout coordinates. When \code{aspect=1} then
the highest of x-axis and y-axis ranges is used for the relative
\code{x} and \code{y} adjustment. Note that \code{igraph::plot()} does not maintain
aspect ratio 1 by default, but \code{jam_igraph()} does maintain
aspect ratio 1 and is preferred.}

\item{debug}{\code{logical} indicating whether to plot the layout before
and after adjustment, drawing arrows to indicate the movement
of particular nodes. This plot is very basic, using base R \code{plot()},
and is only intended as a quick visual review.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}
}
\description{
Nudge igraph layout by node
}
\seealso{
Other jam cnet igraph functions: 
\code{\link{adjust_cnet_nodeset}()},
\code{\link{adjust_cnet_set_relayout_gene}()},
\code{\link{apply_cnet_direction}()},
\code{\link{apply_nodeset_spacing}()},
\code{\link{get_cnet_nodeset}()},
\code{\link{make_cnet_test}()},
\code{\link{rotate_igraph_layout}()}
}
\concept{jam cnet igraph functions}
\concept{jam igraph functions

This function takes an \code{igraph} object that contains a layout
stored in graph attributes, \code{graph_attr(g, "layout")}, and adjusts
the position of one or more nodes based upon relative x- and y-axis
size of the layout.}
