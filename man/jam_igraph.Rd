% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam-igraph-plot.R
\name{jam_igraph}
\alias{jam_igraph}
\title{Jam wrapper to plot igraph}
\usage{
jam_igraph(
  x,
  ...,
  xlim = c(-1, 1),
  ylim = c(-1, 1),
  expand = 0.03,
  rescale = FALSE,
  node_factor = 1,
  node_factor_l = NULL,
  edge_factor = 1,
  edge_factor_l = NULL,
  label_factor = 1,
  label_factor_l = NULL,
  label_dist_factor = 1,
  label_dist_factor_l = 1,
  use_shadowText = FALSE,
  plot_function = jam_plot_igraph,
  edge_bundling = c("none", "connections", "nodegroups"),
  nodegroups = NULL,
  render_nodes = TRUE,
  render_edges = TRUE,
  render_nodelabels = TRUE,
  render_groups = TRUE,
  vectorized_node_shapes = TRUE,
  verbose = FALSE,
  debug = NULL
)
}
\arguments{
\item{x}{\code{igraph} object to be plotted}

\item{...}{additional arguments are passed to \code{igraph::plot.igraph()}}

\item{xlim, ylim}{default x and y axis limits}

\item{expand}{numeric value used to expand the x and y axis ranges,
where \code{0.03} expands each size \verb{3\%}.}

\item{rescale}{logical indicating whether to rescale the layout
coordinates to \code{c(-1, 1)}. When \code{rescale=FALSE} the original
layout coordinates are used as-is without change.}

\item{node_factor}{numeric value multiplied by \code{V(x)$size} to adjust
the relative size of all nodes by a common numeric scalar value.}

\item{node_factor_l, label_factor_l, label_dist_factor_l}{\code{list}
of vectors, where the names of the \code{list} are attribute
names, and the names of each vector are attributes values.
The vector values are used as scalar multipliers, analogous to
\code{node_factor}. The purpose is to apply scalar values to different
subsets of nodes. For example, consider:
\verb{node_factor_l=list(nodeType=c(Gene=1, Set=2)}. The list name
\code{"nodeType"} says to look at \code{vertex_attr(x, "nodeType")}. Nodes
where \code{nodeType="Gene"} will use \code{1}, and where \code{nodeType="Set"}
will use \code{2} as the scalar value.}

\item{edge_factor}{numeric value multiplied by \code{E(x)$width} to adjust
the relative width of all edges by a common numeric scalar value.}

\item{label_factor}{numeric value multiplied by \code{V(x)$label.cex}
and \code{E(x)$label.cex} to adjust the relative size of all labels on
nodes and edges by a common numeric scalar value.}

\item{label_dist_factor}{numeric value multiplied by \code{V(x)$label.dist}
to adjust the relative distance of all nodes labels from the node center
by a common numeric scalar value.}
}
\description{
Jam wrapper to plot igraph
}
\details{
This function is a lightweight wrapper around \code{igraph::plot.igraph()}
intended to handle \code{rescale=FALSE} properly, which is not done
in the former function (as of January 2020). The desired outcome
is for the \code{xlim} and \code{ylim} defaults to be scaled according to the
\code{igraph} layout. Similarly, the \code{vertex.size} and \code{vertex.label.dist}
parameters should also be scaled proportionally.

You can use argument \code{label_factor_l=list(nodeType=c(Gene=0.01, Set=1))}
to hide labels for Gene nodes, and display labels for Set nodes.
Note that due to a quirk in \code{igraph}, setting \code{label.cex=0}
will revert the font to default size, and will not hide the
label.
}
\examples{
## example showing how to use the list form
## This form resizes nodes where V(g)$nodeType \%in\% "Gene" by 2x,
## and resizes nodes where V(g)$nodeType \%in\% "Set" by 3x.
node_factor_l <- list(nodeType=c(Gene=2, Set=3));

## This form multiplies label.dist for nodeType="Gene" nodes by 1,
## and multiplies label.dist for nodeType="Set" nodes by 0.5
label_dist_factor_l <- list(nodeType=c(Gene=1, Set=0.5))

# jam_igraph(g, node_factor_l=node_factor_l, label_dist_factor_l=label_dist_factor_l);

# Example using edge bundling by community detection
g <- igraph::make_graph("Zachary");
gcom <- igraph::cluster_leading_eigen(g);

jam_igraph(g,
   layout=layout_with_qfr,
   edge_bundling="nodegroups",
   nodegroups=gcom,
   vertex.color=colorjam::group2colors(gcom$membership))

cfuncs <- list(cluster_leading_eigen=igraph::cluster_leading_eigen,
   cluster_edge_betweenness=igraph::cluster_edge_betweenness,
   cluster_fast_greedy=igraph::cluster_fast_greedy,
   cluster_spinglass=igraph::cluster_spinglass)
opar <- par("mfrow"=c(2, 2));
for (i in seq_along(cfuncs)) {
   cfunc <- cfuncs[[i]];
   gcom <- cfunc(g);
   jam_igraph(g,
      layout=layout_with_qfr,
      edge_bundling="nodegroups",
      nodegroups=gcom,
      mark.groups=gcom,
      mark.expand=60,
      vertex.color=colorjam::group2colors(gcom$membership))
   title(main=names(cfuncs)[i]);
}
par(opar);

# fancy example showing mark.groups and colorizing
# edges using node colors
gcom <- igraph::cluster_spinglass(g);
igraph::V(g)$color <- colorjam::group2colors(gcom$membership);
g <- color_edges_by_nodes(g);
jam_igraph(g,
   layout=layout_with_qfr,
   edge_bundling="nodegroups",
   nodegroups=gcom,
   mark.groups=gcom)

# same but adjust midpoint of edge bundles
jam_igraph(g,
   layout=layout_with_qfr,
   edge_bundling="nodegroups",
   nodegroups=gcom,
   mark.groups=gcom,
   midpoint=c(0.4, 0.6),
   detail=14)

# same but using node connections
jam_igraph(g,
   layout=layout_with_qfr,
   edge_bundling="connections",
   nodegroups=gcom,
   mark.groups=gcom)

}
\seealso{
Other jam igraph functions: 
\code{\link{cnet2df}()},
\code{\link{cnet2im}()},
\code{\link{cnetplotJam}()},
\code{\link{cnetplot_internalJam}()},
\code{\link{color_edges_by_nodes_deprecated}()},
\code{\link{color_edges_by_nodes}()},
\code{\link{drawEllipse}()},
\code{\link{edge_bundle_nodegroups}()},
\code{\link{enrichMapJam}()},
\code{\link{fixSetLabels}()},
\code{\link{get_bipartite_nodeset}()},
\code{\link{igraph2pieGraph}()},
\code{\link{jam_plot_igraph}()},
\code{\link{layout_with_qfrf}()},
\code{\link{layout_with_qfr}()},
\code{\link{memIM2cnet}()},
\code{\link{mem_multienrichplot}()},
\code{\link{rectifyPiegraph}()},
\code{\link{relayout_with_qfr}()},
\code{\link{removeIgraphBlanks}()},
\code{\link{removeIgraphSinglets}()},
\code{\link{reorderIgraphNodes}()},
\code{\link{rotate_igraph_layout}()},
\code{\link{spread_igraph_labels}()},
\code{\link{subgraph_jam}()},
\code{\link{subsetCnetIgraph}()},
\code{\link{subset_igraph_components}()},
\code{\link{with_qfr}()}

Other jam plot functions: 
\code{\link{grid_with_title}()},
\code{\link{mem_enrichment_heatmap}()},
\code{\link{mem_gene_path_heatmap}()},
\code{\link{mem_legend}()},
\code{\link{mem_multienrichplot}()},
\code{\link{mem_plot_folio}()}
}
\concept{jam igraph functions}
\concept{jam plot functions}
