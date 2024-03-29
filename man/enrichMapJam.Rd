% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamenrich-base.r
\name{enrichMapJam}
\alias{enrichMapJam}
\title{Create enrichMap igraph object}
\usage{
enrichMapJam(
  x,
  n = 50,
  vertex.label.font = 1,
  vertex.label.cex = 1,
  keyColname = "ID",
  nodeLabel = c("Name", "Description", "ID"),
  descriptionColname = "Description",
  nodeLabelFunc = NULL,
  overlapThreshold = 0.2,
  msigdbGmtT = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{either \code{enrichResult} or \code{data.frame} containing
enrichment results, specifically expecting colnames to
contain one of \code{c("ID","Description","Name")}
to represent the node name, and \code{c("Description")} to represent
the description, if present.}

\item{n}{numeric value indicating the maximum number of nodes to
include in the final network.}

\item{vertex.label.font, vertex.label.cex}{attributes to define the
default node label font and size.}

\item{keyColname, nodeLabel, descriptionColname}{character vectors
indicating the colname to use for the node name and label.}

\item{nodeLabelFunc}{optional function to apply to \code{V(g)$name} in
order to create \code{V(g)$label}. One suggestion is \code{fixSetLabels()}
which applies word wrap, and optional max character length.}

\item{overlapThreshold}{numeric value indicating the minimum
Jaccard overlap, where edges with lower values are deleted from
the \code{igraph} object.}

\item{...}{additional arguments are passed to \code{enrichDF2enrichResult()}
when the input \code{x} is a \code{data.frame}.}
}
\value{
\code{igraph} object, whose nodes represent each enriched pathway,
and are sized based upon the number of genes involved in the
enrichment, and are colored based upon the \code{log10(Pvalue)}
using \code{colorjam::vals2colorLevels()}, a function that applies
a color gradient to a numeric range.
Each edge has attributes: \code{overlap} containing Jaccard overlap,
\code{overlap_count} with the number of genes in common between
the two nodes, and \code{overlap_max_pct} with the maximum percent
overlap between two nodes (overlap count)/(smaller node size).
}
\description{
Create enrichMap igraph object from enrichResult.
}
\details{
This function could also be called \code{enrichResult2emap()}.

This function is a minor extension to the original function
DOSE::enrichMap() which is now rewritten in the source package
to \code{enrichplot::emapplot()}. The major differences:
\itemize{
\item This function returns an \code{igraph} object, which can be manipulated
using network-related functions.
\item This function calculates overlap using \code{dist(...,method="binary")}
which is a much faster method for calculating the Jaccard overlap.
\item This function also calculates the overlap count, another helpful
measure for filtering network connections, for example to remove
links with only one gene, even if they overlap is above the
required threshold. Many spurious network connections are removed
with this filter, and it appears to be a helpful option.
}
}
\seealso{
Other jam conversion functions: 
\code{\link{cnet2df}()},
\code{\link{cnet2im}()},
\code{\link{enrichDF2enrichResult}()},
\code{\link{enrichList2IM}()},
\code{\link{enrichList2df}()}

Other jam igraph functions: 
\code{\link{cnet2df}()},
\code{\link{cnet2im}()},
\code{\link{cnetplotJam}()},
\code{\link{cnetplot_internalJam}()},
\code{\link{color_edges_by_nodegroups}()},
\code{\link{color_edges_by_nodes_deprecated}()},
\code{\link{color_edges_by_nodes}()},
\code{\link{color_nodes_by_nodegroups}()},
\code{\link{communities2nodegroups}()},
\code{\link{drawEllipse}()},
\code{\link{edge_bundle_bipartite}()},
\code{\link{edge_bundle_nodegroups}()},
\code{\link{fixSetLabels}()},
\code{\link{flip_edges}()},
\code{\link{get_bipartite_nodeset}()},
\code{\link{igraph2pieGraph}()},
\code{\link{jam_igraph}()},
\code{\link{jam_plot_igraph}()},
\code{\link{label_communities}()},
\code{\link{layout_with_qfrf}()},
\code{\link{layout_with_qfr}()},
\code{\link{mem2emap}()},
\code{\link{memIM2cnet}()},
\code{\link{mem_multienrichplot}()},
\code{\link{nodegroups2communities}()},
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
\code{\link{sync_igraph_communities}()},
\code{\link{with_qfr}()}
}
\concept{jam conversion functions}
\concept{jam igraph functions}
