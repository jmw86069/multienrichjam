% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamenrich-utils.R
\name{enrichList2geneHitList}
\alias{enrichList2geneHitList}
\title{Extract gene hit list from list of enrichResult}
\usage{
enrichList2geneHitList(
  enrichList,
  geneColname,
  geneDelim = "[,/ ]",
  make_unique = TRUE,
  verbose = FALSE,
  ...
)
}
\value{
\code{list} of character vectors, containing the unique
set of genes involved in each enrichment.
}
\description{
Extract gene hit list from list of enrichResult
}
\details{
This function is mainly for internal use in multienrichjam,
it takes a list of \code{enrichResult} objects, and determines
the full set of genes involved in each \code{enrichResult}.

This function also works with \code{ComplexHeatmap::HeatmapList}
objects.
}
\seealso{
Other jam utility functions: 
\code{\link{avg_angles}()},
\code{\link{avg_colors_by_list}()},
\code{\link{call_fn_ellipsis_deprecated}()},
\code{\link{cell_fun_bivariate}()},
\code{\link{collapse_mem_clusters}()},
\code{\link{colorRamp2D}()},
\code{\link{deconcat_df2}()},
\code{\link{display_colorRamp2D}()},
\code{\link{filter_mem_genes}()},
\code{\link{filter_mem_sets}()},
\code{\link{find_colname}()},
\code{\link{get_igraph_layout}()},
\code{\link{gsubs_remove}()},
\code{\link{handle_igraph_param_list}()},
\code{\link{isColorBlank}()},
\code{\link{make_legend_bivariate}()},
\code{\link{make_point_hull}()},
\code{\link{mem_find_overlap}()},
\code{\link{order_colors}()},
\code{\link{rank_mem_clusters}()},
\code{\link{rotate_coordinates}()},
\code{\link{subgraph_jam}()},
\code{\link{subset_mem}()},
\code{\link{summarize_node_spacing}()},
\code{\link{xyAngle}()}
}
\concept{jam utility functions}
