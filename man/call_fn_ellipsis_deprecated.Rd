% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam-call-fn-ellipsis.R
\name{call_fn_ellipsis_deprecated}
\alias{call_fn_ellipsis_deprecated}
\title{Call function using safe ellipsis arguments}
\usage{
call_fn_ellipsis_deprecated(FUN, ...)
}
\arguments{
\item{FUN}{\code{function} that should be called with arguments in \code{...}}

\item{...}{arguments are passed to \code{FUN()} in safe manner.}
}
\value{
output from \code{FUN()} when called with relevant named arguments
from ellipsis \code{...}
}
\description{
Call function using safe ellipsis arguments
}
\details{
This function is deprecated, instead please use \code{jamba::call_fn_ellipsis()}.

This function is a wrapper function intended to help
pass ellipsis arguments \code{...} from a parent function
to an external function in a safe way. It will only
include arguments from \code{...} that are recognized by
the external function.

When the external function \code{FUN} arguments \code{formals()} includes
ellipsis \code{...}, then the \code{...} will be passed as-is without
change.

When the external function \code{FUN} arguments \code{formals()} does not
include ellipsis \code{...}, then only named arguments in \code{...} that
are recognized by \code{FUN} will be passed, as defined by
\code{names(formals(FUN))}.

Note that arguments must be named.
}
\examples{
new_mean <- function(x, trim=0, na.rm=FALSE) {
   mean(x, trim=trim, na.rm=na.rm)
}
x <- c(1, 3, 5, NA);
new_mean(x, na.rm=TRUE);
# throws an error as expected (below)
# new_mean(x, na.rm=TRUE, color="red");

call_fn_ellipsis_deprecated(new_mean, x=x, na.rm=TRUE, color="red")
# throws an error as expected (below)
# call_fn_ellipsis_deprecated(new_mean, x=x, color="red")

}
\seealso{
Other jam utility functions: 
\code{\link{avg_angles}()},
\code{\link{avg_colors_by_list}()},
\code{\link{cell_fun_bivariate}()},
\code{\link{collapse_mem_clusters}()},
\code{\link{colorRamp2D}()},
\code{\link{deconcat_df2}()},
\code{\link{display_colorRamp2D}()},
\code{\link{enrichList2geneHitList}()},
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
