
#' Obtain or create layout for igraph object
#'
#' Obtain or create layout for igraph object
#'
#' This function is a simple helper function intended to retrieve
#' the node layout for an `igraph` object.
#' 
#' The layout is defined with the following priority:
#' 
#' 1. When `layout` is supplied as an argument, it is used.
#' When it is a `function` it is applied to `g` to produce
#' numeric `matrix`; otherwise it should be a numeric `matrix`
#' and is used as-is.
#' 2. When graph attribute 'layout' is defined, it is used
#' as described for argument `layout` above, accepting either
#' `function` or `matrix` values.
#' 3. When vertex attributes 'x' and 'y' are defined, optionally 'z',
#' their values are used to produce a numeric `matrix`.
#' 4. When `default_layout` is supplied as a `function` it is applied
#' to graph `g` to produce a numeric `matrix`.
#' 5. Finally, when `default_layout` is NULL, this function returns NULL.
#' This fallback is intended only when it is desirable not to apply
#' a new layout function, useful for large graphs.
#' 
#' ## Additional rules
#' 
#' * When `layout` is defined as a matrix with rownames, the rownames
#' are matched to vertex attribute 'name' if it exists,
#' using `igraph::V(g)$name`. This step is intended to help ensure nodes
#' the layout can be supplied in any order without regard to the
#' order defined in graph `g`.
#' 
#'    * When the `layout` rownames do not match vertex names, this function
#'    will `stop()`.
#' 
#' * When `layout` is defined as a `function`, or when any layout function
#' is applied as relevant, it is expected to return a `numeric` `matrix`,
#' or data which can be coerced to a `matrix` with `as.matrix()` or
#' `as(x, "matrix")`.
#' 
#'    * The matrix rownames are matched to vertex names as described above.
#'    * Note that `data.frame` rownames are only retained at this step when
#'    they were already explicitly defined before coersion to `matrix`.
#' 
#' ## Suggested Usage
#' 
#' * Define `layout` as a function in order to force the use
#' of that function to produce layout coordinates.
#' This step would always ignore pre-existing layout coordinates in
#' graph `g`.
#' * Define `layout` as NULL, and `default_layout` as a function,
#' to use an existing layout stored in graph `g`, then to apply the
#' default layout function only if no layout already existed in graph `g`.
#' * Define `layout` as NULL, and `default_layout` as NULL, to
#' return an existing layout stored in graph `g`, otherwise to return NULL
#' without applying any layout. This option would avoid computationally
#' expensive layout for large graphs, for example.
#' 
#' @family jam utility functions
#'
#' @returns `get_igraph_layout()` returns a `numeric` matrix when:
#'    the input graph `g` contains layout as a graph attribute as either
#'    numeric matrix or function,
#'    or coordinates as x,y,z vertex attributes.
#'    The layout will contain colnames that begin 'x', 'y', 'z',
#'    for consistency with downstream use.
#'    * When  there is no layout defined in `g` and
#'    `default_layout` is NULL, it returns NULL.  
#'    This logic is intended when it is preferable to avoid returning
#'    a layout if it does not already exist.
#'    * When `matrix` is returned, the number of rows
#'    matches the input graph `g` using `igraph::vcount(g)`,
#'    and in that order.
#'    All `rownames()` are defined to match vertex name when it exists,
#'    using `igraph::V(g)$name`.
#'
#' @param g `igraph` object
#' @param layout is always applied when not NULL, even when layout
#'    exists in `g`.  
#'    Input should be one of:
#'    * `numeric` matrix of layout coordinates, with `nrow(layout)`
#'    equal to the number of nodes `igraph::vcount(g)`.
#'    * `function` that takes `igraph` input, and returns `numeric` matrix
#'    of layout coordinates.
#'    * `NULL`: default, uses the graph attribute 'layout' if it exists,
#'    `igraph::graph_attr(g, "layout")` if it exists.
#'    If it does not exist, it follows `default_layout`.
#' @param default_layout is only applied when `layout` is NULL, and no
#'    layout is defined in graph `g`.  
#'    Input should be `function` or NULL:
#'    * Default `igraph::layout_nicely()` is used
#'    for consistency with igraph conventions.
#'    This algorithm should be safe for very large graphs,
#'    which may be extremely inefficient with `layout_with_qfr()` for example.
#'    * When NULL, it will return NULL unless layout is defined
#'    in the input `g` graph.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to any layout function called.
#' 
#' @examples
#' g <- make_cnet_test(2, c(12, 5))
#' gl <- get_igraph_layout(g, verbose=TRUE)
#' jam_igraph(g)
#' 
#' # apply repulse=4
#' gl2 <- get_igraph_layout(g, layout=layout_with_qfrf(repulse=4), verbose=TRUE)
#' jam_igraph(g, layout=gl2)
#' 
#' igraph::graph_attr(g, "layout") <- layout_with_qfrf(repulse=4)
#' gl3 <- get_igraph_layout(g, verbose=TRUE)
#' identical(gl2, gl3)
#' #> TRUE
#' 
#' g2 <- set_igraph_layout(g, spread_labels=TRUE)
#' jam_igraph(g2)
#' 
#' @export
get_igraph_layout <- function
(g,
 layout=NULL,
 default_layout=igraph::layout_nicely,
 verbose=FALSE,
 ...)
{
   #
   # obtain layout
   xy <- NULL;
   if (length(layout) == 0) {
      # igraph layout
      if ("layout" %in% igraph::graph_attr_names(g) &&
            length(igraph::graph_attr(g, "layout")) > 0) {
      	# graph attribute: layout
         if (verbose) {
            jamba::printDebug("get_igraph_layout(): ",
               "using graph_attr layout.");
         }
         xy <- igraph::graph_attr(g, "layout");
         #
         # if function, apply the function
         if (is.function(xy)) {
         	if (verbose) {
         		jamba::printDebug("get_igraph_layout(): ",
         			"applying graph_attr layout()");
         	}
         	# apply the layout function to produce numeric coordinates
         	xy <- xy(g, ...);
         }
      } else if (all(c("x", "y") %in% igraph::vertex_attr_names(g)) &&
      		is.numeric(igraph::vertex_attr(g, "x")) &&
      		is.numeric(igraph::vertex_attr(g, "y"))) {
      	# vertex attributes: x,y,z
      	xy <- cbind(
      		igraph::vertex_attr(g, "x"),
      		igraph::vertex_attr(g, "y"));
      	vattrs <- c("x", "y");
      	if ("z" %in% igraph::vertex_attr_names(g) &&
      			is.numeric(igraph::vertex_attr(g, "z"))) {
      		if (verbose) {
      			jamba::printDebug("get_igraph_layout(): ",
      				"using vertex_attr: ", c("x", "y", "z"));
      		}
      		vattrs <- c("x", "y", "z");
      		xy <- cbind(xy,
      			igraph::vertex_attr(g, "z"));
      	}
   		if (verbose) {
   			jamba::printDebug("get_igraph_layout(): ",
   				"using vertex_attr: ", vattrs);
   		}
      	#
      } else {
      	# apply default_layout
      	if (!is.function(default_layout)) {
      		if (verbose) {
      			jamba::printDebug("get_igraph_layout(): ",
      				"no graph layout, layout, default_layout, returning ",
      				"NULL");
      		}
      		# exit strategy, when NULL do not apply layout
      		return(NULL);
      	}
      	if (verbose) {
      		jamba::printDebug("get_igraph_layout(): ",
      			"applying default_layout()");
      	}
      	xy <- default_layout(g, ...);
      }
   } else if (is.function(layout)) {
      #
      if (verbose) {
         jamba::printDebug("get_igraph_layout(): ",
            "creating layout using ",
            "layout()");
      }
      xy <- layout(g, ...)
   } else {
   	if (verbose) {
   		jamba::printDebug("get_igraph_layout(): ",
   			"using layout as supplied.");
   	}
   	xy <- layout;
   }
   
   # coerce to matrix if possible
   if (!"matrix" %in% class(xy)) {
      # try two common methods to coerce to matrix
      xy <- tryCatch({
         as.matrix(xy);
      }, error=function(e){
         as(xy, "matrix");
      });
   }
   if (!"matrix" %in% class(xy)) {
   	stop_msg <- paste0("layout could not be coerced to matrix.")
   	stop(stop_msg)
   }

   # validate nrow versus input g
   if (!nrow(xy) == igraph::vcount(g)) {
      stop_msg <- paste0(
         "layout nrow (", nrow(xy),
         ") does not equal vcount (", igraph::vcount(g), ").");
      stop(stop_msg);
   }
   if (length(rownames(xy)) > 0) {
   	if ("name" %in% igraph::vertex_attr_names(g)) {
	      if (!all(rownames(xy) %in% igraph::V(g)$name)) {
	         stop("rownames(layout) do not match V(g)$name.")
	      }
   	}
   	# check for mis-ordered rows
   	if (!all(rownames(xy) == igraph::V(g)$name)) {
	   	nmatch <- match(igraph::V(g)$name, rownames(xy));
	   	if (any(is.na(nmatch))) {
	   		stop_msg <- paste0("Not all rownames(layout) are present in ",
	   			"V(g)$name vertex names.");
	   		stop(stop_msg);
	   	}
	   	if (verbose) {
	   		jamba::printDebug("get_igraph_layout(): ",
	   			"layout rows were re-ordered to match vertex names in V(g)$name.");
	   	}
	      xy <- xy[nmatch, , drop=FALSE];
   	}
   } else {
      rownames(xy) <- igraph::V(g)$name;
   }

   # confirm colnames begin with: 'x', 'y', 'z', then cycle the alphabet
   # - bit overkill, but what can you do
   use_colnames <- head(unique(c("x", "y", "z",
   	tolower(jamba::colNum2excelName(seq_len(ncol(xy)))))),
   	ncol(xy))
   colnames(xy) <- use_colnames;

   if (verbose > 1) {
      jamba::printDebug("get_igraph_layout(): ",
         "head(xy, 3):");
      print(head(xy, 3));
   }
   return(xy);
}

#' Set the node layout for an igraph object
#'
#' Set the node layout for an igraph object
#'
#' This function is a simple wrapper to `get_igraph_layout()` which
#' also defines the resulting layout in the graph `g`.
#' 
#' @returns `set_igraph_layout()` returns `igraph` object with
#'    layout defined per function arguments.
#' 
#' @rdname get_igraph_layout
#' @param prefer `character` vector with preferred method of storage:
#'    * 'graph_attr': store the `matrix` in graph attribute 'layout'.
#'    * 'vertex_attr': store coordinates in vertex attributes 'x', 'y',
#'    and optionally 'z' when defined.
#' @param spread_labels `logical` default FALSE, whether to call
#'    `spread_igraph_labels()` to re-position node labels radially
#'    away from incoming edges.
#'    * Note that `spread_igraph_labels()` also by default `do_reorder=TRUE`
#'    which will re-order nodes by color, border, label, and name.
#'    * It is used to use TRUE when node labels were previously spread,
#'    so that the angle of offset is updated per the new layout coordinates.
#' 
#' @export
set_igraph_layout <- function
(g,
 layout=NULL,
 default_layout=igraph::layout_nicely,
 verbose=FALSE,
 prefer=c("graph_attr"),
 spread_labels=FALSE,
 ...)
{
	#
	use_layout <- get_igraph_layout(g=g,
		layout=layout,
		default_layout=default_layout,
		verbose=verbose,
		...)
	
	# Todo: un-set the layout when NULL?
	# Todo: Consider applying spread_igraph_labels()
	# Todo: Consider storing previous layout in versioned state?
	
	# Define into the igraph object
	if ("graph_attr" %in% prefer) {
		#
		igraph::graph_attr(g, "layout") <- use_layout;
	}
	
	# spread labels
	if (isTRUE(spread_labels)) {
		g <- spread_igraph_labels(g=g,
			force_layout=FALSE,
			...)
	}
	
	return(g);
}
