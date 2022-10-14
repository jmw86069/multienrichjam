
#' Convert MultiEnrichment incidence matrix to Cnet plot
#'
#' Convert MultiEnrichment incidence matrix to Cnet plot
#'
#' This function takes `list` output from `multiEnrichMap()`
#' and uses three elements to create a Cnet plot:
#'
#' 1. `"memIM"` gene-pathway incidence matrix
#' 2. `"geneIM"` the gene incidence matrix
#' 3. `"enrichIM"` the pathway enrichment matrix
#'
#' @return `igraph` object with Concept network data, containing
#'    pathways connected to genes. Each node has attribute `"nodeType"`
#'    of either `"Set"` or `"Gene"`.
#'
#' @param memIM `numeric` matrix, or `mem` output in `list` format. When
#'    `mem` format is supplied, relevant arguments which are empty will
#'    use corresponding data from `mem`, for example `geneIM`, `geneIMcolors`,
#'    `enrichIM`, `enrichIMcolors`.
#' @param categoryShape,geneShape `character` string indicating the
#'    preferred default node shape. In general, `pie` shapes with one segment
#'    are converted to `circle` for simplicity during plot functions,
#'    although when using `jam_igraph()` to plot, it will treat single-segment
#'    pie nodes as circle anyway, and vectorizes pie node plotting to
#'    make the overall rendering time substantially faster.
#' @param categoryColor,geneColor `character` R color for default node colors,
#'    used when `geneIMcolors`,`enrichIMcolors` is not supplied, respectively.
#' @param categoryLabelColor,geneLabelColor `character` R color used as
#'    default node label color.
#' @param categorySize,geneSize `numeric` default node size.
#' @param categoryCex,geneCex `numeric` adjustment to default node label
#'    font size.
#' @param frame_darkFactor `numeric` passed to `jamba::makeColorDarker()`
#'    so the frame color is slightly darker than the node fill color.
#' @param geneIM,geneIMcolors,geneIMdirection,enrichIM,enrichIMcolors,enrichIMdirection
#'    `matrix` objects typically associated with `mem` output from
#'    `multiEnrichMap()`, however they are optional so this function
#'    can be applied broadly to any incidence matrix.
#'    * `geneIMcolors` is used to define gene (row) node colors
#'    * `enrichIMcolors` is used to define set (column) node colors
#'    * `geneIMdirection`,`enrichIMdirection` is used to define optional
#'    border colors defined by the direction, where -1 is down, 0 is no change, and
#'    +1 is up. Use `direction_col` to define a custom color function.
#' @param coloredrect_nrow,coloredrect_ncol,coloredrect_byrow arguments
#'    used when `geneShape="coloredrectangle"`, to define layout and
#'    placement of colors across columns in `geneIMcolors`. By default,
#'    one row of colors is used.
#' @param colorV `character` optional vector of R colors, named by enrichment
#'    names that appear in `colnames(enrichIM)` when supplied. When
#'    defined, these colors override those defined in `enrichIMcolors`.
#' @param remove_blanks `logical` indicating whether to remove blank
#'    subsections from each node, by calling `removeIgraphBlanks()`.
#' @param spread_labels `logical` indicating whether to call
#'    `spread_igraph_labels()` to orient labels radially away from incoming
#'    edges.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to downstream functions:
#'    * when `remove_blanks=TRUE` it is passed to `removeIgraphBlanks()`
#'    * when `spread_labels=TRUE` it is passed to `spread_igraph_labels()`
#'
#' @family jam igraph functions
#'
#' @export
memIM2cnet <- function
(memIM,
 categoryShape=c("pie","coloredrectangle", "circle", "ellipse"),
 geneShape=c("pie","coloredrectangle", "circle", "ellipse"),
 categoryColor="#E5C494",
 geneColor="#B3B3B3",
 categoryLabelColor=c("blue3"),
 geneLabelColor="grey40",
 categorySize=5,
 geneSize=2.5,
 categoryCex=0.9,
 geneCex=0.7,
 frame_darkFactor=1.4,
 geneIM=NULL,
 geneIMcolors=NULL,
 geneIMdirection=NULL,
 enrichIM=NULL,
 enrichIMcolors=NULL,
 enrichIMdirection=NULL,
 coloredrect_nrow=1,
 coloredrect_ncol=NULL,
 coloredrect_byrow=TRUE,
 colorV=NULL,
 direction_col=colorjam::col_div_xf(1.2),
 remove_blanks=TRUE,
 spread_labels=FALSE,
 verbose=FALSE,
 ...)
{
   categoryShape <- match.arg(categoryShape);
   geneShape <- match.arg(geneShape);
   geneLabelColor <- head(geneLabelColor, 1);
   categoryLabelColor <- head(categoryLabelColor, 1);

   ## Accept memIM as list mem output from multiEnrichMap()
   if (is.list(memIM) &&
         "memIM" %in% names(memIM)) {
      if (length(geneIM) == 0) {
         geneIM <- memIM[["geneIM"]];
      }
      if (length(geneIMcolors) == 0) {
         geneIMcolors <- memIM[["geneIMcolors"]];
      }
      if (length(geneIMdirection) == 0) {
         geneIMdirection <- memIM[["geneIMdirection"]];
      }
      if (length(enrichIM) == 0) {
         enrichIM <- memIM[["enrichIM"]];
      }
      if (length(enrichIMcolors) == 0) {
         enrichIMcolors <- memIM[["enrichIMcolors"]];
      }
      if (length(enrichIMdirection) == 0) {
         enrichIMdirection <- memIM[["enrichIMdirection"]];
      }
      if (length(colorV) == 0) {
         colorV <- memIM[["colorV"]];
      }
      memIM <- memIM[["memIM"]];
   }

   ## Convert to igraph
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "igraph::graph_from_incidence_matrix()");
   }
   g <- igraph::graph_from_incidence_matrix(memIM != 0);

   ## Adjust aesthetics
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "basic aesthetics");
   }
   igraph::V(g)$nodeType <- "Set";
   igraph::V(g)$nodeType[match(rownames(memIM), igraph::V(g)$name)] <- "Gene";
   # table(igraph::V(g)$nodeType);
   isset <- igraph::V(g)$nodeType %in% "Set";
   igraph::V(g)$size <- ifelse(isset, categorySize, geneSize);
   igraph::V(g)$label.cex <- ifelse(isset, categoryCex, geneCex);
   igraph::V(g)$label.color <- ifelse(isset, categoryLabelColor, geneLabelColor);
   igraph::V(g)$color <- ifelse(isset, categoryColor, geneColor);
   if (length(frame_darkFactor) > 0 && frame_darkFactor != 1) {
      igraph::V(g)$frame.color <- jamba::makeColorDarker(
         igraph::V(g)$color,
         darkFactor=frame_darkFactor,
         ...);
   }

   ## Optionally apply gene node coloring
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "applying node colors");
   }
   if (length(geneIM) > 0 && length(geneIMcolors) > 0) {
      geneIM <- subset(geneIM, rownames(geneIM) %in% igraph::V(g)$name[!isset]);
      gene_match <- match(rownames(geneIM),
         igraph::V(g)$name[!isset]);
      gene_which <- which(!isset)[gene_match];
      for (j in c("pie", "pie.value", "pie.color", "pie.names",
         "coloredrect.color", "coloredrect.nrow",
         "coloredrect.byrow", "coloredrect.value",
         "coloredrect.names")) {
         if (!j %in% igraph::list.vertex.attributes(g)) {
            if (jamba::igrepHas("color$", j)) {
               igraph::vertex_attr(g, j) <- as.list(igraph::V(g)$color);
            } else if (jamba::igrepHas("byrow", j)) {
               igraph::vertex_attr(g, j) <- rep(unlist(coloredrect_byrow),
                  length.out=igraph::vcount(g));
            } else if (jamba::igrepHas("ncol", j) && length(coloredrect_ncol) > 0) {
               igraph::vertex_attr(g, j) <- rep(unlist(coloredrect_ncol),
                  length.out=igraph::vcount(g));
            } else if (jamba::igrepHas("nrow", j) && length(coloredrect_nrow) > 0) {
               igraph::vertex_attr(g, j) <- rep(unlist(coloredrect_nrow),
                  length.out=igraph::vcount(g));
            } else {
               igraph::vertex_attr(g, j) <- as.list(rep(1, igraph::vcount(g)));
            }
         }
      }
      igraph::V(g)$pie[gene_which] <- lapply(gene_which, function(i){
         rep(1, ncol(geneIM));
      });
      igraph::V(g)$pie.value[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         geneIM[i,,drop=FALSE];
      });
      igraph::V(g)$pie.color[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         jamba::nameVector(geneIMcolors[i,],
            colnames(geneIMcolors));
      });
      igraph::V(g)$pie.names[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         colnames(geneIM);
      });
      ## Now do the same for coloredrectangle node shapes
      igraph::V(g)$coloredrect.color[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         jamba::nameVector(geneIMcolors[i,],
            colnames(geneIMcolors));
      });
      igraph::V(g)$coloredrect.value[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         geneIM[i,,drop=FALSE];
      });
      if (length(coloredrect_byrow) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_byrow, length.out=length(gene_which));
      }
      if (length(coloredrect_nrow) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_nrow, length.out=length(gene_which));
      }
      if (length(coloredrect_ncol) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_ncol, length.out=length(gene_which));
      }
      igraph::V(g)$coloredrect.names[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         colnames(geneIM);
      });
      igraph::V(g)$shape[gene_which] <- geneShape;
   }
   ## Optionally apply category/set node coloring
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "applying category/set node colors");
   }
   if (length(enrichIM) > 0 && length(enrichIMcolors) > 0) {
      enrichIM <- subset(enrichIM, rownames(enrichIM) %in% igraph::V(g)$name[isset]);
      enrich_match <- match(rownames(enrichIM),
         igraph::V(g)$name[isset]);
      enrich_which <- which(isset)[enrich_match];
      for (j in c("pie", "pie.value", #"pie.color",
         "coloredrect.color", "coloredrect.nrow",
         "coloredrect.byrow", "coloredrect.value")) {
         if (!j %in% list.vertex.attributes(g)) {
            if (jamba::igrepHas("color$", j)) {
               igraph::vertex_attr(g, j) <- as.list(igraph::V(g)$color);
            } else if (jamba::igrepHas("byrow", j)) {
               igraph::vertex_attr(g, j) <- as.list(rep(coloredrect_byrow,
                  length.out=igraph::vcount(g)));
            } else if (jamba::igrepHas("ncol", j) && length(coloredrect_ncol) > 0) {
               igraph::vertex_attr(g, j) <- as.list(rep(coloredrect_ncol,
                  length.out=igraph::vcount(g)));
            } else if (jamba::igrepHas("nrow", j) && length(coloredrect_nrow) > 0) {
               igraph::vertex_attr(g, j) <- as.list(rep(coloredrect_nrow,
                  length.out=igraph::vcount(g)));
            } else {
               igraph::vertex_attr(g, j) <- as.list(rep(1, igraph::vcount(g)));
            }
         }
      }
      igraph::V(g)$pie[enrich_which] <- lapply(enrich_which, function(i){
         rep(1, ncol(enrichIM));
      });
      igraph::V(g)$pie.value[enrich_which] <- lapply(jamba::nameVector(rownames(enrichIM)), function(i){
         enrichIM[i,,drop=FALSE];
      });
      igraph::V(g)$pie.color[enrich_which] <- lapply(jamba::nameVector(rownames(enrichIM)), function(i){
         jamba::nameVector(enrichIMcolors[i,],
            colnames(enrichIMcolors));
      });
      ## Now do the same for coloredrectangle node shapes
      igraph::V(g)$coloredrect.color[enrich_which] <- lapply(jamba::nameVector(rownames(enrichIM)), function(i){
         jamba::nameVector(enrichIMcolors[i,],
            colnames(enrichIMcolors));
      });
      igraph::V(g)$coloredrect.value[enrich_which] <- lapply(jamba::nameVector(rownames(enrichIM)), function(i){
         enrichIM[i,,drop=FALSE];
      });
      if (length(coloredrect_byrow) > 0) {
         igraph::V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_byrow,
            length.out=length(enrich_which));
      }
      if (length(coloredrect_nrow) > 0) {
         igraph::V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_nrow,
            length.out=length(enrich_which));
      }
      if (length(coloredrect_ncol) > 0) {
         igraph::V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_ncol,
            length.out=length(enrich_which));
      }
      igraph::V(g)$shape[enrich_which] <- categoryShape;
   }
   # optionally remove igraph blanks
   if (remove_blanks) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying removeIgraphBlanks()");
      }
      g <- removeIgraphBlanks(g, ...);
   }

   # Freshen pie.color,coloredrect.color by using colorV by name
   if (length(colorV) > 0) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying colorV");
      }
      for (attr_name in c("pie.color", "coloredrect.color")) {
         igraph::vertex_attr(g, attr_name) <- lapply(igraph::vertex_attr(g, attr_name), function(i){
            j <- ifelse(names(i) %in% names(colorV) & !isColorBlank(i),
               colorV[names(i)],
               i);
         });
      }
      # also update node color when pie has only one color, so they are in sync
      pie_singlets <- which(lengths(igraph::vertex_attr(g, "pie.color")) == 1)
      if (length(pie_singlets) > 0) {
         igraph::vertex_attr(g, "color", index=pie_singlets) <- unname(unlist(
            igraph::vertex_attr(g, "pie.color", index=pie_singlets)))
      }
   }

   # optionally apply direction as a border color
   if (length(geneIMdirection) > 0 && !all(geneIMdirection %in% c(0, NA))) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying geneIMdirection");
      }
      g <- apply_cnet_direction(cnet=g,
         hitim=geneIMdirection,
         col=direction_col,
         ...);
   }
   if (length(enrichIMdirection) > 0 && !all(enrichIMdirection %in% c(0, NA))) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying enrichIMdirection");
      }
      g <- apply_cnet_direction(cnet=g,
         hitim=enrichIMdirection,
         col=direction_col,
         ...);
   }

   # optionally orient labels radially away from incoming edges
   if (spread_labels) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying spread_igraph_labels()");
      }
      g <- spread_igraph_labels(g,
         do_reorder=FALSE,
         # y_bias=y_bias,
         # repulse=repulse,
         ...);
   }

   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "complete");
   }

   return(g);
}

#' @rdname memIM2cnet
#' @export
mem2cnet <- memIM2cnet
