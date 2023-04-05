
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
 categoryShape=c("pie",
    "coloredrectangle",
    "circle",
    "ellipse"),
 geneShape=c("pie",
    "coloredrectangle",
    "circle",
    "ellipse"),
 categoryColor="#E5C494",
 geneColor="#B3B3B3",
 categoryLabelColor="darkblue",
 geneLabelColor="grey25",
 categorySize=6,
 geneSize=3,
 categoryCex=0.8,
 geneCex=0.5,
 # frame_darkFactor=1.4,
 frame_darkFactor=NULL,
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
 hide_solo_pie=TRUE,
 remove_blanks=TRUE,
 spread_labels=FALSE,
 repulse=3.5,
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

   # Confirm rownames and colnames do not contain duplicate values.
   # If so, make them unique.
   make_unique_im_dims <- function
   (memIM,
    warn=TRUE)
   {
      #
      if (length(rownames(memIM)) == 0) {
         rownames(memIM) <- paste0("row", seq_len(nrow(memIM)));
      }
      if (length(colnames(memIM)) == 0) {
         colnames(memIM) <- paste0("set", seq_len(ncol(memIM)));
      }
      rcnames <- c(rownames(memIM),
         colnames(memIM));
      if (any(duplicated(rcnames))) {
         # duplicated rownames
         if (any(duplicated(rownames(memIM)))) {
            rownames(memIM) <- jamba::makeNames(rownames(memIM),
               suffix="_r");
         }
         # duplicated colnames
         if (any(duplicated(rownames(memIM)))) {
            colnames(memIM) <- jamba::makeNames(colnames(memIM),
               suffix="_c");
         }
         # duplicated across rownames and colnames
         rcnames <- c(rownames(memIM),
            colnames(memIM));
         if (any(duplicated(rcnames))) {
            new_rc <- jamba::makeNames(rcnames,
               suffix="_v");
            rownames(memIM) <- head(rcnames,
               nrow(memIM));
            colnames(memIM) <- tail(rcnames,
               -nrow(memIM));
         }
         # very rare change the versioned rownames/colnames are duplicated
         rcnames <- c(rownames(memIM),
            colnames(memIM));
         if (any(duplicated(rcnames))) {
            memIM <- make_unique_im_dims(memIM,
               warn=FALSE)
         }
         # optionally print a warning
         if (TRUE %in% warn) {
            warning(paste0("Note from mem2cnet(): ",
               "The input incidence matrix contained some",
               " duplicated rownames/colnames. They were",
               " made unique."))
         }
      }
      return(memIM);
   }
   memIMdimnames <- dimnames(memIM);
   memIM <- make_unique_im_dims(memIM,
      warn=TRUE);
   if (!identical(memIMdimnames, dimnames(memIM))) {
      if (length(geneIM) > 0) {
         gmatch <- match(rownames(geneIM),
            memIMdimnames[[1]]);
         rownames(geneIM) <- memIMdimnames[[1]][gmatch]
      }
      if (length(geneIMcolors) > 0) {
         gmatch <- match(rownames(geneIMcolors),
            memIMdimnames[[1]]);
         rownames(geneIMcolors) <- memIMdimnames[[1]][gmatch]
      }
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
   # nodeType
   isset <- igraph::V(g)$name %in% colnames(memIM);
   igraph::V(g)$nodeType <- ifelse(isset, "Set", "Gene");
   igraph::V(g)$size <- ifelse(isset, categorySize, geneSize);
   igraph::V(g)$label.cex <- ifelse(isset, categoryCex, geneCex);
   igraph::V(g)$label.color <- ifelse(isset,
      categoryLabelColor,
      geneLabelColor);
   igraph::V(g)$color <- ifelse(isset, categoryColor, geneColor);
   igraph::V(g)$frame.lwd <- 1;
   if (length(frame_darkFactor) > 0 && frame_darkFactor != 1) {
      igraph::V(g)$frame.color <- jamba::makeColorDarker(
         igraph::V(g)$color,
         darkFactor=frame_darkFactor,
         ...);
   } else {
      igraph::V(g)$frame.color <- default_igraph_values()$vertex$frame.color;
   }

   ## Optionally apply gene node coloring
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "applying node colors");
   }
   if (length(geneIM) > 0) {
      if (length(geneIMcolors) == 0) {
         if (length(colorV) < ncol(geneIM)) {
            colorV <- colorjam::rainbowJam(ncol(geneIM),
               Crange=c(70, 110))
            names(colorV) <- colnames(geneIM);
         }
         if (!all(names(colorV) %in% colnames(geneIM))) {
            names(colorV) <- colnames(geneIM);
         } else {
            colorV <- colorV[colnames(geneIM)];
         }
         # make gene color incidence matrix
         geneIMcolors <- do.call(rbind, lapply(rownames(geneIM), function(i){
            j <- geneIM[match(i, rownames(geneIM)),];
            ifelse(j > 0,
               colorV,
               "#FFFFFF")
         }))
         rownames(geneIMcolors) <- rownames(geneIM);
         colnames(geneIMcolors) <- colnames(geneIM);

      }
      geneIM <- subset(geneIM, rownames(geneIM) %in% igraph::V(g)$name[!isset]);
      geneIMcolors <- subset(geneIMcolors, rownames(geneIMcolors) %in% igraph::V(g)$name[!isset]);

      # # enrichIM
      # if (length(enrichIM) == 0) {
      #    # derive from geneIM
      #    lapply(colnames(geneIM), function(i){
      #       j <- rownames(geneIM)[abs(geneIM[,i]) > 0]
      #    })
      # }
      gene_match <- match(rownames(geneIM),
         igraph::V(g)$name[!isset]);
      gene_which <- which(!isset)[gene_match];
      vseq <- seq_len(igraph::vcount(g));

      # iterate node attributes
      i_node_attributes <- c(
         "pie",
         "pie.value",
         "pie.color",
         "pie.names",
         "coloredrect.color",
         "coloredrect.nrow",
         "coloredrect.byrow",
         "coloredrect.value",
         "coloredrect.names")

      for (j in i_node_attributes) {
         if (!j %in% igraph::list.vertex.attributes(g)) {
            if (verbose) {
               jamba::printDebug("mem2cnet(): ",
                  "Adding vertex attribute '", j, "'");
            }
            if (jamba::igrepHas("[.]color$", j)) {
               igraph::vertex_attr(g, j) <- as.list(igraph::V(g)$color);
            } else if (jamba::igrepHas("byrow", j)) {
               igraph::vertex_attr(g, j) <- rep(unlist(coloredrect_byrow),
                  length.out=igraph::vcount(g));
            } else if (jamba::igrepHas("[.]names$", j)) {
               igraph::vertex_attr(g, j) <- rep(j,
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
      # propagate gene-related attributes
      igraph::V(g)$pie[gene_which] <- lapply(gene_which, function(i){
         rep(1, ncol(geneIM));
      });
      igraph::V(g)$pie.value[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         geneIM[match(i, rownames(geneIM)),,drop=FALSE];
      });
      igraph::V(g)$pie.color[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         jamba::nameVector(geneIMcolors[match(i, rownames(geneIMcolors)),],
            colnames(geneIMcolors));
      });
      igraph::V(g)$pie.names[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         colnames(geneIM);
      });
      ## Now do the same for coloredrectangle node shapes
      igraph::V(g)$coloredrect.color[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         jamba::nameVector(geneIMcolors[match(i, rownames(geneIMcolors)),],
            colnames(geneIMcolors));
      });
      igraph::V(g)$coloredrect.value[gene_which] <- lapply(jamba::nameVector(rownames(geneIM)), function(i){
         geneIM[match(i, rownames(geneIM)),,drop=FALSE];
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
      igraph::V(g)$shape <- ifelse(isset, categoryShape, geneShape);
   }

   # Optionally apply category/set node coloring
   if (verbose) {
      jamba::printDebug("memIM2cnet(): ",
         "applying category/set node colors");
   }
   if (length(enrichIM) > 0 && length(enrichIMcolors) > 0) {
      enrichIM <- subset(enrichIM, rownames(enrichIM) %in% igraph::V(g)$name[isset]);
      enrich_match <- match(rownames(enrichIM),
         igraph::V(g)$name[isset]);
      enrich_which <- which(isset)[enrich_match];

      # iterate set-relevant attributes
      i_set_attributes <- c("pie",
         "pie.value",
         #"pie.color",
         "coloredrect.color",
         "coloredrect.nrow",
         "coloredrect.byrow",
         "coloredrect.value")

      for (j in i_set_attributes) {
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
         hide_solo_pie=hide_solo_pie,
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
         hide_solo_pie=hide_solo_pie,
         ...);
   }

   # optionally orient labels radially away from incoming edges
   if (spread_labels) {
      if (verbose) {
         jamba::printDebug("memIM2cnet(): ",
            "applying spread_igraph_labels()");
      }
      g <- spread_igraph_labels(g,
         do_reorder=TRUE,
         # y_bias=y_bias,
         repulse=repulse,
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
