
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
#' @returns `igraph` object with Concept network data, containing
#'    pathways connected to genes. Each node has attribute `"nodeType"`
#'    of either `"Set"` or `"Gene"`.
#'
#' @param memIM one of the following:
#'    * `Mem` S4 object, preferred. In this case, the arguments regarding
#'    'geneIM' and 'enrichIM' data are taken directly from 'mem'.
#'    * legacy `list` mem object, for backward compatibility,
#'    * `numeric` matrix in the form of a `memIM` gene-pathway incidence
#'    matrix. In this case, other arguments involving `geneIM` and `enrichIM`
#'    matrices are required for correct behavior of this function.
#'    When `mem` format is supplied, relevant arguments which are empty will
#'    use corresponding data from `mem`, for example `geneIM`, `geneIMcolors`,
#'    `enrichIM`, `enrichIMcolors`.
#' @param categoryShape,geneShape `character` string with node shape,
#'    default 'pie' uses pie nodes.
#'    Note that 'pie' shapes with only one segment are converted to
#'    'circle' for convenience.
#'    When using `jam_igraph()` to plot, it already treats single-segment
#'    pie nodes as circle anyway, and vectorizes pie node plotting to
#'    make the overall rendering time substantially faster.
#' @param forceColors `logical` default FALSE, whether to force node colors
#'    to use `categoryColor` and `geneColor` arguments, otherwise colors
#'    are defined using `enrichIMcolors()` and `geneIMcolors()`,
#'    respectively.
#'    Enabling this option will skip all pie and directional styling.
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
#'    `matrix` data used only when the input 'mem' is not an `Mem` or legacy
#'    `list` mem object which already contains these data.
#'    When input 'mem' is supplied as a `matrix`, these value enable
#'    this function to operate on nearly any custom data.
#'    * `geneIMcolors` is used to define gene (row) node colors
#'    * `enrichIMcolors` is used to define set (column) node colors
#'    * `geneIMdirection`,`enrichIMdirection` is used to define optional
#'    border colors defined by the direction, where -1 is down,
#'    0 is no change, and +1 is up.
#'    Use `direction_col_fn` to define a custom color function, however the
#'    default uses the reversed "RdBu" Brewer color ramp with blue (down),
#'    white (no change) and red (up).
#' @param coloredrect_nrow,coloredrect_ncol,coloredrect_byrow arguments
#'    used when `geneShape="coloredrectangle"`, to define layout and
#'    placement of colors across columns in `geneIMcolors`. By default,
#'    one row of colors is used.
#' @param colorV `character` optional vector of R colors, taken from 'mem'
#'    for `Mem` or legacy `list` mem objects.
#'    When supplied, it should be a vector named to match
#'    `colnames(enrichIM)`.
#'    When defined, these colors override `enrichIMcolors`.
#' @param direction_col_fn `function` used to colorize 'Set' node borders
#'    via enrichIMdirection. The default uses `colorjam::col_div_xf()`
#'    which applies reverse Brewer "RdBu" for blue (down), white (no change),
#'    and red (up) with maximum color at `2.0`.
#'    When not supplied, `direction_cutoff` and `direction_max` are used
#'    with `colorjam::col_div_xf().`
#' @param direction_cutoff `numeric` default is taken from `memIM` when
#'    provided as `Mem`, typical default is `1.0` to require directional
#'    z-score at least `1.0` in order to apply any directional border
#'    color.
#' @param direction_max `numeric` default 2, the numeric value to apply
#'    maximum color from the color gradient function.
#' @param gene_direction_col_fn `function` used to colorize 'Gene' node borders
#'    via enrichIMdirection. The default uses `colorjam::col_div_xf()`
#'    which applies reverse Brewer "RdBu" for blue (down), white (no change),
#'    and red (up) with maximum color at `1.2`.
#'    *When not supplied, `gene_direction_cutoff` and `gene_direction_max`
#'    are used with `colorjam::col_div_xf().`
#' @param gene_direction_cutoff `numeric` default 0, the minimum value
#'    in geneIMdirection to apply up/down color. Using 0 will colorize
#'    any non-zero color.
#' @param gene_direction_max `numeric` default 1.2, the numeric value to apply
#'    maximum color from the color gradient function. Typical values
#'    are `-1` and `1`, so using `1.2` will apply nearly the maximum color.
#' @param hide_solo_pie `logical` default TRUE, passed to
#'    `apply_cnet_direction()` to determine whether to display border
#'    only as one outer frame color when all colors are identical.
#'    When FALSE, all pie wedges are individually colored.
#' @param remove_blanks `logical` default TRUE, whether to remove blank
#'    color subsections from each node, using `removeIgraphBlanks()`.
#'    This argument is useful for 'pie', 'jampie', or 'coloredrectangle'
#'    node shapes, so they will only indicate the relevant color.
#' @param remove_singlet_genes `logical` default TRUE, whether to remove
#'    singlet genes, which are genes (rows) not represented in any pathway
#'    gene sets (columns).
#' @param do_reorder `logical` default TRUE used to re-order nodes
#'    in "equivalent positions" by border color, color, etc. by calling
#'    `reorder_igraph_nodes()`. Equivalent nodes are common in Cnet
#'    (bipartite) networks, for example numerous 'Gene' nodes may
#'    be connected to the same 2 pathway 'Set' nodes, and therefore
#'    the node positions are interchangeable. Sorting by node attributes
#'    (notably color, then alphabetically by label) helps visual
#'    review.
#' @param spread_labels `logical` default FALSE, whether to spread node
#'    labels away from incoming edges, also adding label distance via
#'    vertex attribute 'label.dist'.
#'    This step calls `spread_igraph_labels()`, and can be customized
#'    further by passing arguments through '...' ellipses.
#' @param repulse `numeric` value passed to `layout_with_qfr()`
#'    when either `do_reorder` or `spread_labels` is TRUE. Otherwise
#'    there is no layout applied.
#'    The 'repulse' value effectively defines node spacing, usually
#'    with values between 3 and 4, which 3 giving broad spacing,
#'    and 4 giving very close spacing of nodes in equivalent network
#'    positions, and broad spacing otherwise. Higher values tend
#'    to "clump" nodes closer together.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to downstream functions:
#'    * when `remove_blanks=TRUE` it is passed to `removeIgraphBlanks()`
#'    * when `spread_labels=TRUE` it is passed to `spread_igraph_labels()`
#'
#' @family jam Mem utilities
#' @family jam igraph functions
#' 
#' @examples
#' use_sets <- c("eNOS Signaling",
#'    "Growth Hormone Signaling",
#'    "mTOR Signaling")
#' jam_igraph(mem2cnet(Memtest[, use_sets, ]),
#'    use_shadowText=TRUE)
#' 
#' @export
mem2cnet <- function
(memIM,
 categoryShape=c("pie",
    "coloredrectangle",
    "circle",
    "ellipse"),
 geneShape=c("pie",
    "coloredrectangle",
    "circle",
    "ellipse"),
 forceColors=FALSE,
 categoryColor="#E5C494",
 geneColor="#B3B3B3",
 categoryLabelColor="darkblue",
 geneLabelColor="grey25",
 categorySize=12,
 geneSize=6,
 categoryCex=1,
 geneCex=0.8,
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
 direction_col_fn=NULL,
 direction_cutoff=NULL,
 direction_max=2,
 gene_direction_cutoff=0,
 gene_direction_max=1.2,
 gene_direction_col_fn=NULL,
 hide_solo_pie=TRUE,
 remove_blanks=TRUE,
 remove_singlet_genes=TRUE,
 do_reorder=TRUE,
 spread_labels=FALSE,
 repulse=3.5,
 verbose=FALSE,
 ...)
{
   categoryShape <- match.arg(categoryShape);
   geneShape <- match.arg(geneShape);
   geneLabelColor <- head(geneLabelColor, 1);
   categoryLabelColor <- head(categoryLabelColor, 1);

   if (is.list(memIM) &&
         "memIM" %in% names(memIM)) {
      # convert legacy mem to Mem
      memIM <- tryCatch({
         list_to_Mem(memIM)
      }, error=function(e) {
         memIM
      })
   }
   if (inherits(memIM, "Mem")) {
      Mem <- memIM;
      memIM <- memIM(Mem);
      geneIM <- geneIM(Mem);
      geneIMcolors <- geneIMcolors(Mem);
      geneIMdirection <- geneIMdirection(Mem);
      enrichIM <- enrichIM(Mem);
      enrichIMcolors <- enrichIMcolors(Mem);
      enrichIMdirection <- enrichIMdirection(Mem);
      colorV <- colorV(Mem);
      if (length(direction_cutoff) == 0) {
         direction_cutoff <- thresholds(Mem)$direction_cutoff;
      }
   } else if (is.list(memIM) &&
         "memIM" %in% names(memIM)) {
      ## Accept memIM as list mem output from multiEnrichMap()
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
   
   if (TRUE %in% remove_singlet_genes) {
      gene_cts <- rowSums(memIM != 0);
      if (any(gene_cts == 0)) {
         keep_genes <- rownames(memIM)[gene_cts > 0]
         if (length(keep_genes) == 0) {
            stop_msg <- paste0("No genes are present in any pathways.");
            stop(stop_msg);
         }
         memIM <- memIM[keep_genes, , drop=FALSE]
         if (length(geneIM) > 0) {
            geneIM <- geneIM[keep_genes, , drop=FALSE]
         }
         if (length(geneIMcolors) > 0) {
            geneIMcolors <- geneIMcolors[keep_genes, , drop=FALSE]
         }
         if (length(geneIMdirection) > 0) {
            geneIMdirection <- geneIMdirection[keep_genes, , drop=FALSE]
         }
      }
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
      jamba::printDebug("mem2cnet(): ",
         "igraph::graph_from_biadjacency_matrix()");
   }
   g <- igraph::graph_from_biadjacency_matrix(memIM != 0);

   ## Adjust aesthetics
   if (verbose) {
      jamba::printDebug("mem2cnet(): ",
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
      jamba::printDebug("mem2cnet(): ",
         "applying node colors");
   }
   if (!isTRUE(forceColors) && length(geneIM) > 0) {
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
         if (!j %in% igraph::vertex_attr_names(g)) {
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
      geneIMrn <- jamba::nameVector(rownames(geneIM));
      igraph::V(g)$pie.value[gene_which] <- lapply(geneIMrn, function(i){
         geneIM[match(i, rownames(geneIM)),,drop=FALSE];
      });
      igraph::V(g)$pie.color[gene_which] <- lapply(geneIMrn, function(i){
         jamba::nameVector(geneIMcolors[match(i, rownames(geneIMcolors)),],
            colnames(geneIMcolors));
      });
      igraph::V(g)$pie.names[gene_which] <- lapply(geneIMrn, function(i){
         colnames(geneIM);
      });
      ## Now do the same for coloredrectangle node shapes
      igraph::V(g)$coloredrect.color[gene_which] <- lapply(geneIMrn,
         function(i){
         jamba::nameVector(geneIMcolors[match(i, rownames(geneIMcolors)),],
            colnames(geneIMcolors));
      });
      igraph::V(g)$coloredrect.value[gene_which] <- lapply(geneIMrn,
         function(i){
         geneIM[match(i, rownames(geneIM)),,drop=FALSE];
      });
      if (length(coloredrect_byrow) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_byrow,
            length.out=length(gene_which));
      }
      if (length(coloredrect_nrow) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_nrow,
            length.out=length(gene_which));
      }
      if (length(coloredrect_ncol) > 0) {
         igraph::V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_ncol,
            length.out=length(gene_which));
      }
      igraph::V(g)$coloredrect.names[gene_which] <- lapply(geneIMrn,
         function(i){
         colnames(geneIM);
      });
      igraph::V(g)$shape <- ifelse(isset, categoryShape, geneShape);
   }

   # Optionally apply category/set node coloring
   if (verbose) {
      jamba::printDebug("mem2cnet(): ",
         "applying category/set node colors");
   }
   if (!isTRUE(forceColors) &&
         length(enrichIM) > 0 &&
         length(enrichIMcolors) > 0) {
      enrichIM <- subset(enrichIM,
         rownames(enrichIM) %in% igraph::V(g)$name[isset]);
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
         if (!j %in% igraph::vertex_attr_names(g)) {
            if (jamba::igrepHas("color$", j)) {
               igraph::vertex_attr(g, j) <- as.list(igraph::V(g)$color);
            } else if (jamba::igrepHas("byrow", j)) {
               igraph::vertex_attr(g, j) <- as.list(rep(coloredrect_byrow,
                  length.out=igraph::vcount(g)));
            } else if (jamba::igrepHas("ncol", j) &&
                  length(coloredrect_ncol) > 0) {
               igraph::vertex_attr(g, j) <- as.list(rep(coloredrect_ncol,
                  length.out=igraph::vcount(g)));
            } else if (jamba::igrepHas("nrow", j) &&
                  length(coloredrect_nrow) > 0) {
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
      enrichIMrn <- jamba::nameVector(rownames(enrichIM));
      igraph::V(g)$pie.value[enrich_which] <- lapply(enrichIMrn, function(i){
         enrichIM[i,,drop=FALSE];
      });
      igraph::V(g)$pie.color[enrich_which] <- lapply(enrichIMrn, function(i){
         jamba::nameVector(enrichIMcolors[i,],
            colnames(enrichIMcolors));
      });
      ## Now do the same for coloredrectangle node shapes
      igraph::V(g)$coloredrect.color[enrich_which] <- lapply(enrichIMrn,
         function(i){
         jamba::nameVector(enrichIMcolors[i,],
            colnames(enrichIMcolors));
      });
      igraph::V(g)$coloredrect.value[enrich_which] <- lapply(enrichIMrn,
         function(i){
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
         jamba::printDebug("mem2cnet(): ",
            "applying removeIgraphBlanks()");
      }
      g <- removeIgraphBlanks(g, ...);
   }

   # Freshen pie.color,coloredrect.color by using colorV by name
   if (!isTRUE(forceColors) && length(colorV) > 0) {
      if (verbose) {
         jamba::printDebug("mem2cnet(): ",
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
   if (!isTRUE(forceColors) &&
         length(geneIMdirection) > 0 &&
         !all(geneIMdirection %in% c(0, NA))) {
      if (verbose) {
         jamba::printDebug("mem2cnet(): ",
            "applying geneIMdirection");
      }
      if (length(gene_direction_col_fn) == 0) {
         if (length(gene_direction_max) == 0 ||
               any(gene_direction_max <= 0)) {
            gene_direction_max <- 1.2;
         }
         if (length(gene_direction_cutoff) > 0 &&
               any(gene_direction_cutoff >= gene_direction_max)) {
            gene_direction_max <- max(gene_direction_cutoff, na.rm=TRUE);
            gene_direction_cutoff <- gene_direction_max - 
               (gene_direction_max / 1e5);
         }
         gene_direction_col_fn <- colorjam::col_div_xf(gene_direction_max,
            floor=gene_direction_cutoff)
      }
      g <- apply_cnet_direction(cnet=g,
         hitim=geneIMdirection,
         col=gene_direction_col_fn,
         hide_solo_pie=hide_solo_pie,
         ...);
   }
   if (!isTRUE(forceColors) &&
         length(enrichIMdirection) > 0 &&
         !all(enrichIMdirection %in% c(0, NA))) {
      if (verbose) {
         jamba::printDebug("mem2cnet(): ",
            "applying enrichIMdirection");
      }
      if (length(direction_col_fn) == 0) {
         if (length(direction_max) == 0 ||
               any(direction_max <= 0)) {
            direction_max <- 2;
         }
         if (length(direction_cutoff) > 0 &&
               any(direction_cutoff >= direction_max)) {
            direction_max <- max(direction_cutoff, na.rm=TRUE);
            direction_cutoff <- direction_max - (direction_max / 1e5);
         }
         direction_col_fn <- colorjam::col_div_xf(direction_max,
            floor=direction_cutoff)
      }
      g <- apply_cnet_direction(cnet=g,
         hitim=enrichIMdirection,
         col=direction_col_fn,
         hide_solo_pie=hide_solo_pie,
         ...);
   }

   # optionally orient labels radially away from incoming edges
   if (spread_labels) {
      if (verbose) {
         jamba::printDebug("mem2cnet(): ",
            "applying spread_igraph_labels()");
      }
      g <- spread_igraph_labels(g,
         do_reorder=do_reorder,
         # y_bias=y_bias,
         repulse=repulse,
         ...);
   } else if (isTRUE(do_reorder)) {
      if (!"layout" %in% igraph::graph_attr_names(g)) {
         g <- relayout_with_qfr(g,
            spread_labels=FALSE,
            repulse=repulse,
            ...);
      }
      g <- reorder_igraph_nodes(g,
         colorV=colorV,
         ...);
   }

   if (verbose) {
      jamba::printDebug("mem2cnet(): ",
         "complete");
   }

   return(g);
}

#' @rdname mem2cnet
#' @export
memIM2cnet <- mem2cnet
