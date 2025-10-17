# jamenrich-igraph.r
# functions related to igraph manipulations

#' igraph layout using qgraph Fruchterman-Reingold
#'
#' igraph layout using qgraph Fruchterman-Reingold
#'
#' This function provides Fruchterman-Reingold layout for an igraph
#' object using the implementation from the qgraph package, which provides
#' important configuration options deprecated in the igraph implementation.
#' Notably, the repulse.rad parameter is helpful in adjusting the relative
#' spacing of vertices, where higher values cause tighter packing of
#' vertices, and lower values allows greater spacing between vertices.
#'
#' @param g igraph object
#' @param repulse exponent power used to scale the radius effect around
#'    each vertex. The default is slightly higher than the cube of
#'    the number of vertices, but as the number of vertices increases,
#'    values from 3.5 to 4 and higher are more effective for layout.
#' @param area The area of the plot, default is the square of the number of
#'    vertices times 8. Changes to plot area will also affect values of
#'    \code{repulse} and \code{repulse.rad}.
#' @param repulse.rad Repulse radius, defaults to the the number of
#'    vertices raised to the \code{repulse} power.
#' @param constraints optional two-column matrix with the coordinates of
#'    nodes which should not be modified, and NA values for nodes where
#'    the position can be modified.
#' @param constrain `character` optional vector of node names that should
#'    be constrained. This argument is a convenient shortcut for defining
#'    `constraints`, which is a layout coordinate matrix with `NA` values
#'    on each row where the coordinate is free to move, and `numeric`
#'    values where the coordinate is fixed. For graph `g` that contains layout
#'    in `igraph::graph_attr(g, "layout")`, the `init` can be defined with
#'    this layout, then `constraints` is defined using `constrain`.
#' @param ... other arguments are sent to
#'    `qgraph::qgraph.layout.fruchtermanreingold()`
#'
#' @return two-column numeric matrix with coordinates for each vertex.
#'
#' @seealso `qgraph::qgraph.layout.fruchtermanreingold()`
#'
#' @family jam igraph functions
#'
#' @examples
#' if (suppressPackageStartupMessages(require(igraph))) {
#' g  <- make_graph( ~ A-B-C-D-A, E-A:B:C:D,
#'    F-G-H-I-F, J-F:G:H:I,
#'    K-L-M-N-K, O-K:L:M:N,
#'    P-Q-R-S-P, T-P:Q:R:S,
#'    B-F, E-J, C-I, L-T, O-T, M-S,
#'    C-P, C-L, I-L, I-P)
#'
#' par("mfrow"=c(2,2));
#' plot(g, main="default layout\n(igraph)");
#'
#' plot(g, main="layout_with_fr\n(igraph)", layout=layout_with_fr);
#'
#' plot(g, main="layout_with_qfr\n(qgraph default)", layout=layout_with_qfr);
#'
#' plot(g, main="layout_with_qfr, repulse=6\n(qgraph custom)",
#'    layout=function(g)layout_with_qfr(g, repulse=6));
#' }
#'
#' @export
layout_with_qfr <- function
(g,
 repulse=3.5,
 area=8*(igraph::vcount(g)^2),
 repulse.rad=(igraph::vcount(g)^repulse),
 constraints=NULL,
 constrain=NULL,
 seed=123,
 weights=NULL,
 niter=NULL,
 max.delta=NULL,
 cool.exp=NULL,
 init=NULL,
 groups=NULL,
 rotation=NULL,
 layout.control=0.5,
 round=TRUE,
 digits=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to apply Fruchterman-Reingold layout from the qgraph
   ## package, which allows tuning some parameters which are no longer
   ## available from the igraph package.
   ##
   ## It also handles the changes made to igraph which produce a character
   ## edgelist instead of numeric edgelist
   e <- igraph::get.edgelist(g, names=FALSE);
   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }

   ## Handle weights from E(g)$weight if supplied
   if (length(weights) == 0 && "weight" %in% igraph::list.edge.attributes(g)) {
      if (verbose) {
         jamba::printDebug("layout_with_qfr(): ",
            "Using E(g)$weight to define weights during layout.");
      }
      weights <- igraph::E(g)$weight;
   }

   # constrain
   if (length(constrain) > 0 && any(constrain %in% igraph::V(g)$name)) {
      if (length(init) == 0) {
         if (verbose) {
            jamba::printDebug("layout_with_qfr(): ",
               "Defining init from graph_attr layout.");
         }
         init <- igraph::graph_attr(g, "layout");
      } else {
         if (verbose) {
            jamba::printDebug("layout_with_qfr(): ",
               "Using init supplied.");
         }
      }
      if (length(constraints) == 0) {
         if (verbose) {
            jamba::printDebug("layout_with_qfr(): ",
               "Defining constraints from init.");
         }
         constraints <- init;
         constraints[] <- NA;
      } else {
         if (verbose) {
            jamba::printDebug("layout_with_qfr(): ",
               "Using init supplied.");
         }
      }
      rownames(init) <- igraph::V(g)$name;
      rownames(constraints) <- igraph::V(g)$name;
      conmatch <- jamba::rmNA(match(constrain, rownames(constraints)));
      if (verbose) {
         jamba::printDebug("layout_with_qfr(): ",
            "Applied constrain to constraints.");
      }
      constraints[conmatch,] <- init[conmatch,];
   }

   frL <- qgraph::qgraph.layout.fruchtermanreingold(
      e,
      vcount=igraph::vcount(g),
      area=area,
      repulse.rad=repulse.rad,
      constraints=constraints,
      weights=weights,
      niter=niter,
      max.delta=max.delta,
      cool.exp=cool.exp,
      init=init,
      groups=groups,
      rotation=rotation,
      layout.control=layout.control,
      round=round,
      digits=digits);
   rownames(frL) <- igraph::V(g)$name;
   frL;
}

#' igraph layout function using qgraph Fruchterman-Reingold
#'
#' igraph layout function using qgraph Fruchterman-Reingold
#'
#' This function returns a layout function, which can be convenient
#' when calling `igraph::plot.igraph()`, in order to set
#' layout parameters in the same call.
#'
#' @return function used to calculate layout coordinates of
#'    an `igraph` object.
#'
#' @family jam igraph functions
#'
#' @param repulse numeric value typically between 3 and 5, passed
#'    to `layout_with_qfr()`, which in turn is passed to
#'    `qgraph::qgraph.layout.fruchtermanreingold()`.
#' @param seed numeric value used to set the R random seed, in order
#'    to make layouts consistent.
#' @param ... additional arguments are passed to `layout_with_qfr()`.
#'
#' @export
layout_with_qfrf <- function
(repulse=3.5,
 seed=123,
 ...)
{
   ## Purpose is to wrap layout_with_qfr() into a function to make
   ## it easier to modify arguments on the fly
   function(g) {
      l <- layout_with_qfr(g=g,
         repulse=repulse,
         seed=seed,
         ...);
      l;
   }
}

#' igraph re-layout using qgraph Fruchterman-Reingold
#'
#' igraph re-layout using qgraph Fruchterman-Reingold
#'
#' This function extends `layout_with_qfr()` by applying the layout
#' to the `igraph` object itself, while also calling
#' `spread_igraph_labels()` to adjust label positions accordingly.
#'
#' The main benefit to using this function is to update the layout
#' and node label positions in one step,
#' while also returning the `igraph` object ready to be plotted as-is.
#'
#' @return `igraph` object, with layout coordinates stored in
#'    graph attribute `"layout"`, accessible for example with
#'    `graph$layout` or `graph_attr(graph, "layout")`.
#'    When `spread_labels=TRUE`,
#'    `V(g)$label.degree` and `V(g)$label.dist` are updated
#'    by calling `spread_igraph_labels()`.
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object
#' @param repulse exponent power used to scale the radius effect around
#'    each vertex. The default is slightly higher than the cube of
#'    the number of vertices, but as the number of vertices increases,
#'    values from 3.5 to 4 and higher are more effective for layout.
#' @param spread_labels logical indicating whether to call
#'    `spread_igraph_labels()`, which places node labels at an angle offset
#'    from the node, in order to improve default label positions.
#' @param ... additional arguments are passed to `layout_with_qfr()` and
#'    `spread_igraph_labels()` as needed.
#'
#' @export
relayout_with_qfr <- function
(g,
 repulse=3.5,
 spread_labels=TRUE,
 seed=123,
 init=NULL,
 constrain=NULL,
 constraints=NULL,
 verbose=FALSE,
 ...)
{
   # if layout exists, use that for init
   if (length(init) == 0 && "layout" %in% igraph::list.graph.attributes(g)) {
      init <- igraph::graph_attr(g, "layout");
      rownames(init) <- igraph::V(g)$name;
      if (verbose) {
         jamba::printDebug("relayout_with_qfr(): ",
            "head(init):");
         print(head(init));
      }
   }
   layout_xy <- layout_with_qfr(g=g,
      repulse=repulse,
      seed=seed,
      init=init,
      constrain=constrain,
      constraints=constraints,
      verbose=verbose,
      ...);
   rownames(layout_xy) <- igraph::V(g)$name;

   if (verbose) {
      jamba::printDebug("relayout_with_qfr(): ",
         "head(layout_xy):");
      print(head(layout_xy));
   }
   g <- igraph::set_graph_attr(g,
      "layout",
      layout_xy);
   if (spread_labels) {
      g <- spread_igraph_labels(g,
         # layout=layout_xy,
         # verbose=verbose,
         ...);
   }
   rownames(igraph::graph_attr(g, "layout")) <- igraph::V(g)$name;
   return(g);
}

#' Create a cnetplot igraph object
#'
#' Create a cnetplot igraph object
#'
#' The purpose of this function is to mimic the steps in `DOSE:::cnetplot()`
#' except not plot the output, and provide some customizations.
#'
#' This function calls `cnetplot_internalJam()`, which among other things
#' adds a node attribute to the resulting `igraph`, `"nodeType"`,
#' where `nodeType="Gene"` identifies gene nodes, and `nodeType="Set"`
#' identifies pathway/gene set nodes.
#'
#' @param x enrichResults object
#' @param showCategory integer number of categories to include in the
#'    resulting Cnet plot.
#' @param categorySize character value indicating how to size the pathway
#'    nodes, where `"geneNum"` sizes nodes by the number of genes in that
#'    pathway, and `"pvalue"` sizes nodes by the enrichment P-value using
#'    the format `-log10(pvalue)`.
#' @param nodeLabel character value indicating which colname in
#'    `as.data.frame(x)` to use as a node label. Depending upon the source
#'    of data, there may be alternative colnames that are more suitable
#'    as node labels.
#' @param foldChange numeric vector named by gene, or NULL. When supplied,
#'    the vector names must use the same nomenclature as the `x` input
#'    object, which can be inspected with `print(head(x@gene))`.
#' @param fixed optional argument passed to `netplot`.
#' @param doPlot logical indicating whether to plot the result.
#' @param categoryColor,geneColor character color, used to colorize
#'    category (pathway) nodes, or gene nodes, respectively.
#' @param normalizeGeneSize logical indicating whether to re-scale the
#'    gene node sizes so the mean gene node size is no larger than the
#'    median category node size. This option is intended to help gene
#'    and category node sizes to be in relative proportion.
#' @param labelCex numeric value to re-scale label text in the igraph
#'    object, and is applied directly to the igraph object.
#' @param colorSub character vector of valid R colors, whose names
#'    are compared to node names, for example `V(g)$name %in% names(colorSub)`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam igraph functions
#'
#' @export
cnetplotJam <- function
(x,
 showCategory=5,
 categorySize="geneNum",
 nodeLabel=c("Name", "Description", "ID"),
 foldChange=NULL,
 fixed=TRUE,
 doPlot=FALSE,
 categoryColor="#E5C494",
 geneColor="#B3B3B3",
 normalizeGeneSize=TRUE,
 labelCex=0.45,
 colorSub=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to run DOSE::cnetplot() but not create a plot.
   ##
   ## categoryColor if supplied is an alternative color for categories.
   ## Note that colorSub can be used to override category colors subsequently.
   ##
   ## geneColor if supplied is an alternative color for categories.
   ## Note that colorSub can be used to override category colors subsequently.
   ##
   ## colorSub if supplied, is a vector of colors, whose names match
   ## V(g)$name. Matching nodes will have their node colors adjusted
   ## based upon colorSub.
   ##
   if (suppressPackageStartupMessages(!require(enrichplot))) {
      stop("cnetplotJam() requires the enrichplot package.");
   }

   y <- as.data.frame(x);
   nodeLabel <- intersect(nodeLabel, colnames(y));
   dColname <- head(nodeLabel, 1);
   if (verbose) {
      jamba::printDebug("cnetplotJam(): ",
         "dColname:",
         dColname);
   }
   keepColnames <- intersect(c(nodeLabel, "pvalue"),
      names(y));
   if (is(x, "enrichResult") || is(x, "gseaResult")) {
      gc <- DOSE::geneInCategory(x);
      names(gc) <- y[[dColname]];
      if (verbose) {
         jamba::printDebug("cnetplotJam(): ",
            "assigned gc <- DOSE::geneInCategory(x), length(gc):",
            length(gc),
            ", head(gc, 2):");
         print(head(gc, 2));
      }
   } else {
      stop("x should be an 'enrichResult' or 'gseaResult' object...")
   }
   if (verbose) {
      jamba::printDebug("cnetplotJam(): ",
         "keepColnames:",
         keepColnames);
   }


   if (is.numeric(showCategory) && (showCategory > length(gc))) {
      showCategory <- length(gc);
   }
   if (categorySize == "pvalue") {
      pvalue <- y$pvalue;
      names(pvalue) <- y[[dColname]];
   } else {
      pvalue <- NULL;
   }
   readable <- x@readable;
   organism <- x@organism;
   if (readable & (!is.null(foldChange))) {
      gid <- names(foldChange);
      if (length(x@gene2Symbol) > 0) {
         if (is(x, "gseaResult")) {
            ii <- gid %in% names(x@geneList);
         } else {
            ii <- gid %in% x@gene;
         }
         gid[ii] <- x@gene2Symbol[gid[ii]];
         names(foldChange) <- gid;
      }
   }

   ## Convert to igraph
   g <- cnetplot_internalJam(inputList=gc,
      showCategory=showCategory,
      categorySize=categorySize,
      pvalue=pvalue,
      foldChange=foldChange,
      fixed=fixed,
      doPlot=doPlot,
      ...);
   igraph::V(g)$frame.color <- jamba::makeColorDarker(igraph::V(g)$color,
      darkFactor=1.5,
      alpha=0.5);
   igraph::V(g)$label.cex <- labelCex;

   ## Match category and gene color
   categoryColorMatch <- "#E5C494";
   geneColorMatch <- "#B3B3B3";
   iWhichCat <- which(igraph::V(g)$color %in% categoryColorMatch);
   iWhichGene <- which(igraph::V(g)$color %in% geneColorMatch);

   #V(g)$color <- ifelse(V(g)$color %in% categoryColorMatch,
   #   categoryColor,
   #   ifelse(V(g)$color %in% geneColorMatch,
   #      geneColor,
   #      V(g)$color));
   if (!all(categoryColor %in% categoryColorMatch) &&
         length(iWhichCat) > 0) {
      igraph::V(g)[iWhichCat]$color <- rep(categoryColor,
         length.out=length(iWhichCat));
   }
   if (!all(geneColor %in% geneColorMatch) &&
         length(iWhichGene) > 0) {
      igraph::V(g)[iWhichGene]$color <- rep(geneColor,
         length.out=length(iWhichGene));
   }

   ## Optionally custom color nodes using colorSub
   if (!is.null(colorSub)) {
      if (any(names(colorSub) %in% igraph::V(g)$name)) {
         iWhich <- which(igraph::V(g)$name %in% names(colorSub));
         if (length(iWhich) > 0) {
            igraph::V(g)[iWhich]$color <- colorSub[igraph::V(g)[iWhich]$name];
         }
      }
   }

   ## Normalize gene and category node sizes
   if (normalizeGeneSize && length(iWhichCat) && length(iWhichGene)) {
      geneSize <- mean(igraph::V(g)[iWhichGene]$size);
      catSizes <- igraph::V(g)[iWhichCat]$size;
      degreeCat <- igraph::degree(g)[iWhichCat];
      #catSize <- sqrt(degreeCat/pi)/2;
      catSize <- sqrt(degreeCat/pi);
      igraph::V(g)[iWhichCat]$size <- catSize;
      if (geneSize > median(catSize)) {
         if (verbose) {
            jamba::printDebug("cnetplotJam(): ",
               "Shrinking gene nodes to median category node size.");
         }
         geneSize <- median(catSize);
         igraph::V(g)[iWhichGene]$size <- geneSize;
      }
   }

   invisible(g);
}

#' cnetplot internal function
#'
#' cnetplot internal function
#'
#' This function is intended to mimic the `DOSE:::cnetplot_internal()`
#' function to support `cnetplotJam()` customizations, including
#' not plotting the output, and including additional custom igraph
#' attributes.
#'
#' @family jam igraph functions
#'
#' @export
cnetplot_internalJam <- function
(inputList,
 categorySize="geneNum",
 showCategory=5,
 pvalue=NULL,
 foldChange=NULL,
 fixed=TRUE,
 DE.foldChange=NULL,
 categoryColor="#E5C494",
 geneColor="#B3B3B3",
 colorRamp="RdBu_r",
 ...)
{
   ## Purpose is to customize DOSE:::cnetplot_internal() to allow
   ## optionally not creating a plot.
   ##
   ## also allow custom categoryColor.
   ## also use vals2colorLevels() instead of DOSE:::get.col.scale() which
   ## does colors nodes with zero fold change using the positive fold change
   ## color gradient, making the positive and negative gradients not
   ## symmetric. It also does not allow specifying the color ramp.
   ##
   ## colorRamp="RdBu_r" uses brewer.pal("RdBu") in reverse,
   ## so blue is low (cold) and red is high (hot).
   ##
   ## categorySize can be numeric vector, or "geneNum" or "pvalue"
   ##
   #categorySize <- match.arg(categorySize);
   if (is.numeric(showCategory)) {
      inputList <- inputList[1:showCategory];
      if (!is.null(pvalue)) {
         pvalue <- pvalue[1:showCategory];
      }
   } else {
      inputList <- inputList[showCategory];
      if (!is.null(pvalue)) {
         pvalue <- pvalue[showCategory];
      }
   }
   g <- list2graph_ggt(inputList);
   #g <- DOSE::setting.graph.attributes(g);
   lengthOfCategory <- length(inputList);

   ## Color gene nodes by fold change if supplied
   if (!is.null(foldChange)) {
      node.idx <- (lengthOfCategory + 1):igraph::vcount(g);
      fcColors <- colorjam::vals2colorLevels(foldChange,
         col=colorRamp,
         divergent=TRUE,
         ...);
      igraph::V(g)[node.idx]$color <- fcColors;
      g <- scaleNodeColor(g, foldChange, node.idx, DE.foldChange);
   }

   igraph::V(g)$size <- 5;
   igraph::V(g)$color <- geneColor;
   igraph::V(g)[seq_len(lengthOfCategory)]$size <- 30;
   igraph::V(g)[seq_len(lengthOfCategory)]$color <- categoryColor;

   ## 0.0.39.900 - update to add nodeType "Set" or "Gene"
   igraph::V(g)$nodeType <- "Gene";
   igraph::V(g)[seq_len(lengthOfCategory)]$nodeType <- "Set";

   ## Size category nodes
   if (is.numeric(categorySize)) {
      ## If supplied a numeric vector, size categories directly
      igraph::V(g)[1:lengthOfCategory]$size <- categorySize;
   } else {
      if (categorySize == "geneNum") {
         n <- igraph::degree(g)[1:lengthOfCategory];
         igraph::V(g)[1:lengthOfCategory]$size <- n/sum(n) * 100;
      } else if (categorySize == "pvalue") {
         if (is.null(pvalue) || any(is.na(pvalue))) {
            stop("pvalue must not be NULL or contain NA values.");
         }
         pScore <- -log10(pvalue);
         igraph::V(g)[1:lengthOfCategory]$size <- pScore/sum(pScore) * 100;
      }
   }
   invisible(g);
}

#' Convert igraph to use pie node shapes
#'
#' Convert igraph to use pie node shapes
#'
#' This function converts an igraph to use pie node shapes,
#' where pie wedges are colored using values derived from a
#' numeric matrix `valueIM` or pre-defined in a character matrix
#' containing colors `valueIMcolors`.
#'
#' Note that pie wedge sizes are equally-sized and do not vary
#' by score, instead the color intensity is applied to each pie
#' wedge.
#'
#' Node names using `V(g)$name` matching `rownames(valueIMcolors)`
#' are colorized and the node shape is converted to pie. All other
#' nodes are not modified.
#'
#' When `valueIMcolors` is not defined, it is derived from `valueIM`
#' using `colorjam::matrix2heatColors()`. In that case, `colorV` defines the color
#' used for numeric values in each column, and other options are passed
#' to `colorjam::matrix2heatColors()` via `...` arguments.
#'
#' @family jam igraph functions
#'
#' @export
igraph2pieGraph <- function
(g,
 valueIM=NULL,
 valueIMcolors=NULL,
 colorV=NULL,
 updateLabels=FALSE,
 maxNchar=62,
 backgroundColor="white",
 seed=123,
 defineLayout=FALSE,
 repulse=3.6,
 removeNA=FALSE,
 NAvalues=c(NA,"transparent"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert an igraph to one using pie nodes, where
   ## wedges are colored using values in an incident matrix valueIM.
   ##
   ## rownames(valueIM) are expected to match V(g)$name
   ## valueIM is used if valueIMcolors is not supplied, it is
   ## given to df2groupColors() to produce valueIMcolors.
   ##
   ## valueIMcolors can be used instead of valueIM, which
   ## will directly assign colors
   ##
   ## TODO:
   ## maintain the names from colnames(valueIM) in the color vectors
   ##
   #V(g)$pie.value <- list(c(1));
   #V(g)$pie.color <- V(g)$color;

   ## Convert valueIM to colors
   if (is.null(valueIMcolors)) {
      if (verbose) {
         jamba::printDebug("igraph2pieGraph(): ",
            "Calling colorjam::matrix2heatColors().");
      }
      valueIMcolors <- colorjam::matrix2heatColors(x=valueIM,
         colorV=colorV,
         verbose=verbose,
         ...);
   }
   if (verbose) {
      jamba::printDebug("igraph2pieGraph(): ",
         "valueIMcolors:");
      print(head(valueIMcolors));
   }

   ## Determine matching node names
   iNodes <- which(toupper(igraph::V(g)$name) %in% toupper(rownames(valueIMcolors)));
   if (length(iNodes) == 0) {
      if (verbose) {
         jamba::printDebug("igraph2pieGraph(): ",
            "No nodes were changed, returning input graph.");
      }
      return(g);
   }
   if (verbose) {
      jamba::printDebug("igraph2pieGraph(): ",
         "found ", format(length(iNodes), big.mark=","),
         " matching nodes, head(iNodes):",
         head(iNodes, 10));
   }

   Vshapes <- igraph::V(g)$shape;
   if (length(Vshapes) == 0) {
      igraph::V(g)$shape <- "circle";
   }

   ## Change above logic to use lapply() then try to assign to igraph
   ## in batch steps
   if (verbose) {
      jamba::printDebug("igraph2pieGraph(): ",
         "iNodeParamsL");
   }
   iNodeParamsL <- lapply(iNodes, function(i){
      iNode <- igraph::V(g)$name[i];
      iV <- match(toupper(iNode), toupper(rownames(valueIMcolors)));

      iPieColor <- valueIMcolors[iV,];
      names(iPieColor) <- colnames(valueIMcolors);
      if (removeNA) {
         iPieColor <- iPieColor[!iPieColor %in% NAvalues];
      }

      iPieValue <- rep(1, length(iPieColor));
      #names(iPieValue) <- names(iPieColor);
      names(iPieValue) <- NULL;

      retVals <- list();
      retVals$pie.color <- iPieColor;
      retVals$pie.value <- iPieValue;
      retVals$pie.names <- colnames(valueIMcolors);
      retVals$shape <- "pie";
      retVals;
   });
   iPieValueL <- lapply(iNodeParamsL, function(i){
      i$pie.value;
   });
   iPieColorL <- lapply(iNodeParamsL, function(i){
      i$pie.color;
   });
   iPieNamesL <- lapply(iNodeParamsL, function(i){
      i$pie.names;
   });
   igraph::V(g)[iNodes]$shape <- "pie";
   igraph::V(g)[iNodes]$pie.value <- iPieValueL;
   igraph::V(g)[iNodes]$pie <- iPieValueL;
   igraph::V(g)[iNodes]$pie.color <- iPieColorL;
   igraph::V(g)[iNodes]$pie.names <- iPieNamesL;

   ## Define stable layout
   if (defineLayout) {
      if (verbose) {
         jamba::printDebug("igraph2pieGraph(): ",
            "layout_with_qfr(g, repulse=",
            repulse,
            ")");
      }
      if (length(seed) > 0) {
         set.seed(head(seed, 1));
      }
      g <- relayout_with_qfr(g,
         repulse=repulse,
         ...);
      #layoutG <- layout_with_qfr(g, repulse=repulse);
      #V(g)$x <- layoutG[,1];
      #V(g)$y <- layoutG[,2];
   }

   ## Optionally update labels for maximum characters
   if (updateLabels) {
      if (verbose) {
         jamba::printDebug("igraph2pieGraph(): ",
            "Updating node labels.");
      }
      if (is.null(igraph::V(g)$label)) {
         igraph::V(g)$label <- igraph::V(g)$name;
      }
      if (!is.null(maxNchar)) {
         igraph::V(g)$label <- substr(igraph::V(g)$label, 1, maxNchar);
      }
      igraph::V(g)$label <- jamba::ucfirst(tolower(igraph::V(g)$label));
   }
   return(g);
}


#' Get angle from origin to vector of x,y coordinates
#'
#' Get angle from origin to vector of x,y coordinates
#'
#' This function gets the angle from origin to x,y coordinates,
#' allowing for vectorized input and output.
#'
#' @param x numeric vector or two-column matrix with columns
#'    representing x,y coordinates when y is `NULL`.
#' @param y numeric vector or `NULL`.
#' @param directed logical indicating whether to return angles around
#'    the full circle, or only half circle. For example, in degrees
#'    `c(1,1)` indicates 45 degrees, `c(-1,-1)` indicates -135 degrees.
#'    When `directed=FALSE` then `c(-1,-1)` indicates 45 degrees.
#' @param deg logical indicating whether to return degrees, or when
#'    `deg=FALSE` it returns radians.
#' @param origin.x,origin.y numeric input defining the coordinates
#'    to use as the origin. When non-zero it implies the first point
#'    of each segment.
#' @param ... additional arguments are ignored.
#'
#' @family jam utility functions
#'
#' @examples
#' # by default output is in degrees
#' xyAngle(1, 1);
#'
#' # output in radians
#' xyAngle(1, 1, deg=FALSE);
#'
#' # optionally different origin
#' xyAngle(1, 1, origin.x=1, origin.y=0);
#'
#' @export
xyAngle <- function
(x,
 y=NULL,
 directed=FALSE,
 deg=TRUE,
 origin.x=0,
 origin.y=0,
 ...)
{
   ## Get angle from zero to given x,y coordinates
   if (length(y) == 0) {
      y <- x[,2];
      x <- x[,1];
   }
   if (length(ncol(origin.x)) > 0) {
      origin.y <- origin.x[,2];
      origin.x <- origin.x[,1];
   }
   out <- base::atan2(y - origin.y,
      x - origin.x);
   if (!directed) {
      out <- out %% pi;
   }
   if (deg) {
      out <- out * 180 / pi;
   }
   out;
}

#' Draw ellipse
#'
#' Draw ellipse
#'
#' This function draws an ellipse centered on the given coordinates,
#' rotated the given degrees relative to the center point, with give
#' x- and y-axis radius values.
#'
#' @return invisible list of x,y coordinates
#'
#' @param x,y `numeric` coordinates, where x can be a two-column numeric
#'    matrix of x,y coordinates.
#' @param a,b `numeric` values indicating x- and y-axis radius, before
#'    rotation if `angle` is non-zero.
#' @param angle `numeric` value indicating the rotation of ellipse.
#' @param segment NULL or `numeric` vector of two values indicating the
#'    start and end angles for the ellipse, prior to rotation.
#' @param arc.only `logical` indicating whether to draw the ellipse
#'    arc without connecting to the center of the ellipse. Set
#'    `arc.only=FALSE` when segment does not include the full circle,
#'    to draw only the wedge.
#' @param nv `numeric` the number of vertices around the center to draw.
#' @param deg `logical` indicating whether input `angle` and `segment`
#'    values are in degrees, or `deg=FALSE` for radians.
#' @param border,col,lty,lwd arguments passed to `graphics::polygon()`.
#' @param draw `logical` indicating whether to draw the ellipse.
#' @param ... additional arguments are passed to `graphics::polygon()`
#'    when `draw=TRUE`.
#'
#' @family jam igraph functions
#'
#' @examples
#' par("mar"=c(2, 2, 2, 2));
#' plot(NULL,
#'    type="n",
#'    xlim=c(-5, 20),
#'    ylim=c(-5, 18),
#'    ylab="", xlab="", bty="L",
#'    asp=1);
#' xy <- drawEllipse(
#'    x=c(1, 11, 11, 11),
#'    y=c(1, 11, 11, 11),
#'    a=c(5, 5, 5*1.5, 5),
#'    b=c(2, 2, 2*1.5, 2),
#'    angle=c(20, -15, -15, -15),
#'    segment=c(0, 360, 0, 120, 120, 240, 240, 360),
#'    arc.only=c(TRUE, FALSE, FALSE, TRUE),
#'    col=jamba::alpha2col(c("red", "gold", "dodgerblue", "darkorchid"), alpha=0.5),
#'    border=c("red", "gold", "dodgerblue", "darkorchid"),
#'    lwd=1,
#'    nv=99)
#' points(x=c(1, 11), y=c(1, 11), pch=20, cex=2)
#' jamba::drawLabels(x=c(12, 3, 13, 5),
#'    y=c(14, 10, 9, 2),
#'    labelCex=0.7,
#'    drawBox=FALSE,
#'    adjPreset=c("topright", "left", "bottomright", "top"),
#'    txt=c("0-120 degrees,\nangle=-15,\narc.only=TRUE",
#'       "120-240 degrees,\nangle=-15,\narc.only=TRUE,\nlarger radius",
#'       "240-360 degrees,\nangle=-15,\narc.only=FALSE",
#'       "angle=20"))
#'
#' @export
drawEllipse <- function
(x,
 y,
 a=1,
 b=1,
 angle=0,
 segment=NULL,
 arc.only=TRUE,
 nv=100,
 deg=TRUE,
 border=NULL,
 col=NA,
 lty=1,
 lwd=1,
 draw=TRUE,
 ...)
{
   ## Purpose is to draw an ellipse
   if (length(deg) > 0 && any(deg %in% TRUE)) {
      deg <- TRUE
   } else {
      deg <- FALSE
   }

   if (length(segment) == 0) {
      if (length(deg) > 0 && any(deg %in% TRUE)) {
         segment <- c(0, 360);
      } else {
         segment <- c(0, pi*2);
      }
   }
   if (length(nv) == 0) {
      nv <- 100;
   } else if (length(nv) > 1) {
      nv <- head(nv, 1)
   }

   ## Fix various vector lengths
   y <- rep(y, length.out=length(x));
   a <- rep(a, length.out=length(x));
   b <- rep(b, length.out=length(x));
   col <- rep(col, length.out=length(x));
   border <- rep(border, length.out=length(x));

   ## if input is in degrees
   if (deg) {
      angle <- angle * pi/180;
      segment <- segment * pi/180;
   }
   segment <- rep(segment,
      length.out=length(x) * 2);
   segment_seq <- seq(from=1, to=length(segment), by=2);
   segment1 <- segment[segment_seq];
   segment2 <- segment[segment_seq + 1];

   angle <- rep(angle,
      length.out=length(segment1));
   if (length(arc.only) == 0) {
      arc.only <- TRUE;
   }
   arc.only <- rep(arc.only,
      length.out=length(segment1));
   if (length(segment1) == 1) {
      z <- seq(from=segment[1],
         to=segment[2],
         length=nv + 1);
      if (!arc.only) {
         z <- c(NA, z, NA, NA);
      }
      z_idx <- rep(1, length(z));
      z_angle <- rep(angle, length.out=length(z));
      z_cumsum <- length(z);
      z_lengths <- length(z);
   } else {
      z_list <- lapply(seq_along(segment1), function(i){
         j <- c(seq(from=segment1[i],
            to=segment2[i],
            length.out=nv + 1), NA);
         if (!arc.only[i]) {
            j <- c(NA, j, NA);
         }
         j
      })
      # jamba::printDebug("z_list:");print(z_list)
      z <- unlist(z_list);
      z_lengths <- lengths(z_list);
      z_idx <- rep(seq_along(z_list), z_lengths);
      z_angle <- rep(angle, z_lengths);
      z_cumsum <- cumsum(z_lengths);
   }
   xx <- a[z_idx] * cos(z);
   yy <- b[z_idx] * sin(z);
   alpha <- xyAngle(xx,
      yy,
      directed=TRUE,
      deg=FALSE);
   rad <- sqrt(xx^2 + yy^2)
   xp <- rad * cos(alpha - z_angle) + x[z_idx];
   yp <- rad * sin(alpha - z_angle) + y[z_idx];
   if (any(!arc.only)) {
      which_wedge <- which(!arc.only);
      # jamba::printDebug("which_wedge: ", which_wedge);
      wedge_x <- x[which_wedge];
      wedge_y <- y[which_wedge];
      wedge_idx1 <- z_cumsum[which_wedge] - z_lengths[which_wedge] + 1;
      wedge_idx2 <- z_cumsum[which_wedge] - 1;
      # jamba::printDebug("wedge_idx1: ", wedge_idx1);
      # jamba::printDebug("wedge_idx2: ", wedge_idx2);
      xp[wedge_idx1] <- wedge_x;
      xp[wedge_idx2] <- wedge_x;
      yp[wedge_idx1] <- wedge_y;
      yp[wedge_idx2] <- wedge_y;
   }
   if (draw) {
      polygon(xp,
         yp,
         border=border,
         col=col,
         lty=lty,
         lwd=lwd,
         ...);
   }
   return(invisible(list(
      x=xp,
      y=yp,
      z=z
   )));
   invisible(list(x=xp,
      y=yp));
}

#' Summarize Cnet igraph as a data.frame
#'
#' Summarize Cnet igraph as a data.frame
#'
#' This function provides a data.frame summary of an igraph object
#' containing "Cnet" data, including vertex attribute `"nodeType"`
#' with values `"Set"` and `"Gene"`, and where `"Set"` nodes are
#' only connected to `"Gene"` nodes.
#'
#' The data.frame is intended to provide a convenient method for
#' subsetting nodes, typically based upon a connected cluster,
#' or the minimum number of edges per node. For example, filter
#' for the connected component containing a node of interest, or
#' filter for `"Set"` nodes with more than one `"Gene"`.
#'
#' @return data.frame with the node name, label, degree (number of
#'    edges), membership (based upon connected component), and
#'    if `getNeighbors=TRUE` it includes comma-delimited names
#'    of neighboring nodes.
#'
#' @family jam igraph functions
#' @family jam conversion functions
#'
#' @param g igraph object containing Cnet data, specifically vertex
#'    attribute name "nodeType" with values "Set" and "Gene", and
#'    where "Set" nodes are only connected to "Gene" nodes.
#' @param getNeighbors logical indicating whether to include
#'    the connected neighbor node names.
#' @param checkSubsets logical indicating whether to test "Set"
#'    nodes to determine if the neighbors are all represented by
#'    another "Set" node.
#' @param ... additional arguments are ignored.
#'
#' @export
cnet2df <- function
(g,
 getNeighbors=TRUE,
 checkSubsets=getNeighbors,
 ...)
{
   ## Purpose is to summarize an igraph object by connectivity,
   ## connected components, and neighbors
   df <- data.frame(nodeType=igraph::V(g)$nodeType,
      name=igraph::V(g)$name,
      label=igraph::V(g)$label,
      degree=igraph::degree(g),
      membership=igraph::components(g)$membership);
   if (getNeighbors || checkSubsets) {
      df$neighbors <- jamba::cPaste(
         lapply(seq_len(igraph::vcount(g)), function(i){
            igraph::neighbors(g, i)$name;
         }),
         doSort=TRUE,
         makeUnique=TRUE);
   }
   if (checkSubsets) {
      im <- cnet2im(df=df)
      ## determine if neighbors for a Set node are completely contained
      ## in another Set node
      imSet <- (t(im) %*% im);
      # isSubset <- (rowSums(imSet >= rowMaxs(imSet)) > 1);
      isSubset <- (rowSums(imSet >= apply(imSet, 1, max, na.rm=TRUE)) > 1);
      df$isSubset <- FALSE;
      if (any(isSubset)) {
         df[match(rownames(imSet), df$name),"isSubset"] <- isSubset;
      }
   }
   df;
}

#' Convert Cnet igraph to incidence matrix
#'
#' Convert Cnet igraph to incidence matrix
#'
#' This function takes igraph object containing "Cnet" data,
#' including vertex attribute `"nodeType"` with values `"Set"`
#' and `"Gene"`, and where `"Set"` nodes are only connected
#' to `"Gene"` nodes. It returns an incidence matrix whose
#' columns are "Set" node names and whose rows are "Gene"
#' node names.
#'
#' @return numeric matrix with colnames defined by `"Set"` node
#'    names, and rownames defined by `"Gene"` node names.
#'
#' @family jam igraph functions
#' @family jam conversion functions
#'
#' @param g igraph object containing Cnet data, specifically vertex
#'    attribute name "nodeType" with values "Set" and "Gene", and
#'    where "Set" nodes are only connected to "Gene" nodes.
#' @param df data.frame as optional input instead of `g`, usually
#'    the result of a previous call to `cnet2df()`.
#' @param ... additional arguments are ignored.
#'
#' @export
cnet2im <- function
(g=NULL,
 df=NULL,
 ...)
{
   ## Purpose is to convert a Cnet igraph to an incidence matrix
   if (length(g) > 0 && "data.frame" %in% class(g)) {
      df <- g;
   } else if (length(df) == 0) {
      df <- cnet2df(g,
         getNeighbors=TRUE,
         checkSubsets=FALSE);
   }
   dfV <- jamba::nameVector(subset(df, nodeType %in% "Set")[,c("neighbors","name")]);
   dfL <- strsplit(
      as.character(dfV),
      ",");
   im <- list2im(dfL);
   im;
}


#' Remove igraph blank wedges
#'
#' Remove igraph blank wedges
#'
#' This function is intended to affect nodes with shapes `"pie"` or
#' `"coloredrectangle"`, and evaluates the vertex attributes
#' `"coloredrect.color"` and `"pie.color"`. For each node, any colors
#' considered blank are removed, along with corresponding values in
#' related vertex attributes, including `"pie","pie.value","pie.names"`,
#' `"coloredrect.names","coloredrect.nrow","coloredrect.ncol","coloredrect.byrow"`.
#'
#' This function calls `isColorBlank()` to determine which colors are
#' blank.
#'
#' This function is originally intended to follow `igraph2pieGraph()` which
#' assigns colors to pie and coloredrectangle attributes, where missing
#' values or values of zero are often given a "blank" color. To enhance the
#' resulting node coloration, these blank colors can be removed in order to
#' make the remaining colors more visibly distinct.
#'
#' @family jam igraph functions
#'
#' @param g igraph object containing one or more attributes from
#'    `"pie.color"` or `"coloredrect.color"`.
#' @inheritParams isColorBlank
#' @param constrain character value indicating for node shape
#'    `"coloredrectangle"` whether to constrain the `"coloredrect.nrow"`
#'    or `"coloredrect.ncol"` values. When `"none"` the nrow is usually
#'    dropped to nrow=1 whenever colors are removed.
#' @param resizeNodes logical indicating whether to resize the resulting
#'    nodes to maintain roughly proportional size to the number of
#'    colored wedges.
#' @param applyToPie logical indicating whether to apply the logic to
#'    nodes with shape `"pie"`.
#' @param pie_to_circle logical indicating whether node shapes for
#'    single-color nodes should be changed from `"pie"` to `"circle"`
#'    in order to remove the small wedge line in each pie node.
#' @param pieAttrs character vector of `vertex.attributes` from `g`
#'    to be adjusted when `applyToPie=TRUE`. Note that `"pie.color"`
#'    is required, and other attributes are only adjusted when
#'    they are present in the input graph `g`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `isColorBlank()`.
#'
#' @examples
#' require(igraph);
#' g <- make_full_graph(n=3);
#' V(g)$name <- c("nodeA", "nodeB", "nodeC");
#' V(g)$shape <- "coloredrectangle";
#' V(g)$coloredrect.names <- split(
#'    rep(c("up","no", "dn"), 7),
#'    rep(V(g)$name, c(2,3,2)*3));
#' V(g)$coloredrect.byrow <- FALSE;
#' V(g)$coloredrect.nrow <- rep(3, 3);
#' V(g)$coloredrect.ncol <- c(2,3,2);
#' V(g)$label.degree <- pi*3/2;
#' V(g)$label.dist <- 3;
#' V(g)$size2 <- c(3, 3, 3);
#'
#' color_v <- rep("white", 21);
#' k <- c(1, 3, 7, 9, 15, 19, 20, 21);
#' color_v[k] <- rep(colorjam::rainbowJam(5), length.out=length(k))
#' V(g)$coloredrect.color <- split(
#'    color_v,
#'    rep(V(g)$name, c(2,3,2)*3));
#' withr::with_par(list(mfrow=c(2, 2)), {
#' lg <- layout_nicely(g);
#' jam_igraph(g, layout=lg, use_shadowText=TRUE);
#'
#' g2 <- removeIgraphBlanks(g, constrain="none");
#' V(g2)$size2 <- V(g2)$size2 / 3;
#' jam_igraph(g2, layout=lg, use_shadowText=TRUE,
#'    main="constrain='none'");
#'
#' g3 <- removeIgraphBlanks(g, constrain="nrow");
#' jam_igraph(g3, layout=lg, use_shadowText=TRUE,
#'    main="constrain='nrow'");
#'
#' g4 <- removeIgraphBlanks(g, constrain="ncol");
#' jam_igraph(g4, layout=lg, use_shadowText=TRUE,
#'    main="constrain='ncol'");
#' })
#' 
#' #
#' g7 <- make_full_graph(n=7);
#' V(g7)$coloredrect.color <- lapply(c(1,2,3,4,2,3,4),
#'    function(i){colorjam::rainbowJam(i)});
#' V(g7)$coloredrect.ncol <- c(1,1,1,1,2,3,4);
#' V(g7)$coloredrect.nrow <- c(1,2,3,4,1,1,1);
#' V(g7)$coloredrect.names <- V(g7)$coloredrect.color;
#' V(g7)$shape <- "coloredrectangle";
#' V(g7)$size <- 10;
#' V(g7)$size2 <- V(g7)$coloredrect.ncol * 1;
#' lg7 <- layout_nicely(g7);
#' jam_igraph(g7, layout=lg7,
#'    use_shadowText=TRUE,
#'    vertex.size2=5);
#'
#' @export
removeIgraphBlanks <- function
(g,
 blankColor=c("#FFFFFF","#FFFFFFFF","transparent"),
 c_max=7,
 l_min=95,
 alpha_max=0.1,
 constrain=c("nrow","ncol","none"),
 resizeNodes=TRUE,
 applyToPie=TRUE,
 pie_to_circle=FALSE,
 pieAttrs=c("pie", "pie.value", "pie.names", "pie.color"),
 verbose=FALSE,
 ...)
{
   ## Remove white from Cnet multinodes
   ##
   ## resizeNodes will proportionally resize nodes based upon the
   ## resulting ncol and nrow.
   ##
   ## 14jun2018: changed to use isColorBlank() helper function,
   ## which helps encapsulate logic regarding nearly-white colors,
   ## and almost fully transparent colors, both of which are intended
   ## to be considered blank for the purposes of this function
   ##
   ## TODO: iterate pie nodes

   constrain <- match.arg(constrain);
   #ixV <- which(V(g)$shape %in% "coloredrectangle");
   ixV <- which(lengths(igraph::V(g)$coloredrect.color) > 0);

   if ("coloredrect.color" %in% igraph::list.vertex.attributes(g)) {
      if (verbose) {
         jamba::printDebug("removeIgraphBlanks(): ",
            "Adjusting coloredrect nodes.");
      }
      ## Rewrote code to use vectorized logic.

      ## Determine the coloredrect.ncol to use in resizing
      ncolVbefore <- igraph::get.vertex.attribute(g, "coloredrect.ncol");
      if (length(ncolVbefore) == 0) {
         ncolVbefore <- lengths(igraph::get.vertex.attribute(g, "coloredrect.color"));
      }
      nrowVbefore <- igraph::get.vertex.attribute(g, "coloredrect.nrow");

      ## Determine which pie wedges are blank
      iCrColorL <- igraph::get.vertex.attribute(g, "coloredrect.color");
      crBlanksL <- isColorBlank(iCrColorL,
         blankColor=blankColor,
         c_max=c_max,
         l_min=l_min,
         alpha_max=alpha_max);

      ## Check for all blanks
      # all_blank <- sapply(crBlanksL, all);
      # if (any(all_blank)) {
      #    #
      # }

      crLengths <- lengths(iCrColorL);
      crSplitV <- rep(factor(seq_len(igraph::vcount(g))), crLengths);
      ## Vector of TRUE,FALSE
      crBlanksV <- unlist(unname(crBlanksL));
      ## Iterate each attribute
      crAttrs <- intersect(c("coloredrect.color", "coloredrect.names"),
         igraph::list.vertex.attributes(g));
      crAttr <- "coloredrect.color";
      crName <- "coloredrect.names";

      crAttrL <- igraph::get.vertex.attribute(g, crAttr);
      crNameL <- igraph::get.vertex.attribute(g, crName);
      ## Confirm each attribute has the same lengths() as pieColorL
      if (!all(lengths(crAttrL) == crLengths)) {
         # Skip this crAttr since its values are
         # not in sync with "coloredrect.color"
         if (verbose) {
            jamba::printDebug("removeIgraphBlanks(): ",
               "Skipped crAttr attribute '",
               crAttr,
               "' because its lengths() were not consistent with ",
               "'pie.color'");
         }
      } else {
         ##
         nrowV <- jamba::rmNULL(igraph::get.vertex.attribute(g, "coloredrect.nrow"),
            nullValue=1);
         ncolV <- jamba::rmNULL(igraph::get.vertex.attribute(g, "coloredrect.ncol"),
            nullValue=1);
         byrowV <- jamba::rmNULL(igraph::get.vertex.attribute(g, "coloredrect.byrow")*1,
            nullValue=TRUE);
         nprodV <- nrowV * ncolV;
         if (any(crLengths != nprodV)) {
            if (verbose) {
               jamba::printDebug("removeIgraphBlanks(): ",
                  "changing constrain to ",
                  "'none'",
                  " since some ncol,nrow were incorrect.");
            }
            constrain <- "none";
         }

         ## Check for shortcuts from different constraints and ncol,nrow
         if ("nrow" %in% constrain && all(nrowV %in% 1)) {
            if (verbose) {
               jamba::printDebug("removeIgraphBlanks(): ",
                  "changing constrain to ",
                  "'none'",
                  " since all nrow=1");
            }
            constrain <- "none";
         }
         if ("ncol" %in% constrain && all(ncolV %in% 1)) {
            if (verbose) {
               jamba::printDebug("removeIgraphBlanks(): ",
                  "changing constrain to ",
                  "'none'",
                  " since all ncol=1");
            }
            constrain <- "none";
         }

         #########################################
         ## Handle each constraint properly
         if ("none" %in% constrain) {
            if (verbose) {
               jamba::printDebug("removeIgraphBlanks(): ",
                  "Applying constrain ",
                  "'none'");
            }
            ########################################
            ## constrain "none"
            crL <- unname(split(unlist(crAttrL)[!crBlanksV],
               crSplitV[!crBlanksV]));
            crLengthsNew <- lengths(crL);
            crChanged <- (crLengths != crLengthsNew);
            ncolV <- ifelse(nrowV == 1 | ncolV > 1, crLengthsNew, ncolV);
            nrowV <- ifelse(nrowV == 1 | ncolV > 1, 1, crLengthsNew);

            g <- igraph::set.vertex.attribute(g,
               name="coloredrect.nrow",
               value=nrowV);
            g <- igraph::set.vertex.attribute(g,
               name="coloredrect.ncol",
               value=ncolV);
            ## TODO: only update nodes that change
            g <- igraph::set.vertex.attribute(g,
               name=crAttr,
               value=crL);
         } else if (any(c("ncol","nrow") %in% constrain)) {
            ########################################
            ## constrain "nrow" or "ncol"
            #
            # for "nrow":
            # constrain the nrow by making a giant wide matrix,
            # find which columns are completely blank and remove
            # only those columns.
            #
            # for "ncol":
            # constrain the ncol by making a giant tall matrix,
            # find which rows are completely blank and remove
            # only those rows.
            #
            # Note that it needs to iterate each unique coloredrect.nrow
            # in order to keep the dimensions correct.
            nrowNcolByrowAll <- paste0(nrowV, "_", ncolV, "_", byrowV);
            nrowNcolByrowU <- unique(nrowNcolByrowAll);
            for (nrowNcolByrowI in nrowNcolByrowU) {
               nrowNcolByrowV <- as.numeric(strsplit(nrowNcolByrowI, "_")[[1]]);
               iUse <- (nrowNcolByrowAll %in% nrowNcolByrowI);
               ## Create the extended matrix for each of four conditions:
               ## constrain="nrow",byrow=TRUE;
               ## constrain="nrow",byrow=FALSE
               ## constrain="ncol",byrow=TRUE;
               ## constrain="ncol",byrow=FALSE
               if ("nrow" %in% constrain) {
                  if (nrowNcolByrowV[3] == 0) {
                     # coloredrect.byrow == FALSE
                     iM <- matrix(nrow=nrowNcolByrowV[1],
                        unlist(unname(crBlanksL[iUse]))*1);
                     iMvals <- matrix(nrow=nrowNcolByrowV[1],
                        unlist(unname(crAttrL[iUse])));
                     iMnames <- matrix(nrow=nrowNcolByrowV[1],
                        unlist(unname(crNameL[iUse])));
                  } else {
                     iM <- do.call(cbind,
                        lapply(crBlanksL[iUse], function(k){
                           matrix(nrow=nrowNcolByrowV[1],
                              byrow=TRUE,
                              k);
                        }));
                     iMvals <- do.call(cbind,
                        lapply(crAttrL[iUse], function(k){
                           matrix(nrow=nrowNcolByrowV[1],
                              byrow=TRUE,
                              k);
                        }));
                     iMnames <- do.call(cbind,
                        lapply(crNameL[iUse], function(k){
                           matrix(nrow=nrowNcolByrowV[1],
                              byrow=TRUE,
                              k);
                        }));
                  }
                  ## keep track of which columns belong to which node
                  iMcol <- factor(rep(which(iUse), crLengths[iUse]/nrowNcolByrowV[1]));
                  ## Find columns where the colMin is 1, meaning all are blank
                  iMblank <- (matrixStats::colMins(iM) == 1);
                  # Subset for non-blank columns
                  iMvalsM <- iMvals[,!iMblank,drop=FALSE];
                  iMnamesM <- iMnames[,!iMblank,drop=FALSE];
                  iMvalsL <- split(as.vector(iMvalsM),
                     rep(iMcol[!iMblank], each=nrow(iMvalsM)));
                  iMnamesL <- split(as.vector(iMnamesM),
                     rep(iMcol[!iMblank], each=nrow(iMnamesM)));
                  iMncol <- lengths(split(iMcol[!iMblank], iMcol[!iMblank]));
                  iMnrow <- nrowNcolByrowV[1];
               } else {
                  ## constrain "ncol"
                  if (nrowNcolByrowV[3] == 0) {
                     byrow <- FALSE;
                  } else {
                     byrow <- TRUE;
                  }
                  iM <- do.call(rbind,
                     lapply(crBlanksL[iUse], function(k){
                        matrix(ncol=nrowNcolByrowV[2],
                           byrow=byrow,
                           k);
                     }));
                  iMvals <- do.call(rbind,
                     lapply(crAttrL[iUse], function(k){
                        matrix(ncol=nrowNcolByrowV[2],
                           byrow=byrow,
                           k);
                     }));
                  iMnames <- do.call(rbind,
                     lapply(crNameL[iUse], function(k){
                        matrix(ncol=nrowNcolByrowV[2],
                           byrow=byrow,
                           k);
                     }));
                  ## keep track of which columns belong to which node
                  iMrow <- factor(rep(which(iUse), crLengths[iUse]/nrowNcolByrowV[2]));
                  ## Find columns where the colMin is 1, meaning all are blank
                  iMblank <- (rowSums(iM) == ncol(iM));
                  # Subset for non-blank rows
                  iMvalsM <- iMvals[!iMblank,,drop=FALSE];
                  iMnamesM <- iMnames[!iMblank,,drop=FALSE];
                  if (byrow) {
                     iMrowSplit <- rep(iMrow[!iMblank], each=ncol(iMvalsM));
                     iMvalsL <- split(as.vector(t(iMvalsM)),
                        iMrowSplit);
                     iMnamesL <- split(as.vector(t(iMnamesM)),
                        iMrowSplit);
                  } else {
                     iMrowSplit <- rep(iMrow[!iMblank], ncol(iMvalsM));
                     iMvalsL <- split(as.vector(iMvalsM),
                        iMrowSplit);
                     iMnamesL <- split(as.vector(iMnamesM),
                        iMrowSplit);
                  }
                  iMnrow <- lengths(split(iMrow[!iMblank], iMrow[!iMblank]));
                  iMncol <- rep(nrowNcolByrowV[2], length(iMnrow));
               }
               iSet <- as.integer(names(iMvalsL));
               g <- igraph::set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.ncol",
                  value=iMncol);
               g <- igraph::set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.nrow",
                  value=iMnrow);
               g <- igraph::set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.color",
                  value=iMvalsL);
               g <- igraph::set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.names",
                  value=iMnamesL);
               #
            }
         }
      }

      ## Now resize coloredrectangle size2 values
      ## so each square is constant size relative to
      ## its expected node size
      if (resizeNodes) {
         if (verbose) {
            jamba::printDebug("removeIgraphBlanks(): ",
               "Resizing coloredrect nodes.");
         }
         ## Make multi-segment gene nodes wider
         ncolVafter <- igraph::get.vertex.attribute(g, "coloredrect.ncol");
         nrowVafter <- igraph::get.vertex.attribute(g, "coloredrect.nrow");
         resizeWhich <- (ncolVbefore != ncolVafter) |  (nrowVbefore != nrowVafter);
         if (any(resizeWhich)) {
            new_size2 <- nrowVbefore / nrowVafter *
               jamba::rmNULL(igraph::get.vertex.attribute(g, name="size2"),
                  nullValue=default_igraph_values()$vertex$size2)
            if (length(new_size2) == 0) {
               new_size2 <- nrowVbefore / nrowVafter *
                  default_igraph_values()$vertex$size2;
            }
            if (verbose > 1) {
               print(data.frame(ncolVbefore,
                  ncolVafter,
                  size2=jamba::rmNULL(
                     igraph::V(g)$size2,
                     nullValue=default_igraph_values()$vertex$size2),
                  new_size2));
            }
            g <- igraph::set.vertex.attribute(g,
               name="size2",
               value=new_size2[resizeWhich],
               index=which(resizeWhich));
         }
      }
   }

   ## Iterate pie nodes
   if (applyToPie) {
      ## Adjust several pie attributes depending upon what is present
      pieAttrs <- intersect(c("pie", "pie.value", "pie.names", "pie.color"),
         igraph::list.vertex.attributes(g));

      if ("pie.color" %in% pieAttrs) {
         if (verbose) {
            jamba::printDebug("removeIgraphBlanks(): ",
               "Iterating pie nodes.");
         }

         ## Determine which pie wedges are blank
         iPieColorL <- igraph::get.vertex.attribute(g, "pie.color");
         pieBlanksL <- isColorBlank(iPieColorL,
            blankColor=blankColor,
            c_max=c_max,
            l_min=l_min,
            alpha_max=alpha_max,
            ...);
         pieLengths <- lengths(iPieColorL);
         pieSplitV <- rep(seq_len(igraph::vcount(g)), pieLengths);
         ## Vector of TRUE,FALSE
         pieBlanksV <- unlist(unname(pieBlanksL));
         ## Iterate each pie attribute
         for (pieAttr in pieAttrs) {
            pieAttrL <- igraph::get.vertex.attribute(g, pieAttr);
            if (length(pieAttrL) > 0) {
               ## Confirm each attribute has the same lengths() as pieColorL
               if (!all(lengths(pieAttrL) == pieLengths)) {
                  # Skip this pieAttr since its values are not in sync with "pie.color"
                  if (verbose) {
                     jamba::printDebug("removeIgraphBlanks(): ",
                        "Skipped pie attribute '",
                        pieAttr,
                        "' because its lengths() were not consistent with ",
                        "'pie.color'");
                  }
               } else {
                  pieL <- split(unlist(pieAttrL)[!pieBlanksV],
                     pieSplitV[!pieBlanksV]);
                  if (verbose > 1) {
                     jamba::printDebug("removeIgraphBlanks(): ",
                        "length(pieL):", length(pieL),
                        ", vcount(g):", vcount(g),
                        ", length(unique(pieSplitV)):", length(unique(pieSplitV)),
                        ", length(unique(pieSplitV[!pieBlanksV])):", length(unique(pieSplitV[!pieBlanksV])),
                        ", sum(!pieBlanksV):", sum(!pieBlanksV)); # debug
                  }
                  if (length(unique(pieSplitV[!pieBlanksV])) < igraph::vcount(g)) {
                     missing_x <- setdiff(seq_len(igraph::vcount(g)), unique(pieSplitV[!pieBlanksV]));
                     if (verbose > 1) {
                        print(missing_x);
                     }
                     # print(head(igraph::V(g)$name[missing_x]))
                     # print(head(igraph::V(g)$pie.color[missing_x]))
                     # print(head(igraph::V(g)$pie.border[missing_x]))
                  }
                  ## TODO: only update nodes that change
                  g <- igraph::set.vertex.attribute(g,
                     index=unique(pieSplitV[!pieBlanksV]),
                     name=pieAttr,
                     value=pieL);
               }
            }
         }
      }
   }
   if (pie_to_circle) {
      is_pie <- igraph::V(g)$shape %in% "pie";
      if (any(is_pie)) {
         is_single_color <- lengths(igraph::V(g)[is_pie]$pie.color) == 1;
         if (any(is_single_color)) {
            switch_nodes <- which(is_pie)[is_single_color];
            i_colors <- unname(unlist(igraph::V(g)[switch_nodes]$pie.color));
            igraph::V(g)[switch_nodes]$color <- i_colors;
            igraph::V(g)[switch_nodes]$shape <- "circle";
         }
      }
   }

   return(g);
}

#' Convert pie igraph node shapes to coloredrectangle
#'
#' Convert pie igraph node shapes to coloredrectangle
#'
#' This function simply converts an igraph network with `"pie"`
#' node shapes, to use the `"coloredrectangle"` node shape
#' provided by the multienrichjam package.
#'
#' In the process, it transfers related node attributes:
#'
#' * `"pie.color"` are copied to `"coloredrect.color"`
#' * `"pie.names"` are copied to `"coloredrect.names"`. The
#' `"coloredrect.names"` can be used to label a color key.
#' * `"size"` is converted to `"size2"` after applying
#' `sqrt(size) * 1.5`. The `"size2"` value is used to
#' define the size of coloredrectangle nodes.
#'
#' @return igraph object where node shapes were changed
#'    from `"pie"` to `"coloredrectangle"`.
#'
#' @family jam igraph functions
#'
#' @param g igraph object, expected to contain one or more
#'    nodes with shape `"pie"`.
#' @param nrow,ncol integer values indicating the default
#'    number of rows and columns to use when displaying
#'    the colors for each node.
#' @param byrow logical indicating whether each vector
#'    of node colors should fill the nrow,ncol matrix
#'    by each row, similar to how values are filled
#'    in `base::matrix()` with argument `byrow`.
#' @param whichNodes integer vector of nodes in `g`
#'    which should be considered. Only nodes with shape
#'    `"pie"` will be converted which are also within
#'    the `whichNodes` vector. By default, all nodes
#'    are converted, but `whichNodes` allows converting
#'    only a subset of nodes.
#' @param ... additional arguments are ignored.
#'
#' @export
rectifyPiegraph <- function
(g,
 nrow=2,
 ncol=5,
 byrow=TRUE,
 whichNodes=seq_len(igraph::vcount(g)),
 ...)
{
   ## Purpose is to convert a piegraph igraph into coloredrectangles,
   ## applying igraph vertex parameters as relevant.
   ##
   ## whichNodes is an integer vector of nodes to use, which is
   ## further filtered for only nodes with V(g)$shape == "pie"
   if (!jamba::igrepHas("igraph", class(g))) {
      stop("Input g must be an igraph object.");
   }
   whichPie <- intersect(whichNodes,
      which(igraph::V(g)$shape %in% "pie"));

   if (length(whichPie) == 0) {
      return(g);
   }

   igraph::V(g)[whichPie]$coloredrect.color <- igraph::V(g)[whichPie]$pie.color;
   igraph::V(g)[whichPie]$coloredrect.names <- igraph::V(g)[whichPie]$pie.names;
   igraph::V(g)[whichPie]$coloredrect.nrow <- nrow;
   igraph::V(g)[whichPie]$coloredrect.ncol <- ncol;
   igraph::V(g)[whichPie]$coloredrect.byrow <- byrow;
   igraph::V(g)[whichPie]$shape <- "coloredrectangle";
   igraph::V(g)[whichPie]$size2 <- sqrt(igraph::V(g)[whichPie]$size)*1.5;
   g;
}

#' Re-order igraph nodes
#'
#' Re-order igraph nodes
#'
#' This function takes an igraph and a layout in the
#' form of coordinates, or a function used to produce
#' coordinates. It repositions nodes within equivalent
#' positions, ordering nodes by color along either the
#' `"x"` or `"y"` direction.
#'
#' Equivalent node positions are those with the same
#' neighboring nodes. For example if node `"A"` and
#' node `"B"` both have neighbors `c("D", "E", "F")`
#' then nodes `"A"` and `"B"` are considered equivalent,
#' and will be reordered by their color.
#'
#' This function is particularly effective with concept
#' network (Cnet) graphs, where multiple terms may
#' be connnected to the same concept. For MultiEnrichmap,
#' it typically works when multiple genes are connected
#' to the same pathways. When this happens, the genes
#' are sorted to group the colors.
#'
#' @return igraph with nodes positioned to order
#' nodes by color. The layout coordinates are stored in
#' the graph attribute `"layout"`, accessible with
#' `g$layout` or `graph_attr(g, "layout")`.
#' When there are not multiple nodes sharing
#' the same neighbors, the original igraph object is
#' returned, with the addition of layout coordinates.
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object, typically expected to have a fixed
#'    graph layout stored as `igraph::graph_attr(g, "layout")`,
#'    or supplied via `layout` argument.
#' @param sortAttributes `character` vector of node attribute
#'    names, to be applied in order when sorting nodes.
#' @param nodeSortBy `character` vector containing `"x"` and
#'    `"y"` indicating the primary axis used to sort nodes.
#'    Note that sort order can be reversed by prepending "-",
#'    for example `"-x"` or `"-y"`.
#' @param orderByAspect `logical` indicating whether the aspect ratio
#'    of each nodeset determines the sort order. When `orderByAspect=TRUE`
#'    it ignores `nodeSortBy` and uses each nodeset aspect ratio
#'    to determine the order. Note that the original values for `nodeSortOrder`
#'    are retained, specifically the directional sign for each axis,
#'    such that `"x"` and `"-y"` will either become `c("x", "-y")` or
#'    `c("-y", "x")`. This way, other locales for which bottom-right
#'    to top-left ordering is more natural, use `nodeSortBy=c("-x", "y")`
#'    to preserve that order. It will become `c("y", "-x")` for tall
#'    aspect nodesets. See `aspectThreshold` below.
#' @param aspectThreshold `numeric` indicating the aspect ratio
#'    threshold below which nodes are sorted `c("-y", "x")` top to bottom.
#'    All aspect ratios (x/y) above this threshold are sorted
#'    `c("x", "-y")` left to right. Aspect ratios with x < y are
#'    internally converted to negative ratios (-y/x) such that -1.25
#'    is equivalent to x/y=0.8 for example.
#' @param layout `numeric` matrix of node coordinates, or
#'    function used to produce layout coordinates. When layout
#'    is `NULL`, this function tries to use graph attribute
#'    `igraph::graph_attr(g, "layout")`, otherwise
#'    the `relayout_with_qfr()` is called.
#' @param nodesets `character` with optional subset of nodesets to
#'    apply re-ordering. Each value must match names generated
#'    by `get_cnet_nodeset()`, otherwise it will be ignored.
#' @param colorV optional `character` vector that contains R colors,
#'    used to order the colors in attributes such as `"pie.color"`
#'    and `"coloredrect.color"`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' require(igraph)
#'    c3 <- c("red", "gold", "blue");
#'    c3l <- list(c3[1], c3[2], c3[3],
#'       c3[c(1,2)], c3[c(1,3)], c3[c(2,3)],
#'       c3[c(1,2,3)]);
#'    set.seed(123);
#'    pc <- c(c3l[1], sample(rep(c3l, c(6,5,5, 4, 1, 4, 4))))
#'    x <- lapply(pc, function(i){
#'       jamba::nameVector(i, paste0("group_", i))
#'    })
#'    g2 <- igraph::graph_from_edgelist(directed=FALSE,
#'       as.matrix(data.frame(
#'          node1=rep("Pathway", length(x)),
#'          node2=paste0("Gene", jamba::colNum2excelName(seq_along(x))))));
#'    V(g2)$pie.color <- x[c(1,seq_along(pc))];
#'    V(g2)$shape <- "pie";
#'    V(g2)$pie <- lapply(lengths(V(g2)$pie.color), function(i){
#'       rep(1, i)
#'    });
#'    V(g2)$frame.color <- "grey80";
#'    V(g2)$pie.border <- NA;
#'    V(g2)$color <- lapply(V(g2)$pie.color, colorjam::blend_colors)
#'
#'    g2 <- relayout_with_qfr(g2, repulse=7, do_reorder=FALSE);
#'    g2b <- spread_igraph_labels(g2, do_reorder=FALSE)
#'    igraph::V(g2b)$label.family <- "Arial"
#'
#'    withr::with_par(list(mar=c(1, 1, 4, 1), xpd=TRUE), {
#'    jam_igraph(g2b,
#'       main="Unordered",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE)
#'    jam_igraph(reorderIgraphNodes(g2b),
#'       main="reorder_igraph_nodes()",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'    jam_igraph(
#'       reorderIgraphNodes(g2b, nodeSortBy=c("-y","x"), orderByAspect=FALSE),
#'       main='reorderIgraphNodes(nodeSortBy=c("-y", "x"))',
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'
#'    jam_igraph(
#'       reorderIgraphNodes(g2b,
#'          nodeSortBy=c("-y", "x"), orderByAspect=FALSE,
#'          sortAttributes=c("-pie.color.length", "pie.color",
#'             "color", "label", "name")),
#'       main="reorder_igraph_nodes() by pie.color.length",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'    })
#'
#'    g2c <- g2b;
#'    set.seed(12)
#'    V(g2c)$frame.color <- sample(c("firebrick3", "#DDDDDD", "dodgerblue3"),
#'       replace=TRUE, size=igraph::vcount(g2c))
#'    V(g2c)$frame.width <- 3;
#'    V(g2c)$pie.lwd <- 3;
#'    withr::with_par(list("lwd"=4, mar=c(1, 1, 4, 1), xpd=TRUE), {
#'    jam_igraph(reorderIgraphNodes(g2c,
#'       nodeSortBy=c("-y", "x")),
#'       main="reorder_igraph_nodes() including frame.color",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'    })
#'
#'    g2d <- reorderIgraphNodes(g2b);
#'    set.seed(12)
#'    mn <- (lengths(V(g2d)$pie.color) > 1);
#'    V(g2d)[!mn]$frame.color <- sample(c("firebrick3", "#DDDDDD", "dodgerblue3"),
#'       replace=TRUE, size=sum(!mn))
#'    V(g2d)$pie.border <- rep(list(character(0)), vcount(g2d))
#'    V(g2d)[mn]$pie.border <- lapply(which(mn), function(i){
#'       jamba::nameVector(
#'          sample(c("firebrick3", "#DDDDDD", "dodgerblue3"),
#'             replace=TRUE, size=lengths(V(g2d)[i]$pie.color)),
#'          names(V(g2d)[i]$pie.color[[1]]))
#'    })
#'    g2e <- reorderIgraphNodes(g2d,
#'       nodeSortBy=c("-y", "x"));
#'    withr::with_par(list("lwd"=4, mar=c(1, 1, 4, 1), xpd=TRUE), {
#'    withr::with_options(list("inner_pie_border"=TRUE), {
#'    jam_igraph(g2e,
#'       main="reorder_igraph_nodes() including frame.color",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'    })})
#'
#'    # frame color sorting
#'    g2f <- g2e;
#'    igraph::V(g2f)["GeneV"]$frame.color <- "green";
#'    igraph::V(g2f)["GeneE"]$frame.color <- "green";
#'    withr::with_par(list("lwd"=5, mar=c(1, 1, 4, 1), xpd=TRUE), {
#'    withr::with_options(list("inner_pie_border"=TRUE), {
#'    jam_igraph(g2f,
#'       main="reorder_igraph_nodes() including frame.color",
#'       label_dist_factor=3,
#'       label_factor=0.7,
#'       node_factor=1.2,
#'       use_shadowText=TRUE);
#'    })})
#'
#'
#' @export
reorderIgraphNodes <- function
(g,
 sortAttributes=c("pie.color",
    "pie.color.length",
    "pie.border",
    "pie.border.length",
    "coloredrect.color",
    "coloredrect.border",
    "color",
    "frame.color",
    "label",
    "name"),
 nodeSortBy=c("x",
    "-y"),
 orderByAspect=TRUE,
 aspectThreshold=-1.25,
 layout=NULL,
 nodesets=NULL,
 colorV=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to reorder nodes based upon some sortable metric.
   ## Logic is as follows:
   ## - all nodes having identical edges are grouped, e.g. if 20
   ##   nodes all have an edge to node "K" and no other edges, they
   ##   are in the same group.
   ##   Similarly, all nodes having an edge only to nodes "K" and "L"
   ##   will be grouped.
   ## - once nodes are grouped, they are re-ordered within that group
   ##   using something like top-to-bottom coordinate, based upon the
   ##   sort metric nodeSortBy=c("y","x").
   ## The desired result for example, if a set of nodes are colored
   ## red or blue, they should be visibly grouped together by that color.
   ##

   # use wrapper function to determine layout
   layout <- get_igraph_layout(g,
      layout=layout,
      verbose=verbose,
      ...)
   # ensure layout is stored in the resulting igraph object
   g <- igraph::set_graph_attr(graph=g,
      name="layout",
      value=layout);
   if (verbose) {
      jamba::printDebug("head(layout, 10):");
      print(head(layout, 10));
   }

   ## comma-delimited neighboring nodes for each node
   g_nodesets <- get_cnet_nodeset(g, filter_set_only=FALSE);
   # names(g_nodesets) <- jamba::makeNames(substr(names(g_nodesets), 1, 25));
   g_nodesets_v <- jamba::nameVector(
      rep(names(g_nodesets), lengths(g_nodesets)),
      unlist(g_nodesets));
   neighborG <- g_nodesets_v[match(igraph::V(g)$name, names(g_nodesets_v))]
   # names(neighborG) <- seq_len(igraph::vcount(g));

   ## Determine which edge groups are present multiple times
   neighborGct <- jamba::tcount(neighborG, minCount=2);
   if (length(neighborGct) == 0) {
      if (verbose) {
         jamba::printDebug("reorderIgraphNodes(): ",
            "found no edge groups, returning input graph unchanged.");
      }
      return(g);
   }

   # single-color attributes
   color_attrs <- c("color",
      "frame.color")
   # multi-color attributes
   multicolor_attrs <- c("coloredrect.color",
      "pie.color",
      "coloredrect.border",
      "pie.border")
   # attributes that may have length
   length_suffices <- c("length",
      "len",
      "n");
   length_attrs <- paste0(
      rep(multicolor_attrs,
         each=length(length_suffices)),
      ".", length_suffices);

   # Use one or more vertex attributes for grouping
   v_attrs <- igraph::vertex_attr_names(g);
   v_attrs_length <- intersect(v_attrs, multicolor_attrs);
   if (length(v_attrs_length) > 0) {
      v_attrs_length_use <- paste0(
         rep(v_attrs_length,
            each=length(length_suffices)),
         ".", length_suffices);
      v_attrs <- unique(c(v_attrs, v_attrs_length_use));
   }

   # validate sortAttributes, also get reverse
   sortOrders <- ifelse(grepl("^[-]", sortAttributes), TRUE, FALSE);
   sortAttributes <- gsub("^[-]", "", sortAttributes);
   keep_attrs <- (!duplicated(sortAttributes) & sortAttributes %in% v_attrs);
   sortOrders <- sortOrders[keep_attrs];
   sortAttributes <- sortAttributes[keep_attrs];
   if (length(sortAttributes) == 0) {
      if (verbose) {
         jamba::printDebug("reorderIgraphNodes(): ",
            "No sortAttributes matched the igraph object, returning g.");
      }
      return(g);
   }
   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "Applying sort to each of ", length(sortAttributes), " sortAttributes.");
   }

   neighborA_df <- do.call(cbind, lapply(sortAttributes,
      function(sortAttribute){
         if (verbose) {
            jamba::printDebug("reorderIgraphNodes(): ",
               "sortAttribute: ", sortAttribute);
         }
         sortOrder <- sortOrders[sortAttribute];
         if (sortAttribute %in% length_attrs) {
            # length attributes convert values to count before sorting
            length_pattern <- paste0("[.](", paste(length_suffices, collapse="|"), ")$");
            use_sortAttribute <- gsub(length_pattern, "", sortAttribute);
            j <- jamba::padInteger(lengths(
               igraph::vertex_attr(g, use_sortAttribute)));
         } else {
            j <- jamba::rmNULL(igraph::vertex_attr(g, sortAttribute),
               nullValue="#555555");
         }
         names(j) <- seq_len(igraph::vcount(g));

         if (sortAttribute %in% length_attrs) {
            jString <- j;
         } else if (sortAttribute %in% c(color_attrs, multicolor_attrs)) {
            # color or multi-color attributes
            j_colors <- jamba::rmNULL(
               igraph::vertex_attr(g, sortAttribute),
               nullValue="#555555");
            # convert to list to handle non-list attributes
            if (!is.list(j_colors)) {
               j_colors <- as.list(j_colors)
            }
            j_colors_u <- jamba::rgb2col(col2rgb(
               unique(unlist(unname(j_colors)))))
            # if only one value is present, return dummy column
            if (length(j_colors_u) == 1) {
               if (verbose) {
                  jamba::printDebug("reorderIgraphNodes(): ",
                     "Only one value for sortAttribute:",
                     sortAttribute);
               }
               jdf <- data.frame(check.names=FALSE,
                  rep(1, length(j_colors)));
               colnames(jdf) <- sortAttribute;
               return(jdf);
            }
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  "new color logic for sortAttribute:",
                  sortAttribute);
            }
            colorVhex <- NULL;
            if (length(colorV) > 0) {
               colorVhex <- jamba::rgb2col(col2rgb(colorV));
            }
            if (!all(j_colors_u %in% colorVhex)) {
               # if not all colors are defined in colorV
               # define colorV using colors_from_list()
               colorV <- colors_from_list(j_colors,
                  verbose=verbose);
            } else if (!all(j_colors_u %in% colorV)) {
               # not they do not match without converting to hex
               # then convert both to hex upfront
               j_colors <- lapply(j_colors, function(ji){
                  jamba::rgb2col(col2rgb(ji))
               })
               colorV[] <- colorVhex;
            } else {
               if (verbose) {
                  jamba::printDebug("reorderIgraphNodes(): ",
                     "Using the supplied colorV:",
                     names(colorV),
                     fgText=list("darkorange1", "dodgerblue", NA),
                     bgText=list(NA, NA, colorV));
               }
            }
            # convert to factor, using colorV as factor levels in order
            # 0.0.67.900: use unique(colorV) to allow for reused colors
            colorattrm <- jamba::rbindList(lapply(j_colors, function(ji){
               factor(ji,
                  levels=unique(colorV))
            }), fixBlanks=TRUE);
            colorattrm <- matrix(ncol=ncol(colorattrm),
               as.numeric(colorattrm));
            # order:
            # mean rank of colors
            # - therefore c(1,4) and c(2,3) would be tied
            # length of colors
            # - lowest to highest
            # rank of colors in order
            # - therefore c(1,4) would appear before c(2,3)
            colorattrm2 <- cbind(
               rowMeans=round(digits=2,
                  rowMeans(colorattrm, na.rm=TRUE)),
               lengths=lengths(j_colors),
               colorattrm);
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  "head(colorattrm2):");
               print(head(colorattrm2));
            }
            colorattrlevels <- unique(jamba::pasteByRow(
               jamba::mixedSortDF(colorattrm2, na.last=FALSE)))
            colorattrfactor <- factor(
               jamba::pasteByRow(colorattrm2),
               levels=colorattrlevels)
            jString <- paste0(sortAttribute,
               jamba::padInteger(as.integer(colorattrfactor)));
         } else if (sortAttribute %in% c("color")) {
            # 0.0.67.900 - this whole section is ignored
            # in favor of re-using the multi-color sort order also
            # for single-color sorting
            j_colors <- igraph::vertex_attr(g, sortAttribute);
            j_colors_u <- jamba::rgb2col(col2rgb(unique(unlist(unname(j_colors)))))
            colorVhex <- NULL;
            if (length(colorV) > 0) {
               colorVhex <- jamba::rgb2col(col2rgb(colorV));
            }
            if (!all(j_colors_u %in% colorVhex)) {
               # if not all colors are defined in colorV
               # define colorV using colors_from_list()
               colorV <- colors_from_list(j_colors,
                  verbose=verbose);
            } else if (!all(j_colors_u %in% colorV)) {
               # not they do not match without converting to hex
               # then convert both to hex upfront
               j_colors <- lapply(j_colors, function(ji){
                  jamba::rgb2col(col2rgb(ji))
               })
               colorV[] <- colorVhex;
            } else {
               if (verbose) {
                  jamba::printDebug("reorderIgraphNodes(): ",
                     "Using the supplied colorV:",
                     names(colorV),
                     fgText=list("darkorange1", "dodgerblue", NA),
                     bgText=list(NA, NA, colorV));
               }
            }

            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  "avg_colors_by_list for ",
                  length(j_colors),
                  " colors");
            }
            j_colors_v <- avg_colors_by_list(j_colors);
            j_sorted <- colorjam::sort_colors(j_colors_v,
               byCols=c("H", "-C", "-L"));
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  c("head(j_colors_v):", head(j_colors_v)));
               jamba::printDebug("reorderIgraphNodes(): ",
                  c("head(j_sorted):", head(j_sorted)));
            }
            j_rank <- match(j_colors_v, unique(j_sorted));
            jString <- factor(j_sorted,
               levels=unique(j_sorted));
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  c("head(jString):", head(jString)));
               jamba::printDebug("reorderIgraphNodes(): ",
                  c("head(j_sorted):", head(j_sorted)));
            }
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  "head(jString):",
                  head(jString));
            }
         } else {
            # all other non-color, and non-length sorting here
            if (jamba::igrepHas("list", class(j))) {
               # convert to comma-delimited string
               jString <- jamba::cPaste(j);
            } else if (jamba::igrepHas("factor", class(j))) {
               jString <- j;
            } else if (jamba::igrepHas("numeric|integer|float", class(j))) {
               jString <- round(j, digits=2);
            } else {
               jString <- j;
            }
         }
         # names(jString) <- seq_len(igraph::vcount(g));
         names(jString) <- igraph::V(g)$name;
         jdf <- data.frame(jString);
         colnames(jdf) <- sortAttribute;
         jdf;
      }
   ));
   # use jamba::mixedSortDF()
   neighborA_df_sorted <- jamba::mixedSortDF(neighborA_df,
      byCols=sortAttributes,
      decreasing=sortOrders);
   # new_order <- as.integer(rownames(neighborA_df_sorted));

   neighborA <- jamba::pasteByRow(neighborA_df, sep="_");
   neighborA_sorted <- jamba::pasteByRow(unique(neighborA_df_sorted), sep="_");
   neighborA <- factor(neighborA,
      levels=neighborA_sorted);
   names(neighborA) <- rownames(neighborA_df);

   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborA_df):");
      print(head(neighborA_df));
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborA):");
      print(head(neighborA));
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborG):");
      print(head(neighborG));
   }

   ## data.frame with the attribute sort, and the node sort
   layout_match <- match(names(neighborG),
      rownames(layout));
   neighborDF <- data.frame(
      vertex=names(neighborG),
      #vertex=names(neighborG),
      edgeGroup=neighborG,
      sortAttribute=neighborA[names(neighborG)],
      x=layout[layout_match, 1],
      y=layout[layout_match, 2]);

   # optionally include label.degree
   if ("label.degree" %in% igraph::vertex_attr_names(g)) {
      neighborDF$label.degree <- igraph::vertex_attr(g, "label.degree");
   }
   # optionally include label.dist
   if ("label.dist" %in% igraph::vertex_attr_names(g)) {
      neighborDF$label.dist <- igraph::vertex_attr(g, "label.dist");
   }

   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborDF):");
      print(head(neighborDF));
   }

   ## The following code iterates each edge group and reassigns
   ## layout coordinates by nodeSortBy axis order.
   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "nodeSortBy:",
         nodeSortBy);
      jamba::printDebug("reorderIgraphNodes(): ",
         "names(neighborGct):",
         paste0('"', names(neighborGct), '"'));
   }

   # optional nodesets
   if (length(nodesets) == 0) {
      nodesets <- names(neighborGct);
   }
   if (length(nodesets) > 0) {
      nodesets <- intersect(nodesets,
         names(neighborGct));
      if (length(nodesets) == 0) {
         if (verbose) {
            jamba::printDebug("reorderIgraphNodes(): ",
               "None of the provided nodesets need to be re-ordered, returning g.");
         }
         # no given nodesets match, therefore we have nothing to do
         return(g)
      }
   }
   if (verbose && length(nodesets) < length(neighborGct)) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "applying to subset of nodesets: ",
         paste0('"', nodesets, '"'));
   }

   if (!any(grepl("x", nodeSortBy))) {
      nodeSortBy <- c(nodeSortBy, "x");
   }
   if (!any(grepl("y", nodeSortBy))) {
      nodeSortBy <- c(nodeSortBy, "-y");
   }
   use_nodeSortBy <- nodeSortBy;
   newDF <- jamba::rbindList(lapply(names(neighborGct), function(Gname){
      iDF <- subset(neighborDF, edgeGroup %in% Gname);
      if (!Gname %in% nodesets) {
         return(iDF)
      }
      if (TRUE %in% orderByAspect && length(aspectThreshold) > 0) {
         if (nrow(iDF) <= 2) {
            xyaspect <- NA;
            nodeSortBy <- jamba::provigrep(c("x", "y"),
               use_nodeSortBy)
         } else {
            xyrange <- apply(iDF[,c("x", "y")], 2, range, na.rm=TRUE)
            j <- apply(xyrange, 2, diff, na.rm=TRUE)
            xyaspect <- unname(ifelse(j[1] >= j[2],
               j[1] / j[2],
               -j[2] / j[1]));
            if (xyaspect >= aspectThreshold) {
               nodeSortBy <- jamba::provigrep(c("x", "y", "."),
                  use_nodeSortBy)
               if (verbose) {
                  jamba::printDebug("reorder_igraph_nodes(): ",
                     "short-wide aspect ratio, using nodeSortBy: ",
                     paste0('"', nodeSortBy, '"'));
               }
            } else {
               nodeSortBy <- jamba::provigrep(c("y", "x", "."),
                  use_nodeSortBy)
               if (verbose) {
                  jamba::printDebug("reorder_igraph_nodes(): ",
                     "tall-skinny aspect ratio, using nodeSortBy: ",
                     paste0('"', nodeSortBy, '"'));
               }
            }
         }
         if (verbose) {
            jamba::printDebug("reorderIgraphNodes(): ",
               iDF[1, "edgeGroup"],
               " xyaspect: ",
               format(xyaspect, digits=3),
               ", nodeSortBy: ",
               nodeSortBy);
         }
      }
      xyOrder <- jamba::mixedSortDF(iDF,
         byCols=nodeSortBy);

      nodeOrder <- jamba::mixedSortDF(iDF,
         byCols=match(c("sortAttribute", "vertex"), colnames(iDF)));

      nodeOrder[,c("x", "y")] <- xyOrder[,c("x", "y"), drop=FALSE];
      if ("label.degree" %in% colnames(iDF)) {
         nodeOrder[,"label.degree"] <- xyOrder[,"label.degree"]
      }
      if ("label.dist" %in% colnames(iDF)) {
         nodeOrder[,"label.dist"] <- xyOrder[,"label.dist"]
      }
      # If there are repeated sortAttributes, we use them to place subsets
      # of nodes top to bottom within each group of coordinates
      # 0.0.67.900 - ignore this section for now
      if (FALSE) {
         if (length(jamba::tcount(nodeOrder[,"sortAttribute"], minCount=2)) > 0) {
            nodeOrder <- jamba::rbindList(lapply(split(nodeOrder, nodeOrder[,"sortAttribute"]), function(jDF){
               if (nrow(jDF) > 1) {
                  byCols <- match(rev(nodeSortBy), colnames(jDF));
                  if (nodeSortBy[2] %in% "y") {
                     #byCols <- byCols * c(-1,1);
                     byCols <- byCols * c(-1,-1);
                  } else {
                     byCols <- byCols * c(1,1);
                  }
                  jDFcoord <- jamba::mixedSortDF(jDF,
                     byCols=byCols);
                  jDF[,c("x","y")] <- jDFcoord[,c("x","y")];
               }
               jDF;
            }));
            rownames(nodeOrder) <- nodeOrder$vertex;
         }
      }
      nodeOrder;
   }));
   iMatch <- match(newDF$vertex, neighborDF$vertex);
   neighborDF[iMatch, c("x", "y")] <- newDF[,c("x", "y")];
   if ("label.degree" %in% colnames(newDF)) {
      neighborDF[iMatch, c("label.degree")] <- newDF[,c("label.degree")];
      igraph::vertex_attr(g, "label.degree") <- neighborDF[,"label.degree"];
   }
   if ("label.dist" %in% colnames(newDF)) {
      neighborDF[iMatch, c("label.dist")] <- newDF[,c("label.dist")];
      igraph::vertex_attr(g, "label.dist") <- neighborDF[,"label.dist"];
   }
   # subDF <- subset(neighborDF, edgeGroup %in% neighborDF[1, "edgeGroup"]);
   new_layout <- as.matrix(neighborDF[, c("x", "y"), drop=FALSE]);
   # re-apply node names as rownames
   rownames(new_layout) <- igraph::V(g)$name;

   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(new_layout):");
      print(head(new_layout));
   }

   g <- igraph::set_graph_attr(g, "layout", new_layout);
   return(g);
}


#' @rdname reorderIgraphNodes
#' @export
reorder_igraph_nodes <- reorderIgraphNodes


#' Remove igraph singlet nodes
#'
#' Remove igraph singlet nodes
#'
#' This function is a lightweight method to remove igraph nodes with
#' no connections. In fact, the `min_degree` can be used to require
#' a minimum number of connections, but the intended use is to remove
#' the singlet nodes that have no connections.
#'
#' @family jam igraph functions
#'
#' @param g igraph object
#' @param min_degree numeric threshold with the minimum number of
#'    connections, also known as the "degree", required for each node.
#' @param ... additional arguments are ignored.
#'
#' @export
removeIgraphSinglets <- function
(g,
 min_degree=1,
 ...)
{
   keep_nodes <- (igraph::degree(g) >= min_degree);
   g_new <- subgraph_jam(g,
      which(keep_nodes));
   return(g_new);
}

#' Spread igraph node labels by angle from node center
#'
#' Spread igraph node labels by angle from node center
#'
#' This function uses the igraph vertex attribute
#' `"label.degree"`, which describes the angular offset for
#' each vertex label. (The `"label.degree"` values are
#' represented as radians, not degrees, starting at 0 for
#' right, and proceeding clockwise starting from the right,
#' down, left, top, right.)
#'
#' This function requires a network layout, which must be fixed
#' in order for the vertex labels to be properly oriented.
#' Labels are oriented opposite the most dominant angular mean
#' of edges from each network node. Typically the side of a node
#' with the fewest edges has the most space to place a label.
#' No further checks are performed for overlapping labels.
#'
#' Note that this function only modifies the other important
#' attribute `"label.dist"` when `label_min_dist`` is defined,
#' in order to enforce a minimum label distance from the center
#' of each node. There is no other logic to position small or
#' large labels to avoid overlapping labels.
#'
#' @family jam igraph functions
#'
#' @param g igraph object
#' @param layout numeric matrix representing the x and y
#'    coordinates of each node in `g`, in the same order as `V(g)`.
#'    When `layout` is not supplied, nodes are checked for
#'    attributes `c("x", "y")` which define a fixed internal
#'    layout. When `force_layout=TRUE` these coordinates are ignored.
#'    If that is not supplied, then `layout_with_qfr()`
#'    is called along with the `repulse` argument. Subsequent
#'    coordinates are stored in `V(g)$x` and `V(g)$y` when
#'    argument `update_g_coords=TRUE`.
#' @param y_bias numeric value indicating the tendency to spread
#'    labels on the y-axis rather than symmetrically around each node.
#'    This argument elongates the circle surrounding a node into
#'    an ellipse with this ratio.
#' @param update_g_coords logical indicating whether the layout
#'    coordinates will be stored in `graph_attr(g, "layout")`.
#' @param do_reorder logical indicating whether to call
#'    `reorderIgraphNodes()` which re-distributes equivalent nodes
#'    based upon the node color(s). A node is "equivalent" to another
#'    node if both nodes have identical edges.
#' @param sortAttributes,nodeSortBy arguments passed to
#'    `reorderIgraphNodes()` when `do_reorder=TRUE`.
#' @param repulse argument passed to `layout_with_qfr()` only
#'    when `layout` is not supplied, and the layout is not stored
#'    in `c(V(g)$x, V(g)$y)`.
#' @param force_relayout logical indicating whether the `igraph` layout
#'    should be recalculated, in order to override coordinates that
#'    may be previously stored in the `igraph` object itself.
#'    Note that when `layout` is supplied, it is always used.
#' @param label_min_dist numeric value used to ensure all labels are
#'    at least some distance from the center. These units are defined
#'    by igraph, and are roughly in units of one line height of text.
#' @param ... additional arguments are passed to `layout_with_qfr()`
#'    when needed.
#'
#' @export
spread_igraph_labels <- function
(g,
 layout=NULL,
 y_bias=1,
 update_g_coords=TRUE,
 do_reorder=TRUE,
 sortAttributes=NULL,
 nodeSortBy=c("x", "-y"),
 repulse=3.5,
 force_relayout=FALSE,
 label_min_dist=0.5,
 verbose=FALSE,
 ...)
{
   ##
   if (verbose) {
      jamba::printDebug("spread_igraph_labels(): ",
         "vcount:", igraph::vcount(g));
   }
   if (length(layout) == 0) {
      if (!force_relayout) {
         if ("layout" %in% igraph::list.graph.attributes(g)) {
            if (verbose) {
               jamba::printDebug("spread_igraph_labels(): ",
                  "Using ","layout"," from graph attributes.");
            }
            layout <- g$layout;
         } else if (all(c("x", "y") %in% igraph::list.vertex.attributes(g))) {
            if (verbose) {
               jamba::printDebug("spread_igraph_labels(): ",
                  "Using ","x,y"," from vertex attributes.");
            }
            layout <- cbind(x=igraph::V(g)$x, y=V(g)$y);
         } else {
            layout <- layout_with_qfr(g,
               repulse=repulse,
               ...);
         }
      } else {
         if (verbose) {
            jamba::printDebug("spread_igraph_labels(): ",
               "Calling ","layout_with_qfr()"," for node coordinates.");
         }
         layout <- layout_with_qfr(g,
            repulse=repulse,
            verbose=verbose,
            ...);
      }
   } else if (is.function(layout)) {
      if (verbose) {
         jamba::printDebug("spread_igraph_labels(): ",
            "Calling ","layout()"," for node coordinates.");
      }
      layout <- layout(g);
   }

   if (length(rownames(layout)) == 0) {
      rownames(layout) <- igraph::V(g)$name;
   }
   if (do_reorder) {
      if (verbose) {
         jamba::printDebug("spread_igraph_labels(): ",
            "Calling multienrichjam::reorderIgraphNodes()");
         jamba::printDebug("spread_igraph_labels(): ",
            "head(layout) before:");
         print(head(layout));
      }
      # if sortAttributes is empty, use defaults from reorderIgraphNodes()
      if (length(sortAttributes) == 0) {
         sortAttributes <- eval(formals(reorderIgraphNodes)$sortAttributes);
      }
      # apply node re-ordering step
      g <- reorderIgraphNodes(g,
         layout=layout,
         nodeSortBy=nodeSortBy,
         sortAttributes=sortAttributes,
         verbose=verbose,
         ...);
      layout <- igraph::graph_attr(g, "layout");
      if (verbose) {
         jamba::printDebug("spread_igraph_labels(): ",
            "head(layout) after:");
         print(head(layout));
      }
   }
   g_angle <- jamba::nameVector(sapply(seq_len(igraph::vcount(g)), function(i){
      xy1 <- layout[i,1:2,drop=FALSE];
      xy2 <- layout[as.numeric(igraph::ego(g, nodes=i, mindist=1)[[1]]),1:2,drop=FALSE];
      if (length(xy2) == 0) {
         xy2 <- matrix(ncol=2, c(0,0));
      }
      xymean <- colMeans(xy1[rep(1, nrow(xy2)),,drop=FALSE] - xy2);
      -(xyAngle(xymean[1], xymean[2]*y_bias, directed=TRUE) + 0) %% 360
   }), igraph::V(g)$name);
   if (update_g_coords) {
      g <- igraph::set_graph_attr(g, "layout", layout);
   }
   igraph::V(g)$label.degree <- jamba::deg2rad(g_angle);
   if (!"label.dist" %in% igraph::list.vertex.attributes(g)) {
      igraph::V(g)$label.dist <- label_min_dist;
   } else {
      igraph::V(g)$label.dist <- pmax(igraph::V(g)$label.dist, label_min_dist);
   }
   g;
}


#' Subset igraph by connected components
#'
#' Subset igraph by connected components
#'
#' This function is intended to help drill down into an igraph
#' object that contains multiple connected components.
#'
#' By default, it sorts the components from largest number of nodes,
#' to smallest, which helps choose the largest connected component,
#' or subsequent components in size order.
#'
#' The components can also be filtered to require a minimum number
#' of connected nodes.
#'
#' At its core, this function is a wrapper to `igraph::components()`
#' and `igraph::subgraph()`.
#'
#' @family jam igraph functions
#'
#' @param g igraph object
#' @param keep numeric vector indicating which component or components
#'    to keep in the final output. When `order_by_size=TRUE`, components
#'    are ordered by size, from largest to smallest, in that case
#'    `keep=1` will return only the one largest connected subgraph.
#' @param min_size numeric value indicating the number of nodes required
#'    in all connected components returned. This filter is applied after
#'    the `keep` argument.
#' @param order_by_size logical indicating whether the connected components
#'    are sorted by size, largest to smallest, and therefore re-numbered.
#'    Otherwise, the components are somewhat randomly labeled based
#'    upon the output of `igraph::components()`.
#' @param ... additional arguments are passed to `igraph::components()`.
#'
#' @export
subset_igraph_components <- function
(g,
 keep=NULL,
 min_size=1,
 order_by_size=TRUE,
 ...)
{
   gc <- igraph::components(g,
      ...);
   vnum <- seq_len(igraph::vcount(g));
   gc_list <- split(vnum, membership(gc));
   if (order_by_size) {
      gc_order <- names(rev(sort(lengths(gc_list))));
      gc_list <- gc_list[gc_order];
      names(gc_list) <- seq_along(gc_list);
   }
   if (length(keep) > 0) {
      gc_list <- gc_list[names(gc_list) %in% as.character(keep)];
   }
   if (length(min_size) > 0) {
      gc_list <- gc_list[lengths(gc_list) >= min_size];
   }
   g <- subgraph_jam(g,
      v=sort(unlist(gc_list)));
   return(g);
}

#' Layout specification for Qgraph Fruchterman-Reingold
#'
#' @family jam igraph functions
#'
#' @export
with_qfr <- function (...,repulse=4) {
   layout_qfr <- function(graph,...){layout_with_qfr(graph,repulse=repulse,...)}
   igraph:::layout_spec(layout_qfr, ...)
}

#' Subgraph using Jam extended logic
#'
#' Subgraph using Jam extended logic
#'
#' This function extends the `igraph::subgraph()` function
#' to include proper subset of the graph attribute `"layout"`,
#' which for some unknown reason does not subset the layout
#' matrix consistent with the subset of `igraph` nodes.
#'
#' @family jam utility functions
#' @family jam igraph functions
#'
#' @param graph `igraph` object
#' @param v `integer` or `logical` vector indicating the nodes to
#'    retain in the final `igraph` object.
#'
#' @export
subgraph_jam <- function
(graph,
 v)
{
   if ("layout" %in% igraph::list.graph.attributes(graph)) {
      g_layout <- igraph::graph_attr(graph, "layout");
      if (any(c("numeric","matrix") %in% class(g_layout))) {
         if (ncol(g_layout) == 2) {
            g_layout_new <- g_layout[v,,drop=FALSE];
         } else if (ncol(g_layout) == 3) {
            g_layout_new <- g_layout[v,,,drop=FALSE];
         } else {
            stop("The layout matrix cannot contain more than 3 dimensions.");
         }
      } else {
         stop("The layout must be numeric matrix.");
      }
   }
   graph <- igraph::induced_subgraph(graph=graph,
      vids=v);
   if ("layout" %in% igraph::list.graph.attributes(graph)) {
      graph <- igraph::set_graph_attr(graph, "layout", g_layout_new);
   }
   return(graph);
}

#' Color igraph edges using node colors (deprecated)
#'
#' Color igraph edges using node colors (deprecated)
#'
#' Note: This function is deprecated in favor of
#' `color_edges_by_nodes()`.
#'
#' This function uses the average color for the two nodes
#' involved in each edge, and applies that as the new edge color.
#'
#' The color for each node depends upon the node shape, where
#' shape `"pie"` uses the average color from `"pie.color"`, and
#' shape `"coloredrectangle"` uses the avereage color from
#' `"coloredrect.color"`. Everything else uses `"color"`.
#'
#' This function relies upon `avg_colors_by_list()` to
#' blend multiple colors together.
#'
#' @param g `igraph` object
#' @param alpha `NULL` or numeric vector with value between 0 and 1,
#'    where 0 is transparent and 1 is non-transparent. When supplied,
#'    this value is passed to `jamba::alpha2col()` to apply alpha
#'    transparency to each edge color.
#' @param ... additional arguments are ignored.
#'
#' @family jam igraph functions
#'
#' @export
color_edges_by_nodes_deprecated <- function
(g,
 alpha=NULL,
 ...)
{
   edge_m <- igraph::as_edgelist(g, names=FALSE);
   g_colors <- ifelse(igraph::V(g)$shape %in% "circle",
      igraph::V(g)$color,
      ifelse(igraph::V(g)$shape %in% "pie",
         igraph::V(g)$pie.color,
         ifelse(igraph::V(g)$shape %in% "coloredrectangle",
            igraph::V(g)$coloredrect.color,
            "#FFFFFF00")));
   g_color <- avg_colors_by_list(g_colors);
   edge_m[] <- g_color[edge_m];
   edge_l <- as.list(data.frame(t(edge_m)));
   edge_colors <- avg_colors_by_list(edge_l);
   if (length(alpha) > 0) {
      edge_colors <- jamba::alpha2col(edge_colors,
         alpha=alpha);
   }
   igraph::E(g)$color <- unname(edge_colors);
   return(g);
}


#' Ordered colors from a list of color vectors
#'
#' Ordered colors from a list of color vectors
#'
#' This function takes a list of colors and returns the unique
#' order of colors based upon the order in vectors of the list.
#' It is mainly intended to be called by `reorderIgraphNodes()`,
#' however the function is useful for inferring the proper order
#' of unique colors from a list of various subsets of colors.
#'
#' The basic assumption is that there exists one true order of
#' unique colors, and that each vector in the list contains a
#' subset of those colors which is consistent with this
#' true order of colors.
#'
#' The function uses only vectors that contain two or more
#' colors, and therefore requires that all unique colors are
#' present in the subset of vectos in the list where length >= 2.
#' It then uses vectors with two or more colors, calculates
#' the average observed rank for each color, then uses that
#' average rank to define the overall color order.
#'
#' If not all unique colors are present in vectors with two or
#' more colors, the fallback sort uses `colorjam::sort_colors()`.
#'
#' @family jam list functions
#'
#' @return character vector of unique colors in `x`
#'
#' @param x list of character vectors that contain valid R colors.
#'
#' @export
colors_from_list <- function
(x,
 return_type=c("colors", "order"),
 verbose=FALSE,
 ...)
{
   return_type <- match.arg(return_type);
   if (!is.list(x)) {
      x <- as.list(x);
   }
   pcu <- unique(unlist(x));
   # pcu_names <- lapply(x, names);
   # names(pcu) <- pcu_names;
   pc2u <- unique(unlist(x[lengths(x) > 1]));
   if (all(pcu %in% pc2u)) {
      pc2 <- x[lengths(x) > 1];
      pcdf <- jamba::rbindList(lapply(pc2, function(pc2i){
         data.frame(color=as.character(pc2i),
            name=jamba::rmNULL(nullValue=NA, names(pc2i)),
            rank=seq_along(pc2i))
      }))
      # pcdf_u <- venndir::shrink_df(pcdf, by=c("color", "name"), num_func=mean);
      pcdf_mean <- sapply(split(pcdf$rank, pcdf$color), mean);
      pcdf_u <- subset(pcdf, !duplicated(color))
      pcdf_u$rank <- pcdf_mean[pcdf_u$color]
      pcdf_u[,c("H", "C", "L")] <- colorjam::colors_to_df(pcdf_u$color)[,c("H", "C", "L")];
      pcdf_u_sort <- jamba::mixedSortDF(pcdf_u,
         byCols=c("rank",
            "name",
            "H",
            "C",
            "L",
            "color"))
      colorV <- pcdf_u_sort$color;
      colorVnames <- pcdf_u_sort$name;
      if (all(!is.na(colorVnames))) {
         names(colorV) <- colorVnames;
      } else {
         names(colorV) <- seq_along(colorV);
      }
      if (verbose) {
         jamba::printDebug("colors_from_list(): ",
            "Derived colorV from node color values:",
            names(colorV),
            fgText=list("darkorange1", "dodgerblue", NA),
            bgText=list(NA, NA, colorV));
      }
   } else {
      pcdf1 <- jamba::rbindList(lapply(jamba::rmNULL(unique(x)), function(pc2i){
         if (all(pc2i %in% c(NA))) {
            return(NULL)
         }
         if (length(names(pc2i)) == 0) {
            names(pc2i) <- pc2i;
         }
         data.frame(color=pc2i, name=names(pc2i), rank=seq_along(pc2i))
      }));
      if (length(pcu) > 1) {
         colorV <- colorjam::sort_colors(pcu,
            byCols=c("H", "-C", "-L"));
      } else {
         colorV <- pcu;
      }
      colorVnames <- pcdf1$name[match(colorV, pcdf1$color)];
      if (all(!is.na(colorVnames))) {
         names(colorV) <- colorVnames;
      } else {
         names(colorV) <- seq_along(colorV);
      }
      if (verbose) {
         jamba::printDebug("colors_from_list(): ",
            "Derived colorV by using sort_colors():",
            names(colorV),
            fgText=list("darkorange1", "dodgerblue", NA),
            bgText=list(NA, NA, colorV));
      }
   }
   if ("colors" %in% return_type) {
      return(colorV);
   }
   colorattrdf <- data.frame(jamba::rbindList(x));
   for (cacol in colnames(colorattrdf)) {
      colorattrdf[[cacol]] <- factor(colorattrdf[[cacol]], levels=colorV);
   }
   jString <- jamba::pasteByRowOrdered(colorattrdf);
   return(order(jString));
}



#' Flip direction of igraph edges
#'
#' Flip direction of igraph edges
#'
#' This function simply flips the direction of igraph edges,
#' keeping all other node and edge attributes.
#'
#' Note that this function will flip the order of nodes for each
#' edge defined by `edge_idx`, regardless whether the `igraph`
#' itself is a directed graph.
#'
#' When `edge_idx` is provided as a `character` vector edge sequence,
#' any entries that do not match edges in `g` are ignored. A summary
#' table is printed when `verbose=TRUE`.
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object
#' @param edge_idx `integer` index of edges in the order they are stored
#'    in `igraph::E(g)`, or
#'    what igraph calls an "edge sequence" which is a character name
#'    for each node, defined as "node1|node2". For example "D|A" would
#'    define an edge from node name "D" to node name "A".
#'    When `verbose=TRUE` a summary table is printed out to show which
#'    edges were flipped.
#' @param verbose `logical` indicating whether to print verbose output.
#'    When `verbose=TRUE` a summary table is printed with these columns:
#'    * `edge_seq`: the input edge sequence, for example when `edge_idx`
#'    is provided as a `character` vector, the input vector is printed
#'    here.
#'    * `edge_seq_matched`: edge sequence that matched the `g` object.
#'    For example, when `edge_idx` input is a `character` vector, only
#'    the edges that match the `g` input are included here.
#'    * `edge_idx`: the integer index values of edges flipped.
#'    An `NA` value indicates the edge was not flipped, which should
#'    only happen when input `edge_idx` is provided as a `character`
#'    vector and some edges do not match the `g` input.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' am <- matrix(ncol=5, nrow=5, byrow=TRUE,
#'    data=c(0,0,0,0,0,
#'       1,0,0,0,0,
#'       1,0,0,0,0,
#'       1,0,0,0,0,
#'       1,0,0,0,0),
#'    dimnames=list(head(LETTERS, 5),
#'       head(LETTERS, 5)))
#' am;
#' g1 <- igraph::graph_from_adjacency_matrix(am)
#' plot(g1);
#' g2 <- flip_edges(g1, 3:4);
#' plot(g2);
#'
#' @export
flip_edges <- function
(g,
 edge_idx,
 verbose=FALSE,
 ...)
{
   #
   # validate edge_idx
   g_edge_ids <- igraph::as_ids(igraph::E(g));
   if (is.character(edge_idx)) {
      edge_idx_match <- match(edge_idx, g_edge_ids);
      edge_summary_df <- data.frame(
         edge_seq=edge_idx,
         edge_seq_matched=ifelse(is.na(edge_idx_match),
            "", edge_idx),
         edge_idx=edge_idx_match)
      edge_idx_seq <- edge_idx[!is.na(edge_idx_match)];
      edge_idx <- edge_idx_match[!is.na(edge_idx_match)];
   } else {
      edge_summary_df <- data.frame(
         edge_seq=g_edge_ids[edge_idx],
         edge_idx=edge_idx)
   }
   if (length(edge_idx) == 0) {
      jamba::printDebug("flip_edges(): ",
         "No edge_idx entries matched. Returning g.")
      return(g)
   }
   if (verbose) {
      jamba::printDebug("flip_edges(): ",
         "summary:");
      print(edge_summary_df);
   }

   edgeattrnames <- igraph::list.edge.attributes(g);
   edgeattrs <- lapply(jamba::nameVector(edgeattrnames), function(edgeattrname){
      igraph::edge_attr(g,
         name=edgeattrname,
         index=edge_idx)
   })
   add_el <- igraph::as_edgelist(g, names=FALSE)[edge_idx, , drop=FALSE]
   rm_edgenames <- attr(igraph::E(g)[edge_idx], "vnames")

   g2 <- igraph::add_edges(
      igraph::delete_edges(g, rm_edgenames),
      edges=as.numeric(t(add_el[,2:1, drop=FALSE])),
      attr=edgeattrs)
   return(g2)
}
