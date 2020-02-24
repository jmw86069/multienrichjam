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
 area=8*(vcount(g)^2),
 repulse.rad=(vcount(g)^repulse),
 constraints=NULL,
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
   set.seed(seed);

   ## Handle weights from E(g)$weight if supplied
   if (length(weights) == 0 && "weight" %in% list.edge.attributes(g)) {
      if (verbose) {
         jamba::printDebug("layout_with_qfr(): ",
            "Using E(g)$weight to define weights during layout.");
      }
      weights <- E(g)$weight;
   }

   frL <- qgraph::qgraph.layout.fruchtermanreingold(e,
      vcount=vcount(g),
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
 verbose=FALSE,
 ...)
{
   layout_xy <- layout_with_qfr(g=g,
      repulse=repulse,
      seed=seed,
      verbose=verbose,
      ...);
   if (verbose) {
      jamba::printDebug("relayout_with_qfr(): ",
         "head(layout_xy):");
      print(head(layout_xy));
   }
   g <- set_graph_attr(g,
      "layout",
      layout_xy);
   if (spread_labels) {
      g <- spread_igraph_labels(g,
         #layout=layout_xy,
         verbose=verbose,
         ...);
   }
   return(g);
}

#' Create a cnetplot igraph object
#'
#' Create a cnetplot igraph object
#'
#' The purpose of this function is to mimic the steps in `DOSE:::cnetplot()`
#' except not plot the output, and provide some customizations.
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
   V(g)$frame.color <- makeColorDarker(V(g)$color,
      darkFactor=1.5,
      alpha=0.5);
   V(g)$label.cex <- labelCex;

   ## Match category and gene color
   categoryColorMatch <- "#E5C494";
   geneColorMatch <- "#B3B3B3";
   iWhichCat <- which(V(g)$color %in% categoryColorMatch);
   iWhichGene <- which(V(g)$color %in% geneColorMatch);
   #V(g)$color <- ifelse(V(g)$color %in% categoryColorMatch,
   #   categoryColor,
   #   ifelse(V(g)$color %in% geneColorMatch,
   #      geneColor,
   #      V(g)$color));
   if (!all(categoryColor %in% categoryColorMatch) &&
         length(iWhichCat) > 0) {
      V(g)[iWhichCat]$color <- rep(categoryColor,
         length.out=length(iWhichCat));
   }
   if (!all(geneColor %in% geneColorMatch) &&
         length(iWhichGene) > 0) {
      V(g)[iWhichGene]$color <- rep(geneColor,
         length.out=length(iWhichGene));
   }

   ## Optionally custom color nodes using colorSub
   if (!is.null(colorSub)) {
      if (any(names(colorSub) %in% V(g)$name)) {
         iWhich <- which(V(g)$name %in% names(colorSub));
         if (length(iWhich) > 0) {
            V(g)[iWhich]$color <- colorSub[V(g)[iWhich]$name];
         }
      }
   }

   ## Normalize gene and category node sizes
   if (normalizeGeneSize && length(iWhichCat) && length(iWhichGene)) {
      geneSize <- mean(V(g)[iWhichGene]$size);
      catSizes <- V(g)[iWhichCat]$size;
      degreeCat <- degree(g)[iWhichCat];
      #catSize <- sqrt(degreeCat/pi)/2;
      catSize <- sqrt(degreeCat/pi);
      V(g)[iWhichCat]$size <- catSize;
      if (geneSize > median(catSize)) {
         if (verbose) {
            jamba::printDebug("cnetplotJam(): ",
               "Shrinking gene nodes to median category node size.");
         }
         geneSize <- median(catSize);
         V(g)[iWhichGene]$size <- geneSize;
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
   g <- enrichplot:::list2graph(inputList);
   #g <- DOSE::setting.graph.attributes(g);
   lengthOfCategory <- length(inputList);

   ## Color gene nodes by fold change if supplied
   if (!is.null(foldChange)) {
      node.idx <- (lengthOfCategory + 1):length(V(g));
      fcColors <- vals2colorLevels(foldChange,
         col=colorRamp,
         divergent=TRUE,
         ...);
      V(g)[node.idx]$color <- fcColors;
      g <- scaleNodeColor(g, foldChange, node.idx, DE.foldChange);
   }

   V(g)$size <- 5;
   V(g)$color <- geneColor;
   V(g)[1:lengthOfCategory]$size <- 30;
   V(g)[1:lengthOfCategory]$color <- categoryColor;

   ## Size category nodes
   if (is.numeric(categorySize)) {
      ## If supplied a numeric vector, size categories directly
      V(g)[1:lengthOfCategory]$size <- categorySize;
   } else {
      if (categorySize == "geneNum") {
         n <- degree(g)[1:lengthOfCategory];
         V(g)[1:lengthOfCategory]$size <- n/sum(n) * 100;
      } else if (categorySize == "pvalue") {
         if (is.null(pvalue) || any(is.na(pvalue))) {
            stop("pvalue must not be NULL or contain NA values.");
         }
         pScore <- -log10(pvalue);
         V(g)[1:lengthOfCategory]$size <- pScore/sum(pScore) * 100;
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
   iNodes <- which(toupper(V(g)$name) %in% toupper(rownames(valueIMcolors)));
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

   Vshapes <- V(g)$shape;
   if (length(Vshapes) == 0) {
      V(g)$shape <- "circle";
   }

   ## Change above logic to use lapply() then try to assign to igraph
   ## in batch steps
   if (verbose) {
      jamba::printDebug("igraph2pieGraph(): ",
         "iNodeParamsL");
   }
   iNodeParamsL <- lapply(iNodes, function(i){
      iNode <- V(g)$name[i];
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
   V(g)[iNodes]$shape <- "pie";
   V(g)[iNodes]$pie.value <- iPieValueL;
   V(g)[iNodes]$pie <- iPieValueL;
   V(g)[iNodes]$pie.color <- iPieColorL;
   V(g)[iNodes]$pie.names <- iPieNamesL;

   ## Define stable layout
   if (defineLayout) {
      if (verbose) {
         jamba::printDebug("igraph2pieGraph(): ",
            "layout_with_qfr(g, repulse=",
            repulse,
            ")");
      }
      set.seed(seed);
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
      if (is.null(V(g)$label)) {
         V(g)$label <- V(g)$name;
      }
      if (!is.null(maxNchar)) {
         V(g)$label <- substr(V(g)$label, 1, maxNchar);
      }
      V(g)$label <- jamba::ucfirst(tolower(V(g)$label));
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
#' @param x,y numeric coordinates, where x can be a two-column numeric
#'    matrix of x,y coordinates.
#' @param a,b numeric values indicating x- and y-axis radius, before
#'    rotation if `angle` is non-zero.
#' @param angle numeric value indicating the rotation of ellipse.
#' @param segment NULL or numeric vector of two values indicating the
#'    start and end angles for the ellipse, prior to rotation.
#' @param arc.only logical indicating whether to draw the ellipse
#'    arc without connecting to the center of the ellipse. Set
#'    `arc.only=FALSE` when segment does not include the full circle,
#'    to draw only the wedge.
#' @param nv the number of vertices around the center to draw.
#' @param deg logical indicating whether input `angle` and `segment`
#'    values are in degrees, or `deg=FALSE` for radians.
#' @param border,col,lty,lwd arguments passed to `graphics::polygon()`.
#' @param draw logical indicating whether to draw the ellipse.
#' @param ... additional arguments are passed to `graphics::polygon()`
#'    when `draw=TRUE`.
#'
#' @family jam igraph functions
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
   if (length(segment) == 0) {
      if (deg) {
         segment <- c(0, 360);
      } else {
         segment <- c(0, pi*2);
      }
   }
   ## if input is in degrees
   if (deg) {
      angle <- angle * pi/180;
      segment <- segment * pi/180;
   }
   z <- seq(from=segment[1],
      to=segment[2],
      length=nv + 1);
   xx <- a * cos(z);
   yy <- b * sin(z);
   alpha <- xyAngle(xx,
      yy,
      directed=TRUE,
      deg=FALSE);
   rad <- sqrt(xx^2 + yy^2)
   xp <- rad * cos(alpha + angle) + x;
   yp <- rad * sin(alpha + angle) + y;
   if (!arc.only) {
      xp <- c(x, xp, x);
      yp <- c(y, yp, y);
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
   df <- data.frame(nodeType=V(g)$nodeType,
      name=V(g)$name,
      label=V(g)$label,
      degree=degree(g),
      membership=components(g)$membership);
   if (getNeighbors || checkSubsets) {
      df$neighbors <- jamba::cPaste(
         lapply(seq_len(vcount(g)), function(i){
            neighbors(g, i)$name;
         }),
         doSort=TRUE,
         makeUnique=TRUE);
   }
   if (checkSubsets) {
      im <- cnet2im(df=df)
      ## determine if neighbors for a Set node are completely contained
      ## in another Set node
      imSet <- (t(im) %*% im);
      isSubset <- (rowSums(imSet >= rowMaxs(imSet)) > 1);
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
#' library(igraph);
#' library(multienrichjam);
#' g <- graph.full(n=3);
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
#' V(g)$size2 <- c(20,30,20);
#'
#' color_v <- rep("white", 21);
#' color_v[c(1,3,7,9,15,19,20,21)] <- colorjam::rainbowJam(5);
#' V(g)$coloredrect.color <- split(
#'    color_v,
#'    rep(V(g)$name, c(2,3,2)*3));
#' par("mfrow"=c(2,2));
#' lg <- layout_nicely(g);
#' plot(g, layout=lg);
#'
#' g2 <- removeIgraphBlanks(g, constrain="none");
#' V(g2)$size2 <- V(g2)$size2 / 3;
#' plot(g2, layout=lg, main="constrain='none'");
#'
#' g3 <- removeIgraphBlanks(g, constrain="nrow");
#' plot(g3, layout=lg, main="constrain='nrow'");
#'
#' g4 <- removeIgraphBlanks(g, constrain="ncol");
#' plot(g4, layout=lg, main="constrain='ncol'");
#'
#' #
#' g7 <- graph.full(n=7);
#' V(g7)$coloredrect.color <- lapply(c(1,2,3,4,2,3,4),
#'    function(i){colorjam::rainbowJam(i)});
#' V(g7)$coloredrect.ncol <- c(1,1,1,1,2,3,4);
#' V(g7)$coloredrect.nrow <- c(1,2,3,4,1,1,1);
#' V(g7)$coloredrect.names <- V(g7)$coloredrect.color;
#' V(g7)$shape <- "coloredrectangle";
#' V(g7)$size <- 10;
#' V(g7)$size2 <- V(g7)$coloredrect.ncol * 10;
#' lg7 <- layout_nicely(g7);
#' par("mfrow"=c(2,2));
#' plot(g7, layout=lg7, vertez.size2=10);
#' plot(g7, layout=lg7, vertex.size2=V(g7)$coloredrect.ncol*10);
#' plot(g7, layout=lg7, vertex.size2=V(g7)$coloredrect.nrow*10);
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
 pie_to_circle=TRUE,
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
   ixV <- which(lengths(V(g)$coloredrect.color) > 0);

   if ("coloredrect.color" %in% list.vertex.attributes(g)) {
      if (verbose) {
         jamba::printDebug("removeIgraphBlanks(): ",
            "Adjusting coloredrect nodes.");
      }
      ## Rewrote code to use vectorized logic.

      ## Determine the coloredrect.ncol to use in resizing
      ncolVbefore <- get.vertex.attribute(g, "coloredrect.ncol");
      if (length(ncolVbefore) == 0) {
         ncolVbefore <- lengths(get.vertex.attribute(g, "coloredrect.color"));
      }
      nrowVbefore <- get.vertex.attribute(g, "coloredrect.nrow");

      ## Determine which pie wedges are blank
      iCrColorL <- igraph::get.vertex.attribute(g, "coloredrect.color");
      crBlanksL <- isColorBlank(iCrColorL,
         blankColor=blankColor,
         c_max=c_max,
         l_min=l_min,
         alpha_max=alpha_max);

      ## Check for all blanks
      all_blank <- lapply(crBlanksL, all);
      if (any(all_blank)) {
         #
      }

      crLengths <- lengths(iCrColorL);
      crSplitV <- rep(factor(seq_len(vcount(g))), crLengths);
      ## Vector of TRUE,FALSE
      crBlanksV <- unlist(unname(crBlanksL));
      ## Iterate each attribute
      crAttrs <- intersect(c("coloredrect.color", "coloredrect.names"),
         list.vertex.attributes(g));
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
         nrowV <- rmNULL(get.vertex.attribute(g, "coloredrect.nrow"),
            nullValue=1);
         ncolV <- rmNULL(get.vertex.attribute(g, "coloredrect.ncol"),
            nullValue=1);
         byrowV <- rmNULL(get.vertex.attribute(g, "coloredrect.byrow")*1,
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
            ########################################
            ## constrain "none"
            crL <- unname(split(unlist(crAttrL)[!crBlanksV],
               crSplitV[!crBlanksV]));
            crLengthsNew <- lengths(crL);
            crChanged <- (crLengths != crLengthsNew);
            ncolV <- ifelse(nrowV == 1 | ncolV > 1, crLengthsNew, ncolV);
            nrowV <- ifelse(nrowV == 1 | ncolV > 1, 1, crLengthsNew);

            g <- set.vertex.attribute(g,
               name="coloredrect.nrow",
               value=nrowV);
            g <- set.vertex.attribute(g,
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
                  iMblank <- (colMins(iM) == 1);
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
               g <- set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.ncol",
                  value=iMncol);
               g <- set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.nrow",
                  value=iMnrow);
               g <- set.vertex.attribute(g,
                  index=iSet,
                  name="coloredrect.color",
                  value=iMvalsL);
               g <- set.vertex.attribute(g,
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
         ncolVafter <- get.vertex.attribute(g, "coloredrect.ncol");
         nrowVafter <- get.vertex.attribute(g, "coloredrect.nrow");
         resizeWhich <- (ncolVbefore != ncolVafter) |  (nrowVbefore != nrowVafter);
         if (any(resizeWhich)) {
            new_size2 <-nrowVbefore / nrowVafter *
               get.vertex.attribute(g,
                  name="size2");
            if (verbose) {
               print(data.frame(ncolVbefore,
                  ncolVafter,
                  size2=V(g)$size2,
                  new_size2));
            }
            g <- set.vertex.attribute(g,
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
         list.vertex.attributes(g));

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
         pieSplitV <- rep(seq_len(vcount(g)), pieLengths);
         ## Vector of TRUE,FALSE
         pieBlanksV <- unlist(unname(pieBlanksL));
         ## Iterate each pie attribute
         for (pieAttr in pieAttrs) {
            pieAttrL <- igraph::get.vertex.attribute(g, pieAttr);
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
               ## TODO: only update nodes that change
               g <- igraph::set.vertex.attribute(g,
                  name=pieAttr,
                  value=pieL);
            }
         }
      }
   }
   if (pie_to_circle) {
      is_pie <- V(g)$shape %in% "pie";
      if (any(is_pie)) {
         is_single_color <- lengths(V(g)[is_pie]$pie.color) == 1;
         if (any(is_single_color)) {
            switch_nodes <- which(is_pie)[is_single_color];
            i_colors <- unname(unlist(V(g)[switch_nodes]$pie.color));
            V(g)[switch_nodes]$color <- i_colors;
            V(g)[switch_nodes]$shape <- "circle";
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
 whichNodes=seq_len(vcount(g)),
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
      which(V(g)$shape %in% "pie"));

   if (length(whichPie) == 0) {
      return(g);
   }

   V(g)[whichPie]$coloredrect.color <- V(g)[whichPie]$pie.color;
   V(g)[whichPie]$coloredrect.names <- V(g)[whichPie]$pie.names;
   V(g)[whichPie]$coloredrect.nrow <- nrow;
   V(g)[whichPie]$coloredrect.ncol <- ncol;
   V(g)[whichPie]$coloredrect.byrow <- byrow;
   V(g)[whichPie]$shape <- "coloredrectangle";
   V(g)[whichPie]$size2 <- sqrt(V(g)[whichPie]$size)*1.5;
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
#' @param g igraph object
#' @param sortAttributes character vector of node attribute
#'    names, to be applied in order when sorting nodes.
#' @param nodeSortBy character vector containing `"x"` and
#'    `"y"` indicating the primary axis used to sort nodes.
#' @param layout numeric matrix of node coordinates, or
#'    function used to produce layout coordinates. When layout
#'    is `NULL`, this function tries to use graph attribute
#'    `"layout"`, otherwise
#'    the `relayout_with_qfr` is called.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
reorderIgraphNodes <- function
(g,
 sortAttributes=c("pie.color","coloredrect.color", "color", "label", "name"),
 nodeSortBy=c("x","y"),
 layout=NULL,
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
   if (length(layout) == 0) {
      if (!"layout" %in% list.graph.attributes(g)) {
         if (all(c("x","y") %in% vertex_attr_names(g))) {
            layout <- cbind(V(g)$x, V(g)$y);
            g <- set_graph_attr(g, "layout", layout)
         } else {
            g <- relayout_with_qfr(g, ...);
         }
      }
   } else if (is.function(layout)) {
      g2 <- layout(g);
      if ("igraph" %in% class(g2)) {
         g <- g2;
      } else {
         g <- set_graph_attr(g, "layout", g2);
      }
   } else if (jamba::igrepHas("data.*frame|matrix", class(layout))) {
      g <- set_graph_attr(g, "layout", layout);
   } else {
      stop("reorderIgraphNodes() could not determine the layout from the given input.");
   }

   ## comma-delimited neighboring nodes for each node
   neighborG <- jamba::cPaste(
      lapply(seq_len(vcount(g)), function(v){
         neighbors(g, v, mode="all");
      }),
      doSort=TRUE,
      makeUnique=FALSE);
   names(neighborG) <- seq_len(vcount(g));

   ## Determine which edge groups are present multiple times
   neighborGct <- jamba::tcount(neighborG, minCount=2);
   if (length(neighborGct) == 0) {
      if (verbose) {
         jamba::printDebug("reorderIgraphNodes(): ",
            "found no edge groups, returning input graph unchanged.");
      }
      return(g);
   }

   ## Use one or more vertex attributes for grouping
   sortAttributes <- intersect(sortAttributes, vertex_attr_names(g));
   if (length(sortAttributes) == 0) {
      if (verbose) {
         jamba::printDebug("reorderIgraphNodes(): ",
            "did not find any matching sortAttributes in the igraph object.");
      }
      return(g);
   }
   neighborA <- pasteByRow(do.call(cbind, lapply(sortAttributes,
      function(sortAttribute){
         j <- (get.vertex.attribute(g, sortAttribute));
         names(j) <- seq_len(vcount(g));

         if (sortAttribute %in% c("coloredrect.color","pie.color","color")) {
            j_colors <- vertex_attr(g, sortAttribute);
            j_colors_v <- avg_colors_by_list(j_colors);
            j_sorted <- sort_colors(j_colors_v);
            j_rank <- match(j_colors_v, unique(j_sorted));
            jString <- factor(j_sorted,
               levels=unique(j_sorted));
            if (1 == 1) {
               jString <- sapply(seq_along(j), function(j1){
                  if ("coloredrect.nrow" %in% list.vertex.attributes(g)) {
                     cnrow <- V(g)[j1]$coloredrect.nrow;
                     cbyrow <- jamba::rmNA(rmNULL=TRUE,
                        naValue=TRUE,
                        nullValue=TRUE,
                        V(g)[j1]$coloredrect.byrow);
                  } else {
                     cnrow <- 1;
                     cbyrow <- TRUE;
                  }
                  j2 <- matrix(
                     data=j[[j1]],
                     nrow=cnrow,
                     byrow=cbyrow);
                  j3 <- matrix(
                     data=1*(!isColorBlank(j2, ...)),
                     nrow=cnrow,
                     byrow=cbyrow);
                  j3blendByRow <- tryCatch({
                     apply(j2, 1, function(i){avg_colors_by_list(list(i))});
                  }, error=function(e) {
                     ""
                  });
                  if (verbose && j1 == 1) {
                     jamba::printDebug("head(rowSums(j3)):");print(head(rowSums(j3), 10));
                     jamba::printDebug("head(j3blendByRow):");print(head(j3blendByRow, 10));
                     jamba::printDebug("head(j2):");print(head(j2, 10));
                     jamba::printDebug("head(j3):");print(head(j3, 10));
                  }
                  paste(c(
                     j_rank[j1],
                     -rowSums(j3),
                     j3blendByRow,
                     pasteByRow(j3, sep=""),
                     pasteByRow(j2)),
                     collapse="_")
               });
            }
            names(jString) <- seq_len(vcount(g));
            if (verbose) {
               jamba::printDebug("reorderIgraphNodes(): ",
                  "head(jString):",
                  head(jString));
            }
            jString;
         } else {
            if (jamba::igrepHas("list", class(j))) {
               jString <- jamba::cPaste(j);
            } else if (jamba::igrepHas("numeric|integer|float", class(j))) {
               jString <- round(j, digits=2);
            } else {
               jString <- j;
            }
            names(jString) <- seq_len(vcount(g));
            jString;
         }
      }
   )), sep="_");
   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborA):");
      print(head(neighborA));
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(names(neighborG)):");
      print(head(names(neighborG)));
   }

   ## data.frame with the attribute sort, and the node sort
   neighborDF <- data.frame(
      vertex=names(neighborG),
      #vertex=names(neighborG),
      edgeGroup=neighborG,
      sortAttribute=neighborA[names(neighborG)],
      x=g$layout[,1],
      y=g$layout[,2]);

   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "head(neighborDF):");
      print(head(neighborDF));
   }

   ## The following code iterates each edge group and reassigns
   ## layout coordinates from left-to-right
   ## Bonus points: figure out the coordinates to assign using
   ## left-to-right logic, but then assign those coordinates
   ## top-to-bottom.
   if (verbose) {
      jamba::printDebug("reorderIgraphNodes(): ",
         "nodeSortBy:",
         nodeSortBy);
   }
   newDF <- jamba::rbindList(lapply(names(neighborGct), function(Gname){
      iDF <- subset(neighborDF, edgeGroup %in% Gname);
      xyOrder <- jamba::mixedSortDF(iDF,
         byCols=nodeSortBy);

      nodeOrder <- jamba::mixedSortDF(iDF,
         byCols=match(c("sortAttribute","vertex"), colnames(iDF)));

      nodeOrder[,c("x","y")] <- xyOrder[,c("x","y")];
      ## If there are repeated sortAttributes, we use them to place subsets
      ## of nodes top to bottom within each group of coordinates
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
      nodeOrder;
   }));
   iMatch <- match(newDF$vertex, neighborDF$vertex);
   neighborDF[iMatch,c("x","y")] <- newDF[,c("x","y")];
   subDF <- subset(neighborDF, edgeGroup %in% neighborDF[1,"edgeGroup"]);
   new_layout <- as.matrix(neighborDF[,c("x","y"),drop=FALSE]);
   g <- set_graph_attr(g, "layout", new_layout);
   return(g);
}

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
   keep_nodes <- (degree(g) >= min_degree);
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
#'    coordinates will be stored in `V(g)$x` and `V(g)$y`.
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
 sortAttributes=c("pie.color", "coloredrect.color", "color", "label", "name"),
 nodeSortBy=c("x","y"),
 repulse=3.5,
 force_relayout=FALSE,
 label_min_dist=0.5,
 verbose=FALSE,
 ...)
{
   ##
   if (length(layout) == 0) {
      if (!force_relayout) {
         if ("layout" %in% list.graph.attributes(g)) {
            if (verbose) {
               jamba::printDebug("spread_igraph_labels(): ",
                  "Using ","layout"," from graph attributes.");
            }
            layout <- g$layout;
         } else if (all(c("x", "y") %in% list.vertex.attributes(g))) {
            if (verbose) {
               jamba::printDebug("spread_igraph_labels(): ",
                  "Using ","x,y"," from vertex attributes.");
            }
            layout <- cbind(x=V(g)$x, y=V(g)$y);
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
      rownames(layout) <- V(g)$name;
   }
   if (do_reorder) {
      if (verbose) {
         jamba::printDebug("spread_igraph_labels(): ",
            "Calling multienrichjam::reorderIgraphNodes()");
         jamba::printDebug("spread_igraph_labels(): ",
            "head(layout) before:");
         print(head(layout));
      }
      g <- reorderIgraphNodes(g,
         layout=layout,
         nodeSortBy=nodeSortBy,
         sortAttributes=sortAttributes,
         verbose=verbose);
      layout <- g$layout;
      if (verbose) {
         jamba::printDebug("spread_igraph_labels(): ",
            "head(layout) after:");
         print(head(layout));
      }
   }
   g_angle <- jamba::nameVector(sapply(seq_len(vcount(g)), function(i){
      xy1 <- layout[i,1:2,drop=FALSE];
      xy2 <- layout[as.numeric(ego(g, nodes=i, mindist=1)[[1]]),1:2,drop=FALSE];
      if (length(xy2) == 0) {
         xy2 <- matrix(ncol=2, c(0,0));
      }
      xymean <- colMeans(xy1[rep(1, nrow(xy2)),,drop=FALSE] - xy2);
      -(xyAngle(xymean[1], xymean[2]*y_bias, directed=TRUE) + 0) %% 360
   }), V(g)$name);
   if (update_g_coords) {
      g <- set_graph_attr(g, "layout", layout);
   }
   V(g)$label.degree <- deg2rad(g_angle);
   if (!"label.dist" %in% list.vertex.attributes(g)) {
      V(g)$label.dist <- label_min_dist;
   } else {
      V(g)$label.dist <- pmax(V(g)$label.dist, label_min_dist);
   }
   g;
}

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
 enrichIM=NULL,
 enrichIMcolors=NULL,
 coloredrect_nrow=1,
 coloredrect_ncol=NULL,
 coloredrect_byrow=TRUE,
 ...)
{
   categoryShape <- match.arg(categoryShape);
   geneShape <- match.arg(geneShape);
   geneLabelColor <- head(geneLabelColor, 1);
   categoryLabelColor <- head(categoryLabelColor, 1);

   ## Accept memIM as list mem output from multiEnrichMap()
   if (is.list(memIM) &&
         all(c("memIM", "enrichIM", "geneIM", "geneIMcolors", "enrichIMcolors") %in% names(memIM))) {
      if (length(geneIM) == 0) {
         geneIM <- memIM[["geneIM"]];
      }
      if (length(geneIMcolors) == 0) {
         geneIMcolors <- memIM[["geneIMcolors"]];
      }
      if (length(enrichIM) == 0) {
         enrichIM <- memIM[["enrichIM"]];
      }
      if (length(enrichIMcolors) == 0) {
         enrichIMcolors <- memIM[["enrichIMcolors"]];
      }
      memIM <- memIM[["memIM"]];
   }

   ## Convert to igraph
   g <- igraph::graph_from_incidence_matrix(memIM != 0);

   ## Adjust aesthetics
   V(g)$nodeType <- "Set";
   V(g)$nodeType[match(rownames(memIM), V(g)$name)] <- "Gene";
   table(V(g)$nodeType);
   isset <- V(g)$nodeType %in% "Set";
   V(g)$size <- ifelse(isset, categorySize, geneSize);
   V(g)$label.cex <- ifelse(isset, categoryCex, geneCex);
   V(g)$label.color <- ifelse(isset, categoryLabelColor, geneLabelColor);
   V(g)$color <- ifelse(isset, categoryColor, geneColor);
   V(g)$frame.color <- jamba::makeColorDarker(V(g)$color,
      darkFactor=frame_darkFactor,
      ...);

   ## Optionally apply gene node coloring
   if (length(geneIM) > 0 && length(geneIMcolors) > 0) {
      geneIM <- subset(geneIM, rownames(geneIM) %in% V(g)$name[!isset]);
      gene_match <- match(rownames(geneIM),
         V(g)$name[!isset]);
      gene_which <- which(!isset)[gene_match];
      for (j in c("pie", "pie.value", "pie.color", "pie.names",
         "coloredrect.color", "coloredrect.nrow",
         "coloredrect.byrow", "coloredrect.value",
         "coloredrect.names")) {
         if (!j %in% list.vertex.attributes(g)) {
            if (igrepHas("color$", j)) {
               vertex_attr(g, j) <- as.list(V(g)$color);
            } else if (igrepHas("byrow", j)) {
               vertex_attr(g, j) <- rep(unlist(coloredrect_byrow), length.out=vcount(g));
            } else if (igrepHas("ncol", j) && length(coloredrect_ncol) > 0) {
               vertex_attr(g, j) <- rep(unlist(coloredrect_ncol), length.out=vcount(g));
            } else if (igrepHas("nrow", j) && length(coloredrect_nrow) > 0) {
               vertex_attr(g, j) <- rep(unlist(coloredrect_nrow), length.out=vcount(g));
            } else {
               vertex_attr(g, j) <- as.list(rep(1, vcount(g)));
            }
         }
      }
      V(g)$pie[gene_which] <- lapply(gene_which, function(i){
         rep(1, ncol(geneIM))
      });
      V(g)$pie.value[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         geneIM[i,,drop=FALSE]
      });
      V(g)$pie.color[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         geneIMcolors[i,]
      });
      V(g)$pie.names[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         colnames(geneIM);
      });
      ## Now do the same for coloredrectangle node shapes
      V(g)$coloredrect.color[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         geneIMcolors[i,]
      });
      V(g)$coloredrect.value[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         geneIM[i,,drop=FALSE]
      });
      if (length(coloredrect_byrow) > 0) {
         V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_byrow, length.out=length(gene_which));
      }
      if (length(coloredrect_nrow) > 0) {
         V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_nrow, length.out=length(gene_which));
      }
      if (length(coloredrect_ncol) > 0) {
         V(g)$coloredrect.byrow[gene_which] <- rep(coloredrect_ncol, length.out=length(gene_which));
      }
      V(g)$coloredrect.names[gene_which] <- lapply(nameVector(rownames(geneIM)), function(i){
         colnames(geneIM);
      });
      V(g)$shape[gene_which] <- geneShape;
   }
   ## Optionally apply category/set node coloring
   if (length(enrichIM) > 0 && length(enrichIMcolors) > 0) {
      enrichIM <- subset(enrichIM, rownames(enrichIM) %in% V(g)$name[isset]);
      enrich_match <- match(rownames(enrichIM),
         V(g)$name[isset]);
      enrich_which <- which(isset)[enrich_match];
      for (j in c("pie", "pie.value", "pie.color",
         "coloredrect.color", "coloredrect.nrow",
         "coloredrect.byrow", "coloredrect.value")) {
         if (!j %in% list.vertex.attributes(g)) {
            if (igrepHas("color$", j)) {
               vertex_attr(g, j) <- as.list(V(g)$color);
            } else if (igrepHas("byrow", j)) {
               vertex_attr(g, j) <- as.list(rep(coloredrect_byrow, length.out=vcount(g)));
            } else if (igrepHas("ncol", j) && length(coloredrect_ncol) > 0) {
               vertex_attr(g, j) <- as.list(rep(coloredrect_ncol, length.out=vcount(g)));
            } else if (igrepHas("nrow", j) && length(coloredrect_nrow) > 0) {
               vertex_attr(g, j) <- as.list(rep(coloredrect_nrow, length.out=vcount(g)));
            } else {
               vertex_attr(g, j) <- as.list(rep(1, vcount(g)));
            }
         }
      }
      V(g)$pie[enrich_which] <- lapply(enrich_which, function(i){
         rep(1, ncol(enrichIM))
      });
      V(g)$pie.value[enrich_which] <- lapply(nameVector(rownames(enrichIM)), function(i){
         enrichIM[i,,drop=FALSE]
      });
      V(g)$pie.color[enrich_which] <- lapply(nameVector(rownames(enrichIM)), function(i){
         enrichIMcolors[i,]
      });
      ## Now do the same for coloredrectangle node shapes
      V(g)$coloredrect.color[enrich_which] <- lapply(nameVector(rownames(enrichIM)), function(i){
         enrichIMcolors[i,]
      });
      V(g)$coloredrect.value[enrich_which] <- lapply(nameVector(rownames(enrichIM)), function(i){
         enrichIM[i,,drop=FALSE]
      });
      if (length(coloredrect_byrow) > 0) {
         V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_byrow, length.out=length(enrich_which));
      }
      if (length(coloredrect_nrow) > 0) {
         V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_nrow, length.out=length(enrich_which));
      }
      if (length(coloredrect_ncol) > 0) {
         V(g)$coloredrect.byrow[enrich_which] <- rep(coloredrect_ncol, length.out=length(enrich_which));
      }
      V(g)$shape[enrich_which] <- categoryShape;
   }
   return(g);
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
   vnum <- seq_len(vcount(g));
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
#' @param v integer or logical vector indicating the nodes to
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

#' Color igraph edges using node colors
#'
#' Color igraph edges using node colors
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
#' @export
color_edges_by_nodes <- function
(g,
 alpha=NULL,
 ...)
{
   edge_m <- as_edgelist(g, names=FALSE);
   g_colors <- ifelse(V(g)$shape %in% "circle",
      V(g)$color,
      ifelse(V(g)$shape %in% "pie",
         V(g)$pie.color,
         ifelse(V(g)$shape %in% "coloredrectangle",
            V(g)$coloredrect.color,
            "#FFFFFF00")));
   g_color <- avg_colors_by_list(g_colors);
   edge_m[] <- g_color[edge_m];
   edge_l <- as.list(data.frame(t(edge_m)));
   edge_colors <- avg_colors_by_list(edge_l);
   if (length(alpha) > 0) {
      edge_colors <- alpha2col(edge_colors,
         alpha=alpha);
   }
   E(g)$color <- unname(edge_colors);
   return(g);
}

#' Jam igraph vectorized plot function
#'
#' Jam igraph vectorized plot function
#'
#' This function is a complete copy of `igraph:::plot.igraph()` with
#' two changes to enable vectorized plotting in these situations:
#'
#' 1. When there are multiple vertex `"shape"` attributes, the
#' base plot function draws each individual `igraph` vertex one by one.
#' The `jam_plot_igraph()` plots nodes of each shape in a group. For
#' fairly large `igraph` objects (about 1000 nodes), this change
#' is substantially faster.
#' 2. When there are multiple font families, the default plot function
#' draws each label one by one. The `jam_plot_igraph()` draws
#' labels in groups of font family. This situation is very rare, however
#' when used the speed improvement is substantial.
#'
#' Note that this function is not called by default, and is only called
#' by `multienrichjam::jam_igraph()`.
#'
#' Smaller changes: `mark.col` and `mark.border` default colors
#' now call `colorjam::rainbowJam()` instead of `grDevices::rainbow()`.
#'
#' All arguments are documented in `igraph::plot.igraph()`.
#'
#' @family jam igraph functions
#'
#' @inheritParams igraph::plot.igraph
#' @param pie_to_jampie logical indicating whether to convert
#'    vertex shape `"pie"` to `"jampie"` in order to use vectorized
#'    plotting.
#'
#' @export
jam_plot_igraph <- function
(x,
 axes = FALSE,
 add = FALSE,
 xlim = c(-1, 1),
 ylim = c(-1, 1),
 mark.groups = list(),
 mark.shape = 1/2,
 mark.col = colorjam::rainbowJam(length(mark.groups), alpha = 0.3),
 mark.border = colorjam::rainbowJam(length(mark.groups), alpha = 1),
 mark.expand = 15,
 pie_to_jampie=TRUE,
...)
{
   graph <- x
   if (!igraph::is_igraph(graph)) {
      stop("Not a graph object")
   }
   params <- igraph:::i.parse.plot.params(graph, list(...))
   vertex.size <- 1/200 * params("vertex", "size")
   label.family <- params("vertex", "label.family")
   label.font <- params("vertex", "label.font")
   label.cex <- params("vertex", "label.cex")
   label.degree <- params("vertex", "label.degree")
   label.color <- params("vertex", "label.color")
   label.dist <- params("vertex", "label.dist")
   labels <- params("vertex", "label")

   shape <- igraph:::igraph.check.shapes(params("vertex", "shape"));
   ## Optionally convert shape "pie" to "jampie" for vectorized plotting
   if (pie_to_jampie) {
      shape <- ifelse(shape %in% "pie",
         "jampie",
         shape);
   }

   edge.color <- params("edge", "color")
   edge.width <- params("edge", "width")
   edge.lty <- params("edge", "lty")
   arrow.mode <- params("edge", "arrow.mode")
   edge.labels <- params("edge", "label")
   loop.angle <- params("edge", "loop.angle")
   edge.label.font <- params("edge", "label.font")
   edge.label.family <- params("edge", "label.family")
   edge.label.cex <- params("edge", "label.cex")
   edge.label.color <- params("edge", "label.color")
   elab.x <- params("edge", "label.x")
   elab.y <- params("edge", "label.y")
   arrow.size <- params("edge", "arrow.size")[1]
   arrow.width <- params("edge", "arrow.width")[1]
   curved <- params("edge", "curved")
   if (is.function(curved)) {
      curved <- curved(graph)
   }
   layout <- params("plot", "layout")
   margin <- params("plot", "margin")
   margin <- rep(margin, length = 4)
   rescale <- params("plot", "rescale")
   asp <- params("plot", "asp")
   frame <- params("plot", "frame")
   main <- params("plot", "main")
   sub <- params("plot", "sub")
   xlab <- params("plot", "xlab")
   ylab <- params("plot", "ylab")
   palette <- params("plot", "palette")
   if (!is.null(palette)) {
      old_palette <- palette(palette)
      on.exit(palette(old_palette), add = TRUE)
   }
   arrow.mode <- igraph:::i.get.arrow.mode(graph, arrow.mode)
   maxv <- max(vertex.size)
   if (rescale) {
      layout <- igraph::norm_coords(layout, -1, 1, -1, 1)
      xlim <- c(xlim[1] - margin[2] - maxv,
         xlim[2] + margin[4] + maxv);
      ylim <- c(ylim[1] - margin[1] - maxv,
         ylim[2] + margin[3] + maxv)
   }
   if (!add) {
      plot(0,
         0,
         type = "n",
         xlab = xlab,
         ylab = ylab,
         xlim = xlim,
         ylim = ylim,
         axes = axes,
         frame = frame,
         asp = asp,
         main = main,
         sub = sub)
   }
   if (!is.list(mark.groups) && is.numeric(mark.groups)) {
      mark.groups <- list(mark.groups)
   }
   mark.shape <- rep(mark.shape, length = length(mark.groups))
   mark.border <- rep(mark.border, length = length(mark.groups))
   mark.col <- rep(mark.col, length = length(mark.groups))
   mark.expand <- rep(mark.expand, length = length(mark.groups))
   for (g in seq_along(mark.groups)) {
      v <- igraph::V(graph)[mark.groups[[g]]]
      if (length(vertex.size) == 1) {
         vs <- vertex.size
      } else {
         vs <- rep(vertex.size, length = igraph::vcount(graph))[v]
      }
      igraph:::igraph.polygon(layout[v, , drop = FALSE],
         vertex.size = vs,
         expand.by = mark.expand[g]/200,
         shape = mark.shape[g],
         col = mark.col[g],
         border = mark.border[g])
   }
   el <- igraph::as_edgelist(graph, names = FALSE)
   loops.e <- which(el[, 1] == el[, 2])
   nonloops.e <- which(el[, 1] != el[, 2])
   loops.v <- el[, 1][loops.e]
   loop.labels <- edge.labels[loops.e]
   loop.labx <- if (is.null(elab.x)) {
      rep(NA, length(loops.e))
   } else {
      elab.x[loops.e]
   }
   loop.laby <- if (is.null(elab.y)) {
      rep(NA, length(loops.e))
   } else {
      elab.y[loops.e]
   }
   edge.labels <- edge.labels[nonloops.e]
   elab.x <- if (is.null(elab.x)) {
      NULL
   } else {
      elab.x[nonloops.e]
   }
   elab.y <- if (is.null(elab.y)){
      NULL
   } else {
      elab.y[nonloops.e]
   }
   el <- el[nonloops.e, , drop = FALSE]
   edge.coords <- matrix(0, nrow = nrow(el), ncol = 4)
   edge.coords[, 1] <- layout[, 1][el[, 1]]
   edge.coords[, 2] <- layout[, 2][el[, 1]]
   edge.coords[, 3] <- layout[, 1][el[, 2]]
   edge.coords[, 4] <- layout[, 2][el[, 2]]
   if (length(unique(shape)) == 1) {
      ec <- igraph:::.igraph.shapes[[shape[1]]]$clip(edge.coords,
         el,
         params = params,
         end = "both")
   } else {
      shape <- rep(shape, length = vcount(graph))
      ec <- edge.coords
      ec[, 1:2] <- t(sapply(seq(length = nrow(el)), function(x) {
         igraph:::.igraph.shapes[[shape[el[x, 1]]]]$clip(
            edge.coords[x, , drop = FALSE],
            el[x, , drop = FALSE],
            params = params,
            end = "from")
      }))
      ec[, 3:4] <- t(sapply(seq(length = nrow(el)), function(x) {
         igraph:::.igraph.shapes[[shape[el[x, 2]]]]$clip(
            edge.coords[x, , drop = FALSE],
            el[x, , drop = FALSE],
            params = params,
            end = "to")
      }))
   }
   x0 <- ec[, 1]
   y0 <- ec[, 2]
   x1 <- ec[, 3]
   y1 <- ec[, 4]
   if (length(loops.e) > 0) {
      ec <- edge.color
      if (length(ec) > 1) {
         ec <- ec[loops.e]
      }
      point.on.cubic.bezier <- function(cp, t) {
         c <- 3 * (cp[2, ] - cp[1, ])
         b <- 3 * (cp[3, ] - cp[2, ]) - c
         a <- cp[4, ] - cp[1, ] - c - b
         t2 <- t * t
         t3 <- t * t * t
         a * t3 + b * t2 + c * t + cp[1, ]
      }
      compute.bezier <- function(cp, points) {
         dt <- seq(0, 1, by = 1/(points - 1))
         sapply(dt, function(t) {
            point.on.cubic.bezier(cp, t)
         })
      }
      plot.bezier <- function
      (cp,
         points,
         color,
         width,
         arr,
         lty,
         arrow.size,
         arr.w)
      {
         p <- compute.bezier(cp, points)
         polygon(p[1, ],
            p[2, ],
            border = color,
            lwd = width,
            lty = lty)
         if (arr == 1 || arr == 3) {
            igraph:::igraph.Arrows(p[1, ncol(p) - 1],
               p[2, ncol(p) - 1],
               p[1, ncol(p)],
               p[2, ncol(p)],
               sh.col = color,
               h.col = color,
               size = arrow.size,
               sh.lwd = width,
               h.lwd = width,
               open = FALSE,
               code = 2,
               width = arr.w)
         }
         if (arr == 2 || arr == 3) {
            igraph:::igraph.Arrows(p[1, 2],
               p[2, 2],
               p[1, 1], p[2, 1],
               sh.col = color,
               h.col = color,
               size = arrow.size,
               sh.lwd = width,
               h.lwd = width,
               open = FALSE,
               code = 2,
               width = arr.w)
         }
      }
      loop <- function
      (x0,
         y0,
         cx = x0,
         cy = y0,
         color,
         angle = 0,
         label = NA,
         width = 1,
         arr = 2,
         lty = 1,
         arrow.size = arrow.size,
         arr.w = arr.w,
         lab.x, lab.y)
      {
         rad <- angle
         center <- c(cx, cy)
         cp <- matrix(
            c(x0,
               y0,
               x0 + 0.4,
               y0 + 0.2,
               x0 + 0.4,
               y0 - 0.2,
               x0,
               y0),
            ncol = 2,
            byrow = TRUE)
         phi <- atan2(cp[, 2] - center[2],
            cp[, 1] - center[1])
         r <- sqrt((cp[, 1] - center[1])^2 + (cp[, 2] - center[2])^2)
         phi <- phi + rad
         cp[, 1] <- cx + r * cos(phi)
         cp[, 2] <- cy + r * sin(phi)
         plot.bezier(cp,
            50,
            color,
            width,
            arr = arr,
            lty = lty,
            arrow.size = arrow.size,
            arr.w = arr.w)
         if (is.language(label) || !is.na(label)) {
            lx <- x0 + 0.3
            ly <- y0
            phi <- atan2(ly - center[2],
               lx - center[1])
            r <- sqrt((lx - center[1])^2 + (ly - center[2])^2)
            phi <- phi + rad
            lx <- cx + r * cos(phi)
            ly <- cy + r * sin(phi)
            if (!is.na(lab.x)) {
               lx <- lab.x
            }
            if (!is.na(lab.y)) {
               ly <- lab.y
            }
            text(lx,
               ly,
               label,
               col = edge.label.color,
               font = edge.label.font,
               family = edge.label.family,
               cex = edge.label.cex)
         }
      }
      ec <- edge.color
      if (length(ec) > 1) {
         ec <- ec[loops.e]
      }
      vs <- vertex.size
      if (length(vertex.size) > 1) {
         vs <- vs[loops.v]
      }
      ew <- edge.width
      if (length(edge.width) > 1) {
         ew <- ew[loops.e]
      }
      la <- loop.angle
      if (length(loop.angle) > 1) {
         la <- la[loops.e]
      }
      lty <- edge.lty
      if (length(edge.lty) > 1) {
         lty <- lty[loops.e]
      }
      arr <- arrow.mode
      if (length(arrow.mode) > 1) {
         arr <- arrow.mode[loops.e]
      }
      asize <- arrow.size
      if (length(arrow.size) > 1) {
         asize <- arrow.size[loops.e]
      }
      xx0 <- layout[loops.v, 1] + cos(la) * vs
      yy0 <- layout[loops.v, 2] - sin(la) * vs
      mapply(loop,
         xx0,
         yy0,
         color = ec,
         angle = -la,
         label = loop.labels,
         lty = lty,
         width = ew,
         arr = arr,
         arrow.size = asize,
         arr.w = arrow.width,
         lab.x = loop.labx,
         lab.y = loop.laby)
   }
   if (length(x0) != 0) {
      if (length(edge.color) > 1) {
         edge.color <- edge.color[nonloops.e]
      }
      if (length(edge.width) > 1) {
         edge.width <- edge.width[nonloops.e]
      }
      if (length(edge.lty) > 1) {
         edge.lty <- edge.lty[nonloops.e]
      }
      if (length(arrow.mode) > 1) {
         arrow.mode <- arrow.mode[nonloops.e]
      }
      if (length(arrow.size) > 1) {
         arrow.size <- arrow.size[nonloops.e]
      }
      if (length(curved) > 1) {
         curved <- curved[nonloops.e]
      }
      if (length(unique(arrow.mode)) == 1) {
         lc <- igraph:::igraph.Arrows(x0,
            y0,
            x1,
            y1,
            h.col = edge.color,
            sh.col = edge.color,
            sh.lwd = edge.width,
            h.lwd = 1,
            open = FALSE,
            code = arrow.mode[1],
            sh.lty = edge.lty,
            h.lty = 1,
            size = arrow.size,
            width = arrow.width,
            curved = curved)
         lc.x <- lc$lab.x
         lc.y <- lc$lab.y
      } else {
         curved <- rep(curved, length = ecount(graph))[nonloops.e]
         lc.x <- lc.y <- numeric(length(curved))
         for (code in 0:3) {
            valid <- arrow.mode == code
            if (!any(valid)) {
               next
            }
            ec <- edge.color
            if (length(ec) > 1) {
               ec <- ec[valid]
            }
            ew <- edge.width
            if (length(ew) > 1) {
               ew <- ew[valid]
            }
            el <- edge.lty
            if (length(el) > 1) {
               el <- el[valid]
            }
            lc <- igraph:::igraph.Arrows(x0[valid],
               y0[valid],
               x1[valid],
               y1[valid],
               code = code,
               sh.col = ec,
               h.col = ec,
               sh.lwd = ew,
               h.lwd = 1,
               h.lty = 1,
               sh.lty = el,
               open = FALSE,
               size = arrow.size,
               width = arrow.width,
               curved = curved[valid])
            lc.x[valid] <- lc$lab.x
            lc.y[valid] <- lc$lab.y
         }
      }
      if (!is.null(elab.x)) {
         lc.x <- ifelse(is.na(elab.x),
            lc.x,
            elab.x)
      }
      if (!is.null(elab.y)) {
         lc.y <- ifelse(is.na(elab.y),
            lc.y,
            elab.y)
      }
      text(lc.x,
         lc.y,
         labels = edge.labels,
         col = edge.label.color,
         family = edge.label.family,
         font = edge.label.font,
         cex = edge.label.cex)
   }
   rm(x0, y0, x1, y1)
   if (length(unique(shape)) == 1) {
      igraph:::.igraph.shapes[[shape[1]]]$plot(
         layout,
         params = params)
   } else {
      ## Loop each unique shape here instead of individual nodes
      unique_shapes <- unique(shape);
      nodes_by_shape <- split(seq_len(vcount(graph)), shape);
      if (1 == 1) {
         sapply(nodes_by_shape, function(x){
            shape1 <- shape[x[1]];
            igraph:::.igraph.shapes[[shape1]]$plot(
               layout[x, , drop = FALSE],
               v = x,
               params = params)
         })
      } else {
         sapply(seq(length = vcount(graph)), function(x) {
            igraph:::.igraph.shapes[[shape[x]]]$plot(
               layout[x, , drop = FALSE],
               v = x,
               params = params)
         })
      }
   }
   par(xpd = TRUE)
   x <- layout[, 1] +
      (label.dist *
            cos(-label.degree) *
            (vertex.size + 6 * 8 * log10(2))/200);
   y <- layout[, 2] +
      (label.dist *
            sin(-label.degree) *
            (vertex.size + 6 * 8 * log10(2))/200);
   if (length(label.family) == 1) {
      text(x,
         y,
         labels = labels,
         col = label.color,
         family = label.family,
         font = label.font,
         cex = label.cex)
   } else {
      if1 <- function(vect, idx) if (length(vect) == 1) {
         vect
      } else {
         vect[idx]
      }
      ## Also loop through each family here instead of each node
      unique_family <- unique(label.family);
      label.family <- rep(label.family, length.out=vcount(graph));
      nodes_by_family <- split(seq_len(vcount(graph)), label.family);
      if (1 == 1) {
         sapply(nodes_by_family, function(v){
            family1 <- label.family[v[1]];
            text(x[v],
               y[v],
               labels = if1(labels, v),
               col = if1(label.color, v),
               family = family1,
               font = if1(label.font, v),
               cex = if1(label.cex, v))
         })
      } else {
         sapply(seq_len(vcount(graph)), function(v) {
            text(x[v],
               y[v],
               labels = if1(labels, v),
               col = if1(label.color, v),
               family = if1(label.family, v),
               font = if1(label.font, v),
               cex = if1(label.cex, v))
         })
      }
   }
   rm(x, y)
   invisible(NULL)
}

