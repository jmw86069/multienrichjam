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
#'    \code{\link{qgraph::qgraph.layout.fruchtermanreingold}}.
#'
#' @return two-column numeric matrix with coordinates for each vertex.
#'
#' @seealso `qgraph::qgraph.layout.fruchtermanreingold()`
#'
#' @family jam igraph functions
#'
#' @examples
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
#'
#' @export
layout_with_qfr <- function
(g,
 repulse=3.1,
 area=8*(vcount(g)^2),
 repulse.rad=(vcount(g)^repulse),
 constraints=NULL,
 seed=123,
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
   frL <- qgraph::qgraph.layout.fruchtermanreingold(e,
      vcount=vcount(g),
      area=area,
      repulse.rad=repulse.rad,
      constraints=constraints,
      ...);
   frL;
}

#' igraph layout function using qgraph Fruchterman-Reingold
#'
#' igraph layout function using qgraph Fruchterman-Reingold
#'
#' @export
layout_with_qfrf <- function
(repulse=3.1,
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
      printDebug("cnetplotJam(): ",
         "dColname:",
         dColname);
   }
   keepColnames <- intersect(c(nodeLabel, "pvalue"),
      names(y));
   if (is(x, "enrichResult") || is(x, "gseaResult")) {
      gc <- DOSE::geneInCategory(x);
      names(gc) <- y[[dColname]];
      if (verbose) {
         printDebug("cnetplotJam(): ",
            "assigned gc <- DOSE::geneInCategory(x), length(gc):",
            length(gc),
            ", head(gc, 2):");
         print(head(gc, 2));
      }
   } else {
      stop("x should be an 'enrichResult' or 'gseaResult' object...")
   }
   if (verbose) {
      printDebug("cnetplotJam(): ",
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
         printDebug("Shrinking gene nodes to median category node size.");
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
         printDebug("igraph2pieGraph(): ",
            "No nodes were changed, returning input graph.");
      }
      return(g);
   }
   if (verbose) {
      printDebug("igraph2pieGraph(): ",
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
      printDebug("igraph2pieGraph(): ",
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
         printDebug("igraph2pieGraph(): ",
            "layout_with_qfr(g, repulse=",
            repulse,
            ")");
      }
      set.seed(seed);
      layoutG <- layout_with_qfr(g, repulse=repulse);
      V(g)$x <- layoutG[,1];
      V(g)$y <- layoutG[,2];
   }

   ## Optionally update labels for maximum characters
   if (updateLabels) {
      if (verbose) {
         printDebug("igraph2pieGraph(): ",
            "Updating node labels.");
      }
      if (is.null(V(g)$label)) {
         V(g)$label <- V(g)$name;
      }
      if (!is.null(maxNchar)) {
         V(g)$label <- substr(V(g)$label, 1, maxNchar);
      }
      V(g)$label <- ucfirst(tolower(V(g)$label));
   }
   return(g);
}

#' plot function for igraph vertex shape coloredrectangle
#'
#' plot function for igraph vertex shape coloredrectangle
#'
#' This function defines the plotting function for custom igraph vertex
#' shape coloredrectangle. The coloredrectangle shape is described as:
#'
#' \itemize{
#'    \item{a vertex drawn as a rectangle, filled with squares which are
#'   each individually colored.}
#'   \item{the squares are arrayed into a number of columns and rows,
#'   which are defined for each vertex using `coloredrect.ncol` and
#'   `coloredrect.nrow`, respectively.}
#'   \item{the vector of colors is arrayed as values of a matrix, therefore
#'   `coloredrect.byrow` is logical to indicate whether colors should fill
#'   the matrix by row, similar to how `byrow=TRUE` is used.}
#'   \item{the colors for each vertex are defined by
#'   `coloredrect.color`. When this value does not exist, this function
#'   will attempt to use values from `pie.color` if they exist.}
#'   \item{the size of the rectangle is defined by `size2`, which is
#'   typically the height of a `rectangle` node. The intent was to use
#'   a size parameter which is already convenient, but which does not
#'   conflict with the size used for vertex shapes such as `"circle"`.}
#'   \item{when a vertex has 3 columns, and 2 rows, and only 5 colors,
#'   the colors are not cycled to fill the complete node. Instead, the
#'   last color is `"transparent"` to render no color in that square.}
#'   \item{the frame around each node is described with `frame.color` and
#'   `frame.lwd` which draws one rectangular border around the full
#'   vertex. The border around all coloredrect vertices is drawn before the
#'   squares are drawn for all coloredrect vertices, so that the border will
#'   not overlap the squares. It has the visible effect of showing vertex
#'   borders behind the vertex colors, while allowing for vectorized drawing,
#'   which is substantially more efficient than per-vertex rendering.
#'   The `graphics::symbols()` function is used to render rectangles,
#'   it accepts only one value for `lwd` and `lty`, so rectangles
#'   are split by `lwd` and `lty` and rendered in order. As a
#'   result, if each vertex frame had a different lwd,lty combination, the
#'   rendering would effectively be non-vectorized. For large numbers of
#'   vertices and distinct lwd values, it would be preferred to reduce the
#'   lwd values to limited significant digits (e.g. with
#'   `signif(...,digits=2)`). That said, line width is not an optimal
#'   way to convey a quantitative measurement.}
#'   \item{a colored border can be drawn around each square in the vertex,
#'   defined by `coloredrect.border` and `coloredrect.lwd`,
#'   however it is recommended to use
#'   this color only as a visible break, for example the default is
#'   `"grey30"` which draws a simple border. Note that squares are
#'   rendered after the vertex frame color, so each square color border will
#'   be drawn over (on top of) the frame border.}
#'   \item{when `coloredrect.ncol` does not exist, it will attempt
#'   to use two rows of colors if `coloredrect.color` contains
#'   two or more values.}
#'   \item{when `coloredrect.nrow` does not exist, it will use a value
#'   based upon `coloredrect.ncol` and the length of
#'   `coloredrect.color`.}
#' }
#'
#' Currently, the plotting of coloredrectangle vertices does not define
#' clipping, which means edges are drawn to the center of each vertex,
#' and the coloredrectangle vertex shapes are drawn on top of the edges.
#' Transparent nodes will therefore show edges beneath them.
#'
#' @return invisible `list` of `data.frame` objects which were
#'    used to draw the rectangle objects.
#'    However the purpose of this function is the by-product that it
#'    draws rectangles onto an igraph graph.
#'
#' @seealso `igraph::shapes()`
#'
#' @family jam igraph functions
#'
#' @param coords two-column numeric matrix with x- and y-coordinates,
#'    respectively.
#' @param v optional ids of the vertices to plot. It should match
#'    the number of rows in the `coords` argument.
#' @param params a function object that can be called to query
#'    vertex/edge/plot graphical parameters. The first argument of
#'    the function is `'vertex'`, `'edge'` or `'plot'` to decide
#'    the type of the parameter, the second is a character string
#'    giving the name of the parameter.
#'
#' @export
shape.coloredrectangle.plot <- function
(coords,
 v=NULL,
 params)
{
   ## Purpose is to extend igraph:::.igraph.shape.rectangle.plot
   ## to handle the custom colored rectangle node type.
   ## Shape "coloredrectangle" has these features for each vertex:
   ## * the vertex is a rectangle filled with squares, where each
   ##   square has its own assigned color.
   ## * the squares are arrayed into a number of columns defined
   ##   by "coloredrect.ncol", and a number of rows defined by
   ##   "coloredrect.nrow". Each vertex can have its own number
   ##   of columns and rows, and its own number of colors.
   ## * the size of the overall rectangle is defined by "size2",
   ##   which is normally the height of a rectangle shape. The
   ##   intention was to allow using "size" for shapes such as
   ##   "circle" but "size2" for
   ##
   ## This function combines some logic from
   ## igraph:::.igraph.shape.pie.plot() as well.

   getparam <- function(pname) {
      p <- params("vertex", pname)
      if (length(p) != 1 && !is.null(v)) {
         p <- p[v]
      }
      p;
   }
   vertex.color <- getparam("color");

   vertex.frame.color <- getparam("frame.color");
   if (length(vertex.frame.color) == 0) {
      vertex.frame.color <- "grey30";
   }
   vertex.frame.lwd <- getparam("frame.lwd");
   if (length(vertex.frame.lwd) == 0) {
      vertex.frame.lwd <- 1;
   }

   vertex.coloredrect.color <- getparam("coloredrect.color");
   if (length(vertex.coloredrect.color) == 0) {
      vertex.coloredrect.color <- getparam("pie.color");
   }

   vertex.coloredrect.border <- getparam("coloredrect.border");
   if (length(vertex.coloredrect.border) == 0) {
      vertex.coloredrect.border <- "grey30";
      #vertex.coloredrect.border <- vertex.frame.color;
   }

   vertex.coloredrect.byrow <- getparam("coloredrect.byrow");
   if (length(vertex.coloredrect.byrow) == 0) {
      vertex.coloredrect.byrow <- TRUE;
   }

   vertex.coloredrect.ncol <- getparam("coloredrect.ncol");
   if (length(vertex.coloredrect.ncol) == 0) {
      vertex.coloredrect.ncol <- ceiling(length(vertex.coloredrect)/2);
   }
   vertex.coloredrect.nrow <- getparam("coloredrect.nrow");
   if (length(vertex.coloredrect.nrow) == 0) {
      vertex.coloredrect.nrow <- ceiling(
         length(vertex.coloredrect) /
            vertex.coloredrect.ncol);
   }
   vertex.coloredrect.lty <- getparam("coloredrect.lty");
   if (length(vertex.coloredrect.lty) == 0) {
      vertex.coloredrect.lty <- 1;
   }
   vertex.coloredrect.lwd <- getparam("coloredrect.lwd");
   if (length(vertex.coloredrect.lwd) == 0) {
      vertex.coloredrect.lwd <- 1;
   }

   vertex.size1 <- getparam("size2");
   vertex.size2 <- getparam("size");
   #printDebug("vertex.size1:", head(vertex.size1, 10),
   #   ", vertex.size2:", head(vertex.size2, 10));
   vertex.size1 <- rep(1/200 * vertex.size1, length.out=nrow(coords));
   vertex.size2 <- rep(1/200 * vertex.size2, length.out=nrow(coords));

   ## Use size1 (height) to define the size of each square, then apply that
   ## to calculate size2 (width)
   vertex.size1 <- vertex.size1 * 5 / vertex.coloredrect.nrow;
   vertex.size2 <- (vertex.size1 *
         (vertex.coloredrect.nrow / vertex.coloredrect.ncol));

   #printDebug("vertex.size1:", head(vertex.size1, 10),
   #   ", vertex.size2:", head(vertex.size2, 10));
   vertex.size <- cbind(vertex.size1, vertex.size2);

   ## Define custom function to help vectorize drawing, by creating a
   ## data.frame of coordinates for each square and rectangle
   mycoloredrectangle <- function
   (x,
      y,
      size1,
      size2,
      col,
      ncol,
      nrow,
      border,
      lty,
      lwd=1,
      frame.lwd=1,
      frame.color="grey30",
      byrow=TRUE,
      ...)
   {
      ## Purpose is to draw a rectangle filles with multi-color squares
      nrow <- rep(nrow, length.out=length(x));
      ncol <- rep(ncol, length.out=length(x));
      border <- rep(border, length.out=length(x));
      lty <- rep(lty, length.out=length(x));
      lwd <- rep(lwd, length.out=length(x));
      frame.lwd <- rep(frame.lwd, length.out=length(x));
      frame.color <- rep(frame.color, length.out=length(x));
      byrow <- rep(byrow, length.out=length(x));

      ## Iterate each vertex, create a data.frame describing
      ## frame and square colors, then combine into one large
      ## data.frame for vectorized drawing.
      rectDF <- rbindList(lapply(seq_along(x), function(k){
         xk <- x[[k]];
         yk <- y[[k]];
         size1k <- size1[[k]];
         size2k <- size2[[k]];
         nrowk <- nrow[[k]];
         ncolk <- ncol[[k]];
         byrowk <- byrow[[k]];
         ## Individual rectangles
         x01 <- seq(from=xk-(size1k/2),
            to=xk+(size1k/2),
            length.out=ncolk+1);
         y01 <- seq(from=yk-(size2k/2),
            to=yk+(size2k/2),
            length.out=nrowk+1);

         numk <- ncolk * nrowk;
         size1v <- rep(min(diff(x01)), numk);
         size2v <- rep(min(diff(y01)), numk);
         if (byrowk) {
            x01v <- rep(head(x01, -1), nrowk)+(size1v/2);
            y01v <- rep(rev(tail(y01, -1)), each=ncolk)-(size2v/2);
         } else {
            x01v <- rep(head(x01, -1), each=nrowk)+(size1v/2);
            y01v <- rep(rev(tail(y01, -1)), ncolk)-(size2v/2);
         }
         ## Make one large rectangle border
         x01all <- mean(x01v);
         y01all <- mean(y01v);
         ##
         colk <- rep(col[[k]], length.out=length(x01v));
         borderk <- rep(border[[k]], length.out=length(colk));
         ltyk <- rep(lty[[k]], length.out=length(colk));
         lwdk <- rep(lwd[[k]], length.out=length(colk));
         if (length(getOption("debug"))>0) {
            printDebug("k:", k,
               ", numk:", numk,
               ", byrowk:", byrowk,
               ", ncolk:", ncolk,
               ", nrowk:", nrowk,
               ", xk:", xk,
               ", x01v:", signif(digits=3, x01v),
               ", x01:", signif(digits=3, x01),
               ", yk:", yk,
               ", y01v:", signif(digits=3, y01v),
               ", y01:", signif(digits=3, y01),
               ", colk:", colk,
               ", size1v:", signif(digits=3, size1v),
               ", size2v:", signif(digits=3, size2v)
            );
         }
         kDF <- data.frame(x=c(xk, x01v),
            y=c(yk, y01v),
            bg=c("transparent", colk),
            fg=c(head(frame.color[[k]],1),
               borderk),
            rectx=c(head(size1v,1)*ncolk, size1v),
            recty=c(head(size2v,1)*nrowk, size2v),
            lty=c(head(ltyk,1), ltyk),
            lwd=c(head(frame.lwd[[k]],1)*2.5,
               lwdk/4),
            rect_type=rep(factor(c("frame","square")),
               c(1,length(borderk)))
         );
         kDF;
      }));

      ## Split into a list of data.frames, because symbols()
      ## can only use one value for lwd and lty.
      rectDFL <- split(rectDF,
         pasteByRowOrdered(rectDF[,c("rect_type","lwd","lty")]));

      for (rectDFi in rev(rectDFL)) {
         symbols(x=rectDFi$x,
            y=rectDFi$y,
            bg=rectDFi$bg,
            fg=rectDFi$fg,
            rectangles=as.matrix(rectDFi[,c("rectx","recty")]),
            add=TRUE,
            inches=FALSE,
            lty=rectDFi$lty,
            lwd=rectDFi$lwd);
      }
      return(rectDFL);
   }

   ## For reference, the code below draws a single rectangle
   if (1 == 2) {
      symbols(x=coords[, 1],
         y=coords[, 2],
         bg=vertex.color,
         fg=vertex.frame.color,
         rectangles=2*vertex.size,
         add=TRUE,
         inches=FALSE);
   }

   mcr <- mycoloredrectangle(x=coords[,1],
      y=coords[,2],
      size1=vertex.size1,
      size2=vertex.size2,
      col=vertex.coloredrect.color,
      ncol=vertex.coloredrect.ncol,
      nrow=vertex.coloredrect.nrow,
      border=vertex.coloredrect.border,
      lty=vertex.coloredrect.lty,
      lwd=vertex.coloredrect.lwd,
      frame.lwd=vertex.frame.lwd,
      frame.color=vertex.frame.color,
      byrow=vertex.coloredrect.byrow);

}

#' plot function for igraph vertex shape ellipse
#'
#' plot function for igraph vertex shape ellipse
#'
#' This function defines the plotting function for custom igraph vertex
#' shape ellipse.
#'
#' @export
shape.ellipse.plot <- function
(coords,
 v=NULL,
 params)
{
   vertex.color <- params("vertex", "color");
   if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v];
   }
   vertex.frame.color <- params("vertex", "frame.color");
   if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v];
   }
   vertex.frame.width <- params("vertex", "frame.width");
   if (length(vertex.frame.width) == 0) {
      vertex.frame.width <- 1;
   }
   if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v];
   }
   vertex.size <- 1/200 * params("vertex", "size");
   if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v];
   }
   vertex.ellipse.ratio <- params("vertex", "ellipse.ratio");
   if (length(vertex.ellipse.ratio) == 0) {
      vertex.ellipse.ratio <- 2;
   }
   if (length(vertex.ellipse.ratio) != 1 && !is.null(v)) {
      vertex.ellipse.ratio <- vertex.ellipse.ratio[v];
   }

   drawEllipse(x=coords[,1],
      y=coords[,2],
      a=vertex.size,
      b=vertex.size/vertex.ellipse.ratio,
      col=vertex.color,
      border=vertex.frame.color,
      lwd=vertex.frame.width,
      draw=TRUE);
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
