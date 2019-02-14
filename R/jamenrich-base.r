# jamenrich-base.r

.onLoad <- function
(libname,
 pkgname)
{
   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("coloredrectangle",
      clip=igraph::shape_noclip,
      plot=shape.coloredrectangle.plot)

}

#' plot function for igraph vertex shape coloredrectangle
#'
#' plot function for igraph vertex shape coloredrectangle
#'
#' This function defines the plotting function for custom igraph vertex
#' shape coloredrectangle. The coloredrectangle shape is described as:
#'
#' \itemize{
#'    \item a vertex drawn as a rectangle, filled with squares which are
#'   each individually colored.
#'   \item the squares are arrayed into a number of columns and rows,
#'   which are defined for each vertex using \code{coloredrect.ncol} and
#'   \code{coloredrect.nrow}, respectively.
#'   \item the vector of colors is arrayed as values of a matrix, therefore
#'   \code{coloredrect.byrow} is boolean to indicate whether colors should fill
#'   the matrix by row, similar to how \code{byrow=TRUE} is used.
#'   \item the colors for each vertex are defined by
#'   \code{coloredrect.color}. When this value does not exist, this function
#'   will attempt to use values from \code{pie.color} if they exist.
#'   \item the size of the rectangle is defined by \code{size2}, which is
#'   typically the height of a \code{rectangle} node. The intent was to use
#'   a size parameter which is already convenient, but which does not
#'   conflict with the size used for vertex shapes such as \code{'circle'}.
#'   \item when a vertex has 3 columns, and 2 rows, and only 5 colors,
#'   the colors are not cycled to fill the complete node. Instead, the
#'   last color is \code{"transparent"} to render no color in that square.
#'   \item the frame around each node is described with \code{frame.color} and
#'   \code{frame.lwd} which draws one rectangular border around the full
#'   vertex. The border around all coloredrect vertices is drawn before the
#'   squares are drawn for all coloredrect vertices, so that the border will
#'   not overlap the squares. It has the visible effect of showing vertex
#'   borders behind the vertex colors, while allowing for vectorized drawing,
#'   which is substantially more efficient than per-vertex rendering.
#'   The \code{\link{symbols}} function is used to render rectangles,
#'   it accepts only one value for \code{lwd} and \code{lty}, so rectangles
#'   are split by \code{lwd} and \code{lty} and rendered in order. As a
#'   result, if each vertex frame had a different lwd,lty combination, the
#'   rendering would effectively be non-vectorized. For large numbers of
#'   vertices and distinct lwd values, it would be preferred to reduce the
#'   lwd values to limited significant digits (e.g. with
#'   \code{signif(...,digits=2)}). That said, line width is not an optimal
#'   way to convey a quantitative measurement.
#'   \item a colored border can be drawn around each square in the vertex,
#'   defined by \code{coloredrect.border} and \code{coloredrect.lwd},
#'   however it is recommended to use
#'   this color only as a visible break, for example the default is
#'   \code{"grey30"} which draws a simple border. Note that squares are
#'   rendered after the vertex frame color, so each square color border will
#'   be drawn over (on top of) the frame border.
#'   \item when \code{coloredrect.ncol} does not exist, it will attempt
#'   to use two rows of colors if \code{coloredrect.color} contains
#'   two or more values.
#'   \item when \code{coloredrect.nrow} does not exist, it will use a value
#'   based upon \code{coloredrect.ncol} and the length of
#'   \code{coloredrect.color}.
#' }
#'
#' Currently, the plotting of coloredrectangle vertices does not define
#' clipping, which means edges are drawn to the center of each vertex,
#' and the coloredrectangle vertex shapes are drawn on top of the edges.
#' Transparent nodes will therefore show edges beneath them.
#'
#' @return invisible \code{list} of \code{data.frame} objects which were
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
#'    the number of rows in the \code{coords} argument.
#' @param params a function object that can be called to query
#'    vertex/edge/plot graphical parameters. The first argument of
#'    the function is \code{'vertex'}, \code{'edge'} or \code{'plot'} to decide
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
      vertex.coloredrect.nrow <- ceiling(length(vertex.coloredrect)/vertex.coloredrect.ncol);
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

#' Convert data.frame to enrichResult
#'
#' Convert data.frame to enrichResult
#'
#' This function takes a data.frame containing gene set enrichment
#' results, and converts it to a proper `enrichResult` object,
#' compatible with methods provided by the `clusterProfiler`
#' package.
#'
#' @param enrichDF data.frame representing gene set enrichment
#'    results.
#' @param msigdbGmtT optional gmtT object, a representation of
#'    `arules::transactions-class`.
#' @param pvalueCutoff numeric value range 0 to 1, to define the
#'    P-value threshold for enrichment results to be considered
#'    in downstream processing.
#' @param pAdjustMethod character string to define the P-value
#'    adjustment method, or `"none"` for no additional adjustment.
#'    See `stats::p.adjust()` for valid values.
#' @param keyColname character value of the `colname(enrichDF)`
#'    containing the unique row identifier. It can be a pathway_ID
#'    or any uniquely identifying value.
#' @param geneColname character value of the `colname(enrichDF)`
#'    containing delimiited genes in each pathway.
#' @param pathGenes character or value of the `colname(enrichDF)`
#'    containing the number of genes in each pathway. This value will be
#'    derived from `geneRatioColname` if needed.
#' @param geneHits character value of the `colname(enrichDF)`
#'    containing the integer count of the gene hits in each pathway.
#'    This value will be derived from `geneRatioColname` if needed.
#' @param geneRatioColname character value of the `colname(enrichDF)`
#'    containing the character ratio of gene hits to pathway size,
#'    in format "50/100". This value is used when either `"pathGenes"`
#'    or `"geneHits"` are not supplied.
#' @param geneDelim regular expression pattern used to separate
#'    genes in the `pathGenes` column into a vector of character
#'    values.
#' @param pvalueColname character value of the `colname(enrichDF)`
#'    containing enrichment P-values to use in downstream processing.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @family jam conversion functions
#'
#' @export
enrichDF2enrichResult <- function
(enrichDF=NULL,
 msigdbGmtT=NULL,#msigdbGmtTv50mouse,
 pvalueCutoff=0.15,
 pAdjustMethod="none",
 keyColname="itemsetID",
 pathGenes="pathGenes",
 geneColname="geneNames",
 geneHits="geneHits",
 geneRatioColname="GeneRatio",
 geneDelim="[,/ ]+",
 pvalueColname="P-Value",
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert an enrichment data.frame
   ## into enrichResult class format usable by clusterProfiler
   ## methods, like enrichMap()
   if (suppressPackageStartupMessages(!require(clusterProfiler))) {
      stop("enrichDF2enrichResult() requires the clusterProfiler package.");
   }
   enrichDF2 <- renameColumn(enrichDF,
      from=c(keyColname, pvalueColname, geneColname),
      to=c("ID", "pvalue", "geneID"));
   enrichDF2[,"p.adjust"] <- enrichDF2[,"pvalue"];

   ## Convert gene delimiters all to "/"
   enrichDF2[,"geneID"] <- gsub(geneDelim,
      "/",
      enrichDF2[,"geneID"]);

   ## Validate input colnames
   keyColname <- intersect(keyColname, colnames(enrichDF));
   pathGenes <- intersect(pathGenes, colnames(enrichDF));
   geneColname <- intersect(geneColname, colnames(enrichDF));
   geneHits <- intersect(geneHits, colnames(enrichDF));
   geneRatioColname <- intersect(geneRatioColname, colnames(enrichDF));
   pvalueColname <- intersect(pvalueColname, colnames(enrichDF));
   if (verbose) {
      printDebug("enrichDF2enrichResult(): ",
         "keyColname:", keyColname,
         "\npathGenes:", pathGenes,
         "\ngeneColname:", geneColname,
         "\ngeneHits:", geneHits,
         "\ngeneRatioColname:", geneRatioColname,
         "\npvalueColname:", pvalueColname);
   }

   if (length(geneRatioColname) > 0) {
      if (length(geneHits) == 0 ||
            geneHits == geneRatioColname) {
         if (verbose) {
            jamba::printDebug("enrichDF2enrichResult(): ",
               "deriving ",
               "geneHits",
               " from ",
               "GeneRatio");
         }
         geneHits <- "geneHits";
         enrichDF2[,geneHits] <- as.numeric(gsub("[/].*$", "",
            enrichDF2[,geneRatioColname]));
      }
      if (length(pathGenes) == 0) {
         pathGenes <- "pathGenes";
         if (verbose) {
            jamba::printDebug("enrichDF2enrichResult(): ",
               "deriving pathGenes '",
               pathGenes,
               "' from GeneRatio");
         }
         enrichDF2[,pathGenes] <- as.numeric(gsub("^.*[/]", "",
            enrichDF2[,geneRatioColname]));
      }
   } else {
      if (length(geneHits) == 0 || length(pathGenes) == 0) {
         stop("enrichDF2enrichResult() must have geneHits,pathGenes or geneRatioColname defined.");
      }
      if (verbose) {
         jamba::printDebug("enrichDF2enrichResult(): ",
            "deriving ",
            "'GeneRatio'",
            " from ",
            "geneHits/pathGenes");
      }
      geneRatioColname <- "GeneRatio";
      enrichDF2[,"GeneRatio"] <- pasteByRow(enrichDF2[,c(geneHits,pathGenes),drop=FALSE],
         sep="/");
   }

   enrichDF2 <- renameColumn(enrichDF2,
      from=c(geneRatioColname, pathGenes, geneHits),
      to=c("GeneRatio", "setSize", "Count"));
      #to=c("BgRatio", "setSize", "Count"));

   #enrichDF2[,"setSize"] <- enrichDF2[,pathGenes];
   #enrichDF2[,"Count"] <- enrichDF2[,geneHits];

   ## Re-order columns so "ID" is the first column
   if (verbose) {
      printDebug("enrichList2df(): ",
         "colnames(enrichDF2):", colnames(enrichDF2));
      printDebug("enrichList2df(): ",
         "class(enrichDF2):", class(enrichDF2));
   }
   enrichDF2a <- dplyr::select(enrichDF2,
      matches("^ID$"), everything());
   enrichDF2 <- enrichDF2a;
   if (verbose) {
      printDebug("enrichList2df(): ",
         "Done.");
   }

   gene <- jamba::mixedSort(unique(unlist(
      strsplit(enrichDF2[,"geneID"],
         "[/]+"))));
   if (verbose) {
      jamba::printDebug("enrichDF2enrichResult(): ",
         "identified ",
         jamba::formatInt(length(gene)),
         " total genes.");
   }

   #geneSets <- as(msigdbGmtT[enrichDF[,"itemsetID"],], "list");
   #names(geneSets) <- enrichDF[,"itemsetID"];

   ## Note geneSets is used in downstream methods, to represent the
   ## genes enriched which are present in a pathway, so it would
   ## be most correct not to use the GmtT items, which represents the
   ## full set of genes in a pathway.
   if (1 == 2 && !is.null(msigdbGmtT)) {
      geneSets <- as(msigdbGmtT[enrichDF2[,"ID"],], "list");
      names(geneSets) <- enrichDF2[,"ID"];
      universe <- mixedSort(msigdbGmtT@itemInfo[,1]);
   } else {
      if (verbose) {
         printDebug("enrichDF2enrichResult(): ",
            "Defined geneSets from delimited gene values.");
      }
      geneSets <- strsplit(enrichDF2[,"geneID"],
         "[/]");
      names(geneSets) <- enrichDF2[,"ID"];
      universe <- gene;
   }

   ## gene is list of hit genes tested for enrichment
   x <- new("enrichResult",
      result=enrichDF2,
      pvalueCutoff=pvalueCutoff,
      pAdjustMethod=pAdjustMethod,
      gene=as.character(gene),
      universe=universe,
      geneSets=geneSets,
      organism="UNKNOWN",
      keytype="UNKNOWN",
      ontology="UNKNOWN",
      readable=FALSE);
   x;
}

#' Prepare MultiEnrichMap data from enrichList
#'
#' Prepare MultiEnrichMap data from enrichList
#'
#' @family jam conversion functions
#'
#' @export
multiEnrichMap <- function
(enrichList,
geneHitList=NULL,
colorV=NULL,
nrow=NULL,
ncol=NULL,
byrow=FALSE,
enrichLabels=NULL,
subsetSets=NULL,
overlapThreshold=0.1,
cutoffRowMinP=0.05,
enrichBaseline=1.5,
enrichLens=0,
enrichNumLimit=4,
nEM=500,
topEnrichN=15,
topEnrichSources=c("Category", "Source"),
topEnrichCurateFrom=c("CP:.*"),
topEnrichCurateTo=c("CP"),
topEnrichSourceSubset=c("C2_CP", "C5_BP", "C5_CC", "C5_MF"),
topEnrichDescriptionGrep=NULL,
topEnrichNameGrep=NULL,
GmtTname="msigdbGmtTv50human",
keyColname="itemsetID",
geneColname="geneNames",
pvalueColname="P-value",
descriptionColname="Description",
descriptionCurateFrom=c("^Genes annotated by the GO term "),
descriptionCurateTo=c(""),
nameColname="Name",
pathGenes="pathGenes",
geneHits="geneHits",
geneDelim="[,/ ]+",
msigdbGmtT=msigdbGmtTv50human2,
verbose=FALSE,
...)
{
   ## Purpose is to create a multiEnrichMap object
   ##
   ## enrichList is a list of enrichment data.frames
   ## geneHitList is optionally the list of gene hits used
   ##    for each enrichment test, which would therefore contain
   ##    more genes than may be represented in enrichment results.
   ##    If not supplied, geneHitList is generated using only
   ##    genes represented in enrichment results, from enrichList.
   ## geneColname is used only when geneHitList is not supplied,
   ##    to generate geneHitList.
   ## geneDelim is a regular expression pattern used by strsplit() to
   ##    split the values in geneColname into a list of vectors of genes.
   ## colorV is a vector of colors, whose names match the names
   ##    of enrichList. If colorV is empty, rainbowCat() is called,
   ##    and colors are assigned to names(enrichList) in the same order.
   ##
   ## geneHits,pathGenes,msigdbGmtT,keyColname,geneColname, are used
   ##    by enrichDF2enrichResult()
   ##    to convert the combined enrichment data.frame to enrichResult
   ##    object.
   ##
   ## pathGenes is the colname in each enrichment data.frame which represents
   ##    the number of genes in each pathway, as tested for enrichment
   ## geneHits is the colname in each enrichment data.frame which represents
   ##    the number of gene hits which were present in the pathway.
   ##
   ## topEnrichN optional number, if supplied then topEnrichBySource() is
   ##    called with some sensible defaults intended for MsigDB data.
   ##
   ## subsetSets not currently implemented. Appears to be future work allowing
   ##    specific pathways to be specified directly.
   ##
   ## nrow,ncol,byrow parameters used for coloredrectangle igraph node color
   ##    placement. E.g. if there are 6 enrichList entries, one may define
   ##    nrow=2, ncol=3, byrow=TRUE. The enrichment colors will be applied in
   ##    2 rows with 3 columns, and filled in order, by row.
   ##
   ## cutoffRowMinP filters to ensure all pathways used have at least
   ##    one P-value at or below this threshold.
   ##
   ## msigdbGmtT is optionally the GmtT object used to test for enrichment,
   ##    and is used by enrichList2df() to convert a list of enrichment
   ##    data.frames into one combined data.frame.
   ##
   ## descriptionCurateFrom,descriptionCurateTo are vectors which are applied
   ##    in order to values in the colname defined by descriptionColname,
   ##    with the function gsub(). They are intended to remove lengthy prefix
   ##    labels, one in particular is used by Gene Ontology. However, the
   ##    result can be any valid gsub replacement, which allows potential
   ##    use of abbrevations to help shorten labels.
   ##
   ## enrichBaseline is the baseline enrichment P-value used when colorizing
   ##    the -log10(P-value) matrix of pathways in enrichIM. It is intended
   ##    to ensure values below this threshold are not colorized, for
   ##    example for entries deemed below a significance threshold.
   ##
   ## TODO:
   ## - create pathway-gene incidence matrix
   ##
   mem <- list();

   ## Some data checks
   if (suppressPackageStartupMessages(!require(igraph))) {
      stop("multiEnrichMap() requires the igraph package.")
   }
   if (suppressPackageStartupMessages(!require(DOSE))) {
      stop("multiEnrichMap() requires the DOSE package.")
   }
   if (suppressPackageStartupMessages(!require(IRanges))) {
      stop("multiEnrichMap() requires the IRanges package.")
   }
   if (length(names(enrichList)) == 0) {
      stop("multiEnrichMap() requires names(enrichList).")
   }

   ## Define some default colnames
   #nameColname <- "Name";
   #geneColname <- "geneNames";

   ## Add some basic information
   if (length(enrichLabels) == 0) {
      enrichLabels <- nameVector(names(enrichList));
   } else if (length(names(enrichLabels)) == 0) {
      names(enrichLabels) <- names(enrichList);
   }
   mem$enrichLabels <- enrichLabels;

   #####################################################################
   ## Define valid nrow and ncol for coloredrectangle igraph nodes
   if (length(nrow) == 0) {
      if (length(ncol) == 0) {
         nrow <- 1;
         ncol <- length(enrichList);
      } else {
         nrow <- ceiling(length(enrichList) / ncol);
      }
   } else {
      if (length(ncol) == 0) {
         ncol <- ceiling(length(enrichList) / nrow);
      } else if (ncol*nrow < length(enrichList)) {
         ncol <- ceiling(length(enrichList) / nrow);
      }
   }

   #####################################################################
   ## colors
   if (length(colorV) == 0) {
      colorV <- jamba::nameVector(colorjam::rainbowJam(length(enrichList)),
         names(enrichList));
   } else {
      colorV <- rep(colorV, length.out=length(enrichList));
      if (length(names(colorV)) == 0) {
         names(colorV) <- names(enrichList);
      }
   }
   useCols <- names(enrichList);
   colorV <- colorV[useCols];
   mem$colorV <- colorV;

   ## Get GmtT for use here
   #GmtT <- get(GmtTname, envir=.GlobalEnv);

   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "dim for each enrichList entry:");
      print(sdim(enrichList));
   }

   #####################################################################
   ## Create geneHitList if not supplied
   ## Note: this step occurs before filtering gene sets,
   ## which may be incorrect. However, this order will
   ## ensure the geneIM is comprehensive, beyond just the
   ## genes involved in enrichment of the filtered subset.
   if (length(geneHitList) == 0) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "creating geneHitList from enrichList.");
      }
      geneHitList <- lapply(enrichList, function(iDF){
         ## Split text field of delimited genes into proper vector
         if (!jamba::igrepHas("data.frame", class(iDF))) {
            jamba::mixedSort(unique(unlist(strsplit(
               as.data.frame(iDF)[[geneColname]], geneDelim))));
         } else {
            jamba::mixedSort(unique(unlist(strsplit(
               iDF[[geneColname]], geneDelim))));
         }
      });
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitList): ",
            lengths(geneHitList));
         #print(lengths(geneHitList));
      }
   }

   #####################################################################
   ## Optionally run topEnrichBySource()
   if (length(topEnrichN) > 0 && all(topEnrichN) > 0) {
      if (verbose) {
         printDebug("multiEnrichMap(): ",
            "running topEnrichBySource().");
      }
      enrichList <- lapply(enrichList, function(iDF){
         iDFtop <- topEnrichBySource(iDF,
            sourceColnames=topEnrichSources,
            n=topEnrichN,
            descriptionGrep=topEnrichDescriptionGrep,
            nameGrep=topEnrichNameGrep,
            curateFrom=topEnrichCurateFrom,
            curateTo=topEnrichCurateTo,
            sourceSubset=topEnrichSourceSubset);
         iDFtop;
      });
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "dims after topEnrichBySource():");
         print(sdim(enrichList));
      }
   }

   #####################################################################
   ## gene IM
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "geneIM <- list2im(geneHitList)");
   }
   geneIM <- list2im(geneHitList)[,,drop=FALSE];
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "geneIMcolors <- matrix2heatColors(geneIM)");
   }

   #####################################################################
   ## geneIM colors
   geneIMcolors <- colorjam::matrix2heatColors(x=geneIM,
      transformFunc=c,
      colorV=colorV,
      shareNumLimit=FALSE,
      numLimit=2);
   mem$geneHitList <- geneHitList;
   mem$geneIM <- geneIM;
   mem$geneIMcolors <- geneIMcolors;
   if (verbose && 1 == 2) {
      print(head(geneIM, 10));
      print(head(geneIMcolors, 10));
   }

   #####################################################################
   ## enrichIM incidence matrix using -log10(P-value)
   ##
   ## Note: the rownames use values from descriptionColname
   enrichLsetNames <- unique(unlist(lapply(enrichList, function(iDF){
      if (!jamba::igrepHas("data.frame", class(iDF))) {
         as.data.frame(iDF)[[nameColname]];
      } else {
         iDF[[nameColname]];
      }
   })));
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "enrichIM <- enrichList2IM() with    pvalueColname:",
         pvalueColname);
      printDebug("multiEnrichMap(): ",
         "head(enrichList[[1]]):");
      print(head(enrichList[[1]]));
   }
   enrichIM <- enrichList2IM(enrichList,
      valueColname=pvalueColname,
      keyColname=nameColname,
      #keyColname=keyColname,
      verbose=verbose,
      GmtT=msigdbGmtT)[enrichLsetNames,,drop=FALSE];

   ## Clean the rownames consistent with descriptionColname
   ## used in igraphs later on
   ##
   ## NOTE: we defer renaming the rownames here since it causes problems
   ## in trying to maintain consistent rownames and descriptionColname
   ## values.
   if (1 == 2) {
      rownames(enrichIM) <- memAdjustLabel(
         x=rownames(enrichIM),
         descriptionCurateFrom=descriptionCurateFrom,
         descriptionCurateTo=descriptionCurateTo);
   }

   enrichIMM <- as.matrix(enrichIM[,names(enrichList),drop=FALSE]);
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "head(enrichIM):");
      ch(head(enrichIM));
      #printDebug("multiEnrichMap(): ",
      #   "head(enrichIMM[,useCols]):");
      #ch(head(enrichIMM[,useCols]));
      #print(rowMins(na.rm=TRUE, head(enrichIMM[,useCols])) <= cutoffRowMinP);
   }

   #####################################################################
   ## enrichIM incidence matrix using geneCount
   #enrichLsetNames <- unique(unlist(lapply(enrichList, function(i){
   #   i[[nameColname]];
   #})));
   geneCountsGrep <- c(
      paste0("^", geneHits, "$"),
      geneHits,
      "^geneHits",
      "geneCount",
      "GeneRatio");
   if (!igrepHas("data.frame", class(enrichList[[1]]))) {
      geneCountColname <- head(provigrep(geneCountsGrep,
         colnames(as.data.frame(enrichList[[1]]))), 1);
   } else {
      geneCountColname <- head(provigrep(geneCountsGrep,
         colnames(enrichList[[1]])), 1);
   }
   if (length(geneCountColname) == 0) {
      stop(paste0("geneCountColname could not be found, by default it uses `geneHits`:",
         geneHits));
   }
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "enrichIM <- enrichList2IM() with geneCountColname:",
         geneCountColname);
   }

   enrichIMgeneCount <- enrichList2IM(enrichList,
      #valueColname="geneCount",
      keyColname=nameColname,
      valueColname=geneCountColname,
      verbose=verbose,
      GmtT=msigdbGmtT)[enrichLsetNames,,drop=FALSE];
   mem$enrichIMgeneCount <- as.matrix(enrichIMgeneCount[,names(enrichList),drop=FALSE]);

   #####################################################################
   ## enrichIM colors
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "enrichIMcolors <- matrix2heatColors(enrichIMM)");
      printDebug("multiEnrichMap(): ",
         "colorV:", colorV,
         fgText=list("orange", "dodgerblue", colorV));
      printDebug("multiEnrichMap(): ",
         "head(enrichIMM):");
      print(head(enrichIMM));
   }
   enrichIMcolors <- matrix2heatColors(x=-log10(enrichIMM),
      colorV=colorV,
      #lens=colLensFactorEnrich,
      lens=enrichLens,
      #numLimit=4,
      numLimit=enrichNumLimit,
      baseline=enrichBaseline);
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "head(enrichIMcolors)");
      print(head(enrichIMcolors));
   }

   #####################################################################
   ## Subset for at least one significant enrichment P-value
   i1use <- rownames(enrichIMM)[(rowMins(enrichIMM[,useCols], na.rm=TRUE) <= cutoffRowMinP)];
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "nrow(enrichIM):",
         formatInt(nrow(enrichIM)),
         ", nrow(filtered for minimum P-value):",
         formatInt(length(i1use)));
      #ch(head(enrichIMM));
   }

   #####################################################################
   ## Now make sure enrichList only contains these sets
   if (nrow(enrichIM) > length(i1use)) {
      if (verbose) {
         printDebug("multiEnrichMap(): ",
            "dims before filtering minimum P-value():");
         print(sdim(enrichList));
      }
      enrichList <- lapply(enrichList, function(iDF){
         if (!igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         subset(iDF, iDF[[nameColname]] %in% i1use);
      });
      if (verbose) {
         printDebug("multiEnrichMap(): ",
            "dims after filtering minimum P-value():");
         print(sdim(enrichList), Inf);
      }
      ## Now circle back and subset the enrichIM and enrichIMcolors rows
      enrichIM <- enrichIM[i1use,,drop=FALSE];
      enrichIMM <- enrichIMM[i1use,,drop=FALSE];
      enrichIMcolors <- enrichIMcolors[i1use,,drop=FALSE];
   }
   mem$enrichList <- enrichList;
   mem$enrichIM <- enrichIMM;
   mem$enrichIMcolors <- enrichIMcolors;

   #####################################################################
   ## Create one enrichment data.frame from the list
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "enrichDF <- enrichList2df(enrichList[c(",
         useCols,
         ")])");
      printDebug("multiEnrichMap(): ",
         "geneCountColname:",
         geneCountColname);
   }
   enrichDF <- enrichList2df(enrichList[useCols],
      msigdbGmtT=msigdbGmtT,
      keyColname=keyColname,
      geneColname=geneColname,
      geneCountColname=geneCountColname,
      pvalueColname=pvalueColname,
      verbose=verbose);

   #####################################################################
   ## Some cleaning of Description
   if (length(descriptionColname) == 1 &&
         descriptionColname %in% colnames(enrichDF)) {
      if (verbose) {
         printDebug("multiEnrichMap(): ",
            "cleaning Description column:",
            descriptionColname);
         #print(lengths(enrichList[useCols]));
      }
      descriptionColnameFull <- paste0(descriptionColname, "Full");
      enrichDF[,descriptionColnameFull] <- enrichDF[,descriptionColname];
      enrichDF[,descriptionColname] <- memAdjustLabel(
         x=enrichDF[,descriptionColname],
         descriptionCurateFrom=descriptionCurateFrom,
         descriptionCurateTo=descriptionCurateTo);
   }

   #####################################################################
   ## Convert the combined enrichDF to enrichResult
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "head(enrichDF):");
      print(head(as.data.frame(enrichDF)));
      printDebug("multiEnrichMap(): ",
         "enrichER <- enrichDF2enrichResult(), keyColname:",
         keyColname);
   }
   enrichER <- enrichDF2enrichResult(enrichDF,
      geneHits=geneHits,
      pathGenes=pathGenes,
      keyColname=keyColname,
      geneColname=geneColname,
      geneCountColname=geneCountColname,
      geneDelim=geneDelim,
      pvalueColname=pvalueColname,
      msigdbGmtT=msigdbGmtT,
      verbose=verbose);

   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "head(enrichER):");
      print(head(as.data.frame(enrichER)));
   }

   mem$multiEnrichDF <- enrichDF;
   mem$multiEnrichResult <- enrichER;


   #####################################################################
   ## Incidence matrix of genes and pathways
   memIM <- list2im(
      strsplit(
         nameVector(mem$multiEnrichDF[,c(geneColname,nameColname)]),
         geneDelim));
   mem$memIM <- memIM;


   #####################################################################
   ## Convert enrichResult to enrichMap igraph network
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "converting enrichER to igraph enrichMap with enrichMapJam().");
   }
   enrichEM <- multienrichjam::enrichMapJam(enrichER,
      overlapThreshold=overlapThreshold,
      msigdbGmtT=msigdbGmtT,
      doPlot=FALSE,
      n=nEM,
      keyColname="ID",
      nodeLabel=c(nameColname, descriptionColname, keyColname, "ID"),
      vertex.label.cex=0.5,
      verbose=verbose);
   mem$multiEnrichMap <- enrichEM;

   ## Convert EnrichMap to piegraph
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "running igraph2pieGraph() on enrichMap.");
      printDebug("multiEnrichMap(): ",
         "head(enrichIMcolors)");
      print(head(enrichIMcolors));
   }
   enrichEMpieUse <- igraph2pieGraph(g=enrichEM,
      defineLayout=FALSE,
      valueIMcolors=enrichIMcolors[i1use,useCols,drop=FALSE],
      verbose=verbose);

   ## Use colored rectangles
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "running rectifyPiegraph() on enrichMap.");
   }
   enrichEMpieUseSub2 <- rectifyPiegraph(enrichEMpieUse,
      nrow=nrow,
      ncol=ncol,
      byrow=byrow);
   mem$multiEnrichMap2 <- enrichEMpieUseSub2;

   #######################################################
   ## Create a CnetPlot
   ## Consider omitting this step if it is slow with large
   ## data, and if downstream workflows would typically
   ## only need a Cnet Plot on a subset of pathways and
   ## genes.
   gCt <- nrow(enrichER);
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "creating cnetPlot with cnetplotJam().");
   }
   gCnet <- cnetplotJam(enrichER,
      showCategory=gCt,
      categorySize=geneCountColname,
      doPlot=FALSE,
      nodeLabel=c(nameColname, descriptionColname, keyColname, "ID"),
      verbose=verbose);
   V(gCnet)$nodeType <- "Gene";
   V(gCnet)[seq_len(gCt)]$nodeType <- "Set";
   mem$multiCnetPlot <- gCnet;

   #######################################################
   ## Convert to coloredrectangle
   V(gCnet)[seq_len(gCt)]$name <- toupper(V(gCnet)[seq_len(gCt)]$name);
   ## Enrichment IM colors
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "running igraph2pieGraph(",
         "enrichIMcolors",
         ") on Cnet Plot.");
   }
   gCnetPie1 <- igraph2pieGraph(g=gCnet,
      defineLayout=FALSE,
      valueIMcolors=enrichIMcolors[i1use,useCols,drop=FALSE],
      verbose=verbose);
   mem$multiCnetPlot1 <- gCnetPie1;
   ## Gene IM colors
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "running igraph2pieGraph(",
         "geneIMcolors",
         ").");
   }
   gCnetPie <- igraph2pieGraph(g=gCnetPie1,
      defineLayout=FALSE,
      valueIMcolors=geneIMcolors[,useCols,drop=FALSE],
      verbose=verbose);
   mem$multiCnetPlot1b <- gCnetPie;

   #######################################################
   ## Now convert CnetPlot to use coloredrectangle
   if (verbose) {
      printDebug("multiEnrichMap(): ",
         "running rectifyPiegraph() on Cnet Plot.");
   }
   gCnetPie2 <- rectifyPiegraph(gCnetPie,
      nrow=nrow,
      ncol=ncol,
      byrow=byrow);
   mem$multiCnetPlot2 <- gCnetPie2;

   #######################################################
   ## Add all colnames to the mem object
   colnamesL <- list(
      geneColname=geneColname,
      keyColname=keyColname,
      nameColname=nameColname,
      descriptionColname=descriptionColname,
      pvalueColname=pvalueColname);
   mem$colnames <- colnamesL;

   return(mem);
}

#' Convert enrichList to IM incidence matrix
#'
#' Convert enrichList to IM incidence matrix
#'
#' @family jam conversion functions
#'
#' @export
enrichList2IM <- function
(enrichList,
 GmtT=NULL,
 addAnnotations=TRUE,
 keyColname="Name",
 valueColname="P-value",
 emptyValue=1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a list of enrichment data.frames as produced by
   ## enrichSimpleM(), and the corresponding GmtT object, and produce
   ## an incidence matrix whose values are the enrichment P-value
   ##
   ## addAnnotations=TRUE will add annotation columns from GmtT@itemsetInto
   ## to the resulting data.frame
   ##
   ## Examples:
   ## enrichSubIMP <- enrichList2IM(enrichSubL, msigdbGmtTv50mouseV2);
   ## enrichSubIMP <- enrichList2IM(enrichSubL, msigdbGmtTv50mouseV2, keyColname="Description", valueColname="geneHits", emptyValue=0);
   ##
   enrichIMP <- as.data.frame(list2imSigned(lapply(enrichList, function(iDF){
      if (!jamba::igrepHas("data.frame", class(iDF))) {
         iDF <- as.data.frame(iDF);
      }
      ## If "GeneRatio" then parse out the geneCount value
      if (verbose) {
         printDebug("enrichList2IM(): ",
            "keyColname:", keyColname);
         printDebug("enrichList2IM(): ",
            "valueColname:", valueColname);
      }
      if (igrepHas("GeneRatio", valueColname)) {
         iDF[,valueColname] <- gsub("[/].*$", "", iDF[,valueColname]);
         if (length(grep("^[0-9]*$", iDF[,valueColname])) == nrow(iDF)) {
            iDF[,valueColname] <- as.numeric(iDF[,valueColname]);
         }
      }
      if (verbose) {
         printDebug("enrichList2IM(): ",
            "head(iDF)");
         print(head(iDF));
      }
      jamba::nameVector(iDF[,c(valueColname,keyColname)]);
   })));
   ## Empty values should be 1 instead of 0
   if (emptyValue != 0) {
      enrichIMP[is.na(enrichIMP) | enrichIMP == 0] <- emptyValue;
   }

   if (addAnnotations && length(GmtT) > 0) {
      ## Add information about pathways to the data.frame
      enrichIMPinfo <- GmtT@itemsetInfo[match(rownames(enrichIMP), GmtT@itemsetInfo[,keyColname]),];
      enrichIMP[,colnames(enrichIMPinfo)] <- enrichIMPinfo;
   }
   return(enrichIMP);
}

#' Convert enrichList to data.frame
#'
#' Convert enrichList to data.frame
#'
#' @export
enrichList2df <- function
(enrichList,
 keyColname="itemsetID",
 geneColname="geneNames",
 geneCountColname="geneCount",
 pvalueColname="P-value",
 pvalueFloor=1e-200,
 msigdbGmtT=NULL,
 verbose=FALSE,
 debug=0,
...)
{
   ## Purpose is to combine a list of enrichment data.frames into one data.frame
   if (!suppressPackageStartupMessages(require(matrixStats))) {
      stop("enrichList2df() requires the matrixStats package.");
   }
   if (verbose) {
      printDebug("enrichList2df(): ",
         "keyColname:",
         keyColname);
      printDebug("enrichList2df(): ",
         "geneColname:",
         geneColname);
      printDebug("enrichList2df(): ",
         "geneCountColname:",
         geneCountColname);
      printDebug("enrichList2df(): ",
         "pvalueColname:",
         pvalueColname);
   }
   ## Keep the best result for these columns as the exemplar
   enrichCols <- c(`P-value`="lo",
      pathGenes="hi",
      geneHits="hi",
      geneCount="hi");
   names(enrichCols)[1] <- head(pvalueColname, 1);
   names(enrichCols)[4] <- head(geneCountColname, 1);
   enrichCols <- enrichCols[unique(names(enrichCols))];

   ## Get first non-NULL data.frame from enrichList
   if (!igrepHas("data.frame", class(head(rmNULL(enrichList), 1)))) {
      iDF <- as.data.frame(rmNULL(enrichList)[[1]]);
   } else {
      iDF <- rmNULL(enrichList)[[1]];
   }

   if (verbose) {
      printDebug("enrichList2df(): ",
         "enrichCols (before):");
      print(enrichCols);
      printDebug("enrichList2df(): ",
         "colnames(iDF):", colnames(iDF));
   }
   enrichCols <- enrichCols[names(enrichCols) %in% colnames(iDF)];
   enrichColsHi <- names(enrichCols)[enrichCols %in% "hi"];
   #c("pathGenes","geneHits");
   enrichColsLo <- names(enrichCols)[enrichCols %in% "lo"];
   #enrichColsLo <- c("P-value");
   keepCols <- setdiff(unvigrep("gene", colnames(iDF)),
      c(enrichColsHi, enrichColsLo, keyColname, geneColname));

   ## Create a P-value incidence matrix
   if (verbose) {
      printDebug("enrichList2df(): ",
         "enrichCols (after):");
      print(enrichCols);
      printDebug("enrichList2df(): ",
         "sdim(enrichL):");
      print(sdim(enrichList));
   }
   enrichValuesM <- do.call(cbind, lapply(nameVector(names(enrichCols)), function(iCol){
      useType <- enrichCols[iCol];
      enrichIMP <- list2imSigned(lapply(enrichList, function(iDF){
         if (!igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         if (useType %in% "lo" && any(iDF[,iCol] <= pvalueFloor)) {
            printDebug("Some ", iCol, " values are less than ",
               "pvalueFloor:",
               pvalueFloor);
            print(table(iDF[,iCol] <= pvalueFloor));
            iDF[iDF[,iCol] <= pvalueFloor, iCol] <- pvalueFloor;
         } else if (igrepHas("[/]", iDF[,iCol])) {
            iDF[,iCol] <- as.numeric(gsub("[/].*$", "", iDF[,iCol]));
         }
         if (length(tcount(iDF[[keyColname]], minCount=2)) > 0) {
            stop("enrichList2df(): There are duplicate values in iDF[[keyColname]], please resolve.");
         }
         if (verbose) {
            jamba::printDebug("enrichList2df(): ",
               "head(iDF):");
            print(head(iDF));
         }
         nameVector(iDF[,c(iCol,keyColname)]);
      }));
      if (useType %in% "lo") {
         enrichIMP[enrichIMP == 0] <- 1;
         nameVector(rowMins(enrichIMP), rownames(enrichIMP));
      } else if (useType %in% "hi") {
         nameVector(rowMaxs(enrichIMP), rownames(enrichIMP));
      }
   }));
   if (verbose) {
      printDebug("enrichList2df(): ",
         "dim(enrichValuesM):",
         dim(enrichValuesM));
      print(head(enrichValuesM));
   }
   if (debug == 1) {
      return(list(enrichCols=enrichCols,
         enrichValuesM=enrichValuesM,
         enrichList=enrichList));
   }

   ## Generate list with genes per pathway
   allGenes <- mixedSort(unique(unlist(lapply(enrichList, function(iDF){
      if (!igrepHas("data.frame", class(iDF))) {
         iDF <- as.data.frame(iDF);
      }
      unlist(strsplit(iDF[,geneColname], ","));
   }))));

   ## If GmtT is supplied, use it to determine genes per pathway
   if (length(msigdbGmtT) > 0) {
      enrichGeneL <- as(msigdbGmtT[match(rownames(enrichValuesM), msigdbGmtT@itemsetInfo[,keyColname]),
         rmNA(match(allGenes, msigdbGmtT@itemInfo[,1]))], "list");
      names(enrichGeneL) <- rownames(enrichValuesM);
      enrichGeneVL <- list(cPaste2(enrichGeneL, doSort=FALSE));
      names(enrichGeneVL) <- geneColname;
      enrichGeneLen <- lengths(enrichGeneL);
   } else {
      ## if GmtT is not supplied, use the pathway enrichment data as a substitute
      enrichL1L1 <- lapply(nameVectorN(enrichList), function(iName){
         iDF <- enrichList[[iName]];
         if (!igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         iDF <- renameColumn(iDF,
            from=geneColname,
            to=iName);
         iDF[,c(keyColname,iName)];
      });
      enrichL1L <- mergeAllXY(enrichL1L1);
      enrichL1V <- nameVector(gsub("^[,]+|[,]+$", "",
         pasteByRow(enrichL1L[,-match(keyColname, colnames(enrichL1L)),drop=FALSE],
            sep=",")),
         enrichL1L[,keyColname]);
      enrichGeneL <- as.list(unique(CharacterList(strsplit(enrichL1V, "[,]+"))));
      enrichGeneVL <- list(cPaste2(enrichGeneL, doSort=FALSE));
      names(enrichGeneVL) <- geneColname;
      enrichGeneLen <- lengths(enrichGeneL);
   }

   if (debug == 2) {
      return(list(enrichCols=enrichCols,
         enrichValuesM=enrichValuesM,
         enrichList=enrichList,
         enrichGeneLen=enrichGeneLen,
         enrichGeneL=enrichGeneL));
   }

   ## Create data.frame with annotation columns, only keep the first
   ## occurrence of any non-NA value
   if (verbose) {
      printDebug("enrichList2df(): ",
         "head(enrichL1L):");
      print(str(head(enrichL1L)));
      printDebug("enrichList2df(): ",
         "dim(enrichValuesM):",
         dim(enrichValuesM));
   }
   keepColDF <- renameColumn(
      data.frame(row.names=rownames(enrichValuesM),
         keyColname=rep(NA, nrow(enrichValuesM))),
      from="keyColname",
      to=keyColname);
   for (iName in names(enrichList)) {
      iDF <- enrichList[[iName]];
      keyVals <- iDF[,keyColname];
      keyValsUse <- setdiff(keyVals, rmNA(keepColDF[,keyColname]));
      if (verbose) {
         printDebug("iName:", iName,
            ", length(keyVals):", length(keyVals),
            ", length(keyValsUse):", length(keyValsUse));
      }
      if (length(keyValsUse) > 0) {
         keepColDF[keyValsUse,keyColname] <- keyValsUse;
         for (keepCol in keepCols) {
            keepColDF[keyValsUse,keepCol] <- iDF[match(keyValsUse, iDF[,keyColname]),keepCol];
         }
      }
      rm(iDF);
   }
   if (debug == 3) {
      return(list(enrichCols=enrichCols,
         enrichValuesM=enrichValuesM,
         enrichList=enrichList,
         enrichGeneLen=enrichGeneLen,
         enrichGeneL=enrichGeneL,
         keepColDF=keepColDF));
   }

   enrichDF <- data.frame(check.names=FALSE, stringsAsFactors=FALSE,
      enrichValuesM, keepColDF, as.data.frame(enrichGeneVL));
   whichCol1 <- max(which(colnames(enrichDF) %in% names(enrichCols))) + 1;
   enrichDF <- insertDFcols(enrichDF, colnum=whichCol1,
      insertDF=data.frame(allGeneHits=enrichGeneLen));
   enrichDF;
}

#' Merge a list of data.frames keeping all rows
#'
#' Merge a list of data.frames keeping all rows
#'
#' This function is an extension to the `base::merge.data.frame()`
#' function, with the default arguments `all.x=TRUE, all.y=TRUE` in
#' order to maintain all rows of all input data.frames.
#'
#' @family jam list functions
#'
#' @export
mergeAllXY <- function
(...)
{
   ## Purpose is a simple wrapper to merge(..., all.x=TRUE, all.y=TRUE);
   ##
   ## But detect whether 'x' and 'y' are defined, and if not, then define them
   inList <- list(...);

   ## name inList for any un-named entries
   if (is.null(names(inList))) {
      names(inList) <- makeNames(rep("x", length(inList)));
   } else {
      names(inList) <- makeNames(names(inList));
   }

   ## Filter out parameters for the merge.data.frame() function
   mergeArgs <- unvigrep("^(x|y|[.]{3})$", names(formals(base:::merge.data.frame)));
   inListArgs <- inList[names(inList) %in% mergeArgs];
   inList <- inList[!names(inList) %in% mergeArgs];

   ## Now un-nest the list of entries remaining, so they're all non-list entities
   inList <- unnestList(inList);

   ## Filter for only data.frame or matrix like objects
   inListClass <- sapply(inList, function(i){
      class(i);
   });
   inList <- inList[igrep("data.*frame|matrix", inListClass)];
   inListClass <- sapply(inList, function(i){
      class(i);
   });
   ## Convert anything not a data.frame into such an object (ha!)
   if (any(!inListClass %in% c("data.frame"))) {
      whichNonDF <- which(!inListClass %in% c("data.frame"));
      inList[whichNonDF] <- lapply(inList[whichNonDF], function(i){
         newDF <- as.data.frame(i);
         rownames(newDF) <- rownames(i);
         colnames(newDF) <- colnames(i);
         newDF;
      })
   }

   ## Accept list of data.frames and run iterative merge on them
   ## Usually seen as one '...' element which is a list of data.frames
   if (length(inList) == 1) {
      return(inList[[1]]);
   } else if (length(inList) == 200) {
      ## If two list elements, we just call merge() once using x and y like normal
      x <- inList[[1]];
      y <- inList[[2]];

      if (length(inListArgs) > 0) {
         x <- do.call(merge,
            c(alist(x=x,
               y=y,
               all.x=TRUE,
               all.y=TRUE),
            inListArgs));
      } else {
         x <- merge(x,
            y,
            all.x=TRUE,
            all.y=TRUE);
      }
      return(x);
   } else if (length(inList) >= 2) {
      ## Run iterative merge() methods
      x <- inList[[1]];
      for (i in 2:length(inList)) {
         y <- inList[[i]];
         ## If we have merge.data.frame() arguments to pass,
         ## use do.call() so we can separately pass those function arguments
         if (length(inListArgs) > 0) {
            x <- do.call(merge, c(alist(x=x, y=y, all.x=TRUE, all.y=TRUE), inListArgs));
         } else {
            x <- merge(x, y, all.x=TRUE, all.y=TRUE);
         }
      }
      return(x);
   } else {
      stop("Could not match input to expected types, e.g. list of data.frames or matrices.");
   }
}

#' Un-nest a nested list into a simple list
#'
#' Un-nest a nested list into a simple list
#'
#' This function inspects a list, and unlists each entry
#' resulting in a simple list of non-list entries as a result.
#'
#' @family jam list functions
#'
#' @param x list potentially containing lists
#' @param unnamedBase character value used as a base for naming any
#'    un-named lists, using the format `makeNamesFunc(rep(unnamedBase, n))`.
#' @param sep character delimiter used between nested list names.
#' @param makeNamesFunc function that takes a character vector and returns
#'    non-duplicated character vector of equal length. By default it
#'    uses `jamba::makeNames()`.
#' @param stopClasses vector of classes that should not be un-nested,
#'    useful in case some classes inherit list properties.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' L <- list(A=letters[1:10], B=list(C=LETTERS[3:9], D=letters[4:11]));
#' L;
#' unnestList(L);
#'
#' # inspect the data using str()
#' str(L);
#' str(unnestList(L))
#'
#' @export
unnestList <- function
(x,
 unnamedBase="x",
 parentName=NULL,
 sep=".",
 makeNameFunc=jamba::makeNames,
 stopClasses=c("dendrogram", "data.frame", "matrix"),
 ...)
{
   ## Purpose is to take a list of lists, and un-nest them
   ## into a list with depth=1, but all the non-list elements
   ## contained within it
   newList <- list();

   ## Create default names if they don't exist already
   if (is.null(names(x))) {
      names(x) <- rep(unnamedBase, length(x));
   }
   emptyNamesX <- which(names(x) %in% c("", NA));
   if (any(emptyNamesX)) {
      names(x)[emptyNamesX] <- rep(unnamedBase, length(emptyNamesX));
   }

   ## add prefix if it were specified
   if (!is.null(parentName)) {
      names(x) <- paste(parentName, names(x), sep=sep);
   }

   ## Make sure names are unique
   names(x) <- makeNameFunc(names(x), ...);

   ## Iterate each list until we hit a non-list entry
   if (inherits(x, "list") || is.list(x)) {
      for (j in names(x)) {
         i <- x[[j]];
         ## If we reached a list, and if the class isn't something
         ## we want to allow (e.g. dendrogram) then unnest one layer deeper
         if (inherits(i, "list") || is.list(i) && !class(i) %in% stopClasses) {
            i <- unnestList(x=i, parentName=j, ...);
         } else {
            i <- list(i);
            names(i) <- j;
         }
         newList <- c(newList, i);
      };
   } else {
      newList <- x;
   }
   newList;
}

#' Create enrichMap igraph object
#'
#' Create enrichMap igraph object
#'
#' @export
enrichMapJam <- function
(x,
n=50,
fixed=TRUE,
vertex.label.font=1,
vertex.label.cex=1,
nodeLabel=c("Name","Description","ID"),
descriptionColname="Description",
keyColname="ID",
nodeLabelFunc=function(i){paste(collapse="\n",strwrap(width=30, ucfirst(tolower(gsub("_", " ", i)))))},
overlapThreshold=0.2,
msigdbGmtT=NULL,
method=2,
verbose=FALSE,
...)
{
   ## Purpose is to customize enrichMap() to work with data
   ## generated outside clusterProfiler
   ##
   if (suppressPackageStartupMessages(!require(reshape2))) {
      stop("enrichMapJam() requires the reshape2 package is required for melt().");
   }
   if (suppressPackageStartupMessages(!require(igraph))) {
      stop("enrichMapJam() requires the igraph package.");
   }
   if (suppressPackageStartupMessages(!require(DOSE))) {
      stop("enrichMapJam() requires the DOSE package.");
   }
   if (is.null(nodeLabelFunc)) {
      nodeLabelFunc <- function(i){
         paste(collapse="\n",strwrap(width=30, ucfirst(gsub("_", " ", tolower(i)))));
      }
   }
   if (igrepHas("data.*frame", class(x))) {
      if (verbose) {
         jamba::printDebug("enrichMapJam(): ",
            "calling enrichDF2enrichResult()");
      }
      x <- enrichDF2enrichResult(x,
         msigdbGmtT=msigdbGmtT,
         verbose=verbose);
   }
   if (is(x, "gseaResult")) {
      geneSets <- x@geneSets;
   } else if (is(x, "enrichResult")) {
      geneSets <- x@geneSets;
      #geneSets <- geneInCategory(x);
   }
   y <- as.data.frame(x);

   ## Make sure nodeLabel is a colname
   if (verbose) {
      jamba::printDebug("enrichMapJam(): ",
         "nodeLabel (before):",
         nodeLabel);
   }
   nodeLabel <- head(intersect(nodeLabel, colnames(y)), 1);
   if (verbose) {
      jamba::printDebug("enrichMapJam(): ",
         "nodeLabel (found in y):",
         nodeLabel);
      jamba::printDebug("enrichMapJam(): ",
         "colnames(y):",
         colnames(y));
   }

   if (nrow(y) < n) {
      n <- nrow(y);
   } else {
      y <- y[1:n,,drop=FALSE];
   }
   if (verbose) {
      jamba::printDebug("enrichMapJam(): ",
         "n:",
         n);
      jamba::printDebug("enrichMapJam(): ",
         "head(y):");
      print(head(y));
   }

   if (n == 0) {
      stop("`enrichMapJam()` found no enriched terms.")
   } else if (n == 1) {
      g <- igraph::graph.empty(0, directed=FALSE);
      g <- igraph::add_vertices(g, nv=1);
      V(g)$name <- y[, descriptionColname];
      V(g)$color <- "red";
   } else {
      pvalue <- jamba::nameVector(y$pvalue, y[[nodeLabel]]);

      ## Define the vector of identifiers
      id <- y[,keyColname];
      if (verbose) {
         jamba::printDebug("enrichMapJam(): ",
            "id:",
            id);
      }
      #id <- y[,keyColname];
      geneSets <- geneSets[id];
      n <- nrow(y)

      if (method == 1) {
         ## Manual all-by-all overlap_ratio() method
         w <- matrix(NA, nrow = n, ncol = n)
         colnames(w) <- rownames(w) <- y[[nodeLabel]];
         for (i in 1:n) {
            for (j in i:n) {
               #return(list(x=x, y=y, geneSets=geneSets, id=id, i=i, j=j));
               w[i, j] = enrichplot:::overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
         }
      } else {
         ## Jaccard coefficient is given as output from
         ## 1-dist(method="binary")
         wIM <- list2im(geneSets);
         w <- 1-as.matrix(dist(t(wIM), method="binary"));
         colnames(w) <- rownames(w) <- y[[nodeLabel]][match(colnames(w), y$ID)];
      }
      wd <- reshape2::melt(w);

      if (method == 1) {
         wd <- wd[wd[, 1] != wd[, 2], ];
         wd <- wd[!is.na(wd[, 3]), ];
      } else {
         wd1 <- match(wd[,1], colnames(w));
         wd2 <- match(wd[,2], colnames(w));
         wd <- wd[wd1 > wd2,,drop=FALSE];
      }

      g <- igraph::graph.data.frame(wd[, -3], directed=FALSE);
      igraph::E(g)$width <- sqrt(wd[, 3] * 20);
      igraph::V(g)$pvalue <- pvalue[V(g)$name];

      ## Attempt to merge annotations from the enrichResult object
      iMatch <- match(V(g)$name, y[[nodeLabel]]);
      if (verbose) {
         jamba::printDebug("enrichMapJam(): ",
            "merging annotations from enrichResult objects");
      }
      if (!any(is.na(iMatch))) {
         #printDebug("match() worked with enrichResult data.frame Name colname.");
         iColnames <- unvigrep("^name$", colnames(y));
         for (iY in iColnames) {
            g <- g %>% set_vertex_attr(iY, value=y[iMatch,,drop=FALSE][[iY]]);
         }
      } else {
         ## Attempt to merge additional pathway annotation from GmtT
         #printDebug("match() worked with enrichResult data.frame Name colname.");
         if (length(msigdbGmtT) > 0) {
            iMatch <- match(V(g)$name, msigdbGmtT@itemsetInfo$Name);
            iMatchWhich <- which(!is.na(iMatch));
            if (length(iMatchWhich) > 0) {
               for (iCol1 in setdiff(colnames(msigdbGmtT@itemsetInfo), "Name")) {
                  g <- set_vertex_attr(g, iCol1, V(g)[iMatchWhich], msigdbGmtT@itemsetInfo[iMatch[iMatchWhich],iCol1]);
               }
            }
         }
      }

      ## Delete edges where overlap is below a threshold
      E(g)$overlap <- wd[,3];
      g <- delete.edges(g, E(g)[E(g)$overlap < overlapThreshold]);

      pvalue <- V(g)$pvalue;

      if (method == 1) {
         ## where is color_scale?
         cols <- color_scale("red", "#E5C494");
         V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))];
      } else {
         #nodeColor <- customNumColors(-log10(pvalue), col="Reds",
         #   colLensFactor=5);
         nodeColor <- vals2colorLevels(-log10(pvalue),
            col="Reds",
            numLimit=4,
            baseline=0,
            lens=2);
         V(g)$color <- nodeColor;
         #nodeColor <- customNumColors(rank(-pvalue), col="Reds", colLensFactor=5);
      }
      if (is(x, "gseaResult")) {
         cnt <- y$setSize/10
      } else if (is(x, "enrichResult")) {
         cnt <- jamba::nameVector(y$Count, y[[nodeLabel]]);
      }
      cnt2 <- cnt[V(g)$name]
      V(g)$size <- log10(cnt2) * 10;

      #if (!nodeLabelFunc %in% c(FALSE)) {
      V(g)$label <- sapply(V(g)$name, nodeLabelFunc);
      #}
   }
   invisible(g);
}

#' Subset Cnet igraph
#'
#' Subset Cnet igraph
#'
#' This function produces a subset of a Cnet igraph based upon supplied
#' set names or gene names. This function is intended to be a convenient
#' method of filtering a Cnet igraph to a pre-defined set of "Set"
#' names.
#'
#' The function assumes graph nodes have an attribute `"nodeType"` with
#' values either `"Set"` or `"Gene"` to indicate the type of node.
#'
#' When `includeSets` is supplied, the graph is subsetted to include
#' only nodes with `nodeType="Set"` with matching `V(gCnet)$name` or
#' `V(gCnet)$label`. Then only neighboring nodes are retained, thus
#' removing any nodes with `nodeType="Gene"` that do not connect to
#' any of the given Set nodes.
#' The result is a proper Cnet igraph that only contains
#' Gene nodes connected to the subset of Set nodes.
#'
#' If `includeGenes` is supplied, the graph is subsetted to include
#' only nodes with `nodeType="Gene"` with matching `V(gCnet)$name` or
#' `V(gCnet)$label`.
#'
#' When `removeSinglets=TRUE` then any nodes that have no remaining
#' edges are removed. Especially when supplying `includeGenes`, this
#' option is useful to hide any Set nodes that have no connected Gene
#' nodes.
#'
#' @param gCnet igraph object representing Cnet concept network data
#' @param includeSets character vector, or NULL, containing the set
#'    names or labels to retain.
#' @param includeGenes character vector, or NULL, containing the gene
#'    names or labels to retain.
#' @param removeSinglets logical whether to remove singlet graph nodes,
#'    which are nodes that have no remaining edges.
#' @param minSetDegree integer value indicating the minimum number
#'    of edges each Set node must have to be retained in the resulting
#'    igraph. Use `minSetDegree=2` to retain only Set nodes that
#'    have multiple connected Gene nodes.
#' @param minGeneDegree integer value indicating the minimum number
#'    of edges each Gene node must have to be retained in the resulting
#'    igraph. Use `minGeneDegree=2` to retain only Gene nodes that
#'    connect to multiple Set nodes.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @export
subsetCnetIgraph <- function
(gCnet,
 includeSets=NULL,
 includeGenes=NULL,
 removeSinglets=TRUE,
 minSetDegree=1,
 minGeneDegree=1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take an Cnet igraph object and subset
   ## by set name or gene symbol
   ##########################################
   ## Optionally subset for certain pathways
   if (length(includeSets) > 0) {
      if (length(V(gCnet)$label) == 0) {
         includeV <- which(V(gCnet)$nodeType %in% "Set" &
               (
                  toupper(V(gCnet)$name) %in% toupper(includeSets)
               ));
      } else {
         includeV <- which(V(gCnet)$nodeType %in% "Set" &
               (
                  toupper(V(gCnet)$name) %in% toupper(includeSets) |
                  toupper(V(gCnet)$label) %in% toupper(includeSets)
               ));
      }
      includeV2 <- unique(unlist(lapply(includeV, function(v){
         as.numeric(
            neighbors(gCnet,
               v=v,
               mode="all"));
      })));
      includeVall <- sort(unique(c(includeV, includeV2)));
      if (verbose) {
         printDebug("subsetCnetIgraph(): ",
            "Filtered ",
            formatInt(sum(V(gCnet)$nodeType %in% "Set")),
            " Set nodes using ",
            formatInt(length(includeSets)),
            " includeSets down to ",
            formatInt(length(includeV)),
            " sets and ",
            formatInt(length(includeV2)),
            " genes in the Cnet igraph object.");
         whichNodeSets <- which(V(gCnet)$nodeType %in% "Set");
      }
      gCnet <- igraph::subgraph(gCnet,
         includeVall);
   }

   ##########################################
   ## Optionally subset for certain genes
   if (length(includeGenes) > 0) {
      keepSetNodes <- which(V(gCnet)$nodeType %in% "Set");
      if (length(V(gCnet)$label) == 0) {
         keepGeneNodes <- which(
            V(gCnet)$nodeType %in% "Gene" &
               toupper(V(gCnet)$name) %in% toupper(includeGenes)
         );
      } else {
         keepGeneNodes <- which(
            V(gCnet)$nodeType %in% "Gene" &
               (toupper(V(gCnet)$name) %in% toupper(includeGenes) |
                     toupper(V(gCnet)$label) %in% toupper(includeGenes))
         );
      }
      keepNodes <- sort(unique(c(keepSetNodes, keepGeneNodes)));
      if (verbose) {
         printDebug("subsetCnetIgraph(): ",
            "Filtered ",
            formatInt(length(includeSets)),
            " includeGenes down to ",
            formatInt(length(keepGeneNodes)),
            " genes and ",
            formatInt(length(keepSetNodes)),
            " sets in the Cnet igraph object.");
      }
      gCnet <- igraph::subgraph(gCnet,
         keepNodes);
   }

   #####################################################
   ## Polish the igraph by removing nodes with no edges
   iDegree <- degree(gCnet);
   if (removeSinglets) {
      if (any(iDegree) == 0) {
         if (verbose) {
            printDebug("subsetCnetIgraph(): ",
               "Filtered ",
               formatInt(length(iDegree)),
               " nodes to remove ",
               formatInt(sum(iDegree == 0)),
               " nodes with no connections.");
         }
         gCnet <- igraph::subgraph(gCnet,
            which(iDegree > 0));
      }
   }
   #####################################################
   ## Optionally subset by degree of Set and Gene nodes
   if (length(minSetDegree) > 0) {
      dropSetNodes <- (V(gCnet)$nodeType %in% "Set" &
            iDegree < minSetDegree);
      if (any(dropSetNodes)) {
         if (verbose) {
            jamba::printDebug("subsetCnetIgraph(): ",
               "Dropping ",
               formatInt(sum(dropSetNodes)),
               " set nodes with less than degree:",
               minSetDegree);
         }
         gCnet <- igraph::subgraph(gCnet,
            which(!dropSetNodes));
      }
   }
   if (length(minGeneDegree) > 0) {
      dropGeneNodes <- (V(gCnet)$nodeType %in% "Gene" &
         iDegree < minGeneDegree);
      if (any(dropGeneNodes)) {
         if (verbose) {
            jamba::printDebug("subsetCnetIgraph(): ",
               "Dropping ",
               formatInt(sum(dropGeneNodes)),
               " gene nodes with less than degree:",
               minGeneDegree);
         }
         gCnet <- igraph::subgraph(gCnet,
            which(!dropGeneNodes));
      }
   }

   return(gCnet);
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
#' @param verbose logical indicating whether to print verbose output.
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

   for (ix in ixV) {
      ## make a color matrix with appropriate dimensions
      ##
      ## First check if any colors are blank colors
      vcc1 <- V(g)[ix]$coloredrect.color[[1]];
      vcc1blank <- isColorBlank(vcc1,
         blankColor=blankColor,
         c_max=c_max,
         l_min=l_min,
         alpha_max=alpha_max,
         ...);
      if (any(vcc1blank)) {
         if (constrain %in% "nrow") {
            vcc1m <- matrix(data=vcc1,
               byrow=V(g)[ix]$coloredrect.byrow,
               ncol=V(g)[ix]$coloredrect.ncol);
            #vcc1mBlank <- colMaxs(1*matrix(data=!vcc1 %in% blankColor,
            vcc1mBlank <- colMaxs(1*matrix(data=!vcc1blank,
               byrow=V(g)[ix]$coloredrect.byrow,
               ncol=V(g)[ix]$coloredrect.ncol))>0;
            vcc1m <- vcc1m[,vcc1mBlank,drop=FALSE];
         } else if (contrain %in% "ncol") {
            ## Constrain ncol
            vcc1 <- V(g)[ix]$coloredrect.color[[1]];
            vcc1m <- matrix(data=vcc1,
               byrow=V(g)[ix]$coloredrect.byrow,
               ncol=V(g)[ix]$coloredrect.ncol);
            vcc1mBlank <- rowMaxs(1*matrix(data=!vcc1blank,
               byrow=V(g)[ix]$coloredrect.byrow,
               ncol=V(g)[ix]$coloredrect.ncol))>0;
            vcc1m <- vcc1m[vcc1mBlank,,drop=FALSE];
         }
         if (!constrain %in% "none") {
            V(g)[ix]$coloredrect.ncol <- ncol(vcc1m);
            V(g)[ix]$coloredrect.nrow <- nrow(vcc1m);
            if (V(g)[ix]$coloredrect.byrow) {
               vcc1new <- as.vector(t(vcc1m));
            } else {
               vcc1new <- as.vector(vcc1m);
            }
            V(g)[ix]$coloredrect.color <- list(vcc1new);
         } else {
            #vcc <- setdiff(vcc1, blankColor);
            vcc <- vcc[!vcc1blank];
            if (length(vcc) > 0 && length(vcc) < length(vcc1)) {
               if (verbose) {
                  printDebug("removeCnetBlanks(): ",
                     ix);
               }
               V(g)[ix]$coloredrect.color <- list(vcc);
               if (V(g)[ix]$coloredrect.ncol == 1) {
                  V(g)[ix]$coloredrect.nrow <- length(vcc);
                  V(g)[ix]$coloredrect.ncol <- 1;
               } else {
                  V(g)[ix]$coloredrect.ncol <- length(vcc);
                  V(g)[ix]$coloredrect.nrow <- 1;
               }
            }
         }
      }
   }
   ## Now resize coloredrectangle size2 values
   ## so each square is constant size relative to
   ## its expected node size
   if (resizeNodes) {
      ## Make multi-segment gene nodes wider
      iG2 <- which(V(g)$coloredrect.ncol > 1 |
            V(g)$coloredrect.nrow > 1 &
            V(g)$shape %in% "coloredrectangle");
      maxNcol <- max(rmNA(unlist(V(g)[iG2]$coloredrect.ncol)));
      maxNrow <- max(rmNA(unlist(V(g)[iG2]$coloredrect.nrow)));
      if (length(iG2) > 0) {
         V(g)[iG2]$size2 <- (
            2 *
               V(g)[iG2]$size *
               V(g)[iG2]$coloredrect.ncol /
               maxNcol
         );
      }
   }

   ## TODO: iterate pie nodes
   if (applyToPie) {
      ## Adjust pie values
      V(g)$pie <- lapply(seq_along(V(g)$pie), function(i){
         iPie <- V(g)[[i]]$pie;
         iPieColor <- V(g)[[i]]$pie.color;
         iPie[!isColorBlank(iPieColor,
            blankColor=blankColor,
            c_max=c_max,
            l_min=l_min,
            alpha_max=alpha_max,
            ...)];
      });
      V(g)$pie.value <- lapply(seq_along(V(g)$pie.value), function(i){
         iPie <- V(g)[[i]]$pie.value;
         iPieColor <- V(g)[[i]]$pie.color;
         iPie[!isColorBlank(iPieColor,
            blankColor=blankColor,
            c_max=c_max,
            l_min=l_min,
            alpha_max=alpha_max,
            ...)];
      });
      ## Adjust pie colors
      V(g)$pie.color <- lapply(V(g)$pie.color, function(i){
         removeBlankColors(i,
            blankColor=blankColor,
            c_max=c_max,
            l_min=l_min,
            alpha_max=alpha_max,
            ...);
      });
   }

   return(g);
}

#' Determine if colors are blank colors
#'
#' Determine if colors are blank colors
#'
#' This function takes a vector of colors and determines if each color
#' is considered a "blank color", based upon direct match and the
#' color chroma saturation and luminance. For example, extremely pale
#' colors from `vals2colorLevels()` may be considered "blank" if the
#' color saturation is extremely low. Similarly, colors with
#' extremely high alpha transparency may be considered "blank".
#'
#' @param x character vector of R colors.
#' @param c_max maximum chroma as determined by HCL color space, in
#'    range of no color 0 to maximum color 100.
#' @param l_min numeric minimum luminance required for a color to be
#'    considered blank, combined with the `c_max` argument. This
#'    threshold prevents grey colors from being considered blank,
#'    unless their luminance is above this threshold.
#' @param alpha_max numeric value indicating the alpha transparency
#'    below which a color is considered blank, in range of fully
#'    transparent 0, to fully non-transparent 1.
#' @param blankColor character vector of R colors directly matched to
#'    the input `x` vector. The value `"transparent"` is useful here,
#'    because it is not easily converted to HCL color space.
#' @param ... additional arguments are ignored.
#'
#' @export
isColorBlank <- function
(x,
 c_max=7,
 l_min=95,
 alpha_max=0.1,
 blankColor=c("#FFFFFF","#FFFFFFFF","transparent"),
 ...)
{
   ## Purpose is to take a vector of colors and determine which are
   ## blank in terms of not having any color saturation, or being
   ## almost totally transparent.
   ##
   ## alpha_max is the highest alpha value considered "blank" regardless
   ## of all other color values
   ##
   ## c_max is the maximum HCL chroma value (usual range is 0 to 100)
   ## to be considered a blank color, combined with l_min
   ##
   ## l_min is the minimum HCL luminance (usual range is 0 to 100)
   ## to be considered a blank color, for example l_min=95 requires nearly
   ## white colors in order to be a blank color. To allow any greyscale
   ## color to be considered a blank color, use l_min=0, which imposes no
   ## luminance constraint.
   ##
   ## blankColors is a vector of colors, for example "transparent" is an
   ## allowable R color with alpha=0, but which cannot usually be converted to
   ## another colorspace.
   ##

   ## Handle missing values
   if (length(alpha_max) != 1) {
      alpha_max <- 0;
   }
   if (length(c_max) != 1) {
      c_max <- 0;
   }
   if (length(c_max) != 1) {
      l_min <- 0;
   }

   ## apply logic
   isBlank <- (is.na(x) |
         (tolower(x) %in% tolower(blankColor)) |
         (col2hcl(x)["C",] <= c_max & col2hcl(x)["L",] >= l_min) |
         col2alpha(x) <= alpha_max);
   names(isBlank) <- names(x);
   return(isBlank);
}

#' Fix Set labels for legibility
#'
#' Fix Set labels for legibility
#'
#' This function is a convenient wrapper for several steps that edit
#' gene set and pathways labels to be slightly more legible. It
#' operates on either a character vector, or an igraph object.
#'
#' @return vector or igraph object, to match the input `x`.
#'
#' @param x character vector, or `igraph` object. When an `igraph`
#'    object is supplied, the `V(g)$name` attribute is used as the
#'    basis of generating a label, which is then stored as `V(g)$label`.
#' @param wrap logical indicating whether to apply word wrap, based upon
#'    the supplied `width` argument.
#' @param width integer value used when `wrap=TRUE`, it is sent to
#'    `base::strwrap()`.
#' @param maxNchar numeric value or `Inf` to limit the maximum characters
#'    allowed for each string. This option is preferred when `wrap=TRUE`
#'    is not feasible, for example heatmap labels. When `NULL` or `Inf`
#'    no limit is applied. See `base::nchar()`.
#' @param suffix character value used as a suffix when `maxNchar` is used,
#'    the string is shortened so that the shortened string and suffix
#'    values are `maxNchar` characters long. It serves as an indicator
#'    that the string label has been shortened.
#' @param nodeType character value compared to the vertex attribute
#'    `"nodeType"` when the input `x` is an `igraph` object. This option
#'    is used to restrict label changes to certain nodes. When `NULL` or
#'    `nodeType="any"` then all node labels are updated.
#' @param removeGrep character regular expression pattern used to remove
#'    patterns from the resulting label. For example, `"^KEGG."` will remove
#'    the prefix `"KEGG_"` from all MSigDB KEGG pathways, in order to
#'    shorten the label. Due respect and credit to KEGG for their
#'    wonderful pathway data.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- c("KEGG_INSULIN_SIGNALING_PATHWAY",
#'    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#'    "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY");
#' fixSetLabels(x);
#'
#' jamba::nullPlot();
#' jamba::drawLabels(x, preset=c("top", "center", "bottom"));
#'
#' @export
fixSetLabels <- function
(x,
 wrap=TRUE,
 width=25,
 maxNchar=Inf,
 suffix="...",
 nodeType=c("Set","Gene","any"),
 removeGrep="^KEGG.",
 ...)
{
   if (igrepHas("igraph", class(x))) {
      xPrep <- jamba::ucfirst(
         tolower(
            gsub("_", " ",
               gsub(removeGrep,
                  "",
                  ignore.case=TRUE,
                  V(x)$name))));
   } else {
      xPrep <- ucfirst(
         tolower(
            gsub("_", " ",
               gsub(removeGrep,
                  "",
                  ignore.case=TRUE,
                  x))));
   }
   ## Optionally limit the character length
   if (length(maxNchar) > 0 && maxNchar < Inf) {
      if (length(suffix) == 0) {
         suffix <- "";
      }
      xLong <- (nchar(xPrep) > maxNchar);
      if (any(xLong)) {
         xPrep[xLong] <- paste0(
            substr(xPrep[xLong],
               1,
               maxNchar-nchar(suffix)),
            suffix);
      }
   }
   ## Optionally apply word wrap
   if (wrap) {
      xNew <- cPaste(sep="\n",
         doSort=FALSE,
         lapply(xPrep, function(i){
            strwrap(i, width=25);
         }));
   } else {
      xNew <- x;
   }
   ## Update the proper data to return
   if (igrepHas("igraph", class(x))) {
      if (length(nodeType) > 0 &&
            !"any" %in% nodeType &&
            "nodeType" %in% list.vertex.attributes(x)) {
         xUpdate <- which(V(x)$nodeType %in% nodeType);
      } else {
         xUpdate <- seq_len(vcount(x));
      }
      if (length(xUpdate) > 0) {
         V(x)[xUpdate]$label <- xNew;
      }
   } else {
      x <- xNew;
   }
   return(x);
}
