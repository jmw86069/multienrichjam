# jamenrich-base.r

.onLoad <- function
(libname,
 pkgname)
{
   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("coloredrectangle",
      clip=igraph::shape_noclip,
      plot=shape.coloredrectangle.plot);

   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("ellipse",
      clip=igraph::shape_noclip,
      plot=shape.ellipse.plot);

   ## define new igraph vertex shape "jampie"
   igraph::add_shape("jampie",
      clip=get("pie", envir=igraph:::.igraph.shapes)$clip,
      plot=shape.jampie.plot);
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
#' @param enrichDF `data.frame` representing gene set enrichment
#'    results.
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
#' @param msigdbGmtT optional gmtT object, a representation of
#'    `arules::transactions-class`. (Not currently implemented.)
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam conversion functions
#'
#' @importFrom jamba nameVector
#'
#' @export
enrichDF2enrichResult <- function
(enrichDF=NULL,
 pvalueCutoff=1,
 pAdjustMethod="none",
 keyColname=c("itemsetID", "ID", "Name", "Pathway"),
 pathGenes="pathGenes",
 geneColname=c("geneNames", "geneID", "Gene", "Genes"),
 geneHits="geneHits",
 geneRatioColname=c("GeneRatio", "^Ratio"),
 geneDelim="[,/ ]+",
 pvalueColname=c("P.Value", "Pvalue", "FDR", "adj.P.Val"),
 descriptionColname=c("Description", "Name", "Pathway", "ID"),
 msigdbGmtT=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert an enrichment data.frame
   ## into enrichResult class format usable by clusterProfiler
   ## methods, like enrichMap()
   if (suppressPackageStartupMessages(!require(clusterProfiler))) {
      stop("enrichDF2enrichResult() requires the clusterProfiler package.");
   }
   ## Find each colname in the input data.frame
   keyColname <- find_colname(keyColname, enrichDF);
   pathGenes <- find_colname(pathGenes, enrichDF);
   geneColname <- find_colname(geneColname, enrichDF);
   geneHits <- find_colname(geneHits, enrichDF);
   geneRatioColname <- find_colname(geneRatioColname, enrichDF);
   pvalueColname <- find_colname(pvalueColname, enrichDF);
   descriptionColname <- find_colname(descriptionColname, enrichDF);
   if (verbose) {
      jamba::printDebug("enrichDF2enrichResult(): ",
         "Colnames matched in the input data:",
         "\nkeyColname:", keyColname,
         "\npathGenes:", pathGenes,
         "\ngeneColname:", geneColname,
         "\ngeneHits:", geneHits,
         "\ngeneRatioColname:", geneRatioColname,
         "\npvalueColname:", pvalueColname,
         "\ndescriptionColname:", descriptionColname);
   }
   if (length(c(keyColname, pvalueColname, geneColname)) < 3) {
      stop("Could not find c(keyColname, pvalueColname, geneColname) in the input colnames.");
   }

   ## Confirm Description column
   if (length(descriptionColname) == 0) {
      warning(paste("No 'Description' colname was found in enrichDF,",
         "which prevents this data from being used by",
         "enrichplot and clusterProfiler functions.",
         "The Description column is recommended to contain",
         "the full name of the gene set.",
         sep="\n"));
   } else if (descriptionColname %in% c(keyColname)) {
      # If it is also the keyColname we need to make two columns
      # so both names can co-exist.
      enrichDF$Description <- enrichDF[[descriptionColname]];
   } else {
      ## Otherwise rename to make sure the final colname
      ## is 'Description' to fit expectations of
      ## enrichplot:::fortify.internal()
      enrichDF <- jamba::renameColumn(enrichDF,
         from=descriptionColname,
         to="Description");
      descriptionColname <- "Description";
   }

   enrichDF2 <- jamba::renameColumn(enrichDF,
      from=c(keyColname, pvalueColname, geneColname),
      to=c("ID", "pvalue", "geneID"));
   enrichDF2[,"p.adjust"] <- enrichDF2[,"pvalue"];


   ## Ensure all entries in column "ID" are unique
   ## because these values need to become rownames for compatibility
   ## with enrichplot::cnetplot()
   enrichDF2[["ID"]] <- jamba::makeNames(enrichDF2[["ID"]]);
   rownames(enrichDF2) <- enrichDF2[["ID"]];

   ## Convert gene delimiters all to "/"
   enrichDF2[["geneID"]] <- gsub(geneDelim,
      "/",
      enrichDF2[["geneID"]]);

   ## Validate input colnames
   keyColname <- intersect(keyColname, colnames(enrichDF));
   pathGenes <- intersect(pathGenes, colnames(enrichDF));
   geneColname <- intersect(geneColname, colnames(enrichDF));
   geneHits <- intersect(geneHits, colnames(enrichDF));
   geneRatioColname <- intersect(geneRatioColname, colnames(enrichDF));
   pvalueColname <- intersect(pvalueColname, colnames(enrichDF));

   if (length(geneRatioColname) > 0) {
      if (length(geneHits) > 0 && geneHits == geneRatioColname) {
         if (jamba::igrepHas("/", enrichDF2[[geneRatioColname]])) {
            if (verbose) {
               jamba::printDebug("enrichDF2enrichResult(): ",
                  "deriving ",
                  "geneHits",
                  " from gene ratio ",
                  geneRatioColname);
            }
            geneHits <- "geneHits";
            enrichDF2[[geneHits]] <- as.numeric(gsub("[/].*$", "",
               enrichDF2[,geneRatioColname]));
         } else {
            geneHits <- NULL;
         }
      }
      if (length(geneHits) == 0) {
         if (verbose) {
            jamba::printDebug("enrichDF2enrichResult(): ",
               "deriving ",
               "geneHits",
               " by splitting ",
               "geneID");
         }
         geneHits <- "geneHits";
         enrichDF2[[geneHits]] <- lengths(
            strsplit(
               as.character(enrichDF2[["geneID"]]),
               "/"));
      }
      if (length(pathGenes) == 0) {
         pathGenes <- "pathGenes";
         if (length(geneRatioColname) > 0) {
            if (verbose) {
               jamba::printDebug("enrichDF2enrichResult(): ",
                  "deriving ",
                  "pathGenes",
                  " from gene ratio ",
                  geneRatioColname);
            }
            if (jamba::igrepHas("/", enrichDF2[[geneRatioColname]])) {
               enrichDF2[[pathGenes]] <- as.numeric(gsub("^.*[/]", "",
                  enrichDF2[[geneRatioColname]]));
            } else if (all(jamba::rmNA(enrichDF2[[geneRatioColname]]) <= 1)) {
               enrichDF2[[pathGenes]] <- enrichDF2[[geneHits]] /
                  enrichDF2[[geneRatioColname]];
            } else {
               if (verbose) {
                  jamba::printDebug("enrichDF2enrichResult(): ",
                     "gene ratio did not contain '/' nor values with maximum 1, ",
                     "therefore no pathGenes values were created.");
               }
               pathGenes <- NULL;
            }
         }
      }
   } else {
      if (length(geneHits) == 0) {
         enrichDF2[["geneHits"]] <- lengths(strsplit(enrichDF2[["geneID"]], "/"));
         geneHits <- "geneHits";
      }
      if (length(pathGenes) == 0) {
         if (verbose) {
            jamba::printDebug("enrichDF2enrichResult(): ",
               "Assigning pathGenes == geneHits since no other information is available.");
         }
         enrichDF2[["pathGenes"]] <- enrichDF2[["geneHits"]];
         pathGenes <- "pathGenes";
      }
      if (verbose) {
         jamba::printDebug("enrichDF2enrichResult(): ",
            "deriving ",
            "'GeneRatio'",
            " from ",
            "geneHits/pathGenes");
      }
      geneRatioColname <- "GeneRatio";
      enrichDF2[,"GeneRatio"] <- jamba::pasteByRow(enrichDF2[,c(geneHits,pathGenes),drop=FALSE],
         sep="/");
   }

   enrichDF2 <- jamba::renameColumn(enrichDF2,
      from=c(geneRatioColname, pathGenes, geneHits),
      to=c("GeneRatio", "setSize", "Count"));
      #to=c("BgRatio", "setSize", "Count"));
   ## Make sure setSize contains integer values, in case we inferred the value
   if ("setSize" %in% colnames(enrichDF2)) {
      enrichDF2$setSize <- round(enrichDF2$setSize);
   }

   ## Re-order columns so "ID" is the first column
   if (verbose >= 2) {
      jamba::printDebug("enrichList2df(): ",
         "colnames(enrichDF2):", colnames(enrichDF2));
      jamba::printDebug("enrichList2df(): ",
         "class(enrichDF2):", class(enrichDF2));
      print(head(enrichDF2, 3));
   }
   keepcolids <- match(
      unique(jamba::provigrep(c("^ID$","."), colnames(enrichDF2))),
      colnames(enrichDF2));
   enrichDF2 <- enrichDF2[,keepcolids,drop=FALSE];
   #enrichDF2a <- dplyr::select(enrichDF2,
   #   dplyr::matches("^ID$"), tidyselect::everything());
   #enrichDF2 <- enrichDF2a;
   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
         "Done.");
   }

   gene <- jamba::mixedSort(unique(unlist(
      strsplit(
         as.character(enrichDF2[,"geneID"]),
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
      universe <- jamba::mixedSort(msigdbGmtT@itemInfo[,1]);
   } else {
      if (verbose) {
         jamba::printDebug("enrichDF2enrichResult(): ",
            "Defined geneSets from delimited gene values.");
      }
      geneSets <- strsplit(
         as.character(enrichDF2[["geneID"]]),
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
#' This function performs most of the work of comparing multiple
#' enrichment results. It takes a list of `enrichResult` objects,
#' generates an overall pathway-gene incidence matrix, assembles
#' a pathway-to-Pvalue matrix, creates EnrichMap `igraph` network
#' objects, and CnetPlot `igraph` network objects. It also applies
#' node shapes and colors consistent with the colors used for
#' each enrichment result.
#'
#' By default, each enrichment result table is subsetted for the
#' top `n=20` pathways sorted by pathway source, defined by
#' colnames `c("Source", "Category")`. For data without a source
#' column, the overall enrichment results are sorted to take the
#' top 20. Once the top 20 from each enrichment table are selected,
#' the overall set of pathways are used to retain these pathways
#' from all enrichment tables. In this way, a significant enrichment
#' result from one table will still be compared to a non-significant
#' result from another table.
#'
#' The default values for `topEnrichN` and related arguments
#' are intended when using enrichment results from `MSigDB`,
#' which has colnames `c("Source","Category")` and represents
#' 100 or more combinations of sources and categories. The
#' default values will select the top 20 entries from the
#' canonical pathways, after curating the canonical pathway
#' categories to one `"CP"` source value.
#'
#' To disable the top pathway filtering, set `topEnrichN=0`.
#'
#' Colors can be defined for each enrichment result using the
#' argument `colorV`, otherwise colors are assigned using
#' `colorjam::rainbowJam()`.
#'
#' @return `list` object containing various result formats:
#' * colorV: named vector of colors assigned to each enrichment,
#' where names match the names of each enrichment in `enrichList`.
#'
#' @param enrichList `list` of `enrichResult` objects, whose
#'    names are used in subsequent derived results.
#' @param geneHitList `list` of character vectors, or `NULL`. When
#'    `NULL` the gene hit list for each enrichment result is inferred
#'    from the enrichment results themselves. Therefore, if a gene
#'    does not appear in any enrichment result, it will not be
#'    used in this function.
#' @param colorV vector of colors, length equal to `length(enrichList)`,
#'    used to assign specific colors to each enrichment result.
#' @param nrow,ncol,byrow optional arguments used to customize
#'    `igraph` node shape `"coloredrectangle"`, useful when the
#'    number of `enrichList` results is larger than around 4. It
#'    defines the number of columns and rows used for each node,
#'    to display enrichment result colors, and whether to fill
#'    colors by row when `byrow=TRUE`, or by column when `byrow=FALSE`.
#' @param enrichLabels character vector of enrichment labels to use,
#'    as an optional alternative to `names(enrichList)`.
#' @param subsetSets character vector of optional set names to
#'    use in the analysis, useful to analyze only a specific subset
#'    of known pathways.
#' @param overlapThreshold numeric value between 0 and 1, indicating
#'    the Jaccard overlap score above which two pathways will be linked
#'    in the EnrichMap `igraph` network. By default, pathways whose
#'    genes overlap more than `0.1` will be connected, which is roughly
#'    equivalent to about a 10% overlap. Note that the Jaccard coefficient
#'    is adversely affected when pathway sets differ in size by more than
#'    about 5-fold.
#' @param cutoffRowMinP numeric value between 0 and 1, indicating the
#'    enrichment P-value required by at least one enrichment result, to
#'    be retained in downstream analyses. This P-value can be confirmed
#'    in the returned list element `"enrichIM"`, which is a matrix of
#'    P-values by pathway and enrichment.
#' @param enrichBaseline numeric value indicating the `-log10(P-value)`
#'    at which colors are defined as non-blank in color gradients.
#'    This value is typically derived from `cutoffRowMinP` to ensure
#'    that colors are only applied when a pathway meets this significance
#'    threshold.
#' @param enrichLens numeric value indicating the "lens" to apply to
#'    color gradients, where numbers above 0 make the color ramp more
#'    compressed, so colors are more vivid at lower numeric values.
#' @param enrichNumLimit numeric value indicating the `-log10(P-value)`
#'    above which each color gradient is considered the maximum color,
#'    useful to apply a fixed threshold for each color gradient.
#' @param nEM integer number, to define the maximum number of pathway
#'    nodes to include in the EnrichMap `igraph` network. This argument
#'    is passed to `enrichMapJam()`.
#' @param topEnrichN integer value with the maximum rows to retain
#'    from each `enrichList` table, by source. Set `topEnrichN=0` or
#'    `topEnrichN=NULL` to disable subsetting for the top rows.
#' @param topEnrichSources,topEnrichCurateFrom,topEnrichCurateTo,topEnrichSourceSubset,topEnrichDescriptionGrep,topEnrichNameGrep
#'    arguments passed to `topEnrichListBySource()` when `topEnrichN`
#'    is greater than `0`. The default values are used only when
#'    input data matches these patterns.
#' @param keyColname,nameColname,geneColname,pvalueColname,descriptionColname
#'    `character` vector in each case indicating the colnames
#'    for `key`, `name`, `gene`, `pvalue`, and `description`,
#'    respectively. Each vector is passed to `find_colname()` to determine
#'    the first suitable matching colname for each `data.frame` in
#'    `enrichList`.
#' @param descriptionCurateFrom,descriptionCurateTo `character` vectors
#'    with patterns and replacements, passed to `gsubs()`, intended to
#'    help curate common descriptions to shorter, perhaps more
#'    user-friendly labels. One example is removing the prefix
#'    `"Genes annotated by the GO term "` from Gene Ontology pathways.
#'    These label can be manually curated later, but it is often
#'    more convenient to curate them upfront in order to keep the
#'    different result objects consistent.
#' @param pathGenes,geneHits `character` values indicating the colnames
#'    that contain the number of pathway genes, and the number of gene
#'    hits, respectively.
#' @param geneDelim character pattern used with `strsplit()` to
#'    split multiple gene values into a list of vectors. The default
#'    for `enrichResult` objects is `"/"`, but the default for other
#'    sources is often `","`. The default pattern `"[,/ ]+"` splits
#'    by either `"/"`, `","`, or whitespace `" "`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to various internal
#'    functions.
#'
#' @examples
#' ## See the Vignette for a full walkthrough example
#'
#' @family jam enrichment functions
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
 enrichBaseline=-log10(cutoffRowMinP),
 enrichLens=0,
 enrichNumLimit=4,
 nEM=500,
 min_count=1,
 topEnrichN=20,
 topEnrichSources=c("gs_cat", "gs_subat"),
 topEnrichCurateFrom=NULL,
 topEnrichCurateTo=NULL,
 topEnrichSourceSubset=NULL,
 topEnrichDescriptionGrep=NULL,
 topEnrichNameGrep=NULL,
 keyColname=c("ID",
    "Name",
    "pathway",
    "itemsetID",
    "Description"),
 nameColname=c("Name",
    "pathway",
    "Description",
    "itemsetID",
    "ID"),
 geneColname=c("geneID",
    "geneNames",
    "Genes"),
 countColname=c("gene_count",
    "count",
    "geneHits"),
 pvalueColname=c("padjust",
    "p.adjust",
    "adjp",
    "padj",
    "qvalue",
    "qval",
    "q.value",
    "pvalue",
    "p.value",
    "pval",
    "FDR"),
 descriptionColname=c("Description",
    "Name",
    "Pathway",
    "ID"),
 descriptionCurateFrom=c("^Genes annotated by the GO term "),
 descriptionCurateTo=c(""),
 pathGenes=c("setSize",
    "pathGenes",
    "Count"),
 geneHits=c("Count",
    "geneHits",
    "gene_count"),
 geneDelim="[,/ ]+",
 GmtTname=NULL,
 #GmtTname="msigdbGmtTv50human",
 msigdbGmtT=NULL,
 #msigdbGmtT=msigdbGmtTv50human2,
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
   iDF1 <- head(enrichList[[1]]@result, 30);
   keyColname <- find_colname(keyColname, iDF1);
   geneColname <- find_colname(geneColname, iDF1);
   pvalueColname <- find_colname(pvalueColname, iDF1);
   descriptionColname <- find_colname(descriptionColname, iDF1);
   nameColname <- find_colname(nameColname, iDF1);
   countColname <- find_colname(countColname, iDF1);
   geneHits <- find_colname(geneHits, iDF1);
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "keyColname:", keyColname);
      jamba::printDebug("multiEnrichMap(): ",
         "nameColname:", nameColname);
      jamba::printDebug("multiEnrichMap(): ",
         "geneColname:", geneColname);
      jamba::printDebug("multiEnrichMap(): ",
         "pvalueColname:", pvalueColname);
      jamba::printDebug("multiEnrichMap(): ",
         "descriptionColname:", descriptionColname);
      jamba::printDebug("multiEnrichMap(): ",
         "countColname:", countColname);
   }

   ## Add some basic information
   if (length(enrichLabels) == 0) {
      enrichLabels <- jamba::nameVector(names(enrichList));
   } else if (length(names(enrichLabels)) == 0) {
      names(enrichLabels) <- names(enrichList);
   }
   mem$enrichLabels <- enrichLabels;

   ## Ensure that enrichBaseline is not greater than enrichNumLimit
   if (enrichBaseline >= enrichNumLimit) {
      enrichNumLimit <- enrichBaseline + 5;
   }
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
      print(jamba::sdim(enrichList));
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
      geneHitList <- enrichList2geneHitList(enrichList,
         geneColname=geneColname,
         geneDelim=geneDelim);
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim(geneHitList): ");
         print(jamba::sdim(geneHitList));
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
   ## Optionally run topEnrichBySource()
   if ((length(topEnrichN) > 0 && all(topEnrichN) > 0) ||
         length(subsetSets) > 0 ||
         length(topEnrichNameGrep) > 0 ||
         length(topEnrichDescriptionGrep) > 0) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "running topEnrichBySource().");
      }
      enrichList <- topEnrichListBySource(enrichList,
         n=topEnrichN,
         min_count=min_count,
         p_cutoff=cutoffRowMinP,
         sourceColnames=topEnrichSources,
         pvalueColname=pvalueColname,
         countColname=geneHits,
         curateFrom=topEnrichCurateFrom,
         curateTo=topEnrichCurateTo,
         sourceSubset=topEnrichSourceSubset,
         subsetSets=subsetSets,
         descriptionGrep=topEnrichDescriptionGrep,
         nameGrep=topEnrichNameGrep,
         verbose=verbose);
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim after topEnrichBySource():");
         print(jamba::sdim(enrichList));
      }
      ## Now we subset geneIM to the filtered genes
      origGenes <- jamba::mixedSort(unique(unlist(geneHitList)));
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitList) before:",
            lengths(geneHitList));
         jamba::printDebug("multiEnrichMap(): ",
            "length(origGenes):",
            length(origGenes));
      }
      geneHitListNew <- enrichList2geneHitList(enrichList,
         geneColname=geneColname,
         geneDelim=geneDelim);
      newGenes <- jamba::mixedSort(unique(unlist(geneHitListNew)));
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "length(newGenes):",
            length(newGenes));
         jamba::printDebug("multiEnrichMap(): ",
            "nrow(geneIM) before:",
            nrow(geneIM));
         jamba::printDebug("setdiff(rownames(geneIM), newGenes):", setdiff(rownames(geneIM), newGenes));
         jamba::printDebug("setdiff(newGenes, rownames(geneIM)):", setdiff(newGenes, rownames(geneIM)));
         jamba::printDebug("setdiff(newGenes, origGenes):", setdiff(newGenes, origGenes));
      }
      geneIM <- geneIM[rownames(geneIM) %in% newGenes,,drop=FALSE];
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "nrow(geneIM) after:",
            nrow(geneIM));
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitListNew):",
            lengths(geneHitListNew));
      }
      geneHitList <- lapply(geneHitList, function(i){
         intersect(i, newGenes)
      });
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitList) after:",
            lengths(geneHitList));
      }
   }

   #####################################################################
   ## geneIM colors
   geneIMcolors <- colorjam::matrix2heatColors(x=geneIM,
      transformFunc=c,
      colorV=colorV,
      shareNumLimit=FALSE,
      numLimit=1,
      trimRamp=c(4,3));
   mem$geneHitList <- geneHitList;
   mem$geneIM <- geneIM;
   mem$geneIMcolors <- geneIMcolors;

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
      jamba::printDebug("multiEnrichMap(): ",
         "enrichIM <- enrichList2IM() with    pvalueColname:",
         pvalueColname);
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichList[[1]]):");
      print(head(enrichList[[1]]));
   }
   enrichIM <- enrichList2IM(enrichList,
      valueColname=pvalueColname,
      keyColname=nameColname,
      #keyColname=keyColname,
      verbose=verbose,
      GmtT=msigdbGmtT)[enrichLsetNames,,drop=FALSE];
   if (all(is.na(enrichIM))) {
      jamba::printDebug("all enrichIM is NA.");
      jamba::printDebug("pvalueColname:", pvalueColname);
      jamba::printDebug("nameColname:", nameColname);
      stop("enrichIM is entirely NA.");
   }

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
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichIM):");
      print(head(enrichIM));
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
   if (!jamba::igrepHas("data.frame", class(enrichList[[1]]))) {
      geneCountColname <- head(jamba::provigrep(geneCountsGrep,
         colnames(as.data.frame(enrichList[[1]]))), 1);
   } else {
      geneCountColname <- head(jamba::provigrep(geneCountsGrep,
         colnames(enrichList[[1]])), 1);
   }
   if (length(geneCountColname) == 0) {
      stop(paste0("geneCountColname could not be found, by default it uses `geneHits`:",
         geneHits));
   }
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
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
      jamba::printDebug("multiEnrichMap(): ",
         "enrichIMcolors <- matrix2heatColors(enrichIMM)");
      jamba::printDebug("multiEnrichMap(): ",
         "colorV:", colorV,
         fgText=list("orange", "dodgerblue", colorV));
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichIM):");
      print(head(enrichIM));
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichIMM):");
      print(head(enrichIMM));
   }
   enrichIMcolors <- colorjam::matrix2heatColors(x=-log10(enrichIMM),
      colorV=colorV,
      lens=enrichLens,
      numLimit=enrichNumLimit,
      baseline=enrichBaseline,
      trimRamp=c(2, 2));
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "enrichBaseline:", enrichBaseline, ", colorV:", colorV);
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichIMcolors)");
      print(head(enrichIMcolors));
   }

   #####################################################################
   ## Subset for at least one significant enrichment P-value
   i1use <- rownames(enrichIMM)[(matrixStats::rowMins(enrichIMM[,useCols,drop=FALSE], na.rm=TRUE) <= cutoffRowMinP)];
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "nrow(enrichIM):",
         jamba::formatInt(nrow(enrichIM)),
         ", nrow(filtered for minimum P-value):",
         jamba::formatInt(length(i1use)));
      #ch(head(enrichIMM));
   }

   #####################################################################
   ## Now make sure enrichList only contains these sets
   if (nrow(enrichIM) > length(i1use)) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "dims before filtering minimum P-value():");
         print(jamba::sdim(enrichList));
      }
      enrichList <- lapply(enrichList, function(iDF){
         if (!jamba::igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         subset(iDF, iDF[[nameColname]] %in% i1use);
      });
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "dims after filtering minimum P-value():");
         print(jamba::sdim(enrichList));
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
      jamba::printDebug("multiEnrichMap(): ",
         "enrichDF <- enrichList2df(enrichList[c(",
         useCols,
         ")])");
      jamba::printDebug("multiEnrichMap(): ",
         "geneCountColname:",
         geneCountColname);
   }
   enrichDF <- enrichList2df(enrichList[useCols],
      msigdbGmtT=msigdbGmtT,
      keyColname=keyColname,
      geneColname=geneColname,
      geneCountColname=geneCountColname,
      pvalueColname=pvalueColname,
      geneDelim=geneDelim,
      verbose=verbose);

   #####################################################################
   ## Some cleaning of Description
   if (length(descriptionColname) == 1 &&
         descriptionColname %in% colnames(enrichDF)) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "cleaning Description column:",
            descriptionColname);
         #print(lengths(enrichList[useCols]));
      }
      descriptionColnameFull <- paste0(descriptionColname, "Full");
      enrichDF[,descriptionColnameFull] <- enrichDF[,descriptionColname];
      #enrichDF[,descriptionColname] <- memAdjustLabel(
      #   x=enrichDF[,descriptionColname],
      #   descriptionCurateFrom=descriptionCurateFrom,
      #   descriptionCurateTo=descriptionCurateTo);
   }

   #####################################################################
   ## Convert the combined enrichDF to enrichResult
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichDF):");
      print(head(as.data.frame(enrichDF)));
      jamba::printDebug("multiEnrichMap(): ",
         "enrichER <- enrichDF2enrichResult(), keyColname:",
         keyColname);
   }
   enrichER <- enrichDF2enrichResult(enrichDF,
      geneHits=geneHits,
      pathGenes=pathGenes,
      keyColname=keyColname,
      descriptionColname=descriptionColname,
      geneColname=geneColname,
      geneCountColname=geneCountColname,
      geneDelim=geneDelim,
      pvalueColname=pvalueColname,
      msigdbGmtT=msigdbGmtT,
      verbose=verbose);

   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichER):");
      print(head(as.data.frame(enrichER)));
   }

   mem$multiEnrichDF <- enrichDF;
   mem$multiEnrichResult <- enrichER;


   #####################################################################
   ## Incidence matrix of genes and pathways
   memIM <- list2im(
      keepCounts=TRUE,
      strsplit(
         jamba::nameVector(
            as.character(mem$multiEnrichDF[[geneColname]]),
            mem$multiEnrichDF[[nameColname]]),
         geneDelim));
   mem$memIM <- memIM;


   #####################################################################
   ## Convert enrichResult to enrichMap igraph network
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
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
   ## jamba::normScale(..., low=0) scales range 0 to maximum, into 0 to 1
   ## then add 0.3, then multiple by 8. Final range is 2.4 to 10.4
   igraph::V(enrichEM)$size_orig <- igraph::V(enrichEM)$size;
   igraph::V(enrichEM)$size <- (jamba::normScale(igraph::V(enrichEM)$size, low=0) + 0.3) * 8;
   igraph::E(enrichEM)$color <- "#99999977";

   mem$multiEnrichMap <- enrichEM;

   ## Convert EnrichMap to piegraph
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "running igraph2pieGraph() on enrichMap.");
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichIMcolors)");
      print(head(enrichIMcolors));
   }
   enrichEMpieUse <- igraph2pieGraph(g=enrichEM,
      defineLayout=FALSE,
      valueIMcolors=enrichIMcolors[i1use,useCols,drop=FALSE],
      verbose=verbose);

   ## Use colored rectangles
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "running rectifyPiegraph() on enrichMap.");
   }
   enrichEMpieUseSub2 <- rectifyPiegraph(enrichEMpieUse,
      nrow=nrow,
      ncol=ncol,
      byrow=byrow);
   igraph::V(enrichEMpieUseSub2)$size <- (jamba::normScale(igraph::V(enrichEMpieUseSub2)$size) + 0.3) * 6;
   igraph::V(enrichEMpieUseSub2)$size2 <- igraph::V(enrichEMpieUseSub2)$size2 / 2;
   mem$multiEnrichMap2 <- enrichEMpieUseSub2;

   #######################################################
   ## Create a CnetPlot
   ## Consider omitting this step if it is slow with large
   ## data, and if downstream workflows would typically
   ## only need a Cnet Plot on a subset of pathways and
   ## genes.
   gCt <- nrow(enrichER);
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "creating cnetPlot with cnetplotJam().");
   }
   gCnet <- cnetplotJam(enrichER,
      showCategory=gCt,
      categorySize=geneCountColname,
      doPlot=FALSE,
      nodeLabel=c(nameColname, descriptionColname, keyColname, "ID"),
      verbose=verbose);
   igraph::V(gCnet)$nodeType <- "Gene";
   igraph::V(gCnet)[seq_len(gCt)]$nodeType <- "Set";
   mem$multiCnetPlot <- gCnet;

   #######################################################
   ## Convert to coloredrectangle
   igraph::V(gCnet)[seq_len(gCt)]$name <- toupper(igraph::V(gCnet)[seq_len(gCt)]$name);
   ## Enrichment IM colors
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
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
      jamba::printDebug("multiEnrichMap(): ",
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
      jamba::printDebug("multiEnrichMap(): ",
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

   ## Add p_cutoff to output
   mem$p_cutoff <- cutoffRowMinP;

   return(mem);
}

#' Convert enrichList to IM incidence matrix
#'
#' Convert enrichList to IM incidence matrix
#'
#' This function takes a `list` of `enrichResult` objects
#' and creates an incidence matrix using the value defined
#' by `valueColname`.
#'
#' TODO: Port to use `venndir::list2im_value()`, which itself
#' may be moved to its own proper R package for set and list
#' manipulation, without the dependencies incurred by `venndir`.
#'
#' @family jam conversion functions
#'
#' @param enrichList `list` of `enrichResult` objects
#' @param GmtT not implemented here, experimental gene set object
#'    format that uses `arules::transactions` class, an efficient
#'    object with robust access functions in `arules`.
#' @param addAnnotations `logical` not implemented, this argument
#'    is paired with `GmtT`.
#' @param keyColname `character` used to match colnames, referring
#'    to the unique identifier. Values in this column will become
#'    the column headers in the resulting incidence matrix.
#' @param valueColname `character` used to match colnames to determine
#'    the value to place in each cell of the incidence matrix.
#' @param emptyValue `numeric` value used to fill empty cells in the
#'    incidence matrix. When `NULL` and `valueColname` contains
#'    "gene", "count", "num", "hit" then `emptyValue=0`, otherwise
#'    `emptyValue=1` is used with the assumption that `valueColname`
#'    refers to P-values.
#' @param verbose `logical` indicating whether to print verbose output.
#'
#' @export
enrichList2IM <- function
(enrichList,
 GmtT=NULL,
 addAnnotations=TRUE,
 keyColname=c("ID", "Name", "pathway", "itemsetID"),
 #keyColname="Name",
 valueColname=c("qvalue", "q.value", "pvalue", "p.value"),
 emptyValue=NULL,
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
      if (jamba::igrepHas("enrichResult", class(iDF))) {
         iDF <- iDF@result;
      }
      keyColname <- find_colname(keyColname, iDF);
      valueColname <- find_colname(valueColname, iDF);
      if (verbose) {
         jamba::printDebug("enrichList2IM(): ",
            "keyColname:", keyColname);
         jamba::printDebug("enrichList2IM(): ",
            "valueColname:", valueColname);
      }

      ## If "GeneRatio" then parse out the geneCount value
      if (jamba::igrepHas("GeneRatio", valueColname)) {
         iDF[,valueColname] <- gsub("[/].*$", "", iDF[,valueColname]);
         if (length(grep("^[0-9]*$", iDF[,valueColname])) == nrow(iDF)) {
            iDF[,valueColname] <- as.numeric(iDF[,valueColname]);
         }
      }
      if (verbose) {
         jamba::printDebug("enrichList2IM(): ",
            "head(iDF)");
         print(head(iDF));
      }
      jamba::nameVector(iDF[,c(valueColname,keyColname),drop=FALSE]);
   })));
   ## Empty values should be 1 instead of 0
   valueColname <- find_colname(valueColname, enrichList[[1]]);
   if (length(emptyValue) == 0) {
      if (jamba::igrepHas("gene|count|hits|num", valueColname)) {
         emptyValue <- 0;
      } else {
         emptyValue <- 1;
      }
   }
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
#' @family jam conversion functions
#'
#' @export
enrichList2df <- function
(enrichList,
 keyColname=c("itemsetID","ID","Name"),
 geneColname=c("geneNames", "Genes"),
 geneCountColname=c("geneCount", "geneHits", "Count"),
 pvalueColname=c("P-value", "pvalue", "Pval"),
 pvalueFloor=1e-200,
 msigdbGmtT=NULL,
 geneDelim="[,/ ]",
 verbose=FALSE,
 debug=0,
...)
{
   ## Purpose is to combine a list of enrichment data.frames into one data.frame
   if (!suppressPackageStartupMessages(require(matrixStats))) {
      stop("enrichList2df() requires the matrixStats package.");
   }
   iDF1 <- head(as.data.frame(enrichList[[1]]), 3);
   keyColname <- find_colname(keyColname, iDF1);
   geneColname <- find_colname(geneColname, iDF1);
   geneCountColname <- find_colname(geneCountColname, iDF1);
   pvalueColname <- find_colname(pvalueColname, iDF1);
   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
         "colnames(iDF1):",
         colnames(iDF1));
      jamba::printDebug("enrichList2df(): ",
         "keyColname:",
         keyColname);
      jamba::printDebug("enrichList2df(): ",
         "geneColname:",
         geneColname);
      jamba::printDebug("enrichList2df(): ",
         "geneCountColname:",
         geneCountColname);
      jamba::printDebug("enrichList2df(): ",
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
   if (!jamba::igrepHas("data.frame", class(head(jamba::rmNULL(enrichList), 1)))) {
      iDF <- as.data.frame(jamba::rmNULL(enrichList)[[1]]);
   } else {
      iDF <- jamba::rmNULL(enrichList)[[1]];
   }

   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
         "enrichCols (before):");
      print(enrichCols);
      jamba::printDebug("enrichList2df(): ",
         "colnames(iDF):", colnames(iDF));
   }
   enrichCols <- enrichCols[names(enrichCols) %in% colnames(iDF)];
   enrichColsHi <- names(enrichCols)[enrichCols %in% "hi"];
   #c("pathGenes","geneHits");
   enrichColsLo <- names(enrichCols)[enrichCols %in% "lo"];
   #enrichColsLo <- c("P-value");
   keepCols <- setdiff(jamba::unvigrep("gene", colnames(iDF)),
      c(enrichColsHi, enrichColsLo, keyColname, geneColname));

   ## Create a P-value incidence matrix
   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
         "enrichCols (after):");
      print(enrichCols);
      jamba::printDebug("enrichList2df(): ",
         "jamba::sdim(enrichL):");
      print(jamba::sdim(enrichList));
   }
   enrichValuesM <- do.call(cbind, lapply(jamba::nameVector(names(enrichCols)), function(iCol){
      useType <- enrichCols[iCol];
      enrichIMP <- list2imSigned(lapply(enrichList, function(iDF){
         if (!jamba::igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         if (useType %in% "lo" && any(iDF[,iCol] <= pvalueFloor)) {
            if (verbose) {
               jamba::printDebug("enrichList2df(): ",
                  "Some ", iCol, " values are less than ",
                  "pvalueFloor:",
                  pvalueFloor);
               print(table(iDF[,iCol] <= pvalueFloor));
            }
            iDF[iDF[,iCol] <= pvalueFloor, iCol] <- pvalueFloor;
         } else if (jamba::igrepHas("[/]", iDF[,iCol])) {
            iDF[,iCol] <- as.numeric(gsub("[/].*$", "", iDF[,iCol]));
         }
         if (length(jamba::tcount(iDF[[keyColname]], minCount=2)) > 0) {
            stop("enrichList2df(): There are duplicate values in iDF[[keyColname]], please resolve.");
         }
         if (verbose) {
            jamba::printDebug("enrichList2df(): ",
               "head(iDF):");
            print(head(iDF));
         }
         jamba::nameVector(iDF[,c(iCol,keyColname)]);
      }));
      if (useType %in% "lo") {
         enrichIMP[enrichIMP == 0] <- 1;
         jamba::nameVector(matrixStats::rowMins(enrichIMP), rownames(enrichIMP));
      } else if (useType %in% "hi") {
         jamba::nameVector(matrixStats::rowMaxs(enrichIMP), rownames(enrichIMP));
      }
   }));
   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
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
   allGenes <- jamba::mixedSort(unique(unlist(lapply(enrichList, function(iDF){
      if (!jamba::igrepHas("data.frame", class(iDF))) {
         iDF <- as.data.frame(iDF);
      }
      unlist(
         strsplit(
            as.character(iDF[,geneColname]),
            geneDelim));
   }))));

   ## If GmtT is supplied, use it to determine genes per pathway
   if (length(msigdbGmtT) > 0) {
      enrichGeneL <- as(msigdbGmtT[
         match(rownames(enrichValuesM), msigdbGmtT@itemsetInfo[,keyColname]),
         jamba::rmNA(match(allGenes, msigdbGmtT@itemInfo[,1]))], "list");
      names(enrichGeneL) <- rownames(enrichValuesM);
      enrichGeneVL <- list(jamba::cPaste(enrichGeneL, doSort=FALSE));
      names(enrichGeneVL) <- geneColname;
      enrichGeneLen <- lengths(enrichGeneL);
   } else {
      ## if GmtT is not supplied, use the pathway enrichment data as a substitute
      enrichL1L1 <- lapply(jamba::nameVectorN(enrichList), function(iName){
         iDF <- enrichList[[iName]];
         if (!jamba::igrepHas("data.frame", class(iDF))) {
            iDF <- as.data.frame(iDF);
         }
         iDF <- jamba::renameColumn(iDF,
            from=geneColname,
            to=iName);
         iDF[,c(keyColname,iName)];
      });
      enrichL1L <- jamba::mergeAllXY(enrichL1L1);
      enrichL1V <- jamba::nameVector(
         gsub("^[,/ ]+|[,/ ]+$",
            "",
            jamba::pasteByRow(
               enrichL1L[,-match(keyColname, colnames(enrichL1L)),drop=FALSE],
               sep=",")),
         enrichL1L[[keyColname]]);
      #enrichGeneL <- as.list(unique(
      #   CharacterList(
      #      strsplit(
      #         enrichL1V,
      #         geneDelim))));
      enrichGeneL <- strsplit(
         enrichL1V,
         geneDelim);
      enrichGeneVL <- list(jamba::cPaste(enrichGeneL, doSort=FALSE));
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
      jamba::printDebug("enrichList2df(): ",
         "head(enrichL1L):");
      print(str(head(enrichL1L)));
      jamba::printDebug("enrichList2df(): ",
         "dim(enrichValuesM):",
         dim(enrichValuesM));
   }
   keepColDF <- jamba::renameColumn(
      data.frame(row.names=rownames(enrichValuesM),
         keyColname=rep(NA, nrow(enrichValuesM))),
      from="keyColname",
      to=keyColname);
   for (iName in names(enrichList)) {
      iDF <- enrichList[[iName]];
      if (verbose) {
         jamba::printDebug("iName:", iName);
         jamba::printDebug("dim(iDF):", dim(iDF));
      }
      keyVals <- iDF[,keyColname];
      keyValsUse <- setdiff(keyVals, jamba::rmNA(keepColDF[,keyColname]));
      if (verbose) {
         jamba::printDebug("iName:", iName,
            ", length(keyVals):", length(keyVals),
            ", length(keyValsUse):", length(keyValsUse));
         jamba::printDebug("head(iDF):");
         print(head(iDF));
         jamba::printDebug("head(keepColDF):");
         print(head(keepColDF));
      }
      if (length(keyValsUse) > 0) {
         keepColDF[keyValsUse,keyColname] <- keyValsUse;
         for (keepCol in keepCols) {
            keyNew <- iDF[match(keyValsUse, iDF[,keyColname]),keepCol];
            if (length(keyNew) > 0) {
               keepColDF[keyValsUse,keepCol] <- keyNew;
            }
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
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "Creating enrichDF.");
   }
   enrichDF <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      enrichValuesM,
      keepColDF,
      as.data.frame(enrichGeneVL));
   enrichDF[["allGeneHits"]] <- enrichGeneLen;
   #whichCol1 <- max(which(colnames(enrichDF) %in% names(enrichCols))) + 1;
   #enrichDF <- insertDFcols(enrichDF, colnum=whichCol1,
   #   insertDF=data.frame(allGeneHits=enrichGeneLen));
   enrichDF;
}

#' Create enrichMap igraph object
#'
#' Create enrichMap igraph object
#'
#' This function is a minor extension to the original function
#' DOSE::enrichMap() which is now rewritten in the source package
#' to `enrichplot::emapplot()`. The major differences:
#'
#' * This function returns an `igraph` object, which can be manipulated
#' using network-related functions.
#' * This function calculates overlap using `dist(...,method="binary")`
#' which is a much faster method for calculating the Jaccard overlap.
#' * This function also calculates the overlap count, another helpful
#' measure for filtering network connections, for example to remove
#' links with only one gene, even if they overlap is above the
#' required threshold. Many spurious network connections are removed
#' with this filter, and it appears to be a helpful option.
#'
#' @return `igraph` object, whose nodes represent each enriched pathway,
#'    and are sized based upon the number of genes involved in the
#'    enrichment, and are colored based upon the `log10(Pvalue)`
#'    using `colorjam::vals2colorLevels()`, a function that applies
#'    a color gradient to a numeric range.
#'    Each edge has attributes: `overlap` containing Jaccard overlap,
#'    `overlap_count` with the number of genes in common between
#'    the two nodes, and `overlap_max_pct` with the maximum percent
#'    overlap between two nodes (overlap count)/(smaller node size).
#'
#' @param x either `enrichResult` or `data.frame` containing
#'    enrichment results, specifically expecting colnames to
#'    contain one of `c("ID","Description","Name")`
#'    to represent the node name, and `c("Description")` to represent
#'    the description, if present.
#' @param n numeric value indicating the maximum number of nodes to
#'    include in the final network.
#' @param vertex.label.font,vertex.label.cex attributes to define the
#'    default node label font and size.
#' @param keyColname,nodeLabel,descriptionColname character vectors
#'    indicating the colname to use for the node name and label.
#' @param nodeLabelFunc optional function to apply to `V(g)$name` in
#'    order to create `V(g)$label`. One suggestion is `fixSetLabels()`
#'    which applies word wrap, and optional max character length.
#' @param overlapThreshold numeric value indicating the minimum
#'    Jaccard overlap, where edges with lower values are deleted from
#'    the `igraph` object.
#' @param ... additional arguments are passed to `enrichDF2enrichResult()`
#'    when the input `x` is a `data.frame`.
#'
#' @family jam conversion functions
#' @family jam igraph functions
#'
#' @export
enrichMapJam <- function
(x,
 n=50,
 vertex.label.font=1,
 vertex.label.cex=1,
 keyColname="ID",
 nodeLabel=c("Name","Description","ID"),
 descriptionColname="Description",
 nodeLabelFunc=NULL,
 overlapThreshold=0.2,
 msigdbGmtT=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to customize enrichMap() to work with data
   ## generated outside clusterProfiler
   ##
   if (suppressPackageStartupMessages(!require(reshape2))) {
      stop("enrichMapJam() requires the reshape2 package is required for melt().");
   }
   #if (suppressPackageStartupMessages(!require(igraph))) {
   #   stop("enrichMapJam() requires the igraph package.");
   #}
   if (suppressPackageStartupMessages(!require(DOSE))) {
      stop("enrichMapJam() requires the DOSE package.");
   }
   if (is.null(nodeLabelFunc)) {
      nodeLabelFunc <- function(i){
         paste(collapse="\n",strwrap(width=30, jamba::ucfirst(gsub("_", " ", tolower(i)))));
      }
   }
   if (jamba::igrepHas("data.*frame", class(x))) {
      if (verbose) {
         jamba::printDebug("enrichMapJam(): ",
            "calling enrichDF2enrichResult()");
      }
      x <- enrichDF2enrichResult(x,
         msigdbGmtT=msigdbGmtT,
         verbose=verbose,
         ...);
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
      y <- head(y, n);
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
      igraph::V(g)$name <- y[, descriptionColname];
      igraph::V(g)$color <- "red";
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

      ## Jaccard coefficient is given as output from
      ## 1-dist(method="binary")
      wIM <- list2im(geneSets);
      w <- 1-as.matrix(dist(t(wIM), method="binary"));
      ## overlap counts, use sign() to count each gene only once
      wct <- t(sign(wIM)) %*% sign(wIM);
      ## min counts per cell
      wctmin <- wct;
      wctmin[] <- pmin(rep(diag(wct), ncol(wct)), rep(diag(wct), each=ncol(wct)));
      ## highest pct overlap
      wctmaxpct <- wct / wctmin;
      colnames(w) <- rownames(w) <- y[[nodeLabel]][match(colnames(w), y$ID)];

      wd <- reshape2::melt(w);
      wctd <- reshape2::melt(wct);
      wctmaxpctd <- reshape2::melt(wctmaxpct);

      wd1 <- match(wd[,1], colnames(w));
      wd2 <- match(wd[,2], colnames(w));
      w_keep <- (wd1 > wd2);
      wd <- wd[w_keep,,drop=FALSE];
      wctd <- wctd[w_keep,,drop=FALSE];
      wctmaxpctd <- wctmaxpctd[w_keep,,drop=FALSE];

      g <- igraph::graph.data.frame(wd[, -3], directed=FALSE);
      igraph::E(g)$width <- sqrt(wd[, 3] * 20);
      igraph::E(g)$overlap <- wd[,3];
      igraph::E(g)$overlap_count <- wctd[, 3];
      igraph::E(g)$overlap_max_pct <- wctmaxpctd[, 3];
      igraph::V(g)$pvalue <- pvalue[V(g)$name];

      ## Attempt to merge annotations from the enrichResult object
      iMatch <- match(igraph::V(g)$name, y[[nodeLabel]]);
      if (verbose) {
         jamba::printDebug("enrichMapJam(): ",
            "merging annotations from enrichResult objects");
      }
      if (!any(is.na(iMatch))) {
         #printDebug("match() worked with enrichResult data.frame Name colname.");
         iColnames <- jamba::unvigrep("^name$", colnames(y));
         for (iY in iColnames) {
            g <- g %>% igraph::set_vertex_attr(iY, value=y[iMatch,,drop=FALSE][[iY]]);
         }
      } else {
         ## Attempt to merge additional pathway annotation from GmtT
         #printDebug("match() worked with enrichResult data.frame Name colname.");
         if (length(msigdbGmtT) > 0) {
            iMatch <- match(igraph::V(g)$name, msigdbGmtT@itemsetInfo$Name);
            iMatchWhich <- which(!is.na(iMatch));
            if (length(iMatchWhich) > 0) {
               for (iCol1 in setdiff(colnames(msigdbGmtT@itemsetInfo), "Name")) {
                  g <- igraph::set_vertex_attr(g,
                     iCol1,
                     igraph::V(g)[iMatchWhich],
                     msigdbGmtT@itemsetInfo[iMatch[iMatchWhich],iCol1]);
               }
            }
         }
      }

      ## Apply optional node attributes
      if (length(vertex.label.font) > 0) {
         igraph::V(g)$label.font <- vertex.label.font;
      }
      if (length(vertex.label.cex) > 0) {
         igraph::V(g)$label.cex <- vertex.label.cex;
      }

      ## Delete edges where overlap is below a threshold
      g <- igraph::delete.edges(g, igraph::E(g)[igraph::E(g)$overlap < overlapThreshold]);

      pvalue <- igraph::V(g)$pvalue;

      nodeColor <- colorjam::vals2colorLevels(-log10(pvalue),
         col="Reds",
         numLimit=4,
         baseline=0,
         lens=2);
      igraph::V(g)$color <- nodeColor;

      if (is(x, "gseaResult")) {
         cnt <- jamba::nameVector(y$setSize, y[[nodeLabel]]);
      } else if (is(x, "enrichResult")) {
         if ("allGeneHits" %in% colnames(y)) {
            cnt <- jamba::nameVector(y$allGeneHits, y[[nodeLabel]]);
         } else {
            cnt <- jamba::nameVector(y$Count, y[[nodeLabel]]);
         }
      }
      cnt2 <- cnt[igraph::V(g)$name];
      ## Scale between 0 and 1, add 0.2 then multiply by 10
      ## largest node size will be 12 (1.2*10)
      ## smallest node size will be 2 dependent upon the difference between
      ##   largest and smallest gene count
      node_size <- (jamba::normScale(log10(cnt2+1) * 10, low=0) + 0.2) * 10;
      igraph::V(g)$size <- node_size;

      if (length(nodeLabelFunc) > 0 && is.function(nodeLabelFunc)) {
         igraph::V(g)$label <- sapply(igraph::V(g)$name, nodeLabelFunc);
      }
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
#' @family jam igraph functions
#'
#' @param gCnet igraph object representing Cnet concept network data
#' @param includeSets character vector, or NULL, containing the set
#'    names or labels to retain.
#' @param includeGenes character vector, or NULL, containing the gene
#'    names or labels to retain.
#' @param remove_singlets logical whether to remove singlet graph nodes,
#'    which are nodes that have no remaining edges. Note that
#'    argument `"removeSinglets"` is deprecated but will be recognized
#'    with preference, support will be removed in future versions.
#' @param minSetDegree integer value indicating the minimum number
#'    of edges each Set node must have to be retained in the resulting
#'    igraph. Use `minSetDegree=2` to retain only Set nodes that
#'    have multiple connected Gene nodes.
#' @param minGeneDegree integer value indicating the minimum number
#'    of edges each Gene node must have to be retained in the resulting
#'    igraph. Use `minGeneDegree=2` to retain only Gene nodes that
#'    connect to multiple Set nodes.
#' @param remove_blanks logical indicating whether to call
#'    `removeIgraphBlanks()`, which removes blank colors from
#'    `"pie.color"` and `"coloredrect.color"` attributes.
#' @param spread_labels logical indicating whether to call
#'    `spread_igraph_labels()`, which re-orients the node labels
#'    after the `igraph` object has been subsetted. When `TRUE`,
#'    the arguments `force_relayout`, and `do_reorder` are passed
#'    to that function.
#' @param force_relayout logical indicating whether to re-apply
#'    the `igraph` layout function.
#' @param do_reorder logical indicating whether to call
#'    `reorderIgraphNodes()`, which sorts equivalent nodes by
#'    color then label, intended when there are large numbers of
#'    nodes with the same edges, typically most common in a cnet
#'    plot where many gene nodes may be connected to the same
#'    pathway set nodes.
#' @param layout function that takes `igraph` object and returns a
#'    numeric matrix of node coordinates. This function is only
#'    called when `force_relayout=TRUE`, and must be supplied as
#'    a function in order to be applied properly to the subset
#'    Cnet `igraph`. To apply layout before the subset operation,
#'    do so with `igraph::set_graph_attr(g, "layout", layout)`.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @export
subsetCnetIgraph <- function
(gCnet,
 includeSets=NULL,
 includeGenes=NULL,
 remove_singlets=TRUE,
 minSetDegree=1,
 minGeneDegree=1,
 remove_blanks=FALSE,
 spread_labels=TRUE,
 force_relayout=TRUE,
 do_reorder=TRUE,
 repulse=4,
 layout=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take an Cnet igraph object and subset
   ## by set name or gene symbol
   ##########################################
   ## Check for deprecated arguments
   dots <- list(...);
   if (length(dots) > 0) {
      dep_args <- c(removeSinglets="remove_singlets");
      for (i in names(dep_args)) {
         if (i %in% names(dots)) {
            jamba::printDebug("subsetCnetIgraph(): ",
               "Deprecated argument '", i,
               "', please use '", dep_args[i], "'");
            assign(dep_args[i],
               dots[[i]]);
         }
      }
   }
   ##########################################
   ## Optionally subset for certain pathways
   if (length(includeSets) > 0) {
      if (length(igraph::V(gCnet)$label) == 0) {
         includeV <- which(igraph::V(gCnet)$nodeType %in% "Set" &
               (
                  toupper(igraph::V(gCnet)$name) %in% toupper(includeSets)
               ));
      } else {
         includeV <- which(igraph::V(gCnet)$nodeType %in% "Set" &
               (
                  toupper(igraph::V(gCnet)$name) %in% toupper(includeSets) |
                  toupper(igraph::V(gCnet)$label) %in% toupper(includeSets)
               ));
      }
      includeV2 <- unique(unlist(lapply(includeV, function(v){
         which(igraph::V(gCnet)$nodeType %in% "Gene" &
               igraph::V(gCnet)$name %in%
               names(igraph::neighbors(gCnet,
                  v=v,
                  mode="all")));
      })));
      includeVall <- sort(unique(c(includeV, includeV2)));
      if (verbose) {
         jamba::printDebug("subsetCnetIgraph(): ",
            "Filtered ",
            jamba::formatInt(sum(igraph::V(gCnet)$nodeType %in% "Set")),
            " Set nodes using ",
            jamba::formatInt(length(includeSets)),
            " includeSets down to ",
            jamba::formatInt(length(includeV)),
            " sets and ",
            jamba::formatInt(length(includeV2)),
            " genes in the Cnet igraph object.");
         whichNodeSets <- which(igraph::V(gCnet)$nodeType %in% "Set");
      }
      gCnet <- subgraph_jam(gCnet,
         includeVall);
   }

   ##########################################
   ## Optionally subset for certain genes
   if (length(includeGenes) > 0) {
      keepSetNodes <- which(igraph::V(gCnet)$nodeType %in% "Set");
      if (length(igraph::V(gCnet)$label) == 0) {
         keepGeneNodes <- which(
            igraph::V(gCnet)$nodeType %in% "Gene" &
               toupper(igraph::V(gCnet)$name) %in% toupper(includeGenes)
         );
      } else {
         keepGeneNodes <- which(
            igraph::V(gCnet)$nodeType %in% "Gene" &
               (toupper(igraph::V(gCnet)$name) %in% toupper(includeGenes) |
                     toupper(igraph::V(gCnet)$label) %in% toupper(includeGenes))
         );
      }
      keepNodes <- sort(unique(c(keepSetNodes, keepGeneNodes)));
      if (verbose) {
         jamba::printDebug("subsetCnetIgraph(): ",
            "Filtered ",
            jamba::formatInt(length(includeSets)),
            " includeGenes down to ",
            jamba::formatInt(length(keepGeneNodes)),
            " genes and ",
            jamba::formatInt(length(keepSetNodes)),
            " sets in the Cnet igraph object.");
      }
      gCnet <- subgraph_jam(gCnet,
         keepNodes);
   }

   #####################################################
   ## Optionally subset by degree of Set and Gene nodes
   if (length(minSetDegree) > 0) {
      dropSetNodes <- (igraph::V(gCnet)$nodeType %in% "Set" &
            igraph::degree(gCnet) < minSetDegree);
   } else {
      dropSetNodes <- rep(FALSE, igraph::vcount(gCnet));
   }
   if (length(minGeneDegree) > 0) {
      dropGeneNodes <- (igraph::V(gCnet)$nodeType %in% "Gene" &
            igraph::degree(gCnet) < minGeneDegree);
   } else {
      dropGeneNodes <- rep(FALSE, igraph::vcount(gCnet));
   }
   dropNodes <- (dropSetNodes | dropGeneNodes);
   if (any(dropNodes)) {
      if (verbose) {
         jamba::printDebug("subsetCnetIgraph(): ",
            "Dropping ",
            jamba::formatInt(sum(dropSetNodes)),
            " sets with fewer than ",
            minSetDegree,
            " connected genes, and ",
            jamba::formatInt(sum(dropGeneNodes)),
            " gene nodes with fewer than ",
            minGeneDegree,
            " connected sets.");
      }
      gCnet <- subgraph_jam(gCnet,
         which(!dropNodes));
   }

   #####################################################
   ## Polish the igraph by removing nodes with no edges
   if (remove_singlets) {
      gCnet <- removeIgraphSinglets(gCnet,
         min_degree=1,
         verbose=verbose);
   }

   #####################################################
   ## Optionally apply layout-related functions
   if (remove_blanks) {
      gCnet <- removeIgraphBlanks(gCnet,
         ...);
   }
   if (spread_labels) {
      gCnet <- spread_igraph_labels(gCnet,
         force_relayout=force_relayout,
         do_reorder=do_reorder,
         repulse=repulse,
         verbose=verbose,
         ...);
   } else {
      if (force_relayout) {
         if (length(layout) == 0) {
            gCnet <- relayout_with_qfr(gCnet,
               repulse=repulse,
               spread_labels=spread_labels,
               ...);
            layout <- igraph::graph_attr(gCnet, "layout");
            rownames(layout) <- igraph::V(gCnet)$name;
         } else if (is.function(layout)) {
            layout <- layout(gCnet);
         } else {
            layout <- layout[!dropNodes,,drop=FALSE];
         }
         if (nrow(layout) != igraph::vcount(gCnet)) {
            stop("layout dimensions do not match the subset cnet igraph.");
         }
         gCnet <- igraph::set_graph_attr(gCnet,
            "layout",
            layout);
      }
      if (do_reorder) {
         gCnet <- reorderIgraphNodes(gCnet,
            layout=NULL,
            verbose=verbose,
            ...);
      }
   }

   return(gCnet);
}

#' Determine if colors are blank colors
#'
#' Determine if colors are blank colors
#'
#' This function takes a vector of colors and determines if each color
#' is considered a "blank color", based upon direct match and the
#' color chroma saturation and luminance. For example, extremely pale
#' colors from `colorjam::vals2colorLevels()` may be considered "blank" if the
#' color saturation is extremely low. Similarly, colors with
#' extremely high alpha transparency may be considered "blank".
#'
#' @family jam utility functions
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
   ## handle x input as list
   x_names <- names(x);
   is_list <- FALSE;
   if ("list" %in% class(x)) {
      is_list <- TRUE;
      x_len <- lengths(x);
      ## Note unname(x) is used to drop parent names,
      ## and maintain child vector names if they exist.
      x <- unlist(unname(x));
   }

   ## apply logic
   isBlank <- (is.na(x) |
         (tolower(x) %in% tolower(blankColor)) |
         (jamba::col2hcl(x)["C",] <= c_max & jamba::col2hcl(x)["L",] >= l_min) |
         jamba::col2alpha(x) <= alpha_max);

   ## handle list input
   if (is_list) {
      isBlank <- split(isBlank,
         rep(seq_along(x_len), x_len));
   }

   names(isBlank) <- x_names;
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
#' @family jam igraph functions
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
#' @param adjustCase logical indicating whether to adjust the uppercase
#'    and lowercase lettering.
#' @param removeGrep character regular expression pattern used to remove
#'    patterns from the resulting label. The default values remove the
#'    prefix used in MsigDB canonical pathway names, which is a prefix
#'    indicating the source of each pathway.
#' @param words_from,words_to character vectors of words to match
#'    in case-insensitive manner, to be replaced with fixed-case
#'    alternatives. It uses perl-based regular expression matching
#'    in `base::gsub()`, and the `\\b` expression to enforce a
#'    word boundary, either via delimiter, whitespace, or the end
#'    of the string.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- c("KEGG_INSULIN_SIGNALING_PATHWAY",
#'    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#'    "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY");
#' fixSetLabels(x);
#'
#' jamba::nullPlot();
#' jamba::drawLabels(txt=x,
#'    preset=c("top", "center", "bottom"));
#'
#' @export
fixSetLabels <- function
(x,
 wrap=TRUE,
 width=25,
 maxNchar=Inf,
 suffix="...",
 nodeType=c("Set","Gene","any"),
 adjustCase=TRUE,
 removeGrep="^(KEGG|PID|REACTOME|BIOCARTA|NABA|SA|SIG|ST)[_.]",
 words_from=c("als", "ii", "iii", "iv", "v", "tgf",
    "nfkb", "trna", "rrna",
    "mirna", "mrna", "snrna", "snorna",
    "scrna", "lincrna"),
 words_to=c("ALS", "II", "III", "IV", "V", "TGF",
    "NFKB", "tRNA", "rRNA",
    "miRNA", "mRNA", "snRNA", "snoRNA",
    "scRNA", "lincRNA"),
 ...)
{
   nodeType <- match.arg(nodeType);
   if (jamba::igrepHas("igraph", class(x))) {
      if ("any" %in% nodeType) {
         which_nodes <- seq_len(igraph::vcount(x));
      } else {
         which_nodes <- which(igraph::V(x)$nodeType %in% nodeType);
      }
      xPrep <- gsub("_", " ",
         gsub(removeGrep,
            "",
            ignore.case=TRUE,
            igraph::V(x)$name[which_nodes]));
   } else {
      which_nodes <- seq_along(x);
      xPrep <- gsub("_", " ",
         gsub(removeGrep,
            "",
            ignore.case=TRUE,
            x));
   }
   if (adjustCase) {
      xPrep <- jamba::ucfirst(tolower(xPrep));
   }
   ## Optionally replace certain words with fixed capitalization
   if (length(words_from) > 0 && length(words_to) == length(words_from)) {
      for (i in seq_along(words_from)) {
         xPrep <- gsub(paste0("\\b", words_from[i], "\\b"),
            words_to[i],
            ignore.case=TRUE,
            perl=TRUE,
            xPrep);
      }
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
      xNew <- jamba::cPaste(sep="\n",
         doSort=FALSE,
         lapply(xPrep, function(i){
            strwrap(i, width=width);
         }));
   } else {
      xNew <- xPrep;
   }
   ## Update the proper data to return
   if (jamba::igrepHas("igraph", class(x))) {
      if (!"label" %in% igraph::list.vertex.attributes(x)) {
         igraph::V(x)$label <- igraph::V(x)$name;
      }
      igraph::V(x)[which_nodes]$label <- xNew;
   } else {
      x <- xNew;
   }
   return(x);
}



#' Subset enrichList for top enrichment results by source
#'
#' Subset enrichList for top enrichment results by source
#'
#' This function takes one `enrichResult` object, or
#' a `data.frame` of enrichment results, and determines the
#' top `n` number of pathways sorted by P-values, within
#' each pathway source. This function may optionally require
#' `min_count` genes in each pathway, and `p_cutoff` maximum
#' enrichment P-value, prior to taking the top `topEnrichN`
#' entries. The default arguments do not apply filters
#' to `min_count` and `p_cutoff`.
#'
#' When the enrichment data represents pathways from multiple
#' sources, the filtering and sorting is applied to each source
#' independently. The intent is to retain the top entries from
#' each source, as a method of representing each source
#' consistently even when one source may contain many more
#' pathways, and importantly where the range of enrichment P-values
#' may be very different for each source. For example, a database
#' of small canonical pathways would generally provide less
#' statistically significant P-values than a database of
#' dysregulated genes from gene expression experiments, where
#' each set contains a large number of genes.
#'
#' This function can optionally apply basic curation of pathway
#' source names, and can optionally be applied to multiple
#' source columns. This feature is intended for sources like
#' MSigDB (see http://software.broadinstitute.org/gsea/msigdb/index.jsp)
#' which contains columns `"Source"` and `"Category"`,
#' and where canonical pathways are either represented with `"CP"`
#' or a prefix `"CP:"`. The default parameters recognize this
#' case and curates all prefix `"CP:.*"` down to just `"CP"`
#' so that all canonical pathways are considered to be the
#' same source. For MSigDB there are also numerous other sources,
#' which are each independently filtered and sorted to the
#' top `topEnrichN` entries.
#'
#' Finally, this function is useful to subset enrichment results
#' by name, using `descriptionGrep` or `nameGrep`.
#'
#' @return `data.frame` subset up to `topEnrichN` rows, after
#'    applying optional `min_count` and `p_cutoff` filters.
#'
#' @family jam enrichment functions
#'
#' @param enrichDF `data.frame` containing enrichment results.
#' @param n `integer` maximum number of pathways to retain,
#'    after applying `min_count` and `p_cutoff` thresholds
#'    if relevant.
#' @param min_count `integer` minimum number of genes involved
#'    in an enrichment result to be retained, based upon values
#'    in `countColname`.
#' @param p_cutoff `numeric` value indicating the enrichment
#'    P-value threshold, pathways with enrichment P-value at
#'    or below this threshold are retained, based upon values
#'    in `pvalueColname`.
#' @param sourceColnames character vector of colnames in
#'    `enrichDF` to consider as the `"Source"`. Multiple
#'    columns will be combined using delimiter argument
#'    `sourceSep`. When `sourceColnames` is NULL or
#'    contains no `colnames(enrichDF)`, then data
#'    is considered `"All"`.
#' @param sortColname character vector indicating the colnames
#'    to use to sort data, prior to selecting the top `n`
#'    results by source. This argument is passed to
#'    `jamba::mixedSortDF(x, byCols=sortColname)`. Columns
#'    can be sorted in reverse order by using the prefix `"-"`,
#'    as described in `jamba::mixedSortDF()`.
#' @param countColname `character` vector of possible colnames
#'    in `enrichDF` that should contain the `integer` number
#'    of genes involved in enrichment. This vector is
#'    passed to `find_colname()` to find an appropriate
#'    matching colname in `enrichDF`.
#' @param pvalueColname `character` vector of possible colnames
#'    in `enrichDF` that should contain the enrichment P-value
#'    used for filtering by `p_cutoff`.
#' @param newColname new column name to use when `sourceColname`
#'    matches multiple colnames in `enrichDF`. Values for each
#'    row are combined using `jamba::pasteByRow()`.
#' @param curateFrom,curateTo character vectors with
#'    pattern,replacement values, passed to `gsubs()`
#'    to allow some editing of values. The default values
#'    convert MSigDB canonical pathways from the prefix `"CP:"`
#'    to use `"CP"` which has the effect of combining all
#'    canonical pathways before selecting the top `n` results.
#' @param sourceSubset character vector with a subset of
#'    sources to retain. If there are multiple colnames in
#'    `sourceColnames`, then column values are combined
#'    using `jamba::pasteByRow()` and delimiter `sourceSep`,
#'    prior to filtering.
#' @param sourceSep character string used as a delimiter
#'    when `sourceColnames` contains multiple colnames.
#' @param descriptionColname,nameColname character vectors
#'    indicating the colnames to consider description and name,
#'    as returned from `find_colname()`. These arguments are
#'    used only when `descriptionGrep` or `nameGrep` are
#'    supplied.
#' @param descriptionGrep,nameGrep character vector of patterns, used
#'    to filter pathways to those matching one or more patterns.
#'    This argument is used to help extract a specific subset
#'    of pathways of interest using keywords.
#'    The `descriptionGrep` argument searches only `descriptionColname`;
#'    the `nameGrep` argument searches only `nameColname`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
topEnrichBySource <- function
(enrichDF,
 n=15,
 min_count=1,
 p_cutoff=1,
 sourceColnames=c("gs_cat", "gs_subcat"),
 sortColname=c("P-value", "pvalue", "qvalue", "padjust", "-GeneRatio", "-Count", "-geneHits"),
 countColname=c("gene_count", "count", "geneHits"),
 pvalueColname=c("P.Value", "pvalue", "FDR", "adj.P.Val", "qvalue"),
 newColname="EnrichGroup",
 curateFrom=NULL,
 curateTo=NULL,
 sourceSubset=NULL,
 sourceSep="_",
 subsetSets=NULL,
 descriptionColname=c("Description", "Name", "Pathway"),
 nameColname=c("ID", "Name"),
 descriptionGrep=NULL,
 nameGrep=NULL,
 verbose=FALSE,
...)
{
   ## Purpose is to take a data.frame of pathway enrichment, split the results
   ## by values in sourceColnames, and return the top n rows from each
   ## group.
   ##
   ## Data is sorted using sortColname, prior to selecting the top n rows.
   ## To reverse the sort order, use a prefix "-" to the colname, e.g.
   ## sortColname <- c("P-value","-geneHits")
   ##
   ## If sourceSubset is supplied, it will return only results from
   ## source groups with these names, using names derived from:
   ## jamba::pasteByRow(sourceColnames, sep=sourceSep).
   ##
   ## If supplied, curateFrom and curateTo are a vector of gsub(from, to, ...)
   ## arguments, which are processed sequentially. The defaults will change
   ## "CP:REACTOME" and "CP:KEGG" to "CP". These substitutions should allow
   ## changing the sourceColnames values so they can be grouped in an
   ## appropriate way.
   ##
   ## newColname is the colname of a new column to add to the data.frame,
   ## indicating the sourceColname value used for grouping rows.
   ## If newColname is NULL, no new column is added to the data.frame.
   ##
   ## Results are returned as one data.frame, ordered by the sourceColnames
   ## groups, in order they appear in sourceSubset, or the order they appear
   ## after the split() function. That is, you must sort by P-value downstream
   ## if you wish data to be ordered by P-value.
   ##
   ## descriptionGrep is optionally used to pull out a subset of pathways
   ## whose name matches the grep pattern
   ##
   ## nameGrep is optionally used to pull out a subset of pathways
   ## whose name matches the grep pattern
   ##
   ## Tip: To see the available sourceColnames values, supply a false
   ## sourceSubset, which will cause an error message that prints the
   ## possible values.
   ##

   ## First convert enrichResult class to data.frame if needed
   enrichR <- NULL;
   if (jamba::igrepHas("enrichResult", class(enrichDF))) {
      if (length(rownames(enrichDF@result)) == 0) {
         rownames(enrichDF@result) <- paste0("row",
            jamba::padInteger(seq_len(nrow(enrichDF@result))));
      } else {
         rownames(enrichDF@result) <- jamba::makeNames(rownames(enrichDF@result));
      }
      nameColname <- find_colname(nameColname, enrichDF@result);
      if (!"Name" %in% nameColname) {
         enrichDF@result[,"Name"] <- enrichDF@result[[nameColname]];
      }
      enrichR <- enrichDF;
      enrichDF <- enrichDF@result;
   }

   ## Validate colnames
   sourceColnames <- find_colname(sourceColnames, enrichDF, max=Inf);
   descriptionColname <- find_colname(descriptionColname, enrichDF);
   nameColname <- find_colname(nameColname, enrichDF);
   countColname <- find_colname(countColname, enrichDF);
   pvalueColname <- find_colname(pvalueColname, enrichDF);
   if (!"Name" %in% nameColname) {
      enrichDF[,"Name"] <- enrichDF[[nameColname]];
   }

   if (verbose) {
      jamba::printDebug("topEnrichBySource(): ",
         "sourceColnames:", sourceColnames);
      jamba::printDebug("topEnrichBySource(): ",
         "descriptionColname:", descriptionColname);
      jamba::printDebug("topEnrichBySource(): ",
         "nameColname:", nameColname);
      jamba::printDebug("topEnrichBySource(): ",
         "sortColname:", sortColname);
      jamba::printDebug("topEnrichBySource(): ",
         "countColname:", countColname);
      jamba::printDebug("topEnrichBySource(): ",
         "pvalueColname:", pvalueColname);
   }

   ## Optionally determine column sort order
   if (length(sortColname) > 0) {
      if (!any(sortColname %in% colnames(enrichDF)) &&
            !any(gsub("^-", "", sortColname) %in% colnames(enrichDF))) {
         jamba::printDebug("topEnrichBySource(): ",
            "Warning: sortColname does not match colnames(enrichDF).");
         jamba::printDebug("topEnrichBySource(): ",
            "sortColname:", sortColname);
         jamba::printDebug("topEnrichBySource(): ",
            "colnames(enrichDF):", colnames(enrichDF));
      }
      enrichDF <- jamba::mixedSortDF(enrichDF,
         byCols=sortColname,
         ...);
   }

   ## Apply gene count filter
   if (length(countColname) > 0 && length(min_count) > 0 && min_count > 1) {
      if (verbose) {
         jamba::printDebug("topEnrichBySource(): ",
            c("Applying filter for min_count >= ",
               min_count,
               " retaining ",
               jamba::formatInt(sum(enrichDF[[countColname]] >= min_count)),
               " out of ",
               jamba::formatInt(nrow(enrichDF)),
               " total entries."),
            sep="");
      }
      enrichDF <- subset(enrichDF, enrichDF[[countColname]] >= min_count);
   }

   ## Apply pvalue filter
   if (length(pvalueColname) > 0 && length(p_cutoff) > 0 && p_cutoff < 1) {
      if (verbose) {
         jamba::printDebug("topEnrichBySource(): ",
            c("Applying filter for p_cutoff <= ",
               p_cutoff,
               " retaining ",
               jamba::formatInt(sum(enrichDF[[pvalueColname]] <= p_cutoff)),
               " out of ",
               jamba::formatInt(nrow(enrichDF)),
               " total entries."),
            sep="");
      }
      enrichDF <- subset(enrichDF, enrichDF[[pvalueColname]] <= p_cutoff);
   }

   if (length(sourceColnames) == 0) {
      iDFsplit <- rep("All", length.out=nrow(enrichDF));
      sourceSubset <- NULL;
   } else if (length(sourceColnames) == 1) {
      iDFsplit <- enrichDF[[sourceColnames]];
   } else {
      iDFsplit <- jamba::pasteByRow(sep=sourceSep,
         enrichDF[,sourceColnames,drop=FALSE]);
   }
   if (length(curateFrom) > 0) {
      iDFsplit <- gsubs(curateFrom,
         curateTo,
         iDFsplit,
         ...);
   }

   ## Optionally add grouping value to the input data.frame
   if (length(newColname) > 0) {
      newColname <- head(newColname, 1);
      enrichDF[[newColname]] <- iDFsplit;
   }

   if (verbose) {
      jamba::printDebug("topEnrichBySource(): ",
         "table(iDFsplit) before subsetting:");
      print(table(iDFsplit));
   }

   ## Split the data.frame in a list
   if (length(sourceSubset) == 0 || all(nchar(sourceSubset) == 0)) {
      sourceSubset <- levels(factor(iDFsplit));
   } else {
      if (!any(sourceSubset %in% iDFsplit)) {
         jamba::printDebug("topEnrichBySource(): ",
            "sourceSubset:",
            paste0("'", sourceSubset, "'"),
            " does not match values from sourceColnames:\n   ",
            paste0("'", levels(factor(iDFsplit)), "'"),
            sep="\n   ");
         print(sourceSubset);
         stop("topEnrichBySource() failed due to sourceSubset mismatch.");
      }
      sourceSubset <- intersect(sourceSubset, iDFsplit);
      ikeep <- (iDFsplit %in% sourceSubset);
      if (verbose) {
         jamba::printDebug("topEnrichBySource(): ",
            "sourceSubset: ", sourceSubset);
         jamba::printDebug("topEnrichBySource(): ",
            "table(ikeep): ");
         print(table(ikeep));
      }
      enrichDF <- subset(enrichDF, ikeep);
      iDFsplit <- iDFsplit[ikeep];
   }
   iDFsplitL <- split(enrichDF, iDFsplit);
   if (verbose) {
      jamba::printDebug("topEnrichBySource(): ",
         "sdim(iDFsplitL): ");
      print(jamba::sdim(iDFsplitL));
   }

   iDFtopL <- lapply(jamba::nameVectorN(iDFsplitL), function(iSubset){
      iDFsub <- iDFsplitL[[iSubset]];
      ## descriptionGrep
      if (length(descriptionColname) > 0 && length(descriptionGrep) > 0 && nrow(iDFsub) > 0) {
         descr_keep_vals <- jamba::provigrep(descriptionGrep,
            iDFsub[[descriptionColname]]);
         if (length(descr_keep_vals) > 0) {
            descr_keep <- (iDFsub[[descriptionColname]] %in% descr_keep_vals);
         } else {
            descr_keep <- rep(FALSE, nrow(iDFsub))
         }
      } else {
         descr_keep <- rep(TRUE, nrow(iDFsub))
      }
      ## nameGrep
      if (length(nameColname) > 0 && length(nameGrep) > 0 && nrow(iDFsub) > 0) {
         name_keep_vals <- jamba::provigrep(nameGrep,
            iDFsub[[nameColname]]);
         if (length(descr_keep_vals) > 0) {
            name_keep <- (iDFsub[[nameColname]] %in% name_keep_vals);
         } else {
            name_keep <- rep(FALSE, nrow(iDFsub))
         }
      } else {
         name_keep <- rep(TRUE, nrow(iDFsub))
      }
      ## subsetSets
      if (length(nameColname) > 0 && length(subsetSets) > 0) {
         subset_keep <- (tolower(iDFsub[[nameColname]]) %in% tolower(subsetSets));
      } else {
         subset_keep <- rep(TRUE, nrow(iDFsub))
      }
      rows_keep <- (descr_keep & name_keep & subset_keep);
      if (any(!rows_keep)) {
         iDFsub <- subset(iDFsub, rows_keep);
      }
      iDFtop <- head(iDFsub, n);
      iDFtop;
   });
   if (verbose) {
      jamba::printDebug("topEnrichBySource(): ",
         "sdim(iDFtopL): ");
      print(jamba::sdim(iDFtopL));
      print(head(iDFtopL[[1]], 5));
   }
   iDFnew <- data.frame(jamba::rbindList(unname(iDFtopL)),
      check.names=FALSE,
      stringsAsFactors=FALSE);
   if (verbose) {
      jamba::printDebug("topEnrichBySource(): ",
         "dim(iDFnew): ");
      print(dim(iDFnew));
      print(head(iDFnew, 5));
   }
   if (length(enrichR) > 0) {
      #ematch <- match(rownames(iDFnew), rownames(enrichR@result));
      ematch <- match(iDFnew$Name, rownames(enrichR@result));
      enrichR@result <- enrichR@result[ematch,,drop=FALSE];
      enrichR@geneSets <- enrichR@geneSets[ematch];
      return(enrichR)
   }
   iDFnew;
}


#' Subset enrichList for top enrichment results by source
#'
#' Subset enrichList for top enrichment results by source
#'
#' `topEnrichListBySource()` extends `topEnrichBySource()` by applying
#' filters to each `enrichList` entry, then keeping pathways
#' across all `enrichList` that match the filter criteria in any
#' one `enrichList`. It is most useful in the context of
#' `multiEnrichMap()` where a pathway must meet all criteria
#' in at least one enrichment, and that pathway should then
#' be included for all enrichments for the purpose of
#' comparative analysis.
#'
#' @rdname topEnrichBySource
#'
#' @family jam enrichment functions
#'
#' @param enrichList `list` of `enrichDF` entries, each passed
#'    to `topEnrichBySource()`.
#'
#' @export
topEnrichListBySource <- function
(enrichList,
 n=15,
 min_count=1,
 p_cutoff=1,
 sourceColnames=c("gs_cat", "gs_subcat"),
 sortColname=c("P-value", "pvalue", "qvalue", "padjust", "-GeneRatio", "-Count", "-geneHits"),
 countColname=c("gene_count", "count", "geneHits"),
 pvalueColname=c("P.Value", "pvalue", "FDR", "adj.P.Val", "qvalue"),
 newColname="EnrichGroup",
 curateFrom=NULL,
 curateTo=NULL,
 sourceSubset=NULL,
 sourceSep="_",
 subsetSets=NULL,
 descriptionColname=c("Description", "Name", "Pathway"),
 nameColname=c("ID", "Name"),
 descriptionGrep=NULL,
 nameGrep=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to extend topEnrichBySource() in an important way when
   ## dealing with a list of enrichment results, that it does one full
   ## pass through all enrichment results to determine individual
   ## pathways, then a second pass through each source using the shared
   ## consistent set of pathways.

   ## First create the subset for each enrichment result individually
   enrichLsub <- lapply(jamba::nameVectorN(enrichList), function(iName){
      if (verbose) {
         jamba::printDebug("topEnrichListBySource(): ",
            "iName:",
            iName);
      }
      iDF <- enrichList[[iName]];
      iDFsub <- topEnrichBySource(iDF,
         n=n,
         min_count=min_count,
         p_cutoff=p_cutoff,
         sourceColnames=sourceColnames,
         countColname=countColname,
         pvalueColname=pvalueColname,
         sortColname=sortColname,
         newColname=newColname,
         curateFrom=curateFrom,
         curateTo=curateTo,
         sourceSubset=sourceSubset,
         sourceSep=sourceSep,
         descriptionColname=descriptionColname,
         nameColname=nameColname,
         descriptionGrep=descriptionGrep,
         subsetSets=subsetSets,
         nameGrep=nameGrep,
         verbose=verbose,
         ...);
      if ("enrichResult" %in% class(iDFsub)) {
         eName <- iDFsub@result$Name;
      } else {
         eName <- iDFsub$Name;
      }
      if (verbose) {
         jamba::printDebug("topEnrichListBySource(): ",
            "   length(enrichNames):",
            jamba::formatInt(length(eName)));
      }
      eName;
   });
   enrichNames <- unique(unlist(enrichLsub));
   if (verbose) {
      jamba::printDebug("topEnrichListBySource(): ",
         "length(enrichNames):",
         jamba::formatInt(length(enrichNames)));
   }

   ## Step two, combine
   enrichLsubL <- lapply(jamba::nameVectorN(enrichList), function(iName){
      if (verbose) {
         jamba::printDebug("topEnrichListBySource(): ",
            "iName:",
            iName);
      }
      iDF <- enrichList[[iName]];
      if (jamba::igrepHas("enrichResult", class(iDF))) {
         nameColnameUse <- find_colname(nameColname, iDF@result);
         ematch <- jamba::rmNA(match(enrichNames,
            iDF@result[[nameColnameUse]]));
         iDF@result <- iDF@result[ematch,,drop=FALSE];
         iDF@geneSets <- iDF@geneSets[ematch];
      } else {
         nameColnameUse <- find_colname(nameColname, iDF);
         ematch <- match(enrichNames,
            iDF[[nameColnameUse]]);
         iDF <- iDF[ematch,,drop=FALSE];
      }
      iDF;
   });
   enrichLsubL;
}
