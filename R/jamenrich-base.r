# jamenrich-base.r


#' Prepare MultiEnrichMap data from enrichList
#'
#' Prepare MultiEnrichMap data from enrichList
#'
#' This function performs most of the work of comparing multiple
#' enrichment results.
#' This function takes a list of `enrichResult` objects,
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
#' @param enrichList `list` of `enrichResult` or `data.frame` objects.
#'    * The `names(enrichList)` are used in subsequent results.
#'    * Note that `data.frame` are converted to `enrichResult` using
#'    `enrichDF2enrichResult(x, ...)` where the '...' ellipses are
#'    used to recognize colnames.
#'    * Recommendation is to confirm each `data.frame` is properly
#'    converted to `enrichResult` upfront.
#' @param geneHitList `list` of character vectors, or
#'    `list` of `numeric` vectors whose names represent genes, or
#'    or `NULL`. When `NULL` the gene hit list for each enrichment
#'    result is inferred from the enrichment results themselves,
#'    however this option may incompletely represent which genes
#'    were statistical hits.
#'    Note that `geneHitList` and `geneHitIM` serve the same purpose
#'    and either can be supplied.
#' @param geneHitIM `numeric` matrix with gene rows, enrichment columns,
#'    and `numeric` values indicating the presence and/or direction
#'    of change for each gene.
#'    Note that `geneHitList` and `geneHitIM` serve the same purpose
#'    and either can be supplied.
#' @param colorV `character` vector of colors, length
#'    equal to `length(enrichList)`,
#'    used to assign specific colors to each enrichment result.
#' @param nrow,ncol,byrow optional arguments used to customize
#'    `igraph` node shape `"coloredrectangle"`, useful when the
#'    number of `enrichList` results is larger than around 4. It
#'    defines the number of columns and rows used for each node,
#'    to display enrichment result colors, and whether to fill
#'    colors by row when `byrow=TRUE`, or by column when `byrow=FALSE`.
#' @param enrichLabels `character` vector of enrichment labels to use,
#'    as an optional alternative to `names(enrichList)`.
#' @param subsetSets `character` vector of optional set names to
#'    use in the analysis, useful to analyze only a specific subset
#'    of known pathways.
#' @param overlapThreshold `numeric` (deprecated), value between 0 and 1,
#'    used to define the Enrichment Map, which is only created for legacy
#'    output enabled with `returnType="list"`.
#'    * To create a Multi-Enrichment Map `igraph` object, see `mem2emap()`.
#'    * overlapThreshold is the Jaccard overlap score above which
#'    two pathways will be linked in the resulting network.
#' @param p_cutoff `numeric` value between 0 and 1, default 0.05,
#'    enrichment P-value required by at least one enrichment result to
#'    be retained in downstream analyses.
#'    This P-value can be reviewed in `enrichIM(Mem)` of the output,
#'    which is a matrix of P-values by pathway and enrichment.
#'    The column header assigned as P-value is stored in
#'    `headers(Mem)$pvalueColname`.
#' @param cutoffRowMinP (deprecated in favor of 'p_cutoff').
#'    When it is non-NULL, it will be used for backward compatibility.
#' @param enrichBaseline `numeric` value, default uses `-log10(p_cutoff)`,
#'    the -log10 P-value threshold required for the color gradient to assign
#'    a color associated with enrichment.
#'    * In other words, P-values that do not meet this
#'    threshold are not colored in the color gradient for enrichment
#'    P-values.
#'    * To color all P-values use `enrichBaseline=0`.
#' @param enrichLens `numeric` value, default 0, indicating the "lens"
#'    to enhance the intensity of color gradients. Numbers above 0 make
#'    the color ramp more compressed, and more vivid at lower numeric values.
#' @param enrichNumLimit `numeric` value indicating the `-log10(P-value)`
#'    above which each color gradient is considered the maximum color,
#'    useful to apply a fixed threshold for each color gradient.
#' @param nEM `integer` (deprecated) maximum pathways to include in the
#'    Enrichment Map, which is only created for legacy
#'    output enabled with `returnType="list"`.
#'    * To create a Multi-Enrichment Map `igraph` object, see `mem2emap()`.
#' @param topEnrichN `integer` maximum rows to retain from each
#'    `enrichResult` in 'enrichList', for each source when supplied.
#'    Set `topEnrichN=0` or `topEnrichN=NULL` to retain all rows.
#'    Only rows where the pathway met other filtering criteria in at
#'    least one enrichment in 'enrichList' will be retained.
#' @param topEnrichSources,topEnrichCurateFrom,topEnrichCurateTo,topEnrichSourceSubset,topEnrichDescriptionGrep,topEnrichNameGrep
#'    arguments passed to `topEnrichListBySource()` when `topEnrichN`
#'    is greater than `0`. The default values are used only when
#'    input data matches these patterns.
#' @param keyColname,nameColname,geneColname,pvalueColname,descriptionColname
#'    `character` vector in each case with text strings or patterns to
#'    use when matching or prioritizing colnames to assign to each type:
#'    * `key`: The primary unique key for each pathway
#'    * `name`: The short name for each pathway
#'    * `gene`: column with delimited gene symbols or identifiers tested
#'    for enrichment of each pathway.
#'    * `pvalue`: column with enrichment P-value, typically prioritizing
#'    either 'qvalue' or 'p.adjust'. Per `enrichResult` convention,
#'    the selected column is renamed to 'pvalue'.
#'    * `description`: longer description associated with the pathway.
#'    
#'    Each vector is passed to `find_colname()` to find
#'    a suitable matching colname for each entry in `enrichList`.
#'    That function prioritizes full colname matches, then leading or
#'    trailing matches, then substring.
#' @param descriptionCurateFrom,descriptionCurateTo `character` vectors
#'    with patterns and replacements, passed to `gsubs()`, intended to
#'    help curate common descriptions to shorter, perhaps more
#'    user-friendly labels. One example is removing the prefix
#'    `"Genes annotated by the GO term "` from Gene Ontology pathways.
#'    These label can be manually curated later in the `Mem-class`
#'    methods, specifically `sets()<-` will allows assignment of custom
#'    names.
#' @param pathGenes,geneHits `character` values indicating the colnames
#'    that contain the number of pathway genes, and the number of gene
#'    hits, respectively. These values are optional, and not specifically
#'    used by multienrichjam.
#' @param geneDelim `character` pattern used with `strsplit()` to
#'    split multiple gene values into a list of vectors.
#'    The default (required) delimiter for `enrichResult` objects
#'    is `'/'`. A common alternative is `','` (comma-delimited).
#' @param verbose `logical` indicating whether to print verbose output.
#'    For `verbose` to cascade to internal functions, use `verbose=2`.
#' @param ... additional arguments are passed to various internal
#'    functions.
#'
#' @examples
#' ## See the Vignette for a full walkthrough example
#'
#' @family multienrichjam core functions
#' 
#' @importFrom lifecycle deprecated
#' 
#' @export
multiEnrichMap <- function
(enrichList,
 geneHitList=NULL,
 geneHitIM=NULL,
 colorV=NULL,
 nrow=NULL,
 ncol=NULL,
 byrow=FALSE,
 enrichLabels=NULL,
 subsetSets=NULL,
 overlapThreshold=deprecated(), # deprecated
 p_cutoff=0.05,
 cutoffRowMinP=deprecated(),
 enrichBaseline=-log10(p_cutoff),
 enrichLens=0,
 enrichNumLimit=4,
 nEM=500,
 min_count=3,
 topEnrichN=20,
 topEnrichSources=c("gs_cat", "gs_subcat"),
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
 directionColname=c("activation.z.{0,1}score",
    "z.{0,1}score"),
 direction_cutoff=0,
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
 returnType=c("Mem", "list"),
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

   returnType <- match.arg(returnType);

   ## Some data checks
   # if (suppressPackageStartupMessages(!require(igraph))) {
   #    stop("multiEnrichMap() requires the igraph package.")
   # }
   # if (suppressPackageStartupMessages(!require(DOSE))) {
   #    stop("multiEnrichMap() requires the DOSE package.")
   # }
   # if (suppressPackageStartupMessages(!require(IRanges))) {
   #    stop("multiEnrichMap() requires the IRanges package.")
   # }
   if (length(names(enrichList)) == 0) {
      stop("multiEnrichMap() requires names(enrichList).")
   }

   ## Convert list of data.frame to enrichResult
   enrichList <- lapply(jamba::nameVectorN(enrichList), function(iname){
      ier <- enrichList[[iname]];
      if (inherits(ier, "enrichResult")) {
         return(ier)
      }
      if (inherits(ier, "data.frame")) {
         if (verbose) {
            jamba::printDebug("multiEnrichMap(): ",
               "Converting '", iname, "' from data.frame")
         }
         ier <- enrichDF2enrichResult(ier,
            ...)
      }
   })

   ## Define some default colnames
   #nameColname <- "Name";
   #geneColname <- "geneNames";
   iDF1 <- head(enrichList[[1]]@result, 30);
   # version 0.0.56.900: change to use enrichList for any one
   # result to be non-NA, not just testing the first result
   keyColname <- find_colname(keyColname, enrichList);
   geneColname <- find_colname(geneColname, enrichList);
   pvalueColname <- find_colname(pvalueColname, enrichList);
   descriptionColname <- find_colname(descriptionColname, enrichList);
   nameColname <- find_colname(nameColname, enrichList);
   countColname <- find_colname(countColname, enrichList);
   directionColname <- find_colname(directionColname, enrichList);
   geneHits <- find_colname(geneHits, enrichList);
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
      jamba::printDebug("multiEnrichMap(): ",
         "directionColname:", directionColname);
   }

   ## Add some basic information
   if (length(enrichLabels) == 0) {
      enrichLabels <- jamba::nameVector(names(enrichList));
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
      print(jamba::sdim(enrichList));
   }

   #####################################################################
   ## Create geneHitList if not supplied
   ## Note: this step occurs before filtering gene sets,
   ## which may be incorrect. However, this order will
   ## ensure the geneIM is comprehensive, beyond just the
   ## genes involved in enrichment of the filtered subset.
   if (length(geneHitIM) > 0) {
      keep_hitim_colnames <- intersect(names(enrichList),
         colnames(geneHitIM));
      if (length(keep_hitim_colnames) == 0) {
         if (verbose) {
            jamba::printDebug("multiEnrichMap(): ",
               "geneHitIM had no matching colnames, setting", " NULL");
         }
         geneHitIM <- NULL;
      } else {
         geneHitIM <- geneHitIM[, keep_hitim_colnames, drop=FALSE];
      }
      # now populate geneHitList if needed
      if (length(geneHitList) == 0) {
         if (verbose) {
            jamba::printDebug("multiEnrichMap(): ",
               "Creating geneHitList from geneHitIM.");
         }
         # # Note: it saves character values not the direction here
         # geneHitList <- lapply(imSigned2list(geneHitIM), names);
         ## 0.0.104.900: keep the directionality
         geneHitList <- imSigned2list(geneHitIM)
      }
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim(geneHitList): ");
         print(jamba::sdim(geneHitList));
      }
   } else if (length(geneHitList) > 0) {
      # no geneHitIM but geneHitList is supplied
      if (!all(names(enrichList) %in% names(geneHitList))) {
         stop("names(enrichList) must all be present in names(geneHitList)");
      }
      keep_hitlist_names <- intersect(names(enrichList),
         names(geneHitList));
      geneHitList <- geneHitList[keep_hitlist_names];

      # populate geneHitIM if empty
      if (length(geneHitList) > 0 && length(geneHitIM) == 0) {
         if (is.numeric(geneHitList[[1]])) {
            # 0.0.92.900
            geneHitIM <- jamba::rmNA(naValue=0,
               list2imSigned(geneHitList, emptyValue=0))
         } else {
            geneHitIM <- list2im(geneHitList);
         }
      }
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim(geneHitList): ");
         print(jamba::sdim(geneHitList));
      }
   }
   # Finally if no geneHitList is available, derive from enrichment results
   if (length(geneHitList) == 0) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "creating geneHitList from enrichList.");
      }
      geneHitList <- enrichList2geneHitList(enrichList,
         geneColname=geneColname,
         geneDelim=geneDelim,
         verbose=(verbose - 1) > 0,
         make_unique=TRUE);
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim(geneHitList): ");
         print(jamba::sdim(geneHitList));
      }
      # also create geneHitIM for consistency
      geneHitIM <- list2im(geneHitList,
         sort_rows=jamba::mixedSort);
   }
   # At this point geneHitList and geneHitIM should both exist.
   # Note: There is no strict guarantee that geneHitIM contains all entries
   # in the enrichment results.
   # Might want to confirm in future.

   #####################################################################
   ## gene IM
   #
   # confirm no NA values exist, convert to zero 0
   if (any(is.na(geneHitIM))) {
      geneHitIM[is.na(geneHitIM)] <- 0;
   }
   geneIM <- geneHitIM;

   #####################################################################
   ## deprecated arguments
   if (lifecycle::is_present(cutoffRowMinP)) {
      lifecycle::deprecate_warn("0.0.101.900",
         "multienrichjam::multiEnrichMap(cutoffRowMinP)",
         "multienrichjam::multiEnrichMap(p_cutoff)");
      p_cutoff <- cutoffRowMinP;
   }
   if ("list" %in% returnType) {
      if (lifecycle::is_present(overlapThreshold)) {
         lifecycle::deprecate_warn("0.0.101.900",
            "multienrichjam::multiEnrichMap(overlapThreshold)");
      } else {
         overlapThreshold <- 0.1;
      }
   }
   ## Ensure that enrichBaseline is not greater than enrichNumLimit
   if (enrichBaseline >= enrichNumLimit) {
      enrichNumLimit <- enrichBaseline + 5;
   }

   #####################################################################
   ## store thresholds for reference
   thresholds <- list(
      min_count=min_count,
      p_cutoff=p_cutoff,
      direction_cutoff=direction_cutoff,
      topEnrichN=topEnrichN,
      topEnrichSources=topEnrichSources,
      topEnrichNameGrep=topEnrichNameGrep)
   if ("list" %in% returnType) {
      thresholds$overlapThreshold <- overlapThreshold;
   }

   #####################################################################
   ## Optionally run topEnrichBySource()
   if ((length(topEnrichN) > 0 && all(topEnrichN > 0)) ||
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
         p_cutoff=p_cutoff,
         sourceColnames=topEnrichSources,
         pvalueColname=pvalueColname,
         directionColname=directionColname,
         direction_cutoff=direction_cutoff,
         countColname=geneHits,
         curateFrom=topEnrichCurateFrom,
         curateTo=topEnrichCurateTo,
         sourceSubset=topEnrichSourceSubset,
         subsetSets=subsetSets,
         descriptionGrep=topEnrichDescriptionGrep,
         nameGrep=topEnrichNameGrep,
         verbose=(verbose - 1) > 0,
         ...);
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "sdim after topEnrichBySource():");
         print(jamba::sdim(enrichList));
      }
      ## Now we subset geneIM to the filtered genes
      origGenes <- jamba::mixedSort(unique(rownames(geneHitIM)));
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitList) before:",
            lengths(geneHitList));
         jamba::printDebug("multiEnrichMap(): ",
            "length(origGenes):",
            length(origGenes));
         jamba::printDebug("multiEnrichMap(): ",
            "head(origGenes):",
            head(origGenes));
      }
      geneHitListNew <- enrichList2geneHitList(enrichList,
         geneColname=geneColname,
         geneDelim=geneDelim,
         verbose=(verbose - 1) > 0,
         make_unique=TRUE);
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
      geneIM <- subset(geneIM, rownames(geneIM) %in% newGenes);
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "nrow(geneIM) after:",
            nrow(geneIM));
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitListNew):",
            lengths(geneHitListNew));
      }
      # 0.0.104.900: no need to subset geneHitList
      # geneHitList <- lapply(geneHitList, function(i){
      #    if (is.numeric(geneHit)) {
      #       i[names(i) %in% newGenes]
      #    } else {
      #       intersect(i, newGenes)
      #    }
      # });
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "lengths(geneHitList) after:",
            lengths(geneHitList));
      }
   }
   
   #####################################################################
   ## geneIMdirection (if geneHitIM is supplied)
   
   # 0.0.104.900: always prepare
   geneIMdirection <- (abs(geneIM) > 0) * 1;
   shared_rows <- intersect(rownames(geneHitIM),
      rownames(geneIMdirection));
   shared_columns <- intersect(colnames(geneHitIM),
      colnames(geneIMdirection));
   geneIMdirection[shared_rows, shared_columns] <- geneHitIM[shared_rows,
      shared_columns, drop=FALSE];
   #
   # check for missing values
   if (!all(rownames(geneIMdirection) %in% shared_rows)) {
      missing_rows <- setdiff(rownames(geneIMdirection), shared_rows);
      jamba::printDebug("multiEnrichMap(): ",
         "geneHitIM does not contain ",
         jamba::formatInt(length(missing_rows)),
         " rows present in geneIM, default values will use 1.");
      print(missing_rows);
   }
   if (!all(colnames(geneIMdirection) %in% shared_columns)) {
      missing_columns <- setdiff(colnames(geneIMdirection), shared_columns);
      jamba::printDebug("multiEnrichMap(): ",
         "geneHitIM does not contain ",
         jamba::formatInt(length(missing_columns)),
         " columns present in geneIM, default values will use 1.");
   }
   #
   # confirm geneIM only contains c(0,1)
   if (length(geneIMdirection) > 0) {
      geneIM <- (abs(geneIM) > 0) * 1;
   }

   #####################################################################
   ## geneIM colors
   ## Todo: Move this logic into geneIMcolors() Mem accessor
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "creating geneIMcolors");
   }
   # change this code to use actual colorV colors, no gradient necessary
   geneIMcolors <- do.call(cbind,
      lapply(jamba::nameVector(colnames(geneIM)), function(icolname){
         ifelse(geneIM[,icolname] == 0,
            "#FFFFFF",
            rep(colorV[icolname], nrow(geneIM)))
      }))
   # just to be sure, set rownames and colnames to match geneIM
   rownames(geneIMcolors) <- rownames(geneIM);
   colnames(geneIMcolors) <- colnames(geneIM);

   # assign to mem
   mem$geneHitList <- geneHitList;
   geneHitIMsorted <- jamba::mixedSort(rownames(geneHitIM));
   mem$geneHitIM <- geneHitIM[geneHitIMsorted, , drop=FALSE];
   geneIMsorted <- jamba::mixedSort(rownames(geneIM));
   mem$geneIM <- geneIM[geneIMsorted, , drop=FALSE];
   mem$geneIMcolors <- geneIMcolors[geneIMsorted, , drop=FALSE];
   mem$geneIMdirection <- geneIMdirection[geneIMsorted, , drop=FALSE];

   #####################################################################
   ## make nameColname values unique even when duplicated
   ##
   ## (Pathway Commons has several "Glycolysis" and "Gluconeogenesis")
   ## 0.0.101.900 - check for duplicated name
   ## - if it can be made unique with id:name then use that
   ##   to define a unique suffix number
   ## - if duplicated with no other ID then create unique suffix number
   ##   anyway.
   ##
   check_suffix <- " jam";
   namechecks <- lapply(enrichList, function(ier){
      dfcheck <- data.frame(
         key=ier@result[[keyColname]],
         name=ier@result[[nameColname]])
      dfcheck$keyname <- paste0(dfcheck$key,
         rep(":", length.out=nrow(dfcheck)),
         dfcheck$name)
      dupe_keynames <- unique(dfcheck$keyname[duplicated(dfcheck$keyname)]);
      if (length(dupe_keynames) > 0) {
         k <- dfcheck$keyname %in% dupe_keynames;
         dfcheck$name[k] <- jamba::makeNames(dfcheck$name[k],
            suffix=check_suffix);
      }
      dfcheck
   })
   # make one data.frame
   namecheck <- unique(jamba::rbindList(namechecks));
   # check unversioned duplicate name
   check_suffix_re <- paste0(check_suffix, "[0-9]+$")
   if (any(grepl(check_suffix_re, namecheck$name))) {
      # unversioned stem of versioned names
      versioned_names <- unique(gsub(check_suffix_re, "",
         grepv(check_suffix_re, namecheck$name)))
      k <- (namecheck$name %in% versioned_names);
      if (any(k)) {
         namecheck$name[k] <- paste0(namecheck$name[k],
            check_suffix, "1");
      }
      namecheck <- unique(namecheck);
   }
   
   if (any(duplicated(namecheck$name))) {
      # isolate duplicates to rename
      dupenames <- unique(namecheck$name[duplicated(namecheck$name)])
      dupedf <- jamba::mixedSortDF(
         subset(namecheck, name %in% dupenames),
         byCols=c("name", "key"))
      # create unique names
      dupedf$newname <- jamba::makeNames(dupedf$name,
         suffix=check_suffix)
      # dupedf$keyname <- paste0(dupedf$key, ":", dupedf$name)

      # apply names to enrichList
      enrichList <- lapply(enrichList, function(ier){
         if (nrow(ier@result) == 0) {
            return(ier)
         }
         ier_keyname <- paste0(ier@result[[keyColname]],
            ":",
            ier@result[[nameColname]])
         k <- match(ier_keyname, dupedf$keyname)
         if (any(!is.na(k))) {
            ier_newname <- ifelse(is.na(k),
               ier@result[[nameColname]],
               dupedf$newname[k])
            ier@result[[nameColname]] <- ier_newname;
         }
         ier
      })
   }
   
   #####################################################################
   ## enrichIM incidence matrix using -log10(P-value)
   ##
   ## Note: the rownames use values from descriptionColname
   ##
   enrichLsetNames <- unique(unlist(lapply(enrichList, function(iDF){
      if (!inherits(iDF, "data.frame")) {
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
      verbose=(verbose - 1) > 0,
      emptyValue=1,
      GmtT=msigdbGmtT);
   
   match1 <- match(enrichLsetNames, rownames(enrichIM));
   match2 <- match(names(enrichList), colnames(enrichIM));
   enrichIM <- enrichIM[match1, match2, drop=FALSE];
   rownames(enrichIM) <- enrichLsetNames;
   colnames(enrichIM) <- names(enrichList);

   if (all(is.na(enrichIM))) {
      if (verbose) {
         jamba::printDebug("all enrichIM is NA.");
         jamba::printDebug("pvalueColname:", pvalueColname);
         jamba::printDebug("nameColname:", nameColname);
      }
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
      keyColname=nameColname,
      valueColname=geneCountColname,
      emptyValue=0,
      verbose=(verbose - 1) > 0,
      GmtT=msigdbGmtT);
   match1 <- match(enrichLsetNames, rownames(enrichIMgeneCount));
   match2 <- match(names(enrichList), colnames(enrichIMgeneCount));
   enrichIMgeneCount <- enrichIMgeneCount[match1, match2, drop=FALSE];
   rownames(enrichIMgeneCount) <- enrichLsetNames;
   colnames(enrichIMgeneCount) <- names(enrichList);
   mem$enrichIMgeneCount <- as.matrix(enrichIMgeneCount);

   #####################################################################
   ## enrichIM incidence matrix using direction
   if (length(directionColname) > 0) {
      enrichIMdirection <- enrichList2IM(enrichList,
         keyColname=nameColname,
         valueColname=directionColname,
         emptyValue=0,
         verbose=(verbose - 1) > 0)[enrichLsetNames,,drop=FALSE];
      match1 <- match(enrichLsetNames, rownames(enrichIMdirection));
      match2 <- match(names(enrichList), colnames(enrichIMdirection));
      enrichIMdirection <- enrichIMdirection[match1, match2, drop=FALSE];
      rownames(enrichIMdirection) <- enrichLsetNames;
      colnames(enrichIMdirection) <- names(enrichList);
   } else {
      enrichIMdirection <- enrichIMgeneCount;
      enrichIMdirection[] <- 1;
   }
   mem$enrichIMdirection <- as.matrix(enrichIMdirection);

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
   # i1use <- rownames(enrichIMM)[(matrixStats::rowMins(enrichIMM[,useCols,drop=FALSE], na.rm=TRUE) <= cutoffRowMinP)];
   i1use <- rownames(enrichIMM)[(apply(enrichIMM[,useCols,drop=FALSE], 1, min, na.rm=TRUE) <= p_cutoff)];
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "nrow(enrichIM):",
         jamba::formatInt(nrow(enrichIM)),
         ", nrow(filtered for minimum P-value):",
         jamba::formatInt(length(i1use)));
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
      verbose=(verbose - 1) > 0);

   #####################################################################
   ## Some cleaning of Description
   # 0.0.101.900 - not used for now, will use sets(Mem)<- in future
   if (FALSE &&
         length(descriptionColname) == 1 &&
         descriptionColname %in% colnames(enrichDF)) {
      if (verbose) {
         jamba::printDebug("multiEnrichMap(): ",
            "cleaning Description column:",
            descriptionColname);
         #print(lengths(enrichList[useCols]));
      }
      descriptionColnameFull <- paste0(descriptionColname, "Full");
      enrichDF[,descriptionColnameFull] <- enrichDF[,descriptionColname];
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
      verbose=(verbose - 1) > 0);

   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "head(enrichER):");
      print(head(as.data.frame(enrichER)));
   }

   mem$multiEnrichDF <- enrichDF;
   mem$multiEnrichResult <- enrichER;


   #####################################################################
   ## Incidence matrix of genes and pathways
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "preparing memIM");
   }
   # memIM <- venndir::list2im_opt(do_sparse=FALSE,
   memIM <- list2im(
      keepCounts=TRUE,
      strsplit(
         jamba::nameVector(
            as.character(mem$multiEnrichDF[[geneColname]]),
            mem$multiEnrichDF[[nameColname]]),
         geneDelim),
      sort_rows=jamba::mixedSort);
   # 0.0.93.900 - enforce consistent order
   memIMsorted <- jamba::mixedSort(rownames(memIM));
   memIM <- memIM[memIMsorted, enrichLsetNames, drop=FALSE];
   mem$memIM <- memIM;

   #####################################################################
   ## Convert enrichResult to enrichMap igraph network
   # 0.0.101.900 - remove unless output is legacy mem list
   # in favor of using mem_plot_folio() or mem2emap()
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "converting enrichER to igraph enrichMap with enrichMapJam().");
   }
   if ("list" %in% returnType) {
      enrichEM <- enrichMapJam(enrichER,
         overlapThreshold=overlapThreshold,
         msigdbGmtT=msigdbGmtT,
         doPlot=FALSE,
         n=nEM,
         # keyColname="ID",
         keyColname=keyColname,
         nodeLabel=c(nameColname, descriptionColname, keyColname, "ID"),
         vertex.label.cex=0.5,
         verbose=verbose)
         # verbose=(verbose - 1) > 0);
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
         verbose=(verbose - 1) > 0);
   
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
   }

   #######################################################
   ## Create a CnetPlot
   ## Consider omitting this step if it is slow with large
   ## data, and if downstream workflows would typically
   ## only need a Cnet Plot on a subset of pathways and
   ## genes.
   #
   # 0.0.101.900 - omit this step unless output is legacy list
   # in favor of using mem_plot_folio() or mem2cnet()
   if ("list" %in% returnType) {
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
         verbose=(verbose - 1) > 0);
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
         verbose=(verbose - 1) > 0);
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
         verbose=(verbose - 1) > 0);
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
   }

   #######################################################
   ## Add all colnames to the mem object
   colnamesL <- list(
      keyColname=keyColname,
      nameColname=nameColname,
      geneColname=geneColname,
      countColname=countColname,
      pvalueColname=pvalueColname,
      descriptionColname=descriptionColname,
      directionColname=directionColname,
      pathGenes=pathGenes,
      geneHits=geneHits)
   mem$colnames <- colnamesL;

   ## Add p_cutoff to output
   mem$thresholds <- thresholds;
   # mem$p_cutoff <- cutoffRowMinP;

   if ("Mem" %in% returnType) {
      mem <- list_to_Mem(mem);
   }
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
#' TODO: Consider using `venndir::list2im_value()`, which itself
#' may be moved to its own proper R package for set and list
#' manipulation, without the dependencies incurred by `venndir`.
#'
#' @family jam conversion functions
#' 
#' @returns incidence `matrix` with `numeric` values using the
#'    the enrichment P-value (or FDR or Q-value as defined).
#'    Technically it is a "value incidence matrix" where non-NA values
#'    are derived from the source data, and NA represents absense of
#'    a value in the source data.
#'    It is possible to have a value '0' (pure zero) which represents
#'    presence of a value, where that value was '0', if the P-value
#'    reported by an enrichment tool rounded to zero.
#'
#' @param enrichList `list` of `enrichResult` objects
#' @param addAnnotations `logical` not implemented, this argument
#'    is paired with `GmtT`.
#' @param keyColname `character` used to match colnames, referring
#'    to the unique identifier. Values in this column will become
#'    the column headers in the resulting incidence matrix.
#' @param valueColname `character` used to match colnames to determine
#'    the value to place in each cell of the incidence matrix.
#' @param emptyValue `numeric` value used to fill empty cells in the
#'    incidence matrix, default NA.
#'    * NA is used by default, indicating absence of a value.
#'    * When `NULL` and `valueColname` contains
#'    "gene", "count", "num", "hit" then `emptyValue=0`.
#'    * Otherwise `emptyValue=1` is used with the assumption that
#'    `valueColname` refers to P-values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param GmtT (not currently implemented), alternative gene set object.
#'    In future it may enable passing additional annotations associated
#'    with each gene set, including description, source, category, etc.
#' @param ... additional arguments are ignored.
#' 
#' @export
enrichList2IM <- function
(enrichList,
 addAnnotations=TRUE,
 keyColname=c("ID",
    "Name",
    "pathway",
    "itemsetID"),
 valueColname=c("qvalue",
    "q.value",
    "padj",
    "pvalue",
    "p.value"),
 emptyValue=NA,
 verbose=FALSE,
 GmtT=NULL,
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
   ## Empty values should be 1 instead of 0
   valueColname1 <- find_colname(pattern=valueColname,
      x=enrichList,
      require_non_na=FALSE);
      # data.frame(check.names=FALSE,
       # enrichList[[1]]));
   if (length(emptyValue) == 0) {
      if (jamba::igrepHas("gene|count|hits|num|score|total|sum",
         valueColname1)) {
         emptyValue <- 0;
      } else {
         emptyValue <- 1;
      }
   }

   enrichIMP1 <- lapply(enrichList, function(iDF){
      if (inherits(iDF, "enrichResult")) {
         iDF <- iDF@result;
      }
      keyColname <- find_colname(keyColname,
         x=iDF,
         require_non_na=FALSE);
      valueColname <- find_colname(valueColname1,
         x=iDF,
         require_non_na=FALSE);
      if (verbose) {
         jamba::printDebug("enrichList2IM(): ",
            "keyColname:", keyColname,
            ", valueColname:", valueColname);
      }
      if (length(valueColname) == 0) {
         valueColname <- valueColname1;
         iDF[,valueColname] <- emptyValue;
      }

      ## If "GeneRatio" then parse out the geneCount value
      if (jamba::igrepHas("GeneRatio", valueColname)) {
         iDF[,valueColname] <- gsub("[/].*$", "", iDF[,valueColname]);
         if (length(grep("^[0-9]*$", iDF[,valueColname])) == nrow(iDF)) {
            iDF[,valueColname] <- as.numeric(iDF[,valueColname]);
         }
      }
      if (verbose > 1) {
         jamba::printDebug("enrichList2IM(): ",
            "head(iDF)");
         print(head(iDF[,unvigrep("geneID", colnames(iDF)), drop=FALSE]));
      }
      x <- jamba::nameVector(iDF[,c(valueColname,keyColname),drop=FALSE]);
      # if (any(is.na(x))) {
      # if (length(emptyValue) > 0) {
      #    x_blank <- (x %in% c(NA, ""));
      #    if (any(x_blank)) {
      #       x[x_blank] <- emptyValue;
      #    }
      # }
      x;
   });
   enrichIMP <- data.frame(check.names=FALSE,
      list2imSigned(enrichIMP1,
         emptyValue=emptyValue,
         sort_rows=jamba::mixedSort));

   if (nrow(enrichIMP) == 0) {
      return(enrichIMP);
   }

   # version 0.0.62.900: code below is no longer necessary
   # assuming there is nothing different in matrix format
   # that is not covered in the enrichIMP1 <- lapply() above
   #
   # if (length(emptyValue) > 0 && !emptyValue %in% 0) {
   #    enrichIMP[is.na(enrichIMP) | enrichIMP == 0] <- emptyValue;
   # }

   ## optional method not yet implemented for using GmtT objects
   # if (addAnnotations && length(GmtT) > 0) {
   #    ## Add information about pathways to the data.frame
   #    ematch <- match(rownames(enrichIMP), GmtT@itemsetInfo[,keyColname])
   #    enrichIMPinfo <- GmtT@itemsetInfo[ematch, , drop=FALSE];
   #    enrichIMP[,colnames(enrichIMPinfo)] <- enrichIMPinfo;
   # }
   return(enrichIMP);
}


#' Create enrichMap igraph object from enrichResult, deprecated
#'
#' Create enrichMap igraph object from enrichResult, deprecated,
#' use `mem2emap()`.
#'
#' This function is a minor customization to `enrichplot::emapplot()`,
#' which takes a single `enrichResult` and produces an Enrichment map
#' network. For the equivalent function using `Mem` as input, see
#' `mem2emap()`.
#' 
#' The major differences compared with `enrichplot`:
#'
#' * An `igraph` object is returned instead of `ggplot` object, since the
#' `igraph` object can be manipulated and is useful to review.
#' * This function calculates overlap using `dist(...,method="binary")`
#' which is a much faster method for calculating the Jaccard overlap.
#' * This function also calculates the overlap count, another helpful
#' measure for filtering network connections, for example to remove
#' links with only one gene, even if they overlap is above the
#' required threshold. Many spurious network connections are removed
#' with this filter, and it appears to be a helpful option.
#'
#' @returns `igraph` object with nodes representing each pathway,
#'    sized based upon the number of genes involved in enrichment, and
#'    colored based upon the `-log10(Pvalue)`
#'    using `colorjam::col_linear_xf()`, a function that applies
#'    a color gradient to a numeric range.
#'    Each edge has attributes: `overlap` containing Jaccard overlap,
#'    `overlap_count` with the number of genes in common between
#'    the two nodes, and `overlap_max_pct` with the maximum percent
#'    overlap between two nodes (overlap count)/(smaller node size).
#'
#' @param x `enrichResult` or `data.frame` containing
#'    enrichment results, specifically expecting colnames to
#'    contain one of `c("ID","Description","Name")`
#'    to represent the node name, and `c("Description")` to represent
#'    the description, if present.
#' @param n `numeric` value indicating the maximum number of nodes to
#'    include in the final network.
#' @param vertex.label.font,vertex.label.cex `numeric` values
#'    to define the default node label font and size.
#' @param keyColname,nodeLabel,descriptionColname `character` vectors
#'    indicating the colname to use for the node name and label.
#' @param nodeLabelFunc `function`, default NULL uses
#'    `fixSetLabels(x, width=30, ...)`. Use 'FALSE' to keep labels unchanged.
#'    Note that the `igraph` node name remains the same as input,
#'    and label is used as the visual label on the graph.
#'    * When providing a custom function, note that each node is called
#'    individually and therefore the call is not vectorized.
#' @param overlapThreshold `numeric` value indicating the minimum
#'    Jaccard overlap, where edges with lower values are deleted from
#'    the `igraph` object.
#' @param msigdbGmtT not currently implemented.
#' @param ... additional arguments are passed to `enrichDF2enrichResult()`
#'    when the input `x` is a `data.frame`.
#'
#' @family jam deprecated functions
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
   if (is.null(nodeLabelFunc)) {
      # 0.0.101.900 - use fixSetLabels for consistency, not custom function
      nodeLabelFunc <- function(x, width=30, ...)
         fixSetLabels(x, width=width, ...);
      # nodeLabelFunc <- function(i){
      #    paste(collapse="\n",
      #       strwrap(width=30,
      #          jamba::ucfirst(gsub("_", " ", tolower(i)))));
      # }
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
   if (inherits(x, "gseaResult")) {
      geneSets <- x@geneSets;
   } else if (inherits(x, "enrichResult")) {
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
      g <- igraph::make_empty_graph(0, directed=FALSE);
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
      wIM <- list2im(geneSets,
         sort_rows=jamba::mixedSort);
      # wIM <- venndir::list2im_opt(geneSets, do_sparse=FALSE);
      w <- (1 - as.matrix(dist(t(wIM), method="binary")));
      ## overlap counts, use sign() to count each gene only once
      wct <- t(sign(wIM)) %*% sign(wIM);
      ## min counts per cell
      wctmin <- wct;
      wctmin[] <- pmin(rep(diag(wct), ncol(wct)), rep(diag(wct), each=ncol(wct)));
      ## highest pct overlap
      wctmaxpct <- wct / wctmin;
      colnames(w) <- rownames(w) <- y[[nodeLabel]][match(colnames(w), y$ID)];

      ## 0.0.101.900 - replace reshape2 with tidyr?
      tidymelt <- function(x) {
         xdf <- data.frame(check.names=FALSE,
            tidyr::pivot_longer(
               data.frame(check.names=FALSE,
                  set=rownames(x), x),
               cols=colnames(x)))
         colnames(xdf)[1:2] <- c("Var1", "Var2");
         xdf
      }
      # wd <- reshape2::melt(w);
      # wctd <- reshape2::melt(wct);
      # wctmaxpctd <- reshape2::melt(wctmaxpct);
      wd <- tidymelt(w);
      wctd <- tidymelt(wct);
      wctmaxpctd <- tidymelt(wctmaxpct);
      
      wd1 <- match(wd[,1], colnames(w));
      wd2 <- match(wd[,2], colnames(w));
      w_keep <- (wd1 > wd2);
      wd <- subset(wd, w_keep);
      wctd <- subset(wctd, w_keep);
      wctmaxpctd <- subset(wctmaxpctd, w_keep);

      g <- igraph::graph_from_data_frame(wd[, -3, drop=FALSE], directed=FALSE);
      igraph::E(g)$width <- sqrt(wd[, 3] * 20);
      igraph::E(g)$overlap <- wd[, 3];
      igraph::E(g)$overlap_count <- wctd[, 3];
      igraph::E(g)$overlap_max_pct <- wctmaxpctd[, 3];
      igraph::V(g)$pvalue <- pvalue[igraph::V(g)$name];

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
            g <- igraph::set_vertex_attr(graph=g,
               name=iY,
               value=y[iMatch, ,drop=FALSE][[iY]]);
            # g <- g %>% igraph::set_vertex_attr(iY, value=y[iMatch,,drop=FALSE][[iY]]);
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
      g <- igraph::delete_edges(g, igraph::E(g)[igraph::E(g)$overlap < overlapThreshold]);

      pvalue <- igraph::V(g)$pvalue;

      # 0.0.101.900 - use colorjam::col_linear_df() which uses circlize
      # likely to be slightly less intense but more accurate
      nodeColor <- colorjam::col_linear_xf(4,
         lens=2,
         colramp="Reds")(-log10(pvalue))
      # nodeColor <- colorjam::vals2colorLevels(-log10(pvalue),
      #    col="Reds",
      #    numLimit=4,
      #    baseline=0,
      #    lens=2);
      igraph::V(g)$color <- nodeColor;

      if (is(x, "gseaResult")) {
         cnt <- jamba::nameVector(y$setSize, y[[nodeLabel]]);
      } else if (inherits(x, "enrichResult")) {
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

      if (is.function(nodeLabelFunc)) {
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
#' @param ... additional arguments are ignored.
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


