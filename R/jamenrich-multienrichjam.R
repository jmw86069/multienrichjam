

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
#' @param enrichList `list` of `enrichResult` objects, whose
#'    names are used in subsequent derived results.
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
#' @param overlapThreshold `numeric` value between 0 and 1, indicating
#'    the Jaccard overlap score above which two pathways will be linked
#'    in the EnrichMap `igraph` network. By default, pathways whose
#'    genes overlap more than `0.1` will be connected, which is roughly
#'    equivalent to about a 10% overlap. Note that the Jaccard coefficient
#'    is adversely affected when pathway sets differ in size by more than
#'    about 5-fold.
#' @param cutoffRowMinP `numeric` value between 0 and 1, indicating the
#'    enrichment P-value required by at least one enrichment result, to
#'    be retained in downstream analyses. This P-value can be confirmed
#'    in the returned list element `"enrichIM"`, which is a matrix of
#'    P-values by pathway and enrichment.
#' @param enrichBaseline `numeric` value indicating the `-log10(P-value)`
#'    at which colors are defined as non-blank in color gradients.
#'    This value is typically derived from `cutoffRowMinP` to ensure
#'    that colors are only applied when a pathway meets this significance
#'    threshold.
#' @param enrichLens `numeric` value indicating the "lens" to apply to
#'    color gradients, where numbers above 0 make the color ramp more
#'    compressed, so colors are more vivid at lower numeric values.
#' @param enrichNumLimit `numeric` value indicating the `-log10(P-value)`
#'    above which each color gradient is considered the maximum color,
#'    useful to apply a fixed threshold for each color gradient.
#' @param nEM `integer` number, to define the maximum number of pathway
#'    nodes to include in the EnrichMap `igraph` network. This argument
#'    is passed to `enrichMapJam()`.
#' @param topEnrichN `integer` value with the maximum rows to retain
#'    from each `enrichList` table, by source. Set `topEnrichN=0` or
#'    `topEnrichN=NULL` to disable subsetting for the top rows.
#' @param descriptionCurateFn `function` default `fixSetLabels()` used to
#'    curate pathway description to a user-friendly label.
#'    When `NULL` this step is skipped.
#' @param pathGenes,geneHits `character` values indicating the colnames
#'    that contain the number of pathway genes, and the number of gene
#'    hits, respectively.
#' @param geneDelim `character` pattern used with `strsplit()` to
#'    split multiple gene values into a list of vectors. The default
#'    for `enrichResult` objects is `"/"`, but the default for other
#'    sources is often `","`. The default pattern `"[,/ ]+"` splits
#'    by either `"/"`, `","`, or whitespace `" "`.
#' @param returnType `character` default "Mem" output class:
#'    * 'Mem' returns Mem S4 object
#'    * 'list' returns legacy list format (deprecated)
#' @param verbose `logical` indicating whether to print verbose output.
#'    For `verbose` to cascade to internal functions, use `verbose=2`.
#' @param ... additional arguments are passed to internal functions.
#'    * `topEnrichListBySource()` used to take the top `topEnrichN`
#'    pathways from each enrichment, and may also be used to subset
#'    by other criteria such as `descriptionGrep`, `nameGrep`,
#'    `sourceSubset`, `subsetSets`, etc.
#'
#' @examples
#' ## See the Vignette for a full walkthrough example
#'
#' @family jam enrichment functions
#'
#' @export
multienrichjam <- function
(enrichList,
 enrichLabels=NULL,
 p_cutoff=0.05,
 min_count=3,
 top_enrich_n=20,
 # topEnrichN=20,
 # cutoffRowMinP=0.05,
 geneHitList=NULL,
 geneHitIM=NULL,
 colorV=NULL,
 # direction_cutoff=0,   # used for topEnrichListBySource() filtering
 # overlapThreshold=0.1, # used for enrichMapJam(), defer to mem_plot_folio()
 # enrichBaseline=-log10(cutoffRowMinP), # used for color gradient color start
 # enrichLens=0, # used to adjust color gradient intensity
 enrichNumLimit=4, # color gradient max scale
 nEM=500,
 # subsetSets=NULL, ## used for topEnrichListBySource() filtering
 # topEnrichSources=c("gs_cat", "gs_subat"),
 # topEnrichCurateFrom=NULL,
 # topEnrichCurateTo=NULL,
 # topEnrichSourceSubset=NULL,
 # topEnrichDescriptionGrep=NULL,
 # topEnrichNameGrep=NULL,
 descriptionCurateFrom=c("^Genes annotated by the GO term "),
 descriptionCurateTo=c(""),
 descriptionCurateFn=fixSetLabels,
 geneDelim="[,/ ]+",
 returnType=c("Mem",
    "list"),
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

   # assign simple default names
   if (length(names(enrichList)) == 0) {
      ## do not use LETTERS to avoid clashing with pathway clusters later
      # names(enrichList) <- jamba::colNum2excelName(seq_along(enrichList))
      # use enrich01 format
      names(enrichList) <- paste0("enrich",
         jamba::padInteger(seq_along(enrichList)))
   }

   # validate input is list of enrichResult
   enrichList <- lapply(jamba::nameVectorN(enrichList), function(iname){
      ier <- enrichList[[iname]];
      if (inherits(ier, "enrichResult")) {
         ier;
      } else if (inherits(ier, "data.frame")) {
         if (verbose) {
            jamba::printDebug("multienrichjam(): ",
               "Converting '", iname, "' to enrichResult");
         }
         # convert to enrichResult as needed
         enrichDF2enrichResult(enrichDF=ier,
            keyColname=keyColname,
            pathGenes=pathGenes,
            geneColname=geneColname,
            geneHits=geneHits,
            geneDelim=geneDelim,
            pvalueColname=pvalueColname,
            descriptionColname=descriptionColname,
            ...)
      } else {
         # ??
      }
   })

   ## Define some default colnames
   #nameColname <- "Name";
   #geneColname <- "geneNames";
   iDF1 <- head(enrichList[[1]]@result, 30);
   # version 0.0.56.900: change to use enrichList for any one
   # result to be non-NA, not just testing the first result
   enrich_colnames <- find_enrich_colnames(enrichList,
      verbose=verbose,
      ...)
   keyColname <- enrich_colnames$keyColname
   geneColname <- enrich_colnames$geneColname
   pvalueColname <- enrich_colnames$pvalueColname
   descriptionColname <- enrich_colnames$descriptionColname
   nameColname <- enrich_colnames$nameColname
   countColname <- enrich_colnames$countColname
   directionColname <- enrich_colnames$directionColname
   geneHits <- enrich_colnames$geneHits

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
      jamba::printDebug("multienrichjam(): ",
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
            jamba::printDebug("multienrichjam(): ",
               "geneHitIM had no matching colnames, setting", " NULL");
         }
         geneHitIM <- NULL;
      } else {
         geneHitIM <- geneHitIM[, keep_hitim_colnames, drop=FALSE];
      }
   }
   if (length(geneHitList) == 0 && length(geneHitIM) > 0) {
      # generate geneHitList using geneHitIM
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "Creating geneHitList from geneHitIM.");
      }
      # Note: it saves character values not the direction here
      geneHitList <- lapply(imSigned2list(geneHitIM), names);
   }
   if (length(geneHitList) > 0) {
      if (!all(names(enrichList) %in% names(geneHitList))) {
         stop("names(enrichList) must all be present in names(geneHitList)");
      }
      keep_hitlist_names <- intersect(names(enrichList),
         names(geneHitList));
      if (length(keep_hitlist_names) == 0) {
         if (verbose) {
            jamba::printDebug("multienrichjam(): ",
               "geneHitList had no matching names, setting", " NULL");
         }
         geneHitList <- NULL;
      } else {
         geneHitList <- geneHitList[keep_hitlist_names];
      }

      # populate geneHitIM if empty
      if (length(geneHitList) > 0 && length(geneHitIM) == 0) {
         if (is.numeric(geneHitList[[1]])) {
            # 0.0.92.900
            geneHitIM <- jamba::rmNA(naValue=0, list2imSigned(geneHitList))
         } else {
            geneHitIM <- list2im(geneHitList);
         }
      }
   }
   # Finally if no geneHitList is available, use enrichment results
   if (length(geneHitList) == 0) {
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "creating geneHitList from enrichList.");
      }
      geneHitList <- enrichList2geneHitList(enrichList,
         geneColname=geneColname,
         geneDelim=geneDelim,
         verbose=(verbose - 1) > 0,
         make_unique=TRUE);
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "sdim(geneHitList): ");
         print(jamba::sdim(geneHitList));
      }
      # also create geneHitIM for consistency
      geneHitIM <- list2im(geneHitList);
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
   ## store thresholds for reference
   thresholds <- list(
      min_count=min_count,
      p_cutoff=p_cutoff,
      top_enrich_n=top_enrich_n)
   # Note: additional thresholds may be added by other tools later
      # overlapThreshold # enrichMap Jaccard overlap threshold
      # cutoffRowMinP # deprecated for p_cutoff
      # topEnrichN    # deprecated for top_enrich_n
      # topEnrichSources  # deprecated for topEnrichListBySource() via '...'
      # topEnrichNameGrep # deprecated as above

   #####################################################################
   ## Optionally run topEnrichBySource()
   if (length(top_enrich_n) > 0 && all(top_enrich_n > 0)) {
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "running topEnrichBySource().");
      }
      enrichList <- topEnrichListBySource(enrichList,
         n=top_enrich_n,
         min_count=min_count,
         p_cutoff=p_cutoff,
         # sourceColnames=topEnrichSources,
         # pvalueColname=pvalueColname,
         # directionColname=directionColname,
         # direction_cutoff=direction_cutoff,
         # countColname=geneHits,
         # curateFrom=topEnrichCurateFrom,
         # curateTo=topEnrichCurateTo,
         # sourceSubset=topEnrichSourceSubset,
         # subsetSets=subsetSets,
         # descriptionGrep=topEnrichDescriptionGrep,
         # nameGrep=topEnrichNameGrep,
         verbose=(verbose - 1) > 0,
         ...);
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "sdim after topEnrichBySource():");
         print(jamba::sdim(enrichList));
      }
      ## Now we subset geneIM to the filtered genes
      origGenes <- jamba::mixedSort(unique(unlist(geneHitList)));
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "lengths(geneHitList) before:",
            lengths(geneHitList));
         jamba::printDebug("multienrichjam(): ",
            "length(origGenes):",
            length(origGenes));
         jamba::printDebug("multienrichjam(): ",
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
         jamba::printDebug("multienrichjam(): ",
            "length(newGenes):",
            length(newGenes));
         jamba::printDebug("multienrichjam(): ",
            "nrow(geneIM) before:",
            nrow(geneIM));
         jamba::printDebug("setdiff(rownames(geneIM), newGenes):",
            setdiff(rownames(geneIM), newGenes));
         jamba::printDebug("setdiff(newGenes, rownames(geneIM)):",
            setdiff(newGenes, rownames(geneIM)));
         jamba::printDebug("setdiff(newGenes, origGenes):",
            setdiff(newGenes, origGenes));
      }
      geneIM <- subset(geneIM, rownames(geneIM) %in% newGenes);
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "nrow(geneIM) after:",
            nrow(geneIM));
         jamba::printDebug("multienrichjam(): ",
            "lengths(geneHitListNew):",
            lengths(geneHitListNew));
      }
      geneHitList <- lapply(geneHitList, function(i){
         intersect(i, newGenes)
      });
      if (verbose) {
         jamba::printDebug("multienrichjam(): ",
            "lengths(geneHitList) after:",
            lengths(geneHitList));
      }
   }

   #####################################################################
   ## geneIMdirection (if geneHitIM is supplied)
   geneIMdirection <- NULL;
   if (length(geneHitIM) > 0 &&
         (any(geneHitIM < 0) || !all(geneHitIM %in% c(0, 1)))) {
      geneIMdirection <- (abs(geneIM) > 0) * 1;
      shared_rows <- intersect(rownames(geneHitIM),
         rownames(geneIMdirection));
      shared_columns <- intersect(colnames(geneHitIM),
         colnames(geneIMdirection));
      geneIMdirection[shared_rows, shared_columns] <- geneHitIM[shared_rows, shared_columns];
      if (!all(rownames(geneIMdirection) %in% shared_rows)) {
         missing_rows <- setdiff(rownames(geneIMdirection), shared_rows);
         jamba::printDebug("multienrichjam(): ",
            "geneHitIM does not contain ",
            jamba::formatInt(length(missing_rows)),
            " rows present in geneIM, default values will use 1.");
         print(missing_rows);
      }
      if (!all(colnames(geneIMdirection) %in% shared_columns)) {
         missing_columns <- setdiff(colnames(geneIMdirection), shared_columns);
         jamba::printDebug("multienrichjam(): ",
            "geneHitIM does not contain ",
            jamba::formatInt(length(missing_columns)),
            " columns present in geneIM, default values will use 1.");
      }
   }
   # if geneIMdirection exists, confirm geneIM only contains c(0,1)
   if (length(geneIMdirection) > 0) {
      geneIM <- (abs(geneIM) > 0) * 1;
   }

   #####################################################################
   ## geneIM colors
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
   # geneIMcolors <- colorjam::matrix2heatColors(x=geneIM,
   #    transformFunc=c,
   #    colorV=colorV,
   #    shareNumLimit=FALSE,
   #    numLimit=1,
   #    trimRamp=c(4,3));

   # assign to mem
   mem$geneHitList <- geneHitList;
   mem$geneHitIM <- geneHitIM;
   mem$geneIM <- geneIM;
   mem$geneIMcolors <- geneIMcolors;
   mem$geneIMdirection <- geneIMdirection;

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
   # jamba::printDebug("enrichLsetNames:");print(enrichLsetNames);# debug
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
      emptyValue=1);
   match1 <- match(enrichLsetNames, rownames(enrichIM));
   match2 <- match(names(enrichList), colnames(enrichIM));
   enrichIM <- enrichIM[match1, match2, drop=FALSE];
   rownames(enrichIM) <- enrichLsetNames;
   colnames(enrichIM) <- names(enrichList);
   # jamba::printDebug("enrichIM:");print(enrichIM);# debug

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
   i1use <- rownames(enrichIMM)[(apply(enrichIMM[,useCols,drop=FALSE], 1, min, na.rm=TRUE) <= cutoffRowMinP)];
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
      verbose=(verbose - 1) > 0);

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
         geneDelim));
   # 0.0.93.900 - enforce consistent order
   memIM <- memIM[, enrichLsetNames, drop=FALSE];
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
   mem$p_cutoff <- cutoffRowMinP;

   if ("Mem" %in% returnType) {
      mem <- list_to_Mem(mem);
   }
   return(mem);
}
