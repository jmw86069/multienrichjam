
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
   directionColname=c("activation.z.{0,1}score",
      "z.{0,1}score"),
   direction_cutoff=1,
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
   directionColname <- find_colname(directionColname,
      enrichDF,
      require_non_na=FALSE);
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
      jamba::printDebug("topEnrichBySource(): ",
         "directionColname:", directionColname);
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
      # add internal sort_rank for reference.
      enrichDF$sort_rank <- seq_len(nrow(enrichDF));
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

   ## Apply direction filter
   if (length(directionColname) > 0 &&
         length(direction_cutoff) > 0 &&
         !is.na(direction_cutoff) &&
         direction_cutoff > 0) {
      direction_subset <- (!is.na(enrichDF[[directionColname]]) &
         abs(enrichDF[[directionColname]]) >= direction_cutoff);
      if (verbose) {
         jamba::printDebug("topEnrichBySource(): ",
            c("Applying filter for abs(direction) >= ",
               direction_cutoff,
               " retaining ",
               jamba::formatInt(sum(direction_subset)),
               " out of ",
               jamba::formatInt(nrow(enrichDF)),
               " total entries."),
            sep="");
      }
      enrichDF <- subset(enrichDF, direction_subset);
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
   directionColname=c("activation.z.{0,1}score",
      "z.{0,1}score"),
   direction_cutoff=1,
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
         directionColname=directionColname,
         direction_cutoff=direction_cutoff,
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
