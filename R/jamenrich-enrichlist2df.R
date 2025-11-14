
#' Convert enrichList to data.frame
#'
#' Convert enrichList to a single data.frame, merging results
#'
#' This function takes one or more enrichment results, and
#' combines results into one `data.frame`.
#'
#' * Genes are combined, reporting the unique genes in sorted order.
#' * `numeric` columns report the lowest or highest value, as appropriate:
#'
#'    * P-value columns report the lowest observed value, in `pvalueColname`
#'    * Count columns report the highest observed value, except the
#'    primary count column reports the number of unique genes.
#'
#' @returns `data.frame` with combined results across all `enrichList`.
#'
#' @family jam conversion functions
#' 
#'
#' @param enrichList `list` of `enrichResult` or `data.frame` objects.
#' @param keyColname,geneColname,geneCountColname,pvalueColname
#'    `character` string indicating each column to use, or a vector
#'    of column names in order of priority.
#' @param pvalueFloor `numeric` lowest permitted P-value, intended to
#'    prevent using actual zero `0` which may occur sometimes when
#'    upstream processing rounds the value.
#' @param msigdbGmtT not currently implemented.
#' @param geneDelim `character` pattern used to split genes.
#' @param geneSep `character` string used to join genes, default `","`
#'    uses comma-delimited values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `integer` used for debugging, with values `0, 1, 2, 3`
#'    recognized. Value `0`, default, performs no debug.
#' @param ... additional arguments are ignored.
#'
#' @export
enrichList2df <- function
(enrichList,
 keyColname=c("ID",
    "Name",
    "pathway",
    "itemsetID",
    "Description"),
 geneColname=c("geneID",
    "geneNames",
    "Genes"),
 geneCountColname=c("gene_count",
    "Count",
    "geneCount",
    "geneHits"),
 pvalueColname=c("padjust", "p.adjust", "adjp", "padj",
    "qvalue", "qval", "q.value", "pvalue",
    "p.value", "pval", "FDR"),
 pvalueFloor=1e-200,
 msigdbGmtT=NULL,
 geneDelim="[,/ ]+",
 geneSep=",",
 verbose=FALSE,
 debug=0,
 ...)
{
   ## Purpose is to combine a list of enrichment data.frames into one data.frame
   if (!suppressPackageStartupMessages(require(matrixStats))) {
      stop("enrichList2df() requires the matrixStats package.");
   }
   iDF1 <- head(as.data.frame(enrichList[[1]]), 3);
   # version 0.0.56.900: changed to use enrichList which tests all results not
   # just the first in the list
   keyColname <- find_colname(keyColname, enrichList);
   geneColname <- find_colname(geneColname, enrichList);
   geneCountColname <- find_colname(geneCountColname, enrichList);
   pvalueColname <- find_colname(pvalueColname, enrichList);
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
         # jamba::nameVector(matrixStats::rowMins(enrichIMP), rownames(enrichIMP));
         jamba::nameVector(apply(enrichIMP, 1, min, na.rm=TRUE), rownames(enrichIMP));
      } else if (useType %in% "hi") {
         # jamba::nameVector(matrixStats::rowMaxs(enrichIMP), rownames(enrichIMP));
         jamba::nameVector(apply(enrichIMP, 1, max, na.rm=TRUE), rownames(enrichIMP));
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
         iDF[, c(keyColname, iName), drop=FALSE];
      });
      enrichL1L <- jamba::mergeAllXY(enrichL1L1);
      enrichL1V <- jamba::nameVector(
         gsub("^[,/ ]+|[,/ ]+$",
            "",
            jamba::pasteByRow(
               enrichL1L[, -match(keyColname, colnames(enrichL1L)), drop=FALSE],
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
      # enrichGeneVL <- list(jamba::cPaste(enrichGeneL, doSort=FALSE));
      # 0.0.92.900 - add sort,unique
      enrichGeneVL <- list(jamba::cPasteSU(enrichGeneL,
         sep=geneSep));
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
      if (length(keyValsUse) > 0) {
         keepColDF[keyValsUse,keyColname] <- keyValsUse;
         for (keepCol in keepCols) {
            keyNew <- iDF[match(keyValsUse, iDF[,keyColname]),keepCol];
            if (length(keyNew) > 0) {
               keepColDF[keyValsUse,keepCol] <- keyNew;
            }
         }
      }
      if (verbose) {
         jamba::printDebug("iName:", iName,
            ", length(keyVals):", length(keyVals),
            ", length(keyValsUse):", length(keyValsUse));
         jamba::printDebug("head(iDF):");
         print(head(iDF));
         jamba::printDebug("head(keepColDF):");
         print(head(keepColDF));
      }
      rm(iDF);
   }
   if (debug == 3) {
      return(list(enrichCols=enrichCols,
         enrichValuesM=enrichValuesM,
         enrichList=enrichList,
         enrichGeneLen=enrichGeneLen,
         enrichGeneL=enrichGeneL,
         enrichGeneVL=enrichGeneVL,
         keepColDF=keepColDF));
   }
   if (verbose) {
      jamba::printDebug("multiEnrichMap(): ",
         "Creating enrichDF.");
   }
   enrichDF <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      enrichValuesM,
      keepColDF[match(rownames(enrichValuesM), rownames(keepColDF)),
         , drop=FALSE],
      as.data.frame(enrichGeneVL)[rownames(enrichValuesM), , drop=FALSE]);
   enrichDF[["allGeneHits"]] <- enrichGeneLen;
   #whichCol1 <- max(which(colnames(enrichDF) %in% names(enrichCols))) + 1;
   #enrichDF <- insertDFcols(enrichDF, colnum=whichCol1,
   #   insertDF=data.frame(allGeneHits=enrichGeneLen));
   enrichDF;
}
