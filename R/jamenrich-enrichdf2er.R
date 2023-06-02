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
#' @param pvalueCutoff `numeric` value range 0 to 1, to define the
#'    P-value threshold for enrichment results to be considered
#'    in downstream processing.
#' @param pAdjustMethod `character` string to define the P-value
#'    adjustment method, or `"none"` for no additional adjustment.
#'    See `stats::p.adjust()` for valid values.
#' @param keyColname `character` value of the `colname(enrichDF)`
#'    containing the unique row identifier. It can be a pathway_ID
#'    or any uniquely identifying value.
#' @param geneColname `character` value of the `colname(enrichDF)`
#'    containing delimiited genes in each pathway.
#' @param pathGenes `character` or value of the `colname(enrichDF)`
#'    containing the number of genes in each pathway. This value will be
#'    derived from `geneRatioColname` if needed.
#' @param geneHits `character` value of the `colname(enrichDF)`
#'    containing the integer count of the gene hits in each pathway.
#'    This value will be derived from `geneRatioColname` if needed.
#' @param geneRatioColname `character` value of the `colname(enrichDF)`
#'    containing the character ratio of gene hits to pathway size,
#'    in format "50/100". This value is used when either `"pathGenes"`
#'    or `"geneHits"` are not supplied.
#' @param geneDelim `character` regular expression pattern used to separate
#'    genes in the `pathGenes` column into a vector of character
#'    values.
#' @param pvalueColname `character` value of the `colname(enrichDF)`
#'    containing enrichment P-values to use in downstream processing.
#' @param msigdbGmtT optional gmtT object, a representation of
#'    `arules::transactions-class`. (Not currently implemented.)
#' @param verbose `logical` indicating whether to print verbose output.
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
 geneSep=",",
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
