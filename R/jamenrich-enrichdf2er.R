#' Convert data.frame to enrichResult
#'
#' Convert data.frame to enrichResult
#'
#' This function takes a `data.frame` containing gene set enrichment
#' results, and converts it to a proper `enrichResult`` object
#' defined in `DOSE` for R-4.5 and older, or `enrichit` for
#' R-4.6 and newer.
#' 
#' This object is supported by other functions in the `clusterProfiler`
#' suite of tools.
#' 
#' @returns `enrichResult` object.
#'
#' @param enrichDF `data.frame` representing gene set enrichment
#'    results.
#' @param pvalueCutoff `numeric` value range 0 to 1, to define the
#'    P-value threshold for enrichment results to be considered
#'    in downstream processing.
#' @param pAdjustMethod `character` string to define the P-value
#'    adjustment method, or `"none"` for no additional adjustment.
#'    See `stats::p.adjust()` for valid values.
#' @param keyColname `character` value of the `colnames(enrichDF)`
#'    containing the unique row identifier. It can be a pathway_ID
#'    or any uniquely identifying value.
#' @param geneColname `character` value of the `colnames(enrichDF)`
#'    containing delimiited genes in each pathway.
#' @param pathGenes `character` or value of the `colnames(enrichDF)`
#'    containing the number of genes in each pathway. This value will be
#'    derived from `geneRatioColname` if needed.
#' @param geneHits `character` value of the `colnames(enrichDF)`
#'    containing the integer count of the gene hits in each pathway.
#'    This value will be derived from `geneRatioColname` if needed.
#' @param geneRatioColname `character` value of the `colnames(enrichDF)`
#'    containing the character ratio of gene hits to pathway size,
#'    in format "50/100". This value is used when either `"pathGenes"`
#'    or `"geneHits"` are not supplied.
#' @param geneDelim `character` regular expression pattern used to separate
#'    genes in the `pathGenes` column into a vector of character
#'    values.
#' @param geneSep `character` string with separator to use for output
#'    genes. Only `'/'` is supported for most `enrichResult` methods.
#' @param pvalueColname,padjustColname,qvalueColname `character` vector
#'    to search `colnames(enrichDF)`, with enrichment P-value,
#'    adjusted P-value, and Q-value, respectively.
#'    * The P-value is stored in 'pvalue' in the output object, to match
#'    `enrichResult` convention.
#'    * When there is no adjusted P-value column matched, by default it
#'    uses the P-value without change, since some downstream methods
#'    require 'p.adjust' exists.
#'    * When there is no Q-value column matched, by default it uses
#'    the adjusted P-value column value 'p.adjust' without change.
#'    Note that the 'p.adjust' may also contain the value from 'pvalue'
#'    as described above.
#' @param descriptionColname `character` vector
#'    to search `colnames(enrichDF)`, for pathway or gene set description.
#' @param readable `logical` default NULL, sets the 'readable' flag for
#'    the resulting `enrichResult` object.
#' @param msigdbGmtT optional 'GmtT' object (not currently implemented)
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam import functions
#' @family jam conversion functions
#'
#' @export
enrichDF2enrichResult <- function(
   enrichDF=NULL,
   pvalueCutoff=1,
   pAdjustMethod="none",
   keyColname=c("itemsetID", "ID", "Name", "Pathway"),
   pathGenes="pathGenes",
   geneColname=c("geneNames", "geneID", "Gene", "Genes"),
   geneHits="geneHits",
   geneRatioColname=c("GeneRatio", "^Ratio"),
   geneDelim="[,/ ]+",
   geneSep="/",
   pvalueColname=c("P.Value", "Pvalue", "Pval", "FDR", "adj.P.Val"),
   padjustColname=c("p.adjust", "p.adjusted", "padjust", "adjp", "adj.P.Val", "FDR"),
   qvalueColname=c("Q.Value", "Qvalue", "qval", "FDR"),
   descriptionColname=c("Description", "Name", "Pathway", "ID"),
   readable=NULL,
   msigdbGmtT=NULL,
   verbose=FALSE,
   ...
) {
   ## Purpose is to convert an enrichment data.frame
   ## into enrichResult class format usable by clusterProfiler
   ## methods, like enrichMap()

   ## Find each colname in the input data.frame
   keyColname <- find_colname(keyColname, enrichDF);
   pathGenes <- find_colname(pathGenes, enrichDF);
   geneColname <- find_colname(geneColname, enrichDF);
   geneHits <- find_colname(geneHits, enrichDF);
   geneRatioColname <- find_colname(geneRatioColname, enrichDF);
   pvalueColname <- find_colname(pvalueColname, enrichDF);
   padjustColname <- find_colname(padjustColname, enrichDF);
   qvalueColname <- find_colname(qvalueColname, enrichDF);
   descriptionColname <- find_colname(descriptionColname, enrichDF);
   if (verbose) {
      cli::cli_inform(c(
         "{.pkg enrichDF2enrichResult} colnames matched: ",
         "keyColname: {.val {keyColname}}",
         "pathGenes: {.val {pathGenes}}",
         "geneColname: {.val {geneColname}}",
         "geneHits: {.val {geneHits}}",
         "geneRatioColname: {.val {geneRatioColname}}",
         "pvalueColname: {.val {pvalueColname}}",
         "padjustColname: {.val {padjustColname}}",
         "qvalueColname: {.val {qvalueColname}}",
         "descriptionColname: {.val {descriptionColname}}"
      ));
   }
   if (length(c(keyColname, pvalueColname, geneColname)) < 3) {
      cli::cli_abort(paste0(
         "Could not find {.var keyColname}, {.var pvalueColname}, ",
         "or {.var geneColname} in the input colnames."));
   }

   ## Confirm Description column
   if (length(descriptionColname) == 0) {
      cli::cli_warn(c(
         "No {.var descriptionColname} colname was found in {.var enrichDF},",
         "which prevents this data from being used by",
         "{.pkg enrichplot}, {.pkg clusterProfiler} functions.",
         "The Description column is recommended to contain",
         "the full name of the gene set."));
   } else if (descriptionColname %in% c(keyColname)) {
      # If it is also the keyColname we need to make two columns
      # so both names can co-exist.
      if (!"Description" %in% descriptionColname) {
         enrichDF$Description <- enrichDF[[descriptionColname]];
      }
   } else {
      ## Otherwise rename to make sure the final colname
      ## is 'Description' to fit expectations of
      ## enrichplot:::fortify.internal()
      enrichDF <- jamba::renameColumn(enrichDF,
         from=descriptionColname,
         to="Description");
   }

   enrichDF2 <- jamba::renameColumn(enrichDF,
      from=c(keyColname, pvalueColname, geneColname),
      to=c("ID", "pvalue", "geneID"));
   # p.adjust
   if (length(padjustColname) == 1) {
      enrichDF2[["p.adjust"]] <- enrichDF2[[padjustColname]];
   } else {
      enrichDF2[["p.adjust"]] <- enrichDF2[["pvalue"]];
   }
   # qvalue
   if (length(qvalueColname) == 1) {
      enrichDF2[["qvalue"]] <- enrichDF2[[qvalueColname]];
   } else {
      enrichDF2[["qvalue"]] <- enrichDF2[["p.adjust"]];
   }


   ## Ensure all entries in column "ID" are unique
   ## because these values need to become rownames for compatibility
   ## with enrichplot::cnetplot()
   enrichDF2[["ID"]] <- jamba::makeNames(enrichDF2[["ID"]]);
   rownames(enrichDF2) <- enrichDF2[["ID"]];

   ## Convert gene delimiters all to "/"
   enrichDF2[["geneID"]] <- gsub(geneDelim,
      geneSep,
      enrichDF2[["geneID"]]);
   
   ## set readable when needed
   if (length(readable) == 0) {
      readable <- FALSE;
      # if any alphabetic or -_ we consider it "readable" and not ENTREZID
      if (jamba::igrepHas("[-_a-zA-Z]", enrichDF2[["geneID"]])) {
         readable <- TRUE;
      }
   } else {
      readable <- head(as.logical(readable), 1);
   }

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
      unique(jamba::provigrep(c("^ID$", "."), colnames(enrichDF2))),
      colnames(enrichDF2));
   enrichDF2 <- enrichDF2[, keepcolids, drop=FALSE];
   #enrichDF2a <- dplyr::select(enrichDF2,
   #   dplyr::matches("^ID$"), tidyselect::everything());
   #enrichDF2 <- enrichDF2a;
   if (verbose) {
      jamba::printDebug("enrichList2df(): ",
         "Done.");
   }

   gene <- jamba::mixedSort(unique(unlist(
      strsplit(
         as.character(enrichDF2[["geneID"]]),
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
   return(
      methods::new(
         "enrichResult",
         result=enrichDF2,
         pvalueCutoff=pvalueCutoff,
         pAdjustMethod=pAdjustMethod,
         gene=as.character(gene),
         universe=universe,
         geneSets=geneSets,
         organism="UNKNOWN",
         keytype="UNKNOWN",
         ontology="UNKNOWN",
         readable=readable
      ))
}
