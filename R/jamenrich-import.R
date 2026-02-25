
#' Import Ingenuity Pathway Analysis 'IPA' results
#'
#' Import Ingenuity Pathway Analysis 'IPA' results, by default
#' reverting IPA symbols to input values
#'
#' This function parses Ingenuity Pathway Analysis ('IPA')
#' enrichment data into a `list` of `data.frame`
#' objects for downstream analysis.
#' Each `data.frame` represents the results of one Ingenuity IPA test,
#' however not all sections contain gene set enrichment results.
#' 
#' ## Batch processing
#' 
#' When importing multiple files, argument `ipaFile` can be a vector
#' of '.xlsx' or '.txt' files. This workflow also calls
#' `IPAlist_to_hits()` to generate a gene hit matrix, stored
#' as `attr(ipalist, "geneHitIM")` to use in `multiEnrichMap()`.
#' 
#' ## IPA Gene Xref Data
#' 
#' By default, the argument `revert_ipa_xref=TRUE` will convert the
#' IPA gene symbol values back to the original identifier.
#' In most cases* this behavior is desirable, with caveats:
#'
#' * When using platform data which
#' use a non-gene identifier, including microarray probesets, or
#' RefSeq transcript ID, or protein "UniProt" accession numbers,
#' it is recommended to use `revert_ipa_xref=FALSE`.
#' In these cases, the IPA gene symbol is expected to be
#' more user-friendly, and therefore more useful.
#'
#'    * This option would be helpful to view the IPA gene symbols
#'    as they appear in the IPA report - even if the symbols
#'    sometimes do not match the input row identifiers.
#'
#' * It is helpful to use `revert_ipa_xref=TRUE` when the identifier
#' will also be used to compare to the source data, for example
#' if trying to make an expression heatmap of the genes involved
#' in enrichment results.
#'
#'    * This option would convert IPA gene symbol back to Affymetrix
#'    probeset ID, for example, if the probeset ID values were
#'    used as the primary identifier for each measurement.
#'    The probeset ID might be convenient to align with the input
#'    data matrix.
#'    * The primary reason for this option is when providing
#'    gene symbols as input to IPA, some will be renamed to IPA
#'    preferred gene symbols, which would therefore be difficult
#'    to match with the gene symbols provided to IPA.
#'
#' In any case, the output `list` should contain an entry
#' "Analysis Ready Molecules" with the full IPA data table used
#' for the analysis. This `data.frame` will also contain any
#' statistical columns, if provided to IPA upfront.
#' 
#' See `IPAlist_to_hits()` for an automated way to create a gene hit
#' matrix from the 'Analysis Ready Molecules' returned by IPA.
#'
#' ## Motivation
#'
#' 1. Separate multiple IPA enrichment tables.
#' 2. Rename colnames to be consistent, compatible with
#' `enrichDF2enrichResult()`.
#' 3. Revert IPA gene aliases to original user input, default but optional.
#' 4. Generate `geneHitIM` hit matrix when processing multiple files.
#'
#' ## Input format
#'
#' 1. `ipaFile` can be one or more text `.txt` files,
#' where the text file contains all IPA enrichment data in
#' tall format. This format is most common.
#' 2. `ipaFile` can be one or more Excel `.xlsx` files,
#' which contains all IPA enrichment data in
#' one tall worksheet tab.
#' 3. `ipaFile` can be one or more Excel `.xlsx` files,
#' where each type of IPA enrichment appears on a separate
#' Excel worksheet tab.
#' 4. `ipaFile` can be a list of `data.frame` objects.
#' This option is intended when the IPA data has already
#' been imported into R as separate `data.frame` objects.
#'
#' ## Notes
#'
#' When using `"Export All"` from 'IPA', the default text format
#' includes multiple enrichment tables concatenated together in one
#' file. Each enrichment table contains its own unique column
#' headers, with descriptive text in the line preceding the
#' column headers. This function is intended to separate the
#' enrichment tables into a list of `data.frame` objects, and
#' retain the descriptive text as names of the list.
#'
#' ## Troubleshooting
#'
#' * A common error occurs when reverting IPA gene symbols
#' to the original user-supplied identifier, by default
#' `revert_ipa_xref=TRUE`.
#' For errors during this step, consider `revert_ipa_xref=FALSE`
#' which will retain the gene symbol as recognized by IPA.
#' The downside of this approach is that it may be more difficult
#' to equate to the input identifier.
#' In that case look at the "Analysis Ready Molecules" `data.frame`
#' which should contain the
#' user-provided values as "ID"; the IPA recognized symbol as "Name",
#' and optionally a column "Symbol" which is edited by multienrichjam.
#'
#' @returns `list` of `data.frame` objects, where each `data.frame`
#'    contains enrichment data for one of the Ingenuity IPA
#'    enrichment tests.
#'
#' @family jam import functions
#'
#' @param ipaFile one of the four input types described above:
#'    a character vector of text file names; a character vector of
#'    Excel `.xlsx` file names; a list of `data.frame` objects.
#' @param headerGrep regular expression pattern used to
#'    recognize header columns found in Ingenuity IPA enrichment data.
#' @param ipaNameGrep vector of regular expression patterns
#'    used to recognize the name of the enriched entity,
#'    for example the biological pathway, or network, or
#'    disease category, etc.
#' @param geneGrep regular expression pattern used to recognize
#'    the column containing genes, or the molecules tested for
#'    enrichment which were found in the enriched entity.
#' @param geneCurateFrom,geneCurateTo vector of patterns and
#'    replacements, respectively, used to curate values in
#'    the gene column. These replacement rules are used to
#'    ensure that genes are delimited consistently, with no
#'    leading or trailing delimiters.
#' @param method integer value indicating the method used to
#'    import data from a text file, where: `method=1` uses
#'    `data.table::read.table()` and the `textConnection` argument;
#'    `method=2` uses `readr::read_tsv()`. The motivation to use
#'    `data.table::read.table()` is it performed better in the
#'    presence of UTF-8 characters such as the alpha symbol.
#' @param sheet integer value used only when `ipaFile` is
#'    a vector of Excel `.xlsx` files, and when the Excel
#'    format includes multiple worksheets. This value will
#'    extract enrichment data only from one worksheet from
#'    each Excel file.
#' @param sep character string used when `ipaFile` is a vector
#'    of text files, to split fields into columns. The default
#'    will split fields by the tab character.
#' @param xlsxMultiSheet logical indicating whether input
#'    Excel `.xlsx` files contain multiple worksheets.
#' @param useXlsxSheetNames logicl indicating whether to use the
#'    Excel worksheet name for each imported enrichment table,
#'    when importing from `.xlsx` files, and when
#'    `xlsxMultiSheet=FALSE`. When `xlsxMultiSheet=TRUE` the
#'    name is derived from the value matched using `ipaNameGrep`,
#'    because in this case, there are expected to me multiple
#'    enrichment tables in one worksheet.
#' @param remove_blank_colnames `logical` indicating whether to drop
#'    `colnames()` where all values are contained in `c(NA, "")`.
#'    This option may be preferable `remove_blank_colnames=FALSE`
#'    when all values in some column like `zScore` are `NA`, but
#'    you would still like to retain the column for consistency
#'    with other data. We found that IPA does not report `zScore`
#'    values when there are only 4 or fewer genes involved in
#'    each enrichment result.
#' @param convert_ipa_slash `logical` indicating whether to convert
#'    IPA gene naming conventions, currently some genes are considered
#'    one entity in the IPA system, for example `"HSPA1A/HSPA1B"` is
#'    considered one gene, even though two Entrez gene entries
#'    `"HSPA1A"` and `"HSPA1B"` can be represented. Regardless whether
#'    one or both genes are provided to IPA, it considers it one
#'    entity for the purpose of pathway enrichment hypergeometric testing.
#'    Unfortunately, the forward slash `"/"` is also used by
#'    `clusterProfiler` object `enrichResult` as gene delimiter, and
#'    is hard-coded and cannot be changed. So it will automatically
#'    consider `"HSPA1A/HSPA1B"` as two genes, causing a mismatch with
#'    the IPA results.
#'    When `convert_ipa_slash=TRUE` by default, it converts the
#'    forward slash `"/"` to the value of argument `ipa_slash_sep`.
#' @param ipa_slash_sep `character` string used as a delimited when
#'    `convert_ipa_slash=TRUE`, used to replace genes that contain
#'    forward slash `"/"` to use another character.
#' @param revert_ipa_xref `logical` indicating whether to revert the
#'    IPA gene symbols reported, which requires that the IPA data
#'    contains a section `"Analysis Ready Molecules"`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' ipaFile <- system.file(package="multienrichjam", "extdata",
#'    c("Newborns-IPA.txt", "OlderChildren-IPA.txt"));
#' ipalist <- importIPAenrichment(ipaFile)
#' 
#' @export
importIPAenrichment <- function
(ipaFile,
 headerGrep="(^|\t)((expr.|-log.|)p-value|Pvalue|Score($|\t)|Symbol($|\t)|Ratio($|\t)|Consistency.Score|Master.Regulator($|\t))",
 ipaNameGrep=c("Pathway",
    "Regulator$",
    "Regulators",
    "Regulator", #"Master.Regulator",
    "Disease",
    "Toxicity",
    "Category",
    "Categories",
    "Function",
    "Symbol$",
    "^ID$",
    "My.(Lists|Pathways)"),
 geneGrep=c("Molecules in Network",
    "Target molecules",
    "Molecules",
    "Symbol"),
 geneCurateFrom=c("[ ]*[(](complex|includes others)[)][ ]*",
    "^[, ]+|[, ]+$"),
 geneCurateTo=c("",
     ""),
 signColname=c("Expr Fold Change", "fold.*change",
    "log.*ratio", "log.*fold", "log.*fc",
    "lfc", "ratio", "fold", "fc"),
 signThreshold=0,
 method=1,
 sheet=1,
 sep="\t",
 xlsxMultiSheet=TRUE,
 useXlsxSheetNames=FALSE,
 remove_blank_colnames=TRUE,
 convert_ipa_slash=TRUE,
 ipa_slash_sep=":",
 revert_ipa_xref=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to encapsulate several steps required to import the
   ## IPA Ingenuity Pathway Analysis file created with "Export All."

   ## Import the data, only taking rows containing tab-delimited
   ## columns.
   ## Note: this step preserved UTF-8 characters, since IPA often
   ## encodes symbols like alpha and beta as the proper UTF-8
   ## character.
   ##
   ## xlsxMultiSheet=TRUE is intended for Excel import where the Excel file
   ## contains separate worksheets, each with its own properly defined
   ## summary table.
   ##
   ## xlsxMultiSheet=FALSE is intended for Excel import of IPA data where
   ## the Excel table contains all summary tables appended one after another.
   if (jamba::igrepHas("[.]xlsx$", ipaFile)) {
      if (!requireNamespace("openxlsx", quietly=TRUE)) {
         stop_msg <- paste0("importIPAenrichment() requires the openxlsx ",
            "package for Excel import.");
         stop(stop_msg);
      }
   }

   if (jamba::igrepHas("list", class(ipaFile)) &&
         inherits(ipaFile[[1]], c("data.frame", "matrix"))) {
      ## Process list of data.frames as input
      ipaDFL <- lapply(jamba::nameVectorN(ipaFile), function(iSheet){
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "from data.frame list name:", iSheet);
         }
         iDF <- ipaFile[[iSheet]];
         jDF <- curateIPAcolnames(iDF,
            ipaNameGrep=ipaNameGrep,
            geneGrep=geneGrep,
            geneCurateFrom=geneCurateFrom,
            geneCurateTo=geneCurateTo,
            convert_ipa_slash=convert_ipa_slash,
            ipa_slash_sep=ipa_slash_sep,
            verbose=verbose);
         jDF;
      });
   } else if (jamba::igrepHas("[.](txt|xlsx)$", ipaFile) && length(ipaFile) > 1) {
      ## Process multiple files by calling this function on each file
      ipaDFLL <- lapply(ipaFile, function(ipaFile_i){
         jDF <- importIPAenrichment(ipaFile=ipaFile_i,
            headerGrep=headerGrep,
            ipaNameGrep=ipaNameGrep,
            geneGrep=geneGrep,
            geneCurateFrom=geneCurateFrom,
            geneCurateTo=geneCurateTo,
         	signColname=signColname,
         	signThreshold=signThreshold,
            method=method,
            sheet=sheet,
            sep=sep,
            xlsxMultiSheet=xlsxMultiSheet,
         	useXlsxSheetNames=useXlsxSheetNames,
         	remove_blank_colnames=remove_blank_colnames,
            convert_ipa_slash=convert_ipa_slash,
            ipa_slash_sep=ipa_slash_sep,
            revert_ipa_xref=revert_ipa_xref,
            verbose=verbose,
            ...);
      });
      # extract IPA gene hit matrix?
      tryCatch({
         ipahitim <- IPAlist_to_hits(ipaDFLL,
            signColname=signColname,
            signThreshold=signThreshold,
            returnType="matrix",
            ...);
         attr(ipaDFLL, "geneHitIM") <- ipahitim;
      }, error=function(e){
         # do nothing for now
      });
      return(ipaDFLL);
   } else if (jamba::igrepHas("[.]xlsx$", ipaFile) && xlsxMultiSheet) {
      ## Import IPA data from Excel xlsx file, multiple worksheets
      ## Process Excel import instead of CSV
      sheetNames <- openxlsx::getSheetNames(ipaFile);
      sheetNamesUse <- jamba::nameVector(sheet, sheetNames[sheet]);
      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "importing xlsx as from multiple worksheets:",
            ipaFile);
         jamba::printDebug("importIPAenrichment(): ",
            "sheetNamesUse:");
         print(sheetNamesUse);
      }
      ipaDFL <- lapply(sheetNamesUse, function(iSheet){
         if (verbose) {
            jamba::printDebug("   importIPAenrichment(): ",
               "iSheet:", iSheet);
         }
         if (1 == 2) {
            iDF <- read.xlsx(xlsxFile=ipaFile,
               sheet=iSheet);
            jDF <- curateIPAcolnames(iDF,
               ipaNameGrep=ipaNameGrep,
               geneGrep=geneGrep,
               geneCurateFrom=geneCurateFrom,
               geneCurateTo=geneCurateTo,
               convert_ipa_slash=convert_ipa_slash,
               ipa_slash_sep=ipa_slash_sep,
               verbose=verbose);
            jDF;
         } else {
            jDF <- importIPAenrichment(ipaFile=ipaFile,
               headerGrep=headerGrep,
               ipaNameGrep=ipaNameGrep,
               geneGrep=geneGrep,
               geneCurateFrom=geneCurateFrom,
               geneCurateTo=geneCurateTo,
            	signColname=signColname,
            	signThreshold=signThreshold,
               method=method,
               sheet=sheet,
               sep=sep,
               xlsxMultiSheet=FALSE,
            	useXlsxSheetNames=useXlsxSheetNames,
            	remove_blank_colnames=remove_blank_colnames,
               convert_ipa_slash=convert_ipa_slash,
               ipa_slash_sep=ipa_slash_sep,
            	revert_ipa_xref=revert_ipa_xref,
               verbose=verbose,
               ...);
         }
      });
      ipaDFL <- unlist(recursive=FALSE, unname(ipaDFL));
   } else {
      ## Import IPA data from Excel xlsx or txt file, single worksheet
      if (jamba::igrepHas("[.]xlsx$", ipaFile) && !xlsxMultiSheet) {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "importing xlsx as a single worksheet:",
               ipaFile);
         }
         iDF <- openxlsx::read.xlsx(ipaFile,
            sheet=1,
            colNames=FALSE);
         i <- gsub("\t+$",
            "",
            jamba::pasteByRow(iDF,
               sep="\t",
               condenseBlanks=FALSE));
         #i <- jamba::vigrep("\t.*\t", i);
      } else {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "importing text using grep from ipaFile:",
               ipaFile);
         }
         #i <- system(intern=TRUE,
         #      paste0("grep '\t.*\t' ", ipaFile));
         i <- jamba::vigrep(".", readLines(ipaFile));
         #i <- system(intern=TRUE,
         #   paste0("grep '.*' ", ipaFile));
         ## Clean up trailing newlines
         i <- gsub("[\r\n]+$",
            "",
            i);
         #i;
      }
      ## Often IPA output has descriptor lines with 'for ... ->' in them
      ## which indicates a header for each subtable. We will parse
      ## that header for their correct table name where possible.
      is_desc <- grep("( for ).*(->)", i);

      ## Infer which columns area headers, then set up
      ## ranges of rows to be assigned to each header.
      #i1 <- jamba::igrep(headerGrep, i);
      is_header <- jamba::igrep(headerGrep, i);
      is_footer <- c(tail(is_header, -1) - 2, length(i));

      ## Note that a "valid" condition must have
      ## - descLines one row above i1
      ## - i1 must be contained in (descLines + 1)
      ## - descLines must not be contained in (descLines + 2)
      valid_desc <- ((is_desc + 1) %in% is_header & !(is_desc + 2) %in% is_desc);
      i_desc <- gsub("( for |->).*$", "", i[is_desc]);
      valid_df <- data.frame(is_desc,
         is_header=is_header[match(is_desc, is_header - 1)],
         is_footer=is_footer[match(is_desc, is_header - 1)],
         i_desc,
         valid_desc);
      rownames(valid_df) <- jamba::makeNames(i_desc);
      use_df <- valid_df[valid_desc,,drop=FALSE];
      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "summary of IPA row validation:");
         print(valid_df);
         jamba::printDebug("importIPAenrichment(): ",
            "Rows used during import:");
         print(use_df);
      }


      ## Iterate each set of results and create a data.frame.
      ## Note that each type of result has its own number of columns,
      ## and unique colnames.
      ipaDFL <- lapply(jamba::nameVector(rownames(use_df)), function(j){
         j1 <- use_df[j,"is_header"];
         j2 <- use_df[j,"is_footer"];
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "creating data.frame:",
               j,
               ", j1:", j1, ", j2:", j2);
         }
         ## Create a sequence that removes the descLines as relevant
         jseq <- seq(from=j1, to=j2, by=1);
         if (any(jseq %in% is_desc)) {
            j2 <- min(intersect(jseq, is_desc)) - 1;
            jseq <- seq(from=j1, to=j2, by=1);
         }
         ## Note it is very important to use textConnection() here
         ## because it properly preserves character encodings, such
         ## as alpha, beta characters in some pathway names from IPA.
         if (method == 1) {
            if (verbose) {
               jamba::printDebug("importIPAenrichment(): ",
                  "using read.table().");
            }
            jDF <- read.table(
               textConnection(gsub(paste0(sep, "+$"),
                  "",
                  i[jseq])),
               sep=sep,
               stringsAsFactors=FALSE,
               quote="\"",
               check.names=FALSE,
               as.is=TRUE,
               header=TRUE,
               encoding="UTF-8",
               comment.char="",
               fill=TRUE);
         } else {
            if (!requireNamespace("readr", quietly=TRUE)) {
               stop("importIPAenrichment() requires the readr package.");
            }
            if (verbose) {
               jamba::printDebug("importIPAenrichment(): ",
                  "using readr::read_tsv().");
            }
            jDF <- readr::read_tsv(
               #paste0(gsub("\t$", "", i[j1:j2]),
               paste0(gsub(paste0(sep, "+$"),
                  "",
                  i[jseq]),
                  collapse="\n"),
               quote="\"",
               col_names=TRUE,
               #encoding="UTF-8",
               comment="");
         }
         ## Remove empty colnames
         if (TRUE %in% remove_blank_colnames) {
            non_na_cols <- sapply(colnames(jDF), function(j){
               !all(jDF[[j]] %in% c(NA, ""))
            });
            jDF <- jDF[, non_na_cols, drop=FALSE];
         }

         ## Curate IPA colnames
         jDF <- curateIPAcolnames(jDF,
            ipaNameGrep=ipaNameGrep,
            geneGrep=geneGrep,
            geneCurateFrom=geneCurateFrom,
            geneCurateTo=geneCurateTo,
            convert_ipa_slash=convert_ipa_slash,
            ipa_slash_sep=ipa_slash_sep,
            verbose=verbose);
         jDF;
      });
   }

   ## Optionally revert IPA gene symbols where possible
   if (TRUE %in% revert_ipa_xref) {
      arm <- head(jamba::vigrep("Analysis.*Ready.*Molecules", names(ipaDFL)), 1)
      if (length(arm) == 0) {
         jamba::printDebug("IPA xref data was not found, expecting section: '",
            "Analysis Ready Molecules",
            "', handling as: ", "revert_ipa_xref=FALSE",
            fgText=c("darkorange", "firebrick"))
      } else {
         ## Now curate gene symbols
         # - ID: contains the input gene
         # - Symbol: contains the gene as updated by this function:
         #    - replaced "/" with ":"
         #    - removed things like "(and others)", and "(complex)"
         # - Name: contains the original IPA label
         ipaxref <- ipaDFL[[arm]]
         if (nrow(ipaxref) > 0) {
            to_updates <- jamba::nameVector(setdiff(names(ipaDFL), arm));
            ipaDFL2 <- lapply(to_updates, function(to_update){
               ipadf <- ipaDFL[[to_update]]
               use_genecol <- jamba::provigrep(
                  unique(c(
                     "geneNames",
                     "participating regulators")),
                  colnames(ipadf))
               # remove colnames with ".ipa" suffix
               use_genecol <- jamba::unvigrep("[.]ipa$", use_genecol);
               if (verbose) {
                  jamba::printDebug("Updating gene xrefs for to_update:",
                     to_update);
                  jamba::printDebug("Updating gene xrefs for use_genecol:",
                     use_genecol);
                  jamba::printDebug("head(ipadf):");# debug
                  print(head(ipadf, 5));# debug
               }
               if (length(use_genecol) == 0) {
                  if (verbose) {
                     jamba::printDebug("No changes made to gene xrefs.");
                  }
                  return(ipaDFL[[to_update]]);
               }
               ## version 0.0.97.900 - iterate multiple values for use_genecol
               #
               # Note this step may not be complete, it only converts 1-to-1.
               # Some IPA entries such as "MAPK (and others)" appear to be
               # used by IPA as a single entity to represent one or more
               # actual measured genes/probes.
               # For example, HSPA1A and HSPA1B are provided at input,
               # they are combined into one row "HSPA1A/HSPA1B",
               # however only "HSPA1B" is stored as the exemplar,
               # apparently because it is the most significant of multiple
               # entries.

               # iterate each gene column and revert gene values if possible
               for (i_genecol in use_genecol) {
                  if (verbose && length(use_genecol) > 1) {
                     jamba::printDebug("Updating i_genecol:",
                        i_genecol);
                  }
                  # split by delimiter
                  ipadf_gl <- strsplit(ipadf[[i_genecol]], ",")
                  ipadf_gll <- lengths(ipadf_gl);

                  #
                  ipadf_gldf <- data.frame(
                     Num=factor(
                        rep(seq_len(nrow(ipadf)), lengths(ipadf_gl)),
                        levels=seq_len(nrow(ipadf))),
                     IPA_Symbol=unlist(ipadf_gl))
                  #
                  match_ipa <- match(ipadf_gldf$IPA_Symbol,
                     ipaxref$Symbol)
                  if (any(!is.na(match_ipa))) {
                     ipadf_gldf$User_Symbol <- ifelse(!is.na(match_ipa),
                        ipaxref$ID[match_ipa],
                        ipadf_gldf$IPA_Symbol)
                     new_geneNames <- jamba::cPasteSU(
                        split(ipadf_gldf$User_Symbol, ipadf_gldf$Num))
                     # optional: make backup column with original values
                     ipadf[[paste0(i_genecol, ".ipa")]] <- ipadf[[i_genecol]];
                     # update values
                     ipadf[[i_genecol]] <- new_geneNames;
                  }
               }
               if (verbose) {
                  jamba::printDebug("head(ipadf) after changes:");# debug
                  print(head(ipadf, 5));# debug
               }
               ipadf
            })
            ipaDFL[to_updates] <- ipaDFL2[to_updates];
            if (verbose) {
               jamba::printDebug("importIPAenrichment(): ",
                  "Updated gene symbols to user input values.");
            }
         }
      }
   }
   ipaDFL;
}

#' Curate Ingenuity IPA colnames
#'
#' Curate Ingenuity IPA colnames
#'
#' This function is intended to help curate colnames observed
#' in Ingenuity IPA enrichment data. The IPA enrichment data
#' includes multiple types of enrichment tests, each with slightly
#' different column headers. This function is intended to
#' make the colnames more consistent.
#'
#' This function will rename the first recognized gene colname
#' to `"geneNames"` for consistency with downstream analyses.
#'
#' The values in the recognized gene colname are curated
#' using `geneCurateFrom,geneCurateTo` for multiple
#' pattern-replacement substitutions. This mechanism is
#' used to ensure consistent delimiters and values used
#' for each enrichment table.
#'
#' Any colname matching `"-log.*p.value"` is considered
#' -log10 P-value, and is converted to normal P-values
#' for consistency with downstream analyses.
#'
#' Any recognized P-value column is renamed to `"P-value"`
#' for consistency with downstream analyses.
#'
#' When the recognized P-value column contains a range,
#' for example `"0.00017-0.0023"`, the lower P-value is
#' chosen. In that case, the higher P-value is stored in
#' a new column `"max P-value"`. P-value ranges
#' are reported in the disease category analysis by
#' Ingenuity IPA, after collating individual pathways
#' by disease category and storing the range of enrichment
#' P-values.
#'
#' @family jam utility functions
#'
#' @param jDF data.frame from one Ingenuity IPA enrichment test.
#' @param ipaNameGrep vector of regular expression patterns
#'    used to recognize the name of the enriched entity,
#'    for example the biological pathway, or network, or
#'    disease category, etc.
#' @param geneGrep regular expression pattern used to recognize
#'    the column containing genes, or the molecules tested for
#'    enrichment which were found in the enriched entity.
#' @param geneCurateFrom,geneCurateTo vector of patterns and
#'    replacements, respectively, used to curate values in
#'    the gene column. These replacement rules are used to
#'    ensure that genes are delimited consistently, with no
#'    leading or trailing delimiters.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
curateIPAcolnames <- function
(jDF,
 ipaNameGrep=c("^Name$",
    "^ID$",
    "Canonical Pathways",
    "Upstream Regulator",
    "Diseases or Functions Annotation",
    "Diseases . Functions",
    "My Lists",
    "Ingenuity Toxicity Lists",
    "My Pathways"),
 geneGrep=c("Molecules in Network",
    "Target molecules",
    "Molecules",
    "Symbol"),
 geneCurateFrom=c(" [(](complex|includes others)[)]",
    "^[,]+|[,]+$"),
 geneCurateTo=c("",
    ""),
 convert_ipa_slash=TRUE,
 ipa_slash_sep=":",
 verbose=TRUE,
 ...)
{
   ## Purpose is to curate some colnames in IPA enrichment data
   ## Determine a column to become the "Name"
   nameCol <- head(jamba::provigrep(ipaNameGrep, colnames(jDF)), 1);
   if (length(nameCol) == 1) {
      jDF[["Name"]] <- jDF[[nameCol]];
      ## Skip rename, in order to maintain the original colname
      #jDF <- jamba::renameColumn(jDF,
      #   from=nameCol,
      #   to="Name");
      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "created Name column from:",
            nameCol);
      }
   }

   ## Determine which column contains the genes of interest
   geneCol <- head(jamba::provigrep(geneGrep, colnames(jDF)), 1);
   if (length(geneCol) == 1) {
      if ("Symbol" == geneCol) {
         geneNamesColumn <- "Symbol";
      } else {
         jDF <- jamba::renameColumn(jDF,
            from=geneCol,
            to="geneNames");
         geneNamesColumn <- "geneNames";
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "created ", "'geneNames'", " column from:",
               geneCol);
         }
      }
      ## Curate gene column values
      jDF[[geneNamesColumn]] <- jamba::gsubs(geneCurateFrom,
         geneCurateTo,
         jDF[[geneNamesColumn]]);

      ## Optionally convert IPA delimiter forward slash "/"
      if (TRUE %in% convert_ipa_slash) {
         if (length(ipa_slash_sep) == 0) {
            ipa_slash_sep <- ""
         }
         if (any(grepl("/", fixed=TRUE, jDF[[geneNamesColumn]]))) {
            if (verbose) {
               jamba::printDebug("importIPAenrichment(): ",
                  "converting IPA multi-gene delimiter from '", "/",
                  "' to: '", ipa_slash_sep, "'.")
            }
            jDF[[geneNamesColumn]] <- gsub("/",
               ipa_slash_sep,
               fixed=TRUE,
               jDF[[geneNamesColumn]])
         }
      }

   }

   ## Convert -log(p-value) columns to P-value for compatibility
   ## with the enrichResult object class.
   logpCol <- head(jamba::vigrep("-log.*p.value", colnames(jDF)), 1);
   if (length(logpCol) == 1) {
      if ("p-value" %in% tolower(colnames(jDF))) {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "Used existing P-value column, left existing column as-is:",
               logpCol);
            #print(head(jDF[,c("P-value",logpCol),drop=FALSE]));
         }
      } else {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "Created P-value column from:",
               logpCol);
            #print(head(jDF[,logpCol,drop=FALSE]));
         }
         jDF[["P-value"]] <- 10^(-as.numeric(jDF[[logpCol]]));
      }
   }
   pCol <- head(setdiff(jamba::vigrep("^p.value", colnames(jDF)), "P-value"), 1);
   if (length(pCol) == 1) {
      if (length(logpCol) != 1) {
         jDF <- jamba::renameColumn(jDF,
            from=pCol,
            to="P-value");
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "Renamed P-value column from:",
               pCol);
         }
      } else {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "Did not rename P-value column from:",
               pCol);
         }
      }
      ## Check for multiple P-values
      if (jamba::igrepHas("-[0-9]+-", head(jDF[["P-value"]], 20))) {
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "Splitting P-value range");
         }
         pvals <- gsub("([0-9]+)-([0-9]+)",
            "\\1!\\2",
            jDF[["P-value"]]);
         pvalsM1 <- jamba::rbindList(strsplit(pvals, "!"));
         pvalsM <- matrix(as.numeric(pvalsM1),
            ncol=ncol(pvalsM1));
         maxPcolnames <- jamba::makeNames(rep("max P-value", ncol(pvalsM)-1));
         colnames(pvalsM) <- c("P-value", maxPcolnames);
         jDF[,colnames(pvalsM)] <- pvalsM;
      }
   }
   return(jDF);
}

#' Pattern replacement with multiple patterns
#'
#' Pattern replacement with multiple patterns
#'
#' This function is a simple wrapper around `base::gsub()`
#' when considering a series of pattern-replacement
#' combinations. It applies each pattern match and replacement
#' in order and is therefore not vectorized.
#'
#' @family jam utility functions
#'
#' @param pattern character vector of patterns
#' @param replacement character vector of replacements
#' @param x character vector with input data to be curated
#' @param ignore.case logical indicating whether to perform
#'    pattern matching in case-insensitive manner, where
#'    `ignore.case=TRUE` will ignore the uppercase/lowercase
#'    distinction.
#' @param replace_multiple logical vector indicating whether to perform
#'    global substitution, where `replace_multiple=FALSE` will
#'    only replace the first occurrence of the pattern, using
#'    `base::sub()`. Note that this vector can refer to individual
#'    entries in `pattern`.
#' @param ... additional arguments are passed to `base::gsub()`
#'    or `base::sub()`.
#'
#' @export
gsubs_remove <- function
(pattern,
 replacement,
 x,
 ignore.case=TRUE,
 replaceMultiple=rep(TRUE, length(pattern)),
 ...)
{
   ## Purpose is to curate a text field using a series of gsub()
   ## commands, operating on a vector of from,to vectors.
   ## 'pattern' is expected to be a vector of regular expression patterns
   ## used by gsub()
   ##
   ## 'replacement' is expected to be a vector of replacement patterns, as
   ## used by gsub(), including relevant regular expression operators.
   ## If 'replacement' is empty, the "" is used, thereby replacing patterns with
   ## empty characters.
   ##
   ## replace_multiple is a logical vector indicating whether each pattern
   ## replacement should use gsub() if replaceMultiple==TRUE, or sub()
   ## if replaceMultiple==FALSE. The default is TRUE, which uses gsub().
   ## One would use replaceMultiple=FALSE in order to replace only the
   ## first occurrence of a pattern, like replacing the first tab character
   ## only.
   ##
   ## This function allows the patterns and replacements to be defined
   ## upfront, then applied to any relevant character vectors consistently,
   ## for example across columns of a data.frame.
   if (length(x) == 0 || length(pattern) == 0) {
      return(x);
   }
   if (length(replaceMultiple) == 0) {
      replaceMultiple <- TRUE;
   }
   replaceMultiple <- rep(replaceMultiple, length.out=length(pattern));
   if (length(replacement) == 0) {
      replacement <- "";
   }
   replacement <- rep(replacement, length.out=length(pattern));
   for (i in seq_along(pattern)) {
      if (replaceMultiple[[i]]) {
         x <- gsub(pattern=pattern[i],
            replacement=replacement[i],
            x=x,
            ignore.case=ignore.case,
            ...);
      } else {
         x <- sub(pattern=pattern[i],
            replacement=replacement[i],
            x=x,
            ignore.case=ignore.case,
            ...);
      }
   }
   return(x);
}

#' Find colname by character string or pattern matching
#'
#' Find colname by character string or pattern matching
#'
#' This function is intended to help find a colname given
#' a character vector of expected values or regular expression
#' patterns. By default it returns the first matching value,
#' but can return multiple if `max=Inf`.
#'
#' If there are no `colnames(x)` then `NULL` is returned.
#'
#' The order of operations:
#'
#' 1. Match exact string.
#' 2. Match exact string in case-insensitive manner.
#' 3. Match the start of each string using `jamba::provigrep()`.
#' 4. Match the end of each string using `jamba::provigrep()`.
#' 5. Match each string using `jamba::provigrep()`.
#'
#' The results from the first successful operation is returned.
#'
#' When there are duplicate `colnames(x)` only the first
#' unique name is returned.
#'
#' @family jam utility functions
#'
#' @return character vector with length `max`, or if no pattern
#'    match is found it returns `NULL`. Also if there are no
#'    `colnames(x)` then it returns `NULL`.
#'
#' @param pattern character vector containing text strings or
#'    regular expression patterns.
#' @param x input `data.frame` or other R object that
#'    contains colnames.
#' @param max integer maximum number of results to return.
#' @param index logical indicating whether to return the column
#'    index as an integer vector. When `index=FALSE` it returns
#'    the matching `colnames(x)`; when `index=TRUE` it returns
#'    the matching column numbers as an integer vector.
#' @param require_non_na logical indicating whether to require the
#'    column to contain non-NA values, default is TRUE. The intent
#'    of this function is to find colnames whose data will match
#'    expectations, and when require_non_na is TRUE, this
#'    function will continue until it finds a column with non-NA
#'    values.
#' @param ... additional arguments are passed to `jamba::provigrep()`.
#'
#' @export
find_colname <- function
(pattern,
 x,
 max=1,
 index=FALSE,
 require_non_na=TRUE,
 verbose=FALSE,
 ...)
{
   # handle list of data.frame or enrichResult objects
   # by iterating each objects and applying criteria
   if ("list" %in% class(x)) {
      x_colnames <- lapply(x, function(df) {
         find_colname(pattern=pattern,
            x=df,
            max=1,
            index=index,
            require_non_na=require_non_na,
            verbose=verbose)
      })
      x_vals <- head(unique(unlist(x_colnames)), max);
      return(x_vals);
   }
   #
   if ("enrichResult" %in% class(x)) {
      x <- x@result;
   }
   x_colnames <- colnames(x);
   ## require_non_na
   if (require_non_na) {
      x_colnames <- x_colnames[sapply(x_colnames, function(icol){
         any(!is.na(x[[icol]]))
      })]
   }
   ## if no colnames remain, return NULL
   if (length(x_colnames) == 0) {
      return(x_colnames);
   }

   start_pattern <- paste0("^", pattern);
   end_pattern <- paste0(pattern, "$");

   if (any(pattern %in% x_colnames)) {
      ## 1. max exact colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning exact match.");
      }
      x_vals <- intersect(pattern, x_colnames);
   } else if (any(tolower(pattern) %in% tolower(x_colnames))) {
      ## 2. max exact colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning exact case-insensitive match.");
      }
      x_match <- jamba::rmNA(match(tolower(pattern), tolower(x_colnames)));
      x_vals <- x_colnames[x_match];
      return(head(x_vals, max));
   } else if (jamba::igrepHas(paste(collapse="|", start_pattern), x_colnames)) {
      ## 3. match start of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to colname start.");
      }
      x_vals <- unique(jamba::provigrep(start_pattern, x_colnames));
      return(head(x_vals, max));
   } else if (jamba::igrepHas(paste(collapse="|", end_pattern), x_colnames)) {
      ## 4. match end of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to colname end.");
      }
      x_vals <- unique(jamba::provigrep(end_pattern, x_colnames));
      return(head(x_vals, max));
   } else if (jamba::igrepHas(paste(collapse="|", pattern), x_colnames)) {
      ## 5. match any part of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to part of colname.");
      }
      x_vals <- unique(jamba::provigrep(pattern, x_colnames));
      return(head(x_vals, max));
   } else {
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "No match found.");
      }
      x_vals <- NULL;
   }
   if (index && length(x_vals) > 0) {
      x_vals <- unique(match(x_vals, x_colnames));
   }
   return(head(x_vals, max));
}


#' Convert IPA list to a gene hit list or matrix
#' 
#' Convert IPA list to a gene hit list or matrix
#' 
#' Given a `list` of IPA data produced by `importIPAenrichment()`
#' produce a single `list` or incidence `matrix` representing
#' all gene hits tested by IPA.
#' 
#' @family jam import functions
#' 
#' @returns `list` by default, or `matrix` with `returnType="matrix"`.
#' 
#' @param IPAlist `list` of output from `importIPAenrichment()`
#' @param signColname `character` vector of patterns to match the colname
#'    in the 'Analysis Ready Molecules' IPA xref table which should
#'    indicate the direction of change where present.
#' @param signThreshold `numeric` minimum absolute threshold, default 0,
#'    for a value in 'signColname' to be considered a non-zero change.
#' @param geneColname `character` vector of patterns to match the colname
#'    in the 'Analysis Ready Molecules' IPA xref table which should
#'    indicate the gene involved. This column is validated by matching
#'    with the 'Canonical Pathways' geneID or geneNames column to ensure
#'    all values involved in pathway enrichment are also present in the
#'    xref column.
#' @param emptyValue `numeric` default 0, passed to `list2imSigned()` to
#'    use when `returnType="matrix"` in the matrix for non-hits.
#' @param returnType `character` string
#'    * 'list' (default): returns a signed gene list, with `integer` vectors
#'    named by gene.
#'    * 'matrix': returns a numeric `matrix` with gene rownames.
#' @param ... additional arguments are passed to internal functions
#' 
#' @export
IPAlist_to_hits <- function
(IPAlist,
 signColname=c("Expr Fold Change", "fold.*change",
    "log.*ratio", "log.*fold", "log.*fc", "lfc", "ratio", "fold", "fc"),
 signThreshold=0,
 geneColname=c("name", "symbol", "ID"),
 emptyValue=0,
 returnType=c("list", "matrix"),
 verbose=FALSE,
 ...)
{
   #
   returnType <- match.arg(returnType);
   
   # expand signThreshold
   signThreshold <- rep(signThreshold, length.out=length(IPAlist));
   names(signThreshold) <- names(IPAlist);
   
   #
   ipanames <- jamba::nameVectorN(IPAlist);
   
   # iterate each xref
   ipahitlist <- lapply(ipanames, function(iname){
      arm <- head(jamba::vigrep("analysis.*ready.*molecules",
         names(IPAlist[[iname]])), 1);
      if (length(arm) == 0) {
         stop_msg <- paste0("No 'Analysis Ready Molecules' found for ",
            iname, ".");
         stop(stop_msg);
      }
      can <- head(jamba::vigrep("canonical", names(IPAlist[[iname]])), 1);
      if (length(can) == 0) {
         stop_msg <- paste0("No 'Canonical Pathways' found for ",
            iname, ".");
         stop(stop_msg);
      }
      
      # define genes present in Canonical Pathways
      # to confirm the gene colname has appropriate values
      candf <- IPAlist[[iname]][[can]];
      canGeneColname <- find_colname(c("geneNames", "geneID"), candf);
      candf_genes <- unique(unlist(strsplit(candf[[canGeneColname]], "[/,]+")));
      
      # determine useGeneColname from IPA xref
      idf <- IPAlist[[iname]][[arm]];
      geneColnames <- jamba::provigrep(geneColname, colnames(idf));
      useGeneColname <- NULL;
      for (igeneColname in geneColnames) {
         testGenes <- gsub("[ ]*[(]includes others[)].*", "",
            gsub("/", ":", idf[[igeneColname]]));
         if (all(candf_genes %in% testGenes)) {
            useGeneColname <- igeneColname;
            break;
         }
      }
      if (length(useGeneColname) == 0) {
         stop("Genes in Canonical Pathways not present in IPA xref.");
      }
      
      # determine useSignColname with fold change column to use
      useSignColname <- find_colname(signColname, idf);
      idf_sign <- NULL;
      if (length(useSignColname) == 0) {
         warn_msg <- paste0("signColname did not match any IPA xref columns ",
            "for ", iname, ".");
         warning(warn_msg);
         # assume all are '1'
         idf_gene <- gsub(" [(]includes others.*", "",
            gsub("/", ":", idf[[useGeneColname]]));
         idf_sign <- rep(1, nrow(idf));
      } else {
         # optionally subset data
         if (signThreshold[[iname]] > 0) {
            idf <- subset(idf,
               abs(idf[[useSignColname]]) >= signThreshold[[iname]])
            if (nrow(idf) == 0) {
               return(NULL)
            }
         }
         idf_sign <- sign(idf[[useSignColname]]);
      }

      # generate signed list
      idf_gene <- gsub(" [(]includes others.*", "",
         gsub("/", ":", idf[[useGeneColname]]));
      jamba::nameVector(idf_sign, idf_gene)
   })
   if ("matrix" %in% returnType) {
      ipahitim <- list2imSigned(ipahitlist,
         emptyValue=emptyValue,
         ...)
      return(ipahitim)
   }
   return(ipahitlist)
}
