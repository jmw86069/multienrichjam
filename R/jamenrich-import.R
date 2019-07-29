
#' Import Ingenuity IPA enrichment results
#'
#' Import Ingenuity IPA enrichment results
#'
#' This function parses Ingenuity IPA enrichment data into
#' a form usable as a list of  enrichment `data.frame`
#' objects for downstream analysis. Each `data.frame`
#' will represent the results of one Ingenuity IPA
#' enrichment test.
#'
#' The input data can be one of four forms:
#'
#' 1. `ipaFile` can be a text `.txt` file,
#' where the text file contains all IPA enrichment data in
#' tall format. This format is most common.
#' 2. `ipaFile` can be an Excel `.xlsx` file,
#' which contains all IPA enrichment data in
#' one tall worksheet tab.
#' 3. `ipaFile` can be an Excel `.xlsx` file,
#' where each type of IPA enrichment appears on a separate
#' Excel worksheet tab.
#' 4. `ipaFile` can be a list of `data.frame` objects.
#' This option is intended when the IPA data has already
#' been imported into R as separate `data.frame` objects.
#'
#' The basic motivation for this function is two-fold:
#'
#' 1. Separate multiple IPA enrichment tables.
#' 2. Rename colnames to be consistent.
#'
#' When using `"Export All"` from IPA, the default text format
#' includes multiple enrichment tables concatenated together in one
#' file. Each enrichment table contains its own unique column
#' headers, with descriptive text in the line preceding the
#' column headers. This function is intended to separate the
#' enrichment tables into a list of `data.frame` objects, and
#' retain the descriptive text as names of the list.
#'
#' @return list of `data.frame` objects, where each `data.frame`
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
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
importIPAenrichment <- function
(ipaFile,
 headerGrep="(^|\t)((expr.|-log.|)p-value|Pvalue|Score\t|Symbol\t|Ratio\t|Consistency.Score|Master.Regulator\t)",
 ipaNameGrep=c("Pathway", "Regulator$", "Regulators", "Regulator", #"Master.Regulator",
    "Disease", "Toxicity",
    "Category", "Categories",
    "Function", "Symbol$",
    "^ID$", "My.(Lists|Pathways)"),
 geneGrep=c("Molecules in Network", "Molecules"),
 geneCurateFrom=c("^[,]+|[,]+$"),
 geneCurateTo=c(""),
 method=1,
 sheet=1,
 sep="\t",
 xlsxMultiSheet=TRUE,
 useXlsxSheetNames=FALSE,
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
      if (suppressPackageStartupMessages(!require(openxlsx))) {
         stop("importIPAenrichment() requires the openxlsx package for Excel import.");
      }
   }
   if (suppressPackageStartupMessages(!require(jamba))) {
      stop("importIPAenrichment() requires the jamba package, devtools::install_github('jmw86069/jamba')");
   }

   if (jamba::igrepHas("list", class(ipaFile)) &&
         jamba::igrepHas("data.frame|tibbletbl|matrix|dataframe", class(ipaFile[[1]]))) {
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
            verbose=verbose);
         jDF;
      });
   } else if (jamba::igrepHas("[.]xlsx$", ipaFile) && length(ipaFile) > 1) {
      ## Process multiple files by calling this function on each file
      ipaDFLL <- lapply(ipaFile, function(ipaFile_i){
         jDF <- importIPAenrichment(ipaFile=ipaFile_i,
            headerGrep=headerGrep,
            ipaNameGrep=ipaNameGrep,
            geneGrep=geneGrep,
            geneCurateFrom=geneCurateFrom,
            geneCurateTo=geneCurateTo,
            method=method,
            sheet=sheet,
            sep=sep,
            xlsxMultiSheet=xlsxMultiSheet,
            verbose=verbose,
            ...);
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
               verbose=verbose);
            jDF;
         } else {
            jDF <- importIPAenrichment(ipaFile=ipaFile,
               headerGrep=headerGrep,
               ipaNameGrep=ipaNameGrep,
               geneGrep=geneGrep,
               geneCurateFrom=geneCurateFrom,
               geneCurateTo=geneCurateTo,
               method=method,
               sheet=sheet,
               sep=sep,
               xlsxMultiSheet=FALSE,
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
            pasteByRow(iDF,
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
         i <- system(intern=TRUE,
            paste0("grep '.*' ", ipaFile));
         ## Clean up trailing newlines
         i <- gsub("[\r\n]+$",
            "",
            i);
         i;
      }
      ## Often IPA output has descriptor lines with 'for ... ->' in them
      ## which indicates a header for each subtable. We will parse
      ## that header for their correct table name where possible.
      descLines <- grep("( for ).*(->)", i);
      if (length(descLines) > 0) {
         iDesc <- gsub("( for |->).*$", "", i[descLines]);
      }

      ## Infer which columns area headers, then set up
      ## ranges of rows to be assigned to each header.
      i1 <- jamba::igrep(headerGrep, i);
      if (length(i1) == 0) {
         stop("The worksheet does not have rows matching headerGrep.");
      }
      i2 <- c(tail(i1, -1) - 1, length(i));
      i12 <- cbind(i1, i2);
      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "i1:", i1);
         jamba::printDebug("importIPAenrichment(): ",
            "i2:", i2);
      }

      ## Pull out reasonable rownames based upon IPA naming conventions
      i12prevRow <- (i1-1);
      i12prevRowText <- i[i1-1];
      i12names1 <- jamba::makeNames(
         sapply(i[i1], function(ix){
            head(jamba::provigrep(
               ipaNameGrep,
               unlist(strsplit(ix, sep))), 1);
         })
      );
      i12names2 <- jamba::makeNames(gsub(
         paste0(".*(^|", sep, ")",
            "([^", sep, "]*",
            "(Pathways|Regulator[s]*|Diseases|Lists|Consistency.Score|Symbol)",
            "[^", sep, "]*)",
            "($|", sep, ").*"),
         "\\2",
         i[i12[,1]]));
      if (any(i12prevRow %in% descLines)) {
         iWhich <- which(i12prevRow %in% descLines);
         iWhichMatch <- match(i12prevRow[iWhich], descLines);
         i12names1[iWhich] <- iDesc[iWhichMatch];
      }
      if (verbose) {
         print(head(as.data.frame(do.call(cbind, list(i12names1=i12names1, i12names2=i12names2))), Inf));
         jamba::printDebug("i12names1:");
         print(i12names1);
         jamba::printDebug("i12names2:");
         print(i12names2);
      }
      rownames(i12) <- i12names1;
      i12;
      ## Remove any headers which have no rows of results
      i12 <- i12[which(i1 != i2),,drop=FALSE];

      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "detected these ranges of rows to import:");
         print(i12);
      }


      ## Iterate each set of results and create a data.frame.
      ## Note that each type of result has its own number of columns,
      ## and unique colnames.
      ipaDFL <- lapply(jamba::nameVector(rownames(i12)), function(j){
         if (verbose) {
            jamba::printDebug("importIPAenrichment(): ",
               "creating data.frame:",
               j);
         }
         j1 <- i12[j,1];
         j2 <- i12[j,2];
         ## Create a sequence that removes the descLines as relevant
         jseq <- setdiff(seq(from=j1, to=j2, by=1),
            descLines);
         ## Note it is very important to use textConnection() here
         ## because it properly preserves character encodings, such
         ## as alpha, beta characters in some pathway names from IPA.
         if (method == 1) {
            if (verbose) {
               jamba::printDebug("importIPAenrichment(): ",
                  "using read.table().");
            }
            jDF <- read.table(
               #textConnection(gsub("\t+$", "", i[j1:j2])),
               #textConnection(gsub("\t+$", "", i[jseq])),
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
            if (suppressPackageStartupMessages(!require(readr))) {
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

         jDF <- curateIPAcolnames(jDF,
            ipaNameGrep=ipaNameGrep,
            geneGrep=geneGrep,
            geneCurateFrom=geneCurateFrom,
            geneCurateTo=geneCurateTo,
            verbose=verbose);
         jDF;
      });
   }
   ## Return list of data.frames
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
#' @family jam import functions
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
 ipaNameGrep,
 geneGrep,
 geneCurateFrom,
 geneCurateTo,
 verbose=TRUE,
 ...)
{
   ## Purpose is to curate some colnames in IPA enrichment data
   ## Determine a column to become the "Name"
   if (suppressPackageStartupMessages(!require(jamba))) {
      stop("curateIPAcolnames() requires the jamba package.");
   }
   nameCol <- head(jamba::provigrep(ipaNameGrep, colnames(jDF)), 1);
   if (length(nameCol) == 1) {
      jDF[["Name"]] <- jDF[[nameCol]];
      ## Skip rename, in order to maintain the original colname
      #jDF <- renameColumn(jDF,
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
      jDF <- renameColumn(jDF,
         from=geneCol,
         to="geneNames");
      if (verbose) {
         jamba::printDebug("importIPAenrichment(): ",
            "created ", "'geneNames'", " column from:",
            geneCol);
      }
      ## Curate gene column values
      jDF[["geneNames"]] <- gsubs(geneCurateFrom,
         geneCurateTo,
         jDF[["geneNames"]]);
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
         jDF <- renameColumn(jDF,
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
gsubs <- function
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
#' @param ... additional arguments are passed to `jamba::provigrep()`.
#'
#' @export
find_colname <- function
(pattern,
 x,
 max=1,
 index=FALSE,
 verbose=FALSE,
 ...)
{
   ##
   x_colnames <- colnames(x);
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
