

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
