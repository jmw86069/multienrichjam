

#' Fix Set or pathway labels for legibility
#'
#' Fix Set or pathway labels for legibility
#'
#' This function is a convenient wrapper for several steps that edit
#' gene set and pathways labels to be slightly more legible.
#' It operates on:
#' * `character` vector, returning `character` vector
#' * `igraph` object, where is uses 'name' to update 'label'
#' * `Mem` object, where it updates `sets(Mem)`
#' 
#' The arguments have extensive default values encoded, which are
#' represented in data `multienrichjam::words` for basic word replacement,
#' and `multienrichjam::abbrev` for abbreviations used only
#' when `do_abbreviations=TRUE`.
#' 
#' Summary of typical changes:
#' 
#' * The vast majority of changes are custom biological terms which
#' are expected to have certain capitalization, for example
#' 'Mapk' is usually written 'MAPK'.
#' * Some changes are motivated to fix common artifacts in public data,
#' for example `'PI3kakt'` refers to `'PI3K/AKT'`.
#' * To use your own replacements, supply `words_from` as a two-column
#' `data.frame`, or two vectors `words_from` and `words_to`.
#' * To add custom effects to default, supply `abbrev_from` as a two-column
#' `data.frame`, or use two vectors `abbrev_from` and `abbrev_to`.
#' 
#' For `igraph` input, the vertex 'name' is used as the starting point.
#' To revert changes, use `igraph::V(x)$label <- igraph::V(x)$name`.
#' 
#' For `Mem` input, the `sets(x)` are updated, with no immediate way
#' to revert changes. It may become useful to do so in future, however.
#' 
#' @returns object whose class matches input 'x':
#'    * `character` vector
#'    * `igraph` object with vertex 'label' updated from 'name'
#'    * `Mem` object with updated `sets()`
#'
#' @family jam Mem utilities
#' @family jam igraph functions
#'
#' @param x any of the following objects:
#'    * `character` vector
#'    * `igraph` object. The `igraph::V(g)$name` attribute is used as input,
#'    and the resulting label is then stored as `V(g)$label`.
#'    When `nodeType` is also defined, and nodes have attribute 'nodeType',
#'    only nodes with that attribute value will be edited.
#'    The default is `nodeType="Set"`.
#'    * `Mem` object. The `sets()` are adjusted by this function.
#' @param wrap `logical` indicating whether to apply word wrap, based upon
#'    the supplied `width` argument.
#' @param width integer value used when `wrap=TRUE`, it is sent to
#'    `base::strwrap()`.
#' @param maxNchar `numeric` value or `Inf` to limit the maximum characters
#'    allowed for each string. This option is preferred when `wrap=TRUE`
#'    is not feasible, for example heatmap labels. When `NULL` or `Inf`
#'    no limit is applied. See `base::nchar()`.
#' @param suffix `character` value, default `"..."`, used
#'    when `maxNchar` is below `Inf`.
#'    When a string is shortened to `maxNchar`, the `suffix` helps indicate
#'    that there was additional text.
#' @param nodeType `character` string ussed when `x` is `igraph`,
#'    to limit changes to nodes by attribute values in `"nodeType"`.
#'    Use `"any"` or `NULL` to affect all nodes.
#' @param do_abbreviations `logical`, default TRUE, whether to apply
#'    `abbrev_from,abbrev_to`. These patterns are intended specifically
#'    to help shorten a long phrase, possibly removing words, or using
#'    common abbreviations.
#' @param adjustCase `logical`, default TRUE,  indicating whether to
#'    adjust the uppercase and lowercase lettering by calling
#'    `jamba::ucfirst()`. The default sets all characters to lowercase,
#'    then applies uppercase to the first letter of each word.
#' @param lowercaseAll `logical` used only when `adjustCase=TRUE`,
#'    passed to `jamba::ucfirst()`
#' @param removeGrep `character` regular expression pattern used to remove
#'    patterns from the resulting label.
#'    * The default removes common canonical pathway source prefix terms
#'    use in MSigDB data, for example KEGG, BIOCARTA, PID, etc.
#'    Use `""` or `NULL` to skip this step.
#'    * Multiple values can be defined, they are applied in order.
#' @param words_from,words_to `character` default NULL uses internal data
#'    `words` with pattern and replacement.
#'    * Input `words_from` can be a two-column `data.frame` expected to have
#'    'from' and 'to' in order.
#'    * Supplied as vectors, the 'words_from' are regular expression patterns,
#'    replaced with 'words_to' in order, applied case-sensitive.
#'    It does use Perl regular expression in `base::gsub()`, which is useful
#'    to use with 'backslash-b' to enforce a word boundary for example.
#' @param add_from,add_to `character` vectors used in addition to
#'    `words_from`,`words_to`.
#'    * 'add_from' can be supplied as a two-column `data.frame` as described
#'    for 'words_from'.
#'    * These values are applied after `words_from`,`words_to`, so that
#'    user-defined replacements have priority.
#' @param abbrev_from,abbrev_to `character` default NULL uses internal data
#'    `abbrev` with pattern and replacement.
#'    Intended to apply a specific abbreviation, and only applied when
#'    `do_abbreviations=TRUE`.
#'    * 'abbrev_from' can be supplied as a `data.frame` as described for
#'    'words_from'.
#'    * The abbreviations are "opinionated" in that they may remove words
#'    or shorten common phrases which do not seem critical to
#'    understanding the meaning of most biological pathways.
#'    
#'    Examples:
#'    * "Extracellular Matrix" becomes "ECM"
#'    * "Mitochondrial" becomes "Mito"
#'    * " Pathway" at the end of a phrase is removed, as it is not
#'    required to understand the rest of the label.
#'    * "Signaling by " at the start of a phrase is removed, as
#'    it also is not typically necessary to understand the label.
#' @param perl `logical` default TRUE, passed to `gsub()` for pattern matching.
#'    When Perl-mode is enabled, it also enforces word boundaries before and
#'    after each pattern. When Perl-mode is not enabled, there are no word
#'    boundary conditions applied.
#' @param ... additional arguments are passed to `jamba::ucfirst(x, ...)`,
#'    for example `firstWordOnly=TRUE` will capitalize only the first word.
#'
#' @examples
#' x <- c("KEGG_INSULIN_SIGNALING_PATHWAY",
#'    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#'    "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY");
#' fixSetLabels(x);
#' fixSetLabels(x, do_abbreviations=FALSE);
#'
#' jamba::nullPlot();
#' jamba::drawLabels(txt=x,
#'    preset=c("top", "center", "bottom"));
#'
#' @export
fixSetLabels <- function
(x,
 wrap=TRUE,
 width=40,
 maxNchar=Inf,
 suffix="...",
 nodeType=c("Set",
    "Gene",
    "any"),
 do_abbreviations=TRUE,
 adjustCase=TRUE,
 lowercaseAll=TRUE,
 removeGrep="^(KEGG(_MEDICUS|)(_REFERENCE|_VARIANT|)|PID|REACTOME|BIOCARTA|NABA|SA|SIG|ST|WP|HALLMARK)[_. ]",
 words_from=NULL,
 words_to=NULL,
 add_from=NULL,
 add_to=NULL,
 abbrev_from=NULL,
 abbrev_to=NULL,
 perl=TRUE,
 ...)
{
   # validate nodeType
   nodeType <- match.arg(nodeType);
   
   # use default words
   if (length(words_from) == 0) {
      words_from <- words$from;
      words_to <- words$to;
   }
   # use default abbrev
   if (length(abbrev_from) == 0) {
      abbrev_from <- multienrichjam::abbrev$from;
      abbrev_to <- multienrichjam::abbrev$to;
      # jamba::printDebug("abbrev:");print(data.frame(abbrev_from, abbrev_to));# debug
   }
   # accept data.frame
   if (inherits(words_from, "data.frame")) {
      if (ncol(words_from) != 2) {
         stop("words_from as data.frame must have two columns.")
      }
      if (all(c("from", "to") %in% colnames(words_from))) {
         words_to <- words_from[["to"]];
         words_from <- words_from[["from"]];
      } else {
         words_to <- words_from[[2]];
         words_from <- words_from[[1]];
      }
   }
   # accept data.frame
   if (inherits(add_from, "data.frame")) {
      if (ncol(add_from) != 2) {
         stop("add_from as data.frame must have two columns.")
      }
      if (all(c("from", "to") %in% colnames(add_from))) {
         add_to <- add_from[["to"]];
         add_from <- add_from[["from"]];
      } else {
         add_to <- add_from[[2]];
         add_from <- add_from[[1]];
      }
   }
   # accept data.frame
   if (inherits(abbrev_from, "data.frame")) {
      if (ncol(abbrev_from) != 2) {
         stop("abbrev_from as data.frame must have two columns.")
      }
      if (all(c("from", "to") %in% colnames(abbrev_from))) {
         abbrev_to <- abbrev_from[["to"]];
         abbrev_from <- abbrev_from[["from"]];
      } else {
         abbrev_to <- abbrev_from[[2]];
         abbrev_from <- abbrev_from[[1]];
      }
   }
   
   # validate length for from,to pairs
   if (length(words_from) != length(words_to)) {
     stop(paste0(
       "length(words_from) does not equal length(words_to).",
       "Please remedy."));
   }
   if (length(add_from) != length(add_to)) {
     stop(paste0(
       "length(add_from) does not equal length(add_to).",
       "Please remedy."));
   }
   if (length(abbrev_from) != length(abbrev_to)) {
     stop(paste0(
       "length(abbrev_from) does not equal length(abbrev_to).",
       "Please remedy."));
   }

   xMem <- NULL;
   if (inherits(x, "Mem")) {
      xMem <- x;
      x <- sets(xMem);
      which_nodes <- seq_along(x);
      for (i in seq_along(removeGrep)) {
         xPrep <- gsub("[_ ]+", " ",
            gsub(removeGrep[[i]],
               "",
               ignore.case=TRUE,
               x));
      }
   } else if (inherits(x, "igraph")) {
      if ("any" %in% nodeType) {
         which_nodes <- seq_len(igraph::vcount(x));
      } else {
         which_nodes <- which(igraph::V(x)$nodeType %in% nodeType);
      }
      for (i in seq_along(removeGrep)) {
         xPrep <- gsub("[_ ]+", " ",
            gsub(removeGrep[[i]],
               "",
               ignore.case=TRUE,
               igraph::V(x)$name[which_nodes]));
      }
   } else if (is.atomic(x)) {
      which_nodes <- seq_along(x);
      for (i in seq_along(removeGrep)) {
         xPrep <- gsub("[_ ]+", " ",
            gsub(removeGrep[[i]],
               "",
               ignore.case=TRUE,
               x));
      }
   } else {
      stop("Input must be atomic vector, or 'igraph', or 'Mem'.")
   }
   if (TRUE %in% adjustCase) {
      xPrep <- jamba::ucfirst(xPrep,
         lowercaseAll=lowercaseAll,
         ...);
   }

   ## Confirm words_from,words_to have identical length
   if (length(words_from) > 0) {
      words_to <- rep(words_to, length.out=length(words_from))
   }
   ## When add_from,add_to are defined
   if (length(add_from) > 0) {
      add_to <- rep(add_to, length.out=length(add_from))
      # append to words_from
      words_from <- c(words_from, add_from)
      words_to <- c(words_to, add_to)
   }
   ## When abbrev_from,abbrev_to are defined
   if (isTRUE(do_abbreviations) && length(abbrev_from) > 0) {
      abbrev_to <- rep(abbrev_to, length.out=length(abbrev_from))
      # trim trailing space in pattern match
      abbrev_from <- gsub("[ ]*$", "", abbrev_from);
      # append to words_from
      words_from <- c(words_from, abbrev_from)
      words_to <- c(words_to, abbrev_to)
   }
   
   # jamba::printDebug("words:");print(data.frame(words_from, words_to));# debug
   
   # Apply replacements, case-insensitive match, case-sensitive replacement
   if (length(words_from) > 0 && length(words_to) == length(words_from)) {
      for (i in seq_along(words_from)) {
      	if (isTRUE(perl)) {
      		xPrep <- gsub(paste0("\\b", words_from[i], "\\b"),
	            words_to[i],
	            ignore.case=TRUE,
	            perl=TRUE,
	            xPrep);
      	} else {
      		xPrep <- gsub(paste0("", words_from[i], ""),
      			words_to[i],
      			ignore.case=TRUE,
      			perl=FALSE,
      			xPrep);
      	}
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
   
   ## Apply base::trimws() to trim leading/trailing whitespace
   xPrep <- base::trimws(xPrep, which="both")
   
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
   if (inherits(xMem, "Mem")) {
      sets(xMem) <- xNew;
      return(xMem);
   } else if (inherits(x, "igraph")) {
      if (!"label" %in% igraph::vertex_attr_names(x)) {
         igraph::V(x)$label <- igraph::V(x)$name;
      }
      igraph::V(x)[which_nodes]$label <- xNew;
   } else {
      x <- xNew;
   }
   return(x);
}
