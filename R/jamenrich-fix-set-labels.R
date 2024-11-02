

#' Fix Set labels for legibility
#'
#' Fix Set labels for legibility
#'
#' This function is a convenient wrapper for several steps that edit
#' gene set and pathways labels to be slightly more legible. It
#' operates on either a character vector, or an igraph object.
#'
#' * To use custom from,to replacements, along with the default replacements,
#' supply the custom replacements with arguments `add_from`,`add_to`.
#' * To use custom from,to replacements, without applying the defaults,
#' supply the custom replacements with arguments `words_from`,`words_to`.
#'
#' @return vector or igraph object, to match the input `x`.
#'
#' @family jam igraph functions
#'
#' @param x any of the following objects:
#'    * `character` vector
#'    * `igraph` object. The `V(g)$name` attribute is used as the input,
#'    and the resulting label is then stored as `V(g)$label`.
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
#' @param words_from,words_to `character` vectors of pattern, replacement,
#'    respectively. The pattern is matched in case-insensitive manner,
#'    with case-sensitive replacements where applicable.
#'    It uses perl-based regular expression matching with
#'    `base::gsub()`, so that the expression `\\b` can be used
#'    to enforce a word boundary, either via delimiter, whitespace, or
#'    the end of the string.
#' @param add_from,add_to `character` vectors used in addition to
#'    `words_from`,`words_to`.
#'    * These values are applied after `words_from`,`words_to`, so that
#'    user-defined replacements have priority.
#' @param abbrev_from,abbrev_to `character` vectors used when
#'    `do_abbreviations=TRUE`. These defaults are "opinionated",
#'    they are intended to shorten common phrases which do not
#'    seem critical to understanding the meaning of most biological
#'    pathways.
#'    Some abbreviations are used for relatively common phrases
#'    and terms, for which the abbreviation seems to be unambiguous
#'    and fairly widely recognized.
#'    Examples:
#'    * "Extracellular Matrix" becomes "ECM"
#'    * "Mitochondrial" becomes "Mito"
#'    * " Pathway" at the end of a phrase is removed, as it is not
#'    required to understand the rest of the label.
#'    * "Signaling by " at the start of a phrase is removed, as
#'    it also is not typically necessary to understand the label.
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
 removeGrep="^(KEGG|PID|REACTOME|BIOCARTA|NABA|SA|SIG|ST|WP|HALLMARK)[_. ]",
 words_from=c("als",
    "ii", "iii", "iv", "v",  "vi", "Vii", "Viii", "ix", "x",
    "trna", "rrna", "rna", "dna", "mirna", "mrna", "snrna", "snorna",
    "scrna", "lincrna",
    "Il", "Ecm", "Nk",
    "Pi3k.Akt", "Akt", "Pi3k", "tgf", "nfkb", "NK.Kappa.B",
    "Pi3kaktmtorsignaling", "Pi3kaktmtor",
    "PI3kakt", "aktmtor", "mtorsignaling",
    "Pi3kci", "Pi3kgamma", "Pi 3k",
    "Ppar(alpha|a)", "Ppar(gamma|g)", "Ppar",
    "Udp N Acetyl Glucosamine",
    "Mtor", "Gpcr", "Gpcrs",
    "Tnfa|Tnfalpha", "Tnfr1", "Tnfr2", "Tnfs", "Tnf", "Tnfsf", "Tnfrelated",
    "Tgf(beta|b)", "Tgfbr", "TGF Beta",
    "Akt1", "Igf1akt",
    "Foxo", "Hdacs", "Hdac", "Hat", "Hats",
    "Nicotinic Acetylcholine Receptors|Acetylcholine Nicotinic Receptors",
    "G CSF", "Gm Csf", "M Csf",
    "IL([0-9]+)[-. _]([0-9]+)pathway",
    "Jak Stat([1-9]*)", "Mkk3 6pathway",
    "MAP([23]*)K([0-9]*)", "P38[ ]*MAPK", "MAPK1 3", "Erk MAPK", "Erk", "Erks",
    "Interferon([a-z]*)",
    "([A-Za-z0-9]+[a-qs-z])mediated",
    "Lncrna|Lincrna", "Microrna([s]*)",
    "125 DIHYDROXYVITAMIN D3",
    "Tweak", "Er",
    "Vegfavegfr2", "Vegfr1 2", "Vegf", "Vegfr([1-9]*)", "Egfr", "Egfrviii",
    "Egfegfr",
    "Smad2 3(pathway|nuclear)", "Smad2 3",
    "Nf Kb", "Nfkappab",
    "C Jun", "C Fos", "Ap 1", "Apoe",
    "([AG])tpase", "Kras", "Uv",
    "Sarscov2", "Sars Cov ([12])", "Sars Cov", "Covid19",
    "Cellspecific", "([a-zA-Z0-9]+)like",
    "([TB]|NK) Cell"),
 words_to=c("ALS",
    "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "tRNA", "rRNA", "RNA", "DNA", "miRNA", "mRNA", "snRNA", "snoRNA",
    "scRNA", "lincRNA",
    "IL", "ECM", "NK",
    "PI3K/AKT", "AKT", "PI3K", "TGF", "NFKB", "NFKB",
    "PI3K/AKT/mTOR Signaling", "PI3K/AKT/mTOR",
    "PI3K/AKT", "AKT/mTOR", "mTOR Signaling",
    "PI3KCI", "PI3Kgamma", "PI3K",
    "PPARalpha", "PPARgamma", "PPAR",
    "UDP-GlcNAc",
    "mTOR", "GPCR", "GPCRs",
    "TNFa", "TNFR1", "TNFR2", "TNFs", "TNF", "TNFSF", "TNF-related",
    "TGFbeta", "TGFbeta-Receptor", "TGFbeta",
    "AKT1", "IGF1/AKT",
    "FOXO", "HDACs", "HDAC", "HAT", "HATs",
    "Nicotinic Acetylcholine Receptors", #"nAChRs",
    "G-CSF", "GM-CSF", "M-CSF",
    "IL-\\1/IL-\\2 Pathway",
    "JAK/STAT\\1", "MKK3/MKK6 Pathway",
    "MAP\\1K\\2", "p38-MAPK", "MAPK1/MAPK3", "ERK-MAPK", "ERK", "ERKs",
    "IFN\\1",
    "\\1-Mediated",
    "lncRNA", "microRNA\\1",
    "1,25-dihydroxyvitamin D3",
    "TWEAK", "ER",
    "VEGFA/VEGFR2", "VEGFR1/VEGFR2", "VEGF", "VEGFR\\1", "EGFR", "EGFRvIII",
    "EGF/EGFR",
    "Smad2/3 \\1", "Smad2/3",
    "NFKB", "NFKB",
    "C-jun", "C-fos", "AP-1", "APOE",
    "\\1TPase", "KRAS", "UV",
    "SARS-CoV-2", "SARS-CoV-\\1", "SARS-CoV", "COVID19",
    "Cell-Specific", "\\1-Like",
    "\\1-cell"),
 add_from=NULL,
 add_to=NULL,
 abbrev_from=c(
    "Extracellular.Matrix",
    "Mitochondri(um|a|al|on)",
    "Interferon",
    "(IL|Interleukin[ ]*)([0-9]+)",
    "Subsequent",
    "Signaling (pathway|system)",
    "Of The",
    "^Signaling by ", # remove leading "Signaling by" as unnecessary
    " Pathway[s]*$", # remove trailing " pathway" as unnecessary
    "Expression",
    "The Role "),
 abbrev_to=c(
    "ECM",
    "Mito",
    "IFN",
    "IL-\\2",
    "",
    "Signaling",
    "Of",
    "",
    "",
    "Expr.",
    "Role "),
 ...)
{
   # validate nodeType
   nodeType <- match.arg(nodeType);
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

   if (inherits(x, "igraph")) {
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
   } else {
      which_nodes <- seq_along(x);
      for (i in seq_along(removeGrep)) {
         xPrep <- gsub("[_ ]+", " ",
            gsub(removeGrep[[i]],
               "",
               ignore.case=TRUE,
               x));
      }
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
   if (TRUE %in% do_abbreviations && length(abbrev_from) > 0) {
      abbrev_to <- rep(abbrev_to, length.out=length(abbrev_from))
      # append to words_from
      words_from <- c(words_from, abbrev_from)
      words_to <- c(words_to, abbrev_to)
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
