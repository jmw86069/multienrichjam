
#' Multi-enrichment analysis with multienrichjam
#' 
#' @details
#' This package aims to enable analysis for one or multiple
#' pathway or functional enrichment datasets, with flexible
#' options for publication of data visualizations as figures.
#' 
#' The core workflow follows these steps:
#' 
#' 1. Import data for use.
#' 2. Call `multiEnrichMap()` using relevant parameters.
#' 3. Visualize with `mem_plot_folio()` for review.
#' 4. Iterate then polish final figures.
#' 
#' The analysis produces:
#' 
#' * `Mem-class`: multi-enrichment data used for visualization.
#' * Data visualizations from `mem_plot_folio()` in `list` format.
#' 
#' Concept network (Cnet) plots can be customized using ShinyCat:
#' 
#' * `launch_shinycat()` - interactive network adjustments
#' 
#' @docType package
#' @keywords internal
"_PACKAGE"
