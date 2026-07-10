
#' Plot function for MemPlotFolio objects
#'
#' Plot function for MemPlotFolio objects
#' 
#' Main purpose is to permit independent plotting
#' for MemPlotFolio objects, which should allow custom
#' settings for things like Rmd or Qmd tabs, or for
#' multi-page PDF output.
#' 
#' Note argument `plot_which` takes `character` vector,
#' to be distinct from `plot_which` which takes `integer`
#' as used historically by `mem_plot_folio()`.
#' 
#' 
#' Todo:
#' 
#' * Future work will allow customization of each plot.
#' 
#' ## Markdown tab output
#' 
#' When `do_md_tabs=TRUE` the default is to autodetect
#' appropriate values for `md_tab_open` and `md_tab_close`.
#' These values are defined if running in knitr context,
#' with HTML output, using recognized input file type
#' 'Rmd', 'rmd', 'Qmd', 'qmd'. In all other cases, it sets
#' tabs to empty character ''.
#' 
#' You can set `do_md_tabs='Rmd'` to force Rmarkdown tabset style,
#' or `do_md_tabs='Qmd'` to force Quarto tabset style.
#' In either case, the tabs are defined using defaults,
#' where: Rmd uses '{.tabset}' and no footer; and
#' Qmd uses '::: {.panel-tabset}' and ':::' footer.
#' 
#' Setting `do_md_tabs=TRUE` then all other values `''` blank
#' will print markdown-style headings for each plot type,
#' without adding tabset tags. This output may be preferred
#' for a "long form" report where each plot has its own
#' distinct header, which may also appear in a
#' table of contents.
#' 
#' Secondary goals:
#' 
#' * Decouple all the sub-plot customizations
#' so they can be applied without re-creating the object.
#' * Permit subsetting the Mpf data as relevant, perhaps by
#' `Clusters()` as a drill-down technique.
#' 
#' Active development:
#' 
#' * Each plot will support more customization options over time,
#' particularly to adjust Cnet network layouts, and other
#' custom labeling options.
#' * Notably, support for `ClusterLabels()` is propagating
#' through various plots, and `ClusterData()` may provide
#' interactive capability, for example detailed mouse-over
#' or event-driven data for each pathway or gene cluster.
#' 
#' @returns `list` plot data, named by type, returned invisibly.
#' 
#' @family custom plot functions
#' 
#' @param Mpf `MemPlotFolio` object as returned by
#'    `prepare_folio()` or `mem_plot_folio()`.
#' @param plot_which `character` vector with one or more plots,
#' default 'all' renders all plots defined in 'Mpf'.
#' Plots will be created in the order provided, except for 'all'
#' which follows the default order.  
#' The following values are recognized:
#'    * 'all': All available plot types are included.
#'    * 'EnrichmentHeatmap' or 'eh': `EnrichmentHeatmap()`
#'    * 'GenePathHeatmap' or 'gp': `GenePathHeatmap()`
#'    * 'CnetCollapsed' or 'cc': `CnetCollapsed()`
#'    * 'CnetExemplar' or 'ce': `CnetExemplar()`
#'    * 'CnetCluster' or 'c': `CnetCluster()`
#'    * 'EnrichmentMap' or 'em': `mem2emap()`. Note this plot
#'    requires either 'Mem' argument, or 'Mem' being included
#'    in the `Mpf@metadata$Mem`.
#' @param do_md_tabs `logical` default FALSE, whether to print
#'    markdown-compatible tab headers before rendering each
#'    plot, done using `cat()` to STDOUT.
#'    To enable specific Rmarkdown or Quarto (Qmd) style,
#'    provide a `character` string, which will override all
#'    other checks, and will define 'md_tab_open' and 'md_tab_close'.
#'    * 'Rmd' or 'rmd' to enable Rmarkdown '{.tabset}' style
#'    * 'Qmd' or 'qmd' or 'Quarto' or 'quarto' to enable
#'    Quarto '::: {.panel-tabset}' style, with footer ':::'.
#' @param md_tab_open `character` string to define the tab open
#'    style. When NULL, it will auto-detect an appropriate style.
#'    * Rmarkdown with HTML output: `''`
#'    * Quarto with HTML output: `'\\n\\n::: {.panel-tabset}\\n\\n'`.
#'    * All others: `''`.
#' @param md_tab_suffix `character` string to define the tab suffix
#'    style. When NULL, it will auto-detect an appropriate style.
#'    * Rmarkdown with HTML output: `' {.tabset}\\n\\n'`
#'    * All others: `''`.
#' @param md_tab_close `character` string to define the tab close
#'    style. When NULL, it will auto-detect an appropriate style.
#'    * Rmarkdown with HTML output: `''`
#'    * Quarto with HTML output: `'\\n\\n:::\\n\\n'`.
#'    * All others: `''`
#' 
#'    In Quarto, to keep the tabset open, define `md_tab_close=''`.
#' @param md_tab_level `integer` heading level to begin, default 2
#'    uses heading '##'.
#' @param md_title `character` default 'Multi-Enrichment Folio' will
#'    print a title heading using `md_tab_level`, then create the tabset
#'    underneath at a lower heading.
#'    * Use md_title=FALSE or md_title=NULL to suppress this header.
#' @param cc_type,ce_num,c_cluster passed to corresponding functions,
#'    and supports multiple values.
#'    * cc_type: passed to CnetCollapsed(Mpf, type=cc_type).
#'    When NULL it will use 'cluster' if `ClusterLabels(Mpf)` is
#'    available, otherwise 'set'. If more than 500 nodes, it appends
#'    '2' to the end, which hides gene node labels.
#'    * ce_num: passed to CnetExemplar(Mpf, num=ce_num). When NULL
#'    it will use ce_num=1.
#'    * c_cluster: passed to CnetCluster(Mpf, cluster=c_cluster).
#'    When NULL it will iterate each cluster in Mpf.
#' @param em_group `character` default NULL uses 
#'    community detection in `mem2emap()`.
#'    * 'community': uses community detection in `mem2emap()`.
#'    * 'clusters': uses `Clusters(Mpf)`.
#'    * 'cluster_labels': uses `ClusterLabels(Mpf)`.
#' @param Mem `Mem` default NULL, only included for backwards support of
#'    Mpf objects which do not have 'Mem' included in metadata.
#'    Experimental.
#' @param params `list` named by shorthand plot type
#'    ('eh','gp','em', 'cc', 'ce', 'c') each containing a list of
#'    optional parameters relevant to each plot.
#'    In each case, parameters in `params` are passed to the
#'    internal function as arguments.For example,
#'    `em=list(repulse=3.5)` is passed as
#'    `EnrichmentMap(..., repulse=3.5)`.
#'    Limited functionality currently, however check here for
#'    recognized arguments as they become available.
#'    * 'em' EnrichmentMap:
#' 
#'       * repulse: apply via `layout_with_qfr(g, repulse=repulse)`,
#'       default 3.5.
#'       * width: `integer` word-wrap character width, default 30.
#'       Word-wrap is applied using `fixSetLabels()`
#'       however it uses `do_abbreviations=FALSE` and `removeGrep=NULL`
#'       so the primary effect is to apply word-wrap, and capitalization.
#' @param do_newpage `logical` default NULL sets TRUE when knitr is running,
#'    whether to call `grid::grid.newpage()` after certain plots to
#'    encourage Rmd and Qmd to recognize end of a plot for the purpose
#'    of placing figures inside the correct tabset regions. Experimental.
#'    In some cases, knitr appears to embed a figure only when it
#'    detects that figure is "complete", and sometimes it does so after
#'    the next markdown tabset has been printed to STDOUT via `cat()`.
#'    Using `grid::grid.newpage()` appears to force knitr to recognize
#'    a new plot has started, therefore embeds the previous figure image.
#'    Also, knitr skips blank plots, so calling newpage, then having
#'    another plot call a newpage is not problematic.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
plot_mpf <- function
(Mpf,
 plot_which=c("all",
    "EnrichmentHeatmap", "eh",
    "GenePathHeatmap", "gp",
    "CnetCollapsed", "cc",
    "CnetExemplar", "ce",
    "CnetCluster", "c",
    "EnrichmentMap", "em"),
 do_md_tabs=FALSE,
 md_tab_open=NULL,
 md_tab_suffix=NULL,
 md_tab_close=NULL,
 md_tab_level=2,
 md_title="Multi-Enrichment Folio",
 cc_type=NULL,
 ce_num=NULL,
 c_cluster=NULL,
 em_group=NULL,
 Mem=NULL,
 params=list(
   eh=list(do_plot=TRUE),
   em=list(repulse=3.5,
      node_factor=1,
      width=30),
   cc=list(repulse=3.5,
      node_factor=1),
   ce=list(node_factor=1),
   "c"=list(node_factor=1)
 ),
 do_newpage=NULL,
 verbose=FALSE,
 ...)
{
   # context
   is_knitting <- isTRUE(getOption("knitr.in.progress"));
   # validate arguments
   params <- modifyList(
      eval(formals(plot_mpf)$params),
      params)
   
   # plot_which <- match.arg(plot_which, several.ok=TRUE)
   if (length(do_newpage) == 0) {
      do_newpage <- isTRUE(is_knitting);
   }
   # verbose remove '##' comments to prevent markdown headings
   if (verbose && (isTRUE(do_md_tabs) || nchar(do_md_tabs) > 0)) {
      withr::local_options(list(
         jam.comment=FALSE));
   }
   plot_list <- list();

   # define recognized values in do_which
   plot_which_set <- c(
      eh="EnrichmentHeatmap",
      EnrichmentHeatmap="EnrichmentHeatmap",
      gp="GenePathHeatmap",
      gphm="GenePathHeatmap",
      GenePathHeatmap="GenePathHeatmap",
      cc="CnetCollapsed",
      CnetCollapsed="CnetCollapsed",
      ce="CnetExemplar",
      CnetExemplar="CnetExemplar",
      "c"="CnetCluster",
      CnetCluster="CnetCluster",
      em="EnrichmentMap",
      EnrichmentMap="EnrichmentMap")
   # em alternative that uses clusters instead of communities

   if ("all" %in% plot_which) {
      # do them all
      plot_which <- unique(plot_which_set[
         c("eh", "em", "gp", "cc", "ce", "c", "em")]);
   } else {
      # convert to recognized plot types
      plot_which <- unique(plot_which_set[plot_which])
   }
   if (verbose) {
      jamba::printDebug("plot_mpf(): ",
         "plot_which: ", plot_which);
   }
   
   # Iterate plot types stored in Mpf
   # Optionally create (or re-create) plot types using stored params,
   # for example keep clustering output constant.

   # Determine type of markdown tabset to use
   # - if running knitr, HTML output, Rmd file
   #   --> '{.tabset}', footer ''
   # - if running knitr, HTML output, Qmd file 
   #   --> '::: {.panel-tabset}', footer ':::'
   if (do_md_tabs %in% c("Rmd", "rmd")) {
      md_tab_open <- "";
      md_tab_suffix <- " {.tabset}";
      md_tab_close <- "";
      do_md_tabs <- TRUE;
   } else if (do_md_tabs %in% c("Qmd", "qmd", "Quarto", "quarto")) {
      md_tab_open <- "::: {.panel-tabset}";
      md_tab_suffix <- "";
      md_tab_close <- ":::";
      do_md_tabs <- TRUE;
   }
   if (isTRUE(do_md_tabs)) {
      is_html_output <- FALSE;
      is_qmd <- FALSE;
      is_rmd <- FALSE;
      if (isTRUE(is_knitting)) {
         is_html_output <- knitr::is_html_output();
         # Returns NULL if not knitting, otherwise the source file path
         input_file <- knitr::current_input()
         if (length(input_file) > 0) {
            is_qmd <- grepl("[.]qmd$", ignore.case=TRUE, input_file)
            is_rmd <- grepl("[.]rmd$", ignore.case=TRUE, input_file)
         }
      }
      if (isTRUE(is_qmd) && isTRUE(is_html_output)) {
         # Quarto output
         if (length(md_tab_open) == 0) {
            md_tab_open <- "::: {.panel-tabset}";
         }
         if (length(md_tab_suffix) == 0) {
            md_tab_suffix <- "";
         }
         if (length(md_tab_close) == 0) {
            md_tab_close <- ":::";
         }
      } else if (isTRUE(is_rmd) && isTRUE(is_html_output)) {
         # Rmarkdown output
         if (length(md_tab_open) == 0) {
            md_tab_open <- "";
         }
         if (length(md_tab_suffix) == 0) {
            md_tab_suffix <- " {.tabset}";
         }
         if (length(md_tab_close) == 0) {
            md_tab_close <- "";
         }
      } else {
         if (length(md_tab_open) == 0) {
            md_tab_open <- "";
         }
         if (length(md_tab_suffix) == 0) {
            md_tab_suffix <- "";
         }
         if (length(md_tab_close) == 0) {
            md_tab_close <- "";
         }
      }
   }
   if (verbose) {
      jamba::printDebug("plot_mpf(): ",
         "do_md_tabs: ", do_md_tabs);
   }

   # Cat the markdown header, optional title, optional tabset open
   # 
   cat_md_header <- function
   (md_tab_level=2,
    heading="Header",
    md_tab_open=NULL,
    md_tab_suffix=NULL,
    ...)
   {
      md_header <- paste0(rep("#", md_tab_level), collapse="")
      if (length(md_tab_suffix) == 0) {
         md_tab_suffix <- "";
      }
      if (length(heading) > 0 && nchar(heading) > 0) {
         md_string <- paste0("\n\n",
            md_header, " ",
            heading, " ",
            md_tab_suffix, "\n\n");
         cat(md_string);
      }
      if (length(md_tab_open) > 0 && nchar(md_tab_open) > 0) {
         md_string <- paste0("\n\n", md_tab_open, "\n\n")
         cat(md_string);
      }
   }

   # Cat the markdown tab, optional tab suffix
   # 
   cat_md_tab <- function
   (md_tab_level=2,
    heading="Header",
    md_tab_suffix=NULL,
    ...)
   {
      md_header <- paste0(rep("#", md_tab_level), collapse="")
      if (length(md_tab_suffix) == 0) {
         md_tab_suffix <- "";
      } else {
         md_tab_suffix <- paste0(" ", md_tab_suffix);
      }
      if (length(heading) > 0 && nchar(heading) > 0) {
         md_string <- paste0("\n\n",
            md_header, " ",
            heading,
            md_tab_suffix, "\n\n");
         cat(md_string);
      }
   }

   # Cat the markdown tab close, optional
   # 
   cat_md_tab_close <- function
   (md_tab_close=NULL,
    ...)
   {
      if (length(md_tab_close) > 0 && nchar(md_tab_close) > 0) {
         md_string <- paste0("\n\n",
            md_tab_close,
            "\n\n");
         cat(md_string);
      }
   }

   # Optional tab header and tab open
   if (isTRUE(do_md_tabs)) {
      if (isFALSE(md_title)) {
         md_title <- "";
      }
      # first the title heading
      cat_md_header(md_tab_level=md_tab_level,
         heading=md_title,
         md_tab_open=md_tab_open,
         md_tab_suffix=md_tab_suffix)
      if (length(md_title) > 0 && nchar(md_title) > 0) {
         md_tab_level <- md_tab_level + 1;
      }
   }


   ###################################################
   # Iterate plot_which
   plot_names <- character();
   for (iplot in plot_which) {
      if (verbose) {
         jamba::printDebug("plot_mpf(): ",
            "Doing plot type: ", iplot);
      }
   
      if (isTRUE(do_md_tabs)) {
         cat_md_tab(md_tab_level=md_tab_level,
            heading=iplot,
            md_tab_suffix=md_tab_suffix)
      }

      #################################
      ## EnrichmentHeatmap
      if (grepl("EnrichmentHeatmap", iplot)) {
         plot_list$EnrichmentHeatmap <- do.call(
            EnrichmentHeatmap,
            c(
               alist(x=Mpf,
                  params=params$eh)))
         # plot_list$EnrichmentHeatmap <- EnrichmentHeatmap(Mpf, ...);
         if (isTRUE(do_newpage)) grid::grid.newpage();
         plot_names <- c(plot_names, "EnrichmentHeatmap")
      }

      #################################
      ## GenePathHeatmap
      if (grepl("GenePathHeatmap", iplot)) {
         plot_list$GenePathHeatmap <- do.call(GenePathHeatmap,
            c(
               alist(x=Mpf),
               params$gp))
         # plot_list$GenePathHeatmap <- GenePathHeatmap(Mpf, ...);
         if (isTRUE(do_newpage)) grid::grid.newpage();
         plot_names <- c(plot_names, "EnrichmentHeatmap")
      }

      #################################
      ## EnrichmentMap
      if (grepl("EnrichmentMap", iplot)) {
         # makeshift mem object?
         if (length(Mem) == 0) {
            Mem <- metadata(Mpf)$Mem;
         }
         if (!inherits(Mem, "Mem")) {
            # Decide what to do when one plot fails.
            stop_msg <- paste0("Input data does not contain Mem.")
            if (verbose) {
               jamba::printDebug("plot_mpf(): ",
                  "Skipping, Mem not available.",
                  fgText=c("darkorange2", "red"));
            }
            warning(stop_msg);
            next;
         }
         em_group <- intersect(em_group,
            c("community", "clusters", "cluster_labels"));
         if (length(em_group) == 0) {
            em_group <- "community";
         }
         if (length(em_group) > 1 && isTRUE(do_md_tabs)) {
            cat_md_header(md_tab_level=md_tab_level,
               heading=NULL,
               md_tab_open=md_tab_open,
               md_tab_suffix=md_tab_suffix)
            md_tab_level <- md_tab_level + 1;
         }
         for (icc in em_group) {
            if (verbose) {
               jamba::printDebug("plot_mpf(): ",
                  "Doing plot subtype: ", c(iplot, icc));
            }
            if (length(em_group) > 1 && isTRUE(do_md_tabs)) {
               cat_md_tab(md_tab_level=md_tab_level,
                  heading=icc,
                  md_tab_suffix=md_tab_suffix)
            }
            # em_repulse <- ifelse(
            #    is.numeric(params$em$repulse),
            #    params$em$repulse, 3.5);
            Emap <- do.call(EnrichmentMap,
               c(
                  alist(x=Mpf,
                     do_plot=TRUE,
                     params=params$em)))
            # Emap <- EnrichmentMap(Mpf,
            #    do_plot=TRUE,
            #    params=params$em,
            #    ...)
            if (isTRUE(do_newpage)) grid::grid.newpage();
            plot_names <- c(plot_names,
               paste0("EnrichmentMap", icc))
            if (length(plot_list$EnrichmentMap) == 0) {
               plot_list$EnrichmentMap <- list();
            }   
            plot_list$EnrichmentMap[[icc]] <- Emap;
         }
         if (length(cc_type) > 1 && isTRUE(do_md_tabs)) {
            md_tab_level <- md_tab_level - 1;
            cat_md_tab_close(md_tab_close=md_tab_close);
         }
         # Todo: figure out how to accept arguments
      }

      #################################
      ## CnetCollapsed
      if (grepl("CnetCollapsed", iplot)) {
         if (length(cc_type) == 0) {
            if (length(ClusterLabels(Mpf)) > 0) {
               cc_type <- "cluster";
            } else {
               cc_type <- "set";
            }
            if (igraph::vcount(Mpf@cnet_collapsed[[1]]) > 500) {
               cc_type <- paste0(cc_type, "2");
            }
         }
         if (length(cc_type) > 1 && isTRUE(do_md_tabs)) {
            cat_md_header(md_tab_level=md_tab_level,
               heading=NULL,
               md_tab_open=md_tab_open,
               md_tab_suffix=md_tab_suffix)
            md_tab_level <- md_tab_level + 1;
         }
         for (icc in cc_type) {
            if (length(cc_type) > 1 && isTRUE(do_md_tabs)) {
               cat_md_tab(md_tab_level=md_tab_level,
                  heading=icc,
                  md_tab_suffix=md_tab_suffix)
            }
            if (length(plot_list$CnetCollapsed) == 0) {
               plot_list$CnetCollapsed <- list();
            }
            plot_list$CnetCollapsed[[icc]] <- do.call(CnetCollapsed,
               c(
                  alist(x=Mpf,
                  type=icc,
                  params=params$cc)))
            # plot_list$CnetCollapsed <- CnetCollapsed(Mpf,
            #    type=icc,
            #    params=params$cc,
            #    ...);
            if (isTRUE(do_newpage)) grid::grid.newpage();
            plot_names <- c(plot_names, paste0("CnetCollapsed ", icc))
         }
         if (length(cc_type) > 1 && isTRUE(do_md_tabs)) {
            md_tab_level <- md_tab_level - 1;
            cat_md_tab_close(md_tab_close=md_tab_close);
         }
         if (isTRUE(do_newpage)) grid::grid.newpage();
      }

      #################################
      ## CnetExemplar
      if (grepl("CnetExemplar", iplot)) {
         if (length(ce_num) == 0) {
            ce_num <- 1;
         }
         if (length(ce_num) > 1 && isTRUE(do_md_tabs)) {
            cat_md_header(md_tab_level=md_tab_level,
               heading=NULL,
               md_tab_open=md_tab_open,
               md_tab_suffix=md_tab_suffix)
            md_tab_level <- md_tab_level + 1;
         }
         for (icc in ce_num) {
            if (length(ce_num) > 1 && isTRUE(do_md_tabs)) {
               cat_md_tab(md_tab_level=md_tab_level,
                  heading=icc,
                  md_tab_suffix=md_tab_suffix)
            }
            if (length(plot_list$CnetExemplar) == 0) {
               plot_list$CnetExemplar <- list();
            }
            plot_list$CnetExemplar[[icc]] <- do.call(CnetExemplar,
               c(
                  alist(x=Mpf,
                  num=icc,
                  params=params$ce)))
            # plot_list$CnetExemplar <- CnetExemplar(Mpf,
            #    num=icc,
            #    params=params$ce,
            #    ...);
            if (isTRUE(do_newpage)) grid::grid.newpage();
            plot_names <- c(plot_names, paste0("CnetExemplar ", icc))
         }
         if (length(ce_num) > 1 && isTRUE(do_md_tabs)) {
            md_tab_level <- md_tab_level - 1;
            cat_md_tab_close(md_tab_close=md_tab_close);
         }
         if (isTRUE(do_newpage)) grid::grid.newpage();
      }

      #################################
      ## CnetCluster
      if (grepl("CnetCluster", iplot)) {
         if (length(c_cluster) == 0) {
            c_cluster <- seq_along(Clusters(Mpf));
         }
         if (length(c_cluster) > 1 && isTRUE(do_md_tabs)) {
            cat_md_header(md_tab_level=md_tab_level,
               heading=NULL,
               md_tab_open=md_tab_open,
               md_tab_suffix=md_tab_suffix)
            md_tab_level <- md_tab_level + 1;
         }
         for (icc in c_cluster) {
            if (length(c_cluster) > 1 && isTRUE(do_md_tabs)) {
               cat_md_tab(md_tab_level=md_tab_level,
                  heading=icc,
                  md_tab_suffix=md_tab_suffix)
            }
            if (length(plot_list$CnetExemplar) == 0) {
               plot_list$CnetCluster <- list();
            }
            plot_list$CnetCluster[[icc]] <- do.call(CnetCluster,
               c(
                  alist(x=Mpf,
                  cluster=icc,
                  params=params[["c"]])))
            # plot_list$CnetCluster <- CnetCluster(Mpf,
            #    cluster=icc,
            #    params=params[["c"]],
            #    ...);
            if (isTRUE(do_newpage)) grid::grid.newpage();
            plot_names <- c(plot_names, paste0("CnetCluster ", icc))
         }
         if (length(c_cluster) > 1 && isTRUE(do_md_tabs)) {
            md_tab_level <- md_tab_level - 1;
            cat_md_tab_close(md_tab_close=md_tab_close);
         }
         if (isTRUE(do_newpage)) grid::grid.newpage();
      }
   }

   # Optionally close md tabs
   if (isTRUE(do_md_tabs)) {
      cat_md_tab_close(md_tab_close=md_tab_close);
      if (length(md_title) > 0 && nchar(md_title) > 0) {
         md_tab_level <- md_tab_level - 1;
      }
   }
   return(invisible(plot_list))
   return(invisible(plot_names))
}

#' Determine whether Mem should call newpage
#' 
#' Determine whether Mem should call newpage
#' 
#' ## Rules
#' 
#' * Running outside Positron, then yes. End.
#' * If running inside knitr, yes. End.
#' * Otherwise no.
#' 
#' @param x not currently implemented, it may be used
#'    in future.
#' @param ... additional arguments are ignored.
#' 
#' @keywords internal
mem_do_newpage <- function
(x,
 ...)
{
   #
   is_positron <- function() {
      identical(Sys.getenv("POSITRON"), "1")
   }
   if (isFALSE(is_positron)) {
      return(TRUE)
   }
   is_knitting <- isTRUE(getOption("knitr.in.progress"));
   if (isTRUE(is_knitting)) {
      return(TRUE)
   }
   return(FALSE)
}
