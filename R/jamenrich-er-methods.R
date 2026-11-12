# jamenrich-er-methods.R

#' Simple wrapper for pathway enrichment methods
#' 
#' Simple wrapper for pathway enrichment methods, to produce
#' enrichment results suitable for use with multienrichment
#' analysis.
#' 
#' ## Workflow
#' 
#' 1. Iterate each 'genedata'
#'    * Validate data format.
#'    * Optional: Update gene annotations using `genejam::freshenGenes()`.
#'    * Produce `data.frame` with numeric `hit` and `sign`.
#' 
#' 2. Iterate each 'pathwaydata'
#'    * Optional: Update gene annotations using `genejam::freshenGenes()`.
#'    * Produce `data.frame` TERM2GENE suitable for 'clusterProfiler'.
#' 
#' Note that gene annotations can be synchronized to ensure the
#' 'genedata' and 'pathwaydata' both use the same version of annotation.
#' 
#' 3. Iterate each 'enrich_method'
#'    * ORA - over-representation analysis, aka hypergeometric enrichment,
#'    Fisher's Exact Test.
#'    * GSEA - gene set enrichment analysis, using fast-GSEA (fgsea)
#' 
#' 4. Collate results
#'    * Convert to `enrichResult` or `gseaResult`.
#'    * Return all results as a nested `list`, indexed by:
#'       * 'genedata'
#'       * 'pathwaydata'
#'       * 'enrich_method'
#' 
#' ## Todo
#' 
#' * Add examples. Add example data.
#' * Add 'testthat' functional unit tests.
#' * Consider alternative enrichment methods:
#'    * `clusterProfiler` enrichment functions
#'    * `enrichR::enrichr()`
#' 
#' * Consider some custom "hooks":
#'    * Custom enrichment function with x,y,z input
#'    * Custom formatter to convert to `enrichResult` or
#'    `gseaResult`.
#' 
#' * Consider storing parameters as attributes:
#'    * 'ann_lib' and its package version for 'genejam'
#'    * `packageVersion('clusterProfiler')`, perhaps `sessionInfo()`
#'    * 'genedata_params'
#'    * 'pathwaydata_params'
#'    * `list` of 'genedata' `data.frame` used for enrichment
#' 
#' @returns nested `list` of `enrichResult` objects.
#'    * `list` by 'genedata', each with:
#'    * `list` by 'pathwaydata', each with
#'    * `list` by 'enrich_method', each with
#'    * `enrichResult` or `gseaResult`
#' 
#' @family jam enrich functions
#' 
#' @param genedata `list` named by the Enrichment condition,
#'    where each list element contains a suitable form of data.
#'    Each 'genedata' will be tested for enrichment independently.  
#' 
#'    Expected formats:
#'    * **`character`** vector of genes, which represent "gene hits".  
#'    Compatible with 'ORA' and cannot be used with 'GSEA'.
#'    * **`integer`** vector of integer signs limited to `1` or `-1`,
#'    named by gene.  
#'    Compatible with 'ORA' and cannot be used with 'GSEA'.  
#'    The sign applies directionality to genes, is used to 
#'    define 'direction' in the pathway, but not used for
#'    enrichment itself.
#'    * **`numeric`** with scalar values, named by gene.  
#'    Compatible with 'GSEA' and cannot be used with 'ORA'.  
#'    The `numeric` value ranks the genes, and
#'    is applied to directionality of the results.
#'    * **`data.frame`** of gene stats, with columns:
#'    'gene', 'pvalue', 'fold' in some conceptual form.  
#'    Compatible with 'ORA' and 'GSEA'.  
#' @param pathwaydata one of the following:
#'    * **Pathways**: `list` of pathways, named by pathway,
#'    whose values are individual genes assigned to each pathway.
#'    * **List of Pathway Sources**: `list` named by pathway source,
#'    each containing a `list` of pathways as described above,
#'    named by pathway and containing pathway genes.  
#'    Each Pathway Source will be analyzed one by one, independently.
#' @param enrich_method `character` string with enrichment method
#'    to perform:
#'    * **'ORA'**: uses `clusterProfiler::enricher()` for a standard
#'    hypergeometric (Fisher's Exact) test of enrichment.
#'    * **'GSEA'**: uses `clusterProfiler::GSEA()` for ranked gene set
#'    enrichment analysis. It calls 'fgsea' for fast-GSEA.
#' @param universe one of the following:
#'    * **`character`** vector of genes to define the background
#'    universe for the enrichment test.  
#'    The same universe is used for all enrichment tests.  
#'    When NULL as by default, it uses all genes defined in 'pathwaydata'.
#'    * `list` named by `names(genedata)` to define a universe
#'    for each individual entry in 'genedata'.
#'    Todo: Consider unique universe by genedata/pathwaydata.
#' @param genedata_params `list` to manage 'genedata'.
#'    * **refresh_genes**: `logical` whether to update gene symbols
#'    using `genejam::freshenGenes()`
#'    * **ann_lib**: `character` with Bioconductor library to use
#'    when 'refresh_genes' is TRUE.
#'    * **handle_multiple**: `character` string in `genejam::freshenGenes()`,
#'    default "best_each" uses the best match for each gene symbol
#'    per row.
#'    * **empty_rule**: `character` how to handle empty gene 'SYMBOL'.
#'    Default 'empty' leaves the value blank.
#'    Other options are 'original' to keep the original value, or
#'    'na' to insert `NA`.
#'    * **pvalue_cutoff**,**fold_cutoff**: `numeric` thresholds used when
#'    `genedata` is supplied as a `data.frame`. It will define gene
#'    statistical hits using these thresholds.  
#'    Default NULL will not apply thresholds, and will not define hits.
#'    * **gene_colnames**,**pvalue_colnames**,**fold_colnames**:
#'    `character` vectors with patterns used to match column names
#'    when 'genedata' contains a `data.frame`.
#'    These columns are required.  
#'    The 'gene_colname' is used as the initial gene, and will be
#'    used when `refresh_genes=TRUE`.
#'    * **tstat_colnames**: `character` (optional) pattern for
#'    t-statistic column name when 'genedata' contains a `data.frame`.  
#'    It is only used when `rank_strategy='t-statistic'`.
#'    * **rank_strategy**: `character` string used when `genedata`
#'    contains a `data.frame`:
#'       * **'signed-P'** (default) uses the `-log10(Pvalue)` multiplied by
#'       the sign of the fold change. It produces the same rank order as
#'       using the t-statistic, useful when the `'t'` column is not available.
#'       * **'t-statistic'** uses the t-statistical column.
#'       * **'P-value'** uses the `-log10(Pvalue)` to rank entries.
#'       * **'log2fold'** uses the fold column to rank entries. Although a
#'       common default, we recommend the signed-P for improved resolution.
#' 
#'    * **try_list**: passed to `genejam::freshenGenes()` to determine which
#'    annotation methods to search when resolving to gene SYMBOL.
#'    They are handled in priority order, there 'best_each' will return the
#'    first (best) result for each input gene.
#'    When input genes are entirely integer values, they are assumed
#'    to be ENTREZID values, which skips this step.
#' @param pathwaydata_params `list` to manage 'pathwaydata'.
#'    * See 'genedata_params' for description of most fields.
#'    * pathway_colnames: `character` vector with patterns for preferred
#'    pathway column names, the first match with `find_colname(pattern, x)`
#'    is used as the pathway name.
#'    For the purpose of this process, only the pathway name is retained,
#'    and all other annotations should be retrieved separately as needed.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to internal functions.
#' 
#' 
#' @examples
#' # To be defined
#' # genedata is data.frame with gene stat results, or vector of gene symbols
#' # pathwaydata is a list of gene symbol vectors, or data.frame with TERM2GENE
#' # or data.frame with gene column with delimited gene symbols.
#' 
#' @export
jam_enricher <- function
(
   genedata,
   pathwaydata,
   enrich_method=c(
      "ORA",
      "GSEA"
   ),
   universe=NULL,
   genedata_params=list(
      ann_lib="org.Hs.eg.db"
   ),
   pathwaydata_params=list(
      ann_lib="org.Hs.eg.db"
   ),
   verbose=FALSE,
   ...
) {
   # Workflow:
   #
   # Iterate each genedata element
   # - determine datatype, compatible enrich methods
   # - apply optional refresh gene step
   #
   # Iterate each pathwaydata element
   # - determine format
   # - apply optional refresh gene step
   #
   # Iterate each enrich_method
   # - apply each compatible enrich_method
   # - convert output to `enrichResult` or `gseaResult`
   
   # update default params with user-supplied values
   #
   enrich_method <- match.arg(enrich_method, several.ok=TRUE)
   enrich_method <- jamba::nameVector(unique(enrich_method));
   # genedata_params
   default_genedata_params <- eval(formals(
      process_genedata)$genedata_params)
   genedata_params <- modifyList(
      default_genedata_params,
      genedata_params,
      keep.null=FALSE
   )
   # pathwaydata_params
   default_pathwaydata_params <- eval(formals(
      process_pathwaydata)$pathwaydata_params)
   pathwaydata_params <- modifyList(
      default_pathwaydata_params,
      pathwaydata_params,
      keep.null=FALSE
   )

   # genedata
   if (!is.list(genedata) || length(names(genedata)) == 0) {
      cli::cli_abort(paste0(
         "{.var genedata} must be a {.cls list} with {.code names}."
      ))
   }

   # process pathwaydata
   if (!is.list(pathwaydata) ||
      (is.list(pathwaydata) && !is.list(pathwaydata[[1]]))
   ) {
      pathwaydata <- list(pathwaydata=pathwaydata)
   }
   pathwaydata <- process_pathwaydata(
      pathwaydata=pathwaydata,
      pathwaydata_params=pathwaydata_params,
      ...
   )
   # names
   gdnames <- jamba::nameVectorN(genedata)
   pwnames <- jamba::nameVectorN(pathwaydata)
   # iterate genedata
   gd_results <- lapply(gdnames, function(gdname){
      igd <- genedata[[gdname]]
      # each genedata as igdf
      igdf <- process_genedata(
         genedata=igd,
         genedata_params=genedata_params,
         ...
      )
      if (verbose) {
         jamba::printDebug(
            "jam_enricher(): ",
            "processing genedata '", gdname, "'"
         )
      }

      # iterate pathwaydata
      pd_results <- lapply(pwnames, function(pwname){
         ipdf <- pathwaydata[[pwname]]
         # each pathwaydata as ipdf
         use_em <- enrich_method;
         if (!any(igdf$hit)) {
            use_em <- jamba::unvigrep("ORA", use_em);
         }
         if (all(igdf$sign %in% c(-1, 1, 0, NA))) {
            use_em <- jamba::unvigrep("gsea", use_em);
         }
         em_results <- lapply(use_em, function(iem){
            # iterate enrichment as iem
            #
            if ("ORA" %in% iem) {
               # ORA
               if (verbose) {
                  jamba::printDebug(
                     "jam_enricher(): ",
                     sep="",
                     c("ORA: ", "clusterProfiler::enricher()")
                  )
               }
               use_igdf <- subset(igdf, hit & !is.na(sign) & sign != 0)
               use_gene <- use_igdf$gene;
               er <- clusterProfiler::enricher(
                  gene=use_gene,
                  pvalueCutoff=1,
                  qvalueCutoff=1,
                  minGSSize=10,
                  maxGSSize=10000,
                  TERM2GENE=ipdf
               )
               # add pathwaydata name as 'gs_cat'
               er@result$gs_cat <- pwname;
               if (any(jamba::rmNA(igdf$sign) < 0)) {
                  if (verbose) {
                     jamba::printDebug(
                        "jam_enricher(): ",
                        sep="",
                        c("ORA: ", "add_pathway_direction()")
                     )
                  }
                  use_gene_hits <- jamba::nameVector(
                     sign(use_igdf$sign),
                     use_igdf$gene
                  )
                  er <- add_pathway_direction(
                     x=er,
                     gene_hits=use_gene_hits,
                     delim="/",
                     dir_colname="direction"
                  )
               }
               er
            } else if ("GSEA" %in% iem) {
               # GSEA
               if (verbose) {
                  jamba::printDebug(
                     "jam_enricher(): ",
                     sep="",
                     c("GSEA: ", "clusterProfiler::GSEA()")
                  )
               }
               use_genes_ordered <- jamba::nameVector(
                  jamba::mixedSortDF(igdf, "-sign")[, c("sign", "gene")]
               )
               gsear <- clusterProfiler::GSEA(
                  geneList=use_genes_ordered,
                  TERM2GENE=ipdf,
                  verbose=FALSE,
                  minGSSize=10,
                  maxGSSize=10000,
                  pvalueCutoff=1,
                  eps=1e-50
               )
               attr(gsear, "genes_ordered") <- use_genes_ordered;
               gsear
            } else {
               cli::cli_abort("fgsea is not yet implemented.")
            }
         })
         em_results;
      })
      pd_results
   })
   return(gd_results);
}

#' Detect gene data format for pathway enrichment
#' 
#' @returns `list` with
#'    * 'found_type'
#'    * 'compatible_method' which may contain
#'       * 'ORA' for over-representation analysis (not ranked), or
#'       * 'ranked' for ranked analysis (including GSEA, fgsea).
#' 
#' @inheritParams jam_enricher
#' 
#' @keywords internal
#' @noRd
detect_genedata_format <- function(
   genedata,
   ...
) {
   #
   # iterate each element of list input
   if (is.list(genedata)) {
      typelist <- lapply(genedata, function(igenedata){
         detect_genedata_format(
            genedata=igenedata,
            ...
         )
      })
      return(typelist);
   }

   # process individual genedata input
   found_type <- "";
   if (inherits(genedata, "character")) {
      found_type <- "hits";
   } else if (inherits(genedata, "numeric")) {
      if (all(genedata %in% c(-1, 1, 0, NA))) {
         found_type <- "signedhits";
      } else {
         found_type <- "scores";
      }
      if (length(names(genedata)) == 0) {
         # it must have names
         cli::cli_abort(paste0(
            "{.var genedata} with {.cls numeric} values ",
            "must have {.code names(genedata)}."
         ))
      }
   } else if (inherits(genedata, c("data.frame", "tbl"))) {
      # validate the colnames here?
      found_type <- "stats";
   }
   # compatible methods
   compatible_methods <- list(
      hits=c("ORA"),
      signedhits=c("ORA"),
      scores=c("ranked"),
      stats=c(
         "ORA",
         "ranked"
      )
   )
   compatible_method <- compatible_methods[[found_type]];
   retlist <- list(
      found_type=found_type,
      compatible_method=compatible_method
   )
   return(retlist)
}


#' Process pathway data for enrichment
#' 
#' @returns `data.frame` with two columns:
#'    * 'TERM': the pathway name
#'    * 'GENE': one gene per row, after updating genes when
#'    `params$refresh_genes=TRUE`.
#' 
#' @keywords internal
#' @noRd
process_pathwaydata <- function(
   pathwaydata,
   pathwaydata_params=list(
      refresh_genes=TRUE,
      ann_lib="org.Hs.eg.db",
      handle_multiple="best_each",
      empty_rule="empty",
      gene_colnames=c(
         "gene_symbol",
         "db_gene_symbol",
         "gene",
         "ncbi_gene"
      ),
      pathway_colnames=c(
         "gs_name",
         "gs_description",
         "name",
         "ID",
         "Pathway",
         "GO[:]*(BP|MF|CC)"
      ),
      try_list=c(
         "SYMBOL2EG", "ACCNUM2EG",
         "ENSEMBL2EG", "ALIAS2EG"
      ),
      sep="/"
   ),
   ...
) {
   # Workflow
   # - refresh gene nomenclature if necessary
   
   # update default params with user-supplied values
   default_pathwaydata_params <- eval(formals(
      process_pathwaydata)$pathwaydata_params)
   params <- modifyList(
      default_pathwaydata_params,
      pathwaydata_params,
      keep.null=FALSE
   )
   
   # iterate list of pathwaydata
   if (is.list(pathwaydata) && is.list(pathwaydata[[1]])) {
      pathwaydata_list <- lapply(pathwaydata, function(ipd){
         process_pathwaydata(
            pathwaydata=ipd,
            pathwaydata_params=pathwaydata_params,
            ...
         )
      })
      return(pathwaydata_list);
   }

   # process individual pathwaydata
   if (inherits(pathwaydata, c("data.frame", "tbl"))) {
      # convert data.frame,tbl to 'list'
      # then process the same as 'list'
      use_gene_colname <- find_colname(
         params$gene_colnames,
         pathwaydata
      )
      # use_gene_colnames <- head(intersect(
      #    tolower(params$gene_colnames),
      #    tolower(colnames(pathwaydata))
      # ), 1)
      # use_gene_colname <- colnames(pathwaydata)[match(
      #    tolower(use_gene_colnames),
      #    tolower(colnames(pathwaydata))
      # )];
      pathwaygenes <- pathwaydata[[use_gene_colname]];
      if (is.list(pathwaygenes)) {
         if (inherits(pathwaygenes, "AsIs")) {
            pathwaygenes <- unclass(pathwaygenes)
         }
      } else {
         pathwaygenes <- strsplit(pathwaygenes, "[;,/ ]+")
      }
      pathway_colname <- find_colname(
         params$pathway_colnames,
         pathwaydata
      )
      # Note: Error will occur if names are duplicated, this is expected
      # Todo: Check first and write useful error message.
      names(pathwaygenes) <- pathwaydata[[pathway_colname]];
      pathwaydata <- pathwaygenes;
      rm(pathwaygenes);
   }
   if (!is.list(pathwaydata)) {
      cli::cli_abort(paste0(
         "Input must be {.cls data.frame} or {.cls list}, and ",
         "could not be processed as intended."
      ))
   }
   TERM2GENE <- data.frame(
      TERM=character(0),
      GENE=character(0)
   )
   if (isTRUE(params$refresh_genes)) {
      # assume list of character vectors
      # - initial strategy: concatenate values per row
      pwdf <- data.frame(
         check.names=FALSE,
         pw_gene=jamba::cPasteSU(pathwaydata, sep=";")
      )
      rownames(pwdf) <- names(pathwaydata)
      if (all(grepl("^[0-9;,/]*$", pwdf$pw_gene))) {
         # ENTREZID format
         pwdf <- jamba::renameColumn(
            pwdf,
            from="pw_gene",
            to="ENTREZID"
         )
      }
      pwdf1 <- genejam::freshenGenes(
         pwdf,
         handle_multiple=params$handle_multiple,
         empty_rule=params$empty_rule,
         ann_lib=params$ann_lib,
         try_list=params$try_list,
         sep=params$sep
      )
      rownames(pwdf1) <- rownames(pwdf)
      TERM2GENEW <- data.frame(
         check.names=FALSE,
         TERM=rownames(pwdf),
         GENE=pwdf1$SYMBOL
      )
      # new_pathwaydata <- strsplit(pwdf1$SYMBOL, ";")
      # names(new_pathwaydata) <- rownames(pwdf)
      # ln0 <- which(lengths(new_pathwaydata) == 0)
      # new_pathwaydata[ln0] <- ""
   } else {
      TERM2GENEW <- data.frame(
         check.names=FALSE,
         TERM=names(pathwaydata),
         GENE=jamba::cPasteSU(
            pathwaydata,
            sep=params$sep
         )
      )
   }
   TERM2GENE <- deconcat_df2(
      TERM2GENEW,
      column="GENE",
      split=params$sep
   )
   # subset for
   # rows with GENE values
   # unique TERM:GENE entries
   TERM2GENE <- unique(
      subset(TERM2GENE, !GENE %in% c(NA, ""))
   )
   return(TERM2GENE);
}

#' Process gene data for enrichment
#' 
#' Process gene data for enrichment
#' 
#' Workflow:detect_genedata_format
#' * Define the gene column
#' * Define the sign, use "1" for no sign
#'    * Calculate sign for "signed-P"
#'    * Use 't' for "t-statistic"
#' 
#' * Refresh gene symbols to 'SYMBOL'
#' * Remove entries with no 'SYMBOL'
#' * Deconcatenate / Expand 'SYMBOL' with multiple genes
#' * Resolve duplicate 'SYMBOL' to one row
#'    * Sort by absolute value of sign descending, then sign descending.
#'    * Retain the first observation per 'SYMBOL'
#' 
#' ## Params recognized
#' 
#' * ann_lib: passed to `genejam::freshenGenes()`
#' * pvalue_cutoff, fold_cutoff: used to define gene hits when
#' 'genedata' is `data.frame`. When not supplied, data can only
#' be used for ranked enrichment.
#' 
#' ## Todo
#' 
#' * Provide some kind of summary for gene resolution
#'    * Total gene values
#'    * Number of distinct input genes
#'    * Percent matched, percent retained
#'    * Number of multi-gene entries
#'    * Number of duplicated genes
#'    * Number of distinct output genes
#' 
#' @returns `data.frame` with columns:
#'    * 'gene': the gene symbol
#'       * `genedata_params$refresh_genes=FALSE` (default) uses
#'       existing gene symbol from 'gene_colnames'
#'       * `genedata_params$refresh_genes=FALSE` uses 'SYMBOL'
#'       after applying `genejam::freshenGenes()`.
#'       **Make sure to use correct 'ann_lib'.**
#' 
#'    * 'sign': the `numeric` sign or scalar value:
#'       * entirely '1' for non-signed gene hits,
#'       cannot be used for ranked enrichment.
#'       * '1' or '-1' for signed gene hits,
#'       cannot be used for ranked enrichment.
#'       * `numeric` scalar values used for ranked enrichment.
#' 
#'    * 'hit': `logical` indicating whether each entry met
#'    the statistical cutoffs `pvalue_cutoff` and `fold_cutoff`
#'    when defined.
#'    When all entries are FALSE there were no cutoffs applied,
#'    and data should only be used for ranked enrichment.
#' 
#'    * The names are assigned using gene symbols, either as supplied
#'    or after resolution when `genedata_params$refresh_genes=TRUE`.
#'    * Finally, attribute 'compatible_method' indicates which
#'    enrichment tests are suitable.
#' 
#' @inheritParams jam_enricher
#' 
#' @keywords internal
#' @noRd
process_genedata <- function(
   genedata,
   genedata_params=list(
      refresh_genes=TRUE,
      ann_lib="org.Hs.eg.db",
      handle_multiple="best_each",
      empty_rule="empty",
      fold_cutoff=NULL,
      pvalue_cutoff=NULL,
      gene_colnames=c(
         "SYMBOL",
         "probes",
         "gene",
         "gene_symbol",
         "geneid",
         "gene_id",
         "db_gene_symbol",
         "ncbi_gene"
      ),
      fold_colnames=c(
         "fold", "logFC", "lfc",
         "log2FC", "log2fold"
      ),
      pvalue_colnames=c(
         "adj.p.val", "adj.*P.*val", "FDR",
         "q.val|qval",
         "adjp", "padjust", "p.adjust",
         "p.value", "pval"),
      tstat_colnames=c("t", "t.stat"),
      rank_strategy=c(
         "signed-P",
         "t-statistic",
         "P-value",
         "log2fold"
      ),
      try_list=c("SYMBOL2EG", "ACCNUM2EG", "ENSEMBL2EG", "ALIAS2EG"),
      sep="/"
   ),
   verbose=FALSE,
   ...
) {
   #
   # update default params with user-supplied values
   default_genedata_params <- eval(formals(
      process_genedata)$genedata_params)
   params <- modifyList(
      default_genedata_params,
      genedata_params,
      keep.null=FALSE
   )

   # process each type of genedata as relevant
   use_genes <- character();
   use_gene_hits <- character();
   use_signs <- numeric();
   if (inherits(genedata, c("character", "factor"))) {
      if (verbose) {
         jamba::printDebug(
            "process_genedata(): ",
            "detected gene hits"
         )
      }
      use_genes <- unique(as.character(genedata));
      use_signs <- rep_len(1, length(use_genes));
      use_gene_hits <- rep_len(TRUE, length(use_genes))
   } else if (inherits(genedata, "numeric")) {
      if (length(names(genedata)) == 0) {
         cli::cli_abort(paste0(
            "{.var genedata} with {.cls numeric} values must also ",
            "have {.code names(genedata)} with gene values."
         ))
      }
      use_genes <- names(genedata);
      if (all(genedata %in% c(1, -1, 0, NA))) {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "detected directional gene hits"
            )
         }
         # directional sign
         #
         # Todo: confirm this step
         # replace NA with 0
         use_signs <- jamba::rmNA(
            genedata,
            naValue=0
         )
         # remove genes with 0
         non_zero <- (!use_signs %in% 0);
         use_signs <- genedata[non_zero];
         use_genes <- use_genes[non_zero];
      } else {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "detected scalar gene hits"
            )
         }
         # directional scalar
         # remove entries with NA, as "not tested"
         not_na <- (!genedata %in% NA)
         use_signs <- genedata[not_na];
         use_genes <- use_genes[not_na];
      }
      use_gene_hits <- rep_len(TRUE, length(use_genes))
   } else if (inherits(genedata, c("data.frame", "tbl"))) {
      # data.frame
      if (verbose) {
         jamba::printDebug(
            "process_genedata(): ",
            "detected gene stat table"
         )
      }
      rank_strategy <- match.arg(
         params$rank_strategy,
         choices=eval(
            formals(process_genedata)$genedata_params$rank_strategy
         )
      )
      fold_colname <- find_colname(params$fold_colnames, genedata)
      pvalue_colname <- find_colname(params$pvalue_colnames, genedata)
      gene_colname <- find_colname(params$gene_colnames, genedata)
      tstat_colname <- find_colname(params$tstat_colnames, genedata)
      # validate colnames
      if ("t-statistic" %in% rank_strategy && length(tstat_colname) != 1) {
         cli::cli_abort(c(
            "{.var rank_strategy} was {.val {rank_strategy}} but ",
            "{.var tstat_colname} was not found."
         ))
      }
      if (!all(
         length(fold_colname) == 1 &&
         length(pvalue_colname) == 1 &&
         length(gene_colname) == 1
      )) {
         cli::cli_abort(c(
            "Not all columns were defined: ",
            "{.var fold_colname} = {.val {fold_colname}},",
            "{.var pvalue_colname} = {.val {pvalue_colname}},",
            "{.var gene_colname} = {.val {gene_colname}}"
         ))
      }
      if (verbose) {
         jamba::printDebug(
            "process_genedata():\n     ",
            sep="",
            c("\n     pvalue_colname:", pvalue_colname),
            c("\n     fold_colname:", fold_colname),
            c("\n     gene_colname:", gene_colname)
         )
      }

      # Define hits
      pvalue_cutoff <- params$pvalue_cutoff;
      fold_cutoff <- params$fold_cutoff;
      if (length(pvalue_cutoff) == 1 && is.numeric(pvalue_cutoff)) {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "apply pvalue_cutoff:",
               pvalue_cutoff
            )
         }
         if (length(use_gene_hits) == 0) {
            use_gene_hits <- rep_len(TRUE, nrow(genedata))
         }
         use_gene_hits <- (
            use_gene_hits &
            !is.na(genedata[[pvalue_colname]]) &
            genedata[[pvalue_colname]] <= pvalue_cutoff
         );
      }
      if (length(fold_cutoff) == 1 && is.numeric(fold_cutoff)) {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "apply fold_cutoff:",
               fold_cutoff
            )
         }
         if (length(use_gene_hits) == 0) {
            use_gene_hits <- rep_len(TRUE, nrow(genedata))
         }
         use_gene_hits <- (
            use_gene_hits &
            !is.na(genedata[[fold_colname]]) &
            abs(genedata[[fold_colname]]) >= fold_cutoff
         );
      }
      if (length(use_gene_hits) == 0) {
         use_gene_hits <- rep_len(FALSE, nrow(genedata));
      } else {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "defined gene stat hits"
            )
         }
      }

      # process rank
      use_genes <- genedata[[gene_colname]];
      if ("t-statistic" %in% rank_strategy) {
         use_signs <- genedata[[tstat_colname]];
      } else if ("signed-P" %in% rank_strategy) {
         # Todo: check and convert to log2fold if needed
         use_signs <- (
            -log10(genedata[[pvalue_colname]]) *
               sign(genedata[[fold_colname]])
         );
      } else if ("P-value" %in% rank_strategy) {
         use_signs <- -log10(genedata[[pvalue_colname]]);
      } else if ("log2fold" %in% rank_strategy) {
         # Todo: check and convert to log2fold if needed
         use_signs <- (genedata[[fold_colname]]);
      } else {
         cli::cli_abort(c(
            "{.var rank_strategy} = {.val {rank_strategy}} is not supported."
         ))
      }
      # Re-order by sign, descending
      new_order <- order(-use_signs);
      use_signs <- use_signs[new_order];
      use_genes <- use_genes[new_order];
      use_gene_hits <- use_gene_hits[new_order];
   } else {
      cli::cli_abort(c(
         "{.var genedata} as {.cls {class(genedata)}} is not supported."
      ))
   }

   # Process gene annotation
   if (isTRUE(params$refresh_genes)) {
      if (verbose) {
         jamba::printDebug(
            "process_genedata(): ",
            "refresh gene SYMBOL"
         )
      }
      # assume list of character vectors
      # - initial strategy: concatenate values per row
      gddf <- data.frame(
         check.names=FALSE,
         gd_gene=use_genes
      )
      if (all(grepl("^[0-9;,/]*$", gddf$gd_gene))) {
         # ENTREZID format
         gddf <- jamba::renameColumn(
            gddf,
            from="gd_gene",
            to="ENTREZID"
         )
      }
      # freshen gene SYMBOL
      gddf1 <- genejam::freshenGenes(
         gddf,
         handle_multiple=params$handle_multiple,
         empty_rule=params$empty_rule,
         ann_lib=params$ann_lib,
         try_list=params$try_list,
         sep=params$sep
      )
      # Add the sign
      gddf1$use_signs <- use_signs;
      gddf1$abs_use_signs <- abs(use_signs);
      gddf1$use_gene_hits <- use_gene_hits;
   } else {
      # Use gene as-is
      gddf1 <- data.frame(
         check.names=FALSE,
         SYMBOL=use_genes,
         use_signs=use_signs,
         abs_use_signs=abs(use_signs),
         use_gene_hits=use_gene_hits
      )

   }

   # remove entries with no SYMBOL
   gddf1 <- subset(gddf1, !SYMBOL %in% c("", NA))

   # handle multiple gene per row
   if (any(grepl(params$sep, gddf1$SYMBOL))) {
      if (verbose) {
         jamba::printDebug(
            "process_genedata(): ",
            "deconcatenate values with multiple gene"
         )
      }
      gddf1 <- deconcat_df2(
         gddf1,
         column="SYMBOL",
         split=params$sep
      )
   }

   # resolve duplicate SYMBOL
   if (anyDuplicated(gddf1$SYMBOL)) {
      if (verbose) {
         jamba::printDebug(
            "process_genedata(): ",
            "resolve genes with multiple values",
            head(which(duplicated(gddf1$SYMBOL)))
         )
         dv <- unique(gddf1$SYMBOL[duplicated(gddf1$SYMBOL)]);
         print(head(subset(gddf1, SYMBOL %in% dv)));
      }
      # verify sort order:
      # - first priority: stat hits
      #   (requires pvalue_cutoff and fold_cutoff)
      # - then highest magnitude sign
      # - then first sign in descending order
      #
      # There is potential to remove highest magnitude sign when
      # it is not a stat hit, and when a lower sign does meet
      # statistical significance.
      # What should be done? Note when it happens?
      gddf1 <- jamba::mixedSort(
         gddf1,
         byCols=c(
            "-use_gene_hits",
            "-abs_use_signs",
            "-use_signs"
         )
      )
      # Retain only the first of each SYMBOL
      gddf1 <- subset(gddf1, !duplicated(SYMBOL));
   }

   # Prepare named numeric vector
   use_signs <- gddf1$use_signs;
   use_gene_hits <- gddf1$use_gene_hits;
   use_genes <- gddf1$SYMBOL;
   # use_genedata <- jamba::nameVector(
   #    use_signs,
   #    use_genes
   # )
   use_genedata_df <- data.frame(
      check.names=FALSE,
      gene=use_genes,
      sign=use_signs,
      hit=use_gene_hits
   )

   # attribute 'compatible_method'
   # genetypelist <- detect_genedata_format(genedata=use_genedata)
   # compatible_method <- genetypelist$compatible_method;
   # attr(use_genedata, "compatible_method") <- compatible_method;
   return(use_genedata_df);
}
