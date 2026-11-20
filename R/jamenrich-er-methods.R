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
#' * Return genedata after processing so it can be reviewed.
#' * Consider preparing hit incidence matrices, or lists.
#'    * ORA would use genedata
#'    * GSEA would use enrichment genes with input sign.
#' 
#' * Add example data and examples.
#' * Add 'testthat' functional unit tests.
#' * Consider multiple cutoffs, as a named list, only for ORA.
#'    * cutoffs=list(
#'         p01_fc1.5=c(pvalue_cutoff=0.01, fold_cutoff=1.5),
#'         p05_fc1.5=c(pvalue_cutoff=0.05, fold_cutoff=1.5)
#'      )
#' 
#' * Consider flattening the nested list output.
#' * Consider alternative enrichment methods:
#'    * 'clusterProfiler' enrichment functions
#'    * 'enrichR::enrichr()'
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
#'    * `list` named by 'genedata', each with:
#'    * `list` named by 'pathwaydata', each with
#'    * `list` named by 'enrich_method', each with
#'    * 'ORA'=`enrichResult` or 'GSEA'=`gseaResult`
#'    
#'    In addition, an attribute 'pathwaydata' is added with the
#'    updated `pathwaydata` after `process_pathwaydata()`
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
      verbose=verbose > 1,
      ...
   )
   # names
   gdnames <- jamba::nameVectorN(genedata)
   pwnames <- jamba::nameVectorN(pathwaydata)
   # iterate genedata
   gd_results <- lapply(gdnames, function(gdname){
      igd <- genedata[[gdname]]
      # each genedata as igdf
      if (verbose) {
         jamba::printDebug(
            "jam_enricher(): ",
            "processing genedata '", gdname, "'"
         )
      }
      igdf <- process_genedata(
         genedata=igd,
         genedata_params=genedata_params,
         verbose=verbose > 1,
         ...
      )
      if (length(universe) == 0) {
         universe <- unique(igdf[[1]]);
         if (verbose) {
            jamba::printDebug(
               "jam_enricher(): ",
               "Defined universe size: ",
               jamba::formatInt(length(universe))
            )
         }
      }
      if (verbose > 1) {
         jamba::printDebug(
            "jam_enricher(): ",
            "middle(igdf):"
         )
         print(jamba::middle(igdf, 20));
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
                     c(
                        "ORA: ",
                        "clusterProfiler::enricher()")
                  )
               }
               use_igdf <- subset(igdf, hit & !is.na(sign) & sign != 0)
               use_gene <- use_igdf$gene;
               er <- clusterProfiler::enricher(
                  gene=use_gene,
                  pvalueCutoff=1,
                  qvalueCutoff=1,
                  universe=universe,
                  minGSSize=10,
                  maxGSSize=10000,
                  TERM2GENE=ipdf
               )
               # add pathwaydata name as 'gs_cat'
               er@result$gs_cat <- pwname;
               use_gene_hits <- jamba::nameVector(
                  sign(use_igdf$sign),
                  use_igdf$gene
               )
               if (any(jamba::rmNA(igdf$sign) < 0)) {
                  if (verbose) {
                     jamba::printDebug(
                        "jam_enricher(): ",
                        sep="",
                        c(
                           "ORA: ",
                           "add_pathway_direction()")
                     )
                  }
                  er <- add_pathway_direction(
                     x=er,
                     gene_hits=use_gene_hits,
                     delim="/",
                     dir_colname="direction"
                  )
               }
               attr(er, "genes_signed") <- use_gene_hits;
               er
            } else if ("GSEA" %in% iem) {
               # GSEA
               if (verbose) {
                  jamba::printDebug(
                     "jam_enricher(): ",
                     sep="",
                     c(
                        "GSEA: ",
                        "clusterProfiler::GSEA()")
                  )
               }
               if (verbose > 1) {
                  jamba::printDebug("jam_enricher(): ",
                     "middle(igdf):");
                  jamba::middle(igdf, 20);
               }
               use_genes_ordered <- jamba::nameVector(
                  jamba::mixedSortDF(igdf, "-sign")[, c("sign", "gene")]
               )
               if (verbose > 1) {
                  jamba::printDebug("jam_enricher(): ",
                     "middle(use_genes_ordered):");
                  jamba::middle(use_genes_ordered, 20);
               }
               gsear <- clusterProfiler::GSEA(
                  geneList=use_genes_ordered,
                  TERM2GENE=ipdf,
                  minGSSize=10,
                  maxGSSize=10000,
                  pvalueCutoff=1,
                  eps=1e-50,
                  verbose=verbose > 1
               )
               # add pathwaydata name as 'gs_cat'
               gsear@result$gs_cat <- pwname;
               attr(gsear, "genes_ordered") <- use_genes_ordered;
               gsear
            } else {
               cli::cli_abort("No other enrichment methods yet implemented.")
            }
         })
         em_results;
      })
      # add genedata output from process_genedata()
      attr(pd_results, "genedata") <- igdf;
      pd_results
   })

   # add pathwaydata as attribute for review
   attr(gd_results, "pathwaydata") <- pathwaydata;
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


#' Process pathway data for enrichment using 'jam_enricher()'
#' 
#' Process pathway data for enrichment using 'jam_enricher()'
#' 
#' This function recognizes several types of pathway data as input,
#' and returns a consistent format as output: a `list` of pathway
#' sources, each containing a `data.frame` with columns
#' 'TERM' and 'GENE'.
#' 
#' ## Custom parameters
#' 
#' Argument `pathwaydata_params` controls additional options:
#' 
#' * `refresh_genes=TRUE` enables `genejam::freshenGenes()` which
#' uses Bioconductor annotations to update the current gene annotation
#' to the current version available (to the active R session).  
#' Several arguments are passed to this function:
#'    * `ann_lib`: the core annotation library,
#'    'org.Hs.eg.db' is default for Human; 'org.Mm.eg.db' is used
#'    for Mouse. This annotation database must be available,
#'    for example install Rat data with
#'    `BiocManager::install('org.Rn.eg.db')`.
#'    Alternatively, `ann_lib` can be a `list` with one or more
#'    `AnnotationDbi` or `OrgDb` objects.
#'    * `handle_multiple='best_each'`: converts each gene to the best
#'    available match, using the priority order defined in `try_list`.
#'    * `ignore.case=TRUE`: Note that Human symbols are UPPERCASE,
#'    so for a moderate speed boost use `ignore.case=FALSE`.
#'    Note that Mouse by convention uses mixed case, for example "Apoe",
#'    and so for mouse and most other species, use `ignore.case=TRUE`.
#'    * `empty_rule='empty'`: converts any unmatched entry to '' or NA
#'    so it is not retained in the output. This decision prevents keeping
#'    a gene or accession which can never be recognized as an official
#'    gene during pathway enrichment analysis.  
#'    Use `empty_rule='original'` to retain the original entry whenever
#'    the entry cannot be matched to any other annotation in `try_list`.
#'    * `try_list=c('SYMBOL2EG', 'ACCNUM2EG', 'ENSEMBL2EG', 'ALIAS2EG')`:
#'    is the recommended order of sources to use when equating a gene
#'    symbol or identifier with Bioconductor gene annotations. They are
#'    attempted in order, to convert each entry to an ENTREZID ('EG').
#'       * When using `handle_multiple='best_each'`, each method in `try_list`
#'       is attempted, in order, for each gene entry, and the first
#'       successful result is retained as one or more `ENTREZID` values
#'       for each input gene.
#'       * Specifically, when `'SYMBOL2EG'` is successful, it is used.
#'       Only when `'SYMBOL2EG'` is not successful, `'ACCNUM2EG'` is
#'       attempted. Only when `'ACCNUM2EG'` is also unsuccessful,
#'       `'ENSEMBL2EG'` is used, and so on.
#'       * The lowest priority is `'ALIAS2EG'` since this association has
#'       the lowest confidence for association, and is most prone to
#'       multiple associations. In general pathway enrichment is somewhat
#'       tolerant to such lower-confidence associations, in return for
#'       having more potential opportunities for the correct association.
#'       * Finally, the matching ENTREZID values are converted to 'SYMBOL'
#'       so the "official gene SYMBOL" is associated with the original
#'       gene entry.
#'       * For each pathway, only the unique `'SYMBOL'` values are retained.
#'       * In the event two gene identifiers in a pathway are converted to
#'       the same `'SYMBOL'`, only the unique `'SYMBOL'` is associated with
#'       the corresponding pathways.
#'       * In the event one gene identifier is converted to two 'SYMBOL' values
#'       of equivalent priority, both 'SYMBOL' values are retained and
#'       associated with the corresponding pathways.  
#'       For this purpose, 'priority' refers to the order in `try_list`.
#'  
#' * `gene_colnames`: vector of columns that could be recognized to
#' contain gene symbols or identifiers.
#' * `pathway_colnames`: vector of columns that may contain pathway
#' names or identifiers. The first matching entry is assumed to contain
#' the pathway identifier associated with each gene set.
#' 
#' ## Handling other species
#' 
#' In the event the `pathwaydata` entries were derived from one species,
#' and you want to test enrichment using gene data from species.
#' In this case, the recommended order of operations:
#' 1. Convert `pathwaydata` to the most current gene symbols in the
#' species used for that pathway data, for example Human.
#' 2. Use another resource such as `babelgene::orthologs()` to convert
#' Human to the corresponding species; or use OrtholoGene to provide
#' a species-to-species mapping to the relevant gene SYMBOL.
#' 
#' 
#' ## The decision to use SYMBOL
#' 
#' The primary reason to use SYMBOL is to focus pathway analyses,
#' and comparison of pathway analyses, on the functional product
#' of gene loci in the cell: a protein product, or active transcript
#' (lncRNA, snRNA, snoRNA, miRNA, tRNA).
#' 
#' There are often multiple genetic copies of a gene locus which
#' produce identical transcript and/or protein products, and
#' the decision to use SYMBOL is motivated to focus on the product
#' and not the genetic source/coordinates of that product.
#' 
#' In Gene Ontology, Biological Process (GOBP), the gene
#' `'RN7SK'` (also known as '7S rRNA') is represented only once,
#' although there are over 200 copies of `'RN7SK'` in the
#' human genome assembly hg38 as of 2026.
#' 
#' In contrast, the 'Transcription Factor Target' (TFT) resource
#' of the Molescular Signature Database (MSigDB version 2025.1.Hs)
#' contains `'RN7SK'` and 97 distinct pseudogenes
#' named `'RN7SKP1'` through `'RN7SKP296'`.
#' The utility of each pseudogene for the purpose of gene set
#' enrichment is unclear, specifically whether enrichment of
#' one particular genetic copy of 7S rRNA is truly indicative
#' of transcription factor control, or indicative of cell type
#' specificity as represented in the 'Cell Type (C8)' signatures
#' also in MSigDB.
#' 
#' See [https://docs.gsea-msigdb.org/#MSigDB/Release_Notes/MSigDB_2025.1.Hs/](Release_Notes:MSigDB_2025.1.Hs)
#' 
#' The intent of 'MultiEnrichJam' is to use one gene symbol entry
#' to represent any and all copies of genetic locus associated with
#' that gene symbol.
#' 
#' As of November 2026, there are three gene SYMBOL values associated
#' to multiple ENTREZID values in human, which appear to be distinctly
#' different gene products:
#' |**SYMBOL** |**ENTREZID**     |**GENENAME**                                                                     |**ALIAS**   |
#' |:----------|:----------------|:--------------------------------------------------------------------------------|:-----------|
#' |TEC        |7006,100124696   |tec protein tyrosine kinase,transient erythroblastopenia of childhood            |PSCTK4,TEC  |
#' |HBD        |3045,100187828   |hemoglobin subunit delta,hypophosphatemic bone disease                           |HBD,HBK     |
#' |MMD2       |221938,100505381 |Miyoshi muscular dystrophy 2,monocyte to macrophage differentiation associated 2 |MMD2,PAQR10 |
#' 
#' @references Subramanian, Tamayo, et al. PNAS 2005, https://www.pnas.org/content/102/43/15545
#' 
#' @returns `data.frame` with two columns:
#'    * 'TERM': the pathway name
#'    * 'GENE': one gene per row, after updating genes when
#'    `params$refresh_genes=TRUE`.
#' 
#' @family jam enrich functions
#' 
#' @inheritParams jam_enricher
#' 
#' @export
process_pathwaydata <- function(
   pathwaydata,
   pathwaydata_params=list(
      refresh_genes=TRUE,
      ann_lib="org.Hs.eg.db",
      handle_multiple="best_each",
      empty_rule="empty",
      ignore.case=TRUE,
      gene_colnames=c(
         "gene_symbol",
         "db_gene_symbol",
         "gene",
         "ncbi_gene",
         "SYMBOL"
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
   verbose=FALSE,
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
            verbose=verbose,
            ...
         )
      })
      return(pathwaydata_list);
   }

   # process individual pathwaydata
   if (inherits(pathwaydata, c("data.frame", "tbl"))) {
      # convert data.frame,tbl to 'list'
      # then process the same as 'list'
      if (verbose) {
         jamba::printDebug("process_pathwaydata(): ",
            "data.frame input.")
      }
      use_gene_colname <- find_colname(
         params$gene_colnames,
         pathwaydata
      )
      if (length(use_gene_colname) == 0) {
         cli::cli_abort(c(
            "Could not assign {.var use_gene_colname} using: ",
            "{.val {params$gene_colnames}}"
         ))
      }
      pathway_colname <- find_colname(
         params$pathway_colnames,
         pathwaydata
      )
      if (length(pathway_colname) == 0) {
         cli::cli_abort(c(
            "Could not assign {.var pathway_colname} using: ",
            "{.val {params$pathway_colnames}}"
         ))
      }
      pathwaygenes <- pathwaydata[[use_gene_colname]];
      pathwaynames <- pathwaydata[[pathway_colname]];
      if (is.list(pathwaygenes)) {
         if (inherits(pathwaygenes, "AsIs")) {
            pathwaygenes <- unclass(pathwaygenes)
         }
         # Todo: Handle duplicated pathwaynames
         names(pathwaygenes) <- pathwaynames;
      } else {
         if (jamba::igrepHas("[;,/]", pathwaygenes)) {
            pathwaygenes <- strsplit(pathwaygenes, "[;,/ ]+")
            names(pathwaygenes) <- pathwaynames;
         } else {
            pathwaynames <- factor(
               pathwaynames,
               levels=unique(pathwaynames)
            )
            pathwaygenes <- split(pathwaygenes, pathwaynames);
         }
      }
      # Note: Error will occur if names are duplicated, this is expected
      # Todo: Check first and write useful error message.
      pathwaydata <- pathwaygenes;
      if (verbose > 1) {
         jamba::printDebug(
            "middle(pathwaydata, 20):"
         )
         print(jamba::middle(pathwaydata));
         print(table(lengths(pathwaydata)));
      }
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
#' ## Custom parameters
#' 
#' Argument `genedata_params` controls additional options:
#' 
#' * `refresh_genes=TRUE` enables `genejam::freshenGenes()` which
#' uses Bioconductor annotations to update the current gene annotation
#' to the current version available (to the active R session).  
#' Several arguments are passed to this function:
#'    * `ann_lib`: the core annotation library,
#'    'org.Hs.eg.db' is default for Human; 'org.Mm.eg.db' is used
#'    for Mouse. This annotation database must be available,
#'    for example install Rat data with
#'    `BiocManager::install('org.Rn.eg.db')`.
#'    Alternatively, `ann_lib` can be a `list` with one or more
#'    `AnnotationDbi` or `OrgDb` objects.
#'    * `handle_multiple='best_each'`: converts each gene to the best
#'    available match, using the priority order defined in `try_list`.
#'    * `ignore.case=TRUE`: Note that Human symbols are UPPERCASE,
#'    so for a moderate speed boost use `ignore.case=FALSE`.
#'    Note that Mouse by convention uses mixed case, for example "Apoe",
#'    and so for mouse and most other species, use `ignore.case=TRUE`.
#'    * `empty_rule='empty'`: converts any unmatched entry to '' or NA
#'    so it is not retained in the output. This decision prevents keeping
#'    a gene or accession which can never be recognized as an official
#'    gene during pathway enrichment analysis.  
#'    Use `empty_rule='original'` to retain the original entry whenever
#'    the entry cannot be matched to any other annotation in `try_list`.
#'    Note that this decision is intended only when the pathway data
#'    also contain similar gene entries which may also not have official
#'    gene symbols.
#'    * `try_list=c('SYMBOL2EG', 'ACCNUM2EG', 'ENSEMBL2EG', 'ALIAS2EG')`:
#'    is the recommended order of sources to use when equating a gene
#'    symbol or identifier with Bioconductor gene annotations. They are
#'    attempted in order, to convert each entry to an ENTREZID ('EG').
#'       * When using `handle_multiple='best_each'`, each method in `try_list`
#'       is attempted, in order, for each gene entry, and the first
#'       successful result is retained as one or more `ENTREZID` values
#'       for each input gene.
#'       * Specifically, when `'SYMBOL2EG'` is successful, it is used.
#'       Only when `'SYMBOL2EG'` is not successful, `'ACCNUM2EG'` is
#'       attempted. Only when `'ACCNUM2EG'` is also unsuccessful,
#'       `'ENSEMBL2EG'` is used, and so on.
#'       * The lowest priority is `'ALIAS2EG'` since this association has
#'       the lowest confidence for association, and is most prone to
#'       multiple associations. In general pathway enrichment is somewhat
#'       tolerant to such lower-confidence associations, in return for
#'       having more potential opportunities for the correct association.
#'       * Finally, the matching ENTREZID values are converted to 'SYMBOL'
#'       so the "official gene SYMBOL" is associated with the original
#'       gene entry.
#'       * For each gene result, only the first unique `'SYMBOL'` is retained,
#'       after applying the rank/sort operation. Typically the most
#'       significant result is retained for each gene.
#'       * In the event two gene identifiers are converted to
#'       the same `'SYMBOL'`, only the first unique `'SYMBOL'` is retained,
#'       after applying the rank/sort operation.
#'       Typically, the most significant result is retained for each gene.
#'       * In the event one gene identifier is converted to two 'SYMBOL' values
#'       of equivalent priority, the both 'SYMBOL' values are retained,
#'       by splitting the row into multiple rows which each contain one
#'       'SYMBOL' entry. Then the `'SYMBOL'` entries overall are sorted
#'       as described above, and the first unique `'SYMBOL'` is
#'       retained.
#'       Typically, the most significant result is retained for each gene.
#'  
#' * `gene_colnames`: vector of columns that could be recognized to
#' contain gene symbols or identifiers.
#' * `pvalue_cutoff`, `fold_cutoff`: used to define gene hits when
#' 'genedata' is `data.frame`. When not supplied, data can only
#' be used for ranked enrichment.
#' * `mgm_cutoff`: applied to max group mean (mgm) values to
#' require each gene hit to involve one group mean of at least
#' this value.
#' * `apply_mgm_universe`: When `TRUE` it restricts the gene universe
#' to genes with max group mean (mgm) of at least this value.
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
#' * Add optional `detected_genes`.
#'    * For now, this step is pushed onto the scientist.
#' 
#' * DONE. Add `mgm_cutoff` (max group mean)
#'    * Detect 'mgm' or 'mean' columns, calculate row max.
#'    * Apply `mgm_cutoff` to define gene hits.
#'    * Optional: Define universe using `mgm_cutoff` genes.
#'    * Optional: When choosing one row per gene, prefer `mgm_cutoff`
#'    before choosing most significant row?
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
#' @family jam enrich functions
#' 
#' @inheritParams jam_enricher
#' 
#' @export
process_genedata <- function(
   genedata,
   genedata_params=list(
      refresh_genes=TRUE,
      ann_lib="org.Hs.eg.db",
      handle_multiple="best_each",
      empty_rule="empty",
      ignore.case=TRUE,
      fold_cutoff=NULL,
      pvalue_cutoff=NULL,
      mgm_cutoff=NULL,
      apply_mgm_universe=TRUE,
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
      mgm_colnames=c(
         "mgm", "mean", "AveExpr"
      ),
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
      apply_mgm_universe <- match.arg(
         params$apply_mgm_universe,
         choices=eval(
            formals(process_genedata)$genedata_params$apply_mgm_universe
         )
      )
      fold_colname <- find_colname(params$fold_colnames, genedata)
      pvalue_colname <- find_colname(params$pvalue_colnames, genedata)
      mgm_colname <- find_colname(params$mgm_colnames, genedata, max=5);
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
      mgm_cutoff <- params$mgm_cutoff;
      # Apply P-value cutoff
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
      # Apply fold_cutoff
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
      # Apply mgm_cutoff
      met_mgm_cutoff <- rep_len(TRUE, nrow(genedata))
      if (length(mgm_cutoff) == 1 &&
         is.numeric(mgm_cutoff) &&
         length(mgm_colname) > 0) {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "apply mgm_cutoff:",
               mgm_cutoff
            )
         }
         if (length(mgm_colname) > 1) {
            mgm_data <- mgm_data
            mgm_values <- apply(mgm_data, 1, function(imgm){
               max(imgm, na.rm=TRUE)
            })
         } else {
            mgm_values <- genedata[[mgm_colname]]
         }
         met_mgm_cutoff <- (
            !is.na(mgm_values) &
            mgm_values >= mgm_cutoff
         )
         use_gene_hits <- (use_gene_hits & met_mgm_cutoff)
      } else {
         mgm_cutoff <- NULL;
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
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "sign by t-statistic: ",
               tstat_colname
            )
         }
         use_signs <- genedata[[tstat_colname]];
      } else if ("signed-P" %in% rank_strategy) {
         # Todo: check and convert to log2fold if needed
         use_signs <- (
            -log10(genedata[[pvalue_colname]]) *
               sign(genedata[[fold_colname]])
         );
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "sign by signed-P: ",
               c(pvalue_colname, fold_colname)
            )
         }
      } else if ("P-value" %in% rank_strategy) {
         use_signs <- -log10(genedata[[pvalue_colname]]);
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "sign by P-value: ",
               c(pvalue_colname)
            )
         }
      } else if ("log2fold" %in% rank_strategy) {
         # Todo: check and convert to log2fold if needed
         use_signs <- (genedata[[fold_colname]]);
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "sign by log2fold: ",
               c(fold_colname)
            )
         }
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
      met_mgm_cutoff <- met_mgm_cutoff[new_order];

      # Optionally filter universe using mgm_cutoff
      if (isTRUE(apply_mgm_universe) && !all(met_mgm_cutoff)) {
         if (verbose) {
            jamba::printDebug(
               "process_genedata(): ",
               "filtering gene universe by mgm_cutoff ",
               mgm_cutoff,
               ", retaining ",
               sep="",
               c(
                  jamba::formatInt(sum(met_mgm_cutoff)),
                  " of ",
                  jamba::formatInt(length(met_mgm_cutoff))
               ),
               " rows."
            )
         }
         use_signs <- use_signs[met_mgm_cutoff];
         use_genes <- use_genes[met_mgm_cutoff];
         use_gene_hits <- use_gene_hits[met_mgm_cutoff];
         met_mgm_cutoff <- met_mgm_cutoff[met_mgm_cutoff];
      }
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

   # remove entries with sign NA
   if (anyNA(gddf1$use_signs)) {
      jamba::printDebug(
         "process_genedata(): ",
         "remove entries with NA sign"
      )
      gddf1 <- subset(gddf1, !is.na(use_signs))
   }

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
