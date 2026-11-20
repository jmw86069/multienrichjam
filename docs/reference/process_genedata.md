# Process gene data for enrichment

Process gene data for enrichment

## Usage

``` r
process_genedata(
  genedata,
  genedata_params = list(refresh_genes = TRUE, ann_lib = "org.Hs.eg.db", handle_multiple
    = "best_each", empty_rule = "empty", ignore.case = TRUE, fold_cutoff = NULL,
    pvalue_cutoff = NULL, mgm_cutoff = NULL, apply_mgm_universe = TRUE, gene_colnames =
    c("SYMBOL", "probes", "gene", "gene_symbol", "geneid", "gene_id", "db_gene_symbol",
    "ncbi_gene"), fold_colnames = c("fold", "logFC", "lfc", "log2FC", "log2fold"),
    pvalue_colnames = c("adj.p.val", "adj.*P.*val", "FDR", "q.val|qval", "adjp",
    "padjust", "p.adjust", "p.value", 
     "pval"), mgm_colnames = c("mgm", "mean",
    "AveExpr"), tstat_colnames = c("t", "t.stat"), rank_strategy = c("signed-P",
    "t-statistic", "P-value", "log2fold"), try_list = c("SYMBOL2EG", "ACCNUM2EG",
    "ENSEMBL2EG", "ALIAS2EG"), sep = "/"),
  verbose = FALSE,
  ...
)
```

## Arguments

- genedata:

  `list` named by the Enrichment condition, where each list element
  contains a suitable form of data. Each 'genedata' will be tested for
  enrichment independently.

  Expected formats:

  - **`character`** vector of genes, which represent "gene hits".
    Compatible with 'ORA' and cannot be used with 'GSEA'.

  - **`integer`** vector of integer signs limited to `1` or `-1`, named
    by gene. Compatible with 'ORA' and cannot be used with 'GSEA'. The
    sign applies directionality to genes, is used to define 'direction'
    in the pathway, but not used for enrichment itself.

  - **`numeric`** with scalar values, named by gene. Compatible with
    'GSEA' and cannot be used with 'ORA'. The `numeric` value ranks the
    genes, and is applied to directionality of the results.

  - **`data.frame`** of gene stats, with columns: 'gene', 'pvalue',
    'fold' in some conceptual form. Compatible with 'ORA' and 'GSEA'.

- genedata_params:

  `list` to manage 'genedata'.

  - **refresh_genes**: `logical` whether to update gene symbols using
    [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html)

  - **ann_lib**: `character` with Bioconductor library to use when
    'refresh_genes' is TRUE.

  - **handle_multiple**: `character` string in
    [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html),
    default "best_each" uses the best match for each gene symbol per
    row.

  - **empty_rule**: `character` how to handle empty gene 'SYMBOL'.
    Default 'empty' leaves the value blank. Other options are 'original'
    to keep the original value, or 'na' to insert `NA`.

  - **pvalue_cutoff**,**fold_cutoff**: `numeric` thresholds used when
    `genedata` is supplied as a `data.frame`. It will define gene
    statistical hits using these thresholds. Default NULL will not apply
    thresholds, and will not define hits.

  - **gene_colnames**,**pvalue_colnames**,**fold_colnames**: `character`
    vectors with patterns used to match column names when 'genedata'
    contains a `data.frame`. These columns are required. The
    'gene_colname' is used as the initial gene, and will be used when
    `refresh_genes=TRUE`.

  - **tstat_colnames**: `character` (optional) pattern for t-statistic
    column name when 'genedata' contains a `data.frame`. It is only used
    when `rank_strategy='t-statistic'`.

  - **rank_strategy**: `character` string used when `genedata` contains
    a `data.frame`:

    - **'signed-P'** (default) uses the `-log10(Pvalue)` multiplied by
      the sign of the fold change. It produces the same rank order as
      using the t-statistic, useful when the `'t'` column is not
      available.

    - **'t-statistic'** uses the t-statistical column.

    - **'P-value'** uses the `-log10(Pvalue)` to rank entries.

    - **'log2fold'** uses the fold column to rank entries. Although a
      common default, we recommend the signed-P for improved resolution.

  - **try_list**: passed to
    [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html)
    to determine which annotation methods to search when resolving to
    gene SYMBOL. They are handled in priority order, there 'best_each'
    will return the first (best) result for each input gene. When input
    genes are entirely integer values, they are assumed to be ENTREZID
    values, which skips this step.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are passed to internal functions.

## Value

`data.frame` with columns:

- 'gene': the gene symbol

  - `genedata_params$refresh_genes=FALSE` (default) uses existing gene
    symbol from 'gene_colnames'

  - `genedata_params$refresh_genes=FALSE` uses 'SYMBOL' after applying
    [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html).
    **Make sure to use correct 'ann_lib'.**

- 'sign': the `numeric` sign or scalar value:

  - entirely '1' for non-signed gene hits, cannot be used for ranked
    enrichment.

  - '1' or '-1' for signed gene hits, cannot be used for ranked
    enrichment.

  - `numeric` scalar values used for ranked enrichment.

- 'hit': `logical` indicating whether each entry met the statistical
  cutoffs `pvalue_cutoff` and `fold_cutoff` when defined. When all
  entries are FALSE there were no cutoffs applied, and data should only
  be used for ranked enrichment.

- The names are assigned using gene symbols, either as supplied or after
  resolution when `genedata_params$refresh_genes=TRUE`.

- Finally, attribute 'compatible_method' indicates which enrichment
  tests are suitable.

## Details

Workflow:detect_genedata_format

- Define the gene column

- Define the sign, use "1" for no sign

  - Calculate sign for "signed-P"

  - Use 't' for "t-statistic"

- Refresh gene symbols to 'SYMBOL'

- Remove entries with no 'SYMBOL'

- Deconcatenate / Expand 'SYMBOL' with multiple genes

- Resolve duplicate 'SYMBOL' to one row

  - Sort by absolute value of sign descending, then sign descending.

  - Retain the first observation per 'SYMBOL'

### Custom parameters

Argument `genedata_params` controls additional options:

- `refresh_genes=TRUE` enables
  [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html)
  which uses Bioconductor annotations to update the current gene
  annotation to the current version available (to the active R session).
  Several arguments are passed to this function:

  - `ann_lib`: the core annotation library, 'org.Hs.eg.db' is default
    for Human; 'org.Mm.eg.db' is used for Mouse. This annotation
    database must be available, for example install Rat data with
    `BiocManager::install('org.Rn.eg.db')`. Alternatively, `ann_lib` can
    be a `list` with one or more `AnnotationDbi` or `OrgDb` objects.

  - `handle_multiple='best_each'`: converts each gene to the best
    available match, using the priority order defined in `try_list`.

  - `ignore.case=TRUE`: Note that Human symbols are UPPERCASE, so for a
    moderate speed boost use `ignore.case=FALSE`. Note that Mouse by
    convention uses mixed case, for example "Apoe", and so for mouse and
    most other species, use `ignore.case=TRUE`.

  - `empty_rule='empty'`: converts any unmatched entry to ” or NA so it
    is not retained in the output. This decision prevents keeping a gene
    or accession which can never be recognized as an official gene
    during pathway enrichment analysis. Use `empty_rule='original'` to
    retain the original entry whenever the entry cannot be matched to
    any other annotation in `try_list`. Note that this decision is
    intended only when the pathway data also contain similar gene
    entries which may also not have official gene symbols.

  - `try_list=c('SYMBOL2EG', 'ACCNUM2EG', 'ENSEMBL2EG', 'ALIAS2EG')`: is
    the recommended order of sources to use when equating a gene symbol
    or identifier with Bioconductor gene annotations. They are attempted
    in order, to convert each entry to an ENTREZID ('EG').

    - When using `handle_multiple='best_each'`, each method in
      `try_list` is attempted, in order, for each gene entry, and the
      first successful result is retained as one or more `ENTREZID`
      values for each input gene.

    - Specifically, when `'SYMBOL2EG'` is successful, it is used. Only
      when `'SYMBOL2EG'` is not successful, `'ACCNUM2EG'` is attempted.
      Only when `'ACCNUM2EG'` is also unsuccessful, `'ENSEMBL2EG'` is
      used, and so on.

    - The lowest priority is `'ALIAS2EG'` since this association has the
      lowest confidence for association, and is most prone to multiple
      associations. In general pathway enrichment is somewhat tolerant
      to such lower-confidence associations, in return for having more
      potential opportunities for the correct association.

    - Finally, the matching ENTREZID values are converted to 'SYMBOL' so
      the "official gene SYMBOL" is associated with the original gene
      entry.

    - For each gene result, only the first unique `'SYMBOL'` is
      retained, after applying the rank/sort operation. Typically the
      most significant result is retained for each gene.

    - In the event two gene identifiers are converted to the same
      `'SYMBOL'`, only the first unique `'SYMBOL'` is retained, after
      applying the rank/sort operation. Typically, the most significant
      result is retained for each gene.

    - In the event one gene identifier is converted to two 'SYMBOL'
      values of equivalent priority, the both 'SYMBOL' values are
      retained, by splitting the row into multiple rows which each
      contain one 'SYMBOL' entry. Then the `'SYMBOL'` entries overall
      are sorted as described above, and the first unique `'SYMBOL'` is
      retained. Typically, the most significant result is retained for
      each gene.

- `gene_colnames`: vector of columns that could be recognized to contain
  gene symbols or identifiers.

- `pvalue_cutoff`, `fold_cutoff`: used to define gene hits when
  'genedata' is `data.frame`. When not supplied, data can only be used
  for ranked enrichment.

- `mgm_cutoff`: applied to max group mean (mgm) values to require each
  gene hit to involve one group mean of at least this value.

- `apply_mgm_universe`: When `TRUE` it restricts the gene universe to
  genes with max group mean (mgm) of at least this value.

### Todo

- Provide some kind of summary for gene resolution

  - Total gene values

  - Number of distinct input genes

  - Percent matched, percent retained

  - Number of multi-gene entries

  - Number of duplicated genes

  - Number of distinct output genes

- Add optional `detected_genes`.

  - For now, this step is pushed onto the scientist.

- DONE. Add `mgm_cutoff` (max group mean)

  - Detect 'mgm' or 'mean' columns, calculate row max.

  - Apply `mgm_cutoff` to define gene hits.

  - Optional: Define universe using `mgm_cutoff` genes.

  - Optional: When choosing one row per gene, prefer `mgm_cutoff` before
    choosing most significant row?

## See also

Other jam enrich functions:
[`jam_enricher()`](https://jmw86069.github.io/multienrichjam/reference/jam_enricher.md),
[`process_pathwaydata()`](https://jmw86069.github.io/multienrichjam/reference/process_pathwaydata.md)
