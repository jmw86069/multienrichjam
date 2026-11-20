# Simple wrapper for pathway enrichment methods

Simple wrapper for pathway enrichment methods, to produce enrichment
results suitable for use with multienrichment analysis.

## Usage

``` r
jam_enricher(
  genedata,
  pathwaydata,
  enrich_method = c("ORA", "GSEA"),
  universe = NULL,
  genedata_params = list(ann_lib = "org.Hs.eg.db"),
  pathwaydata_params = list(ann_lib = "org.Hs.eg.db"),
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

- pathwaydata:

  one of the following:

  - **Pathways**: `list` of pathways, named by pathway, whose values are
    individual genes assigned to each pathway.

  - **List of Pathway Sources**: `list` named by pathway source, each
    containing a `list` of pathways as described above, named by pathway
    and containing pathway genes. Each Pathway Source will be analyzed
    one by one, independently.

- enrich_method:

  `character` string with enrichment method to perform:

  - **'ORA'**: uses
    [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
    for a standard hypergeometric (Fisher's Exact) test of enrichment.

  - **'GSEA'**: uses
    [`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html)
    for ranked gene set enrichment analysis. It calls 'fgsea' for
    fast-GSEA.

- universe:

  one of the following:

  - **`character`** vector of genes to define the background universe
    for the enrichment test. The same universe is used for all
    enrichment tests. When NULL as by default, it uses all genes defined
    in 'pathwaydata'.

  - `list` named by `names(genedata)` to define a universe for each
    individual entry in 'genedata'. Todo: Consider unique universe by
    genedata/pathwaydata.

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

- pathwaydata_params:

  `list` to manage 'pathwaydata'.

  - See 'genedata_params' for description of most fields.

  - pathway_colnames: `character` vector with patterns for preferred
    pathway column names, the first match with
    `find_colname(pattern, x)` is used as the pathway name. For the
    purpose of this process, only the pathway name is retained, and all
    other annotations should be retrieved separately as needed.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are passed to internal functions.

## Value

nested `list` of `enrichResult` objects.

- `list` named by 'genedata', each with:

- `list` named by 'pathwaydata', each with

- `list` named by 'enrich_method', each with

- 'ORA'=`enrichResult` or 'GSEA'=`gseaResult`

In addition, an attribute 'pathwaydata' is added with the updated
`pathwaydata` after
[`process_pathwaydata()`](https://jmw86069.github.io/multienrichjam/reference/process_pathwaydata.md)

## Details

### Workflow

1.  Iterate each 'genedata'

    - Validate data format.

    - Optional: Update gene annotations using
      [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html).

    - Produce `data.frame` with numeric `hit` and `sign`.

2.  Iterate each 'pathwaydata'

    - Optional: Update gene annotations using
      [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html).

    - Produce `data.frame` TERM2GENE suitable for 'clusterProfiler'.

Note that gene annotations can be synchronized to ensure the 'genedata'
and 'pathwaydata' both use the same version of annotation.

1.  Iterate each 'enrich_method'

    - ORA - over-representation analysis, aka hypergeometric enrichment,
      Fisher's Exact Test.

    - GSEA - gene set enrichment analysis, using fast-GSEA (fgsea)

2.  Collate results

    - Convert to `enrichResult` or `gseaResult`.

    - Return all results as a nested `list`, indexed by:

      - 'genedata'

      - 'pathwaydata'

      - 'enrich_method'

### Todo

- Return genedata after processing so it can be reviewed.

- Consider preparing hit incidence matrices, or lists.

  - ORA would use genedata

  - GSEA would use enrichment genes with input sign.

- Add example data and examples.

- Add 'testthat' functional unit tests.

- Consider multiple cutoffs, as a named list, only for ORA.

  - cutoffs=list( p01_fc1.5=c(pvalue_cutoff=0.01, fold_cutoff=1.5),
    p05_fc1.5=c(pvalue_cutoff=0.05, fold_cutoff=1.5) )

- Consider flattening the nested list output.

- Consider alternative enrichment methods:

  - 'clusterProfiler' enrichment functions

  - 'enrichR::enrichr()'

- Consider some custom "hooks":

  - Custom enrichment function with x,y,z input

  - Custom formatter to convert to `enrichResult` or `gseaResult`.

- Consider storing parameters as attributes:

  - 'ann_lib' and its package version for 'genejam'

  - `packageVersion('clusterProfiler')`, perhaps
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)

  - 'genedata_params'

  - 'pathwaydata_params'

  - `list` of 'genedata' `data.frame` used for enrichment

## See also

Other jam enrich functions:
[`process_genedata()`](https://jmw86069.github.io/multienrichjam/reference/process_genedata.md),
[`process_pathwaydata()`](https://jmw86069.github.io/multienrichjam/reference/process_pathwaydata.md)

## Examples

``` r
# To be defined
# genedata is data.frame with gene stat results, or vector of gene symbols
# pathwaydata is a list of gene symbol vectors, or data.frame with TERM2GENE
# or data.frame with gene column with delimited gene symbols.
```
