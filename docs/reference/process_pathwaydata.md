# Process pathway data for enrichment using 'jam_enricher()'

Process pathway data for enrichment using 'jam_enricher()'

## Usage

``` r
process_pathwaydata(
  pathwaydata,
  pathwaydata_params = list(refresh_genes = TRUE, ann_lib = "org.Hs.eg.db",
    handle_multiple = "best_each", empty_rule = "empty", ignore.case = TRUE,
    gene_colnames = c("gene_symbol", "db_gene_symbol", "gene", "ncbi_gene", "SYMBOL"),
    pathway_colnames = c("gs_name", "gs_description", "name", "ID", "Pathway",
    "GO[:]*(BP|MF|CC)"), try_list = c("SYMBOL2EG", "ACCNUM2EG", "ENSEMBL2EG",
    "ALIAS2EG"), sep = "/"),
  verbose = FALSE,
  ...
)
```

## Arguments

- pathwaydata:

  one of the following:

  - **Pathways**: `list` of pathways, named by pathway, whose values are
    individual genes assigned to each pathway.

  - **List of Pathway Sources**: `list` named by pathway source, each
    containing a `list` of pathways as described above, named by pathway
    and containing pathway genes. Each Pathway Source will be analyzed
    one by one, independently.

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

`data.frame` with two columns:

- 'TERM': the pathway name

- 'GENE': one gene per row, after updating genes when
  `params$refresh_genes=TRUE`.

## Details

This function recognizes several types of pathway data as input, and
returns a consistent format as output: a `list` of pathway sources, each
containing a `data.frame` with columns 'TERM' and 'GENE'.

### Custom parameters

Argument `pathwaydata_params` controls additional options:

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
    any other annotation in `try_list`.

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

    - For each pathway, only the unique `'SYMBOL'` values are retained.

    - In the event two gene identifiers in a pathway are converted to
      the same `'SYMBOL'`, only the unique `'SYMBOL'` is associated with
      the corresponding pathways.

    - In the event one gene identifier is converted to two 'SYMBOL'
      values of equivalent priority, both 'SYMBOL' values are retained
      and associated with the corresponding pathways. For this purpose,
      'priority' refers to the order in `try_list`.

- `gene_colnames`: vector of columns that could be recognized to contain
  gene symbols or identifiers.

- `pathway_colnames`: vector of columns that may contain pathway names
  or identifiers. The first matching entry is assumed to contain the
  pathway identifier associated with each gene set.

### Handling other species

In the event the `pathwaydata` entries were derived from one species,
and you want to test enrichment using gene data from species. In this
case, the recommended order of operations:

1.  Convert `pathwaydata` to the most current gene symbols in the
    species used for that pathway data, for example Human.

2.  Use another resource such as
    [`babelgene::orthologs()`](https://igordot.github.io/babelgene/reference/orthologs.html)
    to convert Human to the corresponding species; or use OrtholoGene to
    provide a species-to-species mapping to the relevant gene SYMBOL.

### The decision to use SYMBOL

The primary reason to use SYMBOL is to focus pathway analyses, and
comparison of pathway analyses, on the functional product of gene loci
in the cell: a protein product, or active transcript (lncRNA, snRNA,
snoRNA, miRNA, tRNA).

There are often multiple genetic copies of a gene locus which produce
identical transcript and/or protein products, and the decision to use
SYMBOL is motivated to focus on the product and not the genetic
source/coordinates of that product.

In Gene Ontology, Biological Process (GOBP), the gene `'RN7SK'` (also
known as '7S rRNA') is represented only once, although there are over
200 copies of `'RN7SK'` in the human genome assembly hg38 as of 2026.

In contrast, the 'Transcription Factor Target' (TFT) resource of the
Molescular Signature Database (MSigDB version 2025.1.Hs) contains
`'RN7SK'` and 97 distinct pseudogenes named `'RN7SKP1'` through
`'RN7SKP296'`. The utility of each pseudogene for the purpose of gene
set enrichment is unclear, specifically whether enrichment of one
particular genetic copy of 7S rRNA is truly indicative of transcription
factor control, or indicative of cell type specificity as represented in
the 'Cell Type (C8)' signatures also in MSigDB.

See
[https://docs.gsea-msigdb.org/#MSigDB/Release_Notes/MSigDB_2025.1.Hs/](NA)

The intent of 'MultiEnrichJam' is to use one gene symbol entry to
represent any and all copies of genetic locus associated with that gene
symbol.

As of November 2026, there are three gene SYMBOL values associated to
multiple ENTREZID values in human, which appear to be distinctly
different gene products:

|  |  |  |  |
|----|----|----|----|
| **SYMBOL** | **ENTREZID** | **GENENAME** | **ALIAS** |
| TEC | 7006,100124696 | tec protein tyrosine kinase,transient erythroblastopenia of childhood | PSCTK4,TEC |
| HBD | 3045,100187828 | hemoglobin subunit delta,hypophosphatemic bone disease | HBD,HBK |
| MMD2 | 221938,100505381 | Miyoshi muscular dystrophy 2,monocyte to macrophage differentiation associated 2 | MMD2,PAQR10 |

## References

Subramanian, Tamayo, et al. PNAS 2005,
https://www.pnas.org/content/102/43/15545

## See also

Other jam enrich functions:
[`jam_enricher()`](https://jmw86069.github.io/multienrichjam/reference/jam_enricher.md),
[`process_genedata()`](https://jmw86069.github.io/multienrichjam/reference/process_genedata.md)
