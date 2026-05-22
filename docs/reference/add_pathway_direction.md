# Add directionality to pathway enrichment

Add directionality to pathway enrichment

## Usage

``` r
add_pathway_direction(
  x,
  gene_hits,
  delim = "/",
  dir_colname = "z-score",
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `enrichResult` object as returned by
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)

- gene_hits:

  `numeric` vector, with values `c(-1, 1)` to indicate up or down,
  respectively, and whose `names(gene_hits)` should match the gene
  symbols (or values) in `x@result$geneID`. If no `names(gene_hits)`
  match values in the `geneID` column, the direction score will
  ultimately be populated with `NA` values.

- delim:

  `character` string with delimiter used in the `enrichResult` object,
  default: `"/"`

- dir_colname:

  `character` string indicating the column name to create, default
  `"z-score"`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`enrichResult` object, with directional values populated into the result
column name defined by `dir_colname`. When no genes in a gene set are
present in `gene_hits` the expected directional value is `NA`.

## Details

This function uses pathway enrichment data provided as an `enrichResult`
object. Data is added as a new column (default "z-score") with `numeric`
values indicating the direction of change.

It requires `gene_hits` with the direction of change for each gene,
which is used to determine the up- and down-regulated genes in each gene
set.

The default equation uses the IPA z-score of activation, since the
z-score has beneficial properties in downstream visualization functions.
It is symmetric around zero, and implies reasonable thresholds.

Note that when no genes involved in enrichment are present in the
`gene_hits` data, the resulting directional value is expected to be
`NA`, which will be treated as "non-directional" in downstream
visualization functions.

### IPA z-score activation calculation

The IPA Activation `"z-score"` uses the following formula:

    z <- (N_genes_up - N_genes_down) / sqrt(N_genes_up + N_genes_down)

Citation for IPA formula: [Kramer
2014](https://doi.org/10.1093/bioinformatics/btt703)
https://doi.org/10.1093/bioinformatics/btt703, as referenced: [IPA FAQ -
Statistical
Calculations](https://qiagen.my.salesforce-sites.com/KnowledgeBase/KnowledgeNavigatorPage?id=kA41i000000L5nQCAS&categoryName=BioX)

### Direction calculation (deprecated)

The deprecated calculation for `"direction"` is shown below:

    direction_1 <- (N_genes_up - N_genes_down) / (N_genes);
    # then log2(1.2 + x) transformed
    direction <- log2(1.2 + abs(direction_1)) * sign(direction_1)

## See also

Other jam enrichment functions:
[`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
