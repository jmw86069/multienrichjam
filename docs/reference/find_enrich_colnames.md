# Find enrichment colnames

Find enrichment colnames

## Usage

``` r
find_enrich_colnames(
  x,
  keyColname = c("ID", "Name", "pathway", "itemsetID", "Description"),
  nameColname = c("Name", "pathway", "Description", "itemsetID", "ID"),
  descriptionColname = c("Description", "Name", "Pathway", "ID"),
  geneColname = c("geneID", "geneNames", "Genes"),
  countColname = c("gene_count", "count", "geneHits"),
  geneRatioColname = c("GeneRatio", "^Ratio"),
  pvalueColname = c("padjust", "p.adjust", "adjp", "padj", "qvalue", "qval", "q.value",
    "pvalue", "p.value", "pval", "FDR"),
  directionColname = c("activation.z.{0,1}score", "NES", "direction", "z.{0,1}score"),
  pathGenes = c("setSize", "pathGenes", "Count"),
  geneHits = c("Count", "geneHits", "gene_count"),
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `data.frame`, `enrichList`, `Mem`, or `list` of `data.frame` objects.

- keyColname:

  `character` default 'ID' indicating the primary identifier for each
  set. This column may be a numeric identifier.

- nameColname:

  `character` default 'Name', with the set name, typically a short name
  for each set.

- descriptionColname:

  `character` default 'Description' with the longer set description. It
  will use `nameColname` or `keyColname` as needed.

- geneColname:

  `character` default 'geneID' containing delimited genes associated
  with each enrichment result.

- countColname:

  `character` default 'Count' with the number of genes in the
  `geneColname` column. It will be calculated as needed.

- geneRatioColname:

  `character` default 'GeneRatio' with the numeric (decimal) ratio of
  test genes to pathway genes, or a character indication in the form
  '6/24' with 'tested/pathway' gene counts.

- pvalueColname:

  `character` default 'padjust' with the best available column to use
  for statistical significance of enrichment.

- directionColname:

  `character` default 'zscore' with a directional score, typically a
  z-score or some other reasonably scaled numeric value where the sign
  indicates directionality, with '+' meaning activated and '-' meaning
  suppressed.

- pathGenes:

  `character` default 'setSize' indicating the number of genes in each
  set as tested for enrichment. This number is not always reported,
  however it is not used by 'multienrichjam', but is required by some
  `clusterProfiler` functions.

- geneHits:

  `character` default 'Count' indicating the number of genes in the
  `geneColname` column. It will be calculated as needed.

- verbose:

  `logical` default FALSE, whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`character` of recognized colnames named by the type of column, using
`NA` for any column types not found.

`list` named by each column name argument, with one `character` value or
`NULL` in each entry.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md),
[`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md),
[`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md),
[`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md),
[`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
[`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md),
[`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md),
[`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md),
[`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md),
[`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md),
[`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md),
[`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md),
[`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md),
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
newborn_txt <- system.file("extdata",
   "Newborns-IPA.txt",
   package="multienrichjam");
ipa_dfs <- importIPAenrichment(newborn_txt);
find_enrich_colnames(ipa_dfs[[1]])
#> $keyColname
#> [1] "Name"
#> 
#> $nameColname
#> [1] "Name"
#> 
#> $descriptionColname
#> [1] "Name"
#> 
#> $geneColname
#> [1] "geneNames"
#> 
#> $countColname
#> NULL
#> 
#> $geneRatioColname
#> [1] "Ratio"
#> 
#> $pvalueColname
#> [1] "P-value"
#> 
#> $directionColname
#> NULL
#> 
#> $pathGenes
#> NULL
#> 
#> $geneHits
#> NULL
#> 

er <- enrichDF2enrichResult(ipa_dfs[[1]])
find_enrich_colnames(er)
#> $keyColname
#> [1] "ID"
#> 
#> $nameColname
#> [1] "Description"
#> 
#> $descriptionColname
#> [1] "Description"
#> 
#> $geneColname
#> [1] "geneID"
#> 
#> $countColname
#> [1] "Count"
#> 
#> $geneRatioColname
#> [1] "GeneRatio"
#> 
#> $pvalueColname
#> [1] "p.adjust"
#> 
#> $directionColname
#> NULL
#> 
#> $pathGenes
#> [1] "setSize"
#> 
#> $geneHits
#> [1] "Count"
#> 
```
