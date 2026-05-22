# Curate Ingenuity IPA colnames

Curate Ingenuity IPA colnames

## Usage

``` r
curateIPAcolnames(
  jDF,
  ipaNameGrep = c("^Name$", "^ID$", "Canonical Pathways", "Upstream Regulator",
    "Diseases or Functions Annotation", "Diseases . Functions", "My Lists",
    "Ingenuity Toxicity Lists", "My Pathways"),
  geneGrep = c("Molecules in Network", "Target molecules", "Molecules", "Symbol"),
  geneCurateFrom = c(" [(](complex|includes others)[)]", "^[,]+|[,]+$"),
  geneCurateTo = c("", ""),
  convert_ipa_slash = TRUE,
  ipa_slash_sep = ":",
  verbose = TRUE,
  ...
)
```

## Arguments

- jDF:

  data.frame from one Ingenuity IPA enrichment test.

- ipaNameGrep:

  vector of regular expression patterns used to recognize the name of
  the enriched entity, for example the biological pathway, or network,
  or disease category, etc.

- geneGrep:

  regular expression pattern used to recognize the column containing
  genes, or the molecules tested for enrichment which were found in the
  enriched entity.

- geneCurateFrom, geneCurateTo:

  vector of patterns and replacements, respectively, used to curate
  values in the gene column. These replacement rules are used to ensure
  that genes are delimited consistently, with no leading or trailing
  delimiters.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

This function is intended to help curate colnames observed in Ingenuity
IPA enrichment data. The IPA enrichment data includes multiple types of
enrichment tests, each with slightly different column headers. This
function is intended to make the colnames more consistent.

This function will rename the first recognized gene colname to
`"geneNames"` for consistency with downstream analyses.

The values in the recognized gene colname are curated using
`geneCurateFrom,geneCurateTo` for multiple pattern-replacement
substitutions. This mechanism is used to ensure consistent delimiters
and values used for each enrichment table.

Any colname matching `"-log.*p.value"` is considered -log10 P-value, and
is converted to normal P-values for consistency with downstream
analyses.

Any recognized P-value column is renamed to `"P-value"` for consistency
with downstream analyses.

When the recognized P-value column contains a range, for example
`"0.00017-0.0023"`, the lower P-value is chosen. In that case, the
higher P-value is stored in a new column `"max P-value"`. P-value ranges
are reported in the disease category analysis by Ingenuity IPA, after
collating individual pathways by disease category and storing the range
of enrichment P-values.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md),
[`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md),
[`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md),
[`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
[`find_enrich_colnames()`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md),
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
