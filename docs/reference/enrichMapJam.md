# Create enrichMap igraph object from enrichResult, deprecated

Create enrichMap igraph object from enrichResult, deprecated, use
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

## Usage

``` r
enrichMapJam(
  x,
  n = 50,
  vertex.label.font = 1,
  vertex.label.cex = 1,
  keyColname = "ID",
  nodeLabel = c("Name", "Description", "ID"),
  descriptionColname = "Description",
  nodeLabelFunc = NULL,
  overlapThreshold = 0.2,
  msigdbGmtT = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `enrichResult` or `data.frame` containing enrichment results,
  specifically expecting colnames to contain one of
  `c("ID","Description","Name")` to represent the node name, and
  `c("Description")` to represent the description, if present.

- n:

  `numeric` value indicating the maximum number of nodes to include in
  the final network.

- vertex.label.font, vertex.label.cex:

  `numeric` values to define the default node label font and size.

- keyColname, nodeLabel, descriptionColname:

  `character` vectors indicating the colname to use for the node name
  and label.

- nodeLabelFunc:

  `function`, default NULL uses `fixSetLabels(x, width=30, ...)`. Use
  'FALSE' to keep labels unchanged. Note that the `igraph` node name
  remains the same as input, and label is used as the visual label on
  the graph.

  - When providing a custom function, note that each node is called
    individually and therefore the call is not vectorized.

- overlapThreshold:

  `numeric` value indicating the minimum Jaccard overlap, where edges
  with lower values are deleted from the `igraph` object.

- msigdbGmtT:

  not currently implemented.

- ...:

  additional arguments are passed to
  [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  when the input `x` is a `data.frame`.

## Value

`igraph` object with nodes representing each pathway, sized based upon
the number of genes involved in enrichment, and colored based upon the
`-log10(Pvalue)` using
[`colorjam::col_linear_xf()`](https://jmw86069.github.io/colorjam/reference/col_linear_xf.html),
a function that applies a color gradient to a numeric range. Each edge
has attributes: `overlap` containing Jaccard overlap, `overlap_count`
with the number of genes in common between the two nodes, and
`overlap_max_pct` with the maximum percent overlap between two nodes
(overlap count)/(smaller node size).

## Details

This function is a minor customization to
[`enrichplot::emapplot()`](https://rdrr.io/pkg/enrichplot/man/emapplot.html),
which takes a single `enrichResult` and produces an Enrichment map
network. For the equivalent function using `Mem` as input, see
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

The major differences compared with `enrichplot`:

- An `igraph` object is returned instead of `ggplot` object, since the
  `igraph` object can be manipulated and is useful to review.

- This function calculates overlap using `dist(...,method="binary")`
  which is a much faster method for calculating the Jaccard overlap.

- This function also calculates the overlap count, another helpful
  measure for filtering network connections, for example to remove links
  with only one gene, even if they overlap is above the required
  threshold. Many spurious network connections are removed with this
  filter, and it appears to be a helpful option.

## See also

Other jam deprecated functions:
[`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md),
[`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md),
[`mem_multienrichplot()`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md),
[`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
