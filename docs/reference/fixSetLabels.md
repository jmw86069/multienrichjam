# Fix Set or pathway labels for legibility

Fix Set or pathway labels for legibility

## Usage

``` r
fixSetLabels(
  x,
  wrap = TRUE,
  width = 40,
  maxNchar = Inf,
  suffix = "...",
  nodeType = c("Set", "Gene", "any"),
  do_abbreviations = TRUE,
  adjustCase = TRUE,
  lowercaseAll = TRUE,
  removeGrep = c("^(KEGG(_MEDICUS|)(_REFERENCE|_VARIANT|))[_. ]",
    "^(PID|REACTOME|BIOCARTA|NABA|SA|SIG|ST|WP|HALLMARK)[_. ]",
    "^((Mmu|Hsa|Rno|Dme)[0-9]+|Wp[0-9]+|^M[0-9]+)[:]",
    "^(R-(Hsa|Mmu|Rno|Dme)-[0-9]+)[:]"),
  words_from = NULL,
  words_to = NULL,
  add_from = NULL,
  add_to = NULL,
  abbrev_from = NULL,
  abbrev_to = NULL,
  perl = TRUE,
  makeUnique = FALSE,
  ...
)

words

abbrev
```

## Format

`words` is a `data.frame` with colnames 'from' and 'to', for pathway
word Perl-compatible regular expression pattern and replacement. It is
used as default in `fixSetLabels(..., words_from)`.

`abbrev` is a `data.frame` with colnames 'from' and 'to', for pathway
word Perl-compatible regular expression pattern and replacement. It is
used as default in `fixSetLabels(..., abbrev_from)`.

## Arguments

- x:

  any of the following objects:

  - `character` vector

  - `igraph` object. The `igraph::V(g)$name` attribute is used as input,
    and the resulting label is then stored as `V(g)$label`. When
    `nodeType` is also defined, and nodes have attribute 'nodeType',
    only nodes with that attribute value will be edited. The default is
    `nodeType="Set"`.

  - `Mem` object. The
    [`sets()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
    are adjusted by this function.

- wrap:

  `logical` indicating whether to apply word wrap, based upon the
  supplied `width` argument.

- width:

  integer value used when `wrap=TRUE`, it is sent to
  [`base::strwrap()`](https://rdrr.io/r/base/strwrap.html).

- maxNchar:

  `numeric` value or `Inf` to limit the maximum characters allowed for
  each string. This option is preferred when `wrap=TRUE` is not
  feasible, for example heatmap labels. When `NULL` or `Inf` no limit is
  applied. See [`base::nchar()`](https://rdrr.io/r/base/nchar.html).

- suffix:

  `character` value, default `"..."`, used when `maxNchar` is below
  `Inf`. When a string is shortened to `maxNchar`, the `suffix` helps
  indicate that there was additional text.

- nodeType:

  `character` string ussed when `x` is `igraph`, to limit changes to
  nodes by attribute values in `"nodeType"`. Use `"any"` or `NULL` to
  affect all nodes.

- do_abbreviations:

  `logical`, default TRUE, whether to apply `abbrev_from,abbrev_to`.
  These patterns are intended specifically to help shorten a long
  phrase, possibly removing words, or using common abbreviations.

- adjustCase:

  `logical`, default TRUE, indicating whether to adjust the uppercase
  and lowercase lettering by calling
  [`jamba::ucfirst()`](https://jmw86069.github.io/jamba/reference/ucfirst.html).
  The default sets all characters to lowercase, then applies uppercase
  to the first letter of each word.

- lowercaseAll:

  `logical` used only when `adjustCase=TRUE`, passed to
  [`jamba::ucfirst()`](https://jmw86069.github.io/jamba/reference/ucfirst.html)

- removeGrep:

  `character` regular expression pattern used to remove patterns from
  the resulting label. If given a vector, it will iterate each value
  individually.

  - The default removes common canonical pathway source prefix terms use
    in MSigDB data, for example KEGG, BIOCARTA, PID, etc. Use `""` or
    `NULL` to skip this step.

  - Multiple values can be defined, they are applied in order.

- words_from, words_to:

  `character` default NULL uses internal data `words` with pattern and
  replacement.

  - Input `words_from` can be a two-column `data.frame` expected to have
    'from' and 'to' in order.

  - Supplied as vectors, the 'words_from' are regular expression
    patterns, replaced with 'words_to' in order, applied case-sensitive.
    It does use Perl regular expression in
    [`base::gsub()`](https://rdrr.io/r/base/grep.html), which is useful
    to use with 'backslash-b' to enforce a word boundary for example.

- add_from, add_to:

  `character` vectors used in addition to `words_from`,`words_to`.

  - 'add_from' can be supplied as a two-column `data.frame` as described
    for 'words_from'.

  - These values are applied after `words_from`,`words_to`, so that
    user-defined replacements have priority.

- abbrev_from, abbrev_to:

  `character` default NULL uses internal data `abbrev` with pattern and
  replacement. Intended to apply a specific abbreviation, and only
  applied when `do_abbreviations=TRUE`.

  - 'abbrev_from' can be supplied as a `data.frame` as described for
    'words_from'.

  - The abbreviations are "opinionated" in that they may remove words or
    shorten common phrases which do not seem critical to understanding
    the meaning of most biological pathways.

  Examples:

  - "Extracellular Matrix" becomes "ECM"

  - "Mitochondrial" becomes "Mito"

  - " Pathway" at the end of a phrase is removed, as it is not required
    to understand the rest of the label.

  - "Signaling by " at the start of a phrase is removed, as it also is
    not typically necessary to understand the label.

- perl:

  `logical` default TRUE, passed to
  [`gsub()`](https://rdrr.io/r/base/grep.html) for pattern matching.
  When Perl-mode is enabled, it also enforces word boundaries before and
  after each pattern. When Perl-mode is not enabled, there are no word
  boundary conditions applied.

- makeUnique:

  `logical` default FALSE, whether to make resulting labels unique,
  using
  [`jamba::makeNames()`](https://jmw86069.github.io/jamba/reference/makeNames.html).
  This option is necessary when the resulting names will be applied back
  to the object `x`, and therefore cannot permit duplicated values.

  - When input 'x' is `Mem` it forces `makeUnique=TRUE`.

  - For `igraph` input, it does not force `makeUnique=TRUE` because the
    result is to populate vertex 'labels' which permits duplicated
    values.

  - For all other 'x' types, `makeUnique` is used as provided.

- ...:

  additional arguments are passed to `jamba::ucfirst(x, ...)`, for
  example `firstWordOnly=TRUE` will capitalize only the first word.

## Value

object whose class matches input 'x':

- `character` vector

- `igraph` object with vertex 'label' updated from 'name'

- `Mem` object with updated
  [`sets()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)

## Details

This function is a convenient wrapper for several steps that edit gene
set and pathways labels to be slightly more legible. It operates on:

- `character` vector, returning `character` vector

- `igraph` object, where is uses 'name' to update 'label'

- `Mem` object, where it updates `sets(Mem)`

The arguments have extensive default values encoded, which are
represented in data `multienrichjam::words` for basic word replacement,
and `multienrichjam::abbrev` for abbreviations used only when
`do_abbreviations=TRUE`.

Summary of typical changes:

- The vast majority of changes are custom biological terms which are
  expected to have certain capitalization, for example 'Mapk' is usually
  written 'MAPK'.

- Some changes are motivated to fix common artifacts in public data, for
  example `'PI3kakt'` refers to `'PI3K/AKT'`.

- To use your own replacements, supply `words_from` as a two-column
  `data.frame`, or two vectors `words_from` and `words_to`.

- To add custom effects to default, supply `abbrev_from` as a two-column
  `data.frame`, or use two vectors `abbrev_from` and `abbrev_to`.

For `igraph` input, the vertex 'name' is used as the starting point. To
revert changes, use `igraph::V(x)$label <- igraph::V(x)$name`.

For `Mem` input, the `sets(x)` are updated, with no immediate way to
revert changes. It may become useful to do so in future, however.

## See also

Other jam Mem utilities:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
[`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md),
[`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md),
[`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)

## Examples

``` r
x <- c("KEGG_INSULIN_SIGNALING_PATHWAY",
   "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
   "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY");
fixSetLabels(x);
#> [1] "Insulin Signaling"         "T-cell Receptor Signaling"
#> [3] "Neurotrophin Signaling"   
fixSetLabels(x, do_abbreviations=FALSE);
#> [1] "Insulin Signaling Pathway"         "T-cell Receptor Signaling Pathway"
#> [3] "Neurotrophin Signaling Pathway"   

jamba::nullPlot();
jamba::drawLabels(txt=x,
   preset=c("top", "center", "bottom"));

```
