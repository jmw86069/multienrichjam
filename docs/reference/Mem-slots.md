# Mem description of slots

Mem description of slots, for reference, slots may change at any time

## Details

Note that the Mem S4 object is intended to be used via the accessor
functions as described in `Mem-class`.

However, for reference, and in case *really necessary* the slots are
described in more detail below. These slots may change at any time, as
the 'Mem' S4 object evolves after more usage, so relying upon specific
slots to be stable may not be ideal.

The descriptions below should be edited whenever the 'Mem' S4 slot
structure changes, but be prepared there may sometimes be some delay
even under the best of intentions.

## Enrichment data input

- **enrichList**: `list` of `enrichResults` named by enrichment test.
  This list represents the input enrichment data.

- **enrichLabels**: `character` vector with custom enrichment labels,
  named by enrichment test.

- **geneHitList**: `list` of vectors named by enrichment test. Vectors
  are either:  
  `character` vector of genes tested as hits for enrichment, or  
  `integer` vector named by genes with values indicating the
  directionality, using `1` or `-1`.

- **geneHitIM**: `numeric` matrix with enrichment tests as colnames, and
  genes as rownames.  
  Values are `1` or `0` to indicate gene hits tested,  
  optionally `1` or `-1` to indicate directionality.

## Incidence matrix data

- **memIM**: `numeric` matrix with gene sets / pathways as colnames, and
  genes as rownames.  
  Values are `1` or `0` indicating genes present in each pathway.

- **geneIM**: `numeric` matrix with enrichment test as colnames, and
  genes as rownames.  
  This matrix differs from 'geneHitIM' in that it only includes genes
  also in 'memIM', and in that order.

- **enrichIM**: `numeric` matrix with enrichment test as colnames, and
  gene sets / pathways as rownames.  
  Matrix values are enrichment P-values, using the appropriate adjusted
  P-value or FDR to represent significance.

## Combined enrichment data

Enrichment data are combined into singular objects, mainly intended to
be used with `clusterProfiler` related functions. Also the `data.frame`
could be used as a supplemental table.

- **multiEnrichDF**: `data.frame` with single combined enrichment
  table.  
  Data include union of genes involved in each pathway, and lowest
  P-value observed across all enrichment tests.

- **multiEnrichResult**: `enrichResult` equivalent of 'multiEnrichDF'.  
  It is intended to facilitate use in `enrichplot` functions.

### Analysis details

- **colorV**: `character` vector of R colors, named by enrichment test.

- **thresholds**: `list` of thresholds used when calling
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

  - **p_cutoff**: `numeric` P-value threshold for significant
    enrichment.

  - **min_count**: `integer` number of genes involved in a pathway for
    the pathway to be considered significant.

- **headers**: `list` to associate actual colnames (headers) in
  'enrichList' with the conceptual names used in multienrichjam.

### Additional incidence matrix data

- **geneIMcolors**: `character` matrix with enrichment test colnames,
  and gene rownames.  
  Values are colors assigned using 'colorV'.

- **geneIMdirection**: `numeric` matrix equivalent to 'geneIM' with
  values `1` and `-1` indicating directionality.  
  When no directionality is provided, all values will be `1` or `0`.

- **enrichIMcolors**: `character` matrix with enrichment test colnames,
  and gene set / pathway rownames.  
  Values are color gradients using the appropriate P-value threshold
  used with
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

- **enrichIMdirection**: `numeric` matrix equivalent to 'enrichIM' with
  values indicating directional score.  
  Directional scores may include a IPA z-score or other directional
  score where available.  
  When not available, all values are `1`.

- **enrichIMgeneCount**: `integer` matrix equivalent to 'enrichIM'.  
  Values are the number of genes involved in the enrichment of each
  pathway.

## See also

Other Mem:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`Mem-data`](https://jmw86069.github.io/multienrichjam/reference/Mem-data.md),
[`check_Mem()`](https://jmw86069.github.io/multienrichjam/reference/check_Mem.md)
