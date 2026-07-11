# multienrichjam — AI Assistant Notes

Project: `multienrichjam` R package (Bioconductor-style)  
Working directory: `c:/Users/james.ward/Documents/Projects/Ward/multienrichjam`

---

## Package overview

Analysis and visualisation of multiple gene-set enrichment results.
Key features: multi-enrichment heatmaps, EnrichmentMap networks, concept-network
(Cnet) plots, pathway clustering, and an S4 `Mem` class that is the central data
container for all downstream analyses.

Load for development: `devtools::load_all("c:/Users/james.ward/Documents/Projects/Ward/multienrichjam")`

---

## Key files

| File | Purpose |
|------|---------|
| `R/AllClasses.R` | `Mem` and `MemPlotFolio` S4 class definitions and validity checks |
| `R/AllGenerics.R` | All `setGeneric()` calls for the package |
| `R/methods-Mem.R` | All `setMethod()` implementations for `Mem` |
| `R/methods-MemPlotFolio.R` | Methods for `MemPlotFolio` |
| `R/jamenrich-get-mem.R` | `multiEnrichMap()` — the main function that builds a `Mem` object |
| `R/jamenrich-base.r` | Core enrichment utilities |
| `TODO.md` | Active task list (most-recent date block at top) |
| `multienrichjam_dev_notes.md` | Developer notes on package internals |
| `NEWS.md` | Changelog |

---

## The `Mem` S4 class

Central data container produced by `multiEnrichMap()`.

### Dimensions

| Symbol | Dimension | Example accessor |
|--------|-----------|-----------------|
| i | genes | `genes(x)` / `rownames(geneIM(x))` |
| j | sets (pathways) | `sets(x)` / `rownames(enrichIM(x))` |
| k | enrichments (comparison groups) | `enrichments(x)` / `names(enrichList(x))` |

`dim(Mem)` returns `c(n_genes, n_sets, n_enrichments)`.

### Key slots

| Slot | Type | Dimensions | Description |
|------|------|-----------|-------------|
| `memIM` | matrix | genes × sets | Gene-pathway incidence; magnitude = enrichment count |
| `enrichIM` | matrix | sets × enrichments | P-values per pathway per enrichment |
| `geneIM` | matrix | genes × enrichments | Gene membership per enrichment |
| `geneHitIM` | matrix | all_tested_genes × enrichments | All tested genes (including those not enriched) |
| `enrichList` | list | length k | Raw `enrichResult` objects, named by enrichment |
| `enrichLabels` | character | length k | Display labels for enrichments |
| `colorV` | character | length k | Named colour vector for enrichments |
| `geneHitList` | list | length k | Named list of character vectors — genes tested per enrichment |
| `enrichIMcolors` | matrix | sets × enrichments | Colour matrix for enrichIM |
| `enrichIMdirection` | matrix | sets × enrichments | Direction scores for enrichIM |
| `enrichIMgeneCount` | matrix | sets × enrichments | Gene counts per pathway per enrichment |
| `geneIMcolors` | matrix | genes × enrichments | Colour matrix for geneIM |
| `geneIMdirection` | matrix | genes × enrichments | Direction scores for geneIM |
| `multiEnrichDF` | data.frame | — | Combined enrichment results as a flat table |
| `multiEnrichResult` | enrichResult | — | Combined enrichResult (from first object) |
| `thresholds` | list | — | Analysis parameters (p_cutoff, topEnrichN, min_count, …) |
| `headers` | list | — | Column name mappings (geneColname, nameColname, pvalueColname, …) |

Slots without a dedicated accessor function must be accessed directly as `x@slotname`:
`multiEnrichDF`, `multiEnrichResult`, `enrichLabels`.

### Validity check

`check_Mem(object)` — returns TRUE/FALSE; also called automatically by `validObject()`.

---

## Implemented `Mem` methods (selected)

| Method | Description |
|--------|-------------|
| `x[i, j, k]` | Subset by genes (i), sets (j), enrichments (k) |
| `show(x)` | Pretty-print summary |
| `dim(x)` | `c(n_genes, n_sets, n_enrichments)` |
| `dimnames(x)` | Named list: genes, sets, enrichments |
| `names(x)` | Enrichment names (same as `enrichments(x)`) |
| `c(Mem1, Mem2, ...)` | Combine along enrichment dimension — see below |
| `as(x, "list")` / `Mem_to_list(x)` | Coerce to legacy list format |
| `as(list, "Mem")` / `list_to_Mem(x)` | Coerce from legacy list format |
| `updateObject(x)` | Update serialised Mem to current class version |
| `EnrichmentMap(x)` | Build/plot EnrichmentMap igraph network |
| `geneInCategory(x)` | List of genes per pathway |
| `setsByGene(x)` | List of pathways per gene |

---

## `c()` for Mem objects — implementation details (added 10 Jul 2026)

### Usage

```r
Mem_combined <- c(Mem1, Mem2, Mem3)   # any number of Mem objects
```

### Rules

- **Enrichments** must be non-overlapping across all inputs → error if any overlap.
- **Genes and sets** are **unioned** across all inputs and sorted with `jamba::mixedSort()`.
- Matrices are expanded to the union dimensions with these fill values:

  | Matrix type | Fill for new rows/cols |
  |-------------|----------------------|
  | Numeric (`geneIM`, `enrichIMgeneCount`, `enrichIMdirection`, `memIM`) | `0` |
  | P-value (`enrichIM`) | `NA` (pathway not enriched in that comparison) |
  | Colour (`geneIMcolors`, `enrichIMcolors`) | `NA_character_` |

- `memIM` is combined by **adding** the two expanded matrices. For unique sets this
  propagates values directly; for shared sets gene membership accumulates correctly.
- `geneHitList` is concatenated (named by enrichment, so no overlap).
- **Metadata** (`thresholds`, `headers`, `multiEnrichDF`, `multiEnrichResult`) is taken
  from the **first** object; other objects' metadata is silently discarded.

### Implementation

- `setMethod("c", "Mem", ...)` in `R/methods-Mem.R`
  (base `c()` is an implicit S4 generic; dispatches on class of first arg)
- Internal helper `combine_mem_objects(x, y)` does the pairwise work
- Internal helper `expand_matrix(mat, new_rows, new_cols, fill)` expands a matrix
- 18 unit tests in `tests/testthat/test-Mem-combine.R`

### Why not `S4Vectors::combine()`?

Base R's `c()` does **not** dispatch to `combine()` — it collects into a list.
The correct approach for S4 is `setMethod("c", "Signature", ...)`.
No new generic needs to be added to `AllGenerics.R` for this.

---

## Dev workflow notes

- Dev notes file: `multienrichjam_dev_notes.md`
- Adding igraph shapes: after `devtools::load_all()`, re-register shapes manually
  (see `zzz.R` and dev notes).
- Adding package data: create object in R session, then `usethis::use_data(obj, overwrite=TRUE)`.
- The deprecated argument `cutoffRowMinP` in `multiEnrichMap()` should be replaced
  with `p_cutoff`.
