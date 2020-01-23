# TODO

This document describes plans for enhancements to the
multienrichjam R package.

## Freshen gene symbols

The `genejam` R package provides `freshenGenes()` to update
gene symbols to the most current nomenclature per gene.

Current recommendation is to run `genejam::freshenGenes()`
before running pathway enrichment. Install from Github:

> `devtools::install_github("jmw86069/genejam")`

* The goal is to ensure each input data is using the same
gene nomenclature, so genes are comparable.
* Running `freshenGenes()` is most useful when combining
pathway enrichment results from different source data,
for example two different platforms, or from two different
studies.

### Why use gene symbol?

There are multiple strategies widely used that attempt to maintain
stable references to genes, such as using an authoritative
gene identifier like NCBI Entrez Gene ID, or EnsEMBL Gene ID.
In these cases, the official gene symbol can be obtained from
resources like HGNC (Human Genome Nomenclature Committee).
No solution is perfect, but the current rationale for using
gene symbol in `multienrichjam`:

* NCBI Entrez Gene ID is an integer numeric value, and sometimes
multiple numbers refer to the same official gene symbol.

   * An integer value is risky as a primary identifier, because it
   is not individually recognizable as a gene identifier, for
   example the number `348` is not itself recognized as the
   identifier for gene `"APOE"`, and in fact the number `348`
   could mean anything by itself. However `"APOE"` is more
   recognizable as a gene symbol.

* Sometimes multiple EnsEMBL Gene ID values refer to the same
official gene symbol. Note: EnsEMBL Gene ID values are individually
unique, using the format `"ENSG00000000"` that is not
represented anywhere in the world except by EnsEMBL Genes.
However, few people in the world recognize `"ENSG00000130203"`
as the gene `"APOE"`.
* There is not a clean 1-to-1 mapping from NCBI Entrez Gene ID
to EnsEMBL Gene ID -- there are one-to-many, many-to-one, and
many-to-many relationships. In most cases, using official
gene symbol resolves these disputes, and reduces the need to
investigate these differences in detail.
* For cases where multiple Gene ID values refer to the same
official gene symbol, conversion to gene symbol is beneficial
to the outcome, because it ensures two alternate possible
input Gene ID values are converted to the same gene
symbol before comparisons.
* In general, gene sets are defined at the level of gene
product (often an expressed protein), and therefore most frequently
at the level of gene symbol. Most exceptions to this rule are
even less specific than one gene, for example KEGG defines its own
concept of a "gene" which may include multiple possible gene symbols.
Ingenuity IPA also converts some gene symbols to a small group
of "equivalent genes" for their purposes.
* Most gene set enrichment tools provide the
list of genes involved in each enrichment. If gene symbols are
not provided, the very next step is almost always to convert
the ID values to gene symbol for human interpretation. However,
many tools already provide gene symbols, which must be
compared back to the source data.


### Desired Outcomes

* The major goal is to ensure GeneA in one dataset is also called GeneA
in another dataset.
* Follow-up goal is to enable tracking from enrichment results
to the source data, for example seeing the fold change of genes
involved in enrichment results, which requires using the
correct gene identifier.


### Ideal Strategy: Run freshenGenes() upfront

* Potential "ideal strategy" is to run `freshenGenes()` on each
source data **before** running pathway enrichment. The major
benefit is that the source data would therefore already be
consistent across comparisons.


### Follow-up Strategy: Run freshenGenes() on each enrichment table

* The strategy would be to convert recognized gene symbols to
the most current annotation, with option to leave unrecognized
gene symbols as-is with no conversion.
* Option to provide a custom set of gene symbol curation events,
for example to reverse known conversions performed by software
such as Ingenuity IPA.

### Specific Challenges

* Some gene symbols will be split or duplicated after conversion.

   * A "split" occurs when one gene symbol refers to two actual genes,
   for example input gene `HSPA1` may ultimately refer to both
   `HSPA1A` and `HSPA1B`.
   * A duplication occurs when two different original gene symbols
   are resolved to the same final gene symbol, for example
   `HIST1H2BC`, and `HIST1H2BD` both refer to current
   official gene symbol `H2BC5`.
   * When a split or duplication occurs, it must be recorded.

* When a gene split or duplication occurs, it will affect
the number of genes returned in the enrichment result, and
is therefore not the same number used in the enrichment
calculation.

   * This effect is another reason to support converting gene
   symbols before running gene set enrichment.
   * However, the enrichment tool should only count each unique
   gene once, so the correct methodology is to accept the outcome
   of the enrichment tool within the caveats and limitations of
   using that tool. The purpose of modifying the gene symbol
   afterward is to facilitate comparison to the source data,
   and comparison across enrichment results. The underlying
   statistical enrichment test should remain unchanged.

* Each enrichment table may have its own distinct original annotation.
Therefore, each enrichment must record the gene update process,
in order to trace back to the source gene symbol.
* The "Use Case": given a gene involved in pathway enrichment,
find the source data from each input dataset.
