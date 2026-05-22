# Import Ingenuity Pathway Analysis 'IPA' results

Import Ingenuity Pathway Analysis 'IPA' results, by default reverting
IPA symbols to input values

## Usage

``` r
importIPAenrichment(
  ipaFile,
  headerGrep =
    "(^|\t)((expr.|-log.|)p-value|Pvalue|Score($|\t)|Symbol($|\t)|Ratio($|\t)|Consistency.Score|Master.Regulator($|\t))",
  ipaNameGrep = c("Pathway", "Regulator$", "Regulators", "Regulator", "Disease",
    "Toxicity", "Category", "Categories", "Function", "Symbol$", "^ID$",
    "My.(Lists|Pathways)"),
  geneGrep = c("Molecules in Network", "Target molecules", "Molecules", "Symbol"),
  geneCurateFrom = c("[ ]*[(](complex|includes others)[)][ ]*", "^[, ]+|[, ]+$"),
  geneCurateTo = c("", ""),
  signColname = c("Expr Fold Change", "fold.*change", "log.*ratio", "log.*fold",
    "log.*fc", "lfc", "ratio", "fold", "fc"),
  signThreshold = 0,
  method = 1,
  sheet = 1,
  sep = "\t",
  xlsxMultiSheet = TRUE,
  useXlsxSheetNames = FALSE,
  remove_blank_colnames = TRUE,
  convert_ipa_slash = TRUE,
  ipa_slash_sep = ":",
  revert_ipa_xref = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- ipaFile:

  one of the four input types described above: a character vector of
  text file names; a character vector of Excel `.xlsx` file names; a
  list of `data.frame` objects.

- headerGrep:

  regular expression pattern used to recognize header columns found in
  Ingenuity IPA enrichment data.

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

- method:

  integer value indicating the method used to import data from a text
  file, where: `method=1` uses `data.table::read.table()` and the
  `textConnection` argument; `method=2` uses
  [`readr::read_tsv()`](https://readr.tidyverse.org/reference/read_delim.html).
  The motivation to use `data.table::read.table()` is it performed
  better in the presence of UTF-8 characters such as the alpha symbol.

- sheet:

  integer value used only when `ipaFile` is a vector of Excel `.xlsx`
  files, and when the Excel format includes multiple worksheets. This
  value will extract enrichment data only from one worksheet from each
  Excel file.

- sep:

  character string used when `ipaFile` is a vector of text files, to
  split fields into columns. The default will split fields by the tab
  character.

- xlsxMultiSheet:

  logical indicating whether input Excel `.xlsx` files contain multiple
  worksheets.

- useXlsxSheetNames:

  logicl indicating whether to use the Excel worksheet name for each
  imported enrichment table, when importing from `.xlsx` files, and when
  `xlsxMultiSheet=FALSE`. When `xlsxMultiSheet=TRUE` the name is derived
  from the value matched using `ipaNameGrep`, because in this case,
  there are expected to me multiple enrichment tables in one worksheet.

- remove_blank_colnames:

  `logical` indicating whether to drop
  [`colnames()`](https://rdrr.io/r/base/colnames.html) where all values
  are contained in `c(NA, "")`. This option may be preferable
  `remove_blank_colnames=FALSE` when all values in some column like
  `zScore` are `NA`, but you would still like to retain the column for
  consistency with other data. We found that IPA does not report
  `zScore` values when there are only 4 or fewer genes involved in each
  enrichment result.

- convert_ipa_slash:

  `logical` indicating whether to convert IPA gene naming conventions,
  currently some genes are considered one entity in the IPA system, for
  example `"HSPA1A/HSPA1B"` is considered one gene, even though two
  Entrez gene entries `"HSPA1A"` and `"HSPA1B"` can be represented.
  Regardless whether one or both genes are provided to IPA, it considers
  it one entity for the purpose of pathway enrichment hypergeometric
  testing. Unfortunately, the forward slash `"/"` is also used by
  `clusterProfiler` object `enrichResult` as gene delimiter, and is
  hard-coded and cannot be changed. So it will automatically consider
  `"HSPA1A/HSPA1B"` as two genes, causing a mismatch with the IPA
  results. When `convert_ipa_slash=TRUE` by default, it converts the
  forward slash `"/"` to the value of argument `ipa_slash_sep`.

- ipa_slash_sep:

  `character` string used as a delimited when `convert_ipa_slash=TRUE`,
  used to replace genes that contain forward slash `"/"` to use another
  character.

- revert_ipa_xref:

  `logical` indicating whether to revert the IPA gene symbols reported,
  which requires that the IPA data contains a section
  `"Analysis Ready Molecules"`.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` of `data.frame` objects, where each `data.frame` contains
enrichment data for one of the Ingenuity IPA enrichment tests.

## Details

This function parses Ingenuity Pathway Analysis ('IPA') enrichment data
into a `list` of `data.frame` objects for downstream analysis. Each
`data.frame` represents the results of one Ingenuity IPA test, however
not all sections contain gene set enrichment results.

### Batch processing

When importing multiple files, argument `ipaFile` can be a vector of
'.xlsx' or '.txt' files. This workflow also calls
[`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md)
to generate a gene hit matrix, stored as `attr(ipalist, "geneHitIM")` to
use in
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

### IPA Gene Xref Data

By default, the argument `revert_ipa_xref=TRUE` will convert the IPA
gene symbol values back to the original identifier. In most cases\* this
behavior is desirable, with caveats:

- When using platform data which use a non-gene identifier, including
  microarray probesets, or RefSeq transcript ID, or protein "UniProt"
  accession numbers, it is recommended to use `revert_ipa_xref=FALSE`.
  In these cases, the IPA gene symbol is expected to be more
  user-friendly, and therefore more useful.

  - This option would be helpful to view the IPA gene symbols as they
    appear in the IPA report - even if the symbols sometimes do not
    match the input row identifiers.

- It is helpful to use `revert_ipa_xref=TRUE` when the identifier will
  also be used to compare to the source data, for example if trying to
  make an expression heatmap of the genes involved in enrichment
  results.

  - This option would convert IPA gene symbol back to Affymetrix
    probeset ID, for example, if the probeset ID values were used as the
    primary identifier for each measurement. The probeset ID might be
    convenient to align with the input data matrix.

  - The primary reason for this option is when providing gene symbols as
    input to IPA, some will be renamed to IPA preferred gene symbols,
    which would therefore be difficult to match with the gene symbols
    provided to IPA.

In any case, the output `list` should contain an entry "Analysis Ready
Molecules" with the full IPA data table used for the analysis. This
`data.frame` will also contain any statistical columns, if provided to
IPA upfront.

See
[`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md)
for an automated way to create a gene hit matrix from the 'Analysis
Ready Molecules' returned by IPA.

### Motivation

1.  Separate multiple IPA enrichment tables.

2.  Rename colnames to be consistent, compatible with
    [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md).

3.  Revert IPA gene aliases to original user input, default but
    optional.

4.  Generate `geneHitIM` hit matrix when processing multiple files.

### Input format

1.  `ipaFile` can be one or more text `.txt` files, where the text file
    contains all IPA enrichment data in tall format. This format is most
    common.

2.  `ipaFile` can be one or more Excel `.xlsx` files, which contains all
    IPA enrichment data in one tall worksheet tab.

3.  `ipaFile` can be one or more Excel `.xlsx` files, where each type of
    IPA enrichment appears on a separate Excel worksheet tab.

4.  `ipaFile` can be a list of `data.frame` objects. This option is
    intended when the IPA data has already been imported into R as
    separate `data.frame` objects.

### Notes

When using `"Export All"` from 'IPA', the default text format includes
multiple enrichment tables concatenated together in one file. Each
enrichment table contains its own unique column headers, with
descriptive text in the line preceding the column headers. This function
is intended to separate the enrichment tables into a list of
`data.frame` objects, and retain the descriptive text as names of the
list.

### Troubleshooting

- A common error occurs when reverting IPA gene symbols to the original
  user-supplied identifier, by default `revert_ipa_xref=TRUE`. For
  errors during this step, consider `revert_ipa_xref=FALSE` which will
  retain the gene symbol as recognized by IPA. The downside of this
  approach is that it may be more difficult to equate to the input
  identifier. In that case look at the "Analysis Ready Molecules"
  `data.frame` which should contain the user-provided values as "ID";
  the IPA recognized symbol as "Name", and optionally a column "Symbol"
  which is edited by multienrichjam.

## See also

Other jam import functions:
[`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md),
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)

## Examples

``` r
ipaFile <- system.file(package="multienrichjam", "extdata",
   c("Newborns-IPA.txt", "OlderChildren-IPA.txt"));
ipalist <- importIPAenrichment(ipaFile)
```
