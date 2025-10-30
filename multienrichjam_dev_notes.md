
# multienrichjam dev notes

* The `.onLoad()` function defined in `zzz.R` adds igraph shapes
'jampie', 'ellipse', and 'coloredrectangle'.

   * When updating the shape functions while using 'devtools' for active
   development, you must also add the shape again to the R session.
   
   ```
   igraph::add_shape("jampie",
      clip=shape.jampie.clip,
      plot=shape.jampie.plot)
   ```

* Adding generic methods for S4 objects

   * Test for existing generic function, the 'getName=TRUE' part tells what
   package provided the generic function
   `isGeneric("score", getName=TRUE)`
   `isGeneric("abbreviate", getName=TRUE)`

   * List all known S3 generic methods:
   `.knownS3Generics` - basic S3 generics
   `utils:::getKnownGenerics()` - more comprehensive S3 generics:
   
      `names(.knownS3Generics)`
      `tools:::.internalGenerics`
      `tools:::.get_S3_primitive_generics()`

* Adding data to include in the R package, and use with R function examples

   * Create the R object in the active R session, for example 'Memtest', then
   `usethis::use_data(Memtest, overwrite=TRUE)`

* Migrating from mem `list` to `Mem` S4 object

* Clean up argument names in `multiEnrichMap()`

   * Change `cutoffRowMinP` to `p_cutoff` 


# Documentation paradigm

* Diataxis: Four types of documentation
https://diataxis.fr/map/

   1. Tutorials      "Can you teach me to...?"

      * learning: a lesson
      * learning a skill, in style of classroom lecture (R vignette)
      * e.g. how to cook

   2. How-to Guides  "How do I...?"

      * goals: series of steps
      * help achieve a particular goal
      * put the skills from the tutorial to work
      * e.g. specific recipe to cook

   3. Reference      "What is...?"

      * information: describe the machinery
      * technical reference used when applying course to practical need (1 and 2)
      * e.g. nutritional information
      
   4. Explanation:   "Why...?"
      
      * understanding: illuminate a topic
      * away from the work, reflect on the practice and knowledge as a whole
      * e.g. culinary social history


# Idea List

* Consider method to add genes into enrichment results `data.frame`

   * Requires: pathway definitions (genes assigned)
   * Requires: gene hits in the enrichment test

* Consider `S4Vectors` object types for S4 object slots

   * Enables things like assignment within a vector or list
   * Some preferred mechanisms for using `metadata()`
   * S4 tips:
   https://bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html#reusing-existing-classes

* Consider tagging with "EDAM" for compatibility with Elixir toolset for example

   * https://blog.bioconductor.org/posts/2025-07-18-edam/
   * https://edamontology.org/page
   * Main goal is to aid discoverability of Bioconductor packages
   * `BiocPkgTools`: manages Bioconductor dependencies, biocViews keywords
   https://github.com/seandavi/BiocPkgTools

* Consider alternative clustering methods to `amap::hcluster()`

   * User can supply `cluster_columns` with custom function, as long as it
   returns 'hclust' compatible output to use with `cutree()`.
   * `cluster::daisy(x, metric='gower')` apparently accepts mixed-type input.

      * It may handle binary 'im' data together with 'enrichIM' -log10(P-value).
      Needs testing, and would need to be told which columns are which type.
      * It scales data to 0-1 using min-max, for directional data probably
      no scaling, and have it "pre-scaled".

   * Directional data is "ternary data" with values: `c(-1, 0, 1)`.
   * "Ternary distance" could be used as a metric for 'geneIMdirection', e.g.:
   `sum(ifelse(x == y, 0, ifelse(x == 0 | y == 0, 1, 2)))` for each pair,
   (effectively `abs(x - y)` but needs to be matrix math).
   Then: `hclust(as.dist(ternary_distance))`
   * Cosine similarity emphasizes directionality more than magnitude.
   * Create two distance matrices, combine then, then run `hclust()`.
   * Consider Multiple Correspondence Analysis (MCA) FactoMineR

# Bioconductor compatibility check

* Bioconductor 'GSEABase'

   * https://bioconductor.org/packages/3.22/bioc/html/GSEABase.html
   * `GeneSet` and `GeneSetCollection` objects

* 'categoryCompare'

   * uses GSEABase::GeneSet
   * https://bioconductor.org/packages/3.22/bioc/vignettes/categoryCompare/inst/doc/categoryCompare_vignette.html
   * Vignette has multiple enrichments as 'enrichLists' that are used downstream

* 'EnrichmentBrowser'

   * wrapper front-end to numerous enrichment approaches
   * https://bioconductor.org/packages/3.22/bioc/vignettes/EnrichmentBrowser/inst/doc/EnrichmentBrowser.html
   * Great discussion/review of points regarding enrichment backgrounds,
   "competitive" versus "self-contained".
   https://bioconductor.org/packages/3.22/bioc/vignettes/EnrichmentBrowser/inst/doc/EnrichmentBrowser.html#125_Underlying_null:_competitive_vs_self-contained
   * great examples for obtaining gene sets from KEGG, Enrichr, others

* 'zenith' - also uses EnrichmentBrowser

   * great examples for obtaining gene sets from MSigDB, GO, etc
   * `zenith_gsa()` returns `data.frame` suitable for input
   * also returns t-statistic with + and - values
   (directionality, or under/over-representation?)
   * uses `limma::camera()` which employs limma-DREAM analysis, accounting
   for correlation within gene expression traits. Need to learn more.

* 'gage'

   * Clunky but might be useful, produces 'greater' and 'less' `data.frame`
   results, however the genes are not included in each row by default.
   * `gageComp()` produces a simplified directional Venn, red/green text labels,
   and no indication of discordance, it only shows concordance.

*'GOstats'

   * original hypergeometric Bioconductor R package - the o.g.
   * returns `GOHyperGResult` object, data.frame-like
   

* Others that produce `data.frame` that should "just work":

   * 'roastgsa' produces data.frame:
   https://bioconductor.org/packages/3.22/bioc/vignettes/roastgsa/inst/doc/roastgsaExample_RNAseq.html

* Bioconductor 'sparrow'

   * `SparrowResult` object
   * https://bioconductor.org/packages/3.22/bioc/vignettes/sparrow/inst/doc/sparrow.html
   * https://bioconductor.org/packages/3.22/bioc/html/sparrow.html

# For reference

Wow, an example of creating a package using a function from the R environment

```
# custom function
say_hello = function(name = "my little friend", punct = ".")
{
  cat( paste0("Hello, ", name, punct), "\n", sep="")
}

# now add to R package
library(devtools)
create("say")
wd("say")
dump("say_hello", file = "R/say_hello.R")
rm(say_hello)  ## remove local copy
load_all()     ## load the package
say_hello()    ## this is the package-version
```
