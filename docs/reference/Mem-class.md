# Mem S4 class, accessors, getters, and setters

Mem class containing results from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
current class version "1.0.0".

## Usage

``` r
# S4 method for class 'Mem,ANY,ANY,ANY'
x[i, j, k, drop = FALSE]

# S4 method for class 'Mem'
show(object)

Mem_to_list(x, ...)

# S4 method for class 'Mem'
as.list(x)

# S4 method for class 'Mem'
names(x)

# S4 method for class 'Mem'
enrichments(x)

# S4 method for class 'Mem'
enrichments(x) <- value

# S4 method for class 'Mem'
sets(x, ...)

# S4 method for class 'Mem'
sets(x) <- value

# S4 method for class 'Mem'
genes(x, ...)

# S4 method for class 'Mem'
genes(x) <- value

# S4 method for class 'Mem'
geneIM(x, ...)

# S4 method for class 'Mem'
geneIMdirection(x, ...)

# S4 method for class 'Mem'
geneIMdirection(x) <- value

# S4 method for class 'Mem'
geneIMcolors(x, ...)

# S4 method for class 'Mem'
geneIMcolors(x) <- value

# S4 method for class 'Mem'
enrichList(x, ...)

# S4 method for class 'Mem'
enrichIM(x, ...)

# S4 method for class 'Mem'
enrichIMcolors(x, ...)

# S4 method for class 'Mem'
enrichIMcolors(x) <- value

# S4 method for class 'Mem'
enrichIMdirection(x, ...)

# S4 method for class 'Mem'
enrichIMdirection(x) <- value

# S4 method for class 'Mem'
enrichIMgeneCount(x, ...)

# S4 method for class 'Mem'
memIM(x, ...)

# S4 method for class 'Mem'
geneHitIM(x, ...)

# S4 method for class 'Mem'
geneHitIM(x) <- value

# S4 method for class 'Mem'
geneHitList(x, ...)

# S4 method for class 'Mem'
geneHitList(x) <- value

# S4 method for class 'Mem'
headers(x)

# S4 method for class 'Mem'
colorV(x, ...)

# S4 method for class 'Mem'
colorV(x) <- value

# S4 method for class 'Mem'
thresholds(x, ...)

# S4 method for class 'Mem'
thresholds(x) <- value

# S4 method for class 'Mem'
dim(x)

# S4 method for class 'Mem'
dimnames(x)

list_to_Mem(mem)

# S4 method for class 'Mem'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'Mem'
geneInCategory(x)

# S4 method for class 'Mem'
setsByGene(x)

# S4 method for class 'Mem'
EnrichmentMap(
  x,
  do_plot,
  legend_x = "bottomleft",
  legend_y = NULL,
  params = list(),
  ...
)
```

## Arguments

- x:

  `Mem` object

- i:

  `character` or `integer` with gene names

- j:

  `character` or `integer` with pathway names

- drop:

  `logical` always set to FALSE for Mem objects.

- ...:

  additional arguments are ignored

- value:

  `character` vector of names to assign to pathway gene sets

- mem:

  `list` output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

## Value

`Mem_to_list()` returns a `list` suitable for other mem functions.

`as.list(Mem)` returns a `list` suitable for other mem functions.

[`names()`](https://rdrr.io/r/base/names.html) returns the `character`
vector of enrichment names

`enrichments()` returns the `character` vector of enrichment names

`sets()` returns the `character` vector of pathway gene sets

`list_to_Mem()` returns a `Mem` S4 object, from 'list' or 'Mem' input.

`updateObject()` returns a `Mem` S4 object with the current class
version, and current valid class definition. It should be able to update
any serialized 'Mem' object from previous class versions.

`geneInCategory()` returns a `list` named by pathway, containing
`character` vectors with genes in each pathway. For the reciprocal, see
`setsByGene()`.

`setsByGene()` returns a `list` named by gene, containing `character`
vectors with associated gene sets. For the reciprocal, see
`geneInCategory()`.

## Functions

- `x[i`: Subset a Mem object, `Mem[i, j, k]`, where `i` are genes, `j`
  are pathways, and `k` are enrichments.

- `show(Mem)`: Show summary of a Mem object, dimensions defined by
  genes, sets, enrichments, showing up to 5 entries of each. Also
  includes a list of analysis parameters: topN, and thresholds.

- `Mem_to_list()`: Coerce S4 `Mem` to legacy `list` mem format

- `as.list(Mem)`: Coerce S4 `Mem` to legacy `list` mem format

- `names(Mem)`: Names for the enrichment tests in a Mem object

- `enrichments(Mem)`: Names for the enrichment tests in a Mem object

- `enrichments(Mem) <- value`: Assign names for the enrichment tests in
  a Mem object

- `sets(Mem)`: List pathway gene set names in a Mem object

- `sets(Mem) <- value`: Assign pathway gene set names to a Mem object

- `genes(Mem)`: List genes represented

- `genes(Mem) <- value`: Assign gene names in a Mem object

- `geneIM(Mem)`: The gene-enrichment matrix of genes represented in
  enrichment tests.

- `geneIMdirection(Mem)`: The matrix of genes tested versus enrichment
  tests, with directionality.

- `geneIMdirection(Mem) <- value`: Assign directional matrix to the
  gene-enrichment matrix

- `geneIMcolors(Mem)`: The matrix of colors indicating genes in each
  enrichment test

- `geneIMcolors(Mem) <- value`: Assign a color matrix to the
  gene-enrichment matrix

- `enrichList(Mem)`: The list of enrichResult data in an Mem object

- `enrichIM(Mem)`: The pathway/P-value matrix

- `enrichIMcolors(Mem)`: The pathway-enrichment directional color matrix

- `enrichIMcolors(Mem) <- value`: Assign the pathway-enrichment
  directional color matrix

- `enrichIMdirection(Mem)`: The pathway-enrichment directional score
  matrix

- `enrichIMdirection(Mem) <- value`: Assign the pathway-enrichment
  directional score matrix

- `enrichIMgeneCount(Mem)`: Pathway-enrichment gene count matrix, genes
  involved in enrichment of each pathway, for each enrichment test.

- `memIM(Mem)`: The gene-pathway incidence matrix

- `geneHitIM(Mem)`: The matrix of genes tested for enrichment, including
  genes not associated with enrichment results.

- `geneHitIM(Mem) <- value`: Assign the matrix of genes tested in a Mem
  object

- `geneHitList(Mem)`: The list of genes tested for enrichment, including
  genes not associated with enrichment results.

- `geneHitList(Mem) <- value`: Assign the list of genes tested for
  enrichment, including genes not associated with enrichment results.

- `headers(Mem)`: Enrichment data column headers associated to
  enrichResult data in a Mem object

- `colorV(Mem)`: List colors assigned to represented

- `colorV(Mem) <- value`: Assign colors to enrichments in a Mem object.
  If supplied 'value' matches the number of 'enrichments(x)', colors are
  assigned in the order supplied. Otherwise, 'names(value)' must be
  present in 'enrichments(x)', and values are assigned in matching
  order.

- `thresholds(Mem)`: List thresholds defined

- `thresholds(Mem) <- value`: thresholds defined in a Mem object. Note
  that thresholds do not trigger any other updates to the Mem object.

- `dim(Mem)`: dimensions in order of genes, sets, and enrichments

- `dimnames(Mem)`: dimension names for each genes, sets, and enrichments

- `list_to_Mem()`: Convert legacy `list` mem format to S4 `Mem`

- `updateObject(Mem)`: Confirm Mem object meets current class
  definition.

- `geneInCategory(Mem)`: Genes in each category (pathway, gene set)

- `setsByGene(Mem)`: Pathway gene sets associated with each gene

- `EnrichmentMap(Mem)`: EnrichmentMap `igraph` network to connect
  pathway gene sets based upon Jaccard overlap between each pathway
  pair. Note '...' arguments are passed to
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  and
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  when `do_plot=TRUE`. Argument `'params'` is a `list` with additional
  arguments: 'repulse', 'width', 'group', 'mark.expand', 'do_legend'.

## See also

Other Mem:
[`Mem-data`](https://jmw86069.github.io/multienrichjam/reference/Mem-data.md),
[`Mem-slots`](https://jmw86069.github.io/multienrichjam/reference/Mem-slots.md),
[`check_Mem()`](https://jmw86069.github.io/multienrichjam/reference/check_Mem.md)

Other jam Mem utilities:
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

## Examples

``` r
# examples for Mem_to_list and as(Mem, "list")
data(Memtest)
mem1 <- Mem_to_list(Memtest)
jamba::sdim(mem1)
#>                   rows cols        class class_v2
#> enrichList           2              list     <NA>
#> enrichLabels         2         character     <NA>
#> colorV               2         character     <NA>
#> geneHitList          2              list     <NA>
#> geneHitIM          203    2       matrix    array
#> memIM               22   16       matrix    array
#> geneIM              22    2       matrix    array
#> enrichIM            16    2       matrix    array
#> multiEnrichDF       16   12   data.frame     <NA>
#> multiEnrichResult   16   14 enrichResult     <NA>
#> thresholds           6              list     <NA>
#> colnames             9              list     <NA>
#> enrichIMcolors      16    2       matrix    array
#> enrichIMdirection   16    2       matrix    array
#> enrichIMgeneCount   16    2       matrix    array
#> geneIMcolors        22    2       matrix    array
#> geneIMdirection     22    2       matrix    array

# list_to_Mem examples
data(Memtest)
mem <- as(Memtest, "list")
Mem <- list_to_Mem(mem)

setsByGene(Memtest)
#> $ADCY3
#> [1] "Role Of Nfat In Cardiac Hypertrophy"        
#> [2] "RAR Activation"                             
#> [3] "Creb Signaling In Neurons"                  
#> [4] "Dopamine-Darpp32 Feedback In Camp Signaling"
#> [5] "Hepatic Cholestasis"                        
#> [6] "eNOS Signaling"                             
#> 
#> $CASP8
#> [1] "Huntington's Disease Signaling" "eNOS Signaling"                
#> 
#> $CITED2
#> [1] "RAR Activation"
#> 
#> $DDIT4
#> [1] "mTOR Signaling"
#> 
#> $IGF1R
#> [1] "Huntington's Disease Signaling"      "Role Of Nfat In Cardiac Hypertrophy"
#> [3] "Synaptic Long Term Depression"       "Glioma Signaling"                   
#> [5] "Growth Hormone Signaling"           
#> 
#> $IL4
#> [1] "Hepatic Cholestasis"                                                  
#> [2] "Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages"
#> [3] "Fc Epsilon RI Signaling"                                              
#> [4] "p70S6K Signaling"                                                     
#> 
#> $MS4A2
#> [1] "Fc Epsilon RI Signaling"
#> 
#> $NAPA
#> [1] "Huntington's Disease Signaling" "Tight Junction Signaling"      
#> 
#> $NR0B2
#> [1] "Hepatic Cholestasis"
#> 
#> $PATJ
#> [1] "Tight Junction Signaling" "HIPPO Signaling"         
#> 
#> $POLR2B
#> [1] "Huntington's Disease Signaling" "Creb Signaling In Neurons"     
#> 
#> $PPP2CA
#> [1] "mTOR Signaling"                                                       
#> [2] "Tight Junction Signaling"                                             
#> [3] "Dopamine-Darpp32 Feedback In Camp Signaling"                          
#> [4] "Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages"
#> [5] "Synaptic Long Term Depression"                                        
#> [6] "HIPPO Signaling"                                                      
#> [7] "p70S6K Signaling"                                                     
#> 
#> $PRKAR2B
#> [1] "Role Of Nfat In Cardiac Hypertrophy"        
#> [2] "RAR Activation"                             
#> [3] "Creb Signaling In Neurons"                  
#> [4] "Tight Junction Signaling"                   
#> [5] "Dopamine-Darpp32 Feedback In Camp Signaling"
#> [6] "Hepatic Cholestasis"                        
#> [7] "eNOS Signaling"                             
#> 
#> $PRKCH
#>  [1] "Huntington's Disease Signaling"                                       
#>  [2] "mTOR Signaling"                                                       
#>  [3] "Role Of Nfat In Cardiac Hypertrophy"                                  
#>  [4] "RAR Activation"                                                       
#>  [5] "Creb Signaling In Neurons"                                            
#>  [6] "Dopamine-Darpp32 Feedback In Camp Signaling"                          
#>  [7] "Hepatic Cholestasis"                                                  
#>  [8] "eNOS Signaling"                                                       
#>  [9] "Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages"
#> [10] "Synaptic Long Term Depression"                                        
#> [11] "Fc Epsilon RI Signaling"                                              
#> [12] "Glioma Signaling"                                                     
#> [13] "Growth Hormone Signaling"                                             
#> [14] "p70S6K Signaling"                                                     
#> 
#> $PRKCZ
#>  [1] "Huntington's Disease Signaling"                                       
#>  [2] "mTOR Signaling"                                                       
#>  [3] "Role Of Nfat In Cardiac Hypertrophy"                                  
#>  [4] "RAR Activation"                                                       
#>  [5] "Creb Signaling In Neurons"                                            
#>  [6] "Tight Junction Signaling"                                             
#>  [7] "Dopamine-Darpp32 Feedback In Camp Signaling"                          
#>  [8] "Hepatic Cholestasis"                                                  
#>  [9] "eNOS Signaling"                                                       
#> [10] "Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages"
#> [11] "Synaptic Long Term Depression"                                        
#> [12] "Fc Epsilon RI Signaling"                                              
#> [13] "Glioma Signaling"                                                     
#> [14] "HIPPO Signaling"                                                      
#> [15] "Growth Hormone Signaling"                                             
#> [16] "p70S6K Signaling"                                                     
#> 
#> $RBL2
#> [1] "Glioma Signaling"
#> 
#> $RPS6KA2
#> [1] "mTOR Signaling"           "Growth Hormone Signaling"
#> 
#> $RPS23
#> [1] "mTOR Signaling"
#> 
#> $RPTOR
#> [1] "mTOR Signaling"
#> 
#> $SLC7A1
#> [1] "eNOS Signaling"
#> 
#> $SMARCD3
#> [1] "RAR Activation"
#> 
#> $YWHAQ
#> [1] "HIPPO Signaling"  "p70S6K Signaling"
#> 
```
