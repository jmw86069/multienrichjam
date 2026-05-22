# Reese 2019 cohort multienrichment data

Data used for multi-enrichment analysis as described in [Reese et al
2019](https://doi.org/10.1016/j.jaci.2018.11.043)
doi:10.1016/j.jaci.2018.11.043, **Epigenome-wide meta-analysis of DNA
methylation and childhood asthma**. It is provided for convenient
testing of multienrichjam.

## Usage

``` r
Memtest

Reese_genes
```

## Format

`Memtest` is a `Mem` multi-enrichment object including two enrichment
tests 'Newborns' and 'OlderChildren'

`Reese_genes` is a `list` with names 'Newborns' and 'OlderChildren',
each with a `character` vector of gene symbols, named by ENTREZID,
representing genes associated with each cohort. For example, the
`names(Reese_genes$Newborns)` may be useful for pathway enrichment
methods that expect Entrez gene ID as input.

## Source

[Reese et al 2019](https://doi.org/10.1016/j.jaci.2018.11.043)
doi:10.1016/j.jaci.2018.11.043, **Epigenome-wide meta-analysis of DNA
methylation and childhood asthma**

## Details

Memtest is a `Mem` object with enrichment results for 'Newborns' and
'OlderChildren'.

A minor change was made to published data:

- Gene symbols were updated to current Entrez gene symbols as of
  November 2025. This process changed 'INADL' to 'PATJ' in the
  enrichment results.

- Two additional genes were updated in 'geneHitIM' for genes not present
  in the IPA pathway results.

- The `Reese_genes` data were also updated for consistency, and to
  ensure that the `names(Reese_genes)` using Entrez ID would be
  consistent with the gene symbols shown. Finally, the genes in
  `Reese_genes` should be usable for other pathway enrichment tests.

## See also

Other Mem:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`Mem-slots`](https://jmw86069.github.io/multienrichjam/reference/Mem-slots.md),
[`check_Mem()`](https://jmw86069.github.io/multienrichjam/reference/check_Mem.md)

## Examples

``` r
# Mem data for Reese 2019 cohorts
data(Memtest)
# enrichment names
names(Memtest)
#> [1] "Newborns"      "OlderChildren"
# genes represented
genes(Memtest)
#>  [1] "ADCY3"   "CASP8"   "CITED2"  "DDIT4"   "IGF1R"   "IL4"     "MS4A2"  
#>  [8] "NAPA"    "NR0B2"   "PATJ"    "POLR2B"  "PPP2CA"  "PRKAR2B" "PRKCH"  
#> [15] "PRKCZ"   "RBL2"    "RPS6KA2" "RPS23"   "RPTOR"   "SLC7A1"  "SMARCD3"
#> [22] "YWHAQ"  
# gene sets, pathways represented
sets(Memtest)
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

# Reese 2019 genes by cohort
data(Reese_genes)
# Number of genes per cohort
lengths(Reese_genes)
#>      Newborns OlderChildren 
#>            41           162 

# first 6 gene symbols for Newborns
head(Reese_genes$Newborns)
#>       6296         87      11174      80216      79572     552889 
#>    "ACSM3"    "ACTN1"  "ADAMTS6"    "ALPK1"  "ATP13A3" "ATXN7L3B" 

# first 6 Entrez gene IDs for OlderChildren
head(names(Reese_genes$OlderChildren))
#> [1] "23457" "11332" "54"    "109"   "285"   "56311"
```
