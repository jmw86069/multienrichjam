# Return the Mem enrichment matrix using column ranks

Return the Mem enrichment matrix using column ranks, applying thresholds
'topEnrichN', 'min_count', and 'p_cutoff'.

## Usage

``` r
enrichIMrank(x, topEnrichN = NULL, p_cutoff = NULL, min_count = NULL, ...)
```

## Arguments

- x:

  `Mem` object

- topEnrichN, p_cutoff, min_count:

  `numeric` values obtained from `thresholds(x)` by default, however
  they can be customized to determine updated ranks using the custom
  thresholds.

- ...:

  additional arguments are ignored.

## Value

`numeric` matrix of `numeric` ranks by column. There should be no ties,
all values should be `integer`.

## Details

Ties are broken first by using `enrichIMgeneCount` then by using the
current row and column order.

## See also

Other jam Mem utilities:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

## Examples

``` r
enrichIMrank(Memtest, min_count=2)
#>                                                                       Newborns
#> Huntington's Disease Signaling                                              NA
#> mTOR Signaling                                                              NA
#> Role Of Nfat In Cardiac Hypertrophy                                         NA
#> RAR Activation                                                              NA
#> Creb Signaling In Neurons                                                   NA
#> Tight Junction Signaling                                                    NA
#> Dopamine-Darpp32 Feedback In Camp Signaling                                 NA
#> Hepatic Cholestasis                                                         NA
#> eNOS Signaling                                                               1
#> Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages       NA
#> Synaptic Long Term Depression                                               NA
#> Fc Epsilon RI Signaling                                                     NA
#> Glioma Signaling                                                            NA
#> HIPPO Signaling                                                             NA
#> Growth Hormone Signaling                                                    NA
#> p70S6K Signaling                                                            NA
#>                                                                       OlderChildren
#> Huntington's Disease Signaling                                                   12
#> mTOR Signaling                                                                    2
#> Role Of Nfat In Cardiac Hypertrophy                                              16
#> RAR Activation                                                                    8
#> Creb Signaling In Neurons                                                        14
#> Tight Junction Signaling                                                         11
#> Dopamine-Darpp32 Feedback In Camp Signaling                                      10
#> Hepatic Cholestasis                                                               5
#> eNOS Signaling                                                                   13
#> Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages            15
#> Synaptic Long Term Depression                                                     9
#> Fc Epsilon RI Signaling                                                           7
#> Glioma Signaling                                                                  6
#> HIPPO Signaling                                                                   4
#> Growth Hormone Signaling                                                          3
#> p70S6K Signaling                                                                  1

enrichIMrank(Memtest, topEnrichN=20)
#>                                                                       Newborns
#> Huntington's Disease Signaling                                              NA
#> mTOR Signaling                                                              NA
#> Role Of Nfat In Cardiac Hypertrophy                                         NA
#> RAR Activation                                                              NA
#> Creb Signaling In Neurons                                                   NA
#> Tight Junction Signaling                                                    NA
#> Dopamine-Darpp32 Feedback In Camp Signaling                                 NA
#> Hepatic Cholestasis                                                         NA
#> eNOS Signaling                                                              NA
#> Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages       NA
#> Synaptic Long Term Depression                                               NA
#> Fc Epsilon RI Signaling                                                     NA
#> Glioma Signaling                                                            NA
#> HIPPO Signaling                                                             NA
#> Growth Hormone Signaling                                                    NA
#> p70S6K Signaling                                                            NA
#>                                                                       OlderChildren
#> Huntington's Disease Signaling                                                   12
#> mTOR Signaling                                                                    2
#> Role Of Nfat In Cardiac Hypertrophy                                              16
#> RAR Activation                                                                    8
#> Creb Signaling In Neurons                                                        14
#> Tight Junction Signaling                                                         11
#> Dopamine-Darpp32 Feedback In Camp Signaling                                      10
#> Hepatic Cholestasis                                                               5
#> eNOS Signaling                                                                   13
#> Production Of Nitric Oxide And Reactive Oxygen Species In Macrophages            15
#> Synaptic Long Term Depression                                                     9
#> Fc Epsilon RI Signaling                                                           7
#> Glioma Signaling                                                                  6
#> HIPPO Signaling                                                                   4
#> Growth Hormone Signaling                                                          3
#> p70S6K Signaling                                                                  1
```
