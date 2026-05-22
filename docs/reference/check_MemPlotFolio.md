# Check MemPlotFolio object

Check whether a MemPlotFolio object is valid.

## Usage

``` r
check_MemPlotFolio(object)
```

## Arguments

- object:

  `MemPlotFolio` object

## Details

It requires:

- slots: 'enrichment_hm', 'gp_hm', 'clusters_mem', 'cnet_clusters',
  'cnet_exemplars', 'cnet_clusters'

General guidance for MemPlotFolio objects:

- The aim is to collect all parameters and associated data together so
  that the various types of plots are all using consistent underlying
  data and settings.

## See also

Other MemPlotFolio:
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`score_gene_path_clusters()`](https://jmw86069.github.io/multienrichjam/reference/score_gene_path_clusters.md)
