# Default igraph parameter values

Default igraph parameter values

## Usage

``` r
default_igraph_values()
```

## Value

`list` of default igraph plotting and data parameters, including
`"plot"`, `"vertex"`, and `"edge"`.

## Details

Note that default values includes some new vertex parameters:

- 'size.scaling': `logical` whether to apply relative size to nodes

- 'relative.size': `numeric` vector length=2, the range of node scaled
  sizes

The 'relative.size' parameter is ignored by
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

Relative node scaling is performed by
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
using a method already implemented in multienrichjam prior to igraph
2.2.0. In future,
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
may implement the equivalent method in igraph 2.2.0+ for consistency.

## See also

Other jam igraph internal functions:
[`get_igraph_arrow_mode()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_arrow_mode.md),
[`jam_igraph_arrows()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph_arrows.md),
[`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)

## Examples

``` r
# summary of default values for graph, vertex, edge
jamba::ssdim(default_igraph_values())
#> $plot
#>         rows     class
#> palette    8 character
#> layout     1  function
#> margin     4   numeric
#> rescale    1   logical
#> asp        1   numeric
#> frame      1   logical
#> main       1  function
#> sub        1 character
#> xlab       1  function
#> ylab       1 character
#> 
#> $vertex
#>                    rows     class
#> color                 1   numeric
#> size                  1   numeric
#> size2                 1   numeric
#> size.scaling          1   logical
#> relative.size         2   numeric
#> label                 1  function
#> label.degree          1   numeric
#> label.color           1 character
#> label.dist            1   numeric
#> label.family          1 character
#> label.font            1   numeric
#> label.cex             1   numeric
#> frame.color           1 character
#> frame.width           1   numeric
#> shape                 1 character
#> pie                   1   numeric
#> pie.color             6      list
#> pie.border            1      list
#> pie.angle             1   numeric
#> pie.density           1   numeric
#> pie.lty               1   numeric
#> pie.lwd               1   numeric
#> coloredrect.lwd       1   numeric
#> coloredrect.border    1 character
#> 
#> $edge
#>              rows     class
#> color           1 character
#> label           1  function
#> lty             1   numeric
#> width           1   numeric
#> loop.angle      1   numeric
#> loopangle2      1   numeric
#> label.family    1 character
#> label.font      1   numeric
#> label.cex       1   numeric
#> label.color     1 character
#> label.x         0      NULL
#> label.y         0      NULL
#> arrow.size      1   numeric
#> arrow.mode      1  function
#> curved          1  function
#> arrow.width     1   numeric
#> 
```
