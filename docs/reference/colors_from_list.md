# Ordered colors from a list of color vectors

Ordered colors from a list of color vectors

## Usage

``` r
colors_from_list(x, return_type = c("colors", "order"), verbose = FALSE, ...)
```

## Arguments

- x:

  list of character vectors that contain valid R colors.

## Value

character vector of unique colors in `x`

## Details

This function takes a list of colors and returns the unique order of
colors based upon the order in vectors of the list. It is mainly
intended to be called by
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
however the function is useful for inferring the proper order of unique
colors from a list of various subsets of colors.

The basic assumption is that there exists one true order of unique
colors, and that each vector in the list contains a subset of those
colors which is consistent with this true order of colors.

The function uses only vectors that contain two or more colors, and
therefore requires that all unique colors are present in the subset of
vectos in the list where length \>= 2. It then uses vectors with two or
more colors, calculates the average observed rank for each color, then
uses that average rank to define the overall color order.

If not all unique colors are present in vectors with two or more colors,
the fallback sort uses
[`colorjam::sort_colors()`](https://jmw86069.github.io/colorjam/reference/sort_colors.html).

## See also

Other jam list functions:
[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md),
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md),
[`list2concordance()`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md),
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md),
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)
