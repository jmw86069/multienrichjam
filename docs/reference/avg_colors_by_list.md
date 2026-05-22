# Average colors by list

Average colors by list

## Usage

``` r
avg_colors_by_list(
  x,
  useWeightedHue = TRUE,
  Cmethod = c("mean", "max", "min"),
  Lmethod = c("mean", "max", "min"),
  c_min = 4,
  grey_hue = 359,
  ...
)
```

## Arguments

- x:

  `list` of character vectors.

- useWeightedHue:

  logical indicating whether to weight the hue wheel using
  [`colorjam::h2hw()`](https://jmw86069.github.io/colorjam/reference/h2hw.html)
  and
  [`colorjam::hw2h()`](https://jmw86069.github.io/colorjam/reference/hw2h.html)
  which effectively converts the RGB angles to RYB (red-yellow-blue),
  and therefore makes additive color blending more sensible.
  Specifically, "yellow and blue makes green".

- ...:

  additional arguments are ignored.

## Value

vector of R colors.

## Details

This is a simple wrapper function intended to provide a rapid average
color, when supplied a list of color vectors in hex or R color name
format.

This function simply converts each color to HCL, determines the color
hue angle (from 0 to 360) then calculates the average angular color hue
using
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
then applies that to the maximum C and L values to determine the new
color. It is deliberately intended to ignore muddiness when averaging
multiple colors.

Colors are only modified for elements with 2 or more entries.

This method also only operates on the unique set of colors, so it should
be substantially more efficient on large lists that contain only a few
unique subsets of colors.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md),
[`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md),
[`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md),
[`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md),
[`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
[`find_enrich_colnames()`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md),
[`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md),
[`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md),
[`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md),
[`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md),
[`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md),
[`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md),
[`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md),
[`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md),
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
x <- list(input1=c(red="red", blue="blue"),
   input2=c(blue="blue", gold="gold"),
   input3=c(red="red", yellow="yellow"));
x_avg <- avg_colors_by_list(x, useWeightedHue=TRUE);
jamba::showColors(list(
   input1=c(x[[1]], x_avg[1]),
   input2=c(x[[2]], x_avg[2]),
   input2=c(x[[3]], x_avg[3])),
   main="With weighted hue")


x_avg <- avg_colors_by_list(x, useWeightedHue=FALSE);
jamba::showColors(list(
   input1=c(x[[1]], x_avg[1]),
   input2=c(x[[2]], x_avg[2]),
   input2=c(x[[3]], x_avg[3])),
   main="Without weighted hue")

```
