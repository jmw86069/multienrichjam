# Draw Heatmap with title and subtitle using grid viewports, deprecated

Draw Heatmap with title and subtitle using grid viewports, deprecated

## Usage

``` r
grid_with_title(
  object,
  title = NULL,
  title_fontsize = 18,
  title_just = "centre",
  caption = NULL,
  caption_fontsize = 12,
  caption_just = "left",
  caption_x = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  an object with class `"Heatmap"`, `"gTree"`, or any object where
  `methods::hasMethod("draw", class(object))` is `TRUE`.

- title:

  character string used as title. When `NULL` or `nchar(title)==0` then
  no title is displayed.

- title_fontsize:

  numeric value indicating the font size, with units `"points"`.

- title_just:

  character string indicating the justification (alignment) of the
  title.

- caption, caption_fontsize, caption_just:

  arguments equivalent to the title\_\* arguments.

- caption_x:

  numeric value in grid units specifying where to position the caption
  text. By default when `caption_just` is `"left"` the `caption_x` is
  defined by `grid::unit(0.2, "npc")`, which positions the caption at
  the left side (20% from the left edge) of the plot, with text
  proceeding to the right of that point.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

The byproduct of this function is to draw a grid visualization that
includes a title and subtitle, then the `object` in the center.

## Details

This function is likely to work only with `grid` or `ComplexHeatmap`
input. Other objects that have `draw()` defined may not work properly.
This function is deprecated, as the equivalent is performed by using
`column_title` in
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

Note: This function may be superceded by the ability to include an
overall title when using
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html),
which was not originally possible.

This function is intended to help make it easier to wrap a heatmap from
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
inside the relevant grid viewports in order to display one large title
at the top of the resulting visualization, and optionally one subtitle
at the bottom of the visualization.

The input can be one heatmap (`Heatmap` object), or object of class
`"gTree"`, or any object with a method `draw()` associated with it,
detected by `methods::hasMethod("draw", class(object))`.

A good example of `"gTree"` object is the result of calling `draw()`
inside [`grid::grid.grabExpr()`](https://rdrr.io/r/grid/grid.grab.html),
for example:

`grid::grid.grabExpr(ComplexHeatmap::draw(...))`

The
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
function has extended capability for arranging one or more heatmaps and
associated annotations.

## See also

Other jam deprecated functions:
[`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md),
[`enrichMapJam()`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md),
[`mem_multienrichplot()`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md),
[`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
