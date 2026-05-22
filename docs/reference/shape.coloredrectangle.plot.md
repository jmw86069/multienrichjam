# custom igraph vertex shape coloredrectangle

custom igraph vertex shape coloredrectangle

## Usage

``` r
shape.coloredrectangle.plot(coords, v = NULL, params)

shape.coloredrectangle.clip(coords, el, params, end = c("both", "from", "to"))
```

## Arguments

- coords:

  two-column numeric matrix with x- and y-coordinates, respectively.

- v:

  optional ids of the vertices to plot. It should match the number of
  rows in the `coords` argument.

- params:

  a function object that can be called to query vertex/edge/plot
  graphical parameters. The first argument of the function is
  `'vertex'`, `'edge'` or `'plot'` to decide the type of the parameter,
  the second is a character string giving the name of the parameter.

## Value

the plot function returns an invisible `list` of `data.frame` objects
which were used to draw the rectangle objects. However the purpose of
this function is the by-product that it draws rectangles onto an igraph
graph. The clip function returns coordinates corresponding to the outer
node edge based upon argument `end`:

- `end="both"` returns 4 columns, the x,y 'edge from' and x,y 'edge to'
  coordinates, respectively

- `end="from"` returns 2 columns, the x,y 'edge from' coordinates

- `end="to"` returns 2 columns, the x,y 'edge to' coordinates

## Details

This function defines the plotting function for custom igraph vertex
shape coloredrectangle. The coloredrectangle shape is described as:

- Each vertex is drawn as a rectangle, filled with squares which are
  each individually colored.

- The squares are arrayed into a number of columns and rows, which are
  defined for each vertex using `coloredrect.ncol` and
  `coloredrect.nrow`, respectively.

- The vector of colors is arrayed as values of a matrix, therefore
  `coloredrect.byrow` is logical to indicate whether colors should fill
  the matrix by row, similar to how `byrow=TRUE` is used.

- The colors for each vertex are defined by `coloredrect.color`, which
  is a `list` of `character` color vectors. When `coloredrect.color`
  does not exist, values from `pie.color` will be used if they exist.

- Any missing colors are displayed as `NA` which applies no color. For
  example, when a vertex has 3 columns, 2 rows, and only 5 colors, the
  colors are not recycled. Instead, the last color is `NA` to render no
  color in that position.

- Each square may also have an optional border, defined by
  `coloredrect.border`, `coloredrect.lwd`, and `coloredrect.lty`. These
  borders are drawn as inner borders so they do not overlap optional
  frame border.

- Each node may have a frame border, defined by `frame.color`,
  `frame.width`, and `frame.lty`. The frame border is drawn as an outer
  border so it does not overlap inner borders, adding some height and
  width to the final node. The frame border is drawn before the node
  squares are drawn.

- The size of each square inside the rectangle is defined by `size2`,
  such that the rectangle width is `coloredrect.ncol * size2` and the
  rectangle height is `coloredrect.nrow * size2`.

- Node sizes may be adjusted by enabling `equalize_sizes` in one of two
  ways:

  1.  `options(coloredrectangle.equalize_sizes=TRUE)` (priority); or

  2.  `vertex.coloredrect.equalize_sizes=TRUE`, which is equivalent to
      `igraph::V(g)$coloredrect.equalize_sizes=TRUE`. in the latter
      case, only the first value is recognized.

- The behavior of `equalize_sizes` is described below:

  - `equalize_sizes=FALSE` (default), nodes are exactly `size2`
    multiples of `coloredrect.ncol` and `coloredrect.nrow`.

  - `equalize_sizes=TRUE` or `1`, the shortest side of each node is
    scaled to `size2`. This option is useful to ensure nodes are never
    smaller than `size2` width or height, but can be larger.

  - `equalize_sizes=2`, the longest side of each node is scaled to
    `size2`. This option is useful to ensure nodes fit insides a square
    with `size2` sides.

- Nodes are drawn using vectorized processes where possible. However,
  the primary function
  [`graphics::symbols()`](https://rdrr.io/r/graphics/symbols.html) only
  permits vectorized plotting for one `lwd`/`lty` combination at a time,
  so rendering is split into unique combinations of `lwd`/`lty`.
  Vectorized rendering is substantially faster than iterative rendering,
  and the majority of circumstances uses only one `lwd`/`lty`
  combination for all nodes. Line width is not an optimal way to convey
  a quantitative measurement, however it can be useful to highlight
  particular nodes of interest.

- When `coloredrect.ncol` does not exist, the `ncol` will be defined by
  allowing up to two rows by default, enough to accommodate the number
  of colors in `coloredrect.color`.

- When `coloredrect.ncol * coloredrect.nrow` will not fit all colors in
  `coloredrect.color`, the `coloredrect.ncol` will be extended to create
  enough positions to display all colors.

- When `coloredrect.nrow` does not exist, it will use a value based upon
  `coloredrect.ncol` and the length of `coloredrect.color`.

The values `coloredrect.color` and other variables described above refer
to `igraph` vertex attributes, and can be accessed for a given `igraph`
object `g` as follows:

- `igraph::V(g)$coloredrect.color`

- `igraph::vertex_attr(g, "coloredrect.color")`

- or during plotting the value can be defined using the syntax
  `igraph::plot(g, vertex.coloredrect.color=list(...))`

Note that blank positions inside coloredrectangle nodes can be removed
via
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
which also has the effect of modifying the `coloredrect.ncol` and
`coloredrect.nrow`, by applying the appropriate logic.

Todo: The clip function should adjust the node boundary to account for
`frame.color` and `frame.lwd` when present, since the frame is drawn as
an outer border and slightly increases the size of the node.

## See also

[`igraph::shapes()`](https://r.igraph.org/reference/shapes.html)

Other jam igraph shapes:
[`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md),
[`shape.ellipse.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.clip.md),
[`shape.ellipse.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.plot.md),
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)

Other jam igraph shapes:
[`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md),
[`shape.ellipse.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.clip.md),
[`shape.ellipse.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.plot.md),
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)

## Examples

``` r
# prepare example igraph object
am <- matrix(ncol=5, nrow=5,
   data=0,
   dimnames=list(LETTERS[1:5], LETTERS[1:5]))
am[2:5, 1] <- 1;
g1 <- igraph::graph_from_adjacency_matrix(am)
igraph::graph_attr(g1, "layout") <- cbind(x=c(0, 1, 1, -1, -1),
   y=c(0, 1, -0.5, 0.5, -1))
colorset <- c("firebrick3", "gold", "deepskyblue",
   "mediumpurple3", "orchid1")
colorset <- c("firebrick3", "dodgerblue3");
vseq <- seq_len(igraph::vcount(g1));
vsizes <- c(3, 2, 2, 2, 1);
set.seed(1);
igraph::V(g1)$coloredrect.border <- lapply(vseq, function(i){
   sample(colorset,
      replace=TRUE,
      size=vsizes[i])
})
igraph::V(g1)$coloredrect.color <- lapply(vseq, function(i){
   jamba::alpha2col(alpha=0.5,
      igraph::V(g1)$coloredrect.border[[i]])
})
igraph::V(g1)$coloredrect.ncol <- c(2, 2, 2, 1, 1);
igraph::V(g1)$coloredrect.nrow <- c(2, 1, 1, 2, 1);
igraph::V(g1)$coloredrect.lwd <- rep(3, igraph::vcount(g1))
igraph::V(g1)$frame.width <- c(2, 1, 1, 1, 1);
igraph::V(g1)$frame.color <- "black"
igraph::V(g1)$size2 <- 10;
igraph::V(g1)$shape <- "coloredrectangle";

jam_igraph(g1, vertex.label="")
#> ##  (11:03:55) 22May2026:   vertex.size: 
#>      vertex.size1 vertex.size2
#> [1,]    0.5555556    0.5555556
#> [2,]    0.5555556    0.2777778
#> [3,]    0.5555556    0.2777778
#> [4,]    0.2777778    0.5555556
#> [5,]    0.2777778    0.2777778
title(font.main=1, line=1.5, main=paste0(
   "Each square is consistent size by vertex.size2\n"))
title(font.main=1, cex.main=1, line=0.5, main=paste0(
   "'vertex.coloredrect.equalize_sizes=FALSE' or\n",
   "'options(coloredrectangle.equalize_sizes=FALSE'"))


# equalize shortest side to size2
jam_igraph(g1, vertex.label="", vertex.coloredrect.equalize_sizes=TRUE)
#> ##  (11:03:55) 22May2026:   vertex.size: 
#>      vertex.size1 vertex.size2
#> [1,]    0.2777778    0.2777778
#> [2,]    0.5555556    0.2777778
#> [3,]    0.5555556    0.2777778
#> [4,]    0.2777778    0.5555556
#> [5,]    0.2777778    0.2777778
title(font.main=1, line=1.5, main=paste0(
   "The shortest side is fixed by vertex.size2\n"))
title(font.main=1, cex.main=1, line=0.5, main=paste0(
   "'vertex.coloredrect.equalize_sizes=TRUE' or\n",
   "'options(coloredrectangle.equalize_sizes=TRUE'"))


# equalize longest side to size2
jam_igraph(g1, vertex.label="", vertex.coloredrect.equalize_sizes=2)
#> ##  (11:03:56) 22May2026:   vertex.size: 
#>      vertex.size1 vertex.size2
#> [1,]    0.2777778    0.2777778
#> [2,]    0.2777778    0.1388889
#> [3,]    0.2777778    0.1388889
#> [4,]    0.1388889    0.2777778
#> [5,]    0.2777778    0.2777778
title(font.main=1, line=1.5, main=paste0(
   "The longest side is fixed by vertex.size2\n"))
title(font.main=1, cex.main=1, line=0.5, main=paste0(
   "'vertex.coloredrect.equalize_sizes=2' or\n",
   "'options(coloredrectangle.equalize_sizes=2'"))

```
