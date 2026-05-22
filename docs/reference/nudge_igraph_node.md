# Nudge igraph layout by node

Nudge igraph layout by node

## Usage

``` r
nudge_igraph_node(
  g,
  nodes = NULL,
  x = 0,
  y = 0,
  nodes_xy = NULL,
  use_grep = FALSE,
  aspect = 1,
  debug = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object that contains layout coordinates stored as graph
  attribute `"layout"`, for example `graph_attr(g, "layout")`.

- nodes:

  `character` vector indicating which nodes in `g` should be nudged.

- x, y:

  `numeric` values indicating the amount to move each node defined by
  `nodes`. These values are relative to the x- and y-axis ranges of the
  layout coordinates, and based upon `aspect` below.

- nodes_xy:

  `list` alternative to using arguments `nodes,x,y`. This argument
  assumes a `list` of `numeric` adjustments to `x` and `y`, where the
  `names(nodes_xy)` are node names. For example:

  - `nodes_xy=list(APOE=c(x=0.05, y=-0.01), GAPDH=c(x=0.1, y=0.0)`
    however the `numeric` vector does not need to contain names.

- use_grep:

  `logical` indicating whether to match values in `nodes` to `V(g)$name`
  and `V(g)$label` using
  [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html),
  which follows case-insensitive
  [`grep()`](https://rdrr.io/r/base/grep.html). When `use_grep=FALSE`
  the values in either `V(g)$name` or `V(g)$label` must be identical to
  `nodes`.

- aspect:

  `numeric` indicating the aspect ratio of the output plot. Any value
  other than `aspect=1` will use the observed x-axis and y-axis range of
  the layout coordinates. When `aspect=1` then the highest of x-axis and
  y-axis ranges is used for the relative `x` and `y` adjustment. Note
  that `igraph::plot()` does not maintain aspect ratio 1 by default, but
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  does maintain aspect ratio 1 and is preferred.

- debug:

  `logical` indicating whether to plot the layout before and after
  adjustment, drawing arrows to indicate the movement of particular
  nodes. This plot is very basic, using base R
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and is only
  intended as a quick visual review.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.
