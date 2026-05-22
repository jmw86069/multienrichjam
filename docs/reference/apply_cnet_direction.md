# Apply Cnet border color by directionality

Apply Cnet border color by directionality

## Usage

``` r
apply_cnet_direction(
  cnet,
  hitim = NULL,
  col = circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c("blue", "grey80",
    "firebrick3")),
  col_l_max = 80,
  hide_solo_pie = TRUE,
  frame_blank = default_igraph_values()$vertex$frame.color,
  frame_lwd_blank = 0.2,
  border_lwd = 2,
  do_reorder = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- cnet:

  `igraph` object with node attribute `"pie.color"` populated as a
  `list` of `character` vectors, named by enrichment. The enrichment
  names should match `colnames(hitim)`.

- hitim:

  `numeric` matrix with values centered at zero for no change, positive
  values `+1` for up-regulation, and negative values `-1` for
  down-regulation. The values are converted to color using color
  function `col`.

- col:

  `function` that takes `numeric` input and assigns a color. The default
  assigns red for positive values, blue for negative values, and white
  for zero, using `colorjam::col_div_xf(1.2)`.

- col_l_max:

  `numeric` maximum HCL Lightness, default 80, for output colors. For
  example, the middle color in `"RdBu_r"` is nearly white, the
  `col_l_max` can be used to apply a darker grey.

- hide_solo_pie:

  `logical` default TRUE, whether a single-color border for a multi-part
  pie node should only apply the color to the overall node with
  'frame.color', and not apply the color to each pie wedge using
  'pie.border'.

  - Default `hide_solo_pie=TRUE`: when all wedges (one or more) have the
    same `pie.border` color, the 'pie.border' is defined as `NA`, and
    'frame.color' is assigned to this color. The effect is to display
    the outline color and not each wedge. Both 'pie.lwd' and
    'frame.width' will be assigned `border_lwd`, and since 'pie.border'
    is `NA` it will not be rendered. Only the 'frame.color' will be
    rendered.

  - When `hide_solo_pie=FALSE` each pie wedge border color is assigned
    to 'pie.border', `frame.color` will be assigned 'frame_blank', and
    'frame.width' will be assigned 'frame_lwd_blank' which is useful for
    displaying a small outer frame for each node.

- frame_blank:

  `character` string to define the color used for 'frame.color' when
  colors are defined in 'pie.border'. The default uses the igraph
  defaults, currently 'black'. In this case, the frame is drawn around
  the inner `pie.border` colors, and only serves to add visual clarity.
  The frame border can be blank `frame_blank="transparent"` or can be a
  thinner line, controlled with `frame_lwd_blank=0.2`.

- frame_lwd_blank:

  `numeric` line width, default 0.2, for nodes that have "blank" frame,
  which also means the 'pie.border' colors must also defined. In this
  case the frame border can be invisible (`frame_lwd_blank=0`) or a very
  thin line (default `frame_lwd_blank=0.2`) to surround the inner
  borders drawn with 'pie.border'.

- border_lwd:

  `numeric` line width, default 2, used when a node matches
  `rownames(hitim)`. When the colors are applied to 'pie.color', the
  border is defined with 'pie.lwd'. When colors are applied to
  'frame.color', the border is defined with 'frame.width'.

- do_reorder:

  `logical` default FALSE, whether to reorder nodes by node attributes
  such as color and border, by calling
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).
  When `do_reorder=TRUE`, other relevant arguments are passed through
  `...` to
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  in particular:

  - `colorV`: to control the expected order of colors. It should be
    supplied if known upfront.

  - `sortAttributes`: to customize the default attribute sort order.

  - `nodeSortBy`: to customize the x-/y- axis arrangement.

  - `orderByAspect`: to enable x-/y- sorting by the aspect ratio of
    nodes in the group. For example, tall-skinny node groups should sort
    by y-axis first, short-wide node groups should sort by x-axis first.

- ...:

  additional arguments are passed to
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  when `do_reorder=TRUE`.

## Details

This function specifically requires the `"pie.color"` node attribute of
the `cnet` object `igraph` is populated with a `list` of colors whose
name is the enrichment, and matches `colnames(hitim)` of the supplied
hit incidence matrix.

The node names `V(cnet)$name` also must match `rownames(hitim)`,
otherwise the `"pie.border"` and `"coloredrect.border"` is assigned its
previous value. In that case `"frame.color"` retains its original value.

This function may be run multiple times with different `hitim`, for
example using `geneIMdirection` and `enrichIMdirection` in two steps.

When there are multiple colors for a given node, they are populated in
node attributes `"pie.border"` in `list` form. For these nodes, the
attribute `"frame.color"` is set to `NULL` so it is not displayed on top
of the pie wedge colors.

When there are no matching `rownames(hitim)`, or all colors are
identical, the `"pie.border"` is populated with `NULL`, so the pie
wedges do not each show a color. Also for single-wedge pie nodes, this
process avoids drawing a small line at the top of each node. Instead the
attribute `"frame.color"` is populated with the one unique border color,
so that only the outer border is colorized.

### Todo

- Consider permitting input to be a simple igraph, then creating
  pie.color, pie.border, pie.border.lwd attributes even if they
  previously did not exist. It would probably optionally also use
  `colorV` to colorize each component of multi-component nodes.

## See also

Other jam cnet igraph functions:
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)
