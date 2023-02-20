# TODO

## 02feb2023

* Low priority visual enhancement, color Cnet edges by Set.

   * Consider coloring Cnet edges by categorical color, by Set
   It may help clarify node groupings, while also reinforcing
   which edges connect to particular nodesets.
   Unclear whether added color would be confusing or beneficial.
   * Referring to Neely et al paper in ACR Open Rheumatology:
   "Gene Expression Meta-Analysis Reveals Concordance in Gene Activation,
   Pathway and Cell-Type Enrichment in Dermatomyositis (DM) Target Tissues"
   * They show Cnet-type plot in Figure 4A, where connections from each set
   were categorically colored

## 26jan2023

* Low priority: It may be useful to create vectorized functions:

   * `polygon()`
   
      * enable multiple line widths for multiple polygons
      split by `NA` coordinates.
      * enable optional inner and outer borders with varying widths
   
   * `text()`: enable multiple family, srt
   * `lines()`: enable multiple col, lty, lwd for split lines, similar
   to `segments()` except enabled for multiple lines split by `NA` coordinates.

## 23jan203

* `jam_igraph()`

   * Return the input `igraph` object with all relevant object attributes
   updated to reflect the plot parameters.
   The returned object could be plotted directly without any customization.
   * Consider storing/using edge coordinates inside the `igraph` object.
   However, whenever layout is re-calculated, edge coordinates would likely
   become invalid. It would be tricky to handle.

* edge bundling and edge clipping integration

   * DONE: Most scenarios described below.
   * Currently, linear edges are clipped at connected node boundaries.
   * When edges are slightly curved, the start and end positions are reasonable.
   * When edges are bundled, and especially when nodes are relatively large,
   the edge curves through the node in a different direction than a linear edge.
   * This situation is not a visualization problem when:
   
      1. nodes are not filled with transparent color, and
      2. edges are not drawn with arrows.
   
   * The situation is only a visible problem when:
   
      1. edges have arrows that would now be partly or fully covered by the
      node, or
      2. nodes that are filled with partial or fully transparent color,
      thus showing the edge underneath.

   * Proper edge clipping would probably be done by calculating the edge
   from the node actual center point, then clipping edge where it exits the
   node shape border.
   * Implementation may benefit from storing the edge coordinates as
   an edge attribute, to be used by the clipping function when present.
   Absence of custom edge coordinates would cause the clipping function to
   use linear edge coordinates.
   
      * The plot function could also use stored edge coordinates as opposed
      to calling edge bundling function; alternatively the edge bundling
      function could simply re-use existing edge coordinates as well.


## 18jan2023

* `jam_igraph()`, proposed drop-in replacement for `igraph:::plot.igraph()`:

   * FIXED: This function does not seem to handle edge arrows,
   nor does it shorten edges based upon node sizes.
   * bonus points for node/edge legend functions
   * When layout is not defined, the xlim/ylim values sometimes do not
   match the dynamic layout calculated on the fly.

* igraph shape "ellipse"

   * DONE: it should have a proper "clip" function, in order for edge arrows
   to appear at the border of each node

* igraph "coloredrectangle" shape

   * DONE: The coloredrect.border should also be capable of adjacent lines that
   do not overlap.

## 30nov2022

* `reorder_igraph_nodes()`

   * DONE: to be fancy, it should also propagate changes to `label.dist`
   and `label.degree`, as created by `spread_igraph_labels()`,
   since switching coordinates for two nodes should also switch
   the `label.degree` and `label.dist` associated with those
   coordinates.

* `jam_igraph()`

   * It shows some lag before plotting all nodes vectorized, check
   if the edge bundling is the slow step, and optimize as needed.
   UPDATE: Yes the edge bundling step is introducing lag; or the
   resulting curved lines are being plotted slowly?
   Another good reason to store edge bundle coordinates in the `igraph`.

* `mem_legend()`

   * consider using `ComplexHeatmap::Legend()` for consistency with
   future Cnet-Heatmap usage, and because that mechanism is really nice.

* Refactor `multiEnrichMap()`

   * likely create new function `multienrich()` so `multiEnrichMap()` can
   be deprecated and remain for backward compatibility.
   * Avoid pre-calculating Cnet and Emap graphs.
   Counterpoint: It could call `mem_plot_folio()` with defaults, and
   store results back into the `multienrichResult`.
   * Consider `multienrichResult` object type?
   
      * could hold what is now a `list` object, with proper slot names.
      * custom `print()`/`summary()` function to display summary info
      about number of genes, pathways, etc.
   
   * Consider `multienrichPlot` object type?
   No, the decision is to re-use `multienrichResult`.
   
      * Benefit is simplicity: this object feeds all downstream plots.
      * The negative is that it requires saving a separate
      `multienrichResult` with alternative clustering if necessary.
      * Main goal is to define the gene-pathway clustering,
      then re-use clusters in subsequent steps without having to repeat
      the clustering the same way. If one uses custom clustering
      in `mem_plot_folio()` to produce gene-pathway heatmap, they have
      to use the same arguments in subsequent steps otherwise the
      clusters will differ, and it is not clearly indicated to be
      a problem.
   
   * `mem_plot_folio()` may likely return the `multienrichResult`
   after updating internally stored data. It should have option to
   use existing results.

* Write a more directed vignette showing at least two common use cases:

   1. Enrichment using `clusterProfiler::enricher()`, then multi-enrichment.
   2. Enrichment using Qiagen Ingenuity IPA (outside of R) then importing
   the files produced using `"Export All"` from within the Ingenuity IPA app.
   3. Enrichment using some other external (non-R) tool, for example DAVID.

* Write a vignette focused solely on Cnet plot custom layout options.
Background:

   * Very common workflow results in Cnet plot to summarize the findings.
   * This Cnet plot has often been included as a figure or supplementary
   figure for a paper, therefore it requires manual adjustments to
   increase legibility and clarity of the figure.
   * Adjust individual nodes: `nudge_igraph_node()`. Note this function
   can be called on a `list` of nodes using `nodes_xy`, or using vectors
   with `nodes`, `x`, `y`.
   * Adjust sets of nodes:
      * `adjust_cnet_nodeset()` - Usually for sub-clusters of gene nodes,
      move the whole group, adjust intra-group node spacing,
      or rotate the group around the group center.
      * `apply_nodeset_spacing()` - Apply a minimum spacing between
      nodes, for each nodeset.
      * `adjust_cnet_set_relayout_gene()` - Adjust "Set" nodes manually,
      then fix them in place but allow "Gene" nodes to move during re-layout,
      usually with `relayout_with_qfr()`
      
   * Adjust all nodes:
   
      * `layout_with_qfr()`, `layout_with_qfrf()`, `relayout_with_qfr()` -
      All are wrappers to `qgraph::qgraph.layout.fruchtermanreingold()`
      with convenient defaults. `layout_with_qfr()` returns coordinates;
      `layout_with_qfrf()` returns a layout function with user-defined
      `repulse` argument; `relayout_with_qfr()` updates layout in-place
      for an `igraph` object, storing in `graph_attr(g, "layout")`.
      * `rotate_igraph_layout()` - rotate the layout coordinates by some
      user-defined degree angle.
      * `spread_igraph_labels()` - positions node labels radially around
      nodes, based upon the average incoming edge angle.
      * `reorder_igraph_nodes()` - within each nodeset, reposition nodes
      in order of node color, border color, then label.
      A nodeset is defined as a set of nodes that all share the same
      connections, which is mostly only useful for bipartite graphs
      such as Cnet plots. This function performs the node re-ordering
      for all nodesets, across the whole `igraph` object.
   
   * plot with `jam_igraph()`
   
      * Vectorized plotting when multiple shapes are used, otherwise
      `igraph::plot()` uses a `for()` to iterate each node.
      * Improved `pie` rendering, also vectorized.
      * Convenience methods to adjust node size, node label font size,
      node label distance
      * Optional shadow text node labels
      * Maintain aspect ratio = 1, so nodes are symmetrically spaced along
      each axis (defined by the layout algorithm used.)
      * optionally render node groups using `mark.groups`
      * Call edge bundling, especially useful for bipartite graphs such
      as Cnet plots.
      * Optionally draw background grid with percent layout units.

* add dev functions `layout_cnet()` and related iterative layout functions.


## 23nov2022

* mechanism to store edge coordinates in `igraph` object

   * to my knowledge, this functionality does not exist in `igraph`,
   nor is it represented in `ggraph` or `tidygraph` objects.
   However `tidygraph` may have capability to supply specific edge
   coordinates if they exist, so it might be the closest to
   implementing this feature.
   * `igraph::plot()` is not equipped to use edge coordinates
   * `ggraph` is not equipped to use edge coordinates, it only creates
   edge coordinates based upon the edge `geom_` being used.
   * The driving use case is to define edge coordinates to handle:
   
      1. edge bundling, a procedure that could be done dynamically,
      but is computationally expensive for large numbers of nodes;
      2. custom edge pathing, specifically for pathway schematics,
      such as those generated when using GraphViz DOT format. While the
      DOT format generates and stores edge coordinates, I could not
      find examples in R that import a fully-described DOT file (with
      edge coordinates embedded) that also imported edge coordinates.
   
   * Reasons to want edge coordinates upfront:

      * pre-calculate edge bundling, saving time during rendering
      * allow custom definition of edges, for example in a pathway
      schematic where edges are positioned to avoid overlaps.
      * calculate better label placement by using the angle of incoming
      edges, not limited to the linear vector from node1 to node2.

   * Issues raised when storing edge coordinates:
   
      * Obviously, whenever the node layout coordinates are changed,
      the edges also need to be changed.
      * `rotate_igraph_layout()` could rotate nodes and all edges together.
      * `nudge_igraph_node()` and `adjust_cnet_nodeset()` must decide
      whether to adjust edges by simple scaling, or simply invalidate
      all edges, then force the user to re-calculate edge coordinates.
   
   * useful helper functions
   
      * `validate_edge_coords()` - Test whether edge coordinates match
      node coordinates. If not, then delete or replace edges with linear
      equivalent.
      * `adjust_edge()` - Wrapper function to adjust the curvature,
      placement, rotation, of edges. Could be called when rotating node
      layout, to rotate edges accordingly.
      * `bundle_edges()` - Wrapper to apply bundling to one set of edges
      together as a group. Optionally define specific coordinate(s)
      through which the bundling loess curve is routed.
      * fancy options like routing edges with preference for vertical/horizontal
      pathing, with slight curvature at right angle turns. Often used in
      schematic diagrams.
      * fancy "subway" style options, where bundled edges are allowed to
      remain visible adjacent to other edges along their path.
      Edges could "snap" to nearby edges in the same bundle.

## 17nov2022

* `get_cnet_nodeset()`

   * This function is called several times by internal functions,
   and could therefore be much faster than currently.
   * Refactor by using `igraph::as_adjacency_matrix()`.
   Subset rows for `nodeType="Gene"` and columns for `nodeType="Set"`.
   Then should be able to convert rapidly by `venndir::im2list()`
   or some `pasteByRow()` magic, produce nodeset per node (row).
   Finally, split node names by nodeset.

* node adjustment to prevent label overlaps

   * Idea: "stretch" out nodeset nodes "to the right", which fixes the
   left edge of the nodeset, then expands the node spacing outward as
   fraction of current range of nodeset nodes. For example,
   expand to the right by 10%, to improve side-to-side spacing, since
   label overlaps typically occur with nodes at the same y-level.
   * Stretch nodes "to the left" would work the same, but fixes the right
   x-coordinate range, and expands the left x-coordinate range.
   * Stretch nodes "to the top", and "to the bottom" work similarly.
   * Future idea to reduce node label overlap is to treat it like
   biscuit dough under a rolling pin. Stretch subset of nodes
   in each direction until the labels no longer overlap. The trick is
   to stretch nodes away from other nodesets, so it does not cause new
   overlaps with other nodes.

* `adjust_cnet_nodeset()`

   * Consider option to restrict "expand" to x- or y-axis expansion.
   Basic idea is to limit expansion to "widen" the node spacing,
   or to make node spacing "taller". The "widen" option is helpful
   to reduce label overlaps for nodes directly beside each other.

* `spread_igraph_labels()`

   * ideal case: somehow take into account the edge bundling to calculate
   input angle to each node, rather than straight vector from node to node.
   * `node_groups` - spread labels relative to node group centroid, so labels
   in this cluster of nodes will be spaced out from each other.
   Bonus points for taking into account the overall average input angle to
   nodes in each group, and applying a fraction of that offset along with the
   node-to-group offset. For example, for a node group in the top right,
   they generally point to the top-right, but are fanned out slightly so the
   bottom-left-most node is not fanned out to bottom-left, but maybe center,
   or only slightly bottom-left of the node.

* new layout functions specific for Cnet plots

   * `iterate_qfr_layout()` - R code version of
   `qgraph::qgraph.layout.fruchtermanreingold()` with custom addition of
   node "shells", and option to call `iterate_node_group_distance()`
   * `iterate_node_group_distance()` - R layout intended only to enforce
   separation across node groups (defined by `get_cnet_nodeset()`) so there
   is additional space between nodesets in the layout.
   * `layout_cnet()` - wrapper function that calls rounds of layouts.
   This series of steps is currently the best default layout to generate
   the most readable Cnet plot possible.
   
      1. initial node placement - qfr layout without node/group distance
      2. node spacing - qfr with node distance
      3. node/group spacing - qfr with node and group distance
      4a. imposed nodeset percent spacing
      4. group spacing - "fix" node layout per nodeset, then layout nodesets
      5. detect best rotation to place first "Set" node top-left, then rotate
      6. spread labels, re-order nodes by color/border.
   
   * Other ideas:
   
      1. Try the new `bubble_force()` calculation for minimum distance,
      instead of the linear desired distance force.
      2. Consider constraining all Gene nodes, then allow only Set nodes
      to move. Goal is to prevent Set nodes from being embedded inside
      Gene nodes due to overall net forces.
      3. Optional fixed coordinate range, to prevent layouts from becoming
      infinitely large, thus endless cycle of imposing minimum precent spacing,
      making the range larger, thus making the spacing smaller, etc. etc.
      It means when a force would push a node/group out of bounds, it gets
      stopped at the boundary/boundaries. Hopefully because forces are applied
      to both node1 and node2, when node1 cannot move, node2 should still move
      albeit at half speed.
      4. Really pie in the sky: Use the node size and shape to calculate
      distance between nodes, nearest polygon distance for example. Probably
      also a processing non-starter since it would add overhead to the
      calculations, however it is the type of thing that could be parallelized
      during each iteration. No idea how efficient it would be to split threads
      during each iteration.
   
   * Code it in Rust, use R package `expandr` to integrate the Rust function
   into the package via an R function. Follow C++ code steps used by
   `qgraph:::qgraph_layout_Cpp()`. It doesn't seem that these layout
   iterations are particularly good for multi-threading, however.
   
      * each iteration could be multi-threaded. Threads could share one
      distance matrix, the extract values when needed.
      Or threads could share node coordinates, then calculate distance
      when needed on the fly. If fast enough, distance matrix (and memory
      allocation) could be avoided, just do the math when needed.
      * Look for Rust libraries for geometric calculations.

* `mem_gene_path_heatmap()`

   * add option to display gene incidence matrix (left annotation)
   using `geneIMdirection` colors, to represent up/down.

* add `plot_cnet_heatmaps()`

   * make default arguments as minimal and close to default settings as possible


## 10nov2022

* `adjust_cnet_nodeset()`

   * COMPLETE: add minimum percent spacing, similar to `apply_nodeset_spacing()` except
   to operate on only one nodeset at a time.

## 04nov2022

* `shape.jampie.plot()` which calls `jam_mypie()`:

   * When drawing borders for multiple pie shapes, the overlapping border
   are overdrawn, hiding all but the last border drawn.
   * Instead, draw border inside the polygon edge, so the borders are
   adjacent and not overlapping.
   * Quick survey revealed no drop-in replacement `polygon()` functions.
   * Workaround involves `sf::st_buffer()` to calculate a buffer inside
   each polygon, creating a "donut" filled with the border color. In
   principle, this issue tends to occur in relatively few nodes for typical
   Cnet plots, however it could be substantial performance hit.
   
      1. Each polygon, convert to corresponding `sf::st_polygon` object
      2. Call `sf::st_buffer()`.
      3. Create polygon from original border, inner buffer border.
      Or split the original polygon at this border, apply color to the outer
      and inner polygons. This way the border and fill colors are drawn
      together so nodes are "complete".


## 26oct2022

* `multiEnrichMap()`

   * accept `p_cutoff` and deprecate `cutoffRowMinP` for future removal.
   `p_cutoff` is used in other package functions, as is `min_count`.

* `mem_enrichment_heatmap()`

   * problem: `style="dotplot"` and `style="heatmap"` have very different
   visual effects, the dotplot de-emphasizes significance in favor of
   gene count (point size) - smaller gene count hides significance.
   The pale bivalent colors are very close to white, while `"Reds"` color
   gradient is much more distinctive.
   * Can try transforming point size, so the minimum size starts out larger?
   Arguments `point_size_min=3` and `point_size_max=10` may help?
   * Another more dramatic option is to create normal heatmap, then just
   draw points on top of the already-filled cells. Note the points would
   be inside the boxes, instead of connected by lines through the center of
   each box. Could test with `style="hybrid"` and implement both.

## 13oct2022

* S4 objects?

   * `mem`:

      * Idea is to have one object type to handle output from `multiEnrichMap()`
      * benefit is mostly to have default behaviors for things like `print()`
      but also for clarity during the workflow. The `mem` object can be
      clear input to other functions.
   
   * `cnet`:
   
      * inherits `igraph` and extends assumptions, again mainly for clarity
      as input argument to other functions

* `mem_gene_path_heatmap()`

   * option to display dot plot format at the top of heatmaps, equivalent
   to calling `mem_enrichment_heatmap()`. Benefit would be to indicate
   direction (color) and number of genes (size) in addition to
   P-value (intensity).

* fix `reorderIgraphNodes()`

   * DONE: add "frame.color":
   `sortAttributes=c("pie.color", "pie.color.length", "coloredrect.color", "color", "pie.border", "frame.color", "label", "name")`
   * future: make sort more intelligent, so it uses appropriate
   color based upon node shape during the sort.

      * `shape="pie"`: use "pie.color", "pie.border", "frame.color", "label", "name"
      * `shape="coloredrect"`: use "coloredrect.color", "coloredrect.border", "frame.color", "label", "name"
      * all others assumed to be `shape="circle"` or similar:
      use "color", "frame.color" and ignore "pie.color", "pie.border",
      "coloredrect.color", "coloredrect.border", "label", "name"

* DONE: Include gene direction in the workflow:

   * DONE: Add `geneIMdirection` to mem object.
   
      * Requires `geneHitIM` input to `multiEnrichMap()`.
      * Alternate backwards compatibility is function to add
      `geneIMdirection` to `mem` object, and/or `cnet` `igraph` object.
      * DONE: `memIM2cnet()` should optionally use `enrichIMdirection` to apply
      border color.
   
   * DONE: `mem2cnet()` as its own function, minor variation and extension to
   `memIM2cnet()` which only uses the `memIM` portion of the data.
   * `jam_igraph()` should recognize `lwd` when plotting nodes. This change
   diverges from default behavior in `igraph` which only uses `par("lwd")`
   as a global adjustment for all lines while plotting nodes, therefore
   does not allow any nodes to have different line widths.
   * DONE: Outline of Cnet gene nodes by direction.
   * Gene-Path heatmap rows should optionally indicate direction,
   possibly another stripe using geneIMdirection.

* `multiEnrichMap()`

   * should call `mem2cnet()` - or just avoid this step altogether since
   the `cnet` process is better done with `mem_plot_folio()`
   * consider removing the `enrichMap` and `cnet` steps altogether.

## 26sep2022

* `mem_plot_folio()`

   * FIXED: `mem_enrichment_heatmap()` does not honor `do_plot=FALSE`.

* `mem_enrichment_heatmap()`

   * DONE: add argument `cluster_rows` to allow `cluster_rows=FALSE`.
   * DONE for manual plot calls: one goal is to draw this heatmap
   using the order from the gene-pathway heatmap clustering.
   * Future idea: allow plotting this data using the same order
   as the gene-pathway heatmap. This process would require
   running the gene-pathway clustering first, determining the
   column order, then using it to order the rows in this heatmap.

## 31aug2022

* Now that `jam_igraph()` and node shapes `"jampie"` render the
`pie.border` and `frame.color`, which can indicate the direction
of change for each gene, in context of each enrichment test, some
other changes should be made to the workflow:

   * Consider adding `pie.lwd` to recognized attributes in
   `shape.jampie.plot()`, but confirm that this argument can be used
   when rendering `"jampie"` node shapes. Currently the line width
   can only be adjusted with `par("lwd")` which is a global setting.
   * The default `frame.color` (`vertex.frame.color`) should probably
   be set to `NA` or `"transparent"` so the frame color is not visible
   for pie nodes. It shows up as a small white line now, which was not
   visible previously.
   * `reorderIgraphNodes()` should include `pie.border` and `frame.color`
   in sensible default locations, so nodes will also be sorted by these
   values when relevant, without the user having to add these columns.
   * New function to populate `frame.color` and `pie.border` for Cnet plots.
   
      * When `pie` has only one value, apply color to `frame.color`
      * When `pie` has multiple values, and different directions, apply
      colors to `pie.border`
      * When `pie` has multiple values, and all have the same direction, apply
      colors to `frame.color`
      * In absence of any directional data, set all `frame.color` to default,
      and set all `pie.border` to `"transparent"` or `NA`.

* Need to include gene direction of change in the workflow:

   * Some easy method to include direction of change in the `multiEnrichMap()`
   workflow, for example argument `geneHitList` currently uses a `character`
   vector, but could accept `integer` vector direction with genes stored as
   `character` labels, same as used with `venndir` signed input lists.
   * When gene direction is available, the `frame.color` and/or
   `pie.border` colors are defined.
   * `mem_gene_path_heatmap()` option to represent direction of change
   in the `gene_im` incidence matrix.

## aug2022

* Edge bundling makes assumptions for bipartite graphs (cnet)
that are difficult to use with normal graphs.

   * `edge_bundle_nodegroups()` probably needs to subset edges
   involved in each nodegroup before bundling all edges connected
   to these nodes. Currently, nodes are assumed to share all
   connections, but the effect should occur for edges where both
   ends of the edge are contained in the nodegroup entries.
   * consider node attribute "`nodegroup`".
   * consider edge attribute `"edgegroup"`, which would probably be the
   optimal approach for edge bundling, apart from implementing another
   technique similar to force-directed, hierarchically-defined, or
   density-directed edge bundling.
   * consider an option to specify the "midpoint" for a nodegroup,
   to allow some control on the spline curvature. Bonus points for
   allowing multiple points in order, to influence a path.

* `jam_igraph()`

   * fancy effect: allow edge colors to have multiple values, then
   interpolate color along each edge. For bundled edges, make
   a gradient equal to the number of line segments. For straight
   edges, break into `detail` number of pieces. The `detail` argument
   is also used by the `edge_bundle_*()` functions.

* subset `igraph` object by nodes, with added benefit that the
`graph_attr(g, "layout")` will also be subset, if present. It is
odd that the default `igraph::subgraph()` function would subset
the graph, leaving the layout which inevitably causes an error.

## 13jul2022

* `mem_gene_pathway_heatmap()`

   * would be useful to have an option to subset enrichments,
   this option would likely be useful to `mem_plot_folio()`, and
   `mem_enrichment_heatmap()`, so the effect could cascade to
   subsequent functions.

* new function idea? `subset_mem()`

   * Could be useful to accomplish the item above, to subset by enrichments
   prior to `mem_plot_folio()` related functions.
   * The challenge is that subsetting only `enrichIM` does not also
   update corresponding `memIM` data. Data to be updated:
   * `geneIM`, `geneIMcolors`, `geneIMdirection` - simple subset by colnames
   * `enrichIM`, `enrichIMcolors`, `enrichIMgeneCount` - simple subset by colnames
   * `memIM` - re-create after updating `geneIM`
   * `enrichList` - simple subset by name
   * Others not necessary for `mem_plot_folio()`:
   
      * `multiEnrichMap`, `multiCnetPlot` should be re-created using
      methods in `multiEnrichMap()`. To be fair, these objects are not
      that useful anymore, since `mem_plot_folio()` is generally preferred.

## 05jul2022

* `mem_gene_pathway_heatmap()` when supplied with custom `column_split`
throws an error when `cluster_columns=TRUE`.

   * internally it uses `amap::hcluster()` to generate a dendrogram/hclust
   * using split and dendrogram together is not allowed by `ComplexHeatmap::Heatmap()`
   * this clustering also uses the incidence matrix combined with the
   pathway enrichment annotation displayed along the top, and these values
   are weighted with `pathway_column_weight`.
   * A proper solution would be to provide a custom function for `cluster_columns`
   that internally combined the incidence matrix and enrichment matrix data
   together prior to clustering. If this function receives a
   numeric matrix with proper colnames, it should work.

* COMPLETE: `mem_plot_folio()` argument `gene_row_title=NULL` is being passed
to `mem_gene_pathway_heatmap()` and is therefore not using the default
`row_title=letters`.

## 23mar2022 issue #7 passing arguments through ... in `mem_plot_folio()`

* User reported an error when calling `mem_plot_folio(mem, node_factor=5, label_factor=1.5)`
* `mem_plot_folio()` is passing `...` to `ComplexHeatmap::Heatmap()` which
does not accept `...` and throws an error when receiving extra arguments.
* I need to limit `...` to arguments accepted by `Heatmap()`. I feel like
there is an R package function to handle this scenario, so I can avoid
writing `do.call(Heatmap, custom_args)`.

## 04feb2022

Workflow that starts with pathway-gene incidence matrix upfront,
bypassing the `multiEnrichMap()` workflow altogether.

* `memIM2cnet()` to convert pathway-to-gene incidence matrix to Cnet `igraph`.

   * Requires `geneIM` which is the enrichment-to-gene incidence matrix.
   In the driving example, each enrichment would represent a disease subgroup,
   with the full set of differential genes associated to each subgroup.
   * Requires `enrichIM` which contains enrichment-to-pathway whose values
   are P-values. When this data is not supplied, it should use 0.001 by default.
   I think these values are only used as an optional gradient color for the
   pathway nodes.
   * Optional `geneIMcolors` which is the same as `geneIM` except the cells
   contain colors to use for each gene. Currently the function does not
   fill these colors, the best method is to use `colorjam::matrix2heatColors()`
   
   ```R
   geneIMcolors <- colorjam::matrix2heatColors(x=geneIM,
      colorV=colorV)
   ```


## 03feb2022

* option to assign pathways to "functional themes" based upon presentation
by Adeline Chin in  Dr. Hanna Kim' group.

   * This step may "collapse" multiple pathways together, similar to
   grouping functional groups:
   union of genes,
   lowest enrichment P-value.


## 25jan2021

* `heatmap_row_order()`, `heatmap_column_order()` should return a flat
vector when there are no row_split or column_split, respectively.
Currently data is returned as a list of one-length vectors.

   * Consider moving into `jamba` package, in case this function
   needs to be re-used, it should not require loading the full
   `multienrichjam` package, which itself requires things like
   `igraph`, `clusterProfiler`, `qgraph`, `DOSE`, `matrixStats`.
   * Import `jamba::heatmap_row_order()` and `jamba::heatmap_column_order()`
   to ensure any functions or packages calling this function will
   succeed without error.


This document describes plans for enhancements to the
multienrichjam R package.

## 05nov2021

Now that directional z-score can also be associated with enrichment
P-values, heatmaps might need to use a bivariate color scale, to
indicate enrichment and directionality. See "stevens.bluered".
For example the `mem_enrichment_heatmap()` colors nodes by enrichment
P-value, more intense is more significant enrichment.
The z-score direction is used to apply red "activated" or blue "inhibited".
However, pathways with no z-score, or z-score below the threshold
are colored red by default. They should use a neutral color.



## 04nov2021

For IPA "Upstream Regulators" it sometimes offers a direction
implied by the `activation z-score`. Design idea is to implement
the directionality so it can be included in downstream analyses.

   * `mem_enrichment_heatmap()` - currently shades by the `-log10(pvalue)`
   however if there is directionality, it could be signed
   `+` for activated,
   `-` for inhibited. Then the heatmap color scale would use blue-white-red
   color gradient.
   * `mem$enrichIMdirection` contains matrix of direction, by default `1`
   means all have same direction.


During `multiEnrichMap()` it filters for `topEnrichN` entries for
each enrichment. It might be useful to retain the rank number for
each enrichment, to review when setting a different `topEnrichN`
threshold.

   * `mem$enrichIMgeneCount` contains matrix of gene counts
   * `mem$enrichIMrank` contains matrix of pathway rank (after filtering gene count)


## 01nov2021

* COMPLETE: `mem_enrichment_heatmap()` - the heatmap circles and legend circles
are not the same size - they should be fixed to the same absolute size.

   * COMPLETE: Optionally label each heatmap cell with the number of genes for
   visual reinforcement.


* COMPLETE: When there are more than 3 enrichments, the color legend on the
gene-pathway heatmap becomes unwieldy - taking over the whole figure.

   * Optionally (and by default) hide color legend for the gene-pathway heatmap.

* COMPLETE: The gene-pathway heatmap use_raster=TRUE causes artifacts in output,
it should be disabled by default. In future, debug why things go wrong.

   * I think the bug is caused by rasterization being done on underlying numbers
   before the color ramp is applied, in which case this problem cannot be
   solved when `colorize_by_gene=TRUE` since categorical values are assigned
   numbers which are not a proper color gradient.



## 20jul2021

* Replace `multiEnrichMap()` with `multienrichjam()` and simplify the arguments:

   * p_cutoff
   * min_gene_count
   * top_enrich_n
   * colnames: id, name, description, pvalue, gene
   * color_sub

* Port `mem_enrichment_heatmap()` argument `colorize_by_gene=TRUE` to
use `ComplexHeatmap::Heatmap()` instead of `jamba::imageByColors()`
for consistency, also so it can support `style="dotplot"`.

## 19jul2021

* `jam_igraph()` with `rescale=TRUE` should also scale igraph vertex
size and igraph label size according to the new axis ranges.
* Debug `jam_igraph()` when `vertex.size` is defined alongside
V(g)$size, and `node_factor`. It appears not to apply the size
properly.
* COMPLETE: `mem_enrichment_heatmap()` new option for dot plot format,
based upon `enrichplot::dotplot()` that sizes each dot by the
number of genes present.

## 24jun2021

* In `multiEnrichMap()` remove default `topEnrichSources`
and `topEnrichSourceSubset` which throw errors when not
using MSigDB data.
* `topEnrichBySource()` and `topEnrichListBySource()` should
be able to accept `enrichResult` as input, and return `enrichResult`
and not `data.frame` which is understandably lossy.
I need to understand the proper method for creating a subset
of an `enrichList` object, including its internal data.
* Streamline the `topEnrichListBySource()` workflow.


## 21jun2021

* COMPLETE: Allow rotating gene-pathway incidence matrix when using
`mem_gene_path_heatmap()`, so the pathway names are rows.
* Add test data object `mem` to be convenient for function examples,
and testthis test suite.


## Smaller usability items

* edge_bundling="connections" should also allow `render_groups=TRUE`
to work, by returning the `"nodegroups"` required for that step.

   * It should probably filter out singlet nodes by default, since
   the typical use case is "nodeset-to-node" for Cnet plots.
   * Main use case is with a Cnet plot, edge bundling by connected nodes
   should allow drawing an optional border around each group of nodes.



## COMPLETE: 10apr2021: implement edge bundling techniques

Priority: high 

Status: Implemented, early active testing of usability and functionality

* "connections" - Cnet edge bundles - special case of nodeset-to-node edge bundling
* "nodegroups": Node group bundles - general case of nodeset-to-nodeset edge bundling

Useful to improve readability/aesthetics of collapsed Cnet plots,
especially Cnet plots with many gene nodes where it is
difficult to tell which pathway clusters are connected
to each gene.


## 10apr2021: enhancements to `jam_igraph()` (must be done inside `jam_plot_igraph()`)

* Consider silently returning the `igraph` object plotted
where the vertex, edge, and graph attributes are updated
to the values used at runtime. For example, if user overrides
any attributes, those attributes will be present in the
object returned, so the next iteration someone could
just call `jam_igraph()` without any custom attributes
and it would produce the same plot again.

* See https://igraph.org/r/doc/plot.common.html
and `igraph::igraph_options()`



## 10apr2021: implement edge bundling in `jam_igraph()` and `jam_plot_igraph()`

Priority: high

Status: Implemented

* This task involves making edge bundling accessible as a simple
option during plotting, to prevent having to run 3 or 4 functions:

   * `jam_igraph()` with hidden edges
   * `bundle_node_edges()` to display edges on top of nodes
   * (optional) `jam_igraph()` with hidden edges, nodes with solid white
   color so they fully cover edges.
   * `jam_igraph()` with hidden edges, to display nodes on top
   of bundled edges while allowing nodes to have alpha transparency

* Note there is a package based on `ggraph` that implements
edge bundling - if `jam_ggraph()` is going to be our future,
maybe we should implement edge bundling using a similar approach
used by that R package.

   * https://github.com/schochastics/edgebundle


## 10apr2021: add `plot_cnet_heatmaps()` to multienrichjam

Priority: high

* `plot_cnet_heatmaps()` is in development, and arranges expression
heatmaps around a central Cnet plot, using genes in each cnet
cluster. It is intended to be used with Cnet clusters from
`collapse_mem_clusters()`, by-product of `mem_plot_folio()`.

This figure seems very useful because it integrates expression changes
alongside gene-pathway connections.

The multienrichment workflow would likely become:

* Prepare enrichment data
* `multiEnrichMap()`
* `mem_plot_folio()`
* `plot_cnet_heatmap()`



## 10apr2021: enable plotting node community data in `jam_igraph()`

Priority: low

* Analogous to `igraph:::plot.communities()` with `mark.groups`
* for edge bundling, this step could optionally display boundaries
around each node group.

Value is reinforcing node groups with a boundary. In testing,
the boundary actually made plots look more complex. And with simple
plots, it does not seem necessary since nodes are already
well-spaced.



## 10apr2021: include geneCount matrix output from `multiEnrichMap()`

Priority: medium

* Purpose is to complement the "enrichIM" matrix, with P-values
for each set. However "dysregulated pathways" may also be required
to have N number of genes, yet this data is not easily available.
* Rows are pathways, columns are enrichments, values are the number
of genes involved in enrichment.

   * Consider a matrix whose cells contain delimited gene symbols

Useful to apply filtering at matrix-level for pathways that
meet enrichment P-value and gene count thresholds.


## 2020: apply node color order based upon shape of node cluster

Priority: medium

`reorderIgraphNodes()`

* when a cluster of nodes is short-wide, the order should be left-to-right
* when a cluster of nodes is tall-skinny, the order should be top-to-bottom
* all else should be sorted left-to-right (or user-defined default order)

   * it is visually confusing when tall-skinny nodes are sorted left-to-right

Useful to automate the ordering of node colors


## 2020: apply subsetSets in `multiEnrichMap()`

Priority: low

* apparently it had never been implemented? User-defined
subset of pathways/sets to include in `multiEnrichMap()`,
as an alternative to using `topEnrichBySource()`.

Users can currently perform this subset step before running
`multiEnrichMap()`.


## COMPLETE: Apply min_count in `multiEnrichMap()`

* Currently `multiEnrichMap()` does not filter by number of
genes involved in enrichment, it only filters by enrichment P-value.
New argument `min_count` is applied only when `topEnrichN` is used,
but nothing else downstream is aware of filtering by `min_count`.
The corresponding argument in `mem_plot_folio()` is
`min_set_ct_each`, which requires a set (pathway) to contain
at least this many entries in at least one enrichment result
which also meets `p_cutoff` criteria for enrichment P-value.


## 2020: Optional highlight genes

Priority: low

* `mem_plot_folio()` and subsequent plots, optional argument
`highlight_genes` which would effectively hide all gene labels
except `highlight_genes` -- to help especially crowded plots.

   * Option 1: User-supplied genes
   * Option 2: Genes defined by one or more pathways/sets to highlight
   * Implement as wrapper to `jam_igraph()` to hide/display relevant node labels
   * `plot_cnet_heatmap()` could label the heatmap rows using `anno_mark()`

Useful in Cnet plots with too many genes to label, it allows
only a subset of genes to be highlighted.
Will become high priority if a manuscript decides to use
this technique.


## Consider a strict output format for `mem_plot_folio()`

Priority: medium

* Goal is to run `mem_plot_folio()` once then be able to plot
any component separately. Would help keep all plots in sync
with the settings used.
* Benefit: Always have ability to see the exact gene-pathway
incidence matrix heatmap, and its clustering, which is used
for a collapsed Cnet plots.
* Design idea: List format

   * Named by type of output:
   
      * enrichment P-values
      * gene-pathway heatmap
      
         * Heatmap object
         * gene clusters
         * pathway clusters
         * filter/cluster settings used
         
      * collapsed cnet igraphs
      
         * with each labeling option implemented
      
      * cnet exemplar igraphs
      
         * 1 exemplar per pathway cluster
         * 2 exemplars per pathway cluster
         * 3 exemplars per pathway cluster
         
      * each cnet cluster igraph
      
         * cluster 1 cnet igraph
         * cluster 2 cnet igraph
         * cluster 3 cnet igraph
         * ...
      
   * Each type of output contains a `list` of relevant components


## 10apr2021: make plot commands word for ggraph

Priority: low

* Goal is to use ggraph for ggplot2-type plotting instead
of using base R for `igraph` plots.
* `jam_ggraph()` ?

   * Test `scatterpie` R package which implements pie
   node shapes for use in ggraph. Unsure the data content
   requirement.


## 10apr2021: document examples for `reorderIgraphNotes()`

* Examples should show how to specify color order, or change
the order if necessary.


## COMPLETE: 2020: Improve `reorderIgraphNotes()`

* `reorderIgraphNodes()` when it encounters attributes with multiple
colors per node, such as `"pie.colors"` and `"coloredrect.colors"`,
calls `avg_colors_by_list()` to generate one blended color per node,
then it sorts those colors using `sort_colors()`. However, the
color order ultimately does not match the order in the color legend,
and other plots such as heatmaps.
* Future idea is to convert node attributes with multiple colors per
node into a `data.frame` where each color is in a separate column.
Then convert each column to a factor whose levels match the color
order from `colorV` (argument to `multiEnrichMap()`). The end result
should sort nodes consistent with the order of colors.
* Note this process assumes that all node colors are from a limited
set, which ideally should match `colorV`. Therefore if any color
gradient is applied to nodes with `nodeType="gene"`, this process
will not work.
