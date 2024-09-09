# multienrichjam 0.0.89.900

## changes to existing functions

Added `%>%` to imported functions, and not all of `dplyr`.

* `mem_enrichment_heatmap()`

   * New argument `cell_size` default `NULL`, with option to specify the
   exact cell size used in the heatmap. The driving use case is to define
   perfectly square heatmap cells for `"dotplot_inverted"` or `"dotplot"`
   so the circles are perfectly centered inside square cells. It forces
   the output figure to be tall and wide enough to accomodate the resulting
   figure, so it needs some user math. Implemented from a user suggestion,
   and it works and looks great. It just requires some upfront work to
   create a figure with large enough canvas.

* `mem_plot_folio()`

   * Added proper pathway cluster labels to `mem_enrichment_heatmap()` output.
   * Finally fixed the pathway order of `mem_enrichment_heatmap()` to be
   identical to `mem_gene_path_heatmap()` without minor re-ordering caused
   by default `ComplexHeatmap` behavior. Now the two formats have identical
   order for direct comparison.
   * Cleaned up the help documentation a little bit.

# multienrichjam 0.0.88.950

## changes to existing functions

* `mem_enrichment_heatmap()`

   * New option `style="dotplot_inverted"` due to a cool recent paper
   using this style. I like it so much, it is the new default.
   Much easier to see the color and the size of the circle, especially
   for cells with very small circles.
   From Jang et al, the Waggoner lab, Nature Genetics 2024:
   https://doi.org/10.1038/s41588-024-01880-x

* `cell_fun_bivariate()`

   * New argument `invert=FALSE` to control whether to draw colored circles,
   or colored cells with white circles on top (`invert=TRUE`).

# multienrichjam 0.0.88.900

## changes to existing functions

* `mem_plot_folio()`

   * Moved workhorse function to its own .R file.
   * Changed default: `verbose=FALSE`.
   * New argument `rotate_heatmap` to handle this workflow properly.
   * New arguments `row_anno_padding`, `column_anno_padding` to control
   the padding between heatmap and row and column annotations, respectively.

* `mem_gene_path_heatmap()`

   * The caption is displayed as a ComplexHeatmap Legend, therefore
   using consistent font and alignment with the color legends,
   rather than being tucked into the corner where it sometimes
   overlapped other heatmap labels.
   * Caption now displays the rows/columns counts first.
   Remaining values are more user-friendly.
   * Default caption font size 10 instead of 6, consistent with
   color legend text.
   * Argument `seed` is utilized.

* `make_point_hull()`

   * New argument default `label_preset=NULL`.
   * Argument `label_preset` is properly recognized.

* Overall: Functions with argument `seed` now properly ignore the seed
when it is `NULL`, thereby allowing random behavior when preferred.

# multienrichjam 0.0.87.900

The next update will likely use `Mem` S4 object instead of `list`,
although the plan is to allow convenient conversion to `list` for
legacy compatibility. The S4 functions and methods should be more
convenient in long term, and will be Bioconductor-compliant.

## new functions

* `add_pathway_direction()`

   * Adds column to `enrichResult` indicating the directionality,
   using the IPA Activation z-score calculation (refs in help doc).
   It requires `gene_hits` as a numeric vector named by gene symbol.

## changes to existing functions

* `topEnrichBySource()`, affecting `topEnrichListBySource()`, `multiEnrichMap()`

   * now applies logic for `descriptionGrep`, `nameGrep`, and `subsetSets`
   using OR logic, so any combination of matching results will be
   retained.
   * Help docs have been updated.


## other changes

* `.onLoad()`

   * Minor change, it now checks whether the new `igraph` shapes already
   exist before adding, mostly useful when reloading this package
   in a live R session, which should be rare.



# multienrichjam 0.0.86.900

## Bug fixes

* `topEnrichBySource()`, `topEnrichListBySource()`, `multiEnrichMap()`

   * Weird rare scenario that appears limited to custom enrichment data
   where input data contains multiple P-value columns, none of which
   match the defaults for argument `sortColnames` in `topEnrichBySource()`.
   * First bug: The `sortColnames` argument was not evaluated with
   `find_colname()`, which is intended to match the patterns to actual
   data. This decision was probably to avoid handling the optional
   prefix `"-"` to reverse sort by colname.
   * Second bug: The `sortColnames` should really only use
   `pvalueColname` then reverse order for `countColname`.
   * End result: When `sortColname` matched no existing columns,
   the data was still sorted using `jamba::mixedSortDF()`.
   It found no matching colnames so by default it sorted
   starting with the first column. This is incorrect, and caused the bug.
   * New argument default: `sortColname=NULL` will use `pvalueColname`
   then reverse order of `countColname`.
   Alternatively `sortColname=FALSE` causes no sort to be performed,
   using data in the order it was provided.
   Finally, if `sortColname` is defined, and its values do not match
   existing colnames, no sort is performed.

## changes to existing functions

* `multiEnrichMap()`

   * Minor change to silence the text output when `enrichIM` is entirely NA,
   which can occur when the input data does not contain any significantly
   enriched pathway results. This outcome causes `stop()` but should not
   otherwise print output unless `verbose=TRUE`.

* `importIPAenrichment()` and `topEnrichBySource()`

   * Both now use `jamba::gsubs()` instead of the "temporary" internal
   function by the same name.

* `gsubs()` is removed, renamed `gsubs_remove()` in favor of `jamba::gsubs()`.

## other changes

Two functions are no longer imported, instead they are both called using
the package prefix. It only affects `grid_with_title()` which is no longer
used by default.

   * `ComplexHeatmap::draw()`
   * `jamba::nameVector()`
   

# multienrichjam 0.0.85.900

## Changes to existing functions

* `importIPAenrichment()`

   * New argument `revert_ipa_xref=TRUE` changes the default behavior
   (for the better) such that the resulting gene symbols associated
   with IPA pathway enrichment will match the input gene symbols,
   instead of using the customized IPA symbols.
   
      * The situation does not have an ideal solution. Ultimately, IPA
      provides results which do not completely represent the input data.
      When two genes are combined to one entity by IPA, they only retain
      one gene symbol in `"Analysis Ready Molecules"`, and so there
      is no recorded association of all entities which were combined.
      Sometimes the combined symbol is `"HSPA1A/HSPA1B"` which can be
      separated... but they do not indicate which symbol(s) the user
      provided, they only record one. Some entities are called
      `"NBPF10 (includes others)"` and which appears to include
      `'NBPF10"` and `"NBPF19"`, possibly others.


# multienrichjam 0.0.84.900

## Bug fixes

* `make_point_hull()`

   * Finally "fixed" the full bug, sometimes causing point hull to fail.
   Sometimes `alphahull::ahull()` would return weird results with small
   `alpha` values: individual points, segments, or multiple disconnected
   polygons.
   * Now `make_point_hull()` calls `get_hull_data()` which performs additional
   validation checks: Confirms there is actually a polygon; confirms each
   point in the `edges` are used exactly twice; confirms that all edge points
   are used in a continuous polygon, not two separate polygons.
   * The default value for `max_iterations=100` is vastly increased, since
   many of these weird cases were caused by having coordinate ranges
   orders of magnitude higher than `alpha=0.1` and so 10 iterations would
   not be enough to avoid these weird situations.
   * Also when `make_point_hull()` does not find a suitable solution after
   `max_iterations` tries, it returns `NULL` instead of proceeding to
   process the inadequate/empty polygon coordinate data. Hopefully this
   will allow weird cases to be skipped rather than throwing an obscure
   error.

# multienrichjam 0.0.83.900

## Bug fixes

* `make_point_hull()`

   * One of the more obscure bugs, caused by creating a point hull with
   only two points, where the two points apparently resided at exactly the
   wrong angle. Only apparently occurred when the two points were within 3
   degrees of some "evil angle" - resolved when points were rotated
   more than 3 degrees from this "evil angle".
   Typically, a point hull is not possible with two points, since it does
   not create a polygon, so the workaround was to add a third point (recycling
   the first point) adding amount of noise to create a polygon. Except for
   some reason even `rnorm()` was not adding *enough* or the right kind
   of noise. No random value should ever exactly, reproducible, reside
   on the line. There must be some rounding that takes place in an
   internal function. Nonetheless, the "fix" for this purpose was
   to add more than one dummy point.

* `list2im()`

   * Another bizarre "bug" was caused by some unconfirmed R function
   that appears to change `options("warn")` to `options("warn"=2)`
   - and then does not change it back! Probably related to RMarkdown
   knitting, since it seems to occur when the knitting is interrupted.
   To be fair, warnings should be resolved in Jam packages, that's true.
   But they should not impose an error.
   * Filed under: "Things that worked just two minutes ago, but now
   cause an error for no reason."
   * Added workaround for random but painful issue when `options("warn"=2)`
   which forces warnings to errors. No idea what made that setting, but
   it caused `list2im()` to fail due to implicit conversion of `empty`
   to whatever datatype was defined in the input matrix. In future this
   function may change, but for now this change keeps it working.

# multienrichjam 0.0.82.900

* `multiEnrichMap()`

   * now `geneHitIM` and `geneHitList` are equivalent as input, either
   or both can be provided. The incidence matrix with values c(0, 1)
   are stored as `geneIM`, and if there are any other numeric values
   those values are stored in `geneIMdirection`.
   It does not (yet) verify that all genes involved in enrichment are
   also present in the `geneIM` and `geneIMdirection` matrices.

* `mem_gene_path_heatmap()`

   * New argument `gene_annotations` allows display of the `mem$geneIM`
   (gene hits per enrichment) and `mem$geneIMdirection` (gene direction
   of change) as annotations alongside the gene axis of the heatmap.
   By default, when `mem$geneIMdirection` is available, both are shown,
   otherwise only `geneIM` is shown. It can now be hidden.
   * New argument `simple_anno_size` to control the size of heatmap annotations.
   * New argument `annotation_suffix` to add optional suffix to the
   gene annotation labels, helpful to indicate the unit. For example
   `im` data has default suffix `"hit"`, while `direction` uses `"dir"`.
   * When `mem$geneIMdirection` is present, it is now used by default
   during clustering. Values are multiplied by `mem$geneIM` and any `NA`
   values in `mem$geneIMdirection` are replaced with `1` in order to
   maintain the `mem$geneIM` value. They are considered `unknown direction`
   in that sense.

# multienrichjam 0.0.81.900

* `jam_igraph()` and `jam_igraph_arrows()`

   * fixed bug when plotting directed edges with arrows, the arrow head
   width (`h.lwd`) was not properly expanded to the number of arrows,
   causing arrow heads to appear "twisted".
   * changed default obscure option `sh.adj=1` which draws straight edges
   to the base of the arrow head, not to the arrow tip. It allows the
   arrow tip to be a point without needing to equal the edge line width.

# multienrichjam 0.0.80.900

## changes to existing functions

* `make_point_hull()`

   * Silenced some default verbose output.


# multienrichjam 0.0.79.900

## changes to existing functions

* `spread_igraph_labels()`

   * argument `nodeSortBy=c("x", "y")` changed to `nodeSortBy=c("x", "-y")`,
   consistent with top-to-bottom sorting on the y-axis.
   This change indirectly affects `mem_plot_folio()` default output.

# multienrichjam 0.0.78.900

## changes to existing functions

* `mem_enrichment_heatmap()`

   * The bivariate color scale was too pale for lower significance P-values,
   so the colors were encoded to have higher color saturation at the low end.
   The intermediate colors are improved, see `colorRamp2D()` changes below.

* `colorRamp2D()`

   * The blue-yellow color blending by default in `circlize::colorRamp2()`
   was still producing grey, despite using the "LAB" (or "LUV", "XYZ")
   color models. This is usually a symptom of using RGB color space,
   blending "blue" with "red/green" (yellow) produces "red/green/blue" (grey),
   and not usually seen when using "LUV" which is a 360-degree hue radial
   color wheel. It is probably still converted to RGB before blending,
   then back to LUV.
   * The default argument was changed to `use_model="sRGB"` which
   produces a somewhat green color when blending blue and gold.
   The red-gold blending is improved as well, producing a more saturated
   orange color. End result: The intermediate directional colors are
   more recognizable as partly up (orange) or partly down (green).
   * TODO: This color legend needs x-axis labels, showing z-score values.

* `make_legend_bivariate()`

      * finally displays the x-axis label and numerical units,
      by default `"z-score"`.

## bug fixes

* `reorder_igraph_nodes()`,`reorderIgraphNodes()`

   * error occurred when attributes in `list` form were entirely `NA`,
   one section did not use `unlist()` properly. The bug was fixed.


# multienrichjam 0.0.77.900

## bug fixes

* `mem_plot_folio()` and `mem_enrichment_heatmap()` threw errors

   * using `pathway_column_split` when there was only one column
   (one enrichment), this error is corrected.
   * Further, passing `cluster_rows` as a `function` caused the
   resulting `Heatmap` object not to store the `obj` with the
   dendrogram/hclust, instead it stored the function.
   This error was also corrected.

# multienrichjam 0.0.76.900

## changes to existing functions

The theme of this update is "Customizing mark.groups labeling".
Useful with `jam_igraph()` when supplying a `list` of node groups
or communities via `mark.groups`, with labels defined in
`names(mark.groups)`. The labels are automatically placed
outside the mark polygon, and can be adjusted in size
with `mark.cex`, and position with `mark.x.nudge`, `mark.y.nudge`.

* `make_point_hull()`

   * new arguments: `label.cex`, `label.x.nudge`, `label.y.nudge`
   to customize the label font size, and label placement, when
   `label` is provided.

* `jam_igraph()` new arguments:

   * `mark.cex` passed to `make_point_hull()` as `label.cex`
   * `mark.x.nudge` passed to `make_point_hull()` as `label.x.nudge`
   * `mark.y.nudge` passed to `make_point_hull()` as `label.y.nudge`


## major changes to existing functions

* `mem_plot_folio()`

   * New default creates the gene-pathway incidence matrix
   heatmap data first, which serves as the basis for all other plots.
   * Impact on `mem_gene_path_heatmap()`, now inherits this
   clustering defined using the gene-pathway incidence matrix
   and only when called by `mem_plot_folio()`.
   This change is the primary motivation for the update, so that
   the enrichment heatmap clustering is informed and driven by
   gene content, and no longer reflects only the enrichment P-values.

## other changes to existing functions

* `mem_enrichment_heatmap()`

   * default `point_size_min=1` changed to `point_size_min=2` so the
   smallest points are still clearly visible, including the fill color.
   * default `p_cutoff=1e-6` changed to `p_cutoff=1e-10`

* `mem_gene_path_heatmap()`

   * default `p_cutoff=1e-6` changed to `p_cutoff=1e-10`

## bug fixes

* `multiEnrichMap()`

   * Edge case where `geneHitList` could be supplied as a `list` of
   signed values, a directional hit list as used in `venndir::venndir()`.
   It did not get recognized, and was converted using the `numeric`
   values rather than using the names of the values.
   * The preferred option is to supply `geneHitIM` for signed data,
   however it now works with a signed hit `list` as well.

# multienrichjam 0.0.75.900

## bug fixes

* `importIPAenrichment()`

   * New arguments to handle forward-slash delimited gene symbols used
   by IPA to indicate two or more genes they consider to be one biomolecule
   for the purpose of pathway enrichment. The forward-slash "/" is also
   the delimiter user in `clusterProfiler` object `enrichResult` which
   causes these genes to break that compatibility. The new arguments
   represent a workaround to handle IPA data, so that downstream
   functions do not require changes.
   * `convert_ipa_slash=TRUE` enables the workaround, which converts
   forward-slash "/" to another delimiter.
   * `ipa_slash_sep=":"` defines the alternate delimiter to use, the default
   `":"` was chosed because it does not interfere with other common delimiters
   used in gene symbols, and does not cause problems with regular expressions,
   which would have been a risk with using `"|"`.

* `removeIgraphBlanks()` was throwing an error when `size2` was not
already defined in `igraph` vertex attributes. When there is no `size2`
it now calls `default_igraph_values()$vertex$size2` to use the appropriate
default value.

   * other related errors are now caught and avoided, relating to steps
   that avoid updating empty entries in `list` vertex attributes,
   now it properly ignores pie attributes which were previously empty,
   and only updates entries with non-zero results.

* `mem_gene_path_heatmap()`

   * Fixed longstanding errors when trying to split the heatmap
   by row or column into more pieces than the data will allow.
   The current workaround is almost complete, covers obvious cases
   where the requested number is higher than the number of columns
   or rows. It does not determine if there can only be one or two clusters
   but there are more columns or rows present. In future.

## other changes

* moved `enrichDF2enrichResult()` to a separate R file, for future
maintenance.


# multienrichjam 0.0.74.900

## updates to existing functions

* `apply_cnet_direction()`

   * argument docs were updated
   * default colors in `col` were defined with more saturated default
   colors: `c("blue", "grey80", "firebrick3")`
   * argument default changed to `col_l_max=80` to accommodate higher
   `"grey80"` middle color.

* `mem_legend()`

   * new argument `pt.lwd=2` to control line width of open circles,
   used when `do_direction=TRUE`. The previous alternative was to define
   `par(lwd=2)` prior to calling `mem_legend()`, then setting
   `par(lwd=1)` afterwards.
   * `directional_colors` use the same default colors used by
   `apply_cnet_direction()`; and now includes `"no change"` as
   a specific legend entry.


# multienrichjam 0.0.73.900

## bug fixes

* Fixed errors rendering `igraph` nodes with `shape="circle"`
only when called by `jam_igraph()`, and only with newer
versions of `igraph` R package that expect the new attribute
`vertex.frame.width` (not `vertex.frame.lwd` as I had hoped).

   * `default_igraph_values()` now defines `vertex.frame.width` and
   `vertex.frame.lwd`.
   * Note: `vertex.frame.lwd` is very likely to be removed from
   this package altogether, for compliance with `igraph`. However, I will
   do proper testing before making the change.

* Bug in `rank_mem_clusters()` which created a `data.frame`,
and did not specify `stringsAsFactos=FALSE`, is fixed. Possibly similar
bugs in other functions.
* Bug in `mem_plot_folio()` and `mem_gene_path_heatmap()` when enrichment
only involves one gene, causing the `row_split` and `row_title` values
to be incorrect. It is not a useful workflow, but the functions should
handle the error edge case.
* `apply_cnet_direction()` now checks input argument `cnet` and `hitim`
to make sure they are non-empty before processing. In theory it should
not happen, but apparently it does when enrichment results are sparse
or possibly empty.
* `jam_plot_igraph()` internal function uses proper `igraph::shapes(shape)`
instead of previously internal `list` object `igraph:::.igraph.shapes[[shape]]`.

# multienrichjam 0.0.72.900

## new functions

* `subset_mem()`

   * convenience function to subset an entire `mem` object by sets (pathways)
   or genes. It will subset all internal incidence matrix objects for
   consistency.
   * This approach will become the preferred approach to display
   a specific set of pathways in a Cnet plot, by calling `subset_mem()`
   then `mem2cnet()`.
   * `mem_gene_path_heatmap()` currently subsets `mem` data internally,
   but for consistency may call `subset_mem()` instead. The arguments
   are designed to be very similar.

* `get_igraph_layout()`

   * simple helper function to retrieve or define layout coordinates
   for an `igraph` object. Because there is some logic to the process,
   it makes sense to put into its own function.

## changes to existing functions

* `jam_igraph()`, `jam_plot_igraph()`

   * recognizes new `igraph` attributes: `vertex.label.fontsize`,
   `edge.label.fontsize`.
   
      * These are not standard `igraph` attributes and will not be honored
      by default `igraph::plot.igraph()` functions.
      * When `vertex.label.fontsize` is specified as a font size in points,
      this font size is used with no modification by `vertex.label.cex`.
      * When any `vertex.label.fontsize` value is `NA`, the default behavior
      is used to calculate font size, which uses `vertex.label.cex`.
      Therefore the `vertex.label.fontsize` can be defined for a single
      node, as long as all other node values are `NA`, and only the one
      node font size will be adjusted to this specific fontsize.
      * The `igraph` labels are drawn using `text()`, and the final exact
      font size is calculated for nodes:
      `par("ps") * par("cex") * vertex.label.cex`
      and for edges:
      `par("ps") * par("cex") * edge.label.cex`
   
   * new argument `label_fontsize_l` used to apply specific
   `vertex.label.fontsize` based upon node attribute values. For example
   `label_fontsize_l=list(nodeType=c(Gene=10, Set=14))` will define
   Gene nodes with fontsize 10, Set nodes with fontsize 14.

* `apply_cnet_direction()`

   * when `frame_blank=NULL` is passed as an argument, it is interpreted
   as `frame_blank=NA` which uses no `frame.color` for blank nodes.

* `reorder_igraph_nodes()`, `reorderIgraphNodes()`

   * The following change was made to default behavior, so any method
   that reorders igraph nodes by color/border/name will be affected.
   Nodes will by default be sorted left-right when the nodes in a nodeset
   have less than 25% the x-span (width) compared to y-span (height),
   otherwise nodes will be sorted top-bottom. The direction along each
   axis respects the original argument `nodeSortBy`, to allow specific
   order based upon the data, or the natural order per the locale.
   * new arguments `orderByAspect` and `aspectThreshold` control when
   to sort left-right or top-bottom.
   * When the nodeset coordinate aspect ratio is taller
   (25% higher y-span than x-span) the nodes are sorted top-bottom,
   otherwise nodes are sorted left-right.

## bug fixes

* `mem2cnet()` threw an error for certain custom input that did not
meet expected constraints. The function was updated to prevent these
errors and to be more robust to this type of issue.
* `jam_mypie()` is called when rendering `shape="jampie"` nodes.

   * Previously when `frame.lwd=0` and `frame.color="black"` a small
   black border was drawn around the node. The new behavior when
   `frame.lwd=0` is to replace the color with `NA` so no outer border
   is drawn.

* `reorder_igraph_nodes()`, `reorderIgraphNodes()`

   * Fixed bug where `nodesets` argument was not always matched due
   to truncating the nodeset label to 25 characters. Now the nodeset
   label is not truncated.

# multienrichjam 0.0.71.900

## changes to existing functions

* `mem_plot_folio()`

   * After gene-pathway heatmap clustering, a collapsed Cnet `igraph`
   is created using pathway column clusters. The pathway cluster nodes
   are colorized based upon the proportion of each enrichment in that
   cluster, however it uses `mem$enrichIMcolors` by default.
   This function now applies statistical thresholds `p_cutoff` and
   `min_set_ct_each` prior to this step so the resulting colors will
   reflect those thresholds.
   * Help documentation was updated to include this information.


# multienrichjam 0.0.70.900

## bug fixes

* `shape.jampie.plot()`

   * Rare scenario using jampie nodes, when vertex.pie.lty is not defined,
   the default was not properly expanded to vcount, causing error
   `"subscript out of bounds"` when referencing `vertex.pie.lty[[i]]`.
   * Error above is caused by missing node attribute `"pie"`, which
   are now filled in with uniform values of 1 based upon
   `lengths(vertex.pie.color)`. This scenario usually occurs when
   trying to create pie nodes outside the "typical" scenarios,
   for example manually assigning attributes and not populating
   all the necessary values.


# multienrichjam 0.0.69.900

## new functions

* `color_edges_by_nodegroups()`
* `communities2nodegroups()`, `nodegroups2communities()`

   * conversion functions that help interconvert between `igraph`
   `communities` objects, and `nodegroups` which is a `list` of
   node names.

* `mem2emap()`

   * replacement for `enrichMapJam()`
   * converts `mem` output to an `igraph` object with multienrichment
   features, such as `pie` nodes with appropriate color fill,
   `pie.border` colored by direction when `enrichIMdirection` is defined.
   * by default, the resulting network has community detection called,
   then visualized with boundaries around the various nodes.

## changes to existing functions

* `mem_gene_path_heatmap()`, `mem_enrichment_heatmap()`

   * Updated to pass `raster_device` to `ComplexHeatmap::Heatmap()` to
   work around temporary error when the `"magick"` package is not
   available, `use_raster=TRUE`, which causes an error during rasterization.
   The error is resolved when changing from default raster device to
   `raster_device="agg_png"`, although this change requires the `"ragg"`
   R package is installed. So the change tests if `"ragg"` is available,
   and if so it passes `raster_device="agg_png"`. The change should not
   affect any other scenarios.

* `edge_bundle_nodegroups()`

   * New argument `bundle_self=FALSE` changes previous default behavior
   by not bundling nodes that connect from and to the same nodegroups.
   Previously, nodes connecting within the same nodegroup would bundle
   through the center point of the cluster, which does minimize the
   busy edge lines, but makes it difficult to follow any paths.

* `make_point_hull()`

   * new arguments `label`, `label_preset`, `label_adj_preset` are ussed
   to define optional label to appear outside the resulting hull.
   The label is intended to be used for network communities, to
   allow a label associated with each community when relevant.
   * Label placement is experimental and could change in future.
   * Labels are placed relative to the center of the layout, using
   the angle from layout center to hull center. Labels are placed outside
   the rectangular bounding box of the point hull, with text aligned
   to the outer edge based upon the nearest 45 degree angle from
   layout center to hull center.

* `apply_cnet_direction()` default `frame_lwd=0.2` changed from `frame_lwd=1`.
* `jam_igraph()`, `jam_plot_igraph()`

   * new argument `bundle_self=FALSE` passed to `jam_plot_igraph()` and
   ultimately to `edge_bundle_nodegroups()`. When TRUE any edges that
   connect from and to the same nodegroup will be bundled through the
   nodegroup center. The default FALSE does not bundle within nodegroup,
   so edges are only bundled when connecting two different nodegroups.

* `shape.jampie.plot()` now handles various combinations of missing
`pie.lwd` and `frame.lwd` more gracefully.
* `mem2cnet()` uses smaller node size by default.

# multienrichjam 0.0.68.910

It turns out that `graphics::polygon()` only properly closes the
polygon when it has a color fill with non-zero alpha transparency.
In this case, the final "closed" corner is correctly extended
to complete the sharp corner edge using line join `par("ljoin")`.
With no color fill, or with completely transparent color fill,
the final corner is not completed, the lines are "ended" using
`par("lend")`, and therefore there is no sharp corner.
The workaround is to apply a color with alpha transparency 1
(on scale of 0 to 255), which causes the border to be drawn
completely. However, some rare graphical output devices do not
support alpha transparency, so there is the chance of rendering
unintended opaque color fill. The situation should only affect
node shapes `"pie"` and `"jampie"`, which
are designed to draw the outer border then inner border, so the
impact should be minimal. Shape `"coloredrectangle"` actually calls
`graphics::symbols()` which avoids this issue.
However, cases where outer border is expected to be drawn after
the inner border, it has small risk of rendering the fill color
fully opaque, covering the inner border.
One day I probably need to replace `graphics::polygon()` with a
custom function `polygon_with_borders()` that can handle inner and
outer border properly.

## updates to existing functions

* `jam_mypie()` which is used to render `igraph` nodes `shape="jampie"`
was modified to use `col="#FFFFFF01"` for color fill of polygon outer
borders.
* `adjust_polygon_border()` examples were modified to show the
effect of using `col=NA` to `col="#FFFFFF01"` on polygon border rendering.

## bug fixes

* `jam_mypie()` and `shape.coloredrectangle.plot()` were updated to force
`stringsAsFactors=FALSE` when creating `data.frame` objects, fixing weird
color glitch in the examples for `shape.jampie.plot()`.

# multienrichjam 0.0.68.900

## Notable release notes

* `igraph` nodes with `shape="pie"` and `shape="jampie"`

   * New attributes `vertex.frame.lwd` and `vertex.pie.lwd`
   to customize the respective line widths of node borders.
   These attributes require using `jam_igraph()` for `shape="pie"`,
   or require using `vertex.shape="jampie"` with `igraph::plot()`.
   * `jam_igraph()` is now recommended as a more complete
   replacement for `igraph::plot.igraph()`.
   * `shape="jampie"` no longer uses `par("lwd")`, which was global
   change that affected all nodes, edges, and plot features.
   * Pie wedges use "inner borders" for each node,
   so adjacent wedge borders will not overlap. Note that inner border
   slightly overlaps the interior node fill color, to maintain
   consistent node sizes.
   * `frame.color` uses "outer borders" for each node, so these borders
   will not overlap inner pie wedge borders. Nodes are slightly adjusted
   smaller to maintain consistent node sizes, by default.
   * `apply_cnet_direction()` has slightly different logic for
   pie and frame colors, and now assigns `pie.lwd` and `frame.lwd`.
   * `removeIgraphBlanks()` no longer changes single-color `shape="pie"`
   nodes to `shape="circle"` by default.
   * **Changes suggested**:
   
      * Prefer `jam_igraph()` over `plot()` or `igraph::plot.igraph()`.
      * Previous use of `par("lwd")` should be removed, and replaced with
      `vertex.pie.lwd` and `frame.lwd`.
      * Previously, single-color `shape="pie"` or `shape="jampie"` nodes
      were changed to `shape="circle"` by `removeIgraphBlanks()`, however
      they now stay `shape="pie"` in order to maintain control of line widths.
      The `igraph` nodes `shape="circle"` do not respect line width `lwd`.
      * Use of `shape="circle"` should be changed to `vertex.shape="pie"`
      or `vertex.shape="jampie"`, this change affects `nodeType="Gene"`
      moreso than `nodeType="Set"`.
   
   * Edges are now clipped to the outer border of nodes, which helps
   when using transparent node fill, or when using edge arrows.
   
      * With transparent nodes, the edge was previously shown connected
      to the center of each node (without clipping); similarly, edge
      arrows connected to the center of each node, effectively invisible
      and certainly not useful. However, cnet plots do not use edge arrows.
      * Most edges appear identical to previous rendering, however the
      control points used during bundling are used to determine where
      the edge connects to each node border. Thus, when edges converge
      on one node, bundled edges appear to "merge together" at one point
      connecting to the node. Previously nodes entered from multiple
      angles, since they connected to the interior center of the node,
      and showed somewhat more spacing based upon the number of edges.
      * Future options may include retaining some edge spacing
      based upon the intermediate edge curvature of each edge, instead of
      using the same control point for all edges, so the
      connection point on the central node will be slightly spaced out.

* `mem2cnet()` and `memIM2cnet()` changed some default sizes.
The new defaults should work well without adjustment.

   * Node sizes are 4x larger than before, because when calling
   downstream functions we found ourselves always scaling nodes 4x larger.
   
      * The defaults: `categorySize=20`, `geneSize=10` are exactly 4x larger.
      * **Changes required**: Change previous `node_factor=4` to
      `node_factor=1` for reproducibility; or when calling `memIM2cnet()`
      use arguments `memIM2cnet(..., categorySize=5, geneSize=2.5)`.
      
   * node label sizes also have new defaults, similarly because we most
   often applied `label_factor=1.3` by default.
   
      * The defaults: `categoryCex=1.2`, `geneCex=0.9` are adjusted from
      previous `categoryCex=0.9`, `geneCex=0.7`.
      * **Changes required**: Change `label_factor` or `label_factor_l` to 1,
      or change the call `memIM2cnet(..., categoryCex=0.9, geneCex=0.7)`.

* igraph defaults were changed when using `jam_igraph()` for plotting:

   * `vertex.label.family` changed from `"serif"` to `"sans"`.
   * `vertex.pie.border` changed to `"grey30"`, the `igraph` default
   values were identical to `vertex.pie.color` and therefore not visible.
   * `vertex.frame.lwd` set to 1, although this value is ignored by
   `igraph::plot.igraph()` since it only uses `par("lwd")` for all lines.


## changes to dependencies

* `arules` was removed as a dependency, as two functions that used
its class `"transactions"` were rewritten to remove the requirement.
The previous version of those functions are internal and renamed
from: `im2list()` to `im2list_dep()`; and `imSigned2list()` to
`imSigned2list_dep()`. The replacement functions should return identical
data.
* `bezier` was added as dependency to generate edge bundling curves,
in head-to-head tests, it out-performed `graphics::xspline()`, and
produced identical results to internal function `ggforce::bezierPath()`.
Since `bezier` has no dependencies, this addition should feel small.
* `jamba` was bumped to version `0.0.88.900` for an important fix to
`mixedSort()`.
* `colorjam` was bumped to version `0.0.23.900` for consistency.

## changes to existing functions

* `adjust_cnet_nodeset()`

   * When nodesets are not matched in the data, it prints a `warning()`
   then returns the input graph without change. Previously it called
   `stop()`, which is problematic for relevant workflows.

* `removeIgraphBlanks()`

   * default argument changed to `pie_to_circle=FALSE` so single-item
   nodes `shape="pie"` are no longer converted to `shape="circle"`, since
   `"jampie"` nodes are rendered much better. Haha.

* `drawEllipse()`

   * now accepts vectorized input and processes accordingly.
   * New examples show the varied features, many of which are not used
   for `igraph` node shape="ellipse", but they could be used in future.

* `shape.jampie.plot()`, `shape.coloredrectangle.plot()`

   * These functions draw custom `igraph` node shapes.
   
      * `shape="jampie"` draws a pie node shape with more features
      than vanilla igraph pie shape. Each pie wedge is color filled
      with `pie.color`, with new optional inner border `pie.border`,
      `pie.lwd`. The pie node overall can have an outer border defined by
      `frame.color` and `frame.lwd`. Also, pie nodes with only one color
      are rendered as a circle, with no tiny internal line.
      * `shape="coloredrectangle"` draws a series of square boxes
      with color fill `coloredrect.color`, and optional inner border
      `coloredrect.border`. The node overall can have an outer border
      defined by `frame.color`, and `frame.lwd`.
      
   * Both functions now use inner and outer borders via
   `adjust_polygon_border()`, so the various borders no longer overlap.
   Previously, pie nodes drew each wedge with default `graphics::polygon()`,
   which allows adjacent borders to overlap 100%.
   * Note that nodes are resized internally so the rendered node size
   is equal across all nodes even when the line widths (lwd) vary.
   * There are some idiosyncracies from calling `graphics::polygon()` to
   render pie wedges, since it does not by default allow vectorized
   plotting of multiple polygons with different line widths.
   Therefore `shape="jampie"` renders line widths in subsets with identical
   line widths for vectorized plotting, which is 10-100x faster than
   plotting each pie node individually as done by `igraph::plot.igraph()`.
   The main potential issue would be seen with partially overlapping nodes,
   with the potential to display inconsistent overlap order.
   * It seems hopeless to evaluate ggraph/tidygraph to render visualizations,
   in part due to visualization and rendering details. Also, ggraph/tidygraph
   does not store nor recognize many `igraph` visualization details in the
   `igraph` object, instead they must be encoded as `ggplot2` visualization
   options and settings. Using that ecosystem would also require creating
   new ggplot2 node geom types, and edge bundling functions.

* `edge_bundle_nodegroups()`

   * This update represents a substantial refactor of logic.
   * Edge bundles can be "invalid", which causes edges to be drawn as
   linear edges between nodes. The criteria were based upon a series of
   test cases which were all so common that they warranted being fixed.
   
      * Bundling usually occurs along the line between two nodegroup
      center points, calculated as mean node coordinates in each nodegroup.
      Sometimes the configuration of the center points, or the line itself,
      cause edge bundling to become unnecessary or ineffective.
      * Co-linear control points: When the edge and control points are
      co-linear (along the same line), the edge is drawn as a line.
      The criteria uses correlation above 0.99 for node and control points
      of each edge.
      This criteria affects indiviual edges, so nodes in the same nodegroup
      may be rendered differently based upon the specific position.
      The problem is clear when one control point appears beyond the path
      between two nodes. For non-linear edges, the path would curve around
      to the far side of the node. For linear edges, it appears as a line
      that extends beyond one node with optional arrow pointing backward.
      The solution effectively draws the same edge, except clips the edge
      at the first node boundary.
      * Both nodegroups contain only one node: with only one node in each
      group, there is nothing to "bundle".
      * Both nodegroups have identical center points, within some small
      tolerance as a small percentage (0.5%) of the overall layout range.
      * When one nodegroup contains only one node, and the other nodegroup
      center point sits inside the node boundary, there is no bundling.
      Note this criteria is dependent upon node size during rendering.
      The problem is the edge spline control point is inside the node,
      so the spline would curve inside the node boundary, then point
      back out to the node border from the inside. Instead, the edge
      is drawn as a straight line to the node boundary from the outside.
      This situation occurs when one nodegroup fully surrounds a
      central node, so edges are drawn directly to the central node.
   
   * Edge bundling "midpoint" represents the position along the line between
   two nodegroups.
   
      * When one nodegroup contains only one node, this line
      is clipped to the node boundary, so the midpoint is defined
      beginning at the outer edge of the node, toward the center of the
      other nodegroup.
      * Note that when both nodegroups contain only one node, edge bundling
      is already invalid (see above).
      * Note than when the other nodegroup center sits inside the node
      boundary, the edge bundling is also invalid.
      * Therefore this situation only occurs when one nodegroup center
      is already outside the border of the single-node nodegroup.
   
   * Edges are properly clipped using the relevant `igraph` shape
   clip function. See `igraph::shapes()`, and `igraph::shapes("circle")$clip`
   for specific examples.
   * Edge labels are rendered along edges as follows:
   Linear edges are encoded with three coordinates: start, middle, end.
   Spline edges are encoded using default `graphics::xspline()` which
   returns 100 points by default.
   Edge labels are placed using the coordinate most distant from the
   start and end node. For linear edges, the middle coordinate is used,
   for splines a point very near the middle of the edge is used.

* `jam_igraph()`

   * This function is intended as an enhanced drop-in replacement for
   `igraph::plot.igraph()`. It was updated to fulfill previously
   un-implemented features, so fulfill the promise of being a replacement.
   * New arguments, formally passed to internal `jam_plot_igraph()`:
   
      * `mark.groups`, `mark.shape`, `mark.col`, `mark.border`, `mark.expand`,
      `mark.lwd`, `mark.lty`, `mark.smooth` - arguments to enable and
      customize the rendering of nodes within clusters or groups.
      * Note the new options: `mark.lwd`, `mark.lty` for each group border;
      `mark.smooth` to control whether the group polygon is smoothed;
      `mark.alpha` to control alpha transparency of fill colors when not
      already defined.
      * `mark.expand` is now expected to be provided as a fraction of plot
      layout range, which is very close to default behavior in `igraph`
      since the default `rescale=TRUE` forced all layout ranges between
      `c(-1, 1)`.
   
   * `edge_bundling` new option `"default"` will try to detect the
   most appropriate bunding method, based upon whether `nodegroups`,
   `mark.groups` are defined, otherwise it chooses `"connections"`.
   * Edge labels are now rendered for straight edges and edge bundled edges.
   * Edge labels can now accept multiple `"edge.family"` values,
   in the unlikely event of multiple fonts on the same plot. This scenario
   will cause an error with `igraph::plot.igraph()`.
   * Edges are properly clipped based upon the `igraph` node shape.
   * Internal adjustments to `node_factor`, `edge_factor`, `node_factor_l`,
   `edge_factor_l`, `label_factor`, `label_factor_l`, `label_dist_factor`,
   and `label_dist_factor_l` were adjusted to be applied more consistently.
   * Undefined layout is now properly calculated dynamically and passed
   to the internal rendering function `jam_plot_igraph()`, so the values
   for `xlim`,`ylim` are now properly calculated.


## new functions

* `make_cnet_test()` to create Cnet plot `igraph` data for testing.
* `adjust_polygon_border()` defines inner and outer borders for polygons.

   * Using inner borders allows adjacent polygons to have their borders
   visible beside each other, without overlap.
   * Using an outer border allows the display of a border around a collection
   of polygons without overlapping their inner borders.
   * Borders can be layered inside or outside existing borders.
   * There are extensive examples showing various combinations of borders.

* `shape.jampie.clip()`, `shape.coloredrectangle.clip()`

   * Invisible to the user, however they are called for igraph shapes
   `"jampie"` and `"coloredrectangle"`, respectively.
   * These functions now properly clip edges to the outer border of each
   node, including optional inner and outer borders.

* `shape.ellipse.clip()`

   * calculates the optional rotation and size of the ellipse and adjusts
   the edge endpoints accordingly.

* `parse_igraph_plot_params()`

   * reproduces an internal `igraph` function, which cannot be called by
   CRAN-approved R packages.

* `default_igraph_values()`

   * reproduces internal `igraph` package data in a function call, for
   CRAN compliance.

* `jam_igraph_arrows()`

   * Mimics and extends internal `igraph:::igraph.Arrows()` for CRAN
   compliance.
   * It also optionally only renders edge arrows, useful when an edge
   bundling function renders edges itself.

* `get_igraph_arrow_mode()`

   * Another mimic of internal `igraph:::i.get.arrow.mode()` for
   CRAN compliance.


# multienrichjam 0.0.67.900

* added packages to Suggests to support new functions for node layout,
and creation of more "correct" alpha hull polygons around points.

   * `alphahull` - best implementation of alpha hull. However, it also
   requires `sp` package, a heavy install which not advised because it
   is being retired in 2023 in favor of `sf`. Slight risk that the
   `alphahull` package is removed from CRAN if it is not updated.

* Added to Depends

   * `sf` - lightweight replacement of `sp` that provides useful
   geometric functions. It is added primarily because it improves
   rendering pie node borders, which are resized by the exact line
   width determined at the time of plotting pie nodes.

## new functions

* `make_point_hull()`

   * This function is in active development, and is not yet used
   in other functions, but will be used in the next version.
   * takes set of points, makes an alpha hull using `alphahull::ashape()`,
   then expands using `sf::st_buffer()`.
   * If `alphahull` is not available, it uses `grDevices::chull()` which
   does not produce the "ideal" shape but has no additional R depedencies.
   It is also identical to output from `igraph` and `sf` packages for
   point hulls, so it has strong precedent.

## changes to existing functions

* `get_cnet_nodeset()`

   * Refactored for greater speed, in test cases with 585 nodes, previous
   output took 1.28 seconds, the new output takes 0.025 seconds,
   50x speed increase. This function is called in numerous places in this
   package, so this improvement will also positively affect all sorts of
   other functions.

* `colors_from_list()`

   * The sort algorithm was improved in cases where the color order were
   ambiguous, but including names in the tiebreak, before using `H,C,L`
   values.
   * Note: This function is used by `reorder_igraph_nodes()` when `colorV`
   is not explicitly provided, by detecting the probable order of colors
   based upon order of colors in multi-color nodes.

* `reorder_igraph_nodes()`

   * The sort order was updated to be more efficient, and to use better
   logic when following `colorV` or when defining `colorV` ad hoc.
   * When sorting by `"color"` the order will be defined by `colorV`
   whenever colors are aligned by `colorV`, otherwise colors are generally
   sorted by hue `"H"` in HCL space.
   * Examples were added to show clearly the different options for ordering,
   including visual examples of the `jam_mypie()` rendering of pie nodes.

* `jam_mypie()`

   * Visual examples are in `reorderIgraphNodes()`, since the `jam_mypie()`
   function is internal to the igraph plotting scheme.
   * This function draws pie nodes in vectorized fashion, properly drawing
   each set of pie polygons per node, then the appropriate borders defined
   by `pie.border` and `frame.color`.
   * `pie.border` has the option to be drawn inside the border of the polygon,
   so that adjacent Venn wedges will have the entire outer border color visible
   without overlapping the adjacent Venn wedge. To enable, set:
   `options("inner_pie_border"=TRUE)`
   * `frame.color` also has the option to be drawn in a manner that does
   not overlap `pie.border`, it will be drawn just outside the `pie.border`.
   * In cases where `frame.color` is not drawn, the `pie.border` radius is
   adjusted to exactly the line width of the `frame.color` border, so nodes
   will always be exactly the same sizes with or without the `frame.border`.
   * In most cases there should either be `pie.border` *or* `frame.color`,
   however it is possible at some point that `frame.color` and `pie.border`
   will both need to be applied, and this function can handle it.
   * Note this process now draws three layers of polygons:
   
      1. each `pie.color` wedge fill color and no border
      2. each `pie.border` wedge border color with no fill
      3. overall `frame.color` border color with no fill


# multienrichjam 0.0.66.900

## changes to existing functions

* changed to remove calls to `matrixStats::rowMaxs()` and
`matrixStats::rowMins()` to use base functions instead.
This change was due to several R crashes that appear to be bugs somewhere
in the upstream packages, that also occurred on Mac OSX, and on linux,
but in both scenarios involved R-3.6.1 which is not likely to gain
support traction with other package authors. Understandable.
* `mem_enrichment_heatmap()` fixed potential bug:

   * When `row_split` is passed via `...` to the underlying
   `ComplexHeatmap::Heatmap()`, it was not aligned to the order of
   `rownames(matrix)` to be displayed in the heatmap, therefore
   the rows were split in the wrong visual order.
   * `row_split` is now a formal argument, and when supplied as
   a vector, the `names(row_split)` are used to align `rownames(matrix)`
   appropriately.

* `reorderIgraphNodes()`

   * new argument `nodesets` to define a subset of nodes for which the
   reordering will be applied, which may be helpful when nodes in
   a nodeset are horizontal or vertical. In future, this option may
   be applied based upon the aspect ratio of nodes in a nodeset,
   or so that the `nodeSortBy` can be defined as a `list` named
   by `nodesets`. It gets complicated.
   * more output when `verbose=TRUE`
   * minor added checks for layout, ensuring matrix input for layout
   will match `rownames(layout)` to `V(g)$name`, just in case the
   order is not identical.

* `spread_igraph_labels()`

   * reverted apparent regression which did not pass `...` to child
   function `reorder_igraph_nodes()`, therefore the `colorV` color
   order was not properly used when called in this manner.

* `rotate_igraph_nodes()`

   * Change default argument from `center="origin"` to `center="median"`.
   Rare change to argument default, justified by the change being fairly
   benign. Also the new default is more consistent with expectations,
   that nodes would be "rotated in place". Practical outcome is the same,
   nodes are rotated exactly as before, but the coordinate range is more
   likely to remain consistent, when input layout is not already centered
   at coordinates `c(0, 0)`.
   * new argument `verbose=FALSE`
   * minor added checks for layout, ensuring matrix input for layout
   will match `rownames(layout)` to `V(g)$name`, just in case the
   order is not identical.

* `rotate_coordinates()`

   * Change default argument from `center="origin"` to `center="median"`.
   Rare change to argument default, justified by the change being fairly
   benign. Also the new default is more consistent with expectations,
   that coordinates would be "rotated in place". Practical outcome is the same,
   points are rotated exactly as before, but the coordinate range is more
   likely to remain consistent, especially when input layout is not
   already centered at coordinates `c(0, 0)`.

# multienrichjam 0.0.65.900

Numerous changes were made to functions in order to improve
the overall Cnet plot layout experience. New functions were
developed offline that focus on layout specific to Cnet plots,
useful for bipartite graphs in general.

## changes to existing functions

* `jam_igraph()` and `jam_plot_igraph()`

   * new argument `plot_grid=FALSE` which optionally plots
   a grey grid in the background, with units equal to "percentage" across
   the layout coordinate range. This option is intended to help when
   manually adjusting node and node_set positions with `nudge_igraph_node()`
   and `adjust_cnet_nodeset()`, both of which take units `x,y` in the form
   of fraction of the overall layout dimensions.

* `nudge_igraph_node()`

   * new argument `nodes_xy` is intended to help enter adjustments when
   many nodes need to be nudged. The `list` is named by node (e.g. by gene),
   and contains x,y coordinate adjustments. It is completely equivalent to
   entering `node`, `x`, `y` as three independent vectors, but may be
   easier to use.
   * argument default changed from `use_grep=TRUE` to `use_grep=FALSE`,
   because nudging a node `"A"` should not also nudge every node that contains
   an `"a"` or `"A"`. This rare change in default argument value seems
   more helpful to avoid erroneous moves by default.

* `layout_with_qfr()`

   * new argument `constrain` which takes a `character` vector of node names,
   then ensures those node coordinate positions are properly configured
   in `constraints` so they do not move during iterative layout.
   * new default behavior is to define `init` using the current graph layout
   stored in `igraph::graph_attr(g, "layout")`, instead of using a random
   circular initial layout.

* `jam_mypie()` internal function `inner_pie_border`

   * new experimental argument `inner_pie_border` is intended to draw
   the pie wedge border on the inside of the pie wedge shape, so it does
   not directly overlap the border of an adjacent pie wedge shape.
   I could not find a polygon function in R that can draw borders on the
   inside edge of the polygon, which is surprising considering the
   considerable effort to draw GIS world maps in R. Any adjacent borders
   are directly overwritten, with no option for borders to be displayed
   side-by-side along the polygon edge. Seems like an opportunity for someone.
   * nonetheless this capability is not yet implemented, but the framework
   is in place to be released soon.

* `mem_legend()`

   * new arguments `do_directional`, `directional_column`, and
   `directional_colors` are intended to add new directional circles
   indicating up- and down-regulation.

## new functions

* `plot_layout_scale()` plots a grey grid background to an igraph plot
indicating percentage units across the range of layout coordinates.


# multienrichjam 0.0.64.900

## changes to existing functions

* `reorderIgraphNodes()` argument `sortAttributes` now includes `"frame.color"`
in order to include `"pie.border"` for pie shape nodes, and `"frame.color"`
for nodes overall.
* `list2im()` and `list2imSigned()`

   * new argument `emptyValue` to control how empty incidence entries should
   appear, either as zero `0`, or as `NA`. This update fixes rare issue
   where missing enrichment P-values were reported as zero `0` instead
   of `1` by default.
   * These functions appear in `venndir` package, however we do not
   want to make this package dependent upon `venndir` just yet, so they
   remain here for now. In `venndir` these functions are named:
   `venndir::list2im_opt()` and `venndir::list2im_value()`.

* `spread_igraph_labels()` now uses `sortAttributes=NULL` default, when
it is NULL it uses defaults from `reorderIgraphNodes()`. Previously there
was inconsistent defaults between the two functions.
* `memIM2cnet()` and `mem2cnet()`

   * new arguments `geneIMdirection`, `enrichIMdirection` are used
   to call new function `apply_cnet_direction()` to colorize node
   borders by default when data is available.

* `mem_plot_folio()` now calls `memIM2cnet()` and no longer edits colors
itself, those steps are performed by `memIM2cnet()`; no longer calls
`removeIgraphBlanks()`.
* `multiEnrichMap()`

   * now applies exact `colorV` colors to `geneIMcolors`
   without calling `colorjam::matrix2heatColors()` since that process
   returned slightly darker colors by default, and for little benefit.
   * new argument `geneHitIM` intended to allow directional hit matrix
   to be supplied, thus enabling other features as described above,
   specifically colorized border on Cnet `igraph` plots.

## new functions

* `apply_cnet_direction()`

   * colorizes node border, pie.border, coloredrect.border based upon
   directionality, using `geneIMdirection` and `enrichIMdirection`
   when available.

* `reorder_igraph_nodes()` is step one in migrating function names
away from camelCase, toward snake_case. Not really a new function, but
a new function name.
* `mem2cnet()` is a similar rename of `memIM2cnet()` which mostly takes
`mem` as input instead of `memIM` anyway.

# multienrichjam 0.0.63.900

## changes to existing functions

* `mem_gene_path_heatmap()`

   * The attribute `"caption"` is formatted more cleanly.
   * New returned attribute `"draw_caption"` which is a function
   that draws the caption in the bottom-right corner of the heatmap,
   mainly because this location is least likely to overlap other
   heatmap labels. The location and style can be customized as needed.

* `reorderIgraphNodes()`

   * argument `sortAttributes` was updated to include `"pie.border"`
   in the default sort order, for future when Cnet nodes also include
   border color with the direction of change.

# multienrichjam 0.0.62.900

## changes

* bumped dependency on `jamba` to 0.0.87.900 to pick up all the
recent updates.
* `mem_enrichment_heatmap()`

   * new argument `cluster_rows` to control row clustering, specifically to
   allow no row clustering.
   * new argument `do_plot=TRUE` to honor `do_plot=FALSE` from `mem_plot_folio()`

* `mem_gene_path_heatmap()`

   * Now generates its own caption that includes relevant clustering and
   filtering parameters used, helpful to reproduce the original result.
   Caption is returned as `attr(hm, "caption")`.

* `mem_plot_folio()`

   * passes `do_plot` to `mem_enrichment_heatmap()`
   * uses `caption` generated by `mem_gene_path_heatmap()` and returns
   both the gene-pathway heatmap `gp_hm` and the caption `gp_hm_caption`.
   The caption is relevant because the parameters define the clusters
   represented in subsequent Cnet plots.

* `importIPAenrichment()`

   * new argument `remove_blank_colnames=TRUE`, same as previous behavior,
   with new option to disable this behavior and keep columns with all
   values in `c(NA, "")`. We found `zScore` is sometimes reported as
   entirely `NA` but may be useful to keep, for consistency with other
   enrichment results.

* `edge_bundle_nodegroups()` no longer includes `...` when calling
`lines()`, which should silence a fair number of harmless but annoying
warning messages. We don't need a warning message for everything, tysm.

## bug fixes

* `reorderIgraphNodes()` somehow broke in a recent `igraph` update,
apparently the return type is no longer coerced to `character` vector,
so needs to be converted directly.

   * `mem_plot_folio()` and `subsetCnetIgraph()` were impacted as well.
   * `mem_plot_folio()` was modified to call `subsetCnetIgraph()`,
   then failing that call, will return the full `cnet_collapsed`
   graph without subset. Messages are printed for review.

* `edge_bundle_nodegroups()` initial work on allowing custom midpoints
within node groups, which would allow defining a custom midpoint
in x,y coordinates, which may or may not be between the two
sets of nodes. Implementation is not in place yet, but in progress.

# multienrichjam 0.0.61.900

## changes to existing functions

* `reorderIgraphNodes()` was updated to improve the consistency of
sorting by different node properties, specifically to allow sorting
by node fill, and node border color(s).
* `shape.jampie.plot()` is the function used to render node `shape="pie"`
by `jam_igraph()`, is an optimized, vectorized method to render nodes
in one shot, rather than drawing each in a `for()` loop by default.

   * Drawing even 20 or more `pie` nodes is substantially faster
   using `shape.jampie.plot()` compared to default `igraph` function.
   * Speed was improved, the method of converting list of polygon coordinates
   to numeric vectors spaced by `NA` was much improved.
   * It can render `pie.border` for each pie wedge of
   each node with `shape="pie"` or `shape="jampie"`. Note that
   the outer line may be covered by subsequent `frame.color`.
   Attribute `pie.border` is expected to be a `list` where
   `lengths(pie.border)` are equal to `lengths(pie)`.
   * It can render `frame.color` around the full circle of
   each node with `shape="pie"` or `shape="jampie"`. This line
   may cover the `pie.border` if also drawn. Attribute `frame.color`
   is expected to have length equal to `igraph::vcount(g)` which
   is the total number of nodes, one `frame.color` value per node.
   * It is recommended to use one style or the other for each node:
   
      1. `pie.border=NA`, and `frame.color="red"`
      2. `pie.border=c("red", "gold")`, and `frame.color=NA`

* Internal function `jam_mypie()` is called by `shape.jampie.plot()`,
and was updated to handle `frame.color`.


# multienrichjam 0.0.60.900

## changes to existing functions

* `mem_plot_folio()` default arguments changed:

   * `min_set_ct=1`, previously was 2
   * `min_gene_ct=1`, previously was 2

## minor bug fixes

* `multiEnrichMap()` was not verifying `names(geneHitList)` also
matched `names(enrichList)`, therefore a mismatch could cause downstream
errors in `mem_plot_folio()` related functions. This issue now only
uses `geneHitList` with names that match.


# multienrichjam 0.0.59.900

* added `amap` package as dependency, it provides `amap::hcluster()`.

## changes to existing functions

* `mem_plot_folio()` verbosity was reduced when `verbose=FALSE`.

## bug fixes, enhancements

* `mem_gene_pathway_heatmap()` now honors `p_floor` when defining the
incidence matrix values to be used in clustering the weighted and
combined enrichment and heatmap gene-pathway incidence matrix data.
* `mem_gene_pathway_heatmap()` was throwing an error when supplying
`column_split` and the default `cluster_columns=TRUE`.

   * The error was caused by creating a dendrogram for `cluster_columns`
   then supplying `ComplexHeatmap::Heatmap()` with a `character` vector
   `column_split`, and a dendrogram/hclust for `cluster_columns`. Instead
   it allows passing a `function` to `cluster_columns`, which also requires
   using custom data, since the data used for clustering is a weighted
   combination of the enrichment P-values across the top of the heatmap,
   and the data inside the heatmap.
   * The new default when `row_split` or `column_split` are `character`
   will be to define `cluster_rows` or `cluster_columns`, respectively,
   to a `function` that calls `amap::hcluster()` on the combined and
   weighted heatmap and respective annotation data matrices.

* `mem_plot_folio()` was not properly defining gene row cluster names
for `mem_gene_pathway_heatmap()`, leaving them empty by default instead
of assigning from `letters`.


# multienrichjam 0.0.58.900

* bumped dependency to `jamba (>= 0.0.84.900)` to retire `call_fn_ellipsis()`

## changes to existing functions

* `topEnrichListBySource()` and `topEnrichBySource()` were updated:

   * `nameColname` is handled properly, without relying upon `"Name"` colname,
   and without relying upon `rownames()` of the enrichment `data.frame`.
   * Now the subset operations use values in the `nameColname`.
   * Also the rows in each enrichment `data.frame` with values in
   `nameColname` are subset using `subset()`, which means in some rare
   cases multiple rows might be returned, if the input enrichment data
   has the same name in `nameColname` for multiple rows. This change
   is intentional, in order to retain all rows with matching names,
   not just the first that occurs.

* `enrichList2geneHitList()`:

   * argument `geneDelim` has a default value, consistent with `multiEnrichMap()`

## bug fixes

* `enrichList2IM()` will return data with zero rows without error, although
it usually occurs because of an error somewhere else (for example input data).


## functions removed

* `call_fn_ellipsis()` was moved to the jamba package version `0.0.84.900`.


# multienrichjam 0.0.57.900

* Issue #6 reported an error when using a series of enrichments
where only a subset contain a z-score column name for `directionColname`.
Related, the current approach ignored the `directionColname` when
all values were NA. Both situations have been corrected, to allow
flexible mish-mash of NA and non-NA values, and presence/absence
of `directionColname` in each enrichment input.
Also related, `mem_plot_folio()` was by default no enabling the
directionality via `mem_enrichment_heatmap(mem, apply_direction=TRUE)`,
therefore there is a new argument to `mem_plot_folio()` default
`apply_direction=NULL` which will auto-detect whether there is
directional non-zero and non-NA values that can be used in the heatmap.
Also, `apply_direction` can be defined on its own to force the issue.


# multienrichjam 0.0.56.900

## bug fixes

* Issue #7 reported an error, traced back to the vignette. The error
was caused by passing arguments to `mem_plot_folio()` by overloading
`...`, when not all arguments were valid in `ComplexHeatmap::Heatmap()`.
The fix:

   * new function `call_fn_ellipsis()` which passes arguments including `...`
   to another function, and when that function arguments do not allow `...`
   then it limits the arguments in `...` to those arguments accepted
   by the other function.
   * Instead of: `x <- some_function(a=1, b=2, ...)`
   * Use: `x <- call_fn_ellipsis(some_function, a=1, b=2, ...)`
   * `mem_enrichment_heatmap()` and `mem_gene_pathway_heatmap()`
   were updated to use `call_fn_ellipsis()`.

* Another error was noted during the vignette workflow, that the
`directionColname` was being populated even when no enrichment data
contained non-NA values, which was inconsistent with `find_colname()`
in function `enrich2IM()`. This situation was corrected by requiring
only one enrichment result to contain a non-NA value in this column.
* `enrichList2df()`, `multiEnrichMap()`, `enrichList2IM()` were
updated to call `find_colname()` with the list of `enrichResult` objects.

## changes to existing functions

* `find_colname()` can accept a `list` object, which is expected to
contain a list of `data.frame` and/or `enrichResult` objects. An
`enrichResult` is converted to `data.frame` by `enrichResult@result`,
all other objects must contain `colnames(x)`. When `require_non_na=TRUE`
it will test each object, and return `max` unique entries that match
`pattern`. The entries should all contain the same matching colnames
for most purposes in `multienrichjam`.


# multienrichjam 0.0.55.900

## bug fixes

* `multiEnrichMap()` was incorrectly populating empty gene counts
with default `1` instead of `0`. The effect is mainly during filtering
by gene count, where the minimum is usually never below `1`, however
it can cause issues in point sizing particularly in
`mem_enrichment_heatmap()`. The previous default `1` is a remnant of
using this function to generate a matrix of enrichment P-values.


## changes to existing functions

* `enrichList2IM()` argument default was changed to
`emptyValue=NA` so in the absence of data to populate into
the incidence matrix, the default cell value will be `NA`
to indicate there is no available data.
* `multiEnrichMap()` was updated to supply a specific `emptyValue`
for all calls to `enrichList2IM()`.
* `cell_fun_bivariate()` new argument `type` is intended to allow
re-using this same function for univariate color functions, so
heatmaps features can be made consistently, specifically for
dotplot or normal heatmap output, and optionally labeling cells
with statistical values.
* `mem_enrichment_heatmap()` is updated to share common heatmap
code for bivariate and univariate color gradients. One by-product
is that output cannot be raster format, which is typically not
an issue for pathway enrichment, since pathways should not represent
more than 1,000 or so pathways. In that event, output should probably
be rasterized (PNG, JPG) instead of vector graphics
(PDF, SVG).
* `make_legend_bivariate()` new argument `digits` used to prevent
displaying weird labels like `5.9999999998` and instead will display `6`.
* `jam_igraph()` default argument was changed to `edge_bundling="connections"`
which will enable edge bundling by default. It can be disabled with
`edge_bundling="none"`, although it is only active when there are edge
connections that can be bundled.


# multienrichjam 0.0.54.900

## updates to existing functions

* `mem_enrichment_heatmap()` was updated to customize directional
information in heatmaps. Cell labels are optionally displayed,
defined by `show`, to display one or more of `z-score`,
`-log10pvalue`, and `gene count`. Argument `sets` was updated
to handle presence of `enrichIM` and `enrichIMgeneCount` and
`enrichIMdirection` - as well as future measurements with prefix
`enrichIM`. New argument `min_count` applies the gene count
filter to the dot plot output heatmap.


# multienrichjam 0.0.53.900

This update mainly focuses on implementing directionality when
available in the pathway enrichment data during `multiEnrichMap()`.

## new functions

* `colorRamp2D()` - implements bivariate color scale, in this case
for enrichment `-log10pvalue` for color intensity, and `z-score`
for color hue directionality: blue "inhibited", gold "neutral",
and red "activated".
* `display_colorRamp2D()` - display the color ramp defined by `colorRamp2D()`
* `cell_fun_bivariate()` - define a heatmap cell function that uses
the color defined by `colorRamp2D()`, optionally plots circle points,
and optional text labels
* `make_legend_bivariate()` - creates a 2-D color ramp legend suitable
for `ComplexHeatmap::draw(..., annotation_legend_list=x)`

## changes to existing functions

* `mem_enrichment_heatmap()` new argument `apply_direction=TRUE` will
enable bivariate colors, where color intensity is defined by the
enrichment -log10 P-value, and color hue is defined by z-score direction:

   * blue = "inhibition" with z-score <= -2
   * gold = "no direction" with z-score > -2 and x score < 2
   * red = "activation" with z-score >= 2

Points are sized by number of genes.



# multienrichjam 0.0.52.900

## changes to existing functions

Moved some functions into their own .R file for better organization:

* `topEnrichListBySource()` and `topEnrichBySource()` new arguments
`directionColname` to define an optional column name that contains
`numeric` values indicating direction of pathway enrichment. These
values are often in the form of an `"Activation z-score"`, as is the
case with IPA "Upstream Regulators". Argument `direction_cutoff`
refers to the absolute value required for direction be given
a "sign" up or down.

* `multiEnrichMap()` new argument `directionColname` used to determine
directionality of pathway enrichment, useful for things like
`"Activation z-score"`, as returned by IPA "Upstream Regulators".
Output included `mem$enrichIMdirection` which contains the `numeric`
values, where `NA` values are substituted with `0` zero.
These values will be used in near future, likely in
`mem_enrichment_heatmap()` to indicate predicted direction
of impact on particular pathways.
Argument `direction_cutoff` is passed to `topEnrichListBySource()`
for optional filtering to require at least one pathway to contain
an absolute direction score at or above this threshold. The IPA
z-score recommends a threshold z=score >= 2 for "activation" or
"inhibition". Note that many pathways have no z-score, so applying
this threshold will remove those pathways from downstream analysis.

* `mem_enrichment_heatmap()` new argument `apply_direction` and
`direction_cutoff` determine whether to indicate direction if it exists,
and optionally applies a threshold. Still in testing currently.
Note that applying this cutoff will hide pathways whose `numeric`
direction is below the threshold.
New argument `gene_count_max` to apply a max gene count threshold
for the point size for the dot plot format.
New argument `legend_height` to control the heatmap legend color bar
height.


# multienrichjam 0.0.51.900

## changes to existing functions

* `mem_gene_path_heatmap()` will now hide all but the main
enrichment colors in color legends when there are more than
8 combinations, to prevent the legend from taking the entire
plot device space and not displaying the heatmap. The threshold
is configurable with `show_heatmap_legend=8`.
New arguments `show_gene_legend`, `show_pathway_legend` are `logical`
and are intended to allow hiding the other color legends. A minimalist
style is to show only the main enrichment colors, which are used for
all other colors anyway.
* `mem_enrichment_heatmap()` was updated so the default dotplot format
has properly controlled point sizes, and legend point sizes.
Previously the two sizes were independent and required manual
adjustment. The current approach ensures the legend point size
exactly match the heatmap dot plot point size. An optional parameter
`cexCellnote` will display labels unless `cexCellnote=0` in which
case gene count labels are hidden.



# multienrichjam 0.0.50.900

## changes to existing functions

* `mem_legend()` new argument `inset` passed to `legend()`.
Also this function uses `tryCatch()` to try to pass `...`
arguments, and if they fail it tries again without `...`.
Fun.


# multienrichjam 0.0.49.900

## extended function help text

* `jam_igraph()` is a drop-in replacement for `igraph::plot.igraph()`,
and its options were described in much more detail. In brief:

* it plots nodes vectorized which is substantially faster
* it changes the default to `rescale=FALSE` and maintains
aspect ratio 1:1, so layout coordinates are rendered without
distorting the x- and y-axis ranges
* it optionally bundles edge connections, which helps particularly
with large bipartite graphs (especially gene-to-pathway graphs.)
* it allows bulk adjustment to node size, label size, and label
distance from node center.

## changes to existing functions

* `mem_plot_folio()` new argument `edge_bundling="connections"`
that by default will bundle edges appropriate for Cnet plots.
Disable with `edge_bundling="none"`.
* `mem_enrichment_heatmap()` new argument `style="dotplot"` will
create a dotplot styled heatmap, whose points are sized proportional
to the number of genes involved in enrichment. This style is
under development and may require additional customization options
such as setting a max gene count for point sizes.


## bug fixes

* `mem_gene_path_heatmap()` fixed bug where enrichment names were
mangled by `data.frame(...)`, avoided by using
`data.frame(check.names=FALSE, ...)`.

Several functions were updated to include proper package prefix
for function calls:

* `matrixStats::colMins()`
* `igraph::V()`
* `igraph::E()`
* `igraph::degree()`
* `igraph::components()`
* `igraph::vcount()`
* `igraph::neighbors()`
* `igraph::vertex_attr_names()`
* `igraph::list.graph.attributes()`
* `igraph::ego()`
* `igraph::set_graph_attr()`
* `igraph::graph_attr()`

affecting functions:

* `subsetCnetIgraph()`
* `apply_nodeset_spacing()`
* `removeIgraphBlanks()`
* `removeIgraphSinglets()`
* `mem_find_overlap()`
* `spread_igraph_labels()`
* `subgraph_jam()`



# multienrichjam 0.0.48.900

## bug fixes

* `mem_gene_path_heatmap()` was updated to handle edge cases:

   * `gene_im_weight` 0 or 1; `enrich_im_weight` 0 or 1
   * user-supplied `cluster_rows` or `cluster_columns` with `colorize_by_gene`
   * The `colorize_by_gene` logic was changed to convert to an integer
   matrix that refers to colors and labels by factor levels, which helps
   for user-defined clustering methods, also helps labeling the color
   legend.
   * new argument `colramp` for user-defined color gradient.
   * `row_title` is no longer assigned when user-defined argument is `NULL`.
   This workaround helps the edge case with `gene_im_weight=1` where
   the `row_split` is forced to have a limited number of values regardless
   what `row_split` integer value is sent, thus causing `row_title`
   mismatch in length.
   * Fixed egregious type in second part of an `if` statement, only
   called when user supplies a custom cluster function.
   

# multienrichjam 0.0.47.900

## bug fixes

* `topEnrichBySource()` was not correctly handling multiple
`sourceColnames` values, instead was only using the first value.
This bug has been corrected.

## changes to existing functions

* * `topEnrichBySource()` and `topEnrichListBySource()` now
return `enrichResult` when supplied with `enrichResult`,
instead of coercing to `data.frame` which would then need
to be converted back to `enrichResult`.
* `topEnrichBySource()` and `topEnrichListBySource()` argument
default `sourceColnames` was changed to reflect
usage with R package `msigdbr` for MSigDB gene set data,
specifically `sourceColnames=c("gs_cat", "gs_subcat")`.
Similarly, default values were removed from `curateFrom`
and `curateTo`, since these defaults imposed a specific
outcome.
* `multiEnrichMap()` default values for `topEnrichCurate*` arguments
were changed to NULL; default argument values changed to
`topEnrichSources=c("gs_cat", "gs_subat")`.


# multienrichjam 0.0.46.900

## changes to existing functions

* `mem_gene_path_heatmap()` new argument `rotate_heatmap=TRUE`
will rotate the heatmap layout so pathway names are displayed
as rows, and genes are displayed as columns.

I still need example data to use for function document examples,
to show utility of rotating the gene-pathway heatmap.


# multienrichjam 0.0.45.900

## changes to existing functions

* `color_edges_by_nodes()` was implemented twice (smh) so the
older function was renamed `color_edges_by_nodes_deprecated()`.
The older function blended colors using a simpler approach
with `avg_colors_by_list()` that took a very fast hue average;
while the new function uses `colorjam::blend_colors()` that
uses red-yellow-blue additive color blending model.
Honestly the function `color_edges_by_nodes()` is not used much yet,
but is likely to be used more with community detection,
and Cnet/community edge bundling.
The new `color_edges_by_nodes()` also applies alpha to the
intermediate blended colors, to retain the relative weight
of colors during the blending step.


# multienrichjam 0.0.44.900

The edge bundling update!

The new functions are under rapid development, and represent potentially
quite useful functions for visualizing complex `igraph` networks,
specifically Cnet plots with Gene and Set node types.

The general technique of bundling edges between two groups
of nodes is relatively stable. The details and extent that edges are
bundled, including visual display of bundled edges, is
under active evaluation and development. So far it doesn't
seem to take much to improve the output figure, compared
to using straight edges.

In the near future, Cnet plots may by default enable
some form of edge bundling, as Cnet plots were in fact
the motivating example.

See pkgdown docs for `edge_bundle_nodegroups()` for visual
examples using the Karate network.

## current issues in edge bundling

`edge_bundle_nodegroups()` is the core function, and it does not
yet implement:

* edge clipping based upon igraph vertex shape, for example
using `igraph:::.igraph.shapes[["circle"]]$clip()`.
* edge arrows, which requires edge clipping for proper usage,
otherwise edge arrows will be underneath the node shape itself.
* edge labels, which could be implemented, except that bundled
edges are by definition much more likely to be too close for
effective labels. I rarely use edge labels, so I will target
the simplest thing that works.


## new igraph edge bundling functions

* `get_bipartite_nodeset()` - is a general use function to return
nodes that are grouped by having identical node neighbors. This
situation usually happens rarely, except with bipartite graphs
where there are often clusters of nodes with the same neighbors,
especially common for Cnet plot data.
* `edge_bundle_bipartite()` is the first edge bundling function,
it is actually a light wrapper around `edge_bundle_nodegroups()`.
Bipartite bundling connects a nodeset (defined as having the
same neighbor nodes) to each neighbor node.
* `edge_bundle_nodegroups()` is a general edge bundling technique
that bundles edges between node groups. These node groups can
be defined by any relevant technique, typically a
community detection algorithm such as `igraph::cluster_walktrap()`
or any other of several `igraph::cluster_*` functions.

## new igraph functions

* `color_edges_by_nodes()` is a simple function that blends the colors
of the two nodes involved in each edge, using `colorjam::blend_colors()`
since it uses RYB color blending. It is optimized so it only
blends unique combinations of colors.


## new utility functions

* `deconcat_df2()` is a necessary utility function that "expands"
a `data.frame` that has one or more columns with multiple
delimited values. It simply expands the rows to represent one
individual value per row for those columns. This function
will very likely be moved into the `"jamba"` package for
much broader use.
* `handle_igraph_param_list()` is a helper function for `jam_igraph()`
which is intended to update node and edge attributes in bulk based
upon another attribute. For example "make all the Gene nodes small,
and all the Set nodes large." Same with labels, colors, etc.


## updates to existing functions

* `jam_igraph()` has new arguments:

* `edge_bundling` - for `"nodegroups"`, `"connections"` and `"none"`
This argument will enable edge bundling by calling `edge_bundle_nodegroups()`
and `get_bipartite_nodeset()` when needed. See examples.
* `render_nodes,render_edges,render_nodelabels,render_nodegroups` are
options to skip various aspects of `igraph` plot features, either to
save time, or to allow more detailed visual layering.
* `vectorized_node_shapes` toggle the vectorized node shape plotting
when there is more than one node shape in the `igraph` object. Note
this feature is *substantially* faster, but changes the order that
nodes are rendered, by design. As a result, nodes are drawn in order
of shape, so each render is a bulk operation. For nodes that overlap,
or partially overlap, this feature will visibly change the ordering
of nodes. In general, speed is still ideal, and reducing node layout
overlaps should be a separate step.


# multienrichjam 0.0.43.900

## changes to existing functions

* `adjust_cnet_nodeset()` logic was updated for `expand`
which was not properly applying negative values to
compress the spacing between nodes in a node set. This
use is relatively rare but feasible.
* `adjust_cnet_nodeset()` argument `set_nodes` now optionally
allows referring to a nodeset using the name of one node
contained in the nodeset, as a convenience.

## new functions

* `apply_nodeset_spacing()` attempts to automate the process
of applying `adjust_cnet_nodeset()` for each subcluster of
nodes in a Cnet plot to enforce a minimum spacing between
nodes. The Fruchterman-Reingold layout algorithm
might offer a minimum distance threshold, but I could
not find it. Incidentally, this new function can also
compress node spacing, which might be useful when
gene labels are not shown.



# multienrichjam 0.0.42.900

## changes to existing functions

* `mem_plot_folio()` has much more detailed help text,
describing more details about clustering parameters,
and describing the returned objects in detail.
* `mem_plot_folio()` new argument `do_plot=TRUE` determines
whether each plot is rendered, or just returned as a plot
object to be reviewed separately.
* `mem_plot_folio()` the returned data now
includes the gene-pathway heatmap caption, which descibes
the gene/set filtering criteria, and the row/column
distance methods used for clustering. Since these arguments
have substantial effect on the pathway clusters, it is
helpful to keep this information readily available.
Similarly, the `cnet` plots as `igraph` objects store
the title as a graph attribute, accessible using
`graph_attr(cnet, "title")`.
* `mem_plot_folio()` argument `colorize_by_gene=TRUE` now
uses `colorjam::blend_colors()`.
* `mem_gene_pathway_heatmap()` is tolerant of `'...'` entries
not valid with `ComplexHeatmap::Heatmap()`, in order to
allow overloading `'...'` for other function arguments.


## new functions

* `nudge_igraph_node()` is intended to help move individual
nodes in an `igraph` layout, useful for adjusting nodes
to reduce label overlaps.
* `get_cnet_nodeset()` is a helper function to get the nodes
in a nodeset, defined as the `"Gene"` nodes that all connect
to the same `"Set"` nodes. It is helpful when looking at
a Cnet plot, and wanting easy access to the `"Gene"` nodes
in a particular cluster of nodes in the plot.
* `adjust_cnet_nodeset()` is intended to manipulate the
layout coordinates for nodes in a nodeset, useful to
expand, shift, or rotate nodes to help with visibility.
* `adjust_cnet_set_relayout_gene()` is a very useful function
with complicated name. It is intended to help move `Set`
nodes in a Cnet `igraph` layout, then re-positions all
the `Gene` nodes while keeping the `Set` nodes in fixed
positions.
* `rotate_igraph_layout()` which allows rotating an `igraph`
layout. It can also reflect layout coordinates across one
or more axes. It also calls other helper functions
as needed, `spread_igraph_labels()` and `reorderIgraphNodes()`.
* `rotate_coordinates()` is the underlying function to
rotate coordinates, it operates on a numeric matrix
and is called by `rotate_igraph_layout()`.

The next version should include `plot_cnet_heatmaps()` which
creates nice Cnet cluster plots where each cluster has a
corresponding expression heatmap displayed at the edge
of the figure.


# multienrichjam 0.0.41.900

## bug fixes

* `mem_plot_folio()` returns `NULL` when there is an error
during the Cnet collapse step, and does not report the error.
Now it at least reports the error. The underlying cause in
this case was ComplexHeatmap not accepting horizontal
color legend orientation in R-3.6.1 (grid package version 3.6.1)
because of some new unit arithmetic only available in
grid 4.0.0+. Hiding the color legend, or using vertical
color legend fixed the issue.

# multienrichjam 0.0.40.900

## bug fixes

* `collapse_mem_clusters()` - fixed error when using `apply()`
on a matrix that had one row or one column.
* `subsetCnetIgraph()` fixed issue when the Cnet igraph object
has no internal layout, it runs `relayout_with_qfr()` by default.

## changes to existing functions

* `mem_plot_folio()` now returns `clusters_mem` in the `list`,
which includes the pathway set names represented in each
cluster shown in the gene-pathway incidence matrix heatmap.

# multienrichjam 0.0.39.900

## changes to existing functions

* `cnetplot_internalJam()` was updated to include `"nodeType"` as
a node attribute, which helps distinguish `"Gene"` and `"Set"`
nodes in downstream operations. This change helps address
#5 to hide the gene labels on Cnet plots.

# multienrichjam 0.0.38.900

## bug fixes

* Update to address issue #4 in `enrichDF2enrichResult()` with
new argument `descriptionColname` which will force the resulting
colname to be `"Description"` to fit expectations of
`enrichplot:::fortify.internal()` which requires `df$Description`.
When `descriptionColname` is not supplied, or not found in the
input `enrichDF` a warning is issued that describes the problem.

## enhancements

* `heatmap_row_order()` and `heatmap_column_order()` now also
work with `HeatmapList` objects. By default they use the first
heatmap in the list, which should be consistent with all other
heatmaps.

# multienrichjam 0.0.37.900

## bug fixed

* `multiEnrichMap()` argument `subsetSets` was never implemented;
existing arguments `descriptionGrep` and `nameGrep` were used for
similar but insufficient purpose. The `subsetSets` argument defines
specific pathway names to retain for analysis, and is implemented
through `topEnrichBySource()` and by proxy `topEnrichListBySource()`.
These functions might better be called `subsetEnrichResult()` and
`subsetEnrichList()`, but will not rename these functions.

## enhancements

* `mem_plot_folio()` includes the number of rows (genes)
and columns (pathways) displayed in the gene-pathway incidence
matrix heatmap.

# multienrichjam 0.0.36.900

## bug fixed

* `enrichDF2enrichResult()` was updated to fix #3, thanks to
@john-lee-johnson for reporting. Issue arose because
`enrichplot::cnetplot()` expected enrichResult rownames to
be equal to values in the `"ID"` column of the enrichment
`data.frame`.

## changes to existing functions

* `multiEnrichMap()` new argument `min_count` requires a pathway
to contain at least `min_count` genes in order to be considered
a "hit". This filter is mostly important when used with `topEnrichN`
to use the top pathways -- it will only sort pathways then take
the top `topEnrichN` number of pathways that also contain at
least `min_count` genes.
* `topEnrichBySource()` and `topEnrichListBySource()` new arguments
`min_count`, `p_cutoff` require pathways to contain at least
`min_count` genes, and have no higher than `p_cutoff` enrichment
P-Value. Previously these functions only took the top pathways,
regardless of these filters (which were applied later). This change
allows the filters to be applied before taking the top `topEnrichN`
pathways, which mostly helps when `min_count` is greater than 1 --
sometimes pathways with only one gene involved in enrichment
have statistically significant P-value, but are not biologically
relevant for interpretation or follow-up experiments.
For `topEnrichListBySource()` the filter is applied to each
enrichment in the list, and any pathways meeting the criteria are
taken for all enrichment lists. So a pathway must 
be present in the top `topEnrichN` entries which meet both
the `p_cutoff` and `min_count` criteria to be retained by these
functions.

# multienrichjam 0.0.35.900

## bug fixes

* `mem_gene_path_heatmap()` was updated to fix a small issue
with `min_set_ct_each` which requires at least one enrichment
to have `min_set_ct_each` genes. However, this filter was
not applied alongside the `p_cutoff` -- therefore some pathways
with enough genes which were not significantly enriched were
fulfilling these criteria, as along as another enrichment was
significant. The new behavior (as expected) is to requires
an pathway to meet both the `min_set_ct_each` and `p_cutoff`
thresholds in the same enrichment in order to be retained
in the gene-pathway incidence matrix.
* `mem_gene_path_heatmap()` fixed edge case where genes present in
multiple enrichments counted more toward `min_set_ct_each`,
but should only count once per gene.


## changes to existing functions

* `mem_gene_path_heatmap()` new arguments `column_title`, and
`row_title` allow custom cluster names, which are also
carried to Cnet cluster names as needed.
* `mem_gene_path_heatmap()` was updated to handle `row_split`
and `column_split` arguments more intuitively, and to allow
`row_split` values `FALSE`,`0`,`1` to inactivate the split
completely. Previously it was not possible to turn off split.
* `mem_plot_folio()` argument was renamed from `pathway_row_split`
to `gene_row_split` to reflect the intent of this argument
more accurately. This change is early in the function lifecycle
and not in broader use yet -- better to change it now. Otherwise,
future argument name changes will not occur without some type
of backward compatibility.
* `mem_gene_path_heatmap()`, `mem_enrichment_heatmap()` and
`mem_plot_folio()` new arguments `row_cex` and `column_cex`
used to adjust row and column heatmap labels, which is helpful
when used with auto-sized axis labels to make minor adjustments.
* `mem_enrichment_heatmap()` was updated to make the color ramp
more consistent with `mem_gene_path_heatmap()`, so the color
scale atop the gene-pathway heatmap more accurately reflects the
color scale used in the enrichment heatmap itself.

# multienrichjam 0.0.34.900

## changes to existing functions

* `mem_gene_path_heatmap()` new argument `colorize_by_gene=TRUE`
will color the heatmap body using blended colors from the
`geneIMcolors` which represents the enrichments in which the
gene is involved. The default `colorize_by_gene=FALSE`
instead colors the heatmap body by the number of enrichments,
which can be confusing if one of the enrichment colors is
also red. The goal is to make it visually apparent when
a gene is involved in one enrichment, by using the color
from that enrichment. This feature is still in development
and testing.
* `mem_gene_path_heatmap()` changed argument name from `enrich_gene_weight`
to `enrich_im_weight`, before this argument is in wider use. Added
new argument `gene_im_weight`. These arguments more accurately reflect
the relative weight between enrichment and incidence matrix for
`enrich_im_weight`.

# multienrichjam 0.0.33.900

## changes to existing functions

* `collapse_mem_clusters()` was updated with new argument
`cluster_color_min_fraction` to help filter the enrichment
colors to include in each resulting cnet cluster. The intent
is not to represent colors where the number of significant
pathways is below this threshold. For example, a cluster of
10 pathways may have only one significant pathway for a
given enrichment set -- therefore this enrichment color
would not be included in the cnet cluster colors.

# multienrichjam 0.0.32.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` cleaned up the heatmap overall:

   * left and top annotations have zero gap between columns and rows
   * left and top annotation legends have discrete color bars; for
   left it is always `c(0,1)`; for the top it uses -log10 integer steps.
   * top annotation legend appends `"-log10P"` to the label, to indicate
   that the color values are based upon `log10()` transform of the
   enrichment P-value.
   * the color legend has discrete color bar steps, indicating the
   number of enrichments where each gene is involved

# multienrichjam 0.0.31.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` new argument `enrich_gene_weight`
used to adjust the relative influence of the enrichment `-log10 P-value`
and the gene incidence matrix on the column clustering. The
effect is to adjust how much the enrichment P-values, or the
gene content, affects the clusters. In principle both should
have similar effects, but sometimes it helps to favor gene
incidence or pathway enrichment.

# multienrichjam 0.0.30.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` was altered to allow filtering
pathways for a minimum number of genes represented by at least
one enrichment result. For example, a pathway may be required to
have at least 4 genes, despite having a statistically significant
enrichment P-value. The new argument `min_set_ct_each` tests
each enrichment to see if any one has enough genes per pathway.
There are two main effects of this filter: The gene-pathway
heatmap will display fewer pathways; and any resulting
Cnet plots will have the same pathways removed.

# multienrichjam 0.0.29.900

## changes to existing functions

* `mem_plot_folio()` new argument `byCols` which is passed along
to `rank_mem_clusters()`, and is used to sort the pathway names
within each Cnet cluster. The default uses `"composite_rank"`,
which sorts by the `floor(-log10(pvalue))` then by the highest
number of genes per pathway. The `floor()` function effectively
sorts by the order of magnitude of the enrichment P-value,
dropping the details. The obvious alternative is `"minp_rank"`
which uses the `-log10(pvalue)` directly, which therefore does
not effectively sort by gene count since enrichment P-values
rarely tie. In our experience, between two pathways with
reasonably similar enrichment P-value (within one order of magnitude
such as 2.4e-5 and 3.5e-5) the pathway with more genes was
usually the more interesting/relevant biological pathway.

# multienrichjam 0.0.28.900

## new functions

`im2list()` and `imSigned2list()` are the reciprocal functions
to `list2im()` and `list2imSigned()`. The new functions convert
incidence matrix to list, or signed incidence matrix to list,
respectively. They're notable because they're blazing fast
in our testing, thanks to efficient methods from the `arules`
R package for interconverting list to compressed logical
matrix, and vice versa.

## changes to existing functions

* `mem_enrichment_heatmap()` was modified to use `p_floor`
as the ceiling for the row hierarchical clustering, previously
it used `ceiling=3` which caused the dendrogram to have `height=0`
when all enrichment results were lower than 0.001.
* `mem_enrichment_heatmap()` was modified to calculate its own
matrix colors when `color_by_column=TRUE`, until the
`colorjam::matrix2heatColors()` can be updated to the improved
method.
* `mem_plot_folio()` now correctly honors `p_floor`.

# multienrichjam 0.0.27.900

Added a new TODO.md file to track some new feature ideas.

## new functions

* `colors_from_list()` infers the proper order of colors from
a list of color vectors. It is called by `reorderIgraphNodes()`
where attributes `"pie.color"` and `"coloredrect.color"` contain
subsets of colors, but in the proper order. This function returns
all colors in their proper order.
* `colors_from_list()` takes a list of colors, and returns the
unique colors, ordered by the overall order inferred from the
list. It is internally called by `reorderIgraphNodes()` when
argument `colorV` is not suppleid, since the attribute `"pie.color"`
can be used to infer the correct order of colors.

## Changes to existing functions

* `collapse_mem_clusters()` now uses `"; "` as delimiter for the
`"set_names"` attribute, to make it easier to distinguish individual
pathway names.
* `reorderIgraphNodes()` was refactored to handle specific color
order defined either by argument `colorV` or by inferring the correct
order of colors from attributes such as `"pie.color"`. This change
should help order node colors more consistent to the original input
colors to `multiEnrichMap()`, specifically the argument `colorV`.
* `reorderIgraphNodes()` now has an example that shows the effects
of ordering pie nodes by color.

# multienrichjam 0.0.26.900

## Changes to existing functions

* `collapse_mem_clusters()` now includes the cluster name
by default with the `set_labels` which are used to display
the top `n` pathway names for each collapsed Cnet node.
For example `"Cluster A: Aryl Hydrocarbons; Granzyme Signaling"`.
* `collapse_mem_clusters()` arguments `max_labels` and `max_char_labels`
accept multiple values, which are recycled and applied to filter each
cluster. For example `max_labels=c(5,2)` will filter the first
cluster to display up to 2 labels, and the second cluster to display
up to 5 labels.

# multienrichjam 0.0.25.900

## Changes to existing functions

* Thanks to `simang5c` the function calls to `renameColumn()` were
corrected to `jamba::renameColumn()`, to resolve issue #1. Also fixed
several other references to `jamba` functions: `nameVector()`, `nameVectorN()`,
`makeNames()`, `printDebug()`, `mixedSortDF()`, `cPaste()`,
`igrepHas()`, `provigrep()`, `vigrep()`, `unvigrep()`, `rbindList()`,
`noiseFloor()`, `rmNA()`, `rmNULL()`, `normScale()`, `getColorRamp()`,
`deg2rad()`. Omg! There were so many more cases than I thought!
(Package-building options to import functions from packages seem
to go too far, the imported functions appear to be contained
in the new package which seems misleading... Ah well.)
* Also made several small changes to handle single-enrichment
input to multiEnrichMap(). It still provides some useful graphical
benefits even with only one enrichment input.

# multienrichjam 0.0.24.900

## Changes to existing functions

* `multiEnrichMap()` was updated to enable analysis with only one
enrichment result. Testing whether the downstream capabilities are
useful in the context of one enrichment result.

# multienrichjam 0.0.23.900

## Changes to existing functions

* Several functions have additional `jamba::printDebug()` output
when `verbose=TRUE`.
* `jam_plot_igraph()` and `jam_igraph()` new argument
`use_shadowText=TRUE` will enable `jamba::shadowText()` to replace
`graphics::text()` for igraph text, which affects node and edge
labels. The goal is to make labels more widely legible when they
are placed on top of light and dark colors, common when using dark
colored nodes on a white background.
* `jam_igraph()` now tries harder to respect a pre-defined xlim and ylim,
before using the range of layout coordinates
* `mem_gene_path_heatmap()` slightly adjusted the number of pathway
clusters to use based upon the number of columns in the gene-pathway
matrix, slightly increasing the number for low number of columns.
* `mem_plot_folio()` new argument pathway_column_split allows setting
a specific number of pathway clusters.
* `mem_plot_folio()` argument `do_which` allows creating only
the requested plots from the full sequence of plots. This function
is likely to become the core part of the analysis workflow:
 
> gene-pathway heatmap -> pathway clusters -> Cnet using pathway clusters

# multienrichjam 0.0.22.900

## New functions

* `shape.jampie.plot()` provides a new `igraph` vertex shape
`"jampie"`. It is a complete clone to the default `igraph`
shape `"pie"` except that it offers vectorized plotting, which
can be substantially faster for large `igraph` objects that
use `"pie"` nodes.
* `jam_plot_igraph()` is a similar clone to `igraph::plot.igraph()`
except that it uses vectorized plotting specifically when there
are multiple shapes in the same `igraph` object. Formerly, each
node is drawn individually which is substantially slower, especially
for `igraph` objects with more than 100 nodes. It also by default
converts shape `"pie"` to `"jampie"` in order to use vectorized
plotting via `shape.jampie.plot()` above.

## Changes to existing functions

* `jam_igraph()` new argument `plot_function` allows use of
custom plot function, specifically to use `jam_plot_igraph()`.

## Other changes

* Functions for igraph vertex shapes were moved into a new .R file,
and into a new function family `"jam igraph shapes"`.

# multienrichjam 0.0.21.900

## New functions

* `grid_with_title()` draws a grid object with title and optional
subtitle, allowing for grid objects such as `ComplexHeatmap::Heatmap`
objects, or even multiple heatmaps. That said, any grid `"gTree"`
object should work.

# multienrichjam 0.0.20.900

## New functions

* `mem_plot_folio()` is a new function that creates the current
recommended folio of multienrichment plots:

1. Enrichment P-value heatmap
2. Gene-pathway heatmap -- importantly with pathway clustering
3. Cnet using the pathway clusters, collapsed by cluster
4. Cnet using the pathway clusters, exemplar pathways per cluster
5. Cnet using each individual pathway cluster

## Changes to existing functions

* `multiEnrichMap()` now includes `p_cutoff` in the output,
which represents the enrichment P-value threshold used for
the analysis. This cutoff is useful in making other color
gradients respect the same threshold required for significant
enrichment results, so P-values that do not meet this threshold
can be colored white (or the background color.)
* `mem_plot_folio()` new argument (minor release) `do_which` to
help produce selected plot pages from a folio of plots.

# multienrichjam 0.0.19.900

## New functions

* `collapse_mem_clusters()` is an intriguing new method of simplifying
the numerous pathways, by using each pathway cluster from
`mem_gene_pathway_heatmap()`. Each cluster is condensed to one
result, combining all genes in each cluster. The results are
surprisingly insightful, especially when numerous pathways
are present per cluster. This function also calls
`rank_mem_clusters()`, so it's possible to pick a handful of the
top pathways per cluster, as relevant.

# multienrichjam 0.0.18.900

## New functions

* `rank_mem_clusters()` is a convenience function which ranks
pathways/sets within clusters, useful with a list of clusters
following `mem_gene_pathway_heatmap()`. It makes it easier to
rank pathways within a cluster, potentially choosing one exemplar
pathway to represent each cluster. This function is part of more
effort to streamline the overall analysis workflow.

## Bug fixes

* Fixed small issue with `mem_gene_path_heatmap()` which was
properly applying `min_gene_ct` and `min_path_ct` however it
did not check the resulting data to remove empty columns and
rows. This change has been made.

# multienrichjam 0.0.17.900

## New functions

* `color_edges_by_nodes()` takes an igraph, determines a
color for each node based upon the shape (using `avg_colors_from_list()`
for shape `"pie"` and `"coloredrectangle"`), then creates an
average color between two nodes, and uses that as the edge color.

## Changes to existing functions

* `jam_igraph()` has more configurable arguments, in the form
of named lists. The driving example is
`label_factor_l=list(nodeType=c(Gene=0.01, Set=1))`, that
will apply `label.cex*0.01` to "Gene" nodes,
and `label.cex*1.5` to "Set" nodes. Similar arguments:
`node_factor_l` to apply to node size, and
`label_dist_factor_l` to apply to label distance from node center.
* `shape.coloredrectangle.plot()` was updated to fix a small bug
in order of rendering colored rectangles, it was somehow drawing the
frame after the fill colors, which caused the frames to overlap
each other and appear transparent.


# multienrichjam 0.0.16.900

The vignette was updated to use the new functions, for
a much cleaner overall workflow.

## Changes to existing functions

* `reorderIgraphNodes()` uses new function `avg_colors_by_list()`
to create one color to represent each node, then `sort_colors()`
to sort nodes by color hue, instead of previous behavior which sorted
by the hex color string. This change affects attributes `"pie.color"`,
`"coloredrect.color"` and `"color"`. Recall that this function
`reorderIgraphNodes()` is intended to help visually organize
groups of nodes by their colors, to make it easier to tell how
many nodes are each color.
* `avg_colors_by_list()` now uses `weighted.mean()` based upon the
`"c"` channel (chroma, color saturation) in order to down-weight
the effect of grey colors, which really should have no hue. The default
weight for grey is 0.1, while the maximum value for fully saturated
colors is 100. Blending `"green"` with `"grey30"` yields slightly less
saturated green, as expected.
* `subsetCnetIgraph()` now includes convenience calls to
`removeIgraphBlanks()`, `spread_igraph_labels()`, `reorderIgraphNodes()`,
and `relayout_with_qfr()`. In almost all cases, creating a subset of
a Cnet necessitates a new layout, which therefore requires new
node ordering (sorting groups of equivalent nodes by color), and
then label orientation (the angle used to offset the node label from
the center of each node).

## New functions

* `sort_colors()` and `order_colors()` uses the `farver::encode_colour()`
function to convert colors to HCL, then sorts by hue, and returns the
original sorted vector or corresponding order, respectively. For
sorting colors inside a `data.frame`, convert colors to a factor
whose levels are `sort_colors(unique(colors))`, in order to maintain
ties during multi-column sorting.
* `apply_color_cap()` imposes a numeric range onto a color channel
for a vector or list of colors, using `jamba::noiseFloor()`. Values
outside the allowed range are forced to the range.
* `subgraph_jam()` is a custom version of `igraph::induced_subgraph()`,
that correctly handles subsetting the graph layout when the graph itself
is subsetted.
* `jam_igraph()` is a custom `igraph::plot()` with proper sizing when
`rescale=FALSE`, allowing better aspect ratios for certain network
layout coordinates.
* `with_qfr()` is a wrapper to `layout_with_qfr()` that returns
a layout specification object, for use with `igraph::add_layout_()`
in a more confusing workflow pattern than I expected.

# multienrichjam 0.0.15.900

## Changes to existing functions

* `multiEnrichMap()` several updates to improve robustness
in the overall workflow.

## new functions

* `heatmap_row_order()` and `heatmap_column_order()` provide enhanced
output from `ComplexHeatmap::row_order()`, mainly that they return
the actual rownames and colnames of data in the heatmap. Very useful
when the data used for the heatmap has been filtered internal
to the function.

# multienrichjam 0.0.14.900

## Changes to existing functions

* `list2im()` removed argument `makeUnique` because the
underlying conversion to matrix no longer requires that step,
and it was a performance hit for extremely large lists.
The argument `keepCounts` is the only requirement to
maintain the count of each entry per list.
* `multiEnrichMap()` now subsets the `geneIM` gene incidence
matrix to match the genes after `topEnrichN` filtering is
applied.
* `layout_with_qfr()` now by default will use edge attribute
`"weight"` as the `weights` argument when calling
`qgraph::qgraph.layout.fruchtermanreingold()`.
* `importIPAenrichment()` was enhanced for more robust import
conditions, specifically for different variations of missing
IPA enrichment results. Also, empty colnames are removed, to
help recognize the proper identifier in each scenario.

## New functions

* `enrichList2geneHitList()` takes a list of `enrichResult`
and returns the list of genes represented in each `enrichResult`.
Intended mainly for internal use by `multiEnrichMap()`.

# multienrichjam 0.0.13.900

## New plotting functions

* `mem_multienrichplot()` allows customized enrichMap-style plotting
of igraph objects. Notably, you can filter by Jaccard overlap,
or by overlap count -- the number of genes involved in the overlap.
There are otherwise a lot of overlaps that involve only one gene,
which is not the best way to build this type of network.
* `relayout_with_qfr()` is a light extension of `layout_with_qfr()`,
but adds default calls to `removeIgraphBlanks()`,
and `spread_igraph_labels()` to help things look pretty.
* `mem_legend()` draws a color legend in the corner of a figure,
using the colors defined in the `multiEnrichMap()` output.

## Other new functions

* `subset_igraph_components()` subsets an igraph based upon
connected components -- i.e. distinct subclusters.
* `memIM2cnet()` takes the pathway-gene incidence matrix and
produces a Cnet plot `igraph` object. If given the output from
`multiEnrichMap()` it will also color nodes using the gene
and enrichment incidence matrix colors.
* `avg_colors_by_list()` and `avg_angles()` are intended for very
rapid color blending.


# multienrichjam 0.0.11.900

## Bug fixed

* `subsetCnetGraph()` was fixed to handle rare cases where
pathway set name is identical to one or more genes, which
happens with IPA pathway analysis, in the "Upstream Regulators"
output.

# multienrichjam 0.0.10.900

## New functions

* `mem_gene_pathway_heatmap()` is a wrapper to `ComplexHeatmap::Heatmap()`
which takes output from `multiEnrichMap()` and produces a
heatmap styled similar to a pathway-gene incidence matrix. It includes
row and column annotations to help interpret the results.
* `mem_enrichment_heatmap()` produces a heatmap of the pathways
and enrichment P-values for each comparison.
* `spread_igraph_labels()` is a general use function for igraph
objects, it takes an igraph network and associated layout, and
arranges igraph nodes at an angel opposite the majority of edges
from each node. The result arranges labels around the outside of
the network, and typically away from other nodes. It isn't perfect,
but is visually a big step in the right direction.
* `removeIgraphSinglets()` is a simple function that removes nodes that
have no edges. I found myself doing it manually enough times I wanted
something quick and easy.

## enhancements

* `multiEnrichMap()` now includes two new edge attributes on MultiEnrichMap
igraphs: overlap_count which contains the number of shared genes
between two sets, and overlap_max_pct which contains the max percent
overlap between two sets -- based upon the # overlapped divided
by the smaller of the two sets. The overlap_count is useful as
an optional edge label, to show how many genes are involved.
* `enrichDF2enrichResult()` now forces `"setSize"` to contain integer
values, useful when pathway size is inferred from something like
gene ratio.

## Bug fixes

* `multiEnrichMap()` was not properly assigning `colnames(memIM)`,
although they were inferred from the `rownames(enrichIM)`. The bug
was because `as.character()` removes names from a character vector,
and the conversion was done to enforce proper handling by `strsplit()`
which gives unexpected results when the list may contain a factor.

# multienrichjam 0.0.9.900

## bug fixes

* `importIPAenrichment()` no longer calls commandline `grep`
to remove blank rows, for now it uses `readLines()` and
`jamba::vigrep()` to select rows with at least one character.

# multienrichjam 0.0.8.900

## bug fixes

* Fixed several calls to `strsplit()` that failed when sent factors,
which only happens when R `options("stringsAsFactors"=TRUE)` which
is the default in base R. Now all calls to `strsplit()` enforce
`as.character()` unless character type has already been enforced.
* Added numerous package prefixes to functions, to avoid importing
all package dependencies.


# multienrichjam 0.0.7.900

## bug fixes

* `topEnrichListBySource()` and `topEnrichBySource()` were updated
to handle more default colnames for P-values, now including the
default colnames for `enrichResult` objects `c("pvalue","padjust")`.
When no sort columns are found, a warning message is printed.

## changes

* Added documentation for `multiEnrichMap()` and other functions.
* Added vignette describing the workflow starting with Ingenuity
IPA enrichment results.

# multienrichjam 0.0.6.900

## new functions

* `importIPAenrichment()` imports Ingenuity IPA enrichment results,
by default splitting each enrichment table into its own `data.frame`.
It curates colnames to be consistent with downstream analyses.
* `curateIPAcolnames()` will curate the colnames of IPA data,
and ensures the values in the gene column are consistently delimited.
* `gsubs()` is similar to `base::gsub()` except that it applies
a vector of pattern-replacement operations in order. This function
may be moved to the `"jamba"` package.
* `find_colname()` is a helper function to find a colname given
a vector of expected values, which is matched directly, then
case-insensitively, then as a vector of patterns to match the
start, end, then any part of the colnames. By default the first
matching value from the first successful method is returned.
Good for matching c("P-Value", "pvalue", "Pval")
* `topEnrichBySource()` and `topEnrichListBySource()` subset the
input pathway enrichment results by taking the top `n` result,
then making sure the overall selected pathways are retained for
all enrichment tables.

# multienrichjam 0.0.5.900

## changes

* `list2imSigned()` is no longer dependent upon the input list `x`
having names, nor having unique names.

## enhancements

* `isColorBlank()` now handles list input, which is helpful when
applied to igraph objects where `"pie"` colors are accessed as a list.
* `removeIgraphBlanks()` was updated to use vectorized logic for `"pie"`
vertex attributes.

# multienrichjam 0.0.4.900

## bug fixed

* Fixed small issues with `fixSetLabels()`.

## changes

* Updated `fixSetLabels()` to handle words the should be kept
uppercase, with some defaults pre-configured, e.g. "mRNA".

# multienrichjam 0.0.3.900

## changes

* removed `mergeAllXY()` and `unnestList()` and moved them to the `jamba`
package. Added corresponding version requirement on jamba.

# multienrichjam 0.0.2.900

## changes

* `fixSetLabels()` was updated to include a wider set of canonical
pathway prefixes used in MsigDB v6.1.

## new functions

* `cnet2im()` and `cnet2df()` are helper functions used to convert
a Cnet igraph to either a incidence matrix, or a table summary useful
for extracting subsets of pathways using various network descriptors.
* `drawEllipse()` and `shape.ellipse.plot()` add an igraph vertex
shape "ellipse" whose shape is controlled by vertex.ellipse.ratio,
where `1` creates a circular node.

# multienrichjam 0.0.1.900

## new functions

* `multiEnrichMap()` the workhorse core function, implementing the full
workflow to pre-process the results for later inspection.
* `enrichDF2enrichResult()` converts a data.frame into an
`enrichResult` object as used in `clusterProfiler` and `DOSE`.
* `cnetplotJam()` custom function to create a cnet plot igraph object.
* `mergeAllXY()` which merges a list of data.frames while keeping
all rows.
* `unnestList()` which un-nests a list of lists, resulting in a flattened
list. It directly supports `mergeAllXY()` in order to provide a simple
list of data.frame objects from potentially nested list of lists of
data.frame objects.
* `enrichList2IM()` converts a list of enrichResult objects (or data.frames)
into an indidence matrix of gene rows and pathway columns.
* `enrichList2df()` combines a list of enrichResult objects (or data.frames)
into a single data.frame, using the union of genes, and the best P-value.
* `fixSetLabels()` to update pathway/gene set labels using some
small set of logic.
* `isColorBlank()` checks if a color is a blank color, either by comparison
to known colors, or transparency, or saturation/brightness.
* `igraph2pieGraph()` converts an igraph into one with pie node shapes,
optionally coloredrectangle node shapes.
* `layout_with_qfrf()` is an extension to `layout_with_qfr()` that returns
a function, therefore convenient to use in `igraph::plot()` with custom
arguments.
* `removeIgraphBlanks()` removes blank colors in either pie or
coloredrectangle node shapes, thereby making non-blank colors easier to
see in plots.
* `subsetCnetIgraph()` takes an igraph with Set and Gene nodeType, and
subsets based upon a fixed list of Set or Gene nodes, removing
singlet disconnected nodes from the output.

## changes

* added `dplyr` to package dependencies.
* added `DOSE` to package dependencies.
