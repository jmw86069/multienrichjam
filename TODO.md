# TODO

## 03dec2025

* `mem2emap()`

	* apply shadowText, bold labels
	* filter words to remove more "stop" words

* FIXED. Fix bug with `CnetCluster()`, `CnetCollapsed()` when drawing
`mem_legend()`, it should include `do_directional=TRUE` when relevant,
and proper fill colors.
* DONE. Consider storing "hasDirection" in MemPlotFolio
* `mem2cnet()` non-directional Cnet plot would not apply a red border.

* Continue adding more abbreviations.

## 01dec2025

* Consider adding `subsetIgraph()` to handle layout and graph attributes
* Consider storing previous graph layouts.
* DONE. Consider `set_igraph_layout()` to use current (evolving) igraph layout
storage recommendations.

## 24nov2025

* FIXED. Bug: Using `mark.groups` seems to break edge bundling? Some not bundled,
and some are bundled with odd spline control points.
* DONE. Add some 'words.txt'

   * "Genes defining", "Genes encoding", "Genes involved in",
   "Genes encoding proteins involved in", "Genes Important For ",
   "Genes Regulated By ",
   * "In response to", "[Geneid=3458]."
   * "Components of ", "Which is ",
* Rename conversions, 'enrichResult_from_*' similar to igraph style.

   * `enrichResult_from_dataframe()`
   * `cnet_to_dataframe()`
   * `im_from_cnet()`
   * `enrichList_to_im()`, `enrichList_to_dataframe()`

* `enrichDF2enrichResult()` - consider keeping pvalueColname intact then
adding 'p.value'.
* Decide where to share im-to-list, list-to-im conversion functions.
* Enhance `list2concordance()`

   * Add argument with 'method', then add methods.
   * Add IPA z-score: `(N_agree - N_disagree) / sqrt(N)`
   * Review previous, and find its source:  
   `dscore <- (N_agree - N_disagree) / N`  
   `log2(1.2 + abs(dscore)) * sign(dscore)`

* `plot_cnet_heatmaps()`, `CnetHeatmaps()`

   * Consider "dev" mode that prints some unit grid onto the plot, to help
   align labels and adjust margins.
   * requires expression data, column groups

* `MemPlotFolio` enhancements

   * Make it so `GenePathHeatmap()` and `EnrichmentHeatmap()` can dynamically
   re-create the heatmap, so customizations can be made without re-running
   `prepare_folio()`?

* Revisit exporters like `Mem_to_dataframe()`, for supplemental tables.
* `rbind()`, `cbind()`, `merge()` to combine multiple Mem objects.
* Evaluate gene clusters, consider pathway/gene cluster Cnet plot,
with some way to label the group of genes in each cluster, e.g. grid table box.

## 18nov2025

* Usability / bugs

   * `mem_plot_folio()`
   
      * FIXED. Bug: `list` for `pathway_column_split` with subset pathways;
      it seems to have clusters as `character` with this approach.
      * FIXED. applies `simple_anno_size` to row and not column in
      gene-path heatmap.
      * DONE. Make 'mpf' argument also use `Clusters()` and `GeneClusters()`.
      May need to review whether this behavior is optimal for genes especially.

   * FIXED. `fixSetLabels()` add option for max width, truncating anything longer
   instead of applying word-wrap. (Already available, haha `maxNchar`.)
   * DONE. `mem2cnet()` consider aesthetic defaults:
   `vertex.label.font=2`, `use_shadowText=TRUE` (store in attributes?)
   * DONE. Change default `mem_gene_path_heatmap()` simple_anno_size to 5mm.
   * UNSURE. Consider new default cluster_column_slices=TRUE, and for rows.
   It extends the dendrograms with connections across slices.

* Consider new function `merge_Mem_sets()`

   * intended to combine "equivalent" pathways
   (TLR signaling v1, TLR signaling v2) together.
   * Bonus points: Pathways with identical names (v1, v2, etc)
   with greater than threshold (75%?) identical genes,
   automatically merge them.

* `CnetCollapsed()` need some way to adjust layout, and apply reorder,
spread labels, etc.
* `CnetExemplar()` consider adding pathway title as prefix,
for example "A: Pleural Mesothelioma", to indicate the cluster for each
exemplar.
* `removeIgraphSinglets()` should subset graph_attr("mark.groups").
* `jam_igraph()`

   * with `mark.groups=TRUE` and non-Cnet, it should try `components()`,
   however using `mark.groups=igraph::components` works currently.

* `make_alpha_hull()`

   * try to split extremely disparate node groups,
   e.g. two sides of the plot, should ideally split into two subcomponents.
   * draw hull labels slightly farther from the edges? It might need
   something helpful from `jamba::drawLabels()`.

* Usability / Workflows

   * Consider workflow to re-Mem a subset of pathway clusters.
   Situation: Pathway clusters often have one "miscellaneous" cluster
   with few genes, no coherence.
   Target: Re-create Mem without this cluster.
   `prepare_folio(Mem, pathway_column_split=Clusters(Mpf)[c(1, 3, 4)])`

* Ideal world: `GenePathHeatmap(Mpf[, clusters])` would use the subset.
* Debug `rank_mem_clusters()` for follow-up with using Mem/Mpf input.

## 14nov2025

* Add vignette for `launch_shinycat()`
* Add `cbind()` and `rbind()` to combine multiple Mem objects.

   * Use case: enrichment with IPA; enrichment with clusterProfiler;
   then `rbind()` them.
   * Use case: enrichment of A,B in IPA; enrichment of C,D in IPA;
   then `cbind()` them.

## 12nov2025

* COMPLETE for now. Continue polishing the workflow:

   * `mem_plot_folio()`
   
      * Accept optional `MemPlotFolio` input, to re-use identical settings.
      * consider new `PrepareFolio()` with `do_plot=FALSE`?  
      Workflow:
      `multiEnrichMap()`
      `PrepareFolio()`
      individual plots: `EnrichmentHeatmap()`, `GenePathHeatmap()`, etc.
      `PlotFolio()`
      * new arg `do_whichbyname` accept `character` string?
      * consider replacing plots 3,4,5 with one Cnet, then move the custom
      label into `CnetCollapsed()`
      
   * `CnetHeatmap()`, `GenePathHeatmap()` could use `Mem` input?
   Secretly call `mem_plot_folio()` behind the scenes.

* COMPLETE. Move `fixSetLabels()` logic to config file which is prepared upon
loading the package.

   * COMPLETE. words_from, words_to; add_from, add_to; abbrev_from, abbrev_to:
   Should they be `data.frame` objects?

## 10nov2025

* COMPLETE. `IPAlist_to_hits()` to extract hit list from IPA xref tables.

   * convert "/" to ":" for genes with two or more symbols
   * remove `" (includes others)"` suffix for proper gene symbol matching

* COMPLETE. `mem_plot_folio()`

   * Replace CnetCollapsed with one `igraph` that updates `V(cnet)$label`
   as needed. The rest should be constant.
   
* `colorRamp2D()` enhancements

   * Use `thresholds(Mem)$direction_cutoff`, add labeling
   * consider continuous color option, also use `direction_cutoff`

* Consider making many functions "internal" to simplify the user experience.

### New ideas

* COMPLETE. Consider some way to score gene clusters by how well they "define" the
cluster? For example, some genes are found in the majority of nodes in
at least one cluster. Other genes are scattered across clusters with
no visible aggregation. The metric could be fraction of genes in at
least half the pathways in one pathway cluster.

## 04nov2025

* PARTIAL. Simplify the function reference categories.

   * DONE. Move core user-facing functions to top categories.
   * Move almost everything into "Details".

* Split vignettes into concise subtopics.
* COMPLETE. Consider `Biobase::Versioned` to version `Mem` class structure.
* COMPLETE. Write `updateObject()`.
* Make some Mem slots optional, consider gradual reduction of slots?

   * Use slot type "matrix_OR_NULL" with
   `setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))`.
   * The `geneIMcolors()` accessor would return geneIMcolors when defined,
   or generate the color matrix dynamically using thresholds and colorV.
   * Candidates:
   geneIMcolors, enrichIMcolors, geneIMdirection
   multiEnrichDF, multiEnrichList

* Add `cbind()` and `rbind()` to combine multiple Mem objects.
* Consider providing some data dynamically, using current thresholds:
`enrichIMcolors`, `geneIMcolors`
* Add `plot_cnet_heatmaps()`. It needs pathway and expression data to test.
* COMPLETE. Add `MemPlotFolio` S4 object

   * Store data in `plotdata` as a `list` or `SimpleList`
   * Store parameters in a slot `parameters` as `list` as well.
   * Accessors: GPHeatmap, EHeatmap, Cnet
   * In theory, all data could be generated via accessors,
   and just store the parameters.

* Consider serialization functions

   * save series of incident matrices, data.frame objects
   * confirm it can be re-loaded directly to `Mem` object

## 29oct2025

* Mem S4 object updates:

   * Consider better mechanism to rename pathway gene sets with `sets()<-`,
   specifically so that the `enrichList` results are also kept consistent.
   Perhaps an alias to change from original name to custom name?
   * Consider auto-updating `Mem` when parameters/thresholds change.
   * Consider removing some slots and replacing the accessor with function:
   
      * `geneIMcolors` and `enrichIMcolors` could be done using
      with dynamic method which applies `colorV`,`p_cutoff`,`min_count` to
      produce the color matrix. Would enable changing thresholds and having
      up-to-date data. Also would remove the slot and reduce S4 complexity.
      * `multiEnrichDF`,`multiEnrichResult`

* Consider S4 classes:

   * `MemPlotFolio`
   * `multiEnrichResult` - simple extension to `enrichResult`?
   * `enrichResultList` - need to determine if it enables something useful,
   other than just making it clear the type of data

## 23oct2025

* Organize vignettes using diataxis philosophy

   * Tutorial: Use IPA; Use clusterProfiler
   * How-to:
   
      * Define pathway-gene clustering parameters
      * Customize Mem Enrichment Heatmap
      * Adjust Cnet plot layout

   * Explanation:
   
      * Determine appropriate parameters for `multiEnrichMap()`
      * Decide Cnet strategy for a visual summary

## 17oct2025

* DONE. Remove all dependencies on arules::transactions
* Clean up the IPA vignette figure sizes, move rambling text to the end.
* Add vignette for `clusterProfiler` and `msigdbr`.
* IN PROGRESS. Convert to Mem S4 object

   * Update all functions to accept 'mem' or 'Mem' input.
   * Update `multiEnrichMap()` to return 'Mem'.
   Consider removing all sections which create cnet/emap objects.

* Add emap to `mem_plot_folio()`, or at least "make it easy".
* Add Cnet-Heatmap and supporting functions.
* Add testthis unit tests for key behaviors.
* ShinyCat

   * Fix bug: singlet node is highlighted, nodesets cannot be highlighted.
   * Consider `shinydashboard::tabBox()` for inputs; or `bsplus::bs_accordion()`
   * Consider adding global layout options
   
      * Relayout with new repulse threshold. Should "nullify" other adjustments.
      * Nodeset spacing overall. (Rather than having to click each nodeset)
   
   * DONE. Consider highlighting edges around the highlighted node/nodeset
   * DONE. Add option to save igraph object to a local file
   * DONE (removed). Add (or remove) data table tab
   * DONE. Document how to retain data

* igraph shapes

   * clean up use of pie.border, pie.lwd, pie.width,
   frame.color, frame.width, frame.lwd
   Official igraph now uses "width" instead of "lwd" so probably best.
   * Consider JamPolygon? Would also require `grid`, which would require
   porting all `igraph` shapes to `grid` equivalents.
   Competes with long-term goal to integrate with ggplot2/tidygraph/ggraph.

* tidygraph/ggraph long-term

   * ggraph requires porting 'jampie' (and potentially 'coloredrectangle')
   vertex shapes into ggplot2/ggraph.
   * ggraph may require porting edge bundling to ggraph.
   * ggraph may require porting 'qfr' layout to include argument 'repulse'.
   * tidygraph workflows may warrant creating examples which use tidy-oriented
   logic instead of base R (igraph-like) syntax.

* reduce R dependencies

   * DONE. alphahull - requires sp, ggplot2, spatstat. Port `alphahull::ashape()`

## 09oct2025

* New function: layout with fixed nodes

   * Purpose is to "make it easy" to apply relayout while holding
   a subset of nodes fixed in place. Simple wrapper around existing
   capability, similar to `adjust_cnet_set_relayout_gene()`.
   Potential to relayout one nodeset, keeping all others static.

* DONE. Support "script-like" Cnet adjustments.

   * Simple text line: node/nodeset:x:y:spacing:rotate

* Performance:

   * Response drawing a plot is noticeably slow in R-shiny context,
   some basic strategies of speeding the process should be applied.
   * Consider storing "edge bundling group" in `igraph` edge attributes,
   since an edge can only belong to one "bundle".
   For Cnet, the default is the Cnet nodeset, however any set
   of edges can be bundled.
   * Consider pre-calculating edge bunding splines, then caching in
   `igraph` attributes for re-use.
   It must be invalidated by: new layout; any change to layout coordinates.
   Rotating the igraph layout could also rotate the bundled splines.
   * Consider pre-calculating mark group alpha hull polygons.
   
* DONE. R-shiny Cnet Adjustment Tool 'ShinyCat' prototype

   * Use `jam_graph()`, associate mouse click x,y to specific node

* DONE. `summarize_igraph_spacing()`

   * Defaults: scaled, within-nodeset, each group
   * call `get_cnet_nodeset()` when not supplied

* Improve igraph node label placement:

   * Consider adding `label.adjx`, `label.adjy` to nudge a label,
   to be applied after `label.dist` and `label.degrees`.
   Unclear if units should be "relative" (by layout) or "absolute",
   "relative" seems more convenient as a default.
   * Reducing label overlaps is a time sink and should be easier,
   or automated, or both.
   Related to whether ggraph/tidygraph rendering will ever become viable.
   If yes then use ggrepel.
   Using ggraph/tidygraph requires porting `jam_igraph()` features:
   edge bundling, node shapes (jampie, coloredrectangle).
   If no, then decide whether to use grid or base R graphics to estimate
   label sizes, then implement non-overlapping labels (somehow).
   Ultimately, that process would use `label.adjx`, `label.adjy`.

## 27mar2025

* Consider approach to store/cache edge bundling for a given layout.
* Consider longer-term use of tidygraph/ggraph environment, which
requires creating pie node type geom, and foregoing edge bundling
until eventually that logic can be ported into that system.

## 13mar2025

* S4 `Mem` object should probably also include Mpf output in optional slots.
* DONE. Fix: `mem_plot_folio()` ignores `byCols` when sorting exemplar pathways.

## 12feb2025

* `importIPAenrichment()`

   * Update to match new variations in the "Export All" text output from IPA.

## 15dec2024

* `mem_enrichment_heatmap()`, and `mem_plot_folio()`

   * Consider option to abbreviate pathway names on the fly.
   Column names are sometimes too ridiculously long,
   and `column_names_max_height` isn't effective enough.
   * Add the clustering method info to legend as in `mem_gene_path_heatmap()`.

* Functions using `nodeset` names

   * Protect against sorting errors, e.g. `nodeset="D,A"` should
   recognize nodeset whose sorted label is: `"A,D"`
   * See `adjust_cnet_nodeset()` and argument `set_nodes`

## 04dec2024

* Fix bug "Error in gp_hm@column_dend_param$fun(t(gp_hm@matrix)) :
  attempt to apply non-function"

   * Seems caused by small number of pathways in the memIM/enrichIM matrix.

* Consider removing `sf`.

   * Potentially replace with direct calls to `polyclip`.
   * Replace with `JamPolygon`, which would suggest replacing with `grid`.
   This would be a "hobby option" - would require updating JamPolygon.

## 12nov2024

* Ideas for S4 helper functions / accessors

   * `exemplars()` - extract exemplar pathways for Gene-Pathway clusters
   * `as_data_frame()` - convert to formal data.frame format
   
      * `Mem`
      * `Mpf` - one section to describe each pathway cluster

## 02nov2024

* S4 objects to consider creating:

   * `Mem` - MultiEnrichMent
   
      * Consider storing MPF clustering params here
      * Consider storing pathway clusters here, to allow custom assignment
      of pathways to clusters in mpf
   
   * `MPF` - Mem Plot Folio
   
      * Consider storing `mem_plot_folio()` output as fixed object.
      Alternative is to add into `Mem` above.
   
   * `Cnet` - Can it inherit from `igraph`?
   (I'm really unsure that it would be more benefit than hassle.)

      * Methodology
      
         * Apparently yes. See `setOldClass()`;
         virtual slot `.Data` for S3 object;
         `as(x, "S3)` to convert to S3 form;
         `extends(class(x))` as equivalent of S3 inheritance;
         see https://stackoverflow.com/questions/64683532/r-as-method-for-converting-custom-s4-class
      
      * Benefits
      
         * Enforce constraints unique to Cnet-type objects: `nodeType`,
         `cnet_nodesets`, etc.
         * Direct the `plot()` function to `jam_igraph()`.
         * Direct other "cnet" functions to use this object type.

* `mem_plot_folio()`

   * Consider adding Mem params to the heatmap legend, e.g.
   
      * `Top Pathways N (topEnrichN)`
      * `Gene min count`
      * `Number of pathway clusters`
      
   * Add Multi-Enrichment Plot to the workflow, `mem2emap()`.

* Add `mem_cnet_heatmap()` - "Cnet-Heatmap"

   * Cnet network plot (center) surrounded by heatmaps for each Cnet cluster.
   

* Consider compatibility with other Bioconductor S4 object classes:

   * `clusterProfiler::enrichResult` - used to some extent already here
   * `clusterProfiler::gseaResult` - not used directly here
   Ref: https://doi.org/10.1016/j.xinn.2021.100141
   * `EGSEA::EGSEAResults`, see `topSets()`, `plotHeatmap()`, `plotPathway()`,
   `plotMethods()`, `plotSummary()`, `plotGOGraph()`, `plotBars()`,
   `limmaTopTable()`

## 01nov2024

* Minor: Consider using `cli` for messaging, for consistency with tidyverse.


## 29oct2024

* General goal is to reduce label overlaps in Cnet plots

   * Many options, including considering `ggraph` ot `tidygraph`, however
   both require re-creating pie/jampie graph node shapes that display both
   fill colors and border outline (without hiding adjacent borders).
   These strategies would not respect edge bundling, so that would either
   be lost or would require porting to ggplot ecosystem somehow. Sigh.
   Then use `ggrepel` for non-overlapping labels.
   * Ideal world: There exists some `grid` tool for non-overlapping labels,
   and not `ggplot2`.
   * A simpler option is to allow direct x/y coordinate label adjustment.
   Current approach requires angle and distance from node center. Difficult
   to move one label "up and left a tiny amount". If given ability to
   adjust a label, it could be possible to adjust all labels to reduce
   overlaps.
   Potential workflow:
   
      * Define label coordinates after calculating angle and distance
      from node center per usual.
      * Apply adjustments.
      * Apply user-provided adjustments. Easiest approach.
      Default nudge matrix is c(0, 0) for all nodes.
      Edit for specific nodes, relative to plot layout coordinates.
      
      * Future: Automated adjustments.
      
         * Use `strwidth()` and `strheight()` to define bounding boxes.
         * Adjust bounding boxes?
         * The `plotrix` approach performs radial search in local space.
         It works best when labels are placed in smart order, e.g. center
         in a cluster, then working out toward edges.
         Define clusters? Define center of each cluster. Sort labels by
         distance from cluster center.
         * Force-directed: calc overlap of two labels, center of overlap,
         center of bounding box, calculate angle, relative "strength",
         and repel along the opposite angle.
      


## 09oct2024

* `fixSetLabels()`

   * Define better default replacements.
   * Consider using default from/to, with option for user-defined additions.

## 09sep2024

* Consider using `mem_enrichment_heatmap()` as the top_annotation
for `mem_gene_path_heatmap()`, including for example the dotplot format.

   * Easiest prototype might be to use `+` to display the
   gene-path heatmap rotated with `rotate_heatmap=TRUE` beside the
   enrichment heatmap,then merge the color legends:
   ```R
   ComplexHeatmap::draw(mpf$enrichment_hm + mpf$gp_hm,
      annotation_legend_list=attr(mpf$enrichment_hm, "annotation_legend_list"),
      merge_legends=TRUE)
   ```

## 04sep2024

* DONE. Fix the `row_title` for the enrichment heatmap when called by
`mem_plot_folio()`, currently it uses numbers instead of `LETTERS`.
* DONE. Consider option for fixed-attribute cells for `mem_enrichment_heatmap()`
when used with dot plot, so each cell is square with the circle centered.
* Consider returning updated `igraph` from `jam_igraph()`

   * Currently returns `invisible(NULL)`.
   * For example, after applying all the updates to `igraph` values
   such as node size, label size, label distance, etc.
   * The end goal is to be able to call `jam_igraph(output)` with no
   additional arguments and have it (mostly) render the identical figure.

## 30aug2024

* Add unit tests.

   * Cover all basic functionality.
   * Examples and vignettes should already cover the core workflow.

* Consider adjustments to `make_point_hull()`

   * font size
   * distance from hull

* Big picture: Consider using `venndir::JamPolygon` instead of `sf` polygons.

   * Would calculate offsets, rescale, transform, using `JamPolygon` functions.
   * Biggest benefit is to reduce dependencies, `sf` is heavy.


## 16aug2024

* `mem_gene_path_heatmap()`

   * Changed default caption to use `ComplexHeatmap::Legend()` format
   so it can be included along with color legends.
   * Increased caption font, it was not legible.
   * Moved the `gene rows, set columns` labels to the top, split in two rows.
   * Return added attribute `"caption_legendlist"`.

* `apply_cnet_direction()`

   * Change default logic to use `pie.border` always instead of
   switching to `frame.color` when all borders are identical.
   I think this logic should belong in the plotting function
   to decide how to render nodes, `jam_igraph()`
   Also more convenient when reviewing the data, user should only
   have to look at `pie.border` and not combination of `pie.border`
   and `frame.color`.


* `mem_enrichment_heatmap()`

   * Consider adding `top_annotation` using the `colorV` colors per
   enrichment, for consistency with `mem_gene_path_heatmap()` and Cnet plots.

* DONE. Add `add_pathway_direction()`

   * DONE. Helper function to add directional z-score column to `enrichResult`
   data, using the formula from QIAGEN IPA:
   `z <- (N_genes_up - N_genes_down) / sqrt(N_genes_up + N_genes_down)`

* Call `add_pathway_direction()` from `multiEnrichMap()` when appropriate.

   * `geneHitList` or `geneHitIM` are provided, and
   * contain positive and negative values, and
   * `"z-score"` is not already defined in each `enrichResult` object.
   (Bonus points for checking the direction colname attribute.)
   
* Review `clusterProfiler::compareClusterResult-class` object
definition. For example `clusterProfiler::merge_result(list(enrichResults))`.

   * Consider some form of integration, if possible, for example conversion
   to/from to call similar functions in `clusterProfiler`.

* S4 object `Mem` to replace `list` returned by `multiEnrichMap()`

   * IT IS COMING
   * Consider `multienrichjam()` to replace `multiEnrichmap()`?
   * **slots**:
   
      * geneIM (im, direction, colors): `matrix` objects
      * enrichIM, (pvalues, direction, geneCount, colors): `matrix` objects
      * memIM: `matrix` object
      * enrichList: `list` of `enrichResult` objects
      * colorV: `character` vector of colors per enrichment
      * colnames: `character` column assignment (consider omitting to enforce
      standardized colnames)

   * **methods**:
   
      * `mem_plot_folio()`, supporting functions: `mem_gene_path_heatmap()`,
      `mem_enrichment_heatmap()`, `mem2cnet()`
      * `enrichList()` - accessor for `mem@enrichList`
      * `mem2dfs()` - create series of `data.frame` summarizing content,
      intended for export to Excel xlsx.
      * `mem2xlsx()` - direct export to Excel xlsx, calling `mem2dfs()`.

   * **behaviors**

      * `multiEnrichMap()` returns `Mem` object instead of `list`
      * `mem_plot_folio()` may store parameters in the `Mem` object
      * `mem_gene_path_heatmap()`, `mem_enrichment_heatmap()` could also
      store/retrieve parameters from the `mem` input object.
      * `mem_plot_folio()` optionally stores plots into `mem`
      to maintian consistent plot attributes

* Consider `Cnet` object that inherits `igraph`?

   * It behaves as an `igraph` object except its class is helpful for
   generic functions: `plot.Cnet()`, `layout.Cnet()`, `relayout.Cnet()`
   * Unclear if S3 object type is preferred since it could inherit
   `igraph` (S3 object) characteristics.

* Consider `subgraph.Cnet()` or `subgraph.igraph()` functions

   * Main purpose is to subset the layout as well as nodes.


## 07aug2024

* Remove all `require()` checks, since they should already be in Dependencies.

## 19jul2024

* Remove `gsubs()` which causes a warning upon loading `multienrichjam`.
It conflicts with `jamba::gsubs()`. They have slightly different logic.
* `mem_plot_folio()`

   * Enrichment heatmap should define `row_title` to match
   `pathway_column_title`, `LETTERS` by default.
   It currently shows numbers.

## 11jul2024

* `mem_plot_folio()` - option to support RMarkdown output

   * Provide optional wrapper for RMarkdown output, specifically to
   print headings/tabs for each plot produced.
   * Slight downside is there isn't an easy way to configure a unique
   figure size for each plot, so plot sizes are at the mercy of the
   Rmd chunk options `fig.height`, `fig.width`.
   * One implementation option is to allow `hook_preplot` and `hook_postplot`
   to allow the user to run a custom function before and after each plot
   is drawn.
   That feels too complicated, when the main driver is just to print
   an RMarkdown header.
   * Investigate whether each Markdown tab can define a new figure size.

## 21jun2024

* `importIPAenrichment()`

   * DONE. Consider handling gene identifiers so that the default behavior
   refers to each enriched gene by the original input gene symbol,
   and not the IPA-generated gene symbol.
   Some considerations:
   
      * The main driver is to associate pathway genes to the input data,
      using the same identifier as the input data. Sometimes IPA
      assigns its own name which will not match the input data.
      * We may consider that data integration (comparison across enrichments)
      may perform better by comparing via IPA Symbol. Consider the example
      `"HSPA1A/HSPA1B"`, one experiment may find `"HSPA1A"` significant,
      another may find `"HSPA1B"` significant. According to IPA, these
      results are equivalent, in which case the IPA symbol `"HSPA1A/HSPA1B"`
      would allow them to use the same identifier.
      * Sometimes the "user input" is a platform ID, such
      as Affymetrix probe. In this case may not be preferable to
      use the "user input". In this case it may be convenient for data
      integration, but less convenient when trying to recognize
      gene symbols as labels.
      * This step requires `"Analysis Ready Molecules"` is available,
      and follows expected convention used by Ingenuity.
      * When multiple genes are combined to "one entity"
      by IPA, and only one input symbol is retained in the
      `"Analysis Ready Molecules"`.
      The driving use case: `"HSPA1A"` and `"HSPA1B"` are combined to
      one gene entry `"HSPA1A/HSPA1B"` by IPA. They are considered one
      gene for the purpose of enrichment. If one or both genes were significant,
      they would appear as `"HSPA1A/HSPA1B"` by IPA.
      The `"Analysis Ready Molecules"` will list one symbol `"HSPA1B"`
      as exemplar, and no entry will appear for `"HSPA1A"`.
      In this case there are two options:
      
         1. Use only the IPA assigned input symbol as provided: `"HSPA1B"`.
         IMPLEMENTED with `revert_ipa_xref=TRUE` (new default).
         2. Split the IPA multi-symbol label into component parts `"HSPA1A"`
         and `"HSPA1B"`. In this case, we need the user to supply the
         actual gene hits, so we retain only the gene hit the user provided.
         NOT IMPLEMENTED.
         3. Leave the entry as `"HSPA1A/HSPA1B"`, however this symbol will
         not match any input gene hit list, or other expression data matrix.
         IMPLEMENTED with `revert_ipa_xref=FALSE`.

   * Consider retaining the header section with analysis details,
   at least for text input.

* `multiEnrichMap()`

   * Sometimes the gene rows in `geneIM` do not match gene rows in `memIM`,
   causing an error downstream. The problem appears to happen when the
   gene hit list does not match all entries in `memIM`, causing `geneIM`
   to have fewer rows.
   * One solution is to reduce `memIM` - not ideal because we do not
   want to lose data.
   * Another option is to expand `geneIM` - which requires inferring the
   incidence matrix. For directional data, it would impose `1` regardless
   of the intended directionality, since no other source is available.



## 16may2024

* edge bundling

   * Between two communities currently calculates the "central path" between
   the group centroids. Consider calculating the "central path" between
   the points with edges between the two communities.

## 10may2024

* `igraph` adjustment scripting language?

   * Combines:
   
      * `nudge_igraph_nodes()`: Fuscobac::x:0.01:y:-0.02
      * `adjust_cnet_nodeset()`: nodeset:A::degrees:45:x0.1:y-0.05:percent_spacing:7
      * `apply_nodeset_spacing()`: nodeset:B::percent_spacing:7

* `adjust_cnet_nodeset()`, `apply_nodeset_spacing()`

   * Debug issue where `label.dist` and `label.degrees` are defined for
   subset of nodes, leaving other node attributes `NA` - which causes an
   error in `jam_igraph()`. Occurs when only one nodeset is adjusted.
   * Debug issue where `percent_spacing` fails in `adjust_cnet_nodeset()`
   when supplying custom nodegroups from community detection, and are
   not proper Cnet nodesets.


## 03may2024

* `mem_gene_path_heatmap()` - lower priority but could be useful

   * Consider interactive view (plotly? InteractiveComplexHeatmap?)
   to enable hover text with the enrichment P-value, gene count, z-score.
   Could be particularly useful with directional output.
   * Consider optionally labeling the significant dots?
   
      * P-value
      * z-score
      * number of genes.

* `mem_legend()` - lower priority, eventually necessary

   * Consider using `ComplexHeatmap::Legend()` for consistency with other
   legends, and to allow combining multiple legends together.

* S4 object `mem` (or `MEM`?) - higher priority - sooner the better

   * streamlined data content:
   
      * geneIM (im, direction, colors): `matrix` objects
      * enrichIM, (pvalues, direction, geneCount, colors): `matrix` objects
      * memIM: `matrix` object
      * enrichList: `list` of `enrichResult` objects
      * colorV: `character` vector of colors per enrichment
      * colnames: `character` column assignment (consider omitting to enforce
      standardized colnames)
      
      * omit: `multiEnrichMap` - in favor of `mem_plot_folio()`, `memIM2cnet()`
      * omit: `multiCnetPlot` - in favor of `mem_multienrichplot()`
      
      * optional: store output from `mem_plot_folio()` to keep a series
      of plots coordinated, using the same options: `pathway_column_split`,
      `gene_row_split`, `enrich_im_weight`, `gene_im_weight`,
      `column_method`, `row_method` (rename `column_method` to `pathway_method`?)
      * optional `enrichment_hm`: `Heatmap` output from
      `mem_enrichment_heatmap()`?
   
   * methods:
   
      * `mem_plot_folio()`, supporting functions: `mem_gene_path_heatmap()`,
      `mem_enrichment_heatmap()`, `mem2cnet()`
      * `enrichList()` - accessor for `mem@enrichList`
      * `mem2dfs()` - create series of `data.frame` summarizing content,
      intended for export to Excel xlsx.
      * `mem2xlsx()` - direct export to Excel xlsx, calling `mem2dfs()`.

   * behaviors

      * `multiEnrichMap()` creates `MEM` object instead of `list` by default
      * `mem_gene_path_heatmap()`, `mem_enrichment_heatmap()` could also
      store/retrieve parameters from the `mem` input object.
      * `mem_plot_folio()` optionally stores plots into `mem`
      to maintian consistent plot attributes


* `mem_plot_folio()`

   * argument `do_which` - consider accepting `character` string terms,
   e.g. `"enrichment_hm"`, `"gp_hm"`, `"cnet_collapse_set"`
   * Consider new argument `clusters_mem` for these uses:
   
      * allow user-defined pathway clusters
      * allow user-defined pathway subsets (missing pathways are dropped)

* `collapse_mem_clusters()`

   * When provided `mpf$clusters_mem` as `list` it may result in singlet
   genes not connected to any pathways - these should (by default) be removed.

* `jam_igraph()` - debug

   * Apparently sometimes with singlet gene nodes it produces an error:
   `"Error in FUN(X[[i]], ...) : !anyNA(x) is not TRUE"`
   traceback pointed to this line:
   `sf::st_polygon(list(polym)) at jamenrich-igraphshapes.R#1635`


## 12mar2024

* `mem_gene_path_heatmap()`

   * Consider option to place the caption elsewhere.
   
      * Problem: The caption sometimes covers part of the color legend
      in the bottom-right corner.
      * Another workaround might be to customize the legend layout,
      so it is not blocked by the caption.

* `mem_plot_folio()`, `mem_gene_path_heatmap()`

   * Consider a workflow to merge pathway clusters, to allow flexibility
   in how pathway clusters are defined.
   
      * Problem: Pathway clusters are sometimes defined inconsistently,
      where clusters `"A"` and `"B"` might be nearly identical.
      * Problem: It might be visually apparent how to sub-divide pathways,
      the user may need a mechanism to define pathways to specific clusters.

## 04mar2024

* `mem_gene_path_heatmap()`

   * consider adding column `top_annotation` with pathway directionality,
   equivalent to the `left_annotation` used for gene directionality.
   * Add row and column annotation padding by default, to help distinguish
   the central heatmap from the left and top annotations.

* `multiEnrichMap()`

   * when supplied with `geneHitIM` or `geneHitList`, calculate the
   `z-score` using the IPA formula:
   `z <- (N_genes_up - N_genes_down) / sqrt(N_genes_up + N_genes_down)`
   as described [https://doi.org/10.1093/bioinformatics/btt703],
   and in their FAQ: [IPA FAQ - Statistical Calculations](https://qiagen.my.salesforce-sites.com/KnowledgeBase/KnowledgeNavigatorPage?id=kA41i000000L5nQCAS&categoryName=BioX)

## 30jan2024

* DONE. `multiEnrichMap()`

   * DONE. `geneHitIM` and `geneHitList` are not behaving as intended, nor
   consistently. They should be interchangeable and equivalent.
   When `geneHitIM` is supplied, it populates `geneIMdirection` but
   not `geneIM`.
   When `geneHitList` is supplied, it populates `geneIM` but not
   `geneIMdirection`.
   When either are supplied, they should populate the relevant row
   in `geneIM` and `geneIMdirection`.

## 19jan2024

* Currently it is cumbersome to edit pathway labels.
It could be done for the `mem` object itself, however the adjustment
might need to be different for different plot outputs: heatmap may not
work well using word-wrap, while Cnet plots might work best with word-wrap.
Publication figures might need an abbreviated label to save plot space.

* `mem_plot_folio()` 

   * Consider argument to enable custom adjustment of pathway labels.

* `mem_gene_path_heatmap()`

   * DONE. Consider argument `gene_annotations` to enable `"geneIM"`,
   `"geneIMdirection"`.
   * DONE. When `mem$geneIMdirection` is present, include directionality
   in the gene clustering step.
   * REWRITTEN ABOVE. Consider option to display pathway z-score
   (`mem$pathwayIMdirection`) similar to display of `mem$geneIMdirection`.
   New argument `pathway_annotations`.

## 06nov2023

* DONE. Debug issue when rendering edge arrow heads, they look wonky.
* Remove dependency on `sf` in `adjust_polygon_border()`

   * use `polyclip` package as it is much more lightweight, without requiring
   RGEOS, LWGEOM, other world globe map-coordinate based libraries
   which are not easily compiled on all computer systems using R.

* Consider porting `deconcat_df2()` to `jamba` for wider re-use by jam
packages that would not otherwise need dependency on `multienrichjam`.
* Consider vectorizing edge arrow size

   * Currently all edge arrows must be the same size (the same limitation
   is present with `igraph::plot()`).
   * This is a niche feature, arrows are rarely used in `multienrichjam`,
   and use of differently-sized arrows does not have a driving use case.
   (This feature is unlikely to move forward until needed.)

* Consider `ggraph` compatibility in future.

   * Many R packages use `ggraph` for `igraph` plotting, so it might make
   sense for consistency to offer this feature. The output `ggplot` objects
   can be combined into larger figures using `patchwork` or `cowplot`
   easier than using base R `plot()` functions.
   * My preference not to use `ggraph` is mainly because the
   returned `ggplot` object is not a useful network graph, it is only the
   instructions for visualization. As such, the layout, node customization,
   is not persistent outside the `ggplot` object. The package
   `clusterProfiler` migrated to `ggraph` and now all the returned objects
   are not useful because they cannot (easily) be customized and analyzed.
   * However `ggraph` does not offer `pie` nor `coloredrectangle` node shapes.
   * Note that `pie` nodes can be emulated with `geom_pie()` but
   (1) is painfully slow to render because it is not vectorized,
   (2) does not use inner borders, which allow adjacent wedge borders to
   be shown beside each other without overlapping.
   * Task is to add new `geom` for `jampie` node shape that accepts
   pie wedge border colors drawn as inner border, with optional outer border
   color that also does not overlap the optional inner border.
   * Another task is to implement edge bundling using `ggraph` compatible
   methods. It involves calculating edge bundles, then rendering the
   curved edges, optionally with arrow heads.
   * Ultimately a utility function `jam_ggraph()` may be the `ggraph`
   equivalent of `jam_igraph()` for plotting `igraph` objects.

## 12oct2023

* Integrate directionality with more steps in the workflow:

   * `mem_gene_path_heatmap()`
   
      * Consider option to display direction of change in row (gene)
      annotations for each enrichment.
      Could display with or without enrichment-colored results?

* Create `bookdown` documentation

   * Should it use a separate Github repository? (Probably yes.)
   It looks like `"jokergoo/circlize_book"` is the repo for the bookdown site.

* Refactor `multiEnrichMap()` - maybe new function `multienrichjam()`?

   * create `mem` S4 object with slots, print, summary functions
   
      * suggested methods:
      
         * `plot()` could default to `mem_gene_path_heatmap()` to show all data
         * `print()` could print summary of content:
         enrichments, pathways, genes
         * consider `as.data.frame()` - convert to wide `data.frame` summary?
         * `memIM()`, `geneIM()`, `enrichIM()` convenient access to slot data
         * convenient way to get `list` format for IMs?
         
      * slots:
      
         * `memIM`: gene/pathway matrix
         * `geneIM`: gene/enrichment matrix
         * `enrichIM`: pathway/enrichment
         * `geneIMdirection`: optional direction per gene/enrichment
         * `geneIMcolors`: colors assigned per gene/enrichment
         * `enrichIMdirection`: optional direction (z-score) per pathway/enrichment
         * `enrichIMcolors`: colors assigned per pathway/enrichment
         * `enrichIMgeneCount`: integer number of genes per pathway/enrichment
   
   * remove steps that create embedded Cnet and Emap `igraph` objects:
   
      * `multiCnetPlot`, `multiCnetPlot1`, `multiCnetPlot1b`, `multiCnetPlot2`
      * `multiEnrichMap`, `multiEnrichMap2`
      * `multiEnrichDF` - consider saving `data.frame` with clear name
      * `multiEnrichResult` - What content is stored here?


## 12jul2023

* New object classes:

   * `"mem"`: store output from `multiEnrichMap()`
   
      * (Essentially a formal replacement for `list` format used currently.)
      * memIM
      * enrichIM, enrichIMcolors, enrichIMdirection
      * geneIM, geneIMcolors, geneIMdirection
      * geneHitList
      * colorV
      * params:
      
         * p_cutoff (from argument `cutoffRowMinP`)
         * min_count
         * topEnrichN
         * pvalueColname
         * directionColname
   
   * `"mem_plots"`: store `mem_plot_folio()` data for re-use.

* New wrapper function `multienrich()` to replace `multiEnrichMap()`?

   * streamlined refactor and replacement for `multiEnrichMap()`
   * avoids defining Cnet and Emap data, pushing into `mem_plot_folio()`
   * reduces arguments by removing all visualization arguments
   * consider storing `memIM` data as `SummarizedExperiment` for convenient
   use with `ComplexHeatmap::Heatmap()` via `jamses::heatmap_se()`.
   * returns `mem` object.

* Update `mem_plot_folio()`

   * input `mem` object
   * add `mem2emap()` plot output.
   * consider using `jamses::heatmap_se()` for heatmap functions
   
      * it adds `jamses` as dependency, along with its dependencies
      * it puts pressure to refactor the jamses contrast stats code
      * should `jamses::heatmap_se()` be moved to new package focused
      only on SummarizedExperiment heatmaps?
      * See Bioconductor package `sechm` which is much less capable, but
      has the same inspiration.
      It provides **row scaling** (ack) as a recommended (!) option.
      (Problematic, sorry to say. For gene expression
      data, the magnitude of change is important, and matters.
      To rescale the numeric range for **consistency** is counter to the
      biology, and to the technology. The technology has measurement
      limitations, for which seeing the **actual** numeric differences is
      important for assessing whether changes are reliable from that platform.
      These differences also transfer into follow-up confirmation assays,
      where changes below a threshold are not feasible to confirm.
      In general, gene (transcript or protein) expression fold changes are
      relatively consistently measured for each gene, so to **enhance**
      the apparent fold change of one gene to fit the fold change of another
      gene is not necessary, certainly not by default.)

* Consider new R package for SummarizedExperiment heatmaps

   * move `jamses::heatmap_se()` into this package
   * move `platjam::design2colors()` into this package
   * minimize R package dependencies

* R packages to review

   * Bioconductor package `"ConsensusClusterPlus"` which can be
   used to determine appropriate `k` values for k-means clustering,
   with metrics to assess consistency of cluster assignment.
   * `"GeneTonic"` function `ggs_graph()` produces a Cnet plot
   which they export to `visIgraph()` for interactivity;
   `enrichment_map()` creates EnrichMap.
   * `"monaLisa"` - motif enrichment, extends HOMER theme; nice heatmap
   Motif labels
   * `"profileplyr"` - coverage heatmap using tidyverse plyr syntax


## 14jun2023

* DONE: `reorder_igraph_nodes()`

   * When `orderByAspect=TRUE`, and it detects tall-skinny aspect ratio,
   it appears to be applying the y-axis sorting bottom-to-top instead of
   top-to-bottom.
   * The culprit was `spread_igraph_labels()` default argument,
   changed from `nodeSortBy=c("x", "y")` to `nodeSortBy=c("x", "-y")`.
   
* `mem_gene_path_heatmap()`

   * Slightly increase the spacing between heatmap body and row/column
   annotations. Currently the gap between heatmap row/column split
   is identical to the gap between heatmap and row/column annotations,
   which makes it harder to distinguish one from the other.
   * The same can be accomplished using `ComplexHeatmap::ht_opts()`
   but the option is hard to remember, and would ideally need to be
   set back to the previous value after drawing the plot.

## 05jun2023

* `mem_plot_folio()`

   * DONE: The enrichment dot plot (or enrichment heatmap) should be created
   after the gene-path heatmap is created, in order to define
   pathway clusters using gene content, rather than defining pathway
   clusters using the `-log10(Pvalue)` matrix.
   
      * The clusters derived from the P-value matrix are sometimes
      not very similar to the gene-pathway clustering result.
      * For this workflow, it makes the most sense to define pathway clusters
      upfront, then share those pathway clusters with all downstream
      visualizations.
      * Consider creating object class "mem_plots".

* consider new function to convert IPA enrichment data to `geneHitList`

   * list element `"Analysis Ready Molecules"` is provided in the IPA output,
   and this data can be used to re-create the directional hit matrix used
   during import (if directionality such as fold change was provided to IPA).

* consider new function to evaluate gene-pathway heatmap output

   * The driving use case is selecting `pathway_column_split=5` upfront,
   but realizing perhaps 3 clusters would be preferred based upon the
   content.
   
      * Criteria for collapsing two pathway clusters together:
      either Jaccard similarity above 0.4, or correlation above 0.6.
      * Definitely requires more testing to determine appropriate
      default thresholds, or whether a reasonable data-driven threshold
      can be defined.
   
   * sometimes two clusters can and perhaps should be merged together
   
      * Create collapsed incidence matrix,
      * Calculate correlation,
      * Any two clusters with correlation >= 0.2 could be merged?
      * Alternatively, if more than 3 clusters would be merged, cancel
      to prevent merging too many clusters together.

* consider new function to edit vertex attributes?

   * use case: existing `igraph` object nodes have attributes to modify:
   
      * sometimes modify only certain nodes: `nodeType="Set"` to edit labels
      * attributes in atomic vector form
      * attributes (`pie.color`, `pie.border`) in `list` form
      * also sometimes want to adjust colors - probably separate function
      * function returns `igraph` object with attribute modified and stored
      in the same state as present originally (e.g. atomic remains atomic;
      list remains list).
   
   * Which format seems most reasonable?
   
      * ```
      gsub_vertex(g,
         pattern_l=list(
            nodeType=c(Set="^(WP|KEGG|BIOCARTA|GO|REACTOME)_")),
         replacement_l=list(
            nodeType=c(Set="")))
      ```
      
      * ```
      gsub_vertex(g,
         subset_attr="nodeType",
         subset_attr_value="Set",
         pattern="^(WP|KEGG|BIOCARTA|GO|REACTOME)_",
         replacement="")
      ```
      
      * ```
      gsub_vertex(g,
         subset_attr_l=list(nodeType=c("Set")),
         pattern="^(WP|KEGG|BIOCARTA|GO|REACTOME)_",
         replacement="")
      ```

* consider new function to adjust `igraph` node colors in all forms

   * modify all fill color attributes `color`, `pie.color`, `coloredrect.color`
   * modify all border attributes `frame.color`, `pie.border`,
   `coloredrect.border`
   * adjust `frame.width`, `pie.lwd`, `coloredrect.lwd` relative to each other.
   For example when `pie.lwd` and `pie.border` are both defined
   (and not transparent), `frame.width=0.1`, otherwise `frame.width=2`.


## 02jun2023

* `mem_plot_folio()`

   * Debug edge cases where `pathway_column_split=4` does not match
   the output number of column split following hierarchical clustering.
   * Debug edge case where there are not enough pathways or genes
   to support the gene-pathway heatmap workflow. E.g. only 1 or 2 pathways,
   or only 1 or 2 genes.

* `multiEnrichMap()`

   * Debug error "multiEnrichMap(): geneHitIM does not contain 5 
   rows present in geneIM, default values will use 1."
   
      * Apparently when some `rownames(mem$geneIM)` are not found in
      `rownames(geneHitIM)`.
      * Ideally change `geneIMdirection` so it does not store "+1"
      and instead shows zero or NA. Genes should be shown but without
      associated direction.
      * The message should describe how to find missing values so the
      user can debug the error.


## 01jun2023

* `reorder_igraph_nodes()`

   * Currently `orderByAspect=TRUE` will order nodes based upon the
   aspect ratio of each nodeset.
   * Ideal world: when nodes are sorted by something like color,
   calculate the aspect ratio of nodes within that color.
   Situation: Assume a tall-skinny nodeset, sorted top-to-bottom by
   different colors. One color has 12 entries, the nodeset is 12 nodes wide,
   so this color appears horizontal among the other nodes. When sorting
   by border color, it goes top-to-bottom, which is not visually intuitive.

* `apply_cnet_direction()`

   * DONE: changed default `col` to use colors:
   `c("blue", "grey80", "firebrick3")` with breaks `c(-1, 0, 1)`.

* change all `frame.lwd` to `frame.width` before `frame.lwd` is widely used.

   * consider backward compatibility: when `frame.lwd` is defined in an
   `igraph` object, copy its values into `frame.width`, then proceed
   using `frame.width` for all other operations.

* `mem_legend()`

   * DONE: new argument `pt.lwd=2` to control the line width used only for point
   borders, useful when `do_direction=TRUE`. The argument `pt.lwd` already
   gets passed to `legend()` however making it a formal argument here
   helps make the option more clear for users.
   * Auto-detect whether to enable `do_direction=TRUE`, by checking if
   any `frame.color` or `pie.border` are defined with red/blue colors.
   Bonus points for using the same colors defined in `frame.color` or
   `pie.border`, however that may be risky if those colors vary based
   upon some fold change value, or vary based upon contrasting with
   the node or pie fill color.

* `jam_igraph()`

   * Consider handling `V(g)$shape="circle"` as shape `"jampie"` during
   rendering. Certain older Cnet plot `igraph` objects appear to break
   the default `shape="circle"` rendering, some cryptic error about
   `names()` not being defined when expected. The error does not appear
   with `igraph::plot.igraph()`, so it is specific to `jam_igraph()`.
   
      * Possible workaround is to sidestep the problem by re-using the
      same `shape="jampie"` rendering method already implemented, which
      would keep all borders consistent when displaying a mixture of
      `shape="circle"` and `shape="jampie"` nodes.
      * Problem with that workaround, if a node has `pie.color` defined
      it will be used when `shape="jampie"` even if the user specified
      `shape="circle"` (which should only use `color`).
      In that case, when a `shape="circle"` node is being rendered internally
      as `shape="jampie"` it should first copy `color` into `pie.color` for
      those nodes; `frame.color` to `pie.border`; and
      `frame.lwd` (`frame.width`) to `pie.lwd`.

## 24may2023

* DONE: Fix bug with node rendering, caused by recent version of `igraph`
adding `vertex.frame.width` (and not `vertex.frame.lwd` ugh).

   * `NULL` or missing `vertex.frame.width` causes an error.
   Ultimately caused by no default value defined in the custom function
   `default_igraph_values()`, which was necessary to create since
   `igraph` does not export that function.
   * FIXED by adding `vertex.frame.width` to default values.
   * Longer term fix is to replace all references to `vertex.frame.lwd`
   with `vertex.frame.width`, before the precedent is set.

* Fix errors caused by `"stringsAsFactors=TRUE"`

   * DONE: `rank_mem_clusters()`

* Fix errors caused when there is only one (or zero) genes.

   * `mem_gene_path_heatmap()` and `mem_plot_folio()`


## 20apr2023

* `mem_enrichment_heatmap()` color legend changes:

   * Show actual P-value `c(1, 0.05, 0.01, 0.001, 10^-4, etc.)` using
   `expression` for labels, and continue using `-log10(p)` for color
   assignment.
   * Add this argument `heatmap_legend_param=list(break_dist=1)` which
   causes the numeric labels to be evenly spaced, instead of having the
   labels at uneven intervals, often with angled lines connecting to
   the color legend.
   * Option for discrete color legend? I.e. Show colors only at the labels,
   and not show the intervening gradient. It is more difficult to show
   abrupt transitions, e.g. it would need to show `c(1, 0.051, 0.05, 0.01)`
   in order to show that `0.051` is not colorized, but `0.05` is colorized.
   The smooth gradient is actually more effective at conveying that effect
   without additional labels.

## 31mar2023

* Nodes with `shape="jampie"` and `frame.lwd=0` are still rendering the
frame color outside the inner borders. When `frame.lwd=0` there should
be no frame drawn even when `frame.color` is defined with a color.

* Big picture musing: Consider replacing base R plotting functions
with corresponding `grid` functions.

   * The `vwline` R package (P. Murrell) is capable of drawing
   internal/external lines.
   * The `gridGraphics` package also provides better methods of
   clipping curved lines to the edge of a node border for example.
   * Major downside, it would likely involve rewriting all the `igraph`
   node shape functions into corresponding `grid` format. It is effectively
   similar to repeating much of `ggraph`, except that this approach
   can be customized. The `ggraph` approach in `ggplot2` is more or
   less untouchable in terms of providing customization. Ugh.
   * Another option may be to figure out how to add custom node shapes
   to `ggraph`, then use `vwline`/`gridGraphics` for rendering. Also
   need to write custom edge bundling function, since those in `ggraph`
   are wholly insufficient. Yeah, this is a no for now.

## 26feb2023

* `reorder_igraph_nodes()`, `reorderIgraphNodes()`

   * DONE: method to specify specific nodes or nodesets to be reordered
   * motivation is to allow sorting based upon relative aspect ratio of nodes,
   so a nodeset whose nodes are "tall-skinny" can be sorted top-to-bottom,
   and nodeset whose nodes are "short-wide" can be sorted left-to-right.
   Frankly, not sure if the inconsistency works for all network layouts,
   but for sure the top-to-bottom is not ideal for "tall-skinny" nodesets,
   it is not visually intuitive.

* consider new igraph shapes, intended to enable inner/outer border, `frame.lwd`

   * Do these make sense?
      * `shape.jamcircle.plot()` - enable custom `frame.lwd` for shape="circle"
      * others: square, csquare, rectangle, vrectangle

* `label_communities()`

   * generalize this method to determine keywords most represented in any
   set of pathway names.

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
