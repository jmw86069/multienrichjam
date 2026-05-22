# Changelog

## multienrichjam 0.0.117.900

### New functions

- `EnrichmentMap()` as proper method for `Mem` and
  [`MemPlotFolio()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  objects. Internally it calls
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  then applies optional arguments to define node groups. When using
  `MemPlotFolio` it allows clustering by `Clusters()`, while `Mem` input
  only supports network communities. Node groups supported: clusters,
  cluster_labels, communities, or none.
  [`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
  now calls `EnrichmentMap()` directly.

### Bug fixes

- edge bundling error fixed when points were co-linear.

## multienrichjam 0.0.116.900

- Hotfix: Added ‘BiocStyle’, ‘kableExtra’ dependencies.
- Hotfix: Removed requirement for directionColname in IPA data, common
  for older format IPA files as included as *sigh* package data.
- Hotfix: Un-hid the initial example code chunk.
- Added dependencies: ‘withr’, ‘ggforce’, ‘circlize’,
  ‘shinydashboardPlus’.
- Support for MemPlotFolio cluster labels, and cluster data. The cluster
  labels are intended to support short summary title or phrase,
  addressing the common “What’s that cluster?”
- New
  [`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
  paradigm for MemPlotFolio, supporting Rmd and Qmd tabbed and
  non-tabbed output.

### Changes

- [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  and
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Quietly adds ‘Mem’ to metadata(Mpf)\$Mem for use in other
    MemPlotFolio related functions. Experimental.

- `MemPlotFolio`

  - New accessor and setter
    [`ClusterLabels()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
    stores custom labels for each pathway cluster.
    [`ClusterData()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
    stores data in `list` format, intended to store summary information
    for each cluster.
  - Slot ‘metadata’ now enforces
  - ‘cluster_labels’ as named `character` vector of labels, named by
    `names(Clusters(Mpf))`. These names are short, usually single-letter
    codes for each cluster, while the labels are words or phrases, which
    may appear in
    [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
    [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
    or potentially the cluster title on
    [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md).
  - ‘cluster_data’ is a named `list` of data, content to be determined.

- [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - ‘type’ argument recognizes: ‘cluster’, ‘set’, ‘title’. A suffix ‘2’
    hides gene node labels, for example `type='set2'`. The ‘cluster’
    uses
    [`ClusterLabels()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
    if present. The ‘set’ uses the top N pathways, shortened to M
    characters.
  - Arguments in ‘…’ are used: ‘layout’, ‘rotate_degrees’, ‘width’
    applies word wrap, ‘maxNchar’ applies label cropping.

- [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

  - Default for argument ‘main’ will use
    [`ClusterLabels()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
    as plot title when present.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - `mark.groups` is a little more robust to subset igraph, and
    cluster_names as vector, list, or NULL.

- [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
  [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  now call `ComplexHeatmap::draw(x, newpage=FALSE)` when running inside
  Positron, and not running inside knitr, to prevent the
  [`grid::grid.newpage()`](https://rdrr.io/r/grid/grid.newpage.html)
  from causing Positron to ignore subsequent figure grobs. More testing
  necessary to confirm no other conditions inside Positron are affected.

### New functions

- [`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)

  - Distinct function from
    [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    to enable plotting one or more plots from the `MemPlotFolio` (Mpf)
    object.
  - Enables optional Rmarkdown (Rmd) or Quarto (Qmd) tabsets, including
    Quarto open/close style, and sub-tabs where relevant.
  - In the next release, this function will likely replace all active
    plotting functions currently done within
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md).

## multienrichjam 0.0.115.950

### Bug Fixes

- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  fixed visual glitch in the color legend causing node fill color to
  appear as fill color for directional borders. Also fixed alignment
  using ‘added-top’ and ‘added-bottom’ for multi-column output.
- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  used to render pie nodes, updated to use proper `frame.lwd` line
  widths for each node.
- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  and
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  now define ‘vertex.size2’ used for igraph node
  shape=‘coloredrectangle’. Rarely used, but an OG that may be preferred
  when analyzing a lot of enrichments.
- [`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)
  fixed bug causing blank cells to re-use recycled border colors.

### Changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  and
  [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  have new option to print Rmarkdown (Rmd) tab headers, using
  `do_rmd_tabs=TRUE`. The original target for
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  was multi-page PDF, however for HTML output the Rmd tabs are a
  convenient option.

- Several functions now use ‘vertex.frame.width’ instead of
  ‘vertex.frame.lwd’ as new attribute name since igraph 1.3.0:
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
  [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md),
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

- [`shape.jampie.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  now clips to the outer border instead of the edge of the inner border.
  Mostly noticeable with transparent outer border, or wide border where
  it was clear the arrow head was too close and was cropped.

- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  node border line width is now scaled by node size, relative to the
  median node size. Otherwise, small nodes only showed the border color.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - now applies `node_factor` to borders ‘frame.width’ and
    ‘pie.border.lwd’ by default. Previously, nodes would be smaller but
    the border maintained absolute width and could fill the node
    completely. To reverse, apply the inverse via `border_factor`.

- [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
  [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
  [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  have new arguments ‘legend_x’, ‘legend_y’ which are passed to
  [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  to control the position of the legend.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Sped Cnet exemplar and cluster generation in large networks by
    skipping the global Cnet step, which is not necessary.

- [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)

  - Changed default ‘byCols’ to favor ‘minp_rank’ instead of
    ‘composite_rank’.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  new argument `border_factor` to adjust frame and border line width.

## multienrichjam 0.0.114.950

- Minor update to improve defaults for
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

### Updates to existing functions

- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

  - New arguments: ‘apply_edge_width=TRUE’, ‘apply_edge_color=TRUE’,
    ‘max_edge_width=25’ apply edge color and width by default.
  - Node size calculation improved to set `median_size` and `max_size`.

## multienrichjam 0.0.114.900

- Changes to support `gseaResult` input to
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
  and to support it alongside `enrichResult` when relevant.

### Bug fixes

- Fixed bug in
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  caused by its use of
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  calling
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).
  When it could not define Cnet nodesets, singlets could be NA in the
  updated layout.

### Changes to existing functions

- [`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
  fixed when input list had different colnames.

- [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)

  - Much more robust in checking layout rownames, and the order of
    rownames in the event the igraph was subset and the layout was not.
    It only works when there are rownames(layout), which is not the
    typical igraph default.

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - With `gseaResult` input, it calls internal
    [`gr2er()`](https://jmw86069.github.io/multienrichjam/reference/gr2er.md)
    which populates the common colnames into the `gseaResult` object for
    convenient use alongside `enrichResult` objects.
  - Arg ‘directionColname’ defaults now include `'NES'` to use GSEA
    normalized enrichment score (NES) for enrichment directionality.
  - Internally, the various ’\*Colnames’ args are applied to each entry
    in `enrichList`, to tolerate different colnames. This change
    embraces situations like having both `enrichResult` and `gseaResult`
    together, and any non-clusterProfiler enrichment table. Previously,
    colnames must be consistent across inputs.
  - When the `enrichList` tables contain different colnames, it now goes
    through each entry and updates colnames to match the first matching
    colname of each type.  
    The updated enrichments can be reviewed with `enrichList(Mem)`.
  - Work is active and ongoing to handle multiple types of enrichment
    data.

## multienrichjam 0.0.113.900

### Changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - When ‘enrichList’ is supplied with `data.frame`, the call to
    [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
    now passes proper arguments in addition to ‘…’, which permits
    columns ‘pvalueColname’ to be respected.

## multienrichjam 0.0.112.900

- Added ‘readr’ to Enhances, it can be used as a backup import option
  for now.
- Removed some dependencies internally. More to be done.

### Updates to existing functions

- [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Added formal arguments to the definition, instead of passing ‘…’ to
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md),
    so they appear with [`args()`](https://rdrr.io/r/base/args.html),
    [`jamba::jargs()`](https://jmw86069.github.io/jamba/reference/jargs.html),
    and R argument auto-fill.

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - Argument ‘removeGrep’ accepts a vector, and will insert delimiter
    ‘\|’.
  - Argument ‘removeGrep’ includes more common pathway prefixes, formats
    including: Wp:Pathway, Wp12345:Pathway, Hsa:Pathway,
    R-Hsa-152:Pathway, M12345:Pathway. Also recognizes some common
    species codes: Hsa\|Mmu\|Rno\|Dme.
  - New argument ‘makeUnique’ optionally calls
    [`jamba::makeNames()`](https://jmw86069.github.io/jamba/reference/makeNames.html)
    to enforce unique names after adjustments are made. It is fairly
    common that several pathway sources provide the same pathway, for
    example “T-cell signaling”.
  - All `Mem` input automatically uses `makeUnique=TRUE` to protect the
    integrity of the output `Mem` object.
  - Added tests to cover new removeGrep, and makeUnique=TRUE.

## multienrichjam 0.0.111.950

### Updates (Quick Fix, still in progress)

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - Argument ‘removeGrep’ now removes the longer KEGG Medicus prefixes:
    “KEGG_MEDICUS\_”, “KEGG_MEDICUS_REFERENCE\_”,
    “KEGG_MEDICUS_VARIANT\_”, “KEGG_MEDICUS_PATHOGEN\_”. It does not
    (yet) remove “KEGG_MEDICUS_ENV_FACTOR” because what is that?
    Environmental factor pathways? We will revisit in future.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  fixed bug where `revert_ipa_xref` was not propagated to multi-sheet
  xlsx import consistently.

## multienrichjam 0.0.111.900

### Updates

- words.txt,abbrev.txt

  - “Signaling by X Y” is only abbreviated to “X Y” with two or more
    words.
  - Added some capitalizations from recent IPA data: BBSome, ID1,
    GABAergic, PDGF, PDGFR, PDGFRB.
  - IL 24 is edited to IL-24; Ifn to IFN; Robo and Robos to ROBO and
    ROBOs; Nmd to NMD; Pml to PML; Hes Hey to HES/HEY; Slits to SLITs;
    Cmyb to c-Myb; Ctcf to CTCF; Gna12 13 to GNA12/GNA13; Rhoa to RhoA;
    Hcmv to HCMV; Pirna to piRNA; Classiii to ClassIII; Ndkdynamin to
    NDK/Dynamin
  - Fixed “Signaling By” so it is only removed when one word follows it,
    e.g. “Signaling By VEGF” is kept as-is, but “Signaling By Erbb2 in
    Cancer” is converted to “Erbb2 in Cancer”.

## multienrichjam 0.0.110.900

### Bug fixes

- [`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md)

  - Fixed error when input `list` did not contain all nodes.
  - Fixed bug causing mis-alignment of colors across nodegroups.
  - Made argument ‘nodegroups’ more robust, calling
    [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
    as needed.
  - Added examples.

- [`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md)
  fixed missing argument ‘nodegroups’.

- [`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md)
  updated to support one-point hulls.

### Changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  passes ‘…’ to
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  when used with `mark.groups=TRUE` to make it easy to pass `min_size`
  to filter nodesets to at least this many nodes.

- [`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md)

  - more robust to argument ‘nodegroups’
  - new adjustments ‘darkFactor’,‘sFactor’ to distinguish from node
    colors.
  - added examples

- [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)
  now always returns colnames ‘x’,‘y’,‘z’ as many as needed.

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - Now produces a point hull for 1 and 2 point data, as originally
    intended.
  - Adds argument ‘min_points’ to allow filtering for 3+ points, for
    legacy behavior.

- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)

  - now tolerates data without ‘nodeType’, but does not filter ‘Set’
    nodes of course.
  - new argument ‘min_size’ to require nodesets to have at least this
    many nodes. Intended mainly for use with
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md).

- [`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md)

  - now tolerates argument nodegroups as NULL, which calls
    [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md).
    Help docs clarify that `filter_set_only=FALSE` will assign a nodeset
    to ‘Set’ nodes, otherwise they are NA.

## multienrichjam 0.0.109.900

### Changes

- Added more to words.txt and abbrev.txt.

## multienrichjam 0.0.108.900

### Changes

- Added more to words.txt and abbrev.txt.

- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

  - New arguments ‘vertex.label.font=2’ and ‘use_shadowText=TRUE’ to
    match
    [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
    default aesthetics.
  - Edge and node attributes include ‘gene_count’ for reference.
  - New argument ‘size_by_genes=TRUE’ and ‘mean_size=5’.

- [`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

  - Now supports other `igraph` objects, such as
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
    output.
  - Supports graph attribute ‘mark.groups’ as populated by
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).
  - Argument ‘nodegroups’ can include a `communities` object, for
    non-Cnet.

- [`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md)
  now supports `list` input which would not need to be converted, as a
  convenient way to confirm proper `list` output.

### Bug fixes

- [`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md)
  added missing prefix to
  [`jamba::tcount()`](https://jmw86069.github.io/jamba/reference/tcount.html).
- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md):
  fixed error when there were no edges, now returns NULL.
- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md):
  fixed error when no edges, now returns NULL.
- [`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md):
  fixed error for smaller Cnet graphs, added tolerance for disconnected
  subsets with larger Jaccard values.

## multienrichjam 0.0.107.900

### Bug fixes

- Visual glitch fixed when using `mark.groups` or `nodegroups` with
  named entries, causing the edge bundling to go awry. Group names often
  clashed with singular node names, causing the group centers to be
  mis-matched. It was usually “nearby” so it was not terrible.
- `headers()` now returns the `list` and not just names of the list.

### Changes

- ‘MemPlotFolio’ metadata now includes `logical` ‘hasDirection’ as
  convenient way to tell if data include directionality.
- [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
  [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
  [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  show direction in the color legend if defined in the ‘MemPlotFolio’
  object. It can be overridden with `do_directional=FALSE`.
- Added more words.txt.

### New functions

- [`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

  - Experimental layout algorithm for Cnet plots, intended to help each
    nodegroup show visual grouping. Still being improved.
  - It helps when two nodegroups are intertwined, with no other force
    helpful in separating the groups. Since nodegroups do not have
    connections among themselves, typical layout algorithms are less
    effective.

- [`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md)

  - Returns a vector of nodegroups for all nodes in a Cnet `igraph`, a
    simple wrapper around
    [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
    which returns a `list`.

## multienrichjam 0.0.106.900

### Changes

- Updated ‘words’ and ‘abbrev’ to remove verbose phrases from gene sets.
  “Genes encoding proteins involved in” can be entirely removed.

### Changes to existing functions

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - New argument ‘perl=TRUE’ to allow disabling the word boundary
    condition.
  - Added testthat conditions to verify behavior.

- [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)

  - Now uses vertex attributes x,y,z when `layout` is not provided, and
    when graph attribute ‘layout’ is not defined, and when x,y,z
    attributes are defined with `numeric` values. The igraph team has
    debated removing support for graph attributes, using only vertex
    attributes, in order to improve subgraph support for layouts.
  - New argument `default_layout` applied only when layout does not
    exist.
  - Argument `make_circular` removed, as it was not used.

### New functions

- [`set_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)

  - complements
    [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)
    by also storing layout in the `igraph`. Intended to future-proof the
    layout, since the igraph team is debating whether to recommend only
    using vertex attributes, or graph attributes.
  - It can optionally call
    [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
    to update label positions.

## multienrichjam 0.0.105.900

### Changes to existing functions

- `pathway_column_split` and `column_split` accept `list` for
  convenience.

- Added default `vertex.label.font=2` (bold labels) for
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  recognizes graph attribute ‘use_shadowText’ as global option.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md),
  [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Passing `pathway_column_split` as `list` is now supported.
  - Default annotation width is 5mm instead of 6mm.
  - Anno size is now correctly used for column annotations.
  - When providing ‘mpf’ it will re-use Clusters() for consistency.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  now applies two new defaults:

  - ‘vertex.label.font=2’ for bold node labels by default
  - ‘use_shadowText=TRUE’ added to graph attributes, for outlined labels

## multienrichjam 0.0.104.900

- Added ‘openxlsx’ to Suggests.
- Removed ‘matrixStats’ dependency.
- Added quick walkthrough to README.Rmd as visual intro of the package.

### New functions

- [`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md):
  creates a gene hit list from IPA data, for
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  argument ‘geneHitIM’.
- [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md):
  for cleaner workflow to create `MemPlotFolio` without directly
  plotting all results.
- [`score_gene_path_clusters()`](https://jmw86069.github.io/multienrichjam/reference/score_gene_path_clusters.md):
  Experimental, intriguing effort to find “hot spots” in the Gene-Path
  Heatmap clustering.
- [`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
  [`local_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md):
  applies withr style options using
  [`ComplexHeatmap::ht_opt()`](https://rdrr.io/pkg/ComplexHeatmap/man/ht_opt.html),
  and added to
  [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

### Changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - argument ‘mark.groups’ now accepts `function` and applies
    dynamically, or `TRUE` which marks Cnet nodesets when applicable.
  - Fixed regression that improperly sized `mark.expand` causing it to
    be 200 times smaller than intended.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - now calls
    [`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md)
    by default.
  - fixed regression where batch import would ignore
    `revert_ipa_xref=TRUE`.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)

  - new argument ‘forceColors’ to use traditional clusterProfiler
    colors, only coloring set nodes as palegoldenrod, and gene nodes as
    grey. It skips all other directional and pie-based styling.
  - supports enrichIMdirection, geneIMdirection, via new arguments
    ‘direction_col_fn’, ‘direction_cutoff’, ‘direction_max’, and by
    default uses
    [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
    to apply divergent color gradient with optional numeric floor. As a
    result, directional scores (often a type of z-score) will only be
    colorized at or above that threshold. Convenient to color only
    z-scores above 1 or below -1.

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)

  - applies border colors with vectorized logic, substantially faster.

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - includes `direction_cutoff` in
    [`thresholds()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.html)
    to use by
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    which calls
    [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md).
    This threshold determines when to colorize set nodes with up/down
    directionality.

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)

  - now calculates aggregate directional scores for
    [`enrichIMdirection()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
    previously it used the aggregate
    [`enrichIM()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
    values.
  - the
    [`enrichIMdirection()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
    values are now used to colorize Set nodes, using `direction_cutoff`
    as a minimum to apply up/down color.

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - Now uses package data `words` and `abbrev` as `data.frame` for
    default replacements, to make it easier to edit over time. Data are
    prepared from files stored in ‘extdata/words.txt’ and
    ‘extdata/abbrev.txt’.

## multienrichjam 0.0.103.900

### Major/Breaking changes

- Added `MemPlotFolio` S4 class as default output from
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md).
  Nobody in theory should have been using
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  output directly. However, you can use `returnType="list"` for legacy
  output.

- There are several accessors, a few of which can plot the data, and
  return the corresponding graphical objects to be customized.

  - [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  - [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  - [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  - [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  - [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

- Vignettes and tests have been updated to use the new MemPlotFolio S4
  functions, a much cleaner approach overall.

### New Mem methods

- [`geneInCategory()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  mimics
  [`DOSE::geneInCategory()`](https://rdrr.io/pkg/DOSE/man/geneInCategory.html)
  by listing genes in each pathway.
- `genesBySet()` lists pathways associated to each gene.

## multienrichjam 0.0.102.900

- Added “Versioned” to ‘Mem’ class, to handle changes over time.
- Added more S4 “setter” methods, to keep rownames geneHitIM consistent
  with geneIM since they may not be identical.
- Updated `Reese_genes` and `Memtest` with Entrez ID and updated
  symbols.
- Re-organized the pkgdown function families, trying to make it simpler.

### General changes

- Added ‘Biobase’ to Imports, to use class ‘Versioned’.
- Added ‘methods’ to Imports, imported several functions.
- S4-related code was moved to prepare for Bioconductor: ‘AllClasses.R’,
  ‘AllGenerics.R’, ‘methods-Mem.R’.
- Added ‘S4Vectors’ dependency, then imported
  [`as.list()`](https://rdrr.io/r/base/list.html) to avoid weird C stack
  error by recursive loop, when unloading multienrichjam. Reproducible
  by unload, but seemed to occur in other rare cases, e.g. unit testing,
  creating vignettes, pkgdown, perhaps others. This approach seems to
  load S4 method appropriately upfront.

### Mem S4 changes

- `Mem-class` now contains `'Versioned'` Biobase class, which assigns a
  class version to help with future compatibility. It should help
  convert old versions of `Mem` objects to the current format, allowing
  me to change the format over time as needed.
- Added `updateObject` generic method for `Mem` S4 objects. It currently
  adds a version to any un-versioned S4 Mem object.
- Updated `genes<-` so it also updates `geneHitIM(x)`, and gives a
  warning if for some reason not all `genes(x)` are present in
  `rownames(geneHitIM(x))`.
- Added `geneHitIM<-` which also confirms all `genes(x)` and
  `enrichments(x)` are correctly represented.
- Similarly added: `geneIMcolors<-`, `geneIMdirection<-`,
  `enrichIMdirection<-`, `enrichIMcolors<-`
- Added ‘Mem-slots’ help topic, to describe the `Mem` slots even though
  it is not intended for end users. The move was motivated to keep the
  `Mem-class` page clean.

### Changes to package data

- `Reese_genes` was updated to assign Entrez gene ID values as
  [`names()`](https://rdrr.io/r/base/names.html). The Entrez ID is
  useful for other pathway enrichment methods, such as clusterProfiler
  and ReactomePA functions. Three genes were updated to use more current
  Entrez gene symbols, as of November 2025: BC034424 -\> HEXA; INADL -\>
  PATJ; SLC24A6 -\> SLC8B1.

- `Memtest` was updated:

  - to use the full `Reese_genes` incidence matrix for
    `geneHitIM(Memtest)`.
  - to apply updated `Reese_genes` to `genes(Memtest)`, and `enrichList`
    ‘geneID’ columns.

## multienrichjam 0.0.101.900

### Breaking changes

- Migrated to the new S4 `Mem-class` object! Most functions accept
  either `Mem` or legacy `list` mem format.
- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  now returns S4 `Mem` object by default, via `returnType="Mem"`. The
  legacy `list` mem format is available with `returnType="list"`.
- [`Mem_to_list()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  converts `Mem-class` to legacy list mem format.
- Removed unofficial functions which were never used:
  filter_mem_genes(), filter_mem_sets(), subset_mem().

### Updates

- Added vignette for using `clusterProfiler`, still work in progress.
- Added test data objects: `Memtest` S4 object, and `Reese_genes`.
- Removed mention of arules transactions objects.
- Moved ‘Mem Plot Concepts’ into its own minimal vignette, still in
  progress.
- Added ‘tidyr’ dependency, removed ‘reshape2’.
- Added ‘vdiffr’ to Enhances, and numerous visual unit tests.
- Added ‘lifecycle’ dependency to help manage deprecated features and
  arguments.

### changes to existing functions

- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)

  - new argument `readable=NULL` which when NULL will set the `logical`
    value accordingly for the output `enrichResult` object, mostly
    useful when using output with `clusterProfiler` related functions.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - Now expects Mem S4 object input, but accepts `list` legacy format.
  - Refactored logic for column_split, row_split, column_title,
    row_title, cluster_columns, cluster_rows.Now accepts integer,
    vector, data.frame, function inputs. Added docs to explain.
  - Now it tests the clustering before applying
    [`cutree()`](https://rdrr.io/r/stats/cutree.html) to confirm the
    data supports the requested number of clusters - because with binary
    matrix data, sometimes the distance is literally zero, preventing a
    cut.
  - Fixed error when number of clusters requested using `column_split`
    or `row_split` did not match the number of title entries.

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Default output is now S4 ‘Mem’ object, set with `returnType="Mem"`.
  - Legacy behavior of creating the total Cnet plot objects, and
    EnrichMap objects, are now only performed when `returnType="list"`.
    These plots were no longer used, in favor of
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    or
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).
  - Added argument `p_cutoff` and deprecated `cutoffRowMinP`.
  - Deprecated `overlapThreshold`, only used for EnrichMap which is
    skipped.

- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

  - Accepts `Mem` input.
  - Applies
    [`igraph::layout_components()`](https://r.igraph.org/reference/merge_coords.html)
    to layout multiple disconnected sub-graphs together.
  - New default `do_plot=FALSE` for broader compliance.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)

  - New argument `remove_singlet_genes=TRUE` will hide singlet genes,
    useful when plotting a subset Cnet with specific pathways.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - It will now apply graph attribute ‘mark.groups’ when defined, to
    help re-use this sub-cluster information without having to keep
    adding it to the plot function. To turn off this feature use
    `mark.groups=FALSE`.

## multienrichjam 0.0.100.900

### Updates

- ShinyCat Cnet Adjustment Tool

  - Now highlights edges of the selected node.

- `Mem`, `Mem-class` methods more fully documented.

- Fixed warnings during jampie,coloredrectangle igraph shape rendering.

- Moved dependency on alphahull to Enhances

### New functions

- [`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md)

  - Given one or more nodes, it colors and expands the width of all
    edges.
  - Optionally de-emphasizes non-highlighted edges with alpha
    transparency.

- [`list_to_Mem()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
  [`Mem_to_list()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  for interconversion with S4 ‘Mem’ object.

- [`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
  and internal functions for alpha hull calculation.

## multienrichjam 0.0.99.900

### Updates

- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)

  - Added buttom “Download as .RData” to save for re-use.
  - Added docs describing how to save the environment.
  - Node and label factors are applied to Gene and Set nodes.

## multienrichjam 0.0.98.900

- Added dependencies: shiny, shinydashboard

### New R-shiny app for Cnet plot adjustments

- R-shiny Cnet Adjustment Tool (ShinyCat)

  - [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)
    takes `g` as input, returns an `environment` containing ‘adj_cnet’
    with the adjusted `igraph` object.
  - Intended to make my life a easier when making a zillion tiny
    adjustments to a Cnet plot as a final figure.
  - Nodeset adjustments: x,y coordinates; percent spacing, rotation.
  - Node adjustments: x,y coordinates; label angle; label distance.

### new functions

- [`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md)

  - takes `data.frame` input for nodesets, and nodes, to apply
    adjustments with just one function call.

### updates to existing functions

- [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
  (the internal function)

  - Removed the `rm(x, y)` to remove warning when y does not exist.
  - Call
    [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
    directly rather than overriding
    [`text()`](https://rdrr.io/r/graphics/text.html).
  - New default `mark.expand=NULL` uses half the median vertex.size.

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - Accepts `lwd=0` and hides the border.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - New default: `mark.expand=NULL` uses half the median vertex.size,
    calculated within
    [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md).

- [`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md)

  - Updated help docs, and output data format to return `list` of
    matrices.
  - Default uses node groups (Cnet nodesets), named by nodeset.

- [`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md)
  now calls
  [`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md)

### In progress

- S4 ‘Mem’ object implementation is in progress, some helper functions
  are included for future use.

## multienrichjam 0.0.97.900

### updates to existing functions

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - Fixed error when importing IPA data where some analysis results
    included two or more gene columns, such as “Causal Networks”.
  - Added documentation to `revert_ipa_xref=TRUE` argument default,
    suggesting `revert_ipa_xref=FALSE` for microarray platforms for
    example.

## multienrichjam 0.0.96.900

- Quick update to restore `list2graph()` removed from non-exported
  remote function. Added `list2graph_ggt()`.

## multienrichjam 0.0.95.900

### S4 objects

- `Mem` - Multi-Enrichment

  - Still todo: Convert most functions to import/export `Mem`

### bug fixes

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  fixed exemplar Cnet plots ignoring the `byCols` when choosing the
  exemplar pathway to display.

### changes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Added `"thresholds"` to `list` output.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  more gracefully handles optional data such as “Networks” which do not
  return traditional gene symbols, but a collection of IPA entries with
  no clearly associated cross-reference information. As a result, these
  columns are skipped.

## multienrichjam 0.0.94.900

### changes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Silenced unintended verbose output.

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - Added extensive corrections for common MSigDB abbreviations. Many
    should probably be moved to a file or improved structure, and not
    the function arguments. For now it solves the problem.

## multienrichjam 0.0.93.900

### changes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - New default `min_count=3` instead of previous `min_count=1`.
  - Updated `mem$colnames` to include all detected colnames, adding
    countColname, directionColname, pathGenes, geneHits. The order now
    matches function argument order.
  - Added new component `mem$thresholds` with `cutoffRowMinP`,
    `min_count`, `overlapThreshold` (may not be necessary), and the
    `topEnrichSource*` arguments, sufficient to reproduce the method.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - Changed ceiling for auto-detected pathway clusters to 8, from 10.

- Added more tests, including much of
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

## multienrichjam 0.0.92.900

### bug fixes

- [`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)

  - Fixed bug which occurred only when `enrichList` had length=1, which
    caused gene/pathway to be out of sync. The bug appeared to be caused
    only with single-enrichment input, and must have been relatively
    recent.
  - Added numerous `testthat` entries to cover this bug, and various
    permutations. Will rapidly add more test cases to cover various core
    behaviors and expectations with
    [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

### changes

- moved `enrichList2DF()` to its own .R file.

## multienrichjam 0.0.91.900

### changes

- Moved `.onLoad()` to `zzz.R`

- Moved
  [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  to its own R file for convenience.

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)

  - default argument change `width=40`, previously 25.
  - arguments `words_from`,`words_to` were edited to add more commonly
    observed patterns.
  - new arguments `add_from`,`add_to` for user-defined replacements, to
    use in addition to the default replacements `word_from`,`word_to`.
  - new argument `do_abbreviations=TRUE` to help shorten common phrass
  - new arguments `abbrev_from`,`abbrev_to` are used with
    `do_abbreviations=TRUE`. These are opinionated changes.

## multienrichjam 0.0.90.950

### bug fixes

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
  fixed typo specific to “right” alignment.

## multienrichjam 0.0.90.900

- Added `shadowtext` to Suggests, optionally adds contrast to
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  when using `show_enrich`.

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - help docs changed to use the correct name
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
    from incorrect `mem_gene_pathway_heatmap()`.
  - The `column_anno_padding`,`row_anno_padding` are now properly
    applied before creating each heatmap, then reverted to previous
    state before returning data.
  - The Enrichment Heatmap default now uses `style="dotplot_inverted"`,
    for previous behavior use `style="dotplot"`.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - Fixed visual bug where caption labels indicated “genes” twice,
    instead of “genes” and “pathways”.
  - When `rotate_heatmap=TRUE`, arguments are more consistently flipped
    from row to column, and help docs are updated to be explicit about
    which arguments are flipped, and which two are not.
  - The heatmap caption uses consistent order to describe rows, columns,
    and reverses the order when `rotate_heatmap=TRUE`.
  - The geneIM legend label now uses title case for consistency: “Gene
    Hit By Enrichment”, rather than “enrichments per gene” which was not
    correct wording.

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - Renamed previously unused arg `show` to `show_enrich` so it is
    passed cleanly via `...` using
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md),
    added help doc. It adds optional label with -log10P, z-score, and/or
    number of genes.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - Fixed issue with `"Analysis Ready Molecules"` not having proper
    header column pattern matching, causing some imports to skip
    applying the user-defined gene symbol.

## multienrichjam 0.0.89.900

### changes to existing functions

Added `%>%` to imported functions, and not all of `dplyr`.

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - New argument `cell_size` default `NULL`, with option to specify the
    exact cell size used in the heatmap. The driving use case is to
    define perfectly square heatmap cells for `"dotplot_inverted"` or
    `"dotplot"` so the circles are perfectly centered inside square
    cells. It forces the output figure to be tall and wide enough to
    accomodate the resulting figure, so it needs some user math.
    Implemented from a user suggestion, and it works and looks great. It
    just requires some upfront work to create a figure with large enough
    canvas.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Added proper pathway cluster labels to
    [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
    output.
  - Finally fixed the pathway order of
    [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
    to be identical to
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
    without minor re-ordering caused by default `ComplexHeatmap`
    behavior. Now the two formats have identical order for direct
    comparison.
  - Cleaned up the help documentation a little bit.

## multienrichjam 0.0.88.950

### changes to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - New option `style="dotplot_inverted"` due to a cool recent paper
    using this style. I like it so much, it is the new default. Much
    easier to see the color and the size of the circle, especially for
    cells with very small circles. From Jang et al, the Waggoner lab,
    Nature Genetics 2024: <https://doi.org/10.1038/s41588-024-01880-x>

- [`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md)

  - New argument `invert=FALSE` to control whether to draw colored
    circles, or colored cells with white circles on top (`invert=TRUE`).

## multienrichjam 0.0.88.900

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - Moved workhorse function to its own .R file.
  - Changed default: `verbose=FALSE`.
  - New argument `rotate_heatmap` to handle this workflow properly.
  - New arguments `row_anno_padding`, `column_anno_padding` to control
    the padding between heatmap and row and column annotations,
    respectively.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - The caption is displayed as a ComplexHeatmap Legend, therefore using
    consistent font and alignment with the color legends, rather than
    being tucked into the corner where it sometimes overlapped other
    heatmap labels.
  - Caption now displays the rows/columns counts first. Remaining values
    are more user-friendly.
  - Default caption font size 10 instead of 6, consistent with color
    legend text.
  - Argument `seed` is utilized.

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - New argument default `label_preset=NULL`.
  - Argument `label_preset` is properly recognized.

- Overall: Functions with argument `seed` now properly ignore the seed
  when it is `NULL`, thereby allowing random behavior when preferred.

## multienrichjam 0.0.87.900

The next update will likely use `Mem` S4 object instead of `list`,
although the plan is to allow convenient conversion to `list` for legacy
compatibility. The S4 functions and methods should be more convenient in
long term, and will be Bioconductor-compliant.

### new functions

- [`add_pathway_direction()`](https://jmw86069.github.io/multienrichjam/reference/add_pathway_direction.md)

  - Adds column to `enrichResult` indicating the directionality, using
    the IPA Activation z-score calculation (refs in help doc). It
    requires `gene_hits` as a numeric vector named by gene symbol.

### changes to existing functions

- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md),
  affecting
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md),
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - now applies logic for `descriptionGrep`, `nameGrep`, and
    `subsetSets` using OR logic, so any combination of matching results
    will be retained.
  - Help docs have been updated.

### other changes

- `.onLoad()`

  - Minor change, it now checks whether the new `igraph` shapes already
    exist before adding, mostly useful when reloading this package in a
    live R session, which should be rare.

## multienrichjam 0.0.86.900

### Bug fixes

- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md),
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md),
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Weird rare scenario that appears limited to custom enrichment data
    where input data contains multiple P-value columns, none of which
    match the defaults for argument `sortColnames` in
    [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md).
  - First bug: The `sortColnames` argument was not evaluated with
    [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
    which is intended to match the patterns to actual data. This
    decision was probably to avoid handling the optional prefix `"-"` to
    reverse sort by colname.
  - Second bug: The `sortColnames` should really only use
    `pvalueColname` then reverse order for `countColname`.
  - End result: When `sortColname` matched no existing columns, the data
    was still sorted using
    [`jamba::mixedSortDF()`](https://jmw86069.github.io/jamba/reference/mixedSortDF.html).
    It found no matching colnames so by default it sorted starting with
    the first column. This is incorrect, and caused the bug.
  - New argument default: `sortColname=NULL` will use `pvalueColname`
    then reverse order of `countColname`. Alternatively
    `sortColname=FALSE` causes no sort to be performed, using data in
    the order it was provided. Finally, if `sortColname` is defined, and
    its values do not match existing colnames, no sort is performed.

### changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Minor change to silence the text output when `enrichIM` is entirely
    NA, which can occur when the input data does not contain any
    significantly enriched pathway results. This outcome causes
    [`stop()`](https://rdrr.io/r/base/stop.html) but should not
    otherwise print output unless `verbose=TRUE`.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  and
  [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)

  - Both now use
    [`jamba::gsubs()`](https://jmw86069.github.io/jamba/reference/gsubs.html)
    instead of the “temporary” internal function by the same name.

- [`gsubs()`](https://jmw86069.github.io/jamba/reference/gsubs.html) is
  removed, renamed
  [`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md)
  in favor of
  [`jamba::gsubs()`](https://jmw86069.github.io/jamba/reference/gsubs.html).

### other changes

Two functions are no longer imported, instead they are both called using
the package prefix. It only affects
[`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md)
which is no longer used by default.

- [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
- [`jamba::nameVector()`](https://jmw86069.github.io/jamba/reference/nameVector.html)

## multienrichjam 0.0.85.900

### Changes to existing functions

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - New argument `revert_ipa_xref=TRUE` changes the default behavior
    (for the better) such that the resulting gene symbols associated
    with IPA pathway enrichment will match the input gene symbols,
    instead of using the customized IPA symbols.

    - The situation does not have an ideal solution. Ultimately, IPA
      provides results which do not completely represent the input data.
      When two genes are combined to one entity by IPA, they only retain
      one gene symbol in `"Analysis Ready Molecules"`, and so there is
      no recorded association of all entities which were combined.
      Sometimes the combined symbol is `"HSPA1A/HSPA1B"` which can be
      separated… but they do not indicate which symbol(s) the user
      provided, they only record one. Some entities are called
      `"NBPF10 (includes others)"` and which appears to include
      `'NBPF10"` and `"NBPF19"`, possibly others.

## multienrichjam 0.0.84.900

### Bug fixes

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - Finally “fixed” the full bug, sometimes causing point hull to fail.
    Sometimes
    [`alphahull::ahull()`](https://rdrr.io/pkg/alphahull/man/ahull.html)
    would return weird results with small `alpha` values: individual
    points, segments, or multiple disconnected polygons.
  - Now
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
    calls
    [`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md)
    which performs additional validation checks: Confirms there is
    actually a polygon; confirms each point in the `edges` are used
    exactly twice; confirms that all edge points are used in a
    continuous polygon, not two separate polygons.
  - The default value for `max_iterations=100` is vastly increased,
    since many of these weird cases were caused by having coordinate
    ranges orders of magnitude higher than `alpha=0.1` and so 10
    iterations would not be enough to avoid these weird situations.
  - Also when
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
    does not find a suitable solution after `max_iterations` tries, it
    returns `NULL` instead of proceeding to process the inadequate/empty
    polygon coordinate data. Hopefully this will allow weird cases to be
    skipped rather than throwing an obscure error.

## multienrichjam 0.0.83.900

### Bug fixes

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - One of the more obscure bugs, caused by creating a point hull with
    only two points, where the two points apparently resided at exactly
    the wrong angle. Only apparently occurred when the two points were
    within 3 degrees of some “evil angle” - resolved when points were
    rotated more than 3 degrees from this “evil angle”. Typically, a
    point hull is not possible with two points, since it does not create
    a polygon, so the workaround was to add a third point (recycling the
    first point) adding amount of noise to create a polygon. Except for
    some reason even [`rnorm()`](https://rdrr.io/r/stats/Normal.html)
    was not adding *enough* or the right kind of noise. No random value
    should ever exactly, reproducible, reside on the line. There must be
    some rounding that takes place in an internal function. Nonetheless,
    the “fix” for this purpose was to add more than one dummy point.

- [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)

  - Another bizarre “bug” was caused by some unconfirmed R function that
    appears to change `options("warn")` to `options("warn"=2)`
  - and then does not change it back! Probably related to RMarkdown
    knitting, since it seems to occur when the knitting is interrupted.
    To be fair, warnings should be resolved in Jam packages, that’s
    true. But they should not impose an error.
  - Filed under: “Things that worked just two minutes ago, but now cause
    an error for no reason.”
  - Added workaround for random but painful issue when
    `options("warn"=2)` which forces warnings to errors. No idea what
    made that setting, but it caused
    [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
    to fail due to implicit conversion of `empty` to whatever datatype
    was defined in the input matrix. In future this function may change,
    but for now this change keeps it working.

## multienrichjam 0.0.82.900

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - now `geneHitIM` and `geneHitList` are equivalent as input, either or
    both can be provided. The incidence matrix with values c(0, 1) are
    stored as `geneIM`, and if there are any other numeric values those
    values are stored in `geneIMdirection`. It does not (yet) verify
    that all genes involved in enrichment are also present in the
    `geneIM` and `geneIMdirection` matrices.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - New argument `gene_annotations` allows display of the `mem$geneIM`
    (gene hits per enrichment) and `mem$geneIMdirection` (gene direction
    of change) as annotations alongside the gene axis of the heatmap. By
    default, when `mem$geneIMdirection` is available, both are shown,
    otherwise only `geneIM` is shown. It can now be hidden.
  - New argument `simple_anno_size` to control the size of heatmap
    annotations.
  - New argument `annotation_suffix` to add optional suffix to the gene
    annotation labels, helpful to indicate the unit. For example `im`
    data has default suffix `"hit"`, while `direction` uses `"dir"`.
  - When `mem$geneIMdirection` is present, it is now used by default
    during clustering. Values are multiplied by `mem$geneIM` and any
    `NA` values in `mem$geneIMdirection` are replaced with `1` in order
    to maintain the `mem$geneIM` value. They are considered
    `unknown direction` in that sense.

## multienrichjam 0.0.81.900

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  and
  [`jam_igraph_arrows()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph_arrows.md)

  - fixed bug when plotting directed edges with arrows, the arrow head
    width (`h.lwd`) was not properly expanded to the number of arrows,
    causing arrow heads to appear “twisted”.
  - changed default obscure option `sh.adj=1` which draws straight edges
    to the base of the arrow head, not to the arrow tip. It allows the
    arrow tip to be a point without needing to equal the edge line
    width.

## multienrichjam 0.0.80.900

### changes to existing functions

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - Silenced some default verbose output.

## multienrichjam 0.0.79.900

### changes to existing functions

- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)

  - argument `nodeSortBy=c("x", "y")` changed to
    `nodeSortBy=c("x", "-y")`, consistent with top-to-bottom sorting on
    the y-axis. This change indirectly affects
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    default output.

## multienrichjam 0.0.78.900

### changes to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - The bivariate color scale was too pale for lower significance
    P-values, so the colors were encoded to have higher color saturation
    at the low end. The intermediate colors are improved, see
    [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)
    changes below.

- [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)

  - The blue-yellow color blending by default in
    [`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
    was still producing grey, despite using the “LAB” (or “LUV”, “XYZ”)
    color models. This is usually a symptom of using RGB color space,
    blending “blue” with “red/green” (yellow) produces “red/green/blue”
    (grey), and not usually seen when using “LUV” which is a 360-degree
    hue radial color wheel. It is probably still converted to RGB before
    blending, then back to LUV.
  - The default argument was changed to `use_model="sRGB"` which
    produces a somewhat green color when blending blue and gold. The
    red-gold blending is improved as well, producing a more saturated
    orange color. End result: The intermediate directional colors are
    more recognizable as partly up (orange) or partly down (green).
  - TODO: This color legend needs x-axis labels, showing z-score values.

- [`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md)

  ``` R
  * finally displays the x-axis label and numerical units,
  by default `"z-score"`.
  ```

### bug fixes

- [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - error occurred when attributes in `list` form were entirely `NA`,
    one section did not use
    [`unlist()`](https://rdrr.io/r/base/unlist.html) properly. The bug
    was fixed.

## multienrichjam 0.0.77.900

### bug fixes

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  and
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  threw errors

  - using `pathway_column_split` when there was only one column (one
    enrichment), this error is corrected.
  - Further, passing `cluster_rows` as a `function` caused the resulting
    `Heatmap` object not to store the `obj` with the dendrogram/hclust,
    instead it stored the function. This error was also corrected.

## multienrichjam 0.0.76.900

### changes to existing functions

The theme of this update is “Customizing mark.groups labeling”. Useful
with
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
when supplying a `list` of node groups or communities via `mark.groups`,
with labels defined in `names(mark.groups)`. The labels are
automatically placed outside the mark polygon, and can be adjusted in
size with `mark.cex`, and position with `mark.x.nudge`, `mark.y.nudge`.

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - new arguments: `label.cex`, `label.x.nudge`, `label.y.nudge` to
    customize the label font size, and label placement, when `label` is
    provided.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  new arguments:

  - `mark.cex` passed to
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
    as `label.cex`
  - `mark.x.nudge` passed to
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
    as `label.x.nudge`
  - `mark.y.nudge` passed to
    [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
    as `label.y.nudge`

### major changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - New default creates the gene-pathway incidence matrix heatmap data
    first, which serves as the basis for all other plots.
  - Impact on
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
    now inherits this clustering defined using the gene-pathway
    incidence matrix and only when called by
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md).
    This change is the primary motivation for the update, so that the
    enrichment heatmap clustering is informed and driven by gene
    content, and no longer reflects only the enrichment P-values.

### other changes to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - default `point_size_min=1` changed to `point_size_min=2` so the
    smallest points are still clearly visible, including the fill color.
  - default `p_cutoff=1e-6` changed to `p_cutoff=1e-10`

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - default `p_cutoff=1e-6` changed to `p_cutoff=1e-10`

### bug fixes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - Edge case where `geneHitList` could be supplied as a `list` of
    signed values, a directional hit list as used in
    [`venndir::venndir()`](https://jmw86069.github.io/venndir/reference/venndir.html).
    It did not get recognized, and was converted using the `numeric`
    values rather than using the names of the values.
  - The preferred option is to supply `geneHitIM` for signed data,
    however it now works with a signed hit `list` as well.

## multienrichjam 0.0.75.900

### bug fixes

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - New arguments to handle forward-slash delimited gene symbols used by
    IPA to indicate two or more genes they consider to be one
    biomolecule for the purpose of pathway enrichment. The forward-slash
    “/” is also the delimiter user in `clusterProfiler` object
    `enrichResult` which causes these genes to break that compatibility.
    The new arguments represent a workaround to handle IPA data, so that
    downstream functions do not require changes.
  - `convert_ipa_slash=TRUE` enables the workaround, which converts
    forward-slash “/” to another delimiter.
  - `ipa_slash_sep=":"` defines the alternate delimiter to use, the
    default `":"` was chosed because it does not interfere with other
    common delimiters used in gene symbols, and does not cause problems
    with regular expressions, which would have been a risk with using
    `"|"`.

- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
  was throwing an error when `size2` was not already defined in `igraph`
  vertex attributes. When there is no `size2` it now calls
  `default_igraph_values()$vertex$size2` to use the appropriate default
  value.

  - other related errors are now caught and avoided, relating to steps
    that avoid updating empty entries in `list` vertex attributes, now
    it properly ignores pie attributes which were previously empty, and
    only updates entries with non-zero results.

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - Fixed longstanding errors when trying to split the heatmap by row or
    column into more pieces than the data will allow. The current
    workaround is almost complete, covers obvious cases where the
    requested number is higher than the number of columns or rows. It
    does not determine if there can only be one or two clusters but
    there are more columns or rows present. In future.

### other changes

- moved
  [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  to a separate R file, for future maintenance.

## multienrichjam 0.0.74.900

### updates to existing functions

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)

  - argument docs were updated
  - default colors in `col` were defined with more saturated default
    colors: `c("blue", "grey80", "firebrick3")`
  - argument default changed to `col_l_max=80` to accommodate higher
    `"grey80"` middle color.

- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)

  - new argument `pt.lwd=2` to control line width of open circles, used
    when `do_direction=TRUE`. The previous alternative was to define
    `par(lwd=2)` prior to calling
    [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md),
    then setting `par(lwd=1)` afterwards.
  - `directional_colors` use the same default colors used by
    [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md);
    and now includes `"no change"` as a specific legend entry.

## multienrichjam 0.0.73.900

### bug fixes

- Fixed errors rendering `igraph` nodes with `shape="circle"` only when
  called by
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
  and only with newer versions of `igraph` R package that expect the new
  attribute `vertex.frame.width` (not `vertex.frame.lwd` as I had
  hoped).

  - [`default_igraph_values()`](https://jmw86069.github.io/multienrichjam/reference/default_igraph_values.md)
    now defines `vertex.frame.width` and `vertex.frame.lwd`.
  - Note: `vertex.frame.lwd` is very likely to be removed from this
    package altogether, for compliance with `igraph`. However, I will do
    proper testing before making the change.

- Bug in
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
  which created a `data.frame`, and did not specify
  `stringsAsFactos=FALSE`, is fixed. Possibly similar bugs in other
  functions.

- Bug in
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  and
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  when enrichment only involves one gene, causing the `row_split` and
  `row_title` values to be incorrect. It is not a useful workflow, but
  the functions should handle the error edge case.

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
  now checks input argument `cnet` and `hitim` to make sure they are
  non-empty before processing. In theory it should not happen, but
  apparently it does when enrichment results are sparse or possibly
  empty.

- [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
  internal function uses proper `igraph::shapes(shape)` instead of
  previously internal `list` object `igraph:::.igraph.shapes[[shape]]`.

## multienrichjam 0.0.72.900

### new functions

- `subset_mem()`

  - convenience function to subset an entire `mem` object by sets
    (pathways) or genes. It will subset all internal incidence matrix
    objects for consistency.
  - This approach will become the preferred approach to display a
    specific set of pathways in a Cnet plot, by calling `subset_mem()`
    then
    [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).
  - [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
    currently subsets `mem` data internally, but for consistency may
    call `subset_mem()` instead. The arguments are designed to be very
    similar.

- [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)

  - simple helper function to retrieve or define layout coordinates for
    an `igraph` object. Because there is some logic to the process, it
    makes sense to put into its own function.

### changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
  [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)

  - recognizes new `igraph` attributes: `vertex.label.fontsize`,
    `edge.label.fontsize`.

    - These are not standard `igraph` attributes and will not be honored
      by default
      [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
      functions.
    - When `vertex.label.fontsize` is specified as a font size in
      points, this font size is used with no modification by
      `vertex.label.cex`.
    - When any `vertex.label.fontsize` value is `NA`, the default
      behavior is used to calculate font size, which uses
      `vertex.label.cex`. Therefore the `vertex.label.fontsize` can be
      defined for a single node, as long as all other node values are
      `NA`, and only the one node font size will be adjusted to this
      specific fontsize.
    - The `igraph` labels are drawn using
      [`text()`](https://rdrr.io/r/graphics/text.html), and the final
      exact font size is calculated for nodes:
      `par("ps") * par("cex") * vertex.label.cex` and for edges:
      `par("ps") * par("cex") * edge.label.cex`

  - new argument `label_fontsize_l` used to apply specific
    `vertex.label.fontsize` based upon node attribute values. For
    example `label_fontsize_l=list(nodeType=c(Gene=10, Set=14))` will
    define Gene nodes with fontsize 10, Set nodes with fontsize 14.

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)

  - when `frame_blank=NULL` is passed as an argument, it is interpreted
    as `frame_blank=NA` which uses no `frame.color` for blank nodes.

- [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - The following change was made to default behavior, so any method
    that reorders igraph nodes by color/border/name will be affected.
    Nodes will by default be sorted left-right when the nodes in a
    nodeset have less than 25% the x-span (width) compared to y-span
    (height), otherwise nodes will be sorted top-bottom. The direction
    along each axis respects the original argument `nodeSortBy`, to
    allow specific order based upon the data, or the natural order per
    the locale.
  - new arguments `orderByAspect` and `aspectThreshold` control when to
    sort left-right or top-bottom.
  - When the nodeset coordinate aspect ratio is taller (25% higher
    y-span than x-span) the nodes are sorted top-bottom, otherwise nodes
    are sorted left-right.

### bug fixes

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  threw an error for certain custom input that did not meet expected
  constraints. The function was updated to prevent these errors and to
  be more robust to this type of issue.

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  is called when rendering `shape="jampie"` nodes.

  - Previously when `frame.lwd=0` and `frame.color="black"` a small
    black border was drawn around the node. The new behavior when
    `frame.lwd=0` is to replace the color with `NA` so no outer border
    is drawn.

- [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - Fixed bug where `nodesets` argument was not always matched due to
    truncating the nodeset label to 25 characters. Now the nodeset label
    is not truncated.

## multienrichjam 0.0.71.900

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - After gene-pathway heatmap clustering, a collapsed Cnet `igraph` is
    created using pathway column clusters. The pathway cluster nodes are
    colorized based upon the proportion of each enrichment in that
    cluster, however it uses `mem$enrichIMcolors` by default. This
    function now applies statistical thresholds `p_cutoff` and
    `min_set_ct_each` prior to this step so the resulting colors will
    reflect those thresholds.
  - Help documentation was updated to include this information.

## multienrichjam 0.0.70.900

### bug fixes

- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)

  - Rare scenario using jampie nodes, when vertex.pie.lty is not
    defined, the default was not properly expanded to vcount, causing
    error `"subscript out of bounds"` when referencing
    `vertex.pie.lty[[i]]`.
  - Error above is caused by missing node attribute `"pie"`, which are
    now filled in with uniform values of 1 based upon
    `lengths(vertex.pie.color)`. This scenario usually occurs when
    trying to create pie nodes outside the “typical” scenarios, for
    example manually assigning attributes and not populating all the
    necessary values.

## multienrichjam 0.0.69.900

### new functions

- [`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md)

- [`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
  [`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md)

  - conversion functions that help interconvert between `igraph`
    `communities` objects, and `nodegroups` which is a `list` of node
    names.

- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

  - replacement for
    [`enrichMapJam()`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md)
  - converts `mem` output to an `igraph` object with multienrichment
    features, such as `pie` nodes with appropriate color fill,
    `pie.border` colored by direction when `enrichIMdirection` is
    defined.
  - by default, the resulting network has community detection called,
    then visualized with boundaries around the various nodes.

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - Updated to pass `raster_device` to
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
    to work around temporary error when the `"magick"` package is not
    available, `use_raster=TRUE`, which causes an error during
    rasterization. The error is resolved when changing from default
    raster device to `raster_device="agg_png"`, although this change
    requires the `"ragg"` R package is installed. So the change tests if
    `"ragg"` is available, and if so it passes
    `raster_device="agg_png"`. The change should not affect any other
    scenarios.

- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)

  - New argument `bundle_self=FALSE` changes previous default behavior
    by not bundling nodes that connect from and to the same nodegroups.
    Previously, nodes connecting within the same nodegroup would bundle
    through the center point of the cluster, which does minimize the
    busy edge lines, but makes it difficult to follow any paths.

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - new arguments `label`, `label_preset`, `label_adj_preset` are ussed
    to define optional label to appear outside the resulting hull. The
    label is intended to be used for network communities, to allow a
    label associated with each community when relevant.
  - Label placement is experimental and could change in future.
  - Labels are placed relative to the center of the layout, using the
    angle from layout center to hull center. Labels are placed outside
    the rectangular bounding box of the point hull, with text aligned to
    the outer edge based upon the nearest 45 degree angle from layout
    center to hull center.

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
  default `frame_lwd=0.2` changed from `frame_lwd=1`.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
  [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)

  - new argument `bundle_self=FALSE` passed to
    [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
    and ultimately to
    [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md).
    When TRUE any edges that connect from and to the same nodegroup will
    be bundled through the nodegroup center. The default FALSE does not
    bundle within nodegroup, so edges are only bundled when connecting
    two different nodegroups.

- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  now handles various combinations of missing `pie.lwd` and `frame.lwd`
  more gracefully.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  uses smaller node size by default.

## multienrichjam 0.0.68.910

It turns out that
[`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html) only
properly closes the polygon when it has a color fill with non-zero alpha
transparency. In this case, the final “closed” corner is correctly
extended to complete the sharp corner edge using line join
`par("ljoin")`. With no color fill, or with completely transparent color
fill, the final corner is not completed, the lines are “ended” using
`par("lend")`, and therefore there is no sharp corner. The workaround is
to apply a color with alpha transparency 1 (on scale of 0 to 255), which
causes the border to be drawn completely. However, some rare graphical
output devices do not support alpha transparency, so there is the chance
of rendering unintended opaque color fill. The situation should only
affect node shapes `"pie"` and `"jampie"`, which are designed to draw
the outer border then inner border, so the impact should be minimal.
Shape `"coloredrectangle"` actually calls
[`graphics::symbols()`](https://rdrr.io/r/graphics/symbols.html) which
avoids this issue. However, cases where outer border is expected to be
drawn after the inner border, it has small risk of rendering the fill
color fully opaque, covering the inner border. One day I probably need
to replace
[`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html) with a
custom function `polygon_with_borders()` that can handle inner and outer
border properly.

### updates to existing functions

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  which is used to render `igraph` nodes `shape="jampie"` was modified
  to use `col="#FFFFFF01"` for color fill of polygon outer borders.
- [`adjust_polygon_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md)
  examples were modified to show the effect of using `col=NA` to
  `col="#FFFFFF01"` on polygon border rendering.

### bug fixes

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  and
  [`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)
  were updated to force `stringsAsFactors=FALSE` when creating
  `data.frame` objects, fixing weird color glitch in the examples for
  [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md).

## multienrichjam 0.0.68.900

### Notable release notes

- `igraph` nodes with `shape="pie"` and `shape="jampie"`

  - New attributes `vertex.frame.lwd` and `vertex.pie.lwd` to customize
    the respective line widths of node borders. These attributes require
    using
    [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
    for `shape="pie"`, or require using `vertex.shape="jampie"` with
    `igraph::plot()`.

  - [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
    is now recommended as a more complete replacement for
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

  - `shape="jampie"` no longer uses `par("lwd")`, which was global
    change that affected all nodes, edges, and plot features.

  - Pie wedges use “inner borders” for each node, so adjacent wedge
    borders will not overlap. Note that inner border slightly overlaps
    the interior node fill color, to maintain consistent node sizes.

  - `frame.color` uses “outer borders” for each node, so these borders
    will not overlap inner pie wedge borders. Nodes are slightly
    adjusted smaller to maintain consistent node sizes, by default.

  - [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
    has slightly different logic for pie and frame colors, and now
    assigns `pie.lwd` and `frame.lwd`.

  - [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
    no longer changes single-color `shape="pie"` nodes to
    `shape="circle"` by default.

  - **Changes suggested**:

    - Prefer
      [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
      over [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
      [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).
    - Previous use of `par("lwd")` should be removed, and replaced with
      `vertex.pie.lwd` and `frame.lwd`.
    - Previously, single-color `shape="pie"` or `shape="jampie"` nodes
      were changed to `shape="circle"` by
      [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
      however they now stay `shape="pie"` in order to maintain control
      of line widths. The `igraph` nodes `shape="circle"` do not respect
      line width `lwd`.
    - Use of `shape="circle"` should be changed to `vertex.shape="pie"`
      or `vertex.shape="jampie"`, this change affects `nodeType="Gene"`
      moreso than `nodeType="Set"`.

  - Edges are now clipped to the outer border of nodes, which helps when
    using transparent node fill, or when using edge arrows.

    - With transparent nodes, the edge was previously shown connected to
      the center of each node (without clipping); similarly, edge arrows
      connected to the center of each node, effectively invisible and
      certainly not useful. However, cnet plots do not use edge arrows.
    - Most edges appear identical to previous rendering, however the
      control points used during bundling are used to determine where
      the edge connects to each node border. Thus, when edges converge
      on one node, bundled edges appear to “merge together” at one point
      connecting to the node. Previously nodes entered from multiple
      angles, since they connected to the interior center of the node,
      and showed somewhat more spacing based upon the number of edges.
    - Future options may include retaining some edge spacing based upon
      the intermediate edge curvature of each edge, instead of using the
      same control point for all edges, so the connection point on the
      central node will be slightly spaced out.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  and
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  changed some default sizes. The new defaults should work well without
  adjustment.

  - Node sizes are 4x larger than before, because when calling
    downstream functions we found ourselves always scaling nodes 4x
    larger.

    - The defaults: `categorySize=20`, `geneSize=10` are exactly 4x
      larger.
    - **Changes required**: Change previous `node_factor=4` to
      `node_factor=1` for reproducibility; or when calling
      [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
      use arguments `memIM2cnet(..., categorySize=5, geneSize=2.5)`.

  - node label sizes also have new defaults, similarly because we most
    often applied `label_factor=1.3` by default.

    - The defaults: `categoryCex=1.2`, `geneCex=0.9` are adjusted from
      previous `categoryCex=0.9`, `geneCex=0.7`.
    - **Changes required**: Change `label_factor` or `label_factor_l` to
      1, or change the call
      `memIM2cnet(..., categoryCex=0.9, geneCex=0.7)`.

- igraph defaults were changed when using
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  for plotting:

  - `vertex.label.family` changed from `"serif"` to `"sans"`.
  - `vertex.pie.border` changed to `"grey30"`, the `igraph` default
    values were identical to `vertex.pie.color` and therefore not
    visible.
  - `vertex.frame.lwd` set to 1, although this value is ignored by
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
    since it only uses `par("lwd")` for all lines.

### changes to dependencies

- `arules` was removed as a dependency, as two functions that used its
  class `"transactions"` were rewritten to remove the requirement. The
  previous version of those functions are internal and renamed from:
  [`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md)
  to `im2list_dep()`; and
  [`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md)
  to `imSigned2list_dep()`. The replacement functions should return
  identical data.
- `bezier` was added as dependency to generate edge bundling curves, in
  head-to-head tests, it out-performed
  [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html), and
  produced identical results to internal function
  `ggforce::bezierPath()`. Since `bezier` has no dependencies, this
  addition should feel small.
- `jamba` was bumped to version `0.0.88.900` for an important fix to
  [`mixedSort()`](https://jmw86069.github.io/jamba/reference/mixedSort.html).
- `colorjam` was bumped to version `0.0.23.900` for consistency.

### changes to existing functions

- [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)

  - When nodesets are not matched in the data, it prints a
    [`warning()`](https://rdrr.io/r/base/warning.html) then returns the
    input graph without change. Previously it called
    [`stop()`](https://rdrr.io/r/base/stop.html), which is problematic
    for relevant workflows.

- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)

  - default argument changed to `pie_to_circle=FALSE` so single-item
    nodes `shape="pie"` are no longer converted to `shape="circle"`,
    since `"jampie"` nodes are rendered much better. Haha.

- [`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md)

  - now accepts vectorized input and processes accordingly.
  - New examples show the varied features, many of which are not used
    for `igraph` node shape=“ellipse”, but they could be used in future.

- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md),
  [`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)

  - These functions draw custom `igraph` node shapes.

    - `shape="jampie"` draws a pie node shape with more features than
      vanilla igraph pie shape. Each pie wedge is color filled with
      `pie.color`, with new optional inner border `pie.border`,
      `pie.lwd`. The pie node overall can have an outer border defined
      by `frame.color` and `frame.lwd`. Also, pie nodes with only one
      color are rendered as a circle, with no tiny internal line.
    - `shape="coloredrectangle"` draws a series of square boxes with
      color fill `coloredrect.color`, and optional inner border
      `coloredrect.border`. The node overall can have an outer border
      defined by `frame.color`, and `frame.lwd`.

  - Both functions now use inner and outer borders via
    [`adjust_polygon_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md),
    so the various borders no longer overlap. Previously, pie nodes drew
    each wedge with default
    [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html),
    which allows adjacent borders to overlap 100%.

  - Note that nodes are resized internally so the rendered node size is
    equal across all nodes even when the line widths (lwd) vary.

  - There are some idiosyncracies from calling
    [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html) to
    render pie wedges, since it does not by default allow vectorized
    plotting of multiple polygons with different line widths. Therefore
    `shape="jampie"` renders line widths in subsets with identical line
    widths for vectorized plotting, which is 10-100x faster than
    plotting each pie node individually as done by
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).
    The main potential issue would be seen with partially overlapping
    nodes, with the potential to display inconsistent overlap order.

  - It seems hopeless to evaluate ggraph/tidygraph to render
    visualizations, in part due to visualization and rendering details.
    Also, ggraph/tidygraph does not store nor recognize many `igraph`
    visualization details in the `igraph` object, instead they must be
    encoded as `ggplot2` visualization options and settings. Using that
    ecosystem would also require creating new ggplot2 node geom types,
    and edge bundling functions.

- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)

  - This update represents a substantial refactor of logic.

  - Edge bundles can be “invalid”, which causes edges to be drawn as
    linear edges between nodes. The criteria were based upon a series of
    test cases which were all so common that they warranted being fixed.

    - Bundling usually occurs along the line between two nodegroup
      center points, calculated as mean node coordinates in each
      nodegroup. Sometimes the configuration of the center points, or
      the line itself, cause edge bundling to become unnecessary or
      ineffective.
    - Co-linear control points: When the edge and control points are
      co-linear (along the same line), the edge is drawn as a line. The
      criteria uses correlation above 0.99 for node and control points
      of each edge. This criteria affects indiviual edges, so nodes in
      the same nodegroup may be rendered differently based upon the
      specific position. The problem is clear when one control point
      appears beyond the path between two nodes. For non-linear edges,
      the path would curve around to the far side of the node. For
      linear edges, it appears as a line that extends beyond one node
      with optional arrow pointing backward. The solution effectively
      draws the same edge, except clips the edge at the first node
      boundary.
    - Both nodegroups contain only one node: with only one node in each
      group, there is nothing to “bundle”.
    - Both nodegroups have identical center points, within some small
      tolerance as a small percentage (0.5%) of the overall layout
      range.
    - When one nodegroup contains only one node, and the other nodegroup
      center point sits inside the node boundary, there is no bundling.
      Note this criteria is dependent upon node size during rendering.
      The problem is the edge spline control point is inside the node,
      so the spline would curve inside the node boundary, then point
      back out to the node border from the inside. Instead, the edge is
      drawn as a straight line to the node boundary from the outside.
      This situation occurs when one nodegroup fully surrounds a central
      node, so edges are drawn directly to the central node.

  - Edge bundling “midpoint” represents the position along the line
    between two nodegroups.

    - When one nodegroup contains only one node, this line is clipped to
      the node boundary, so the midpoint is defined beginning at the
      outer edge of the node, toward the center of the other nodegroup.
    - Note that when both nodegroups contain only one node, edge
      bundling is already invalid (see above).
    - Note than when the other nodegroup center sits inside the node
      boundary, the edge bundling is also invalid.
    - Therefore this situation only occurs when one nodegroup center is
      already outside the border of the single-node nodegroup.

  - Edges are properly clipped using the relevant `igraph` shape clip
    function. See
    [`igraph::shapes()`](https://r.igraph.org/reference/shapes.html),
    and `igraph::shapes("circle")$clip` for specific examples.

  - Edge labels are rendered along edges as follows: Linear edges are
    encoded with three coordinates: start, middle, end. Spline edges are
    encoded using default
    [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html)
    which returns 100 points by default. Edge labels are placed using
    the coordinate most distant from the start and end node. For linear
    edges, the middle coordinate is used, for splines a point very near
    the middle of the edge is used.

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

  - This function is intended as an enhanced drop-in replacement for
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).
    It was updated to fulfill previously un-implemented features, so
    fulfill the promise of being a replacement.

  - New arguments, formally passed to internal
    [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md):

    - `mark.groups`, `mark.shape`, `mark.col`, `mark.border`,
      `mark.expand`, `mark.lwd`, `mark.lty`, `mark.smooth` - arguments
      to enable and customize the rendering of nodes within clusters or
      groups.
    - Note the new options: `mark.lwd`, `mark.lty` for each group
      border; `mark.smooth` to control whether the group polygon is
      smoothed; `mark.alpha` to control alpha transparency of fill
      colors when not already defined.
    - `mark.expand` is now expected to be provided as a fraction of plot
      layout range, which is very close to default behavior in `igraph`
      since the default `rescale=TRUE` forced all layout ranges between
      `c(-1, 1)`.

  - `edge_bundling` new option `"default"` will try to detect the most
    appropriate bunding method, based upon whether `nodegroups`,
    `mark.groups` are defined, otherwise it chooses `"connections"`.

  - Edge labels are now rendered for straight edges and edge bundled
    edges.

  - Edge labels can now accept multiple `"edge.family"` values, in the
    unlikely event of multiple fonts on the same plot. This scenario
    will cause an error with
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

  - Edges are properly clipped based upon the `igraph` node shape.

  - Internal adjustments to `node_factor`, `edge_factor`,
    `node_factor_l`, `edge_factor_l`, `label_factor`, `label_factor_l`,
    `label_dist_factor`, and `label_dist_factor_l` were adjusted to be
    applied more consistently.

  - Undefined layout is now properly calculated dynamically and passed
    to the internal rendering function
    [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md),
    so the values for `xlim`,`ylim` are now properly calculated.

### new functions

- [`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md)
  to create Cnet plot `igraph` data for testing.

- [`adjust_polygon_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md)
  defines inner and outer borders for polygons.

  - Using inner borders allows adjacent polygons to have their borders
    visible beside each other, without overlap.
  - Using an outer border allows the display of a border around a
    collection of polygons without overlapping their inner borders.
  - Borders can be layered inside or outside existing borders.
  - There are extensive examples showing various combinations of
    borders.

- [`shape.jampie.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md),
  [`shape.coloredrectangle.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)

  - Invisible to the user, however they are called for igraph shapes
    `"jampie"` and `"coloredrectangle"`, respectively.
  - These functions now properly clip edges to the outer border of each
    node, including optional inner and outer borders.

- [`shape.ellipse.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.clip.md)

  - calculates the optional rotation and size of the ellipse and adjusts
    the edge endpoints accordingly.

- `parse_igraph_plot_params()`

  - reproduces an internal `igraph` function, which cannot be called by
    CRAN-approved R packages.

- [`default_igraph_values()`](https://jmw86069.github.io/multienrichjam/reference/default_igraph_values.md)

  - reproduces internal `igraph` package data in a function call, for
    CRAN compliance.

- [`jam_igraph_arrows()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph_arrows.md)

  - Mimics and extends internal `igraph:::igraph.Arrows()` for CRAN
    compliance.
  - It also optionally only renders edge arrows, useful when an edge
    bundling function renders edges itself.

- [`get_igraph_arrow_mode()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_arrow_mode.md)

  - Another mimic of internal `igraph:::i.get.arrow.mode()` for CRAN
    compliance.

## multienrichjam 0.0.67.900

- added packages to Suggests to support new functions for node layout,
  and creation of more “correct” alpha hull polygons around points.

  - `alphahull` - best implementation of alpha hull. However, it also
    requires `sp` package, a heavy install which not advised because it
    is being retired in 2023 in favor of `sf`. Slight risk that the
    `alphahull` package is removed from CRAN if it is not updated.

- Added to Depends

  - `sf` - lightweight replacement of `sp` that provides useful
    geometric functions. It is added primarily because it improves
    rendering pie node borders, which are resized by the exact line
    width determined at the time of plotting pie nodes.

### new functions

- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)

  - This function is in active development, and is not yet used in other
    functions, but will be used in the next version.
  - takes set of points, makes an alpha hull using
    [`alphahull::ashape()`](https://rdrr.io/pkg/alphahull/man/ashape.html),
    then expands using
    [`sf::st_buffer()`](https://r-spatial.github.io/sf/reference/geos_unary.html).
  - If `alphahull` is not available, it uses
    [`grDevices::chull()`](https://rdrr.io/r/grDevices/chull.html) which
    does not produce the “ideal” shape but has no additional R
    depedencies. It is also identical to output from `igraph` and `sf`
    packages for point hulls, so it has strong precedent.

### changes to existing functions

- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)

  - Refactored for greater speed, in test cases with 585 nodes, previous
    output took 1.28 seconds, the new output takes 0.025 seconds, 50x
    speed increase. This function is called in numerous places in this
    package, so this improvement will also positively affect all sorts
    of other functions.

- [`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)

  - The sort algorithm was improved in cases where the color order were
    ambiguous, but including names in the tiebreak, before using `H,C,L`
    values.
  - Note: This function is used by
    [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
    when `colorV` is not explicitly provided, by detecting the probable
    order of colors based upon order of colors in multi-color nodes.

- [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - The sort order was updated to be more efficient, and to use better
    logic when following `colorV` or when defining `colorV` ad hoc.
  - When sorting by `"color"` the order will be defined by `colorV`
    whenever colors are aligned by `colorV`, otherwise colors are
    generally sorted by hue `"H"` in HCL space.
  - Examples were added to show clearly the different options for
    ordering, including visual examples of the
    [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
    rendering of pie nodes.

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)

  - Visual examples are in
    [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
    since the
    [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
    function is internal to the igraph plotting scheme.

  - This function draws pie nodes in vectorized fashion, properly
    drawing each set of pie polygons per node, then the appropriate
    borders defined by `pie.border` and `frame.color`.

  - `pie.border` has the option to be drawn inside the border of the
    polygon, so that adjacent Venn wedges will have the entire outer
    border color visible without overlapping the adjacent Venn wedge. To
    enable, set: `options("inner_pie_border"=TRUE)`

  - `frame.color` also has the option to be drawn in a manner that does
    not overlap `pie.border`, it will be drawn just outside the
    `pie.border`.

  - In cases where `frame.color` is not drawn, the `pie.border` radius
    is adjusted to exactly the line width of the `frame.color` border,
    so nodes will always be exactly the same sizes with or without the
    `frame.border`.

  - In most cases there should either be `pie.border` *or*
    `frame.color`, however it is possible at some point that
    `frame.color` and `pie.border` will both need to be applied, and
    this function can handle it.

  - Note this process now draws three layers of polygons:

    1.  each `pie.color` wedge fill color and no border
    2.  each `pie.border` wedge border color with no fill
    3.  overall `frame.color` border color with no fill

## multienrichjam 0.0.66.900

### changes to existing functions

- changed to remove calls to
  [`matrixStats::rowMaxs()`](https://rdrr.io/pkg/matrixStats/man/rowRanges.html)
  and
  [`matrixStats::rowMins()`](https://rdrr.io/pkg/matrixStats/man/rowRanges.html)
  to use base functions instead. This change was due to several R
  crashes that appear to be bugs somewhere in the upstream packages,
  that also occurred on Mac OSX, and on linux, but in both scenarios
  involved R-3.6.1 which is not likely to gain support traction with
  other package authors. Understandable.

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  fixed potential bug:

  - When `row_split` is passed via `...` to the underlying
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
    it was not aligned to the order of `rownames(matrix)` to be
    displayed in the heatmap, therefore the rows were split in the wrong
    visual order.
  - `row_split` is now a formal argument, and when supplied as a vector,
    the `names(row_split)` are used to align `rownames(matrix)`
    appropriately.

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - new argument `nodesets` to define a subset of nodes for which the
    reordering will be applied, which may be helpful when nodes in a
    nodeset are horizontal or vertical. In future, this option may be
    applied based upon the aspect ratio of nodes in a nodeset, or so
    that the `nodeSortBy` can be defined as a `list` named by
    `nodesets`. It gets complicated.
  - more output when `verbose=TRUE`
  - minor added checks for layout, ensuring matrix input for layout will
    match `rownames(layout)` to `V(g)$name`, just in case the order is
    not identical.

- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)

  - reverted apparent regression which did not pass `...` to child
    function
    [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
    therefore the `colorV` color order was not properly used when called
    in this manner.

- `rotate_igraph_nodes()`

  - Change default argument from `center="origin"` to `center="median"`.
    Rare change to argument default, justified by the change being
    fairly benign. Also the new default is more consistent with
    expectations, that nodes would be “rotated in place”. Practical
    outcome is the same, nodes are rotated exactly as before, but the
    coordinate range is more likely to remain consistent, when input
    layout is not already centered at coordinates `c(0, 0)`.
  - new argument `verbose=FALSE`
  - minor added checks for layout, ensuring matrix input for layout will
    match `rownames(layout)` to `V(g)$name`, just in case the order is
    not identical.

- [`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md)

  - Change default argument from `center="origin"` to `center="median"`.
    Rare change to argument default, justified by the change being
    fairly benign. Also the new default is more consistent with
    expectations, that coordinates would be “rotated in place”.
    Practical outcome is the same, points are rotated exactly as before,
    but the coordinate range is more likely to remain consistent,
    especially when input layout is not already centered at coordinates
    `c(0, 0)`.

## multienrichjam 0.0.65.900

Numerous changes were made to functions in order to improve the overall
Cnet plot layout experience. New functions were developed offline that
focus on layout specific to Cnet plots, useful for bipartite graphs in
general.

### changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  and
  [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)

  - new argument `plot_grid=FALSE` which optionally plots a grey grid in
    the background, with units equal to “percentage” across the layout
    coordinate range. This option is intended to help when manually
    adjusting node and node_set positions with
    [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)
    and
    [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
    both of which take units `x,y` in the form of fraction of the
    overall layout dimensions.

- [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)

  - new argument `nodes_xy` is intended to help enter adjustments when
    many nodes need to be nudged. The `list` is named by node (e.g. by
    gene), and contains x,y coordinate adjustments. It is completely
    equivalent to entering `node`, `x`, `y` as three independent
    vectors, but may be easier to use.
  - argument default changed from `use_grep=TRUE` to `use_grep=FALSE`,
    because nudging a node `"A"` should not also nudge every node that
    contains an `"a"` or `"A"`. This rare change in default argument
    value seems more helpful to avoid erroneous moves by default.

- [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)

  - new argument `constrain` which takes a `character` vector of node
    names, then ensures those node coordinate positions are properly
    configured in `constraints` so they do not move during iterative
    layout.
  - new default behavior is to define `init` using the current graph
    layout stored in `igraph::graph_attr(g, "layout")`, instead of using
    a random circular initial layout.

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  internal function `inner_pie_border`

  - new experimental argument `inner_pie_border` is intended to draw the
    pie wedge border on the inside of the pie wedge shape, so it does
    not directly overlap the border of an adjacent pie wedge shape. I
    could not find a polygon function in R that can draw borders on the
    inside edge of the polygon, which is surprising considering the
    considerable effort to draw GIS world maps in R. Any adjacent
    borders are directly overwritten, with no option for borders to be
    displayed side-by-side along the polygon edge. Seems like an
    opportunity for someone.
  - nonetheless this capability is not yet implemented, but the
    framework is in place to be released soon.

- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)

  - new arguments `do_directional`, `directional_column`, and
    `directional_colors` are intended to add new directional circles
    indicating up- and down-regulation.

### new functions

- [`plot_layout_scale()`](https://jmw86069.github.io/multienrichjam/reference/plot_layout_scale.md)
  plots a grey grid background to an igraph plot indicating percentage
  units across the range of layout coordinates.

## multienrichjam 0.0.64.900

### changes to existing functions

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  argument `sortAttributes` now includes `"frame.color"` in order to
  include `"pie.border"` for pie shape nodes, and `"frame.color"` for
  nodes overall.

- [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
  and
  [`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)

  - new argument `emptyValue` to control how empty incidence entries
    should appear, either as zero `0`, or as `NA`. This update fixes
    rare issue where missing enrichment P-values were reported as zero
    `0` instead of `1` by default.
  - These functions appear in `venndir` package, however we do not want
    to make this package dependent upon `venndir` just yet, so they
    remain here for now. In `venndir` these functions are named:
    [`venndir::list2im_opt()`](https://jmw86069.github.io/venndir/reference/list2im_opt.html)
    and
    [`venndir::list2im_value()`](https://jmw86069.github.io/venndir/reference/list2im_value.html).

- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  now uses `sortAttributes=NULL` default, when it is NULL it uses
  defaults from
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).
  Previously there was inconsistent defaults between the two functions.

- [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  and
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)

  - new arguments `geneIMdirection`, `enrichIMdirection` are used to
    call new function
    [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
    to colorize node borders by default when data is available.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  now calls
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  and no longer edits colors itself, those steps are performed by
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md);
  no longer calls
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md).

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

  - now applies exact `colorV` colors to `geneIMcolors` without calling
    [`colorjam::matrix2heatColors()`](https://jmw86069.github.io/colorjam/reference/matrix2heatColors.html)
    since that process returned slightly darker colors by default, and
    for little benefit.
  - new argument `geneHitIM` intended to allow directional hit matrix to
    be supplied, thus enabling other features as described above,
    specifically colorized border on Cnet `igraph` plots.

### new functions

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)

  - colorizes node border, pie.border, coloredrect.border based upon
    directionality, using `geneIMdirection` and `enrichIMdirection` when
    available.

- [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  is step one in migrating function names away from camelCase, toward
  snake_case. Not really a new function, but a new function name.

- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  is a similar rename of
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  which mostly takes `mem` as input instead of `memIM` anyway.

## multienrichjam 0.0.63.900

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - The attribute `"caption"` is formatted more cleanly.
  - New returned attribute `"draw_caption"` which is a function that
    draws the caption in the bottom-right corner of the heatmap, mainly
    because this location is least likely to overlap other heatmap
    labels. The location and style can be customized as needed.

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)

  - argument `sortAttributes` was updated to include `"pie.border"` in
    the default sort order, for future when Cnet nodes also include
    border color with the direction of change.

## multienrichjam 0.0.62.900

### changes

- bumped dependency on `jamba` to 0.0.87.900 to pick up all the recent
  updates.

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)

  - new argument `cluster_rows` to control row clustering, specifically
    to allow no row clustering.
  - new argument `do_plot=TRUE` to honor `do_plot=FALSE` from
    [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)

  - Now generates its own caption that includes relevant clustering and
    filtering parameters used, helpful to reproduce the original result.
    Caption is returned as `attr(hm, "caption")`.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

  - passes `do_plot` to
    [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  - uses `caption` generated by
    [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
    and returns both the gene-pathway heatmap `gp_hm` and the caption
    `gp_hm_caption`. The caption is relevant because the parameters
    define the clusters represented in subsequent Cnet plots.

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

  - new argument `remove_blank_colnames=TRUE`, same as previous
    behavior, with new option to disable this behavior and keep columns
    with all values in `c(NA, "")`. We found `zScore` is sometimes
    reported as entirely `NA` but may be useful to keep, for consistency
    with other enrichment results.

- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
  no longer includes `...` when calling
  [`lines()`](https://rdrr.io/r/graphics/lines.html), which should
  silence a fair number of harmless but annoying warning messages. We
  don’t need a warning message for everything, tysm.

### bug fixes

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  somehow broke in a recent `igraph` update, apparently the return type
  is no longer coerced to `character` vector, so needs to be converted
  directly.

  - [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    and
    [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
    were impacted as well.
  - [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
    was modified to call
    [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
    then failing that call, will return the full `cnet_collapsed` graph
    without subset. Messages are printed for review.

- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
  initial work on allowing custom midpoints within node groups, which
  would allow defining a custom midpoint in x,y coordinates, which may
  or may not be between the two sets of nodes. Implementation is not in
  place yet, but in progress.

## multienrichjam 0.0.61.900

### changes to existing functions

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  was updated to improve the consistency of sorting by different node
  properties, specifically to allow sorting by node fill, and node
  border color(s).

- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  is the function used to render node `shape="pie"` by
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
  is an optimized, vectorized method to render nodes in one shot, rather
  than drawing each in a `for()` loop by default.

  - Drawing even 20 or more `pie` nodes is substantially faster using
    [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
    compared to default `igraph` function.

  - Speed was improved, the method of converting list of polygon
    coordinates to numeric vectors spaced by `NA` was much improved.

  - It can render `pie.border` for each pie wedge of each node with
    `shape="pie"` or `shape="jampie"`. Note that the outer line may be
    covered by subsequent `frame.color`. Attribute `pie.border` is
    expected to be a `list` where `lengths(pie.border)` are equal to
    `lengths(pie)`.

  - It can render `frame.color` around the full circle of each node with
    `shape="pie"` or `shape="jampie"`. This line may cover the
    `pie.border` if also drawn. Attribute `frame.color` is expected to
    have length equal to `igraph::vcount(g)` which is the total number
    of nodes, one `frame.color` value per node.

  - It is recommended to use one style or the other for each node:

    1.  `pie.border=NA`, and `frame.color="red"`
    2.  `pie.border=c("red", "gold")`, and `frame.color=NA`

- Internal function
  [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  is called by
  [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md),
  and was updated to handle `frame.color`.

## multienrichjam 0.0.60.900

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  default arguments changed:

  - `min_set_ct=1`, previously was 2
  - `min_gene_ct=1`, previously was 2

### minor bug fixes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  was not verifying `names(geneHitList)` also matched
  `names(enrichList)`, therefore a mismatch could cause downstream
  errors in
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  related functions. This issue now only uses `geneHitList` with names
  that match.

## multienrichjam 0.0.59.900

- added `amap` package as dependency, it provides
  [`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html).

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  verbosity was reduced when `verbose=FALSE`.

### bug fixes, enhancements

- `mem_gene_pathway_heatmap()` now honors `p_floor` when defining the
  incidence matrix values to be used in clustering the weighted and
  combined enrichment and heatmap gene-pathway incidence matrix data.

- `mem_gene_pathway_heatmap()` was throwing an error when supplying
  `column_split` and the default `cluster_columns=TRUE`.

  - The error was caused by creating a dendrogram for `cluster_columns`
    then supplying
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
    with a `character` vector `column_split`, and a dendrogram/hclust
    for `cluster_columns`. Instead it allows passing a `function` to
    `cluster_columns`, which also requires using custom data, since the
    data used for clustering is a weighted combination of the enrichment
    P-values across the top of the heatmap, and the data inside the
    heatmap.
  - The new default when `row_split` or `column_split` are `character`
    will be to define `cluster_rows` or `cluster_columns`, respectively,
    to a `function` that calls
    [`amap::hcluster()`](https://rdrr.io/pkg/amap/man/hcluster.html) on
    the combined and weighted heatmap and respective annotation data
    matrices.

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  was not properly defining gene row cluster names for
  `mem_gene_pathway_heatmap()`, leaving them empty by default instead of
  assigning from `letters`.

## multienrichjam 0.0.58.900

- bumped dependency to `jamba (>= 0.0.84.900)` to retire
  [`call_fn_ellipsis()`](https://jmw86069.github.io/jamba/reference/call_fn_ellipsis.html)

### changes to existing functions

- [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  were updated:

  - `nameColname` is handled properly, without relying upon `"Name"`
    colname, and without relying upon
    [`rownames()`](https://rdrr.io/r/base/colnames.html) of the
    enrichment `data.frame`.
  - Now the subset operations use values in the `nameColname`.
  - Also the rows in each enrichment `data.frame` with values in
    `nameColname` are subset using
    [`subset()`](https://rdrr.io/r/base/subset.html), which means in
    some rare cases multiple rows might be returned, if the input
    enrichment data has the same name in `nameColname` for multiple
    rows. This change is intentional, in order to retain all rows with
    matching names, not just the first that occurs.

- [`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md):

  - argument `geneDelim` has a default value, consistent with
    [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

### bug fixes

- [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
  will return data with zero rows without error, although it usually
  occurs because of an error somewhere else (for example input data).

### functions removed

- [`call_fn_ellipsis()`](https://jmw86069.github.io/jamba/reference/call_fn_ellipsis.html)
  was moved to the jamba package version `0.0.84.900`.

## multienrichjam 0.0.57.900

- Issue [\#6](http://github.com/jmw86069/multienrichjam/issues/6)
  reported an error when using a series of enrichments where only a
  subset contain a z-score column name for `directionColname`. Related,
  the current approach ignored the `directionColname` when all values
  were NA. Both situations have been corrected, to allow flexible
  mish-mash of NA and non-NA values, and presence/absence of
  `directionColname` in each enrichment input. Also related,
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  was by default no enabling the directionality via
  `mem_enrichment_heatmap(mem, apply_direction=TRUE)`, therefore there
  is a new argument to
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  default `apply_direction=NULL` which will auto-detect whether there is
  directional non-zero and non-NA values that can be used in the
  heatmap. Also, `apply_direction` can be defined on its own to force
  the issue.

## multienrichjam 0.0.56.900

### bug fixes

- Issue [\#7](http://github.com/jmw86069/multienrichjam/issues/7)
  reported an error, traced back to the vignette. The error was caused
  by passing arguments to
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  by overloading `...`, when not all arguments were valid in
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
  The fix:

  - new function
    [`call_fn_ellipsis()`](https://jmw86069.github.io/jamba/reference/call_fn_ellipsis.html)
    which passes arguments including `...` to another function, and when
    that function arguments do not allow `...` then it limits the
    arguments in `...` to those arguments accepted by the other
    function.
  - Instead of: `x <- some_function(a=1, b=2, ...)`
  - Use: `x <- call_fn_ellipsis(some_function, a=1, b=2, ...)`
  - [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
    and `mem_gene_pathway_heatmap()` were updated to use
    [`call_fn_ellipsis()`](https://jmw86069.github.io/jamba/reference/call_fn_ellipsis.html).

- Another error was noted during the vignette workflow, that the
  `directionColname` was being populated even when no enrichment data
  contained non-NA values, which was inconsistent with
  [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  in function `enrich2IM()`. This situation was corrected by requiring
  only one enrichment result to contain a non-NA value in this column.

- [`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md),
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
  [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
  were updated to call
  [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  with the list of `enrichResult` objects.

### changes to existing functions

- [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  can accept a `list` object, which is expected to contain a list of
  `data.frame` and/or `enrichResult` objects. An `enrichResult` is
  converted to `data.frame` by `enrichResult@result`, all other objects
  must contain `colnames(x)`. When `require_non_na=TRUE` it will test
  each object, and return `max` unique entries that match `pattern`. The
  entries should all contain the same matching colnames for most
  purposes in `multienrichjam`.

## multienrichjam 0.0.55.900

### bug fixes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  was incorrectly populating empty gene counts with default `1` instead
  of `0`. The effect is mainly during filtering by gene count, where the
  minimum is usually never below `1`, however it can cause issues in
  point sizing particularly in
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md).
  The previous default `1` is a remnant of using this function to
  generate a matrix of enrichment P-values.

### changes to existing functions

- [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
  argument default was changed to `emptyValue=NA` so in the absence of
  data to populate into the incidence matrix, the default cell value
  will be `NA` to indicate there is no available data.
- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  was updated to supply a specific `emptyValue` for all calls to
  [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md).
- [`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md)
  new argument `type` is intended to allow re-using this same function
  for univariate color functions, so heatmaps features can be made
  consistently, specifically for dotplot or normal heatmap output, and
  optionally labeling cells with statistical values.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  is updated to share common heatmap code for bivariate and univariate
  color gradients. One by-product is that output cannot be raster
  format, which is typically not an issue for pathway enrichment, since
  pathways should not represent more than 1,000 or so pathways. In that
  event, output should probably be rasterized (PNG, JPG) instead of
  vector graphics (PDF, SVG).
- [`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md)
  new argument `digits` used to prevent displaying weird labels like
  `5.9999999998` and instead will display `6`.
- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  default argument was changed to `edge_bundling="connections"` which
  will enable edge bundling by default. It can be disabled with
  `edge_bundling="none"`, although it is only active when there are edge
  connections that can be bundled.

## multienrichjam 0.0.54.900

### updates to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  was updated to customize directional information in heatmaps. Cell
  labels are optionally displayed, defined by `show`, to display one or
  more of `z-score`, `-log10pvalue`, and `gene count`. Argument `sets`
  was updated to handle presence of `enrichIM` and `enrichIMgeneCount`
  and `enrichIMdirection` - as well as future measurements with prefix
  `enrichIM`. New argument `min_count` applies the gene count filter to
  the dot plot output heatmap.

## multienrichjam 0.0.53.900

This update mainly focuses on implementing directionality when available
in the pathway enrichment data during
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

### new functions

- [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md) -
  implements bivariate color scale, in this case for enrichment
  `-log10pvalue` for color intensity, and `z-score` for color hue
  directionality: blue “inhibited”, gold “neutral”, and red “activated”.
- [`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md) -
  display the color ramp defined by
  [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)
- [`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md) -
  define a heatmap cell function that uses the color defined by
  [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
  optionally plots circle points, and optional text labels
- [`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md) -
  creates a 2-D color ramp legend suitable for
  `ComplexHeatmap::draw(..., annotation_legend_list=x)`

### changes to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  new argument `apply_direction=TRUE` will enable bivariate colors,
  where color intensity is defined by the enrichment -log10 P-value, and
  color hue is defined by z-score direction:

  - blue = “inhibition” with z-score \<= -2
  - gold = “no direction” with z-score \> -2 and x score \< 2
  - red = “activation” with z-score \>= 2

Points are sized by number of genes.

## multienrichjam 0.0.52.900

### changes to existing functions

Moved some functions into their own .R file for better organization:

- [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  new arguments `directionColname` to define an optional column name
  that contains `numeric` values indicating direction of pathway
  enrichment. These values are often in the form of an
  `"Activation z-score"`, as is the case with IPA “Upstream Regulators”.
  Argument `direction_cutoff` refers to the absolute value required for
  direction be given a “sign” up or down.

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  new argument `directionColname` used to determine directionality of
  pathway enrichment, useful for things like `"Activation z-score"`, as
  returned by IPA “Upstream Regulators”. Output included
  `mem$enrichIMdirection` which contains the `numeric` values, where
  `NA` values are substituted with `0` zero. These values will be used
  in near future, likely in
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  to indicate predicted direction of impact on particular pathways.
  Argument `direction_cutoff` is passed to
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  for optional filtering to require at least one pathway to contain an
  absolute direction score at or above this threshold. The IPA z-score
  recommends a threshold z=score \>= 2 for “activation” or “inhibition”.
  Note that many pathways have no z-score, so applying this threshold
  will remove those pathways from downstream analysis.

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  new argument `apply_direction` and `direction_cutoff` determine
  whether to indicate direction if it exists, and optionally applies a
  threshold. Still in testing currently. Note that applying this cutoff
  will hide pathways whose `numeric` direction is below the threshold.
  New argument `gene_count_max` to apply a max gene count threshold for
  the point size for the dot plot format. New argument `legend_height`
  to control the heatmap legend color bar height.

## multienrichjam 0.0.51.900

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  will now hide all but the main enrichment colors in color legends when
  there are more than 8 combinations, to prevent the legend from taking
  the entire plot device space and not displaying the heatmap. The
  threshold is configurable with `show_heatmap_legend=8`. New arguments
  `show_gene_legend`, `show_pathway_legend` are `logical` and are
  intended to allow hiding the other color legends. A minimalist style
  is to show only the main enrichment colors, which are used for all
  other colors anyway.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  was updated so the default dotplot format has properly controlled
  point sizes, and legend point sizes. Previously the two sizes were
  independent and required manual adjustment. The current approach
  ensures the legend point size exactly match the heatmap dot plot point
  size. An optional parameter `cexCellnote` will display labels unless
  `cexCellnote=0` in which case gene count labels are hidden.

## multienrichjam 0.0.50.900

### changes to existing functions

- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  new argument `inset` passed to
  [`legend()`](https://rdrr.io/r/graphics/legend.html). Also this
  function uses [`tryCatch()`](https://rdrr.io/r/base/conditions.html)
  to try to pass `...` arguments, and if they fail it tries again
  without `...`. Fun.

## multienrichjam 0.0.49.900

### extended function help text

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  is a drop-in replacement for
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
  and its options were described in much more detail. In brief:

- it plots nodes vectorized which is substantially faster

- it changes the default to `rescale=FALSE` and maintains aspect ratio
  1:1, so layout coordinates are rendered without distorting the x- and
  y-axis ranges

- it optionally bundles edge connections, which helps particularly with
  large bipartite graphs (especially gene-to-pathway graphs.)

- it allows bulk adjustment to node size, label size, and label distance
  from node center.

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new argument `edge_bundling="connections"` that by default will bundle
  edges appropriate for Cnet plots. Disable with `edge_bundling="none"`.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  new argument `style="dotplot"` will create a dotplot styled heatmap,
  whose points are sized proportional to the number of genes involved in
  enrichment. This style is under development and may require additional
  customization options such as setting a max gene count for point
  sizes.

### bug fixes

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  fixed bug where enrichment names were mangled by `data.frame(...)`,
  avoided by using `data.frame(check.names=FALSE, ...)`.

Several functions were updated to include proper package prefix for
function calls:

- [`matrixStats::colMins()`](https://rdrr.io/pkg/matrixStats/man/rowRanges.html)
- [`igraph::V()`](https://r.igraph.org/reference/V.html)
- [`igraph::E()`](https://r.igraph.org/reference/E.html)
- [`igraph::degree()`](https://r.igraph.org/reference/degree.html)
- [`igraph::components()`](https://r.igraph.org/reference/components.html)
- [`igraph::vcount()`](https://r.igraph.org/reference/gorder.html)
- [`igraph::neighbors()`](https://r.igraph.org/reference/neighbors.html)
- [`igraph::vertex_attr_names()`](https://r.igraph.org/reference/vertex_attr_names.html)
- [`igraph::list.graph.attributes()`](https://r.igraph.org/reference/list.graph.attributes.html)
- [`igraph::ego()`](https://r.igraph.org/reference/ego.html)
- [`igraph::set_graph_attr()`](https://r.igraph.org/reference/set_graph_attr.html)
- [`igraph::graph_attr()`](https://r.igraph.org/reference/graph_attr.html)

affecting functions:

- [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
- [`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md)
- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
- [`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md)
- [`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)
- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
- [`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

## multienrichjam 0.0.48.900

### bug fixes

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  was updated to handle edge cases:

  - `gene_im_weight` 0 or 1; `enrich_im_weight` 0 or 1
  - user-supplied `cluster_rows` or `cluster_columns` with
    `colorize_by_gene`
  - The `colorize_by_gene` logic was changed to convert to an integer
    matrix that refers to colors and labels by factor levels, which
    helps for user-defined clustering methods, also helps labeling the
    color legend.
  - new argument `colramp` for user-defined color gradient.
  - `row_title` is no longer assigned when user-defined argument is
    `NULL`. This workaround helps the edge case with `gene_im_weight=1`
    where the `row_split` is forced to have a limited number of values
    regardless what `row_split` integer value is sent, thus causing
    `row_title` mismatch in length.
  - Fixed egregious type in second part of an `if` statement, only
    called when user supplies a custom cluster function.

## multienrichjam 0.0.47.900

### bug fixes

- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  was not correctly handling multiple `sourceColnames` values, instead
  was only using the first value. This bug has been corrected.

### changes to existing functions

- - [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
    and
    [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
    now return `enrichResult` when supplied with `enrichResult`, instead
    of coercing to `data.frame` which would then need to be converted
    back to `enrichResult`.
- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  argument default `sourceColnames` was changed to reflect usage with R
  package `msigdbr` for MSigDB gene set data, specifically
  `sourceColnames=c("gs_cat", "gs_subcat")`. Similarly, default values
  were removed from `curateFrom` and `curateTo`, since these defaults
  imposed a specific outcome.
- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  default values for `topEnrichCurate*` arguments were changed to NULL;
  default argument values changed to
  `topEnrichSources=c("gs_cat", "gs_subat")`.

## multienrichjam 0.0.46.900

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  new argument `rotate_heatmap=TRUE` will rotate the heatmap layout so
  pathway names are displayed as rows, and genes are displayed as
  columns.

I still need example data to use for function document examples, to show
utility of rotating the gene-pathway heatmap.

## multienrichjam 0.0.45.900

### changes to existing functions

- [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  was implemented twice (smh) so the older function was renamed
  `color_edges_by_nodes_deprecated()`. The older function blended colors
  using a simpler approach with
  [`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
  that took a very fast hue average; while the new function uses
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html)
  that uses red-yellow-blue additive color blending model. Honestly the
  function
  [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  is not used much yet, but is likely to be used more with community
  detection, and Cnet/community edge bundling. The new
  [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  also applies alpha to the intermediate blended colors, to retain the
  relative weight of colors during the blending step.

## multienrichjam 0.0.44.900

The edge bundling update!

The new functions are under rapid development, and represent potentially
quite useful functions for visualizing complex `igraph` networks,
specifically Cnet plots with Gene and Set node types.

The general technique of bundling edges between two groups of nodes is
relatively stable. The details and extent that edges are bundled,
including visual display of bundled edges, is under active evaluation
and development. So far it doesn’t seem to take much to improve the
output figure, compared to using straight edges.

In the near future, Cnet plots may by default enable some form of edge
bundling, as Cnet plots were in fact the motivating example.

See pkgdown docs for
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
for visual examples using the Karate network.

### current issues in edge bundling

[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
is the core function, and it does not yet implement:

- edge clipping based upon igraph vertex shape, for example using
  `igraph:::.igraph.shapes[["circle"]]$clip()`.
- edge arrows, which requires edge clipping for proper usage, otherwise
  edge arrows will be underneath the node shape itself.
- edge labels, which could be implemented, except that bundled edges are
  by definition much more likely to be too close for effective labels. I
  rarely use edge labels, so I will target the simplest thing that
  works.

### new igraph edge bundling functions

- [`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md) -
  is a general use function to return nodes that are grouped by having
  identical node neighbors. This situation usually happens rarely,
  except with bipartite graphs where there are often clusters of nodes
  with the same neighbors, especially common for Cnet plot data.
- [`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md)
  is the first edge bundling function, it is actually a light wrapper
  around
  [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md).
  Bipartite bundling connects a nodeset (defined as having the same
  neighbor nodes) to each neighbor node.
- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
  is a general edge bundling technique that bundles edges between node
  groups. These node groups can be defined by any relevant technique,
  typically a community detection algorithm such as
  [`igraph::cluster_walktrap()`](https://r.igraph.org/reference/cluster_walktrap.html)
  or any other of several `igraph::cluster_*` functions.

### new igraph functions

- [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  is a simple function that blends the colors of the two nodes involved
  in each edge, using
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html)
  since it uses RYB color blending. It is optimized so it only blends
  unique combinations of colors.

### new utility functions

- [`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md)
  is a necessary utility function that “expands” a `data.frame` that has
  one or more columns with multiple delimited values. It simply expands
  the rows to represent one individual value per row for those columns.
  This function will very likely be moved into the `"jamba"` package for
  much broader use.
- [`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md)
  is a helper function for
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  which is intended to update node and edge attributes in bulk based
  upon another attribute. For example “make all the Gene nodes small,
  and all the Set nodes large.” Same with labels, colors, etc.

### updates to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  has new arguments:

- `edge_bundling` - for `"nodegroups"`, `"connections"` and `"none"`
  This argument will enable edge bundling by calling
  [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
  and
  [`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md)
  when needed. See examples.

- `render_nodes,render_edges,render_nodelabels,render_nodegroups` are
  options to skip various aspects of `igraph` plot features, either to
  save time, or to allow more detailed visual layering.

- `vectorized_node_shapes` toggle the vectorized node shape plotting
  when there is more than one node shape in the `igraph` object. Note
  this feature is *substantially* faster, but changes the order that
  nodes are rendered, by design. As a result, nodes are drawn in order
  of shape, so each render is a bulk operation. For nodes that overlap,
  or partially overlap, this feature will visibly change the ordering of
  nodes. In general, speed is still ideal, and reducing node layout
  overlaps should be a separate step.

## multienrichjam 0.0.43.900

### changes to existing functions

- [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  logic was updated for `expand` which was not properly applying
  negative values to compress the spacing between nodes in a node set.
  This use is relatively rare but feasible.
- [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  argument `set_nodes` now optionally allows referring to a nodeset
  using the name of one node contained in the nodeset, as a convenience.

### new functions

- [`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md)
  attempts to automate the process of applying
  [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  for each subcluster of nodes in a Cnet plot to enforce a minimum
  spacing between nodes. The Fruchterman-Reingold layout algorithm might
  offer a minimum distance threshold, but I could not find it.
  Incidentally, this new function can also compress node spacing, which
  might be useful when gene labels are not shown.

## multienrichjam 0.0.42.900

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  has much more detailed help text, describing more details about
  clustering parameters, and describing the returned objects in detail.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new argument `do_plot=TRUE` determines whether each plot is rendered,
  or just returned as a plot object to be reviewed separately.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  the returned data now includes the gene-pathway heatmap caption, which
  descibes the gene/set filtering criteria, and the row/column distance
  methods used for clustering. Since these arguments have substantial
  effect on the pathway clusters, it is helpful to keep this information
  readily available. Similarly, the `cnet` plots as `igraph` objects
  store the title as a graph attribute, accessible using
  `graph_attr(cnet, "title")`.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  argument `colorize_by_gene=TRUE` now uses
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html).
- `mem_gene_pathway_heatmap()` is tolerant of `'...'` entries not valid
  with
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  in order to allow overloading `'...'` for other function arguments.

### new functions

- [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)
  is intended to help move individual nodes in an `igraph` layout,
  useful for adjusting nodes to reduce label overlaps.
- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  is a helper function to get the nodes in a nodeset, defined as the
  `"Gene"` nodes that all connect to the same `"Set"` nodes. It is
  helpful when looking at a Cnet plot, and wanting easy access to the
  `"Gene"` nodes in a particular cluster of nodes in the plot.
- [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  is intended to manipulate the layout coordinates for nodes in a
  nodeset, useful to expand, shift, or rotate nodes to help with
  visibility.
- [`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md)
  is a very useful function with complicated name. It is intended to
  help move `Set` nodes in a Cnet `igraph` layout, then re-positions all
  the `Gene` nodes while keeping the `Set` nodes in fixed positions.
- [`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)
  which allows rotating an `igraph` layout. It can also reflect layout
  coordinates across one or more axes. It also calls other helper
  functions as needed,
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  and
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).
- [`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md)
  is the underlying function to rotate coordinates, it operates on a
  numeric matrix and is called by
  [`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md).

The next version should include `plot_cnet_heatmaps()` which creates
nice Cnet cluster plots where each cluster has a corresponding
expression heatmap displayed at the edge of the figure.

## multienrichjam 0.0.41.900

### bug fixes

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  returns `NULL` when there is an error during the Cnet collapse step,
  and does not report the error. Now it at least reports the error. The
  underlying cause in this case was ComplexHeatmap not accepting
  horizontal color legend orientation in R-3.6.1 (grid package version
  3.6.1) because of some new unit arithmetic only available in grid
  4.0.0+. Hiding the color legend, or using vertical color legend fixed
  the issue.

## multienrichjam 0.0.40.900

### bug fixes

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md) -
  fixed error when using [`apply()`](https://rdrr.io/r/base/apply.html)
  on a matrix that had one row or one column.
- [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
  fixed issue when the Cnet igraph object has no internal layout, it
  runs
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)
  by default.

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  now returns `clusters_mem` in the `list`, which includes the pathway
  set names represented in each cluster shown in the gene-pathway
  incidence matrix heatmap.

## multienrichjam 0.0.39.900

### changes to existing functions

- `cnetplot_internalJam()` was updated to include `"nodeType"` as a node
  attribute, which helps distinguish `"Gene"` and `"Set"` nodes in
  downstream operations. This change helps address
  [\#5](http://github.com/jmw86069/multienrichjam/issues/5) to hide the
  gene labels on Cnet plots.

## multienrichjam 0.0.38.900

### bug fixes

- Update to address issue
  [\#4](http://github.com/jmw86069/multienrichjam/issues/4) in
  [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  with new argument `descriptionColname` which will force the resulting
  colname to be `"Description"` to fit expectations of
  `enrichplot:::fortify.internal()` which requires `df$Description`.
  When `descriptionColname` is not supplied, or not found in the input
  `enrichDF` a warning is issued that describes the problem.

### enhancements

- [`heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
  and
  [`heatmap_column_order()`](https://jmw86069.github.io/jamba/reference/heatmap_column_order.html)
  now also work with `HeatmapList` objects. By default they use the
  first heatmap in the list, which should be consistent with all other
  heatmaps.

## multienrichjam 0.0.37.900

### bug fixed

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  argument `subsetSets` was never implemented; existing arguments
  `descriptionGrep` and `nameGrep` were used for similar but
  insufficient purpose. The `subsetSets` argument defines specific
  pathway names to retain for analysis, and is implemented through
  [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and by proxy
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md).
  These functions might better be called `subsetEnrichResult()` and
  `subsetEnrichList()`, but will not rename these functions.

### enhancements

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  includes the number of rows (genes) and columns (pathways) displayed
  in the gene-pathway incidence matrix heatmap.

## multienrichjam 0.0.36.900

### bug fixed

- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  was updated to fix
  [\#3](http://github.com/jmw86069/multienrichjam/issues/3), thanks to
  [@john-lee-johnson](http://github.com/john-lee-johnson) for reporting.
  Issue arose because
  [`enrichplot::cnetplot()`](https://rdrr.io/pkg/ggtangle/man/cnetplot.html)
  expected enrichResult rownames to be equal to values in the `"ID"`
  column of the enrichment `data.frame`.

### changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  new argument `min_count` requires a pathway to contain at least
  `min_count` genes in order to be considered a “hit”. This filter is
  mostly important when used with `topEnrichN` to use the top pathways –
  it will only sort pathways then take the top `topEnrichN` number of
  pathways that also contain at least `min_count` genes.
- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  new arguments `min_count`, `p_cutoff` require pathways to contain at
  least `min_count` genes, and have no higher than `p_cutoff` enrichment
  P-Value. Previously these functions only took the top pathways,
  regardless of these filters (which were applied later). This change
  allows the filters to be applied before taking the top `topEnrichN`
  pathways, which mostly helps when `min_count` is greater than 1 –
  sometimes pathways with only one gene involved in enrichment have
  statistically significant P-value, but are not biologically relevant
  for interpretation or follow-up experiments. For
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  the filter is applied to each enrichment in the list, and any pathways
  meeting the criteria are taken for all enrichment lists. So a pathway
  must be present in the top `topEnrichN` entries which meet both the
  `p_cutoff` and `min_count` criteria to be retained by these functions.

## multienrichjam 0.0.35.900

### bug fixes

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  was updated to fix a small issue with `min_set_ct_each` which requires
  at least one enrichment to have `min_set_ct_each` genes. However, this
  filter was not applied alongside the `p_cutoff` – therefore some
  pathways with enough genes which were not significantly enriched were
  fulfilling these criteria, as along as another enrichment was
  significant. The new behavior (as expected) is to requires an pathway
  to meet both the `min_set_ct_each` and `p_cutoff` thresholds in the
  same enrichment in order to be retained in the gene-pathway incidence
  matrix.
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  fixed edge case where genes present in multiple enrichments counted
  more toward `min_set_ct_each`, but should only count once per gene.

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  new arguments `column_title`, and `row_title` allow custom cluster
  names, which are also carried to Cnet cluster names as needed.
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  was updated to handle `row_split` and `column_split` arguments more
  intuitively, and to allow `row_split` values `FALSE`,`0`,`1` to
  inactivate the split completely. Previously it was not possible to
  turn off split.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  argument was renamed from `pathway_row_split` to `gene_row_split` to
  reflect the intent of this argument more accurately. This change is
  early in the function lifecycle and not in broader use yet – better to
  change it now. Otherwise, future argument name changes will not occur
  without some type of backward compatibility.
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
  [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  and
  [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new arguments `row_cex` and `column_cex` used to adjust row and column
  heatmap labels, which is helpful when used with auto-sized axis labels
  to make minor adjustments.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  was updated to make the color ramp more consistent with
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md),
  so the color scale atop the gene-pathway heatmap more accurately
  reflects the color scale used in the enrichment heatmap itself.

## multienrichjam 0.0.34.900

### changes to existing functions

- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  new argument `colorize_by_gene=TRUE` will color the heatmap body using
  blended colors from the `geneIMcolors` which represents the
  enrichments in which the gene is involved. The default
  `colorize_by_gene=FALSE` instead colors the heatmap body by the number
  of enrichments, which can be confusing if one of the enrichment colors
  is also red. The goal is to make it visually apparent when a gene is
  involved in one enrichment, by using the color from that enrichment.
  This feature is still in development and testing.
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  changed argument name from `enrich_gene_weight` to `enrich_im_weight`,
  before this argument is in wider use. Added new argument
  `gene_im_weight`. These arguments more accurately reflect the relative
  weight between enrichment and incidence matrix for `enrich_im_weight`.

## multienrichjam 0.0.33.900

### changes to existing functions

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  was updated with new argument `cluster_color_min_fraction` to help
  filter the enrichment colors to include in each resulting cnet
  cluster. The intent is not to represent colors where the number of
  significant pathways is below this threshold. For example, a cluster
  of 10 pathways may have only one significant pathway for a given
  enrichment set – therefore this enrichment color would not be included
  in the cnet cluster colors.

## multienrichjam 0.0.32.900

### changes to existing functions

- `mem_gene_pathway_heatmap()` cleaned up the heatmap overall:

  - left and top annotations have zero gap between columns and rows
  - left and top annotation legends have discrete color bars; for left
    it is always `c(0,1)`; for the top it uses -log10 integer steps.
  - top annotation legend appends `"-log10P"` to the label, to indicate
    that the color values are based upon
    [`log10()`](https://rdrr.io/r/base/Log.html) transform of the
    enrichment P-value.
  - the color legend has discrete color bar steps, indicating the number
    of enrichments where each gene is involved

## multienrichjam 0.0.31.900

### changes to existing functions

- `mem_gene_pathway_heatmap()` new argument `enrich_gene_weight` used to
  adjust the relative influence of the enrichment `-log10 P-value` and
  the gene incidence matrix on the column clustering. The effect is to
  adjust how much the enrichment P-values, or the gene content, affects
  the clusters. In principle both should have similar effects, but
  sometimes it helps to favor gene incidence or pathway enrichment.

## multienrichjam 0.0.30.900

### changes to existing functions

- `mem_gene_pathway_heatmap()` was altered to allow filtering pathways
  for a minimum number of genes represented by at least one enrichment
  result. For example, a pathway may be required to have at least 4
  genes, despite having a statistically significant enrichment P-value.
  The new argument `min_set_ct_each` tests each enrichment to see if any
  one has enough genes per pathway. There are two main effects of this
  filter: The gene-pathway heatmap will display fewer pathways; and any
  resulting Cnet plots will have the same pathways removed.

## multienrichjam 0.0.29.900

### changes to existing functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new argument `byCols` which is passed along to
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
  and is used to sort the pathway names within each Cnet cluster. The
  default uses `"composite_rank"`, which sorts by the
  `floor(-log10(pvalue))` then by the highest number of genes per
  pathway. The [`floor()`](https://rdrr.io/r/base/Round.html) function
  effectively sorts by the order of magnitude of the enrichment P-value,
  dropping the details. The obvious alternative is `"minp_rank"` which
  uses the `-log10(pvalue)` directly, which therefore does not
  effectively sort by gene count since enrichment P-values rarely tie.
  In our experience, between two pathways with reasonably similar
  enrichment P-value (within one order of magnitude such as 2.4e-5 and
  3.5e-5) the pathway with more genes was usually the more
  interesting/relevant biological pathway.

## multienrichjam 0.0.28.900

### new functions

[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md)
and
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md)
are the reciprocal functions to
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
and
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md).
The new functions convert incidence matrix to list, or signed incidence
matrix to list, respectively. They’re notable because they’re blazing
fast in our testing, thanks to efficient methods from the `arules` R
package for interconverting list to compressed logical matrix, and vice
versa.

### changes to existing functions

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  was modified to use `p_floor` as the ceiling for the row hierarchical
  clustering, previously it used `ceiling=3` which caused the dendrogram
  to have `height=0` when all enrichment results were lower than 0.001.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  was modified to calculate its own matrix colors when
  `color_by_column=TRUE`, until the
  [`colorjam::matrix2heatColors()`](https://jmw86069.github.io/colorjam/reference/matrix2heatColors.html)
  can be updated to the improved method.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  now correctly honors `p_floor`.

## multienrichjam 0.0.27.900

Added a new TODO.md file to track some new feature ideas.

### new functions

- [`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)
  infers the proper order of colors from a list of color vectors. It is
  called by
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  where attributes `"pie.color"` and `"coloredrect.color"` contain
  subsets of colors, but in the proper order. This function returns all
  colors in their proper order.
- [`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)
  takes a list of colors, and returns the unique colors, ordered by the
  overall order inferred from the list. It is internally called by
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  when argument `colorV` is not suppleid, since the attribute
  `"pie.color"` can be used to infer the correct order of colors.

### Changes to existing functions

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  now uses `"; "` as delimiter for the `"set_names"` attribute, to make
  it easier to distinguish individual pathway names.
- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  was refactored to handle specific color order defined either by
  argument `colorV` or by inferring the correct order of colors from
  attributes such as `"pie.color"`. This change should help order node
  colors more consistent to the original input colors to
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md),
  specifically the argument `colorV`.
- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  now has an example that shows the effects of ordering pie nodes by
  color.

## multienrichjam 0.0.26.900

### Changes to existing functions

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  now includes the cluster name by default with the `set_labels` which
  are used to display the top `n` pathway names for each collapsed Cnet
  node. For example
  `"Cluster A: Aryl Hydrocarbons; Granzyme Signaling"`.
- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  arguments `max_labels` and `max_char_labels` accept multiple values,
  which are recycled and applied to filter each cluster. For example
  `max_labels=c(5,2)` will filter the first cluster to display up to 2
  labels, and the second cluster to display up to 5 labels.

## multienrichjam 0.0.25.900

### Changes to existing functions

- Thanks to `simang5c` the function calls to
  [`renameColumn()`](https://jmw86069.github.io/jamba/reference/renameColumn.html)
  were corrected to
  [`jamba::renameColumn()`](https://jmw86069.github.io/jamba/reference/renameColumn.html),
  to resolve issue
  [\#1](http://github.com/jmw86069/multienrichjam/issues/1). Also fixed
  several other references to `jamba` functions:
  [`nameVector()`](https://jmw86069.github.io/jamba/reference/nameVector.html),
  [`nameVectorN()`](https://jmw86069.github.io/jamba/reference/nameVectorN.html),
  [`makeNames()`](https://jmw86069.github.io/jamba/reference/makeNames.html),
  [`printDebug()`](https://jmw86069.github.io/jamba/reference/printDebug.html),
  [`mixedSortDF()`](https://jmw86069.github.io/jamba/reference/mixedSortDF.html),
  [`cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html),
  [`igrepHas()`](https://jmw86069.github.io/jamba/reference/igrepHas.html),
  [`provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html),
  [`vigrep()`](https://jmw86069.github.io/jamba/reference/vigrep.html),
  [`unvigrep()`](https://jmw86069.github.io/jamba/reference/unvigrep.html),
  [`rbindList()`](https://jmw86069.github.io/jamba/reference/rbindList.html),
  [`noiseFloor()`](https://jmw86069.github.io/jamba/reference/noiseFloor.html),
  [`rmNA()`](https://jmw86069.github.io/jamba/reference/rmNA.html),
  [`rmNULL()`](https://jmw86069.github.io/jamba/reference/rmNULL.html),
  [`normScale()`](https://jmw86069.github.io/jamba/reference/normScale.html),
  [`getColorRamp()`](https://jmw86069.github.io/jamba/reference/getColorRamp.html),
  [`deg2rad()`](https://jmw86069.github.io/jamba/reference/deg2rad.html).
  Omg! There were so many more cases than I thought! (Package-building
  options to import functions from packages seem to go too far, the
  imported functions appear to be contained in the new package which
  seems misleading… Ah well.)
- Also made several small changes to handle single-enrichment input to
  multiEnrichMap(). It still provides some useful graphical benefits
  even with only one enrichment input.

## multienrichjam 0.0.24.900

### Changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  was updated to enable analysis with only one enrichment result.
  Testing whether the downstream capabilities are useful in the context
  of one enrichment result.

## multienrichjam 0.0.23.900

### Changes to existing functions

- Several functions have additional
  [`jamba::printDebug()`](https://jmw86069.github.io/jamba/reference/printDebug.html)
  output when `verbose=TRUE`.
- [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
  and
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  new argument `use_shadowText=TRUE` will enable
  [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
  to replace [`graphics::text()`](https://rdrr.io/r/graphics/text.html)
  for igraph text, which affects node and edge labels. The goal is to
  make labels more widely legible when they are placed on top of light
  and dark colors, common when using dark colored nodes on a white
  background.
- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  now tries harder to respect a pre-defined xlim and ylim, before using
  the range of layout coordinates
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  slightly adjusted the number of pathway clusters to use based upon the
  number of columns in the gene-pathway matrix, slightly increasing the
  number for low number of columns.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new argument pathway_column_split allows setting a specific number of
  pathway clusters.
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  argument `do_which` allows creating only the requested plots from the
  full sequence of plots. This function is likely to become the core
  part of the analysis workflow:

> gene-pathway heatmap -\> pathway clusters -\> Cnet using pathway
> clusters

## multienrichjam 0.0.22.900

### New functions

- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  provides a new `igraph` vertex shape `"jampie"`. It is a complete
  clone to the default `igraph` shape `"pie"` except that it offers
  vectorized plotting, which can be substantially faster for large
  `igraph` objects that use `"pie"` nodes.
- [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
  is a similar clone to
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
  except that it uses vectorized plotting specifically when there are
  multiple shapes in the same `igraph` object. Formerly, each node is
  drawn individually which is substantially slower, especially for
  `igraph` objects with more than 100 nodes. It also by default converts
  shape `"pie"` to `"jampie"` in order to use vectorized plotting via
  [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  above.

### Changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  new argument `plot_function` allows use of custom plot function,
  specifically to use
  [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md).

### Other changes

- Functions for igraph vertex shapes were moved into a new .R file, and
  into a new function family `"jam igraph shapes"`.

## multienrichjam 0.0.21.900

### New functions

- [`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md)
  draws a grid object with title and optional subtitle, allowing for
  grid objects such as
  [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  objects, or even multiple heatmaps. That said, any grid `"gTree"`
  object should work.

## multienrichjam 0.0.20.900

### New functions

- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  is a new function that creates the current recommended folio of
  multienrichment plots:

1.  Enrichment P-value heatmap
2.  Gene-pathway heatmap – importantly with pathway clustering
3.  Cnet using the pathway clusters, collapsed by cluster
4.  Cnet using the pathway clusters, exemplar pathways per cluster
5.  Cnet using each individual pathway cluster

### Changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  now includes `p_cutoff` in the output, which represents the enrichment
  P-value threshold used for the analysis. This cutoff is useful in
  making other color gradients respect the same threshold required for
  significant enrichment results, so P-values that do not meet this
  threshold can be colored white (or the background color.)
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  new argument (minor release) `do_which` to help produce selected plot
  pages from a folio of plots.

## multienrichjam 0.0.19.900

### New functions

- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  is an intriguing new method of simplifying the numerous pathways, by
  using each pathway cluster from `mem_gene_pathway_heatmap()`. Each
  cluster is condensed to one result, combining all genes in each
  cluster. The results are surprisingly insightful, especially when
  numerous pathways are present per cluster. This function also calls
  [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
  so it’s possible to pick a handful of the top pathways per cluster, as
  relevant.

## multienrichjam 0.0.18.900

### New functions

- [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
  is a convenience function which ranks pathways/sets within clusters,
  useful with a list of clusters following `mem_gene_pathway_heatmap()`.
  It makes it easier to rank pathways within a cluster, potentially
  choosing one exemplar pathway to represent each cluster. This function
  is part of more effort to streamline the overall analysis workflow.

### Bug fixes

- Fixed small issue with
  [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  which was properly applying `min_gene_ct` and `min_path_ct` however it
  did not check the resulting data to remove empty columns and rows.
  This change has been made.

## multienrichjam 0.0.17.900

### New functions

- [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  takes an igraph, determines a color for each node based upon the shape
  (using `avg_colors_from_list()` for shape `"pie"` and
  `"coloredrectangle"`), then creates an average color between two
  nodes, and uses that as the edge color.

### Changes to existing functions

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  has more configurable arguments, in the form of named lists. The
  driving example is
  `label_factor_l=list(nodeType=c(Gene=0.01, Set=1))`, that will apply
  `label.cex*0.01` to “Gene” nodes, and `label.cex*1.5` to “Set” nodes.
  Similar arguments: `node_factor_l` to apply to node size, and
  `label_dist_factor_l` to apply to label distance from node center.
- [`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)
  was updated to fix a small bug in order of rendering colored
  rectangles, it was somehow drawing the frame after the fill colors,
  which caused the frames to overlap each other and appear transparent.

## multienrichjam 0.0.16.900

The vignette was updated to use the new functions, for a much cleaner
overall workflow.

### Changes to existing functions

- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  uses new function
  [`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
  to create one color to represent each node, then
  [`sort_colors()`](https://jmw86069.github.io/colorjam/reference/sort_colors.html)
  to sort nodes by color hue, instead of previous behavior which sorted
  by the hex color string. This change affects attributes `"pie.color"`,
  `"coloredrect.color"` and `"color"`. Recall that this function
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  is intended to help visually organize groups of nodes by their colors,
  to make it easier to tell how many nodes are each color.
- [`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
  now uses
  [`weighted.mean()`](https://rdrr.io/r/stats/weighted.mean.html) based
  upon the `"c"` channel (chroma, color saturation) in order to
  down-weight the effect of grey colors, which really should have no
  hue. The default weight for grey is 0.1, while the maximum value for
  fully saturated colors is 100. Blending `"green"` with `"grey30"`
  yields slightly less saturated green, as expected.
- [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
  now includes convenience calls to
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
  and
  [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md).
  In almost all cases, creating a subset of a Cnet necessitates a new
  layout, which therefore requires new node ordering (sorting groups of
  equivalent nodes by color), and then label orientation (the angle used
  to offset the node label from the center of each node).

### New functions

- [`sort_colors()`](https://jmw86069.github.io/colorjam/reference/sort_colors.html)
  and
  [`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md)
  uses the
  [`farver::encode_colour()`](https://farver.data-imaginist.com/reference/encode_colour.html)
  function to convert colors to HCL, then sorts by hue, and returns the
  original sorted vector or corresponding order, respectively. For
  sorting colors inside a `data.frame`, convert colors to a factor whose
  levels are `sort_colors(unique(colors))`, in order to maintain ties
  during multi-column sorting.
- `apply_color_cap()` imposes a numeric range onto a color channel for a
  vector or list of colors, using
  [`jamba::noiseFloor()`](https://jmw86069.github.io/jamba/reference/noiseFloor.html).
  Values outside the allowed range are forced to the range.
- [`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)
  is a custom version of
  [`igraph::induced_subgraph()`](https://r.igraph.org/reference/subgraph.html),
  that correctly handles subsetting the graph layout when the graph
  itself is subsetted.
- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  is a custom `igraph::plot()` with proper sizing when `rescale=FALSE`,
  allowing better aspect ratios for certain network layout coordinates.
- [`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
  is a wrapper to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  that returns a layout specification object, for use with
  [`igraph::add_layout_()`](https://r.igraph.org/reference/add_layout_.html)
  in a more confusing workflow pattern than I expected.

## multienrichjam 0.0.15.900

### Changes to existing functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  several updates to improve robustness in the overall workflow.

### new functions

- [`heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
  and
  [`heatmap_column_order()`](https://jmw86069.github.io/jamba/reference/heatmap_column_order.html)
  provide enhanced output from
  [`ComplexHeatmap::row_order()`](https://rdrr.io/pkg/ComplexHeatmap/man/row_order-dispatch.html),
  mainly that they return the actual rownames and colnames of data in
  the heatmap. Very useful when the data used for the heatmap has been
  filtered internal to the function.

## multienrichjam 0.0.14.900

### Changes to existing functions

- [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
  removed argument `makeUnique` because the underlying conversion to
  matrix no longer requires that step, and it was a performance hit for
  extremely large lists. The argument `keepCounts` is the only
  requirement to maintain the count of each entry per list.
- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  now subsets the `geneIM` gene incidence matrix to match the genes
  after `topEnrichN` filtering is applied.
- [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  now by default will use edge attribute `"weight"` as the `weights`
  argument when calling
  [`qgraph::qgraph.layout.fruchtermanreingold()`](https://rdrr.io/pkg/qgraph/man/qgraph.layout.fruchtermanreingold.html).
- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  was enhanced for more robust import conditions, specifically for
  different variations of missing IPA enrichment results. Also, empty
  colnames are removed, to help recognize the proper identifier in each
  scenario.

### New functions

- [`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md)
  takes a list of `enrichResult` and returns the list of genes
  represented in each `enrichResult`. Intended mainly for internal use
  by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md).

## multienrichjam 0.0.13.900

### New plotting functions

- [`mem_multienrichplot()`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md)
  allows customized enrichMap-style plotting of igraph objects. Notably,
  you can filter by Jaccard overlap, or by overlap count – the number of
  genes involved in the overlap. There are otherwise a lot of overlaps
  that involve only one gene, which is not the best way to build this
  type of network.
- [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)
  is a light extension of
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
  but adds default calls to
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
  and
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  to help things look pretty.
- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  draws a color legend in the corner of a figure, using the colors
  defined in the
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  output.

### Other new functions

- [`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md)
  subsets an igraph based upon connected components – i.e. distinct
  subclusters.
- [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  takes the pathway-gene incidence matrix and produces a Cnet plot
  `igraph` object. If given the output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  it will also color nodes using the gene and enrichment incidence
  matrix colors.
- [`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
  and
  [`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md)
  are intended for very rapid color blending.

## multienrichjam 0.0.11.900

### Bug fixed

- `subsetCnetGraph()` was fixed to handle rare cases where pathway set
  name is identical to one or more genes, which happens with IPA pathway
  analysis, in the “Upstream Regulators” output.

## multienrichjam 0.0.10.900

### New functions

- `mem_gene_pathway_heatmap()` is a wrapper to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  which takes output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  and produces a heatmap styled similar to a pathway-gene incidence
  matrix. It includes row and column annotations to help interpret the
  results.
- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  produces a heatmap of the pathways and enrichment P-values for each
  comparison.
- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  is a general use function for igraph objects, it takes an igraph
  network and associated layout, and arranges igraph nodes at an angel
  opposite the majority of edges from each node. The result arranges
  labels around the outside of the network, and typically away from
  other nodes. It isn’t perfect, but is visually a big step in the right
  direction.
- [`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md)
  is a simple function that removes nodes that have no edges. I found
  myself doing it manually enough times I wanted something quick and
  easy.

### enhancements

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  now includes two new edge attributes on MultiEnrichMap igraphs:
  overlap_count which contains the number of shared genes between two
  sets, and overlap_max_pct which contains the max percent overlap
  between two sets – based upon the \# overlapped divided by the smaller
  of the two sets. The overlap_count is useful as an optional edge
  label, to show how many genes are involved.
- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  now forces `"setSize"` to contain integer values, useful when pathway
  size is inferred from something like gene ratio.

### Bug fixes

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  was not properly assigning `colnames(memIM)`, although they were
  inferred from the `rownames(enrichIM)`. The bug was because
  [`as.character()`](https://rdrr.io/r/base/character.html) removes
  names from a character vector, and the conversion was done to enforce
  proper handling by
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html) which gives
  unexpected results when the list may contain a factor.

## multienrichjam 0.0.9.900

### bug fixes

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  no longer calls commandline `grep` to remove blank rows, for now it
  uses [`readLines()`](https://rdrr.io/r/base/readLines.html) and
  [`jamba::vigrep()`](https://jmw86069.github.io/jamba/reference/vigrep.html)
  to select rows with at least one character.

## multienrichjam 0.0.8.900

### bug fixes

- Fixed several calls to
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html) that failed when
  sent factors, which only happens when R
  `options("stringsAsFactors"=TRUE)` which is the default in base R. Now
  all calls to [`strsplit()`](https://rdrr.io/r/base/strsplit.html)
  enforce [`as.character()`](https://rdrr.io/r/base/character.html)
  unless character type has already been enforced.
- Added numerous package prefixes to functions, to avoid importing all
  package dependencies.

## multienrichjam 0.0.7.900

### bug fixes

- [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  were updated to handle more default colnames for P-values, now
  including the default colnames for `enrichResult` objects
  `c("pvalue","padjust")`. When no sort columns are found, a warning
  message is printed.

### changes

- Added documentation for
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  and other functions.
- Added vignette describing the workflow starting with Ingenuity IPA
  enrichment results.

## multienrichjam 0.0.6.900

### new functions

- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  imports Ingenuity IPA enrichment results, by default splitting each
  enrichment table into its own `data.frame`. It curates colnames to be
  consistent with downstream analyses.
- [`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md)
  will curate the colnames of IPA data, and ensures the values in the
  gene column are consistently delimited.
- [`gsubs()`](https://jmw86069.github.io/jamba/reference/gsubs.html) is
  similar to [`base::gsub()`](https://rdrr.io/r/base/grep.html) except
  that it applies a vector of pattern-replacement operations in order.
  This function may be moved to the `"jamba"` package.
- [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  is a helper function to find a colname given a vector of expected
  values, which is matched directly, then case-insensitively, then as a
  vector of patterns to match the start, end, then any part of the
  colnames. By default the first matching value from the first
  successful method is returned. Good for matching c(“P-Value”,
  “pvalue”, “Pval”)
- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  and
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  subset the input pathway enrichment results by taking the top `n`
  result, then making sure the overall selected pathways are retained
  for all enrichment tables.

## multienrichjam 0.0.5.900

### changes

- [`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)
  is no longer dependent upon the input list `x` having names, nor
  having unique names.

### enhancements

- [`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md)
  now handles list input, which is helpful when applied to igraph
  objects where `"pie"` colors are accessed as a list.
- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
  was updated to use vectorized logic for `"pie"` vertex attributes.

## multienrichjam 0.0.4.900

### bug fixed

- Fixed small issues with
  [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md).

### changes

- Updated
  [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  to handle words the should be kept uppercase, with some defaults
  pre-configured, e.g. “mRNA”.

## multienrichjam 0.0.3.900

### changes

- removed
  [`mergeAllXY()`](https://jmw86069.github.io/jamba/reference/mergeAllXY.html)
  and
  [`unnestList()`](https://jmw86069.github.io/jamba/reference/unnestList.html)
  and moved them to the `jamba` package. Added corresponding version
  requirement on jamba.

## multienrichjam 0.0.2.900

### changes

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  was updated to include a wider set of canonical pathway prefixes used
  in MsigDB v6.1.

### new functions

- [`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md)
  and
  [`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md)
  are helper functions used to convert a Cnet igraph to either a
  incidence matrix, or a table summary useful for extracting subsets of
  pathways using various network descriptors.
- [`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md)
  and
  [`shape.ellipse.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.plot.md)
  add an igraph vertex shape “ellipse” whose shape is controlled by
  vertex.ellipse.ratio, where `1` creates a circular node.

## multienrichjam 0.0.1.900

### new functions

- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  the workhorse core function, implementing the full workflow to
  pre-process the results for later inspection.
- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  converts a data.frame into an `enrichResult` object as used in
  `clusterProfiler` and `DOSE`.
- [`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md)
  custom function to create a cnet plot igraph object.
- [`mergeAllXY()`](https://jmw86069.github.io/jamba/reference/mergeAllXY.html)
  which merges a list of data.frames while keeping all rows.
- [`unnestList()`](https://jmw86069.github.io/jamba/reference/unnestList.html)
  which un-nests a list of lists, resulting in a flattened list. It
  directly supports
  [`mergeAllXY()`](https://jmw86069.github.io/jamba/reference/mergeAllXY.html)
  in order to provide a simple list of data.frame objects from
  potentially nested list of lists of data.frame objects.
- [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
  converts a list of enrichResult objects (or data.frames) into an
  indidence matrix of gene rows and pathway columns.
- [`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
  combines a list of enrichResult objects (or data.frames) into a single
  data.frame, using the union of genes, and the best P-value.
- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  to update pathway/gene set labels using some small set of logic.
- [`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md)
  checks if a color is a blank color, either by comparison to known
  colors, or transparency, or saturation/brightness.
- [`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md)
  converts an igraph into one with pie node shapes, optionally
  coloredrectangle node shapes.
- [`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md)
  is an extension to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  that returns a function, therefore convenient to use in
  `igraph::plot()` with custom arguments.
- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
  removes blank colors in either pie or coloredrectangle node shapes,
  thereby making non-blank colors easier to see in plots.
- [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
  takes an igraph with Set and Gene nodeType, and subsets based upon a
  fixed list of Set or Gene nodes, removing singlet disconnected nodes
  from the output.

### changes

- added `dplyr` to package dependencies.
- added `DOSE` to package dependencies.
