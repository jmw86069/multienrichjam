# Package index

## User-Facing Functions

### Core functions

Functions central to multienrichjam

- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  : Jam custom function to plot an igraph network
- [`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  [`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
  : Multienrichment folio of summary plots
- [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
  : Prepare MultiEnrichMap data from enrichList

### Import functions

Functions to use external data in multienrichjam

- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  : Convert data.frame to enrichResult
- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  : Import Ingenuity Pathway Analysis 'IPA' results
- [`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md)
  : Convert IPA list to a gene hit list or matrix

## Other Useful Functions

### Custom plot functions

Customized plotting for ‘Mem’ results

- [`mem_enrichment_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_enrichment_heatmap.md)
  : MultiEnrichment Heatmap of enrichment P-values
- [`mem_gene_path_heatmap()`](https://jmw86069.github.io/multienrichjam/reference/mem_gene_path_heatmap.md)
  : MultiEnrichment Heatmap of Genes and Pathways
- [`mem_legend()`](https://jmw86069.github.io/multienrichjam/reference/mem_legend.md)
  : MultiEnrichMap color legend
- [`plot_mpf()`](https://jmw86069.github.io/multienrichjam/reference/plot_mpf.md)
  : Plot function for MemPlotFolio objects

### igraph utilities

Utilities for igraph objects

- [`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md)
  : Color edges by nodegroups
- [`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)
  : Colorize igraph edges by nodes
- [`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md)
  : Color edges by nodegroups
- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  : Get Cnet node set by connected Sets
- [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  : Jam custom function to plot an igraph network
- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)
  : Launch ShinyCAT: Cnet Adjustment Tool
- [`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)
  : Subgraph using Jam extended logic

### igraph layout functions

Layout functions for igraph objects

- [`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md)
  : Layout igraph communities
- [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  : igraph layout using qgraph Fruchterman-Reingold
- [`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md)
  : igraph layout function using qgraph Fruchterman-Reingold
- [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)
  : Nudge igraph layout by node
- [`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)
  : igraph re-layout using qgraph Fruchterman-Reingold
- [`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md)
  : Remove igraph singlet nodes
- [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  : Re-order igraph nodes
- [`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)
  : Rotate igraph layout coordinates
- [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  : Spread igraph node labels by angle from node center

### Cnet utilities

Utilities for Cnet igraph objects

- [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  : Adjust Cnet node set
- [`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md)
  : Adjust Set nodes then relayout Gene nodes
- [`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md)
  : Apply minimum node spacing for each Cnet node set
- [`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md)
  : Bulk Cnet plot adjustments
- [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  : Get Cnet node set by connected Sets
- [`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md)
  : Get Cnet nodesets as a named vector
- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)
  : Launch ShinyCAT: Cnet Adjustment Tool
- [`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md)
  : Make Cnet test igraph
- [`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)
  : Relayout each nodegroup in a bipartite (Cnet) graph, experimental

### Mem utilities

Functions useful for ‘Mem’ results

- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  [`words`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  [`abbrev`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  : Fix Set or pathway labels for legibility
- [`` `[`( ``*`<Mem>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`show(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`Mem_to_list()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`as.list(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`names(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichments(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichments<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`sets(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `sets<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`genes(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `genes<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIMdirection(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneIMdirection<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIMcolors(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneIMcolors<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichList(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMcolors(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichIMcolors<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMdirection(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichIMdirection<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMgeneCount(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`memIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneHitIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneHitIM<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneHitList(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneHitList<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`headers(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`colorV(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `colorV<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`thresholds(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `thresholds<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`dim(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`dimnames(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`list_to_Mem()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`updateObject(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneInCategory(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`setsByGene(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`EnrichmentMap(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  : Mem S4 class, accessors, getters, and setters
- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  : Convert MultiEnrichment incidence matrix to Cnet plot
- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  : Convert multiEnrichMap mem output to EnrichmentMap emap
- [`list_to_MemPlotFolio()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`MemPlotFolio_to_list()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`show(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`Clusters(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`GeneClusters(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`ClusterLabels(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`` `ClusterLabels<-`( ``*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`ClusterData(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`` `ClusterData<-`( ``*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`thresholds(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`metadata(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`Caption(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CaptionLegendList(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`EnrichmentHeatmap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`GenePathHeatmap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetCollapsed(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetExemplar(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetCluster(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`plot(`*`<MemPlotFolio>`*`,`*`<ANY>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`EnrichmentMap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  : Mem S4 class, accessors, getters, and setters

## Functions Deep in the Details

### Enrichment-supporting functions

Functions that manipulate enrichment or multiple enrichment results.

- [`add_pathway_direction()`](https://jmw86069.github.io/multienrichjam/reference/add_pathway_direction.md)
  : Add directionality to pathway enrichment
- [`topEnrichBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  : Subset enrichResult for top enrichment results by source

### Plot support functions

Functions to plot MultiEnrichment results

- [`adjust_polygon_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md)
  [`adjust_rect_border()`](https://jmw86069.github.io/multienrichjam/reference/adjust_polygon_border.md)
  : Adjust polygon border to inner or outer edge
- [`plot_layout_scale()`](https://jmw86069.github.io/multienrichjam/reference/plot_layout_scale.md)
  : Plot layout scale by percentage of coordinate range

### Conversion functions

Functions to convert data types

- [`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md)
  : Summarize Cnet igraph as a data.frame
- [`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md)
  : Convert Cnet igraph to incidence matrix
- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  : Convert data.frame to enrichResult
- [`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
  : Convert enrichList to data.frame
- [`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
  : Convert enrichList to IM incidence matrix

### Detailed igraph-related functions

Functions that extend or customize igraph-related features

- [`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md)
  : Convert communities object to nodegroups list format
- [`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md)
  : Draw ellipse
- [`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md)
  : Bundle edges in a bipartite graph
- [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)
  : Bundle edges using node groups
- [`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  [`words`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  [`abbrev`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)
  : Fix Set or pathway labels for legibility
- [`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md)
  : Flip direction of igraph edges
- [`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md)
  : Get partite/connected graph nodesets
- [`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md)
  : Highlight edges connected to a node or nodes
- [`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md)
  : Convert igraph to use pie node shapes
- [`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md)
  : Assign labels to igraph communities
- [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  [`memIM2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
  : Convert MultiEnrichment incidence matrix to Cnet plot
- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
  : Convert multiEnrichMap mem output to EnrichmentMap emap
- [`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md)
  : Convert nodegroups list to communities object
- [`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md)
  : Convert pie igraph node shapes to coloredrectangle
- [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)
  : Remove igraph blank wedges
- [`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)
  : Subset Cnet igraph
- [`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md)
  : Subset igraph by connected components
- [`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
  : Sync igraph nodes and communities

### Cnet igraph-related functions

Functions specific to Cnet igraph objects

- [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
  : Apply Cnet border color by directionality
- [`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)
  : Rotate igraph layout coordinates

### R-shiny functions

Functions supporting ShinyCat

- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)
  : Launch ShinyCAT: Cnet Adjustment Tool
- [`shinycat_server()`](https://jmw86069.github.io/multienrichjam/reference/shinycat_server.md)
  : R-shiny server function for shinycat
- [`shinycat_ui()`](https://jmw86069.github.io/multienrichjam/reference/shinycat_ui.md)
  : R-shiny ui function for shinycat

### Mem functions

Functions supporting Mem S4 object use

- [`check_Mem()`](https://jmw86069.github.io/multienrichjam/reference/check_Mem.md)
  : Check Mem object
- [`` `[`( ``*`<Mem>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`show(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`Mem_to_list()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`as.list(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`names(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichments(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichments<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`sets(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `sets<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`genes(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `genes<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIMdirection(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneIMdirection<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneIMcolors(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneIMcolors<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichList(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMcolors(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichIMcolors<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMdirection(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `enrichIMdirection<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`enrichIMgeneCount(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`memIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneHitIM(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneHitIM<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneHitList(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `geneHitList<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`headers(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`colorV(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `colorV<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`thresholds(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`` `thresholds<-`( ``*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`dim(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`dimnames(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`list_to_Mem()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`updateObject(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`geneInCategory(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`setsByGene(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  [`EnrichmentMap(`*`<Mem>`*`)`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  : Mem S4 class, accessors, getters, and setters
- [`Memtest`](https://jmw86069.github.io/multienrichjam/reference/Mem-data.md)
  [`Reese_genes`](https://jmw86069.github.io/multienrichjam/reference/Mem-data.md)
  : Reese 2019 cohort multienrichment data
- [`Mem-slots`](https://jmw86069.github.io/multienrichjam/reference/Mem-slots.md)
  : Mem description of slots

### MemPlotFolio functions

Functions supporting MemPlotFolio S4 object use

- [`check_MemPlotFolio()`](https://jmw86069.github.io/multienrichjam/reference/check_MemPlotFolio.md)
  : Check MemPlotFolio object
- [`list_to_MemPlotFolio()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`MemPlotFolio_to_list()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`show(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`Clusters(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`GeneClusters(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`ClusterLabels(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`` `ClusterLabels<-`( ``*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`ClusterData(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`` `ClusterData<-`( ``*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`thresholds(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`metadata(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`Caption(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CaptionLegendList(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`EnrichmentHeatmap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`GenePathHeatmap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetCollapsed(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetExemplar(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`CnetCluster(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`plot(`*`<MemPlotFolio>`*`,`*`<ANY>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  [`EnrichmentMap(`*`<MemPlotFolio>`*`)`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
  : Mem S4 class, accessors, getters, and setters
- [`score_gene_path_clusters()`](https://jmw86069.github.io/multienrichjam/reference/score_gene_path_clusters.md)
  : Score Gene-Path Clusters in MemPlotFolio

### igraph vertex shapes

Functions that provide custom igraph vertex shapes

- [`jam_mypie()`](https://jmw86069.github.io/multienrichjam/reference/jam_mypie.md)
  : Vectorized mypie() function for igraph vertex pie polygons
- [`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)
  [`shape.coloredrectangle.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md)
  : custom igraph vertex shape coloredrectangle
- [`shape.ellipse.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.clip.md)
  : clip function for igraph vertex shape ellipse
- [`shape.ellipse.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.plot.md)
  : plot function for igraph vertex shape ellipse
- [`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  [`shape.jampie.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
  : custom igraph vertex shape jampie

### list functions

Functions that operate on list objects

- [`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)
  : Ordered colors from a list of color vectors
- [`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md)
  : convert incidence matrix to list
- [`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md)
  : convert signed incidence matrix to list
- [`list2concordance()`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md)
  : Convert list to concordance matrix
- [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
  : convert list to incidence matrix
- [`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)
  : convert list to directional incidence matrix

### Utility functions

Functions used by other multienrichjam functions

- [`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
  : Alpha shape calculation
- [`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md)
  : Average geometric angles
- [`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
  : Average colors by list
- [`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md)
  : ComplexHeatmap cell function with bivariant color
- [`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)
  : Collapse Multienrichment clusters
- [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)
  : Color ramp for bivariate colors
- [`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md)
  : Curate Ingenuity IPA colnames
- [`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md)
  : Deconcatenate delimited column values in a data.frame
- [`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md)
  : Display colors from bivariate color function
- [`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md)
  : Extract gene hit list from list of enrichResult
- [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  : Find colname by character string or pattern matching
- [`find_enrich_colnames()`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md)
  : Find enrichment colnames
- [`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md)
  : Get data for alpha hull (internal)
- [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)
  [`set_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)
  : Obtain or create layout for igraph object
- [`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md)
  : Pattern replacement with multiple patterns
- [`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md)
  : Handle igraph attribute parameter list
- [`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md)
  : Determine if colors are blank colors
- [`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md)
  : Display colors from bivariate color function
- [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
  : Make alpha hull from points
- [`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)
  : Find recommended overlap threshold for EnrichMap, experimental
- [`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md)
  : Order colors
- [`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)
  : Rank Multienrichment clusters
- [`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md)
  : Rotate numeric coordinates
- [`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md)
  : Summarize spacing between igraph nodes or node groups
- [`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md)
  [`local_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md)
  : Withr mimic of with_options() for ComplexHeatmap options
- [`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)
  : Get angle from origin to vector of x,y coordinates

### Import functions

Functions to import enrichment data

- [`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md)
  : Convert data.frame to enrichResult
- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
  : Import Ingenuity Pathway Analysis 'IPA' results
- [`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md)
  : Convert IPA list to a gene hit list or matrix

### Internal igraph functions

Functions typically for internal use by other Jam functions, or which
were internal igraph functions.

- [`default_igraph_values()`](https://jmw86069.github.io/multienrichjam/reference/default_igraph_values.md)
  : Default igraph parameter values
- [`get_igraph_arrow_mode()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_arrow_mode.md)
  : Get igraph arrow mode
- [`jam_igraph_arrows()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph_arrows.md)
  : Render igraph arrows
- [`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
  : Jam igraph vectorized plot function (internal)

### Deprecated functions

Functions moved or removed

- [`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md)
  : Create a cnetplot igraph object, deprecated
- [`enrichMapJam()`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md)
  : Create enrichMap igraph object from enrichResult, deprecated
- [`grid_with_title()`](https://jmw86069.github.io/multienrichjam/reference/grid_with_title.md)
  : Draw Heatmap with title and subtitle using grid viewports,
  deprecated
- [`mem_multienrichplot()`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md)
  : MultiEnrichMap plot
- [`with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)
  : Layout specification for Qgraph Fruchterman-Reingold, deprecated
