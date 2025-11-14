
# Updated Conceptual Workflow (Nov 2025)

Starting data

   * `enrichResult` from `clusterProfiler`, or
   * `importIPAenrichment()`

Run `multiEnrichMap()`

   - output: `Mem`

Run `prepare_folio()` (without plotting)
or `mem_plot_folio()` (with plotting)

   - output: `MemPlotFolio`

Explore specific plots

   * `GenePathHeatmap()`
   * `EnrichmentHeatmap()`
   * `CnetCollapsed()`
   * `CnetExemplar()`
   * `CnetCluster()`

Create Cnet with Specific Pathways

   * `mem2cnet(Mem)`

Customize Cnet Layout

   * `launch_shinycat()`
