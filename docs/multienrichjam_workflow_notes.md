# Updated Conceptual Workflow (Nov 2025)

Starting data

- `enrichResult` from `clusterProfiler`, or
- [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

Run
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- output: `Mem`

Run
[`prepare_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
(without plotting) or
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)
(with plotting)

- output: `MemPlotFolio`

Explore specific plots

- [`GenePathHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
- [`EnrichmentHeatmap()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
- [`CnetCollapsed()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
- [`CnetExemplar()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)
- [`CnetCluster()`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md)

Create Cnet with Specific Pathways

- `mem2cnet(Mem)`

Customize Cnet Layout

- [`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)

Todo:

- Custom
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)
- [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
- Custom
  [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
