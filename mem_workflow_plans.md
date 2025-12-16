
Need schematic or simple workflow:


Import / enrichResult -->

multiEnrichMap() -->

prepare_folio() -->

- EnrichmentHeatmap()
- GenePathHeatmap()
- CnetCollapsed()
- CnetExemplar()
- CnetCluster()


## New function

* `CnetHeatmaps()`

	* Expression matrix - centered data - rows match `genes(Mem)`
	* Expression column groups - match `enrichments(Mem)`
	* Cnet data: `MemPlotFolio` or cnet `igraph`
	* Lots of settings
	
		* column widths, heatmap placement


## Add examples

* Cnet with custom exemplars

	* Currently:
	`mem2cnet(Mem[, custom_pathways, ])`

* Cnet with custom subset of pathways with some groupings

	* Currently:
	
	```
	pathway_groups <- list(
	   A=pathways_A,
	   B=pathways_B,
	   C=pathways_C
	)
	Mpf <- prepare_folio(Mem,
	   pathway_column_split=pathway_groups);
	CnetCollapsed(Mpf)
	```
