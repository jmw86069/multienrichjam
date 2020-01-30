
## Keep track of what igraph functions call other igraph functions

subsetCnetIgraph()
- removeIgraphSinglets() - optional, default=TRUE
- removeIgraphBlanks()   - optional, default=TRUE
- spread_igraph_labels() - optional, default=TRUE
- reorderIgraphNodes()   - optional, default=TRUE

spread_igraph_labels()
- layout_with_qfr() - called if force_relayout=TRUE, OR (layout=NULL and g$layout=NULL)
   - layout() - only called if force_relayout=TRUE
- reorderIgraphNodes()   - optional, default=TRUE


reorderIgraphNodes()
- returns g which always includes g$layout with layout coordinates
- layout() - called if layout is a function
- relayout_with_qfr() - called if layout=NULL, g$layout=NULL
- isColorBlank() - called with isColorBlank(...)
