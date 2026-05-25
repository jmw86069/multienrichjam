# AllGenerics.R
#
# Todo:
# geneIM()<-
# geneIMcolors()<-
# geneIMdirection()<-
# enrichIM()<-
# enrichIMcolors()<-
# enrichIMdirection()<-
# enrichIMgeneCount()<-

setGeneric("enrichments", signature="x",
   function(x, ...) standardGeneric("enrichments")
)

setGeneric("enrichments<-", signature="x",
   function(x, value) standardGeneric("enrichments<-")
)

setGeneric("sets", signature="x",
   function(x, ...) standardGeneric("sets")
)

setGeneric("sets<-", signature="x",
   function(x, value) standardGeneric("sets<-")
)

setGeneric("genes", signature="x",
   function(x, ...) standardGeneric("genes")
)

setGeneric("genes<-", signature="x",
   function(x, value) standardGeneric("genes<-")
)

setGeneric("geneIM", signature="x",
   function(x, ...) standardGeneric("geneIM")
)

setGeneric("geneIMdirection", signature="x",
   function(x, ...) standardGeneric("geneIMdirection")
)

setGeneric("geneIMdirection<-", signature="x",
   function(x, ...) standardGeneric("geneIMdirection<-")
)

setGeneric("geneIMcolors", signature="x",
   function(x, ...) standardGeneric("geneIMcolors")
)

setGeneric("geneIMcolors<-", signature="x",
   function(x, ...) standardGeneric("geneIMcolors<-")
)

setGeneric("enrichList", signature="x",
   function(x, ...) standardGeneric("enrichList")
)

setGeneric("enrichIM", signature="x",
   function(x, ...) standardGeneric("enrichIM")
)

setGeneric("enrichIMcolors", signature="x",
   function(x, ...) standardGeneric("enrichIMcolors")
)

setGeneric("enrichIMcolors<-", signature="x",
   function(x, ...) standardGeneric("enrichIMcolors<-")
)

setGeneric("enrichIMdirection", signature="x",
   function(x, ...) standardGeneric("enrichIMdirection")
)

setGeneric("enrichIMdirection<-", signature="x",
   function(x, ...) standardGeneric("enrichIMdirection<-")
)

setGeneric("enrichIMgeneCount", signature="x",
   function(x, ...) standardGeneric("enrichIMgeneCount")
)

setGeneric("memIM", signature="x",
   function(x, ...) standardGeneric("memIM")
)

setGeneric("geneHitIM", signature="x",
   function(x, ...) standardGeneric("geneHitIM")
)

setGeneric("geneHitIM<-", signature="x",
   function(x, ...) standardGeneric("geneHitIM<-")
)

setGeneric("geneHitList", signature="x",
   function(x, ...) standardGeneric("geneHitList")
)

setGeneric("geneHitList<-", signature="x",
   function(x, ...) standardGeneric("geneHitList<-")
)

setGeneric("headers", signature="x",
   function(x, ...) standardGeneric("headers")
)

setGeneric("colorV", signature="x",
   function(x, ...) standardGeneric("colorV")
)

setGeneric("colorV<-", signature="x",
   function(x, value) standardGeneric("colorV<-")
)

setGeneric("thresholds", signature="x",
   function(x, ...) standardGeneric("thresholds")
)

setGeneric("thresholds<-", signature="x",
   function(x, ...) standardGeneric("thresholds<-")
)

setGeneric("geneInCategory", signature="x",
   function(x) standardGeneric("geneInCategory")
)

setGeneric("setsByGene", signature="x",
   function(x) standardGeneric("setsByGene")
)

setGeneric("Clusters", signature="x",
   function(x) standardGeneric("Clusters")
)

setGeneric("ClusterLabels", signature="x",
   function(x) standardGeneric("ClusterLabels")
)

setGeneric("ClusterLabels<-", signature="x",
   function(x, value) standardGeneric("ClusterLabels<-")
)

setGeneric("ClusterData", signature="x",
   function(x) standardGeneric("ClusterData")
)

setGeneric("ClusterData<-", signature="x",
   function(x, value) standardGeneric("ClusterData<-")
)

setGeneric("GeneClusters", signature="x",
   function(x) standardGeneric("GeneClusters")
)

setGeneric("metadata", signature="x",
   function(x) standardGeneric("metadata")
)

setGeneric("EnrichmentHeatmap", signature=c("x", "do_plot", "params"),
   function(x, do_plot,
      params=list(width=30),
      ...) standardGeneric("EnrichmentHeatmap")
)

setGeneric("GenePathHeatmap", signature=c("x", "do_plot", "params"),
   function(x, do_plot,
      params=list(width=30),
      ...) standardGeneric("GenePathHeatmap")
)

setGeneric("CnetCollapsed", signature=c("x", "type", "do_plot", "params"),
   function(x, type, do_plot,
      params=list(
         width=30,
         maxNchar=Inf,
         layout=NULL,
         rotate_degrees=0), ...) standardGeneric("CnetCollapsed")
)

setGeneric("CnetExemplar", signature=c("x", "num", "do_plot", "params"),
   function(x, num, do_plot,
      params=list(
         width=30,
         maxNchar=Inf,
         layout=NULL,
         rotate_degrees=0),
      ...) standardGeneric("CnetExemplar")
)

setGeneric("CnetCluster", signature=c("x", "cluster", "do_plot", "params"),
   function(x, cluster, do_plot,
      params=list(
         width=30,
         maxNchar=Inf,
         layout=NULL,
         rotate_degrees=0),
      ...) standardGeneric("CnetCluster")
)

setGeneric("EnrichmentMap", signature=c("x", "do_plot", "params"),
   function(x, do_plot,
      params=list(repulse=3.5,
         width=30,
         group="default",
         mark.expand=4,
         do_legend=TRUE),
      ...) standardGeneric("EnrichmentMap")
)

setGeneric("Caption", signature=c("x"),
   function(x, ...) standardGeneric("Caption")
)

setGeneric("CaptionLegendList", signature=c("x"),
   function(x, ...) standardGeneric("CaptionLegendList")
)

setGeneric("CaptionLegendList<-", signature="x",
   function(x, ...) standardGeneric("CaptionLegendList<-")
)
