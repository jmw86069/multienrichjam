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
