
# 0.0.101.900: no longer skip this step when the shape already exists,
#    presuming that calling this function is intended to overwrite
#    any pre-existing node shapes.

.onLoad <- function
(libname,
 pkgname)
{
   ## Retrieve known igraph shapes, only used when checking
   ## for pre-existing shapes.
   # suppressPackageStartupMessages(igraph_shapes <- igraph::shapes());
   
   # if (!"coloredrectangle" %in% igraph_shapes) {
   # }
   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("coloredrectangle",
      clip=shape.coloredrectangle.clip,
      plot=shape.coloredrectangle.plot);

   ## define new igraph vertex shape "coloredrectangle"
   igraph::add_shape("ellipse",
      clip=shape.ellipse.clip,
      plot=shape.ellipse.plot);

   ## define new igraph vertex shape "jampie"
   igraph::add_shape("jampie",
      clip=shape.jampie.clip,
      plot=shape.jampie.plot);
   
}
