
.onLoad <- function
(libname,
   pkgname)
{
   ## define new igraph vertex shape "coloredrectangle"
   # confirm shape has not already been added
   if (!"coloredrectangle" %in% igraph::shapes()) {
      igraph::add_shape("coloredrectangle",
         # clip=igraph::shape_noclip,
         clip=shape.coloredrectangle.clip,
         plot=shape.coloredrectangle.plot);
   }

   ## define new igraph vertex shape "coloredrectangle"
   if (!"ellipse" %in% igraph::shapes()) {
      igraph::add_shape("ellipse",
         clip=shape.ellipse.clip,
         # clip=igraph::shape_noclip,
         plot=shape.ellipse.plot);
   }

   ## define new igraph vertex shape "jampie"
   if (!"jampie" %in% igraph::shapes()) {
      igraph::add_shape("jampie",
         clip=shape.jampie.clip,
         plot=shape.jampie.plot);
   }
}
