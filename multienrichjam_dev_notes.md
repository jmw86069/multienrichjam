
# multienrichjam dev notes

* The `.onLoad()` function defined in `zzz.R` adds igraph shapes
'jampie', 'ellipse', and 'coloredrectangle'.

   * When updating the shape functions while using 'devtools' for active
   development, you must also add the shape again to the R session.
   
   ```
   igraph::add_shape("jampie",
      clip=shape.jampie.clip,
      plot=shape.jampie.plot)
   ```

* Adding generic methods for S4 objects

   * Test for existing generic function, the 'getName=TRUE' part tells what
   package provided the generic function
   `isGeneric("score", getName=TRUE)`
   `isGeneric("abbreviate", getName=TRUE)`

* Adding data to include in the R package, and use with R function examples

   * Create the R object in the active R session, for example 'Memtest', then
   `usethis::use_data(Memtest, overwrite=TRUE)`

* Migrating from mem `list` to `Mem` S4 object

* Clean up argument names in `multiEnrichMap()`

   * Change `cutoffRowMinP` to `p_cutoff` 
