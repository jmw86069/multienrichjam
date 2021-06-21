
## Keep track of igraph-related functions that call other functions

`subsetCnetIgraph()`

- `removeIgraphSinglets()` - optional, default=TRUE
- `removeIgraphBlanks()`   - optional, default=TRUE
- `spread_igraph_labels()` - optional, default=TRUE
- `reorderIgraphNodes()`   - optional, default=TRUE

`spread_igraph_labels()`

- `layout_with_qfr()` - called if force_relayout=TRUE, OR (layout=NULL and g$layout=NULL)

   - `layout()` - only called if force_relayout=TRUE

- `reorderIgraphNodes()`   - optional, default=TRUE


`reorderIgraphNodes()`

- returns g which always includes `g$layout` with layout coordinates
- `layout()` - called if layout is a function
- `relayout_with_qfr()` - called if `layout=NULL, g$layout=NULL`
- `isColorBlank()` - called with `isColorBlank(...)`


## Re-consider function name style

Re-consider whether renaming functions with some consistent
and intuitive naming scheme could be beneficial.
I actually complain about my own inconsistency so I think
it best to rename functions sooner than later. Old functions
will be deprecated and aliased to the new function names.

### Function Style Option A: `type_verb_subtype()`

For me it is sometimes helpful to see the relevant object type
first in the function name, however this style may not be idea.
For example:

* `input_method_subtype()`
* `igraph_reorder_nodes()` instead of `reorderIgraphNodes()`
* `igraph_spread_labels()` instead of `spread_igraph_labels()`
* `igraph_subset_cnet()` instead of `subsetCnetIgraph()`
* `igraph_remove_singlets()` instead of `removeIgraphSinglets()`
* `igraph_remove_blanks()` instead of `removeIgraphBlanks()`
* `layout_with_qfr()` for consistency with `igraph::layout_with_()` convention.
* `relayout_with_qfr()` consistent with `layout_with_qfr()`
* `color_is_blank()` possible instead of `isColorBlank()`

### Function Style Option B: `verb_type_subtype()`

**This option is preferred because it matches Tidyverse convention**

A counter point, the tidyverse style guide suggests using verbs,
all lowercase with underscore '_' between words (not camelCase).
For example:

* `method_input_typeaffected()`
* `reorder_igraph_nodes()` instead of `reorderIgraphNodes()`
* `spread_igraph_labels()`
* `subset_cnet_igraph()` instead of `subsetCnetIgraph()`
* `remove_igraph_singlets()` instead of `removeIgraphSinglets()`
* `remove_igraph_blanks()` instead of `removeIgraphBlanks()`
* `layout_with_qfr()`
* `relayout_with_qfr()`
* `is_color_blank()` instead of `isColorBlank()`


## New alternate function `subset_igraph()`?

Instead of:

* `remove_igraph_singlets()`
* `remove_igraph_blanks()`

One common function:

* `subset_igraph()`

   *-* or `prune_igraph()` or `filter_igraph()`
   * `remove_singlets=FALSE`
   * `remove_blanks=FALSE`

* `subset_cnet_igraph()`

   * same as above except option to subset
   by specific entries by `V(g)$nodeType`.


## Tidyverse %>% and R-4.1 \> pipe compatibility

Where possible and where relevant the functions should be
compatible with tidyverse pipe `%>%`
and the new R-4.1 operator pipe `\>`

In practice when functions are called for their side-effects,
the function should return the same object supplied after
applying modifications.

* `relayout_with_qfr()`
* `spread_igraph_labels()`
* `subset_cnet_igraph()`
* `remove_igraph_singlets()`
* `remove_igraph_blanks()`


## Finding relevant functions

For convenience one can use `jamba::grepls("igraph.")` to find
functions relevant to `igraph` objects.


