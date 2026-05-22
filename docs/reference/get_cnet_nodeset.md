# Get Cnet node set by connected Sets

Get Cnet node set by connected Sets, using nodeType to delineate the two
types of nodes in the graph.

## Usage

``` r
get_cnet_nodeset(
  g,
  set_nodes = NULL,
  sep = ",",
  filter_set_only = TRUE,
  min_size = 1,
  ...
)
```

## Arguments

- g:

  `igraph` object specifically containing Cnet plot data, with node
  attribute `"nodeType"` that has values `c("Gene", "Set")`.

- set_nodes:

  `character` vector, default NULL, used to subset the returned data to
  nodes which are connected to these 'Set' nodes. of one or more set
  names, which are defined here as the names of the nodes with
  `nodeType="Set"` to which each node with `nodeType="Gene"` is
  connected. For example `set_nodes=c("A", "B")` will return all Gene
  nodes that are connected only to Set nodes `"A"` and `"B"`. If
  `set_nodes=NULL` then this function returns a `list` with all
  set_nodes observed.

- sep:

  `character` string used as a delimited when combining `set_nodes` into
  a name.

- filter_set_only:

  `logical` default TRUE, whether to exclude Set nodes, with vertex
  attribute nodeType='Set', from being members of nodegroups. In
  principle, 'Gene' nodes are typically members of Cnet nodesets, and
  'Set' nodes do not usually share the same 'Gene' nodes with other
  'Set' nodes, and therefore are not included.

- ...:

  additional arguments are ignored.

## Value

`list` named by nodeset, with `character` vectors of node names in each
nodeset. When `filter_set_only=TRUE` by default, the nodes will not
contain 'Set' nodes. The `list` is named by nodes connected by nodes in
each node group, using `sep` as delimiter.

## Details

This function operates on a Cnet `igraph` object, distinguished by node
attribute `'nodeType'` with values: `'Gene'` for Gene nodes, and `'Set'`
for Set nodes. This function returns Gene nodes that are connected to
the Sets given by argument `set_nodes`. It is useful to identify which
genes are connected to Set `"A"` and `"B"` for example, in other words
the Gene nodes in a specific subcluster of the Cnet `igraph`.

When `set_nodes` is `NULL`, this function returns a `list` containing
`character` vectors, where each vector represents Gene nodes that are
connected to each combination of Set nodes. This option is useful for
obtaining all possible Gene subclusters.

This function is also called by
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
in order to manipulate all nodes in a subcluster as a group, for example
calling
[`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)
on the whole set of nodes.

Note that when input data do not contain vertex attribute 'nodeType',
data are still returned but cannot be filtered by 'Set' and 'Gene', and
so will contain all possible node-connection groupings.

## See also

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

## Examples

``` r
cnet1 <- make_cnet_test()
get_cnet_nodeset(cnet1)
#> $SetA
#>  [1] "AB" "AO" "AV" "AW" "BE" "F"  "H"  "J"  "K"  "M"  "O"  "Q"  "R"  "W"  "X" 
#> [16] "Y" 
#> 
#> $`SetA,SetB`
#> [1] "BH" "BT" "BU"
#> 
#> $`SetA,SetB,SetC,SetD`
#> [1] "CL" "CM" "CN" "CO" "CP"
#> 
#> $`SetA,SetB,SetD`
#> [1] "BZ" "CD" "CI"
#> 
#> $`SetA,SetC`
#> [1] "BR" "BS" "BY"
#> 
#> $`SetA,SetC,SetD`
#> [1] "CB" "CC" "CE" "CF" "CG" "CJ" "CK"
#> 
#> $`SetA,SetD`
#> [1] "BK" "BL" "BN" "BQ"
#> 
#> $SetB
#>  [1] "AK" "AL" "AP" "AR" "AS" "AT" "BC" "C"  "D"  "E"  "I"  "P"  "V"  "Z" 
#> 
#> $`SetB,SetC`
#> [1] "BG" "BM" "BO" "BV"
#> 
#> $`SetB,SetC,SetD`
#> [1] "CA" "CH"
#> 
#> $`SetB,SetD`
#> [1] "BF" "BJ" "BP" "BW" "BX"
#> 
#> $SetC
#>  [1] "A"  "AC" "AD" "AE" "AH" "AI" "AN" "AU" "AX" "AY" "AZ" "BA" "BB" "G"  "L" 
#> [16] "N"  "T"  "U" 
#> 
#> $`SetC,SetD`
#> [1] "BI"
#> 
#> $SetD
#> [1] "AA" "AF" "AG" "AJ" "AM" "AQ" "B"  "BD" "S" 
#> 

get_cnet_nodeset(cnet1, c("SetB,SetD"))
#> $`SetB,SetD`
#> [1] "BF" "BJ" "BP" "BW" "BX"
#> 
```
