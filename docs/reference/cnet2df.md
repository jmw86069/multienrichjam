# Summarize Cnet igraph as a data.frame

Summarize Cnet igraph as a data.frame

## Usage

``` r
cnet2df(g, getNeighbors = TRUE, checkSubsets = getNeighbors, ...)
```

## Arguments

- g:

  igraph object containing Cnet data, specifically vertex attribute name
  "nodeType" with values "Set" and "Gene", and where "Set" nodes are
  only connected to "Gene" nodes.

- getNeighbors:

  logical indicating whether to include the connected neighbor node
  names.

- checkSubsets:

  logical indicating whether to test "Set" nodes to determine if the
  neighbors are all represented by another "Set" node.

- ...:

  additional arguments are ignored.

## Value

data.frame with the node name, label, degree (number of edges),
membership (based upon connected component), and if `getNeighbors=TRUE`
it includes comma-delimited names of neighboring nodes.

## Details

This function provides a data.frame summary of an igraph object
containing "Cnet" data, including vertex attribute `"nodeType"` with
values `"Set"` and `"Gene"`, and where `"Set"` nodes are only connected
to `"Gene"` nodes.

The data.frame is intended to provide a convenient method for subsetting
nodes, typically based upon a connected cluster, or the minimum number
of edges per node. For example, filter for the connected component
containing a node of interest, or filter for `"Set"` nodes with more
than one `"Gene"`.

## See also

Other jam conversion functions:
[`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md),
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md),
[`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md),
[`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
