# Launch ShinyCAT: Cnet Adjustment Tool

Launch ShinyCAT: Cnet Adjustment Tool

## Usage

``` r
launch_shinycat(g = NULL, envir = new.env(), ..., options = list(width = 1200))
```

## Arguments

- g:

  `igraph` object, or can be NULL if the variable 'g' is defined in
  'envir'.

- envir:

  `environment` assigned to the R-shiny function space, which may be
  useful when capturing the `igraph` object after manipulation via the
  R-shiny app.

- ...:

  additional arguments are ignored.

- options:

  `list` with additional settings, for example:

  - host: `character` string with hostname or IP address to permit
    incoming connections. Use '0.0.0.0' to accept all incoming host
    names.

  - port: `numerical` port to listen for incoming connections.

  - quiet: `logical` whether to

  - launch.browser: `logical` whether to launch a web browser, default
    TRUE for interactive sessions.

  - display.mode: `character` string with display mode: 'normal' is
    recommended; 'showcase' displays the code beside the app.

## Value

`environment` invisibly, containing

- 'g' the input `igraph`

- 'adj_cnet' the output `igraph` after manipulation in the R-shiny app.
  It also contains two attributes for reference:

  1.  'nodeset_adj': `data.frame` used to adjust nodesets.

  2.  'node_adj': `data.frame` used to adjust nodes.

## Details

This function launches the R-shiny app 'ShinyCat' to manipulate a Cnet
`igraph` object, which is expected to contain node (vertex) attribute
'nodeType' with values 'Gene' and 'Set'.

There are two ways to retain the resulting igraph object:

1.  Capture the output of this function, an `environment` described
    below.

2.  Click "Adjustments" then "Save RData", which stores the 'igraph\`
    object as 'adj_cnet'.

The function returns an `environment` inside which are two objects:

- 'g': the `igraph` data input

- 'adj_cnet': the adjusted Cnet `igraph` output

Additionally, the 'adj_cnet' object contains two attributes:

- 'nodeset_adj', 'node_adj': `data.frame` objects used by
  [`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md)
  to convert 'g' into 'adj_cnet'.

## See also

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

Other jam shiny functions:
[`shinycat_server()`](https://jmw86069.github.io/multienrichjam/reference/shinycat_server.md),
[`shinycat_ui()`](https://jmw86069.github.io/multienrichjam/reference/shinycat_ui.md)

## Examples

``` r
# create Cnet test data
g <- make_cnet_test();
cnetenv <- new.env();

# you must catch the output to use the resulting igraph object
# output_envir <- launch_shinycat();
```
