# Get Cnet nodesets as a named vector

Get Cnet nodesets as a named vector, with nodegroups named by node

## Usage

``` r
get_cnet_nodeset_vector(g, nodegroups = NULL, filter_set_only = TRUE, ...)
```

## Arguments

- g:

  `igraph` Cnet object

- nodegroups:

  default is NULL, which calls
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md).
  Also accepts `list` or `communities`.

- filter_set_only:

  `logical` default TRUE, passed to
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
  to determine whether 'Set' nodes are assigned a nodeset.

  - When TRUE, 'Set' nodes are filtered out, therefore their nodeset
    becomes NA.

  - When FALSE, 'Set' nodes are also assigned to a nodeset.

- ...:

  additional arguments passed to
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)

## Value

`character` vector named by node, whose value is the name of the
nodegroup. All nodes in `g` are returned in order, with `NA` values used
when there is no nodegroup assigned to a particular node.

## Details

This function is a simple conversion of `list` output from
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
that returns a `character` vector alternative.

It is easier to answer: "What nodegroup is this node in?"

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

## Examples

``` r
cnet <- make_cnet_test()
nsv <- get_cnet_nodeset_vector(cnet)
nsv
#>                     A                     B                     C 
#>                "SetC"                "SetD"                "SetB" 
#>                     D                     E                     F 
#>                "SetB"                "SetB"                "SetA" 
#>                     G                     H                     I 
#>                "SetC"                "SetA"                "SetB" 
#>                     J                     K                     L 
#>                "SetA"                "SetA"                "SetC" 
#>                     M                     N                     O 
#>                "SetA"                "SetC"                "SetA" 
#>                     P                     Q                     R 
#>                "SetB"                "SetA"                "SetA" 
#>                     S                     T                     U 
#>                "SetD"                "SetC"                "SetC" 
#>                     V                     W                     X 
#>                "SetB"                "SetA"                "SetA" 
#>                     Y                     Z                    AA 
#>                "SetA"                "SetB"                "SetD" 
#>                    AB                    AC                    AD 
#>                "SetA"                "SetC"                "SetC" 
#>                    AE                    AF                    AG 
#>                "SetC"                "SetD"                "SetD" 
#>                    AH                    AI                    AJ 
#>                "SetC"                "SetC"                "SetD" 
#>                    AK                    AL                    AM 
#>                "SetB"                "SetB"                "SetD" 
#>                    AN                    AO                    AP 
#>                "SetC"                "SetA"                "SetB" 
#>                    AQ                    AR                    AS 
#>                "SetD"                "SetB"                "SetB" 
#>                    AT                    AU                    AV 
#>                "SetB"                "SetC"                "SetA" 
#>                    AW                    AX                    AY 
#>                "SetA"                "SetC"                "SetC" 
#>                    AZ                    BA                    BB 
#>                "SetC"                "SetC"                "SetC" 
#>                    BC                    BD                    BE 
#>                "SetB"                "SetD"                "SetA" 
#>                    BF                    BG                    BH 
#>           "SetB,SetD"           "SetB,SetC"           "SetA,SetB" 
#>                    BI                    BJ                    BK 
#>           "SetC,SetD"           "SetB,SetD"           "SetA,SetD" 
#>                    BL                    BM                    BN 
#>           "SetA,SetD"           "SetB,SetC"           "SetA,SetD" 
#>                    BO                    BP                    BQ 
#>           "SetB,SetC"           "SetB,SetD"           "SetA,SetD" 
#>                    BR                    BS                    BT 
#>           "SetA,SetC"           "SetA,SetC"           "SetA,SetB" 
#>                    BU                    BV                    BW 
#>           "SetA,SetB"           "SetB,SetC"           "SetB,SetD" 
#>                    BX                    BY                    BZ 
#>           "SetB,SetD"           "SetA,SetC"      "SetA,SetB,SetD" 
#>                    CA                    CB                    CC 
#>      "SetB,SetC,SetD"      "SetA,SetC,SetD"      "SetA,SetC,SetD" 
#>                    CD                    CE                    CF 
#>      "SetA,SetB,SetD"      "SetA,SetC,SetD"      "SetA,SetC,SetD" 
#>                    CG                    CH                    CI 
#>      "SetA,SetC,SetD"      "SetB,SetC,SetD"      "SetA,SetB,SetD" 
#>                    CJ                    CK                    CL 
#>      "SetA,SetC,SetD"      "SetA,SetC,SetD" "SetA,SetB,SetC,SetD" 
#>                    CM                    CN                    CO 
#> "SetA,SetB,SetC,SetD" "SetA,SetB,SetC,SetD" "SetA,SetB,SetC,SetD" 
#>                    CP                  SetA                  SetB 
#> "SetA,SetB,SetC,SetD"                    NA                    NA 
#>                  SetC                  SetD 
#>                    NA                    NA 
```
