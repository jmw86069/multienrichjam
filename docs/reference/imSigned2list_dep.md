# convert signed incidence matrix to list

convert signed incidence matrix to list

## Usage

``` r
imSigned2list_dep(x, verbose = FALSE, ...)
```

## Value

`list` of named numeric vectors, where list names are defined by
`colnames(x)`, and vector names are derived from `rownames(x)`. Values
in each vector indicate the signed direction, `c(-1,1)`.

## Details

This function converts an signed incidence `matrix` that contains
positive and negative values, or equivalent `data.frame`, to a list of
named vectors containing values `c(-1, 1)` to indicate signed direction.
The input `matrix` should contain numeric values where positive and
negative values indicate directionality. When the input contains only
logical values `c(TRUE,FALSE)` the direction is assumed to be `+1`
positive.

Values of `NA` are converted to zero `0` and therefore ignored.

This function uses the `transactions` class from the `arules` R package,
which in our testing is substantially faster than similar techniques
from a variety of other R packages.

## See also

Other jam list functions:
[`colors_from_list`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)`()`,
[`im2list`](https://jmw86069.github.io/multienrichjam/reference/im2list.md)`()`,
[`im2list_dep`](https://jmw86069.github.io/multienrichjam/reference/im2list_dep.md)`()`,
[`imSigned2list`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md)`()`,
[`list2concordance`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md)`()`,
[`list2im`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)`()`,
[`list2imSigned`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)`()`
