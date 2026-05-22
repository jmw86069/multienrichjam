# convert incidence matrix to list

convert incidence matrix to list

## Usage

``` r
im2list_dep(x, verbose = FALSE, ...)
```

## Value

`list` of `character vectors`, where list names are defined by
`colnames(x)`, and vectors contain values from `rownames(x)`.

## Details

This function converts an incidence `matrix`, or equivalent
`data.frame`, to a list. The `matrix` should contain either numeric
values such as `c(0, 1)`, or logical values such as `c(TRUE,FALSE)`,
otherwise values are considered either zero == `FALSE`, or non-zero ==
`TRUE`.

The resulting list will be named by `colnames(x)` of the input, and will
contain members named by `rownames(x)` which are either non-zero, or
contain `TRUE`.

Values of `NA` are converted to zero `0` and therefore ignored.

This function uses the `transactions` class from the `arules` R package,
which in our testing is substantially faster than similar techniques
from a variety of other R packages.

## See also

Other jam list functions:
[`colors_from_list`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md)`()`,
[`im2list`](https://jmw86069.github.io/multienrichjam/reference/im2list.md)`()`,
[`imSigned2list`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md)`()`,
[`imSigned2list_dep`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list_dep.md)`()`,
[`list2concordance`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md)`()`,
[`list2im`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)`()`,
[`list2imSigned`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)`()`
