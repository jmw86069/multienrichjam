# convert signed incidence matrix to list

convert signed incidence matrix to list

## Usage

``` r
imSigned2list(x, empty = c(NA, "", 0, FALSE), ...)
```

## Arguments

- x:

  `matrix` or equivalent object with `colnames(x)` indicating list set
  names, and `rownames(x)` indicating list contents.

- empty:

  `character` vector of incidence matrix values that should be
  considered "empty" and therefore do not indicate the row in `x` is
  present for the given column in `x`. All other items are considered to
  be present, and are assigned direction based upon the value in that
  cell of `x`.

- ...:

  additional arguments are ignored.

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

Values that are `logical` with `TRUE` and `FALSE` are converted to
`numeric` before output.

## See also

Other jam list functions:
[`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md),
[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md),
[`list2concordance()`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md),
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md),
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)

## Examples

``` r
im <- matrix(c(0,1,-1,1,1,NA,-1,0,1),
   ncol=3,
   nrow=3,
   dimnames=list(letters[1:3], LETTERS[1:3]))
print(im);
#>    A  B  C
#> a  0  1 -1
#> b  1  1  0
#> c -1 NA  1
# matrix input
im2list(im);
#> $A
#> [1] "b" "c"
#> 
#> $B
#> [1] "a" "b"
#> 
#> $C
#> [1] "a" "c"
#> 
imSigned2list(im);
#> $A
#>  b  c 
#>  1 -1 
#> 
#> $B
#> a b 
#> 1 1 
#> 
#> $C
#>  a  c 
#> -1  1 
#> 
imSigned2list(im != 0);
#> $A
#> b c 
#> 1 1 
#> 
#> $B
#> a b 
#> 1 1 
#> 
#> $C
#> a c 
#> 1 1 
#> 
```
