# convert incidence matrix to list

convert incidence matrix to list

## Usage

``` r
im2list(x, empty = c(NA, "", 0, FALSE), ...)
```

## Arguments

- x:

  `matrix` or equivalent object with `colnames(x)` indicating list set
  names, and `rownames(x)` indicating list contents.

- empty:

  `character` vector of incidence matrix values that should be
  considered "empty" and therefore do not indicate the row in `x` is
  present for the given column in `x`. All other items are considered to
  be present.

- ...:

  additional arguments are ignored.

## Value

`list` of `character vectors`, where list names are defined by
`colnames(x)`, and list elements are vectors that contain values from
`rownames(x)`.

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

## See also

Other jam list functions:
[`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md),
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md),
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

# data.frame
imdf <- data.frame(im);
print(imdf);
#>    A  B  C
#> a  0  1 -1
#> b  1  1  0
#> c -1 NA  1
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

# logical input
imtf <- (!im == 0);
print(imtf);
#>       A    B     C
#> a FALSE TRUE  TRUE
#> b  TRUE TRUE FALSE
#> c  TRUE   NA  TRUE
im2list(imtf);
#> $A
#> [1] "b" "c"
#> 
#> $B
#> [1] "a" "b"
#> 
#> $C
#> [1] "a" "c"
#> 
```
