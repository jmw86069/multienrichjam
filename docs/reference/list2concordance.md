# Convert list to concordance matrix

Convert list to concordance matrix

## Usage

``` r
list2concordance(x, naValue = NA, makeSigned = TRUE, verbose = FALSE, ...)
```

## Arguments

- x:

  `list` of named numerical vectors, where the sign (positive or
  negative sign) indicates directionality, and is used to calculate
  concordance, which is a measure of the agreement of the overall set of
  directions shared between each pair of vectors.

- naValue:

  value passed to
  [`jamba::rmNA()`](https://jmw86069.github.io/jamba/reference/rmNA.html)
  used to replace any `NaN` values in the output matrix. The `NaN`
  values result when a pair of vectors has no shared non-zero named
  entry.

- makeSigned:

  logical indicating whether to force the vectors in the input list `x`
  to contain only values `c(-1,0,1)`, by calling
  [`base::sign()`](https://rdrr.io/r/base/sign.html).

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md).

## Details

This function calculates pairwise concordance using Kruskal concordance
coefficient (ref) using the following equation:

- (number_agree - number_disagree) / (total_shared)

The equation is applied to each pair of named vectors in the input list
`x`, and reflects the degree of agreement in direction (+ or -) between
shared named elements, with +1 being perfect concordance (agreement),
and -1 being perfect discordance (disagreement.) Values of zero indicate
equal agreement and disagreement, and therefore reflect no concordance
nor discordance. Values of `NA` occur when no named entries are shared.

This function calls
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)
to produce a signed incidence matrix, which is then used with
[`base::crossprod()`](https://rdrr.io/r/base/crossprod.html) to
calculate the full matrix of values.

## See also

Other jam list functions:
[`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md),
[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md),
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md),
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md),
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)

## Examples

``` r
set.seed(123);
l123 <- lapply(jamba::nameVector(1:3), function(i){
   jamba::nameVector(
      sample(c(-1,-1,0,1), replace=TRUE, size=15),
      letters[1:15]
   )
});
list2concordance(l123);
#>           1         2         3
#> 1 1.0000000 0.4285714 0.0000000
#> 2 0.4285714 1.0000000 0.3333333
#> 3 0.0000000 0.3333333 1.0000000

# observe the signed incidence matrix
list2imSigned(l123);
#>    1  2  3
#> a  0  0 -1
#> b  0  1 -1
#> c  0 -1  0
#> d -1  0  1
#> e  0  0 -1
#> f -1 -1 -1
#> g -1  1  0
#> h -1 -1  0
#> i  0 -1 -1
#> j -1 -1  1
#> k  1  0  0
#> l -1  1  1
#> m -1 -1 -1
#> n -1  0 -1
#> o -1 -1  0
```
