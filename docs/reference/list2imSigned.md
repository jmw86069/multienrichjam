# convert list to directional incidence matrix

convert list to directional incidence matrix

## Usage

``` r
list2imSigned(x, emptyValue = NA, verbose = FALSE, ...)
```

## Arguments

- x:

  `list` of named vectors, where the names are used to identify each
  element, and become rownames in the output incidence matrix. The
  vector values become values in the incidence matrix.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md).

## Value

`numeric` matrix with rownames defined by vector names from each vector
in the input list. The colnames are defined by names of the input list
if they exist. list. Values in the matrix are values from each vector.

## Details

This function extends
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)
in that it stores the value associated with each element in the list. As
such, the input format is a named vector, where the names of the vector
are the items, and the numeric values are the values to be stored in the
incidence matrix.

A common scenario is to generate a vector of genes, with values
`c(-1, 0, 1)` indicating the direction of gene expression changes, named
by the gene symbol. Each vector in the list represents one statistical
test. Here, `list2imSigned()` will convert this list into a directional
matrix representing the gene changes across the comparisons.

Note that this function currently does not combine multiple values,
instead only the last occurring value is stored in the resulting matrix.
This decision is partly due to efficiency, and partly because there are
multiple possible methods to combine multiple values. For example,
taking the `mean(x)` for a given gene, which has a value `1` and `-1`
would result in `0` and might suggest the gene is not a statistical hit.
Instead, when multiple values are anticipated per named vector entry,
use functions in a package like `data.table` or `dplyr` to apply a
function to combine values.

## See also

Other jam list functions:
[`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md),
[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md),
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md),
[`list2concordance()`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md),
[`list2im()`](https://jmw86069.github.io/multienrichjam/reference/list2im.md)

## Examples

``` r
L1 <- list(A=c("C","A","B","A"),
   D=c("D","E","F","D"),
   A123=c(1:8,3,5),
   T=LETTERS[7:9]);
L1;
#> $A
#> [1] "C" "A" "B" "A"
#> 
#> $D
#> [1] "D" "E" "F" "D"
#> 
#> $A123
#>  [1] 1 2 3 4 5 6 7 8 3 5
#> 
#> $T
#> [1] "G" "H" "I"
#> 

# Convert each vector to a signed vector
set.seed(123);
L2 <- lapply(L1, function(i){
   i <- unique(i);
   jamba::nameVector(sample(c(-1,1), size=length(i), replace=TRUE), i);
});
L2;
#> $A
#>  C  A  B 
#> -1 -1 -1 
#> 
#> $D
#>  D  E  F 
#>  1 -1  1 
#> 
#> $A123
#>  1  2  3  4  5  6  7  8 
#>  1  1 -1 -1  1  1  1 -1 
#> 
#> $T
#>  G  H  I 
#>  1 -1  1 
#> 

# Convert to signed incidence matrix
list2imSigned(L2);
#>    A  D A123  T
#> C -1 NA   NA NA
#> A -1 NA   NA NA
#> B -1 NA   NA NA
#> D NA  1   NA NA
#> E NA -1   NA NA
#> F NA  1   NA NA
#> 1 NA NA    1 NA
#> 2 NA NA    1 NA
#> 3 NA NA   -1 NA
#> 4 NA NA   -1 NA
#> 5 NA NA    1 NA
#> 6 NA NA    1 NA
#> 7 NA NA    1 NA
#> 8 NA NA   -1 NA
#> G NA NA   NA  1
#> H NA NA   NA -1
#> I NA NA   NA  1
```
