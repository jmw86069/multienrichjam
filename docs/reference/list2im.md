# convert list to incidence matrix

convert list to incidence matrix

## Usage

``` r
list2im(
  x,
  empty = 0,
  do_sparse = FALSE,
  keepCounts = FALSE,
  sort_rows = FALSE,
  ...
)
```

## Arguments

- x:

  `list` of item vectors, where vector values are item names

- empty:

  `numeric` used for missing or empty values, default 0.

- do_sparse:

  `logical` default FALSE, whether to return sparse Matrix. When
  `do_sparse=TRUE`, it sets `keepCounts=FALSE` since a sparse Matrix
  does not represent counts.

- keepCounts:

  `logical` default FALSE, whether to include the count of any
  duplicated items.

- sort_rows:

  `logical` default FALSE, whether to sort rows using
  [`jamba::mixedOrder()`](https://jmw86069.github.io/jamba/reference/mixedOrder.html),
  or `function` can be supplied to use another sort algorithm.

- ...:

  additional arguments are ignored.

## Value

`numeric` matrix with values 0 or 1 indicating absense or presence of
each item (row) in each set (column). The sets are derived from
`names(x)`.

## Details

This function converts a list of vectors into an incidence matrix, where
the rows are the vector items and the columns are the list names.

## See also

Other jam list functions:
[`colors_from_list()`](https://jmw86069.github.io/multienrichjam/reference/colors_from_list.md),
[`im2list()`](https://jmw86069.github.io/multienrichjam/reference/im2list.md),
[`imSigned2list()`](https://jmw86069.github.io/multienrichjam/reference/imSigned2list.md),
[`list2concordance()`](https://jmw86069.github.io/multienrichjam/reference/list2concordance.md),
[`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)

## Examples

``` r
L1 <- list(A=c("C","A","B","A"),
   D=c("D","E","F","D"),
   A123=c(1:8,3,5),
   B=c("A1", "A2", "A10"),
   T=LETTERS[7:9]);
# Default behavior is to make items unique
list2im(L1);
#>     A D A123 B T
#> C   1 0    0 0 0
#> A   1 0    0 0 0
#> B   1 0    0 0 0
#> D   0 1    0 0 0
#> E   0 1    0 0 0
#> F   0 1    0 0 0
#> 1   0 0    1 0 0
#> 2   0 0    1 0 0
#> 3   0 0    1 0 0
#> 4   0 0    1 0 0
#> 5   0 0    1 0 0
#> 6   0 0    1 0 0
#> 7   0 0    1 0 0
#> 8   0 0    1 0 0
#> A1  0 0    0 1 0
#> A2  0 0    0 1 0
#> A10 0 0    0 1 0
#> G   0 0    0 0 1
#> H   0 0    0 0 1
#> I   0 0    0 0 1

# Option to report the counts
list2im(L1, keepCounts=TRUE);
#>     A D A123 B T
#> C   1 0    0 0 0
#> A   2 0    0 0 0
#> B   1 0    0 0 0
#> D   0 2    0 0 0
#> E   0 1    0 0 0
#> F   0 1    0 0 0
#> 1   0 0    1 0 0
#> 2   0 0    1 0 0
#> 3   0 0    2 0 0
#> 4   0 0    1 0 0
#> 5   0 0    2 0 0
#> 6   0 0    1 0 0
#> 7   0 0    1 0 0
#> 8   0 0    1 0 0
#> A1  0 0    0 1 0
#> A2  0 0    0 1 0
#> A10 0 0    0 1 0
#> G   0 0    0 0 1
#> H   0 0    0 0 1
#> I   0 0    0 0 1

# Option to sort item rows with jamba::mixedOrder()
list2im(L1, sort_rows=TRUE);
#>     A D A123 B T
#> 1   0 0    1 0 0
#> 2   0 0    1 0 0
#> 3   0 0    1 0 0
#> 4   0 0    1 0 0
#> 5   0 0    1 0 0
#> 6   0 0    1 0 0
#> 7   0 0    1 0 0
#> 8   0 0    1 0 0
#> A1  0 0    0 1 0
#> A2  0 0    0 1 0
#> A10 0 0    0 1 0
#> A   1 0    0 0 0
#> B   1 0    0 0 0
#> C   1 0    0 0 0
#> D   0 1    0 0 0
#> E   0 1    0 0 0
#> F   0 1    0 0 0
#> G   0 0    0 0 1
#> H   0 0    0 0 1
#> I   0 0    0 0 1

# Option to sort item rows with vanilla sort()
list2im(L1, sort_rows=sort);
#>     A D A123 B T
#> 1   0 0    1 0 0
#> 2   0 0    1 0 0
#> 3   0 0    1 0 0
#> 4   0 0    1 0 0
#> 5   0 0    1 0 0
#> 6   0 0    1 0 0
#> 7   0 0    1 0 0
#> 8   0 0    1 0 0
#> A   1 0    0 0 0
#> A1  0 0    0 1 0
#> A10 0 0    0 1 0
#> A2  0 0    0 1 0
#> B   1 0    0 0 0
#> C   1 0    0 0 0
#> D   0 1    0 0 0
#> E   0 1    0 0 0
#> F   0 1    0 0 0
#> G   0 0    0 0 1
#> H   0 0    0 0 1
#> I   0 0    0 0 1
```
