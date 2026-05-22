# Subset Mem object

Extract some elements of potentially long vector for display

## Usage

``` r
# S4 method for class 'Mem,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

some_vector(
  x,
  maxToShow = 5,
  ellipsis = "...",
  ellipsisPos = c("middle", "start", "end"),
  quote = FALSE,
  sep = ", ",
  ...
)
```

## See also

Other Mem:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`check_Mem`](https://jmw86069.github.io/multienrichjam/reference/check_Mem.md)`()`,
[`list_to_Mem`](https://jmw86069.github.io/multienrichjam/reference/list_to_Mem.md)`()`
