# Check Mem object

Check whether a Mem object is valid.

## Usage

``` r
check_Mem(object)
```

## Arguments

- object:

  `Mem` object

## Details

It requires:

- `slotName(object)`: `"memIM", "enrichIM", "geneIM", "enrichList"`

- `object@memIM` to inherit `"matrix"`

- `object@enrichIM` to inherit `"matrix"`

- `object@geneIM` to inherit `"matrix"`

- `colnames(object@memIM)` to equal `rownames(object@enrichIM)`

- `rownames(object@memIM)` to equal `colnames(object@geneIM)`

- `names(enrichList)` to equal `colnames(object@enrichIM)`

General guidance for Mem objects:

- The aim is to support all user-facing tasks using S4 methods, to
  prevent the need to access slots directly.

## See also

Other Mem:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`Mem-data`](https://jmw86069.github.io/multienrichjam/reference/Mem-data.md),
[`Mem-slots`](https://jmw86069.github.io/multienrichjam/reference/Mem-slots.md)
