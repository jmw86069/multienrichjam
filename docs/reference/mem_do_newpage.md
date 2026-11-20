# Determine whether Mem should call newpage

Determine whether Mem should call newpage

## Usage

``` r
mem_do_newpage(x, ...)
```

## Arguments

- x:

  not currently implemented, it may be used in future.

- ...:

  additional arguments are ignored.

## Details

### Rules

- Running outside Positron, then yes. End.

- If running inside knitr, yes. End.

- Otherwise no.
