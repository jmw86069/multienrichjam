# Determine whether Mem should call newpage

Determine whether Mem should call newpage

## Usage

``` r
mem_do_newpage(x, ...)
```

## Details

### Rules

- Running outside Positron, then yes. End.

- If running inside knitr, yes. End.

- Otherwise no.
