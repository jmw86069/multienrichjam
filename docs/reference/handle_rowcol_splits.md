# Handle row and column split parameters for gene-pathway data

Handle row and column split parameters for gene-pathway data

## Usage

``` r
handle_rowcol_splits(
  Mem,
  auto_split = TRUE,
  row_split,
  row_title,
  max_row_split = 12,
  cluster_rows,
  row_method,
  gene_im_weight = 0.5,
  column_split,
  column_title,
  max_column_split = 8,
  cluster_columns,
  column_method,
  enrich_im_weight = 0.3,
  trim_rows = TRUE,
  trim_columns = TRUE,
  p_cutoff = NULL,
  p_floor = NULL,
  seed = 123,
  verbose = FALSE,
  debug = FALSE,
  ...
)
```

## Value

`list` with

- Mem

- row_split

- row_title

- cluster_rows

- column_split

- column_title

- cluster_columns
