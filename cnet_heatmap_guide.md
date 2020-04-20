
## Brief steps to produce a Cnet-Heatmap plot

1. Run multiEnrichMap()

   * Apply appropriate p_cutoff, and topEnrichN to reduce the number
   of pathways per enrichment result, topEnrichN from 10 to 20 for
   two enrichment results; topEnrichN 4 to 10 for four or more enrichment
   results.

2. Run mem_plot_folio()

   * Make sure `do_which` includes `c(2,4)`, which produces cnet_collapsed_set
   * Review the gene-protein heatmap to ensure the clusters are
   representative of the underlying content of the data. Adjust clustering
   parameters as needed, for example:
   
      * `column_method="correlation"` changes column clustering to use
      correlation, instead of euclidean distance, which can emphasize the
      top annotation enrichment colors, slightly more than the gene-protein matrix.
      * `min_gene_ct=1`, `min_sample_ct=1` filters pathways and genes to
      require more than one gene per pathway, or more than one pathway
      per gene. Sometimes this filtering helps the clustering focus on
      genes/pathways which are shared rather than random singlet genes.
      * `pathway_column_split=4` to adjust the number of pathway clusters,
      can be helpful if there are more or fewer clusters visible in the
      gene-pathway heatmap than are being produced by the automatic process.
      For large number of pathways, sometimes more clusters are necessary. 
      * `p_floor` especially when `p_cutoff` is lower than 0.01, the `p_floor`
      should be set appropriately to represent the range of enrichment P-values
      in the top annotation color bar. Note that it may not be useful to 
      show P-values of 1e-100 -- anything less than 1e-10 should in theory
      be equivalent.
   * Note the location of Cnet clusters, `c("A", "B", "C", "D")` because
   they need to be defined manually in the Cnet-Heatmap.

3. Extract `cnet_collapsed_set` as the Cnet plot to use for the Cnet-Heatmap.

   * Optionally adjust the network layout using some helpful functions:
   
      * `rotate_igraph_layout()` - rotate the layout so the clusters can
      be aligned with heatmaps on the outer edges.
      * `spread_igraph_labels()` - spreads out igraph labels based upon
      the angle of incoming edges, it does pretty well at reducing
      overlapping labels. Use `nodeSortBy` to change color-sorting from
      "left-to-right" with `nodeSortBy="x"`, to "top-to-bottom" with
      `nodeSortBy="-y"`. Use `label_min_dist` to increase the distance
      of labels from the node, where `0.5` is roughly half the node diameter,
      which is of course 1x radius. Use `force_relayout` and `repulse` to
      force a new layout with different repulse adjustment, higher repulse
      makes nodes more tightly grouped.

4. Prepare centered expression data, whose rownames are genes that are
present in the gene-pathway data above.

   * Typically data should be log2-transformed, and row-centered,
   which can be done using `jamma::centerGeneData()`.
   
      * Center each group by its associated control group, using
      `centerGroups` and `controlSamples` arguments. It is not
      recommended to scale the data, but that is available with
      `scale="row"`.
      * Optionally omit control groups from the heatmaps to maximize
      plot space. You may want to create a plot with control groups
      to verify that the control groups have low variability, then
      re-create the plot without control groups for the final figure.

   * Define `isamples` using `colnames()` of the expression matrix.
   * Define `igroups` to be a vector of group labels, named by
   `colnames()` of the expression matrix.
   * Define `colorV` whose names are unique values from `igroups`.

      * These `colorV` values are shown at the top of each heatmap
      as annotations to represent the groups in `igroups`. The
      names are displayed below the heatmap, should they should be
      small enough to fit in a limited space. Use `legend_labels`
      below to use a different label for the color legend.
      * Ideally these colors should also match the colors from
      `colorV` produced by `multiEnrichMap()`, but it is not required.
      * Note that when the color is matched between the sample group
      and the gene node, a box is drawn around the heatmap panel to
      indicate these genes were used in the enrichment test.

5. Run `plot_cnet_heatmap()`.

   * There are numerous arguments to this function, intended to
   help configure the many details of the figure.
   * `iexprs_ctr` - the expression data matrix, intended to show centered
   data, but in fact any numeric data is accepted.
   * `allowed_degree=c(1,2)` uses any gene node with 1 or 2 connections
   for each cluster heatmap. When `allowed_degree=1` the heatmap will
   include genes that connect *only* to that pathway cluster, and no
   other pathway clusters. Sometimes it is useful to focus on specific
   genes, to make sure there are no changes visible in the heatmap across
   the other groups. If there are changes across other groups, the
   specificity may be an artifact of choosing too stringent a statistical
   threshold for gene hits prior to running the enrichment test.
   * `legend_labels` are used to provide a custom label for the colors,
   helpful when `colorV` has abbreviated labels.
   * `cnet_mar` - adjust margins `c(bottom, left, top, right)` 
   around the central Cnet plot, to reduce overlaps across labels.
   * `vjust` - vertical justification of pathway labels for each cluster,
   where `1` positions labels at the top of the top layout cell,
   `0` positions labels at the top of the next layout cell lower, and `-1`
   positions labels at the top of the layout cell two cells lower.
   * `strwrap_width` is used to define the word wrap width for each
   pathway cluster. Sometimes it helps to make some sections wider or
   narrower than others, to reduce label overlaps.
   * `set_names` vector of pathway clusters to use, in order. These names
   should reflect the cluster names such as `c("A", "B", "C", "D")`, but
   should be in order they appear in the Cnet layout.
   * Define `set_panel_row` and `set_panel_col` as the range of rows and
   columns for each heatmap. They should be positioned relatively close
   to the Cnet pathway clusters, with height relative to the number of
   genes.
   * `use_shadowText=TRUE` uses a slight outline effect around node labels,
   which makes it easier to read over light/dark background colors.
