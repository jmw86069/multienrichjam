
## Brief steps to produce a Cnet-Heatmap plot

1. Run multiEnrichMap()

   * Apply appropriate p_cutoff, and topEnrichN to reduce the number
   of pathways per enrichment result, topEnrichN from 10 to 20 for
   two enrichment results; topEnrichN 4 to 10 for four or more enrichment
   results.

2. Run mem_plot_folio()

```{r}
cairo_pdf(file="folio_mem_DM_AffySoma_top20_23apr2020.pdf",
   height=12, width=12, pointsize=11, onefile=TRUE);
mem_dm_plots <- mem_plot_folio(mem_dm_affysoma,
   do_which=2:5,
   pathway_column_split=4,
   main="Dermatomyositis Affy-SomaLogic");
dev.off();
```

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

```{r}
ht_opt("legend_border"="black");
cnet_collapsed_set <- mem_dm_plots$cnet_collapsed_set;
E(cnet_collapsed_set)$color <- "#BBBBBB99";
V(cnet_collapsed_set)$frame.color <- "#77777799";
isset <- V(cnet_collapsed_set)$nodeType %in% "Set";
V(cnet_collapsed_set)$label <- ifelse(isset,
   V(cnet_collapsed_set)$name,
   V(cnet_collapsed_set)$label);
V(cnet_collapsed_set)$label.family <- "Arial";
V(cnet_collapsed_set)$label.color <- ifelse(isset,
   "grey85",
   "grey10");
## Manually add distance for label from node
V(cnet_collapsed_set)$label.dist <- 0.8;
V(cnet_collapsed_set)$label.dist[igrep("^PR", V(cnet_collapsed_set)$name)] <- 2;
```

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

```{r}
{cairo_pdf(file="DM_GeneProtein_MultiEnrich_Heatmap_23apr2020.pdf",
   height=9, width=10, pointsize=12);
plot_cnet_heatmaps(mem_folio=cnet_collapsed_set,
   iexprs_ctr=all_pergene_exprs_ctr,
   col=circlize::colorRamp2(breaks=seq(from=-2, to=2, length.out=15),
      colors=getColorRamp("RdBu_r", n=15, lens=2)),
   isamples=all_samples_dm,
   column_title_gp=gpar(fontsize=10),
   igroups=gsub(" ", "\n", all_cohort_type[all_samples_dm]),
   allowed_degree=c(1,2),
   use_shadowText=TRUE,
   layout_widths=c(1, 4, 1),
   cnet_mar=c(0.5, 7, 0, 6.3),
   use_gridBase=FALSE,
   border=FALSE,
   title_prefix=NULL,
   label_factor_l=list(nodeType=c(Gene=0.63, Set=1.8)),
   label_dist_factor_l=list(nodeType=c(Gene=1, Set=0)),
   set_names=c("C","D","A","B"),
   set_panel_row=list(c(1:2), c(1:4), c(5:8), c(6:8)),
   vjust=c(0.9, 0.9, -1, 0),
   strwrap_width=c(36, 64, 46, 46),
   fontsize_set_names=9,
   colorV=c(`DM\nGenes`="firebrick", `DM\nProteins`="gold"),
   ipa_xref_df=ipa_xref_df)
dev.off();}
```

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

6. Save the same data to Excel:


```{r}
dmlupus_df <- find_shared_pathways(ipa_pathways_all,
   filter_dysreg=FALSE,
   hm=mem_dmlupus_plots$gp_hm,
   prefix1="Affy_DM",
   prefix2="Affy_Lupus");
jamba::writeOpenxlsx("DM_Lupus_MultiEnrich_12mar2020.xlsx",
   sheetName="DM_Lupus",
   x=dmlupus_df,
   colorSub=c(Affy_DM="tomato3", Soma_DM="gold", Affy_Lupus="blue3", Soma_Lupus="cyan3", cluster_colors),
   highlightColumns=igrep("Cluster|Name", colnames(dmlupus_df)),
   pvalueColumns=igrep("p.value", colnames(dmlupus_df)));
set_xlsx_colwidths("DM_Lupus_MultiEnrich_12mar2020.xlsx",
   widths=c(11,8,11,11, 30, 11,11,11,11,11,30,15,40))
set_xlsx_rowheights("DM_Lupus_MultiEnrich_12mar2020.xlsx",
   heights=rep(56, nrow(dmlupus_df)))
```
