
#' Make Cnet test igraph
#'
#' Make Cnet test igraph
#'
#' This function simply creates an `igraph` object with attributes
#' expected for a cnet plot object:
#' * node attribute `"nodeType"` with values `c("Gene", "Set")`.
#'
#' It optionally derives random directionality when `add_direction=TRUE`,
#' and calls `apply_cnet_direction()` so node borders are updated
#' appropriately.
#'
#' @family jam cnet utilities
#' 
#' @returns `igraph` object containing Cnet concept network data, specifically
#'    with vertex attribute 'nodeType' with values 'Gene' or 'Set'.
#' 
#' @param num_sets `integer` number of sets, default 4.
#' @param overlap_counts `integer` vector of counts, with length 'num_sets',
#'    default 'c(57, 20, 12, 5)'.
#' @param row_prefix `character` default "" used as a prefix for rows names.
#' @param column_prefix `character` default "Set" used as a prefix for
#'    set names.
#' @param add_direction `logical` default TRUE, whether to include direction.
#' @param set_colors `character` vector of colors per set.
#' @param seed `numeric` used to fix the random seed.
#' @param repulse `numeric` default 3.5, passed to `layout_with_qfr()`.
#' @param hide_solo_pie `logical` default TRUE, passed to `mem2cnet()`.
#' @param ... additional arguments are passed to `mem2cnet()`.
#'
#' @examples
#' # by default, single-border-color pie is shown as circle
#' cnet1 <- make_cnet_test(border_lwd=2)
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet1, use_shadowText=TRUE)
#' })
#'
#' # hide_solo_pie=FALSE shows every pie wedge bordder
#' cnet2 <- make_cnet_test(hide_solo_pie=FALSE, border_lwd=2)
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet2, use_shadowText=TRUE)
#' })
#' 
#' # Set nodes can be adjusted, reorienting the Gene nodes
#' cnet2_adj <- adjust_cnet_set_relayout_gene(cnet2,
#'    nodes=c("SetB", "SetD"),
#'    x=c(-0.1, 0), y=c(0, -0.2),
#'    repulse=3.6);
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet2_adj, use_shadowText=TRUE, label_dist_factor=0)
#' })
#'
#' # nodeset spacing can be enforced
#' cnet3 <- make_cnet_test(num_sets=3)
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3, use_shadowText=TRUE)
#' })
#' cnet3_sp <- apply_nodeset_spacing(cnet3,
#'    percent_spacing=7)
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3_sp, use_shadowText=TRUE)
#' })
#'
#' # a specific nodeset can be individually adjusted
#' cnet3_adj <- adjust_cnet_nodeset(cnet3_sp,
#'    set_nodes=list(c("SetA", "SetB")),
#'    x=c(-0.2), y=c(0.2))
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3_adj, use_shadowText=TRUE)
#' })
#'
#' # several nodesets can be adjusted at once
#' cnet3_adj2 <- adjust_cnet_nodeset(cnet3_sp,
#'    set_nodes=list("SetA,SetB", "SetA,SetC", "SetB,SetC"),
#'    x=c(-0.2, 0.2, 0), y=c(0.2, 0.2, -0.2))
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3_adj2, use_shadowText=TRUE)
#' })
#'
#' # individual nodes can be nudged
#' cnet3_adj2_nudge <- nudge_igraph_node(cnet3_adj2,
#'    nodes=c("T"), x=c(-0.02), y=c(0.1))
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3_adj2_nudge, use_shadowText=TRUE, vertex.label.font=2)
#' })
#'
#' # nodes can be nudged in larger sets using nodes_xy
#' cnet3_adj2_nudge2 <- nudge_igraph_node(cnet3_adj2,
#'    nodes_xy=list(
#'       T=c(-0.02, 0.2),
#'       AK=c(0.02, 0.2),
#'       AG=c(-0.2, 0),
#'       Q=c(0.2, 0)
#' ))
#' withr::with_par(list("mar"=c(0, 0, 0, 0) + 0.5),{
#'    jam_igraph(cnet3_adj2_nudge2, use_shadowText=TRUE, vertex.label.font=2)
#' })
#'
#' @export
make_cnet_test <- function
(num_sets=4,
 overlap_counts=c(57, 20, 12, 5),
 row_prefix="",
 column_prefix="Set",
 add_direction=TRUE,
 set_colors=NULL,
 seed=123,
 repulse=3.5,
 hide_solo_pie=TRUE,
 ...)
{
   # set seed for reproducibility
   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }

   # create example cnet data
   if (length(overlap_counts) < num_sets) {
      overlap_ratio <- mean(sapply(head(seq_along(overlap_counts), -1), function(i){
         overlap_counts[i+1]/overlap_counts[i]
      }));
      for (i in seq(from=length(overlap_counts)+1, to=num_sets)) {
         overlap_counts <- c(overlap_counts,
            ceiling(tail(overlap_counts, 1) * overlap_ratio))
      }
      overlap_counts
   } else {
      overlap_counts <- head(overlap_counts, num_sets)
   }

   # prepare incidence matrix
   overlap_nums <- rep(seq_len(num_sets), overlap_counts);
   cnetim <- do.call(rbind, lapply(seq_len(sum(overlap_counts)), function(i){
      iset <- rep(c(1, 0),
         c(overlap_nums[i],
            num_sets - overlap_nums[i]))
      sample(iset)
   }))
   rownames(cnetim) <- jamba::colNum2excelName(seq_len(sum(overlap_counts)));
   colnames(cnetim) <- paste0("Set", jamba::colNum2excelName(seq_len(num_sets)));

   # set_colors
   if (length(set_colors) < num_sets) {
      set_colors <- colorjam::rainbowJam(num_sets,
         Crange=c(70, 110),
         ...);
   } else {
      set_colors <- head(set_colors, num_sets);
   }
   names(set_colors) <- colnames(cnetim);

   # make gene color incidence matrix
   geneIMcolors <- do.call(rbind, lapply(rownames(cnetim), function(i){
      j <- cnetim[match(i, rownames(cnetim)),];
      ifelse(j > 0,
         set_colors,
         "#FFFFFF")
   }))
   rownames(geneIMcolors) <- rownames(cnetim);
   colnames(geneIMcolors) <- colnames(cnetim);

   # optionally add directionality
   geneIMdirection <- NULL;
   if (TRUE %in% add_direction) {
      geneIMdirection <- cnetim * sample(c(1, -1),
         size=prod(dim(cnetim)),
         replace=TRUE)
   }

   cnet <- mem2cnet(list(
      memIM=cnetim,
      geneIM=cnetim,
      geneIMdirection=geneIMdirection,
      geneIMcolors=geneIMcolors),
      spread_labels=TRUE,
      repulse=repulse,
      hide_solo_pie=hide_solo_pie,
      ...)
   isset <- igraph::V(cnet)$nodeType %in% "Set";
   igraph::V(cnet)$color[isset] <- set_colors[igraph::V(cnet)$name[isset]];

   igraph::V(cnet)$pie.color[isset] <- as.list(set_colors[igraph::V(cnet)$name[isset]]);
   igraph::V(cnet)$coloredrect.color[isset] <- as.list(set_colors[igraph::V(cnet)$name[isset]]);

   igraph::V(cnet)$size <- igraph::V(cnet)$size * 1;
   igraph::V(cnet)$size2 <- igraph::V(cnet)$size;
   igraph::V(cnet)$label.dist <- ifelse(isset, 0,
   	igraph::V(cnet)$label.dist);
   return(cnet);
}
