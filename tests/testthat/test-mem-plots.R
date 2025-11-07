
# test Mem manipulations

testthat::test_that("mem_plots", {
   #
   # test first gp heatmap
   set.seed(123)
   hco <- NULL;
   hro <- NULL;
   gp_hm <- mem_gene_path_heatmap(Memtest, column_split=4,
      row_split=3, column_cex=0.5)
   # withr::local_options(list(warn=-1));
   withr::with_options(list(warn=-1), {
      hco <- jamba::heatmap_column_order(gp_hm);
      hro <- jamba::heatmap_row_order(gp_hm);
   })
   testthat::expect_equal(
      lengths(hro),
      c(a=4, b=3, c=15))
   testthat::expect_equal(
      lengths(hco),
      c(A=6, B=4, C=2, D=4))
   gp_hm_fn <- function() {
      ComplexHeatmap::draw(gp_hm,
         annotation_legend_list=attributes(gp_hm)$caption_legendlist,
         merge_legends=TRUE)
   }
   if (jamba::check_pkg_installed("vdiffr")) {
      vdiffr::expect_doppelganger("GP Heatmap",
         gp_hm_fn)
   }

   # column_split with named vector
   withr::with_options(list(warn=-1), {
      hco <- jamba::heatmap_column_order(gp_hm)
   })
   gcs <- jamba::nameVector(rep(names(hco), lengths(hco)), unlist(hco))
   gcs[] <- gsub("[BC]+", "BC", gcs)
   gcs <- factor(gcs, levels=c("BC", "A", "D"))
   gp_hm_colgrp <- mem_gene_path_heatmap(Memtest,
      column_cex=0.5,
      column_split=gcs, row_split=5)
   gp_hm_fn2 <- function() {
      ComplexHeatmap::draw(gp_hm_colgrp,
         annotation_legend_list=attributes(gp_hm)$caption_legendlist,
         merge_legends=TRUE)
   }
   withr::with_options(list(warn=-1), {
      hco <- jamba::heatmap_column_order(gp_hm_colgrp);
   })
   testthat::expect_equal(
      lengths(hco),
      c(BC=6, A=6, D=4))

   # custom cluster_columns function
   if (requireNamespace("cluster", quietly=TRUE)) {
      gp_hm_clfn <- function() {
         colfn <- function(x){
            withr::with_options(list(warn=-1), {
               stats::hclust(cluster::daisy(x, metric="gower"))
            })
         }
         gp_hm_colfn <- mem_gene_path_heatmap(Memtest,
            column_method="gower",
            column_names_max_height=grid::unit(8, "cm"),
            cluster_columns=colfn, column_split=4, row_split=5)
         ComplexHeatmap::draw(gp_hm_colfn,
            annotation_legend_list=attributes(gp_hm_colfn)$caption_legendlist,
            merge_legends=TRUE)
      }
      if (jamba::check_pkg_installed("vdiffr")) {
         vdiffr::expect_doppelganger("GP Heatmap Column-cluster-function",
            gp_hm_clfn)
      }
   }

   # mem_plot_folio gp heatmap
   set.seed(123)
   mpf <- mem_plot_folio(Memtest, pathway_column_split=4,
      do_which=c(1, 2, 4), do_plot=FALSE,
      row_method="euclidean",
      gene_row_split=3, column_cex=0.5)
   withr::with_options(list(warn=-1), {
      hro <- jamba::heatmap_row_order(mpf$gp_hm);
      hco <- jamba::heatmap_column_order(mpf$gp_hm);
   })
   testthat::expect_equal(
      lengths(hro),
      c(a=18, b=2, c=2))
   testthat::expect_equal(
      lengths(hco),
      c(A=6, B=4, C=2, D=4))
   gp_hm_fn <- function() {
      ComplexHeatmap::draw(mpf$gp_hm,
         annotation_legend_list=attributes(mpf$gp_hm)$caption_legendlist,
         merge_legends=TRUE)
   }
   if (jamba::check_pkg_installed("vdiffr")) {
      vdiffr::expect_doppelganger("mpf GP Heatmap",
         gp_hm_fn)
   }
   
   # cnet_cluster
   cnet <- mpf$cnet_collapsed_set;
   igraph::V(cnet)$label.dist <- ifelse(igraph::V(cnet)$nodeType %in% "Gene",
      igraph::V(cnet)$label.dist,
      0)
   # test Cnet nodesets
   nodesets <- get_cnet_nodeset(cnet);
   testthat::expect_equal(
      lengths(nodesets),
      c(A=5, "A,B"=3, "A,B,C,D"=2, "A,B,D"=1,
         "A,C"=1, "A,D"=1, B=1, "B,C"=1,
         "B,D"=1, C=1, "C,D"=1, D=4))
   
   # visual Cnet plot
   cnet_fn <- function() {
      jam_igraph(cnet,
         label_dist_factor=5,
         use_shadowText=TRUE,
         node_factor_l=list(nodeType=c(Gene=1.5, Set=2)),
         label_factor_l=list(nodeType=c(Gene=1.5, Set=0.9)))
   }
   if (jamba::check_pkg_installed("vdiffr")) {
      vdiffr::expect_doppelganger("mpf Cnet cluster plot",
         cnet_fn)
   }
   
   # Cnet manipulations
   cnet2 <- relayout_with_qfr(cnet, repulse=3)
   cnet2_fn <- function() {
      jam_igraph(cnet2,
         label_dist_factor=5,
         use_shadowText=TRUE,
         node_factor_l=list(nodeType=c(Gene=1.5, Set=2)),
         label_factor_l=list(nodeType=c(Gene=1.5, Set=0.9)))
   }
   if (jamba::check_pkg_installed("vdiffr")) {
      vdiffr::expect_doppelganger("mpf Cnet2 cluster plot",
         cnet2_fn)
   }
   
   # Cnet node spacing
   cnet_spacing <- summarize_node_spacing(cnet)$nearest_within[c("A", "A,B",
      "A,B,C,D", "D"), ]
   testthat::expect_equal(
      round(digits=2, summary(as.vector(cnet_spacing))),
      c("Min."=6.91, "1st Qu."=9.63, Median=11.53,
         Mean=11.23, "3rd Qu."=13.89, "Max."=14.85))
   
   cnet2_spacing <- summarize_node_spacing(cnet2)$nearest_within[c("A", "A,B",
      "A,B,C,D", "D"), ]
   testthat::expect_equal(
      round(digits=1, summary(as.vector(cnet2_spacing))),
      c("Min."=16.4, "1st Qu."=17.2, Median=18.2,
         Mean=19.8, "3rd Qu."=18.6, "Max."=48.9))

   # Todo: Cnet manipulations
   # adjust_cnet_nodeset(), nudge_igraph_node()
})
