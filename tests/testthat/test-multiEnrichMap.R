
test_that("multiEnrichMap", {
   test_enrichdf <- data.frame(check.names=FALSE,
      row.names=paste("Pathway", LETTERS[5:1]),
      ID=paste("Pathway", LETTERS[5:1]),
      Description=paste("Description", LETTERS[5:1]),
      pvalue=c(0.001, 0.005, 0.01, 0.05, 0.10),
      Count=c(7, 5, 4, 3, 2),
      geneID=c(
         "ACTB/FKBP5/GAPDH/MAPK3/MAPK8/PPIA/PTEN",
         "CALM1/COL4A1/MAPK3/PTEN/SGK",
         "CALM1/GAPDH/PPIA/TTN",
         "ACTB/ESR1/ZBTB16",
         "ESR1/TTP"),
      GeneRatio=c("7/100", "5/90", "4/85", "3/110", "2/200")
   )

   erlist <- list(EnrichmentA=enrichDF2enrichResult(test_enrichdf))
   Mem <- multiEnrichMap(erlist)
   mem <- Mem_to_list(Mem)

   pw <- paste("Pathway", LETTERS[2:5]);
   testthat::expect_equal(
      sort(rownames(Mem@multiEnrichDF)),
      pw)
   testthat::expect_equal(
      Mem@multiEnrichDF[pw, c("pvalue", "ID", "Count", "Description")],
      test_enrichdf[pw, c("pvalue", "ID", "Count", "Description")])
   testthat::expect_equal(
      Mem@multiEnrichDF[pw, "geneID"],
      gsub("/", ",", test_enrichdf[pw, "geneID"]))

   # gene constraints
   geneset <- c("ACTB", "CALM1", "COL4A1", "ESR1", "FKBP5", "GAPDH",
      "MAPK3", "MAPK8", "PPIA", "PTEN", "SGK", "TTN", "ZBTB16");
   testthat::expect_setequal(
      rownames(memIM(Mem)),
      geneset)
   testthat::expect_equal(
      rownames(geneIM(Mem)),
      geneset)
   testthat::expect_setequal(
      rownames(geneHitIM(Mem)),
      c(geneset, "TTP"))
   testthat::expect_equal(
      rownames(geneIM(Mem)),
      rownames(geneIMcolors(Mem)))
   testthat::expect_equal(
      rownames(geneIM(Mem)),
      rownames(memIM(Mem)))

   # memIM constraints
   desc <- paste("Description", LETTERS[2:5])
   testthat::expect_setequal(
      colnames(memIM(Mem)),
      desc)
   testthat::expect_equal(ignore_attr=TRUE,
      colSums(memIM(Mem))[desc],
      test_enrichdf[pw, "Count"])

   # memIM-enrichIM consistent order
   testthat::expect_equal(
      colnames(memIM(Mem)),
      rownames(enrichIM(Mem)))
   # memIM-geneIM consistent order
   testthat::expect_equal(
      rownames(memIM(Mem)),
      rownames(geneIM(Mem)))

   # enrichIM constraints
   testthat::expect_equal(
      rownames(enrichIM(Mem)),
      rownames(enrichIMcolors(Mem)))
   testthat::expect_equal(
      rownames(enrichIM(Mem)),
      rownames(enrichIMgeneCount(Mem)))
   testthat::expect_equal(
      rownames(enrichIM(Mem)),
      rownames(enrichIMdirection(Mem)))

   ## mem_plot_folio
   # gene-path heatmap
   mpf2 <- mem_plot_folio(Mem, do_plot=FALSE,
      do_which=2, gene_column_split=1, gene_column_title="A")
   testthat::expect_true(
      inherits(mpf2@gp_hm, "Heatmap"),
      TRUE)
   testthat::expect_contains(
      strsplit(Caption(mpf2), "\n")[[1]],
      c("13 genes (rows)",
         "4 sets (columns)",
         "enrichment P <= 0.05"))
   
   # gene-path heatmap with filtering criteria
   mpf2b <- mem_plot_folio(Mem, do_plot=FALSE,
      min_set_ct_each=4,
      do_which=2, gene_column_split=1, gene_column_title="A")
   testthat::expect_contains(
      strsplit(Caption(mpf2b), "\n")[[1]],
      c("11 genes (rows)",
         "3 sets (columns)",
         "enrichment P <= 0.05"))

   # enrichment heatmap
   mpf1 <- mem_plot_folio(Mem, do_plot=FALSE, do_which=1)
   testthat::expect_true(
      inherits(mpf1@enrichment_hm, "Heatmap"),
      TRUE)
   
   # Cnet plot 1
   mpf3 <- mem_plot_folio(Mem, do_plot=FALSE,
      pathway_column_split=3, do_which=3)
   testthat::expect_true(
      inherits(CnetCollapsed(mpf3, do_plot=FALSE), "igraph"),
      TRUE)
   # verify nodes in the Cnet igraph
   cnet1 <- CnetCollapsed(mpf3, do_plot=FALSE);
   testthat::expect_setequal(
      igraph::V(cnet1)$name,
      c(geneset, LETTERS[1:3]))

   # Cnet plot 1b with filter critera
   mpf3b <- mem_plot_folio(Mem, do_plot=FALSE,
      min_set_ct_each=4,
      pathway_column_split=3, do_which=3)
   # verify nodes in the Cnet igraph
   cnet1b <- CnetCollapsed(mpf3b, do_plot=FALSE);
   testthat::expect_setequal(
      igraph::V(cnet1b)$name,
      c(setdiff(geneset, c("ESR1", "ZBTB16")), LETTERS[1:3]))
   
   # Cnet plot 2
   mpf4 <- mem_plot_folio(Mem, do_plot=FALSE, do_which=4)
   testthat::expect_true(
      inherits(CnetCollapsed(mpf3, do_plot=FALSE, type="set"), "igraph"),
      TRUE)

})
