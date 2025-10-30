
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
   mem <- multiEnrichMap(erlist)

   pw <- paste("Pathway", LETTERS[2:5]);
   testthat::expect_equal(
      sort(rownames(mem$multiEnrichDF)),
      pw)
   testthat::expect_equal(
      mem$multiEnrichDF[pw, c("pvalue", "ID", "Count", "Description")],
      test_enrichdf[pw, c("pvalue", "ID", "Count", "Description")])
   testthat::expect_equal(
      mem$multiEnrichDF[pw, "geneID"],
      gsub("/", ",", test_enrichdf[pw, "geneID"]))

   # gene constraints
   geneset <- c("ACTB", "CALM1", "COL4A1", "ESR1", "FKBP5", "GAPDH",
      "MAPK3", "MAPK8", "PPIA", "PTEN", "SGK", "TTN", "ZBTB16");
   testthat::expect_setequal(
      rownames(mem$memIM),
      geneset)
   testthat::expect_equal(
      rownames(mem$geneIM),
      geneset)
   testthat::expect_setequal(
      rownames(mem$geneHitIM),
      c(geneset, "TTP"))
   testthat::expect_equal(
      rownames(mem$geneIM),
      rownames(mem$geneIMcolors))
   testthat::expect_equal(
      rownames(mem$geneIM),
      rownames(mem$memIM))

   # memIM constraints
   desc <- paste("Description", LETTERS[2:5])
   testthat::expect_setequal(
      colnames(mem$memIM),
      desc)
   testthat::expect_equal(ignore_attr=TRUE,
      colSums(mem$memIM)[desc],
      test_enrichdf[pw, "Count"])

   # memIM-enrichIM consistent order
   testthat::expect_equal(
      colnames(mem$memIM),
      rownames(mem$enrichIM))
   # memIM-geneIM consistent order
   testthat::expect_equal(
      rownames(mem$memIM),
      rownames(mem$geneIM))

   # enrichIM constraints
   testthat::expect_equal(
      rownames(mem$enrichIM),
      rownames(mem$enrichIMcolors))
   testthat::expect_equal(
      rownames(mem$enrichIM),
      rownames(mem$enrichIMgeneCount))
   testthat::expect_equal(
      rownames(mem$enrichIM),
      rownames(mem$enrichIMdirection))

   ## mem_plot_folio
   # gene-path heatmap
   mpf2 <- mem_plot_folio(mem, do_plot=FALSE,
      do_which=2, gene_column_split=1, gene_column_title="A")
   testthat::expect_contains(
      names(mpf2),
      "gp_hm")
   testthat::expect_contains(
      strsplit(mpf2$gp_hm_caption, "\n")[[1]],
      c("11 genes (rows)",
         "3 sets (columns)",
         "enrichment P <= 0.05"))
   # enrichment heatmap
   mpf1 <- mem_plot_folio(mem, do_plot=FALSE, do_which=1)
   testthat::expect_contains(
      names(mpf1),
      "enrichment_hm")
   # Cnet plot 1
   mpf3 <- mem_plot_folio(mem, do_plot=FALSE,
      pathway_column_split=3, do_which=3)
   testthat::expect_contains(
      names(mpf3),
      "cnet_collapsed")
   # verify nodes in the Cnet igraph
   cnet1 <- mpf3$cnet_collapsed;
   testthat::expect_setequal(
      igraph::V(cnet1)$name,
      c(setdiff(geneset, c("ESR1", "ZBTB16")), LETTERS[1:3]))

   # Cnet plot 2
   mpf4 <- mem_plot_folio(mem, do_plot=FALSE, do_which=4)
   testthat::expect_contains(
      names(mpf4),
      "cnet_collapsed_set")

})
