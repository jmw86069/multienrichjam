
test_that("enrichList2df", {
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
   pw <- paste("Pathway", LETTERS[1:5])
   # test single enrichment input
   enrichdf1 <- enrichList2df(list(EnrichmentA=test_enrichdf))
   testthat::expect_equal(
      enrichdf1[pw, c("pvalue", "ID", "Count", "Description")],
      test_enrichdf[pw, c("pvalue", "ID", "Count", "Description")])
   testthat::expect_equal(
      enrichdf1[pw, "geneID"],
      gsub("/", ",", test_enrichdf[pw, "geneID"]))

   # test custom geneSep
   enrichdf1sep <- enrichList2df(list(EnrichmentA=test_enrichdf), geneSep="/")
   testthat::expect_equal(
      enrichdf1sep[pw, c("pvalue", "ID", "Count", "Description")],
      test_enrichdf[pw, c("pvalue", "ID", "Count", "Description")])
   testthat::expect_equal(
      enrichdf1sep[pw, "geneID"],
      test_enrichdf[pw, "geneID"])

   # test two enrichment input
   test_enrichdf2 <- test_enrichdf;
   test_enrichdf2$pvalue <- test_enrichdf$pvalue / 2;
   test_enrichdf2$geneID <- paste0(test_enrichdf$geneID, "/ZZZ");
   test_enrichdf2$Count <- test_enrichdf$Count + 1;
   enrichdf2 <- enrichList2df(list(EnrichmentA=test_enrichdf,
      EnrichmentB=test_enrichdf2))
   testthat::expect_equal(
      enrichdf2[pw, c("pvalue", "ID", "Count", "Description")],
      test_enrichdf2[pw, c("pvalue", "ID", "Count", "Description")])
   testthat::expect_equal(
      enrichdf2[pw, "geneID"],
      gsub("/", ",", test_enrichdf2[pw, "geneID"]))

   # test two enrichment input, one is not significant
   test_enrichdf3 <- test_enrichdf2;
   test_enrichdf3$pvalue <- 1;
   enrichdf3 <- enrichList2df(list(EnrichmentA=test_enrichdf,
      EnrichmentB=test_enrichdf3))
   testthat::expect_equal(
      enrichdf3[pw, c("pvalue", "ID", "Description")],
      test_enrichdf[pw, c("pvalue", "ID", "Description")])
   testthat::expect_equal(
      enrichdf3[pw, c("Count")],
      test_enrichdf3[pw, c("Count")])
   testthat::expect_equal(
      enrichdf3[pw, "geneID"],
      gsub("/", ",", test_enrichdf2[pw, "geneID"]))

   # test two enrichment input, non-overlapping pathways
   test_enrichdf4 <- test_enrichdf2;
   test_enrichdf4$ID <- paste("Pathway", LETTERS[5:9])
   test_enrichdf4$Description <- paste("Description", LETTERS[5:9])
   rownames(test_enrichdf4) <- test_enrichdf4$ID;
   enrichdf4 <- enrichList2df(list(EnrichmentA=test_enrichdf,
      EnrichmentB=test_enrichdf4))
   pw1 <- paste("Pathway", LETTERS[1:4])
   pw4 <- paste("Pathway", LETTERS[5:9])
   testthat::expect_equal(
      enrichdf4[c(pw1, pw4), c("pvalue", "ID", "Description")],
      rbind(test_enrichdf[pw1, c("pvalue", "ID", "Description")],
         test_enrichdf4[pw4, c("pvalue", "ID", "Description")]))
   testthat::expect_equal(
      enrichdf4[c(pw1, pw4), c("Count")],
      c(test_enrichdf[pw1, c("Count")],
         test_enrichdf4[pw4, c("Count")]))
   testthat::expect_equal(
      enrichdf4[c(pw1, pw4), "geneID"],
      c(gsub("/", ",", test_enrichdf[pw1, "geneID"]),
         gsub("/", ",", test_enrichdf4[pw4, "geneID"])))

   # test two enrichment input, one with empty gene
   test_enrichdf5 <- test_enrichdf4;
   test_enrichdf5$geneID[5] <- "";
   test_enrichdf5$Count[5] <- 0;
   test_enrichdf5
   enrichdf5 <- enrichList2df(list(EnrichmentA=test_enrichdf,
      EnrichmentB=test_enrichdf5))
   testthat::expect_equal(
      enrichdf5[c(pw1, pw4), c("pvalue", "ID", "Description")],
      rbind(test_enrichdf[pw1, c("pvalue", "ID", "Description")],
         test_enrichdf5[pw4, c("pvalue", "ID", "Description")]))
   testthat::expect_equal(
      enrichdf5[c(pw1, pw4), c("Count")],
      c(test_enrichdf[pw1, c("Count")],
         test_enrichdf5[pw4, c("Count")]))
   testthat::expect_equal(
      enrichdf5[c(pw1, pw4), "geneID"],
      c(gsub("/", ",", test_enrichdf[pw1, "geneID"]),
         gsub("/", ",", test_enrichdf5[pw4, "geneID"])))

})
