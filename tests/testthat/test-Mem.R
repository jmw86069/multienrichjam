
# test-Mem.R

test_that("Mem", {
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

   # convert to Mem
   Mem <- list_to_Mem(mem)
   testthat::expect_setequal(
      class(Mem),
      "Mem")

   # convert back to list
   new_mem <- as(Mem, "list")
   # geneIMdirection was added during import, so omit during the comparison
   new_names <- setdiff(names(new_mem), "geneIMdirection")
   testthat::expect_setequal(
      new_mem[new_names],
      mem[new_names])

   # accessors
   testthat::expect_setequal(
      dim(geneIM(Mem)),
      c(13, 1))
   testthat::expect_setequal(
      dim(enrichIM(Mem)),
      c(4, 1))
   testthat::expect_setequal(
      dim(memIM(Mem)),
      c(13, 4))


   # gene constraints
   geneset <- c("ACTB", "CALM1", "COL4A1", "ESR1", "FKBP5", "GAPDH",
      "MAPK3", "MAPK8", "PPIA", "PTEN", "SGK", "TTN", "ZBTB16");
   testthat::expect_setequal(
      rownames(Mem@memIM),
      geneset)
   testthat::expect_equal(
      rownames(Mem@geneIM),
      geneset)
   testthat::expect_setequal(
      rownames(Mem@geneHitIM),
      c(geneset, "TTP"))
   testthat::expect_equal(
      rownames(Mem@geneIM),
      rownames(Mem@geneIMcolors))
   testthat::expect_equal(
      rownames(Mem@geneIM),
      rownames(Mem@memIM))


})

test_that("Mem.asList", {

})
