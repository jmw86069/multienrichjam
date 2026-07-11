
# test-Mem-combine.R
#
# Tests for combining Mem objects with c()

# Shared test data factory
make_mem <- function(name, enrichdf) {
   multiEnrichMap(setNames(
      list(enrichDF2enrichResult(enrichdf)),
      name))
}

base_enrichdf <- data.frame(check.names=FALSE,
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

alt_enrichdf <- data.frame(check.names=FALSE,
   row.names=paste("Alt Pathway", LETTERS[3:1]),
   ID=paste("Alt Pathway", LETTERS[3:1]),
   Description=paste("Alt Pathway", LETTERS[3:1]),
   pvalue=c(0.002, 0.006, 0.011),
   Count=c(5, 3, 2),
   geneID=c(
      "ACTB/GAPDH/BRCA1/BRCA2/TP53",
      "BRCA1/PTEN/MAPK3",
      "TP53/ESR1"),
   GeneRatio=c("5/80", "3/80", "2/80")
)


test_that("Mem combine - basic two-object c()", {
   Mem1 <- make_mem("EnrichmentA", base_enrichdf)
   Mem2 <- make_mem("EnrichmentB", alt_enrichdf)

   Mem_combined <- c(Mem1, Mem2)

   testthat::expect_s4_class(Mem_combined, "Mem")
   testthat::expect_equal(enrichments(Mem_combined), c("EnrichmentA", "EnrichmentB"))
   testthat::expect_equal(dim(Mem_combined)[3], 2L)   # 2 enrichments
   testthat::expect_true(check_Mem(Mem_combined))
})


test_that("Mem combine - genes and sets are unioned and sorted", {
   Mem1 <- make_mem("EnrichmentA", base_enrichdf)
   Mem2 <- make_mem("EnrichmentB", alt_enrichdf)

   Mem_combined <- c(Mem1, Mem2)

   expected_genes <- jamba::mixedSort(union(genes(Mem1), genes(Mem2)))
   expected_sets  <- jamba::mixedSort(union(sets(Mem1),  sets(Mem2)))

   testthat::expect_equal(genes(Mem_combined), expected_genes)
   testthat::expect_equal(sets(Mem_combined),  expected_sets)
})


test_that("Mem combine - three-object c()", {
   base2 <- base_enrichdf; base2$pvalue <- base2$pvalue + 0.002
   base3 <- base_enrichdf; base3$pvalue <- base3$pvalue + 0.005

   Mem1 <- make_mem("TestA", base_enrichdf)
   Mem2 <- make_mem("TestB", base2)
   Mem3 <- make_mem("TestC", base3)

   Mem_combined <- c(Mem1, Mem2, Mem3)

   testthat::expect_equal(dim(Mem_combined)[3], 3L)
   testthat::expect_equal(enrichments(Mem_combined), c("TestA", "TestB", "TestC"))
   testthat::expect_true(check_Mem(Mem_combined))
})


test_that("Mem combine - error on overlapping enrichment names", {
   Mem1 <- make_mem("SameName", base_enrichdf)
   Mem2 <- make_mem("SameName", base_enrichdf)

   testthat::expect_error(c(Mem1, Mem2), "overlapping enrichment names")
})


test_that("Mem combine - error on non-Mem argument", {
   Mem1 <- make_mem("EnrichmentA", base_enrichdf)

   testthat::expect_error(c(Mem1, "not_a_mem"), "must be Mem objects")
})


test_that("Mem combine - preserves first object metadata", {
   base2 <- base_enrichdf; base2$pvalue <- base2$pvalue + 0.002

   Mem1 <- multiEnrichMap(
      list(TestA=enrichDF2enrichResult(base_enrichdf)),
      topEnrichN=5, cutoffRowMinP=0.05)
   Mem2 <- multiEnrichMap(
      list(TestB=enrichDF2enrichResult(base2)),
      topEnrichN=10, cutoffRowMinP=0.01)

   Mem_combined <- c(Mem1, Mem2)

   testthat::expect_equal(thresholds(Mem_combined), thresholds(Mem1))
   testthat::expect_equal(headers(Mem_combined),    headers(Mem1))
})


test_that("Mem combine - matrix dimensions are correct", {
   Mem1 <- make_mem("EnrichmentA", base_enrichdf)
   Mem2 <- make_mem("EnrichmentB", alt_enrichdf)

   Mem_combined <- c(Mem1, Mem2)
   n_genes  <- length(jamba::mixedSort(union(genes(Mem1), genes(Mem2))))
   n_sets   <- length(jamba::mixedSort(union(sets(Mem1),  sets(Mem2))))
   n_enrich <- 2L

   testthat::expect_equal(nrow(geneIM(Mem_combined)),       n_genes)
   testthat::expect_equal(ncol(geneIM(Mem_combined)),       n_enrich)
   testthat::expect_equal(nrow(enrichIM(Mem_combined)),     n_sets)
   testthat::expect_equal(ncol(enrichIM(Mem_combined)),     n_enrich)
   testthat::expect_equal(dim(memIM(Mem_combined)),         c(n_genes, n_sets))
})
