# test fixSetLabels()

test_that("fixSetLabels examples", {
	a <- c(
		"Genes involved in Gluconeogenesis",
		"Genes Up-regulated during fasting",
		"Genes regulated by Tnfalpha, TGF beta, Pi3kakt, mtorc"
	)
	fixSetLabels(a, wrap=FALSE)
	
	testthat::expect_equal(
		fixSetLabels(a, wrap=FALSE),
		c("Gluconeogenesis",
			"Up-regulation During Fasting",
			"Regulation By TNFa, TGFbeta, PI3K/AKT, mTORc"))

	testthat::expect_equal(
		fixSetLabels(a),
		c("Gluconeogenesis",
			"Up-regulation During Fasting",
			"Regulation By TNFa, TGFbeta, PI3K/AKT,\nmTORc"))
	
	testthat::expect_equal(
		fixSetLabels(a, do_abbreviations=FALSE),
		c("Genes Involved In Gluconeogenesis",
			"Genes Up-Regulated During Fasting",
			"Genes Regulated By TNFa, TGFbeta,\nPI3K/AKT, mTORc"))
	
	b <- c(
		"Extracellular Matrix effects on lincRNA, mirna, interleukin-1, il-1b",
		"PPAR-alpha signaling pathways",
		"Clear signaling pathway")
	testthat::expect_equal(
		fixSetLabels(b),
		c(
			"ECM Effects On lncRNA, miRNA, IL-1,\nIL-1b",
			"PPARalpha Signaling",
			"CLEAR Signaling"))

})
