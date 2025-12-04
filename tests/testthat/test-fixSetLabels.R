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
			"Regulation By TNFa, TGF-beta, PI3K/AKT, mTORC"))

	testthat::expect_equal(
		fixSetLabels(a),
		c("Gluconeogenesis",
			"Up-regulation During Fasting",
			"Regulation By TNFa, TGF-beta, PI3K/AKT,\nmTORC"))
	
	testthat::expect_equal(
		fixSetLabels(a, do_abbreviations=FALSE),
		c("Genes Involved In Gluconeogenesis",
			"Genes Up-Regulated During Fasting",
			"Genes Regulated By TNFa, TGF-beta,\nPI3K/AKT, mTORC"))
	
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

	# Special abbreviations
	# TNF
	# c-FLIP
	# Interferons / IFN
	ifns <- c("Tumor necrosisFactor alpha",
		"Regulation by c flip",
		"Class A 1 Rhodopsin Like Receptors",
		"Interferon beta 1a",
		"IFN alpha 2", "IFN beta", "IFNb2",
		"IFNgamma", "IFN gamma 2a")
	testthat::expect_equal(
		fixSetLabels(ifns),
		c("TNFa", "Regulation By c-FLIP",
			"Class A 1 Rhodopsin-Like Receptors",
			"IFN-beta-1a", "IFN-alpha-2", "IFN-beta",
			"IFN-beta-2", "IFN-gamma", "IFN-gamma-2a"))
	
})
