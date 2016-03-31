#' Fit models using limma
#'
#' This function will fit the flat and nesed models, and tests contrasts
#'
#' @param GR GRanges-object.
#' @param bed_fname Path to where file will be saved. Must end in ".bed". Defaults to tempfile().
#' @return Path to file.
#' @author Malte Thodberg
#' @import magrittr SummarizedExperiment dplyr
#' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
#' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' @export
fitAPModel <- function(rse, modelMatrix, groupingVariable) {
	### Check input
	message("Performing checks...")
	# Class checks
	# Check on design and contrast matrix
	# Check grouping column exists.



	### Sum before fitting models
	message("Summing counts over genes...")

	# Promoter counts
	EMPromoter <- assay(rse, "counts")

	# Grouping variable
	geneIds <- mcols(rse)[,groupingVariable]

	# Gene counts
	EMGene <- EMPromoter %>%
		tibble::as_data_frame() %>%
		mutate(geneIds=geneIds) %>%
		group_by(geneIds) %>%
		summarise_each(funs(sum)) %>%
		ungroup

	# Resave as matrix
	EMGene <- data.frame(EMGene[,-1],
											 row.names=EMGene$geneIds) %>%
		as.matrix

	### Normalizing
	message("Normalizing and weighting...")

	# Normalize using promoters
	dgePromoter <- EMPromoter %>%
		edgeR::DGEList() %>%
		edgeR::calcNormFactors()

	# Transfer norm factors to genes
	dgeGene <- edgeR::DGEList(EMGene)
	dgeGene$samples$norm.factors <- dgePromoter$samples$norm.factors

	# Estimate weights
	#vPromoter <- voomWithQualityWeights(dgePromoter, design=mod)
	#vGene <- voomWithQualityWeights(dgeGene, design=mod)

	# Estimate only gene weights
	vPromoter <- limma::voom(dgePromoter, design=modelMatrix)
	vGene <- limma::voom(dgeGene, design=modelMatrix)

	### Fit models
	message("Fitting flat and nested models")

	# Fit flat models
	fitPromoter <- vPromoter %>%
		limma::lmFit(object=., design=modelMatrix) %>%
		limma::contrasts.fit(fit=., contrasts=contrastMatrix)

	ebPromoter <- limma::eBayes(fit=fitPromoter, robust=TRUE)

	fitGene <- vGene %>%
		limma::lmFit(object=., design=modelMatrix) %>%
		limma::contrasts.fit(fit=., contrasts=contrastMatrix)

	ebGene <- limma::eBayes(fit=fitGene, robust=TRUE)

	# Fit nested model
	ebNested <- limma::diffSplice(fit=fitPromoter,
												 geneid=geneIds, exonid=rownames(EMPromoter),
												 robust=TRUE, verbose=FALSE)

	### Test desired contrasts
	message("Testing contrasts...")

	#
	# ### Batch correct
	# message("Obtaining corrected expression values...")
	# EM <- removeBatchEffect(x=vPromoter, design=mod[,-c(2:4)], covariates=mod[,2:4])
	#

}


extractContrastRes <- function(contrast, flatPromoter, flatGene, Nested)
# Scaffold to merge on
scaffold <- data_frame(GeneID=geneIds, ExonID=mcols(SE)$name, feat_cat=mcols(SE)$feat_cat)

# Add flat promoter
scaffold <- topTable(ebPromoter, coef=conCol, number=Inf, sort.by="none", confint=TRUE) %>%
	set_colnames(paste0("flatPromoter_", colnames(.))) %>%
	mutate(ExonID=rownames(.)) %>%
	left_join(x=scaffold, y=., by="ExonID")

# Add flat promoter
scaffold <- topTable(ebGene, coef=conCol, number=Inf, sort.by="none", confint=TRUE) %>%
	set_colnames(paste0("flatGene_", colnames(.))) %>%
	mutate(GeneID=rownames(.)) %>%
	left_join(x=scaffold, y=., by="GeneID")

# Add nested gene
simesTest <- topSplice(fit=ebNested, coef=conCol, test="simes", number=Inf) %>%
	set_colnames(c("GeneID", "NExons", "simes.P.Value", "simes.FDR"))

scaffold <- full_join(x=topSplice(fit=ebNested, coef=conCol, test="F", number=Inf) %>%
												set_colnames(c("GeneID", "NExons", "nestedGene_F", "nestedGene_FTest.P.Value", "nestedGene_FTest.FDR")),
											y=topSplice(fit=ebNested, coef=conCol, test="simes", number=Inf) %>%
												set_colnames(c("GeneID", "NExons", "nestedGene_simes.P.Value", "nestedGene_simes.FDR")),
											by=c("GeneID", "NExons")) %>%
	left_join(x=scaffold, y=., by="GeneID")

# Add nested promoter

scaffold <- topSplice(fit=ebNested, coef=conCol, test="t", number=Inf) %>%
	set_colnames(c("ExonID", "GeneID", "nestedPromoter_logFC", "nestedPromoter_t", "nestedPromoter_P.Value", "nestedPromoter_FDR")) %>%
	left_join(x=scaffold, y=., by=c("GeneID", "ExonID"))


#' Fit models using limma
#'
#' This function will fit the flat and nesed models, and tests contrasts
#'
#' @param GR GRanges-object.
#' @param bed_fname Path to where file will be saved. Must end in ".bed". Defaults to tempfile().
#' @return Path to file.
#' @author Malte Thodberg
#' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
#' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' @export
summariseAPModel <- function(){
	1
}
