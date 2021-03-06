#' Fit GLMs at the gene and promoter level
#'
#' Uses edgeR to fit a robust QL GLM to promoters and genes
#'
#' @param promoterCounts matrix: Number of CAGE tags in each promoter.
#' @param modelMatrix matrix: Full rank model matrix describing the experiment.
#' @param geneIds character: Vector marking which genes promoters belong two.
#'
#' @import magrittr
#'
#' @return APfit object
#' @examples
#' # ADD EXAMPLES HERE
#' @export
APfit <- function(promoterCounts, modelMatrix, geneIds){
		### Check input
		message("Performing checks...")
		# Class checks
		# Check on design and contrast matrix
		# Check grouping column exists.

		### Gene strucutre
		geneStructure <- dplyr::data_frame(promoterIds=rownames(promoterCounts),
																geneIds=geneIds)

		### Sum before fitting models
		message("Summing counts over genes...")

		# Gene counts
		geneCounts <- data.frame(geneIds=geneIds, promoterCounts) %>%
			dplyr::group_by(geneIds) %>%
			dplyr::summarise_each(funs(sum)) %>%
			dplyr::ungroup()

		# Resave as matrix
		geneCounts <- data.frame(geneCounts[,-1],
												 row.names=geneCounts$geneIds) %>%
			as.matrix

		### Normalizing
		message("Normalizing...")

		# Normalize using promoters
		dgePromoter <- promoterCounts %>%
			edgeR::DGEList() %>%
			edgeR::calcNormFactors()

		# Transfer norm factors to genes
		dgeGene <- edgeR::DGEList(geneCounts)
		dgeGene$samples$norm.factors <- dgePromoter$samples$norm.factors

		# Get TPM values
		tpmPromoter <- edgeR::cpm(x=dgePromoter, normalized.lib.sizes=TRUE, log=TRUE)
		tpmGene <- edgeR::cpm(x=dgeGene, normalized.lib.sizes=TRUE, log=TRUE)

		### Dispersion
		message("Estimating Dispersions...")
		dispPromoter <- edgeR::estimateDisp(y=dgePromoter, design=modelMatrix, robust=TRUE)
		dispGene <- edgeR::estimateDisp(y=dgeGene, design=modelMatrix, robust=TRUE)

		### Fitting model
		message("Fitting GLMs...")
		fitPromoter <- edgeR::glmQLFit(y=dispPromoter, design=modelMatrix, robust=TRUE)
		fitGene <- edgeR::glmQLFit(y=dispGene, design=modelMatrix, robust=TRUE)

		# Return
		o <- list(genes=geneStructure,
							modelMatrix=modelMatrix,
							dge=list(promoter=dgePromoter, gene=dgeGene),
							tpm=list(promoter=tpmPromoter, gene=tpmGene),
							disp=list(promoter=dispPromoter, gene=dispGene),
							fit=list(promoter=fitPromoter, gene=fitGene))

		class(o) <- c("APfit", class(o))

		o
}

#' Internal function: Test single contrast
#'
#' Called by APtest
#'
#' @param o APfit object
#' @param contrast numeric: Valid contrast.
#'
#' @return RETURN list with info returned by APtest
#' @examples
#' # ADD EXAMPLES HERE
#' @export
singleContrast <- function(o, contrast){
	### Check input
	message("Performing checks...")
	stopifnot(length(contrast) == ncol(o$modelMatrix))

	### Flat tests
	message("Testing flat models...")

	testPromoter <- edgeR::glmQLFTest(glmfit=o$fit$promoter, contrast=contrast)
	testGene <- edgeR::glmQLFTest(glmfit=o$fit$gene, contrast=contrast)

	### Nested test
	testNested <- edgeR::diffSpliceDGE(glmfit=o$fit$promoter,
															contrast=contrast,
															geneid=o$genes$geneIds,
															exonid=o$genes$promoterIds,
															prior.count=0.125,
															verbose=FALSE)

	### Extract info
	message("Extracting info...")

	# From flat
	flatPromoter <- data.frame(edgeR::topTags(object=testPromoter, n=Inf))
	flatGene <- data.frame(edgeR::topTags(object=testGene, n=Inf))

	# From Nested
	nestedPromoter <- edgeR::topSpliceDGE(lrt=testNested, test="exon", number=Inf, FDR=1)
	nestedGene <- edgeR::topSpliceDGE(lrt=testNested, test="gene", number=Inf, FDR=1)
	nestedSimes <- edgeR::topSpliceDGE(lrt=testNested, test="Simes", number=Inf, FDR=1)

	### Extract info
	message("Merging results...")

	# Rename promoter stats
	colnames(flatPromoter) <- paste0("flatPromoter_", colnames(flatPromoter))
	flatPromoter$promoterIds <- rownames(flatPromoter)

	colnames(nestedPromoter) <- paste0("nestedPromoter_", colnames(nestedPromoter))
	colnames(nestedPromoter)[1:2] <- c("promoterIds", "geneIds")

	# Rename gene stats
	colnames(flatGene) <- paste0("flatGene_", colnames(flatGene))
	flatGene$geneIds <- rownames(flatGene)

	colnames(nestedGene) <- paste0("nestedGene_", colnames(nestedGene))
	colnames(nestedGene)[1:2] <- c("geneIds", "nPromoters")

	colnames(nestedSimes) <- paste0("nestedSimes_", colnames(nestedSimes))
	colnames(nestedSimes)[1:2] <- c("geneIds", "nPromoters")

	# Use original order as scaffold
	mergedResults <- o$genes

	mergedResults <- dplyr::left_join(x=mergedResults, y=flatGene, by=c("geneIds"))
	mergedResults <- dplyr::left_join(x=mergedResults, y=nestedGene, by=c("geneIds"))
	mergedResults <- dplyr::left_join(x=mergedResults, y=nestedSimes, by=c("geneIds", "nPromoters"))

	mergedResults <- dplyr::left_join(x=mergedResults, y=flatPromoter, by="promoterIds")
	mergedResults <- dplyr::left_join(x=mergedResults, y=nestedPromoter, by=c("promoterIds", "geneIds"))

	# Return
	o <- list(test=list(promoter=testPromoter,
											gene=testGene,
											nested=testNested),
						stats=list(flatPromoter=flatPromoter,
											 flatGene=flatGene,
											 nestedPromoter=nestedPromoter,
											 nestedGene=nestedGene,
											 nestedSimes=nestedSimes),
						results=mergedResults)

	o
}


#' Perform test for APU
#'
#' Use edgeR's diffspliceDGE to perform flat DE and nested APU analyses
#'
#' @param APfitObject APfit Object.
#' @param contrast numeric or matrix: A single contrast or a contrast matrix to be tested.
#'
#' @return RETURN APtest object
#' @examples
#' # ADD EXAMPLES HERE
#' @export
APtest <- function(APfitObject, contrast){
	if(class(contrast) == "numeric"){
		o <- list(singleContrast=singleContrast(o=APfitObject, contrast=contrast))

		}else if(class(contrast) == "matrix"){
			message("Testing multiple contrasts...")
			o <- pbapply::pbapply(contrast, 2, function(x) suppressMessages(singleContrast(o=fits, contrast=x)))
	}

	# Return
	class(o) <- c("APtest", class(o))

	o
}

#' Internal function: Classify single contrast
#'
#' Called by APclassify
#'
#' @param con numeric: contrast
#' @param p numeric: FDR-corrected p-value threshold.
#' @param promoterTest See APclassify.
#' @param geneTest See APclassify.
#'
#' @import magrittr
#'
#' @return See APclassify
#' @examples
#' # ADD EXAMPLES HERE
#' @export
singleClassify <- function(con, p=0.05, promoterTest, geneTest){
	### Summarise stats based on settings
	huge <- con$results %>%
		dplyr::group_by(geneIds) %>%
		dplyr::summarise(nPromoter=unique(nPromoters),
										 # Flat info
										 geneDE=unique(flatGene_FDR < p),
										 geneDir=ifelse(unique(flatGene_logFC) >= 0, "Up", "Down"),
										 # APU info
										 geneAPUgene=unique(nestedGene_FDR) < p,
										 geneAPUsimes=unique(nestedSimes_FDR) < p,
										 geneAPUboth=unique(nestedGene_FDR) < p &
										 	unique(nestedSimes_FDR) < p,
										 # Promoter info
										 nUpFlat=sum(flatPromoter_FDR < p &
										 							flatPromoter_logFC >= 0),
										 nUpNested=sum(nestedPromoter_FDR < p &
										 								nestedPromoter_logFC >= 0),
										 nUpBoth=sum((nestedPromoter_FDR < p &
										 						 	nestedPromoter_logFC >= 0) &
										 							(flatPromoter_FDR < p &
										 							 	nestedPromoter_logFC >= 0)),
										 nDownFlat=sum(flatPromoter_FDR < p &
										 								flatPromoter_logFC <= 0),
										 nDownNested=sum(nestedPromoter_FDR < p &
										 									nestedPromoter_logFC <= 0),
										 nDownBoth=sum((nestedPromoter_FDR < p &
										 							 	nestedPromoter_logFC <= 0) &
										 								(flatPromoter_FDR < p &
										 								 	nestedPromoter_logFC <= 0)),
										 # Actual gene stats for sorting
										 geneFDR=unique(nestedGene_FDR),
										 simesFDR=unique(nestedSimes_FDR)) %>%
		dplyr::filter(!is.na(nPromoter))

	### Apply selection criteria

	# Extract data needed for all comparisons
	o <- huge %>% dplyr::select(geneIds, nPromoter, geneDE, geneDir, geneFDR, simesFDR)

	# Append correct gene test
	if(geneTest == "simes"){
		o <- dplyr::mutate(o, geneAPU=huge$geneAPUsimes)
	}else if(geneTest == "gene"){
		o <- dplyr::mutate(o, geneAPU=huge$geneAPUgene)
	}else if(geneTest == "both"){
		o <- dplyr::mutate(o, geneAPU=huge$geneAPUsimes & huge$geneAPUgene)
	}else{
		message("Warning: Something went wrong!")
	}

	# Append correct promoter test
	if(promoterTest == "nested"){
		o <- o %>%
			dplyr::mutate(nUp=huge$nUpNested,
										nDown=huge$nDownNested)
	}else if(promoterTest == "flat"){
		o <- o %>%
			dplyr::mutate(nUp=huge$nUpFlat,
										nDown=huge$nDownFlat)
	}else if(promoterTest == "both"){
		o <- o %>%
			dplyr::mutate(nUp=huge$nUpBoth,
										nDown=huge$nDownBoth)
	}else{
		message("Warning: Something went wrong!")
	}

	### Classify genes into APU types.

	# Orient direction relative to gene
	o <- o %>%
		mutate(nConcordant=ifelse(geneDir=="Up", nUp, nDown),
					 nOpposite=ifelse(geneDir=="Up", nDown, nUp))

	# Classify
	o$APUclass <- "NoAPU"

	o <- o %>%
		dplyr::mutate(APUclass="NoAPU") %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & nConcordant == 0 & nOpposite == 0 & geneDE==FALSE, "StableDiffuse", APUclass)) %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & nConcordant == 0 & nOpposite == 0 & geneDE==TRUE, "ComplexDiffuse", APUclass)) %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & (nConcordant > 0 & nOpposite == 0) & geneDE==FALSE, "WeakEmergence", APUclass)) %>%
		dplyr::	mutate(APUclass=ifelse(geneAPU == TRUE & (nConcordant > 0 & nOpposite == 0) & geneDE==TRUE, "StrongEmergence", APUclass)) %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & (nConcordant == 0 & nOpposite > 0) & geneDE==FALSE, "WeakRetention", APUclass)) %>%
		dplyr::	mutate(APUclass=ifelse(geneAPU == TRUE & (nConcordant == 0 & nOpposite > 0) & geneDE==TRUE, "StrongRetention", APUclass)) %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & nConcordant > 0 & nOpposite > 0 & geneDE==FALSE, "StableShift", APUclass)) %>%
		dplyr::mutate(APUclass=ifelse(geneAPU == TRUE & nConcordant > 0 & nOpposite > 0 & geneDE==TRUE, "ComplexShift", APUclass)) %>%
		dplyr::mutate(APUclass=factor(APUclass, levels=c("NoAPU", "StableDiffuse", "ComplexDiffuse",
																										 "WeakRetention", "StrongRetention",
																										 "WeakEmergence", "StrongEmergence",
																										 "StableShift", "ComplexShift")))
	### Return
	o
}

#' Classify genes based on types of APUs
#'
#' Classfiy multipromoter genes into 7 classes based on types of APU.
#'
#' @param APtestObject APtest object.
#' @param significanceThreshold FDR-corrected p-value threshold for significance.
#' @param promoterTest Not yet implemented.
#' @param geneTest Not yet implemented.
#'
#' @return APclassify object
#' @examples
#' # ADD EXAMPLES HERE
#' @export
APclassify <- function(APtestObject, significanceThreshold=0.05, promoterTest="nested", geneTest="simes"){
	message("Classifying APU genes...")
	o <- pbapply::pblapply(APtestObject, function(x) suppressMessages(singleClassify(con=x,
																																								 p=0.05,
																																								 promoterTest=promoterTest,
																																								 geneTest=geneTest)))

	# Return
	class(o) <- c("APclassify")#, class(o))

	o
}
