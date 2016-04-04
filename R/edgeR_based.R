APfit <- function(promoterCounts, modelMatrix, geneIds){
		### Check input
		message("Performing checks...")
		# Class checks
		# Check on design and contrast matrix
		# Check grouping column exists.

		### Gene strucutre
		geneStructure <- data_frame(promoterIds=rownames(promoterCounts),
																geneIds=geneIds)

		### Sum before fitting models
		message("Summing counts over genes...")

		# Gene counts
		geneCounts <- data.frame(geneIds=geneIds, promoterCounts) %>%
			group_by(geneIds) %>%
			summarise_each(funs(sum)) %>%
			ungroup

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
							disp=list(promoter=dispPromoter, gene=dispGene),
							fit=list(promoter=fitPromoter, gene=fitGene))

		class(o) <- c("APfit", class(o))

		o
}

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

	mergedResults <- left_join(x=mergedResults, y=flatGene, by=c("geneIds"))
	mergedResults <- left_join(x=mergedResults, y=nestedGene, by=c("geneIds"))
	mergedResults <- left_join(x=mergedResults, y=nestedSimes, by=c("geneIds", "nPromoters"))

	mergedResults <- left_join(x=mergedResults, y=flatPromoter, by="promoterIds")
	mergedResults <- left_join(x=mergedResults, y=nestedPromoter, by=c("promoterIds", "geneIds"))

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

singleClassify <- function(con, p=0.05, promoterTest="both", geneTest="both"){
	### Summarise stats based on settings
	if(promoterTest == "both" & geneTest == "both"){
		o <- con$results %>%
			group_by(geneIds) %>%
			summarise(nPromoter=unique(nPromoters),
								geneDE=unique(flatGene_FDR < p),
								geneDir=ifelse(unique(flatGene_logFC) >= 0, "Up", "Down"),
								geneAPU=unique(nestedGene_FDR) < p &
									unique(nestedSimes_FDR) < p,
								nUp=sum(nestedPromoter_FDR < p &
													nestedPromoter_logFC >= 0 &
													flatPromoter_FDR < p &
													nestedPromoter_logFC >= 0),
								nDown=sum(nestedPromoter_FDR < p &
														nestedPromoter_logFC <= 0 &
														flatPromoter_FDR < p &
														nestedPromoter_logFC <= 0)) %>%
			filter(!is.na(nPromoter))

	}else{
		warning("NOT YET IMPLEMENTED")
	}

	### Classify
	o$APUclass <- "NoAPU"

	o <- o %>%
		mutate(APUclass=ifelse(geneAPU == TRUE, "WeakDiffuse", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & nUp == 0 & nDown == 0 & geneDE==FALSE, "StableDiffuse", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & nUp == 0 & nDown == 0 & geneDE==TRUE, "ComplexDiffuse", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & (nUp > 0 | nDown > 0) & geneDE==FALSE, "WeakEmergence", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & (nUp > 0 | nDown > 0) & geneDE==TRUE, "StrongEmergence", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & nUp > 0 & nDown > 0 & geneDE==FALSE, "StableShift", APUclass)) %>%
		mutate(APUclass=ifelse(geneAPU == TRUE & nUp > 0 & nDown > 0 & geneDE==TRUE, "ComplexShift", APUclass)) %>%
		mutate(APUclass=factor(APUclass, levels=c("NoAPU", "StableDiffuse", "ComplexDiffuse", "WeakEmergence", "StrongEmergence", "StableShift", "ComplexShift")))

	### Return
	o
}

APclassify <- function(APtestObject, significanceThreshold=0.05, promoterTest="both", geneTest="both"){
	message("Classifying APU genes...")
	o <- pbapply::pblapply(APtestObject, function(x) suppressMessages(singleClassify(con=x,
																																								 p=0.05,
																																								 promoterTest="both",
																																								 geneTest="both")))

	# Return
	class(o) <- c("APclassify")#, class(o))

	o
}
