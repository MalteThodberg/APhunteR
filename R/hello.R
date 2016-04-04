#' #' Fit models using limma
#' #'
#' #' This function will fit the flat and nesed models, and tests contrasts
#' #'
#' #' @param GR GRanges-object.
#' #' @param bed_fname Path to where file will be saved. Must end in ".bed". Defaults to tempfile().
#' #' @return Path to file.
#' #' @author Malte Thodberg
#' #' @import magrittr SummarizedExperiment dplyr
#' #' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
#' #' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' #' @export
#' fitAPModel <- function(rse, modelMatrix, contrastMatrix, groupingVariable) {
#' 	### Check input
#' 	message("Performing checks...")
#' 	# Class checks
#' 	# Check on design and contrast matrix
#' 	# Check grouping column exists.
#'
#'
#'
#' 	### Sum before fitting models
#' 	message("Summing counts over genes...")
#'
#' 	# Promoter counts
#' 	EMPromoter <- assay(rse, "counts")
#'
#' 	# Grouping variable
#' 	geneIds <- mcols(rse)[,groupingVariable] %>% as.character
#'
#' 	# Gene counts
#' 	EMGene <- EMPromoter %>%
#' 		tibble::as_data_frame() %>%
#' 		mutate(geneIds=geneIds) %>%
#' 		group_by(geneIds) %>%
#' 		summarise_each(funs(sum)) %>%
#' 		ungroup
#'
#' 	# Resave as matrix
#' 	EMGene <- data.frame(EMGene[,-1],
#' 											 row.names=EMGene$geneIds) %>%
#' 		as.matrix
#'
#' 	### Normalizing
#' 	message("Normalizing and weighting...")
#'
#' 	# Normalize using promoters
#' 	dgePromoter <- EMPromoter %>%
#' 		edgeR::DGEList() %>%
#' 		edgeR::calcNormFactors()
#'
#' 	# Transfer norm factors to genes
#' 	dgeGene <- edgeR::DGEList(EMGene)
#' 	dgeGene$samples$norm.factors <- dgePromoter$samples$norm.factors
#'
#' 	# Estimate weights
#' 	#vPromoter <- voomWithQualityWeights(dgePromoter, design=mod)
#' 	#vGene <- voomWithQualityWeights(dgeGene, design=mod)
#'
#' 	# Estimate only gene weights
#' 	vPromoter <- limma::voom(dgePromoter, design=modelMatrix)
#' 	vGene <- limma::voom(dgeGene, design=modelMatrix)
#'
#' 	### Fit models
#' 	message("Fitting flat and nested models...")
#'
#' 	# Fit flat models
#' 	fitPromoter <- vPromoter %>%
#' 		limma::lmFit(object=., design=modelMatrix) %>%
#' 		limma::contrasts.fit(fit=., contrasts=contrastMatrix)
#'
#' 	ebPromoter <- limma::eBayes(fit=fitPromoter, robust=TRUE)
#'
#' 	fitGene <- vGene %>%
#' 		limma::lmFit(object=., design=modelMatrix) %>%
#' 		limma::contrasts.fit(fit=., contrasts=contrastMatrix)
#'
#' 	ebGene <- limma::eBayes(fit=fitGene, robust=TRUE)
#'
#' 	# Fit nested model
#' 	ebNested <- limma::diffSplice(fit=fitPromoter,
#' 												 geneid=geneIds, exonid=rownames(EMPromoter),
#' 												 robust=TRUE, verbose=FALSE)
#'
#' 	### Test desired contrasts
#' 	message("Testing contrasts...")
#'
#' 	contrasts <- colnames(contrastMatrix) %>% as.list
#' 	names(contrasts) <- colnames(contrastMatrix)
#'
#' 	mergedModels <- lapply(contrasts, extractContrastRes,
#' 												 flatPromoter=ebPromoter,
#' 												 flatGene=ebGene,
#' 												 nested=ebNested)
#'
#' 	### Output
#'
#' 	o <- list(weights=list(gene=vGene, promoter=vPromoter),
#' 						eb=list(gene=ebGene, promoter=ebPromoter, nested=ebNested),
#' 						results=mergedModels)
#'
#' 	o
#'
#' 	#
#' 	# ### Batch correct
#' 	# message("Obtaining corrected expression values...")
#' 	# EM <- removeBatchEffect(x=vPromoter, design=mod[,-c(2:4)], covariates=mod[,2:4])
#' 	#
#'
#' }
#'
#' extractContrastRes <- function(conCol, flatPromoter, flatGene, nested){
#' 	# Scaffold to merge on
#' 	#scaffold <- data_frame(GeneID=geneIds, ExonID=mcols(SE)$name, feat_cat=mcols(SE)$feat_cat)
#'
#' 	scaffold <- as_data_frame(nested$genes) %>% set_colnames(c("PromoterId", "GeneId"))
#' 	scaffold <- data_frame(promoterId=as.character(nested$genes$ExonID),
#' 												 geneId=as.character(nested$genes$GeneID))# as_data_frame(nested$genes) %>% set_colnames(c("promoterId", "geneId"))
#'
#' 	# Add flat promoter
#' 	scaffold <- limma::topTable(ebPromoter,
#' 											 coef=conCol,
#' 											 number=Inf,
#' 											 sort.by="none",
#' 											 confint=TRUE) %>%
#' 		set_colnames(paste0("flatPromoter_", colnames(.))) %>%
#' 		mutate(promoterId=rownames(.)) %>%
#' 		left_join(x=scaffold, y=., by="promoterId")
#'
#' 	# Add flat promoter
#' 	scaffold <- limma::topTable(ebGene, coef=conCol, number=Inf, sort.by="none", confint=TRUE) %>%
#' 		set_colnames(paste0("flatGene_", colnames(.))) %>%
#' 		mutate(geneId=rownames(.)) %>%
#' 		left_join(x=scaffold, y=., by="geneId")
#'
#' 	# Add nested gene
#' 	simesTest <- limma::topSplice(fit=ebNested, coef=conCol, test="simes", number=Inf) %>%
#' 		set_colnames(c("geneId", "nPromoters", "simes.P.Value", "simes.FDR"))
#'
#' 	scaffold <- full_join(x=limma::topSplice(fit=ebNested, coef=conCol, test="F", number=Inf) %>%
#' 													set_colnames(c("geneId", "nPromoters", "nestedGene_F", "nestedGene_FTest.P.Value", "nestedGene_FTest.FDR")),
#' 												y=limma::topSplice(fit=ebNested, coef=conCol, test="simes", number=Inf) %>%
#' 													set_colnames(c("geneId", "nPromoters", "nestedGene_simes.P.Value", "nestedGene_simes.FDR")),
#' 												by=c("geneId", "nPromoters")) %>%
#' 		left_join(x=scaffold, y=., by="geneId")
#'
#' 	# Add nested promoter
#' 	scaffold <- limma::topSplice(fit=ebNested, coef=conCol, test="t", number=Inf) %>%
#' 		set_colnames(c("promoterId", "geneId", "nestedPromoter_logFC", "nestedPromoter_t", "nestedPromoter_P.Value", "nestedPromoter_FDR")) %>%
#' 		mutate(promoterId=as.character(promoterId)) %>%
#' 		left_join(x=scaffold, y=., by=c("geneId", "promoterId"))
#'
#' 	# Return
#' 	scaffold
#' }
#'
#' mod <- model.matrix(~patient+treatment+time, data=colData(parathyroidExonsSE))
#' con <- con <- data.frame(DPN=c(0,0,0,0,1,0,0),
#' 												 OHT=c(0,0,0,0,0,1,0),
#' 												 time=c(0,0,0,0,0,0,1)) %>% as.matrix
#'
#' res <- fitAPModel(rse=parathyroidExonsSE,
#' 									modelMatrix=mod,
#' 									contrastMatrix=con,
#' 									groupingVariable="gene_id")
#'
#' #' Fit models using limma
#' #'
#' #' This function will fit the flat and nesed models, and tests contrasts
#' #'
#' #' @param GR GRanges-object.
#' #' @param bed_fname Path to where file will be saved. Must end in ".bed". Defaults to tempfile().
#' #' @return Path to file.
#' #' @author Malte Thodberg
#' #' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
#' #' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' #' @export
#' summariseAPModel <- function(){
#' 	1
#' }
