#' Plot basic statistics for a single gene
#'
#' Plots mean expression by a grouping of samples or the BCV from a GLM fit.
#'
#' @param o APfit object.
#' @param gene character: gene name.
#' @param stat character: "expression" or "BCV".
#' @param grouping factor: groups to summarize expression by.
#'
#' @return ggplot2 object
#' @examples
#' # ADD EXAMPLES HERE
#' @export
plot.APfit <- function(o, gene, stat, grouping=NULL){
	### IDEAS
	# Weighted boxplots

	### Checks
	stopifnot(gene %in% o$genes$geneIds)
	stopifnot(stat %in% c("BCV", "expression"))

	### Data used in all plots
	plottingOrder <- o$genes$promoterIds[o$genes$geneIds == gene]

	### Plot depending on selected stat
	if(stat == "BCV"){
		message(paste0(gene, ": Plotting BCVs..."))

		# Merge BCV for promoters and genes
		geneBCV <- sqrt(o$disp$gene$tagwise.dispersion[rownames(o$disp$gene$counts) == gene])

		B <- data.frame(o$genes, Promoter=sqrt(o$disp$promoter$tagwise.dispersion)) %>%
			dplyr::filter(geneIds==gene) %>%
			dplyr::mutate(Gene=geneBCV) %>%
			tidyr::gather(key="levelOfGLM", value="BCV", Promoter, Gene) %>%
			dplyr::mutate(promoterIds=factor(promoterIds, levels=plottingOrder))

		# Plots
		ggp <- ggplot2::ggplot(B, ggplot2::aes(x=promoterIds, y=BCV,
																	group=levelOfGLM, linetype=levelOfGLM)) +
			ggplot2::geom_line() +
			ggplot2::scale_linetype_manual("Level of GLM", values=c("dotted", "dashed")) +
			ggplot2::labs(x="Promoter ID", y="BCV") +
			ggplot2::theme_bw() +
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))
	}else if(stat == "expression"){
		stopifnot(!is.null(grouping))
		message(paste0(gene, ": Plotting log expression in each group..."))

		# Sample info
		sampleGrouping <- data.frame(sample=colnames(o$fit$promoter$fitted.values),
																 grouping=grouping, stringsAsFactors=FALSE)

		# Extract expression values
		E <- data.frame(o$genes, o$tpm$promoter) %>%
			dplyr::filter(geneIds==gene) %>%
			tidyr::gather(key="sample", value="logExpression", -promoterIds, -geneIds) %>%
			dplyr::left_join(y=sampleGrouping, by="sample") %>%
			dplyr::mutate(promoterIds=factor(promoterIds, levels=plottingOrder))

		# Plot
		ggp <- ggplot2::ggplot(E, ggplot2::aes(x=promoterIds, y=logExpression, fill=grouping)) +
			ggplot2::geom_boxplot(alpha=0.75) +
			ggplot2::labs(x="Promoter ID", y="log2(Expression)") +
			ggplot2::scale_fill_brewer("Grouping", palette="Set1") +
			ggplot2::theme_bw() +
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))

	}else{
		message("Something went wrong...")
		ggp <- NULL
	}

	### Return
	ggp
}

#' Plot DE statistics for a single gene
#'
#' Plot either logFCs or p-values for a single gene, at both gene and promoter levels.
#'
#' @param o APtests object
#' @param gene character: gene name
#' @param stat character: Either "logFC" or "pValue"
#'
#' @import magrittr
#'
#' @return ggplot2 object
#' @examples
#' # ADD EXAMPLES HERE
#' @export
plot.APtest <- function(o, gene, stat){
	### IDEAS
	# Allow for subselection of contrasts

	### Checks
	stopifnot(gene %in% o[[1]]$results$geneIds)
	stopifnot(stat %in% c("logFC", "pValue"))

	### Data used in all plots

	# Original order of promoters
	plottingOrder <- o[[1]]$results$promoterIds[o[[1]]$results$geneIds == gene]

	# Complete wide data
	pWide <- lapply(names(o), function(x) dplyr::filter(o[[x]]$results, geneIds==gene) %>%
										dplyr::mutate(contrast=x)) %>%
		do.call(rbind, .)

	### Select stat to be plotted
	if(stat == "logFC"){
		message(paste0(gene, ": Plotting logFCs..."))
		# Tidy up data
		pTidy <- pWide %>%
			tidyr::gather(key="levelOfDE",
										value="logFC",
										flatPromoter_logFC, nestedPromoter_logFC, flatGene_logFC) %>%
			dplyr::select(promoterIds, levelOfDE, logFC, contrast)

		# Ensure correct order of factors
		pTidy <- pTidy %>%
			dplyr::mutate(promoterIds=factor(promoterIds, levels=plottingOrder),
										levelOfDE=gsub(pattern="_logFC", replacement="", x=levelOfDE))

		# Plot using ggplot
		ggp <- ggplot2::ggplot(pTidy, ggplot2::aes(x=promoterIds, y=logFC,
									 group=paste0(contrast, levelOfDE),
									 color=contrast,
									 linetype=levelOfDE)) +
			ggplot2::geom_line() +
			ggplot2::scale_linetype_manual("Level of DE", values=c("dotted", "dashed", "solid")) +
			ggplot2::scale_color_brewer("Contrast", palette="Set1") +
			ggplot2::labs(x="Promoter ID", y="logFC") +
			ggplot2::theme_bw() +
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))

	}else if(stat == "pValue"){
		message(paste0(gene, ": Plotting -log10(pValues)..."))

		# Tidy up data
		pTidy <- pWide %>%
			tidyr::gather(key="levelOfDE",
										value="pval",
										flatPromoter_FDR, nestedPromoter_FDR, flatGene_FDR, nestedGene_FDR, nestedSimes_FDR) %>%
			dplyr::select(promoterIds, levelOfDE, pval, contrast)

		# Ensure correct order of factors
		pTidy <- pTidy %>%
			dplyr::mutate(promoterIds=factor(promoterIds, levels=plottingOrder),
										levelOfDE=gsub(pattern="_FDR", replacement="", x=levelOfDE),
										logpval=-log10(pval))


		ggp <- ggplot2::ggplot(pTidy, ggplot2::aes(x=promoterIds, y=logpval,
														 group=paste0(contrast, levelOfDE),
														 color=contrast,
														 linetype=levelOfDE)) +
			ggplot2::geom_line() +
			ggplot2::scale_linetype_manual("Level of DE", values=c("dotted", "dashed", "dotdash", "solid", "twodash")) +
			ggplot2::scale_color_brewer("Contrast", palette="Set1") +
			ggplot2::labs(x="Promoter ID", y="-log10(pValue)") +
			ggplot2::theme_bw() +
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))

	}else{
		message("Something went wrong!")
		ggp <- NULL
	}

	### Return
	ggp
}

#' Plot APU classification results
#'
#' Plots either the result of a single contrast, or compares the output of multiple contrasts.
#'
#' @param o APclassify object.
#' @param contrast Name of contrast in APclassify object.
#'
#' @return RETURN DESCRIPTION
#' @examples
#' # ADD EXAMPLES HERE
#' @export
plot.APclassify <- function(o, contrast=NULL){
	# Always chose single plot if only one contrast
	if(length(o) == 1){
		contrast <- names(o)
	}

	if(!is.null(contrast)){
		message("Plotting APU classification for a single contrast...")

		# Summarise distribution of APUs
		S <- with(o[[contrast]], table(nUp, nDown, geneDE, geneAPU)) %>%
			as.data.frame %>%
			dplyr::mutate(typeAPU=ifelse(nUp != 0 & nDown != 0, "Multidirectional", "Unidirectional")) %>%
			dplyr::mutate(typeAPU=ifelse(nUp == 0 & nDown == 0, "None", typeAPU))

		# Plot
		ggp <- ggplot2::ggplot(S, ggplot2::aes(x=nUp, y=nDown, fill=log10(Freq), label=Freq, color=typeAPU)) +
			ggplot2::geom_tile() +
			ggplot2::geom_text() +
			ggplot2::facet_grid(geneAPU~geneDE, labeller = ggplot2::label_both) +
			ggplot2::scale_fill_distiller(palette="GnBu", direction="up") +
			ggplot2::scale_color_manual(values=c("red", "blue", "black")) +
			ggplot2::scale_x_discrete(expand=c(0,0)) +
			ggplot2::scale_y_discrete(expand=c(0,0)) +
			ggplot2::ggtitle(paste0("Freq of total multipromoter genes: ", sum(S$Freq))) +
			ggplot2::theme_bw()
	}else{
		message("Plotting APU classification for multiple contrasts...")

		# Count classifications
		M <- sapply(o, function(x) table(x$APUclass))
		M <- data.frame(APUclass=factor(rownames(M),
																		levels=c("NoAPU",
																						 "StableDiffuse", "ComplexDiffuse",
																						 "WeakEmergence", "StrongEmergence",
																						 "StableShift", "ComplexShift")), M)
		M <- dplyr::filter(M, APUclass != "NoAPU")
		M <- tidyr::gather(M, key="contrast", value="count", -APUclass)

		# Plot
		ggp <- ggplot2::ggplot(M, ggplot2::aes(x=APUclass, y=count, fill=contrast)) +
			ggplot2::geom_bar(stat="identity", position="dodge") +
			ggplot2::coord_flip() +
			ggplot2::xlab("APU Class") +
			ggplot2::ylab("Count") +
			ggplot2::theme_bw()
	}

	### Return
	ggp
}
