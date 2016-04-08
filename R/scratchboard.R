# mod <- model.matrix(~patient+treatment+time, data=exampleData$designMatrix)
# con <- as.matrix(data.frame(DPN=c(0,0,0,0,1,0,0),
# 									OHT=c(0,0,0,0,0,1,0),
# 									time=c(0,0,0,0,0,0,1)))
#
# fits <- APfit(exampleData$EM, modelMatrix=mod, geneIds=exampleData$geneIds )
# tests <- APtest(APfitObject=fits, contrast=con)
# apus <- APclassify(tests)

#
#
#
# plot(apus)
#
# P <- singleClassify(tests$time)
#
# P1 <- with(P, table(nUp, nDown, geneDE, geneAPU)) %>%
# 	as.data.frame %>%
# 	tbl_df %>%
# 	#mutate(nUp=as.numeric(as.character(nUp)), nDown=as.numeric(as.character(nDown))) %>%
# 	mutate(typeAPU=ifelse(nUp != 0 & nDown != 0, "Multidirectional", "Unidirectional")) %>%
# 	mutate(typeAPU=ifelse(nUp == 0 & nDown == 0, "None", typeAPU))
#
# ggplot(P1, aes(x=nUp, y=nDown, fill=log10(Freq), label=Freq, color=typeAPU)) +
# 	geom_tile() +
# 	geom_text() +
# 	facet_grid(geneAPU~geneDE, labeller = label_both) +
# 	scale_fill_distiller(palette="GnBu", direction="up") +
# 	scale_color_manual(values=c("red", "blue", "black")) +
# 	scale_x_discrete(expand=c(0,0)) +
# 	scale_y_discrete(expand=c(0,0)) +
# 	ggtitle(paste0("Freq of total multipromoter genes: ", sum(P1$Freq))) +
# 	theme_bw()


###### IBD

# source("http://people.binf.ku.dk/~nzl922/BasicSetup.R")
# source("http://people.binf.ku.dk/~nzl922/BioC.R")
#
# library(APhunteR)
#
# ### NOTES
#
# # Paths to data
# project_dir <- file.path(run_loc(),
# 												 "seqdata/sandelin/projects/CERNADUC")
#
# personal_dir <- file.path(run_loc(),
# 													"seqdata/sandelin/people/nzl922/ibd")
#
# # Load promoters
# SE_fname <- file.path(project_dir,
# 											"CAGE/FREEZE/March2015/R_files/SE_promoters.RData")
#
# load(SE_fname)
#
# ### Subset to only active cases
#
# # Only use active cases
# SE <- subset(SE, select=condition %in% c("CDa", "con", "UCa"))
# colnames(SE)[26] <- "cont091"
#
#
# # Add subgroup
# CD_subgroup <- c("CD1954", "CD1947", "CD1945", "CD1957", "CD1955", "CD1949", "CD1956", "CD1951", "CD1948",
# 								 "CD1950")
#
# SE$subgroup <- SE$condition
# SE$subgroup[str_replace(SE$sample_id, pattern="a", replacement="") %in% CD_subgroup] <- "subgroup"
# SE$subgroup <- factor(SE$subgroup, levels=c("con", "UCa", "CDa", "subgroup"))
#
# ### Model matrix
#
# # Model matrix
# mod <- model.matrix(~batch_in_four + subgroup, data=colData(SE))
#
# # Contrast matrix
# con <- cbind(InflAny=c(0,0,0,0,1,1,1),
# 						 CDvUC=c(0,0,0,0,-1,0.5,0.5),
# 						 SGvCD=c(0,0,0,0,0,-1,1))
#
# geneIds <- mcols(SE)$transcript
# geneIds[geneIds == "."] <- paste0("Novel_", 1:sum(geneIds == "."))

### Actual code
#
# fits <- APfit(promoterCounts=assay(SE, "counts"), modelMatrix=mod, geneIds=geneIds )
# tests <- APtest(APfitObject=fits, contrast=con)
# apus <- APclassify(APtestObject=tests)
#
# geneName <- "HNF4G"
#
# p1 <- plot(fits, gene=geneName, stat="expression", grouping=factor(SE$subgroup))
# p2 <- plot(tests, gene=geneName, stat="logFC")
# p3 <- plot(tests, gene=geneName, stat="pValue")
#
# # Extra functions
# removeXText <- theme(axis.text.x=element_blank(), axis.title.x=element_blank())
#
# stackLegends <- function(n){
# 	guides(col=guide_legend(ncol=n),
# 				 linetype=guide_legend(ncol=n),
# 				 fill=guide_legend(ncol=n))
# }
#
# stackPlots <- function(ggplots, legendColumns=2){
# 	# Remove legends
# 	pl <- lapply(ggplots, function(x) x + stackLegends(legendColumns))
#
# 	# Remove x-axis on all plots
# 	pl[-length(pl)] <- lapply(pl[-length(pl)], function(x) x + removeXText)
#
# 	pl <- lapply(pl, ggplotGrob)
#
# 	grid.draw(do.call(rbind, pl))
# }
#
# stackPlots(list(p1, p2, p3), 2)
#
# grid.draw(rbind(ggplotGrob(p1 + removeXText + stackLegends()),
# 								ggplotGrob(p2 + removeXText + stackLegends()),
# 								ggplotGrob(p3 + stackLegends()),
# 								size="last"))
#
# grid.draw(rbind(ggplotGrob(p1 +
#
# 													 	theme(axis.text.x=element_blank(),
# 																			axis.title.x=element_blank())),
# 								ggplotGrob(p2 +
# 													 	theme(axis.text.x=element_blank(),
# 													 				axis.title.x=element_blank()) +
# 													 	guides(col=guide_legend(ncol=3), linetype=guide_legend(ncol=3), fill=guide_legend(ncol=3))),
# 								ggplotGrob(p3),
# 								size="last"))
#
# o <- tests
#
# gene <- "LPP"
#
# # Original order of promoters
# promoterOrder <- o[[1]]$results$promoterIds[o[[1]]$results$geneIds == gene]
# contrastOrder <- names(o)
#
# # Complete wide data
# pWide <- lapply(names(o), function(x) dplyr::filter(o[[x]]$results, geneIds==gene) %>%
# 									dplyr::mutate(contrast=x)) %>%
# 	do.call(rbind, .)
#
# # Tidy up data
# pTidy <- pWide %>%
# 	tidyr::gather(key="levelOfDE",
# 								value="pval",
# 								flatPromoter_FDR, nestedPromoter_FDR, flatGene_FDR, nestedGene_FDR, nestedSimes_FDR) %>%
# 	dplyr::select(promoterIds, levelOfDE, pval, contrast)
#
# # Ensure correct order of factors
# pTidy <- pTidy %>%
# 	dplyr::mutate(promoterIds=factor(promoterIds, levels=promoterOrder),
# 								contrast=factor(contrast, levels=contrastOrder),
# 								levelOfDE=gsub(pattern="_FDR", replacement="", x=levelOfDE),
# 								logpval=ifelse(levelOfDE %in% c("flatPromoter", "flatGene"), log10(pval), -log10(pval)))
#
# ggplot2::ggplot() +
# 	ggplot2::geom_bar(data=filter(pTidy, levelOfDE %in% c("flatPromoter", "nestedPromoter")),
# 										ggplot2::aes(x=promoterIds, y=logpval,
# 																 fill=contrast),
# 										stat="identity", position="dodge", alpha=0.75) +
# 	ggplot2::geom_line(data=filter(pTidy, !levelOfDE %in% c("flatPromoter", "nestedPromoter")),
# 										ggplot2::aes(x=promoterIds, y=logpval,
# 																 color=contrast, linetype=levelOfDE,
# 																 group=paste0(contrast, levelOfDE))) +
# 	ggplot2::scale_linetype_manual("Level of DE", values=c("dotted", "dashed", "dotdash")) +
# 	ggplot2::scale_color_brewer("Contrast", palette="Set1") +
# 	ggplot2::scale_fill_brewer("Contrast", palette="Set1") +
# 	ggplot2::labs(x="Promoter ID", y="-/+ log10(pValue)\nFlat / Nested") +
# 	ggplot2::geom_hline(yintercept=0, alpha=0.75) +
# 	ggplot2::theme_bw() +
# 	ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))
#
#
# ggplot2::ggplot(pTidy, ggplot2::aes(x=promoterIds, y=logpval,
# 																					 group=paste0(contrast, levelOfDE),
# 																					 color=contrast,
# 																					 linetype=levelOfDE)) +
# 	ggplot2::geom_line() +
# 	ggplot2::scale_linetype_manual("Level of DE", values=c("dotted", "dashed", "dotdash", "solid", "twodash")) +
# 	ggplot2::scale_color_brewer("Contrast", palette="Set1") +
# 	ggplot2::labs(x="Promoter ID", y="-log10(pValue)") +
# 	ggplot2::theme_bw() +
# 	ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))

