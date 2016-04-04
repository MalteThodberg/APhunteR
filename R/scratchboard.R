# fits
# summary(fits)
#
# mod <- model.matrix(~patient+treatment+time, data=exampleData$designMatrix)
# con <- data.frame(DPN=c(0,0,0,0,1,0,0),
# 									OHT=c(0,0,0,0,0,1,0),
# 									time=c(0,0,0,0,0,0,1)) %>%
# 	as.matrix
#
# fits <- APfit(exampleData$EM, modelMatrix=mod, geneIds=exampleData$geneIds )
# tests <- APtest(APfitObject=fits, contrast=con)
#
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
