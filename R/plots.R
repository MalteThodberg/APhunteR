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
			tbl_df %>%
			#mutate(nUp=as.numeric(as.character(nUp)), nDown=as.numeric(as.character(nDown))) %>%
			mutate(typeAPU=ifelse(nUp != 0 & nDown != 0, "Multidirectional", "Unidirectional")) %>%
			mutate(typeAPU=ifelse(nUp == 0 & nDown == 0, "None", typeAPU))

		# Plot
		ggp <- ggplot(S, aes(x=nUp, y=nDown, fill=log10(Freq), label=Freq, color=typeAPU)) +
			geom_tile() +
			geom_text() +
			facet_grid(geneAPU~geneDE, labeller = label_both) +
			scale_fill_distiller(palette="GnBu", direction="up") +
			scale_color_manual(values=c("red", "blue", "black")) +
			scale_x_discrete(expand=c(0,0)) +
			scale_y_discrete(expand=c(0,0)) +
			ggtitle(paste0("Freq of total multipromoter genes: ", sum(P1$Freq))) +
			theme_bw()
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
		ggp <- ggplot(M, aes(x=APUclass, y=count, fill=contrast)) +
			geom_bar(stat="identity") +
			coord_flip() +
			xlab("APU Class") +
			ylab("Count") +
			theme_bw()
	}

	### Return
	ggp
}
