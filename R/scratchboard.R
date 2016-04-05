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

source("http://people.binf.ku.dk/~nzl922/BasicSetup.R")
source("http://people.binf.ku.dk/~nzl922/BioC.R")

### NOTES

# Paths to data
project_dir <- file.path(run_loc(),
												 "seqdata/sandelin/projects/CERNADUC")

personal_dir <- file.path(run_loc(),
													"seqdata/sandelin/people/nzl922/ibd")

# Load promoters
SE_fname <- file.path(project_dir,
											"CAGE/FREEZE/March2015/R_files/SE_promoters.RData")

load(SE_fname)

### Subset to only active cases

# Only use active cases
SE <- subset(SE, select=condition %in% c("CDa", "con", "UCa"))
colnames(SE)[26] <- "cont091"


# Add subgroup
CD_subgroup <- c("CD1954", "CD1947", "CD1945", "CD1957", "CD1955", "CD1949", "CD1956", "CD1951", "CD1948",
								 "CD1950")

SE$subgroup <- SE$condition
SE$subgroup[str_replace(SE$sample_id, pattern="a", replacement="") %in% CD_subgroup] <- "subgroup"
SE$subgroup <- factor(SE$subgroup, levels=c("con", "UCa", "CDa", "subgroup"))

### Model matrix

# Model matrix
mod <- model.matrix(~batch_in_four + subgroup, data=colData(SE))

# Contrast matrix
con <- cbind(InflAny=c(0,0,0,0,1,1,1),
						 CDvUC=c(0,0,0,0,-1,0.5,0.5),
						 SGvCD=c(0,0,0,0,0,-1,1))

geneIds <- mcols(SE)$transcript
geneIds[geneIds == "."] <- paste0("Novel_", 1:sum(geneIds == "."))

### Actual code

fits <- APfit(promoterCounts=assay(SE, "counts"), modelMatrix=mod, geneIds=geneIds )
tests <- APtest(APfitObject=fits, contrast=con)
apus <- APclassify(APtestObject=tests)

geneName <- "PCBD2"

p1 <- plot(fits, gene=geneName, stat="expression", grouping=factor(SE$subgroup))
p2 <- plot(tests, gene=geneName, stat="logFC")
p3 <- plot(tests, gene=geneName, stat="pValue")

grid.draw(rbind(ggplotGrob(p1 + theme(axis.text.x=element_blank(),
																			axis.title.x=element_blank())),
								ggplotGrob(p2 + theme(axis.text.x=element_blank(),
																			axis.title.x=element_blank())),
								ggplotGrob(p3),
								size="last"))

