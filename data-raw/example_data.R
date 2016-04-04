library(coRe)
setup_R()

# Example data
library(parathyroidSE)
data("parathyroidExonsSE")

# Remove technical replicates
parathyroidExonsSE <- subset(parathyroidExonsSE, select=!run%in% c("SRR479078", "SRR479075", "SRR479061","SRR479064"))

# TPM trimming
parathyroidExonsSE <- subset(parathyroidExonsSE,
														 cpm(assay(parathyroidExonsSE)) %>%
														 	rowMeans %>%
														 	is_greater_than(10))

# Study design
designMatrix <- colData(parathyroidExonsSE) %>% as.data.frame
rownames(designMatrix) <- designMatrix$experiment

# Get and label count matrix
EM <- assay(parathyroidExonsSE)

geneIds <- mcols(parathyroidExonsSE)$gene_id %>% unlist

colnames(EM) <- designMatrix$experiment
rownames(EM) <- paste0(geneIds, "_", mcols(parathyroidExonsSE)$exonic_part)

exampleData <- list(designMatrix=designMatrix, EM=EM, geneIds=geneIds)

devtools::use_data(exampleData, overwrite=TRUE)
