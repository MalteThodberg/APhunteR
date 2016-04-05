#' #' @export
#' print.APfit <- function(o){
#' 	message("APfit object from the the APhunteR package")
#' }

#' Print basic info about APfit
#'
#' Short representation fo APfit object
#'
#' @param o APfit object
#'
#' @examples
#' # ADD EXAMPLES HERE
#' @export
print.APfit <- function(o){
	# Calculate some basic statistics
	nPromoter <- nrow(fits$genes)
	nGenes <- length(unique(fits$genes$geneIds))

	# Stats on distribution of promoters
	ns <- table(o$genes$geneIds)
	singlePromoter <- sum(ns == 1)
	maxN <- max(ns)
	meanN <- mean(ns)
	medianN <- median(ns)

	# Message output
	message("Alternative Promoter Usage (APU) model:")
	message("=======================================")
	message(paste0("Number of promoters: ", nPromoter))
	message(paste0("Number of genes: ", nGenes))
	message(paste0("Single promoter genes: ", singlePromoter))
	message(paste0("Multi promoter genes: ", nGenes-singlePromoter))
	message(paste0("Highest promoter count: ", maxN))
	message(paste0("Mean promoter count: ", round(meanN, digits=2)))
	message(paste0("Median promoter count: ", medianN))
	message("=======================================")
}

#' #' @export
#' print.APtest <- function(o){
#' 	message("APtest object from the the APhunteR package")
#' }

#' Print basic info about APtest
#'
#' Short representation fo APtest object
#'
#' @param o APtest object
#'
#' @examples
#' # ADD EXAMPLES HERE
#' @export
print.APtest <- function(o){
	# Calculate some basic statistics
	nContrast <- length(o)

	# Message output
	message("Alternative Promoter Usage (APU) tests:")
	message(paste0("Number of contrasts: ", nContrast))
}

#' Print basic info about APclassify
#'
#' Short representation fo APclassify object
#'
#' @param o APClassify object
#'
#' @examples
#' # ADD EXAMPLES HERE
#' @export
print.APclassify <- function(o){
	# Calculate some basic statistics
	nContrast <- length(o)

	# Message output
	message("Alternative Promoter Usage (APU) classifications:")
	message(paste0("Number of contrasts: ", nContrast))
}
