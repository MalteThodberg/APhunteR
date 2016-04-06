#' Internal function: Plot a ggplot2 legend over multiple colums
#'
#' Called by stackPlots
#'
#' @param n integer: Number of columns.
#'
#' @return ggplot2 object
#'
#' @examples
#' # ADD EXAMPLES HERE
#' @export
stackLegends <- function(n){
	ggplot2::guides(col=ggplot2::guide_legend(ncol=n),
				 linetype=ggplot2::guide_legend(ncol=n),
				 fill=ggplot2::guide_legend(ncol=n))
}

#' Stack multiple single gene plots
#'
#' This function will appropriately stack plots by reformating legends, removing x-axis labels and aligning plot windows.
#'
#' @param ggplots list: list of ggplot2 objects.
#' @param legendColumns integer: number of columns in legend.
#'
#' @return None
#' @examples
#' # ADD EXAMPLES HERE
#' @export
stackPlots <- function(ggplots, legendColumns=2){
	### Checks
	# Whether it's a list of ggplots
	# Whether column is an integer


	# Snippet for removing x axis labels
	removeXText <- ggplot2::theme(axis.text.x=ggplot2::element_blank(),
																axis.title.x=ggplot2::element_blank())

	# Remove legends
	pl <- lapply(ggplots, function(x) x + stackLegends(legendColumns))

	# Remove x-axis on all plots
	pl[-length(pl)] <- lapply(pl[-length(pl)], function(x) x + removeXText)

	# Make plots into grobs
	pl <- lapply(pl, ggplot2::ggplotGrob)

	# Plot
	grid.draw(do.call(rbind, pl))
}
