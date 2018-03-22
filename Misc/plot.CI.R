#' @title Plot confidence intervals
#'
#' @description Plot confidence intervals in a "boxplot" style as lines or polygons
#'
#' @param x a list of vectors
#' @param type \code{"line"} or \code{"polygon"}
#' @param CI the confidence intervals to plot (in percentages, \code{default = c(95, 50)}) 
#' @param cent.tend the central tendency (\code{default = mean})
#' @param col.x a vector of one or more colours to apply to x (default is \code{"black"} for lines and \code{"grey"} for polygons)
#' @param ... Optional arguments to be handled by \code{\link[graphics]{polygon}, \code{\link[graphics]{lines} or \code{\link[graphics]{points} 
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

plot.CI <- function(x, type, CI = c(95, 50), cent.tend = mean, col.x, ...) {

    ## Converts one or more CI into quantile probabilities
    CI.converter <- function(CI) {
        sort(c(50-CI/2, 50+CI/2)/100)
    }

    ## Convert the quantiles
    quantiles <- CI.converter(CI)

    ## Summarising the data
    x_summary <- lapply(x, stats::quantile, probs = quantiles)

    ## Number of plots
    points_n <- length(x)
    quantiles_n <- length(CI)

    ## col.x
    if(missing(col.x)) {
        col.x <- ifelse(type == "line", "black", "grey")
    }
    if(length(col.x) != length(x)) {
        col.x <- rep(col.x, length(x))
    }

    ## Options
    dots <- list(...)

    if(is.null(dots$border)) dots$border <- "black"
    if(is.null(dots$point.col)) dots$point.col <- col.x
    if(is.null(dots$width)) dots$width <- 0.5
    if(is.null(dots$shift)) dots$shift <- 0
    if(is.null(dots$shift)) dots$lwd <- 1.5

    ## Plot polygons
    if(type == "polygon") {
        for (point in 1:points_n) {
            for(cis in 1:quantiles_n) {
                ## Setting X
                x_vals <- c(point-dots$width/(quantiles_n - cis + 1.5), point+dots$width/(quantiles_n - cis + 1.5), point+dots$width/(quantiles_n - cis + 1.5), point-dots$width/(quantiles_n - cis + 1.5)) + dots$shift

                ## Setting Y
                y_vals <- c(rep(x_summary[[point]][cis], 2), rep(x_summary[[point]][quantiles_n*2 - (cis-1)], 2))

                ## Plotting the box
                graphics::polygon(x_vals, y_vals, col = col.x[[point]], border = dots$border, dots$density)
            }
        }
    }

    if(type == "line") {
        for (point in 1:points_n) {
            for(cis in 1:quantiles_n) {

                ## Setting Y
                y_vals <- c(x_summary[[point]][cis], x_summary[[point]][quantiles_n*2 - (cis-1)])

                ## Plotting the line
                graphics::lines(x = rep((point + dots$shift), 2), y = y_vals, lty = (quantiles_n - cis + 1), lwd = cis * dots$lwd, col = col.x[[1]])

            }
        }
    }

    graphics::points(1:points_n + dots$shift, unlist(lapply(x, cent.tend)), col = dots$point.col, pch = 19)

    return(invisible())
}
