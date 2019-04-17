#' @title Phytools branch colours
#'
#' @description Get a trait colour palette as used in the phytools package tree plots
#'
#' @param traits a vector of trait values (\code{numeric}).
#' @param col.fun a function from which to generate the colours (\code{\link[grDevives]{grey.colors}} is used by default).
#' @param contrast a numeric value from which to contrast the colours. The bigger the value, the more contrasted the colours (default is \code{100} to give similar results as \code{phytools}).
#' 
#' @examples
#' set.seed(42)
#' ## A random tree
#' tree <- ape::rcoal(11)
#' ## Some random traits
#' traits <- rnorm(10)
#' 
#' ## Getting the edge colours
#' edge_cols <- phytools.plot.col(traits)
#' 
#' ## Plotting the results (and comparing them to phytools)
#' par(mfrow = c(1,2))
#' ## The phytools results
#' phytools::plotBranchbyTrait(tree, x = traits, mode = "edge", palette = "gray", main = "phytools")
#' ## The ape results
#' ape::plot.phylo(tree, edge.color = edge_cols, edge.width = 3, main = "ape")
#' 
#' @seealso
#' \code{\link[grDevives]{grey.colors}}, \code{\link[grDevives]{rainbow}}, \code{\link[grDevives]{heat.colors}}
#' 
#' @author Thomas Guillerme
#' @export

phytools.plot.col <- function(traits, col.fun = grDevices::grey.colors, contrast = 100) {
    ## Get the trait variables
    n_traits <- length(traits)
    trait_range <- range(traits)
    n_colours <- n_traits*contrast
    ## colours
    col <- rev(col.fun(n_colours))
    # col <- grey(n_colours:1/n_colours)
    ## The colours breaks
    breaks <- 0:n_colours/n_colours * diff(trait_range)+ trait_range[1] * 1e-06

    ## Select each individual colours
    color.selector <- function(one_trait, col, breaks) {
        ## Initialise colour selection
        selection <- 1
        ## Select the right colour for the one trait
        while(one_trait >= breaks[selection] && one_trait > breaks[selection+1]) {
            ## Increment the selection
            selection <- selection + 1
        }
        ## Return the corresponding colour
        return(col[selection])
    }

    ## Select the colour for each traits
    branch_colours <- sapply(traits, color.selector, col, breaks)

    return(branch_colours)
}


