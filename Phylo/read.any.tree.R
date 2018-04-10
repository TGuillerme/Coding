#' @title Read any tree
#'
#' @description Reads any type of tree, whether it's a nexus or newick file
#'
#' @param file a file name specified by either a variable of mode character, or a double-quoted string.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

read.any.tree <- function(file) {

    if(scan(what = "#NEXUS", file = file, nlines = 1, quiet = TRUE) == "#NEXUS") {
        ## The tree is a true nexus
        return(read.nexus(file))
    } else {
        ## The tree is a true newick
        return(read.tree(file))
    }
}