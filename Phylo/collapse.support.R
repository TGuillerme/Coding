#' @title Collapse tree support
#'
#' @description Collapses nodes if support is below a threshold value
#'
#' @param tree a \code{phylo} object
#' @param values a \code{vector} of support to consider. If left missing, \code{tree$node.label} is used.
#' @param tol a \code{numeric} value giving the tolerance below which nodes are collapsed. By default, the mean of \code{values} is used.
#' 
#' @details
#' Note that if the root has a value below the tolerance, it will \textit{not} be collapsed.
#' 
#' @examples
#' ## Generating a random tree
#' tree_a <- rtree(10)
#' 
#' ## Attributing random bootstrap values
#' tree_a$node.label <- sample(1:100, Nnode(tree_a))
#' 
#' ## Default collapsing
#' ## (using the mean bootstrap value as default tolerance)
#' tree_b <- collapse.support(tree_a)
#' 
#' ## Graphical parameters
#' par(mfrow = c(1,2))
#' ## Plot both trees
#' plot(tree_a, main = "uncollapsed")
#' nodelabels(tree_a$node.label)
#' plot(tree_b, main = "collapsed")
#' nodelabels(tree_b$node.label)
#' 
#' set.seed(0)
#' ## Generating a random tree
#' tree_1 <- rtree(10)
#' 
#' ## Attributing random node values for each node
#' node_values <- rnorm(Nnode(tree_1))
#' 
#' ## Collapsing the nodes with positive values
#' tree_2 <- collapse.support(tree_1, node_values, tol = 0)
#' 
#' ## Graphical parameters
#' par(mfrow = c(1,2))
#' ## Plot both trees
#' plot(tree_1, main = "uncollapsed")
#' nodelabels(round(node_values, digits = 1))
#' plot(tree_2, main = "collapsed")
#' nodelabels(round(node_values[-which(node_values < 0)], digits = 1))
#' 
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

collapse.support <- function(tree, values, tol) {

    if(missing(values)) {
        ## Use the node labels (e.g. bs values)
        values <- tree$node.label
        ## Convert into numeric values
        if(class(values) != "numeric") {
            values <- as.numeric(values)
            ## Replace NAs by 100 bs values if the values look like bs (between 1 and 100)
            if(any(is.na(values)) && range(na.omit(values)) <= 100 && range(na.omit(values)) >= 1) {
                values[which(is.na(values))] <- 100
            }
        }
    }

    if(missing(tol)) {
        ## Define the default tolerance
        tol <- mean(values)
    }

    ## Get the nodes to remove
    remove_nodes <- which(values < tol) + Ntip(tree)

    ## Checking if a node is the root (can't collapse the root!)
    remove_root <- which(remove_nodes == Ntip(tree) + 1)
    if(length(remove_root) != 0) {
        remove_nodes <- remove_nodes[-remove_root]
    }

    ## Function for removing a node from a tree (creating a polytomy)
    remove.node <- function(tree, node) {

        ## Find the node row in the edge table
        node_row <- which(tree$edge[,2] == node)

        ## Get the parent node
        parent_node <- tree$edge[node_row, 1]

        ## Get the parent_node/node edge length
        node_edge_length <- tree$edge.length[node_row]

        ## Add the node edge length to the node descendants
        tree$edge.length[which(tree$edge[,1] == node)] <- tree$edge.length[which(tree$edge[,1] == node)] + node_edge_length

        ## Connect the nodes descendants to the parent_node
        tree$edge[which(tree$edge[,1] == node), 1] <- parent_node

        ## Remove the edge between the parent_node and the node
        tree$edge <- tree$edge[-node_row, ]
        tree$edge.length <- tree$edge.length[-node_row]

        ## Updating the edge matrix (renaming nodes and tips)
        tree$edge[which(tree$edge > node)] <- tree$edge[which(tree$edge > node)] - 1

        ## Remove the node labels (if they exist)
        if(!is.null(tree$node.label)) {
            tree$node.label <- tree$node.label[-(node-Ntip(tree))]
        }

        ## Updating the number of nodes
        tree$Nnode <- tree$Nnode - 1

        ## Return the collapsed tree
        return(tree)
    }

    ## Recursive function wrapper to loop through the nodes
    recursive.remove.node <- function(remove_nodes, tree, counter) {

        ## Remove the fist node of the list
        tree <- remove.node(tree, remove_nodes[1])

        ## Update the list (removing the first element and decrementing the node numbers)
        remove_nodes <- remove_nodes[-1]
        remove_nodes <- remove_nodes - 1

        ## Decrementing the counter
        counter <- counter - 1 

        ## If the counter is 0, return the tree, else, continue the recursion
        if(counter == 0) {
            return(tree)
        } else {
            recursive.remove.node(remove_nodes, tree, counter)
        }
    }

    return(recursive.remove.node(remove_nodes, tree, counter = length(remove_nodes)))
}
