library(ape)

context("collapse.support")
test_that("Example correct behaviour", {
    set.seed(42)
    ## Generating a random tree with "bootstrap values"
    tree_a <- rtree(10)
    tree_a$node.label <- sample(1:100, Nnode(tree_a))

    ## Expected errors
    


    ## Default collapsing
    tree_b <- collapse.support(tree_a)

    ## Same attributes
    expect_equal(class(tree_b), class(tree_a))
    ## Same tips numbers
    expect_equal(Ntip(tree_b), Ntip(tree_a))
    ## Number of nodes is the number of nodes in a with a bs > mean(bs)
    expect_equal(Nnode(tree_b), length(which(tree_a$node.label > mean(tree_a$node.label))))


    set.seed(0)
    ## Generating a random tree
    tree_1 <- rtree(10)

    ## Attributing random node values for each node
    node_values <- rnorm(Nnode(tree_1))

    ## Collapsing the nodes with positive values
    tree_2 <- collapse.support(tree_1, node_values, tol = 0)

    ## Graphical parameters
    par(mfrow = c(1,2))
    ## Plot both trees
    plot(tree_1, main = "uncollapsed")
    nodelabels(round(node_values, digits = 1))
    plot(tree_2, main = "collapsed")
    nodelabels(round(node_values[-which(node_values < 0)], digits = 1))

})