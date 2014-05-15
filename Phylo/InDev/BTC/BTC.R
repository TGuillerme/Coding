##########################
#Calculate difference between sets of trees
##########################
#Calculate the difference between sets of trees from Bayesian distributions. Allows to use any tree difference metric function. Allows to compare the trees all to all (i.e. !n comparisons) or one to all (i.e. n comparisons to a specified reference set of trees). If one just want to compare two sets of trees, both options have the same effect.
#v.0.4
#Update: When output function is specified, each comparison is stored in an output file and then overwritten allowing a consistent RAM saving (especially useful when using many trees and/or many metrics)
#Update: The list returned by the function is now named after the tree combinations
##########################
#SYNTAX :
#<tree_sets> the Bayesian trees to be compared (should be a list of at least two multiPhylo objects). If names are provided to the trees, they will be used in the output, otherwise default names are 'tree_set1', 'tree_set2', etc...
#<metrics> the metrics to use for calculating the difference (should be at least one R function - the function should take as input just two trees x and y such as (fun(x,y)) such as 'dist.topo()' {ape}.  If names are provided to the metrics, they will be used in the output, otherwise default names are 'metric1', 'metric2', etc...
#<sample_size> number of random iterations to calculate (default = 1000)
#<verbose> be verbose (default=TRUE)
#<opt_output> any name for the writing an output file.
#<opt_reference> should be one of the tree sets specified in tree_sets. If specified, then this tree set will be used as a reference set.
##########################
#----
#guillert(at)tcd.ie - 19/02/2014
##########################
#Requirements:
#-R 3
#-R package 'ape'
##########################

BTC<-function(tree_sets, metrics, sample_size=1000, opt_reference, opt_output, verbose=TRUE){

warning('In developement')

#HEADER

#Loading the libraries
    require(ape)

#DATA INPUT

#tree_sets
    #Is tree_sets a list?
    if(class(tree_sets) !='list') {
    stop('Tree_sets must be a list of at least 2 multiPhylo objects')
    }

    #Does tree_sets contain at least 2 elements?
    if(length(tree_sets) < 2) {
    stop('Tree_sets must be a list of at least 2 multiPhylo objects')
    }

    #Are all the elements of tree_sets multiPhylo objects?
    for(n in 1:length(tree_sets)) {
        if(class(tree_sets[[n]]) !='multiPhylo') {
            stop('Tree_sets must be a list of at least 2 multiPhylo objects')
        }
    }

    #Are names provided?
    if(length(names(tree_sets)) == 0) {
        for (i in 1:length(tree_sets)) {
            names(tree_sets)[[i]]<-paste('tree_set', i, sep='')
        }
    }

#metrics
    #Is metrics a list?
    if(class(metrics) !='list') {
    stop('Metrics must be a list of at least 1 function object')
    }

    #Does metrics contain at least 1 elements?
    if(length(metrics) < 1) {
    stop('Metrics must be a list of at least 1 function object')
    }

    #Are all the elements of metrics functions objects?
    for(n in 1:length(metrics)) {
        if(class(metrics[[n]]) !='function') {
            stop('Metrics must be a list of at least 1 function object')
        }
    }

    #Are names provided?
    if(length(names(metrics)) == 0) {
        for (i in 1:length(metrics)) {
            names(metrics)[[i]]<-paste(metrics, i, sep='')
        }
    }

#sample_size
    if(class(sample_size) != 'numeric') {
        stop('sample_size must be numeric')
    }

#verbose
    if(class(verbose) != 'logical') {
        stop('verbose must be TRUE or FALSE')
    }

#method (implied by opt_reference or by length(tree_set) == 2)
    #Setting the method
    if(missing(opt_reference)) {

        #Method is one to one if there are only two tree sets provided (tree_sets[[1]] vs. tree_sets[[2]])
        if(length(tree_sets) == 2) {
            method<-'one_to_one'

        #Method is all to all if more than two tree sets are provided and the opt_reference is not provided (tree_sets[[1]] vs. tree_sets[[2]] ; tree_sets[[1]] vs. tree_sets[[n]] ; tree_sets[[2]] vs. tree_sets[[n]] ; !n times) 
        } else {
            method<-'all_to_all'
        }

    } else {

        #Opt_reference must me a multiPhylo object
        if(class(opt_reference) !='multiPhylo') {
            stop('The reference tree must be a multiPhylo object')
        
        #If opt_reference is provided and a multiPhylo object, method is all to ref (opt_reference vs. tree_sets[[1]] ; opt_reference vs. tree_sets[[2]])
        } else {
            method<-'all_to_ref'
        }
    }

#opt_output
    #Is opt_output specified?
    if(missing(opt_output)) {
        warning('Optional output was specified: no output files will be generated which can lead to RAM usage issues.')
        output<-FALSE

    #Is opt_output a string of characters?
    } else {
        output<-TRUE
        if(class(opt_output) !='character') {
            stop('If the optional output is specified, output name must be a string of characters')
        }
    }


#FUNCTIONS

#Tree difference function: calculate the difference(s) between two trees randomly selected from two sets of trees 
    FUN.single.tree.diff<-function(tree_set.x, tree_set.y, metrics){

        #Preparing the output vector (number of metrics + 2 for the selected trees)
        diff.x.y<-rep(NA, (length(metrics)+2))

        #Selecting the tree from the tree_set.x
        diff.x.y[1]<-sample(length(tree_set.x),1)
        tree.x<-tree_set.x[[diff.x.y[1]]]

        #Selecting the tree from the tree_set.y
        diff.x.y[2]<-sample(length(tree_set.y),1)
        tree.y<-tree_set.y[[diff.x.y[2]]]

        #Calculating the metrics
        for (i in 1:length(metrics)){
            diff.x.y[i+2]<-metrics[[i]](tree.x,tree.y)
        }

        return(diff.x.y)
    }

#Measure the difference between the two sets of trees
    FUN.set.tree.diff<-function(tree_set.x, tree_set.y, metrics, sample_size){

        #Build the comparison data.frame
        X.Y<-data.frame(t(replicate(sample_size, FUN.single.tree.diff(tree_set.x[[1]], tree_set.y[[1]], metrics), simplify='matrix')))

        #Renaming the columns
        colnames(X.Y)<-c(names(tree_set.x),names(tree_set.y), names(metrics))
        return(X.Y)
    }

#CALCULATE THE DIFFERENCE BETWEEN SETS OF TREES

    #One to one method
    if(method == 'one_to_one') {

        #Assigning the tree sets
        x<-1
        y<-2

        #Calculating the matrix
        Comp_one_to_one<-FUN.set.tree.diff(tree_sets[x], tree_sets[y], metrics, sample_size)

        #Save when done (conditional)
        if(output == TRUE) {
            write.csv(Comp_one_to_one, file=paste(paste(names(tree_sets[x]), names(tree_sets[y]), opt_output, sep='_'), 'csv', sep='.'))
        }

        return(Comp_one_to_one)
        cat('Analyse completed', '\n')
        cat('Use sumBTC to summarise the comparisons', '\n')

    }

    #All to all method
    if(method=='all_to_all') {

        #Renaming the reference tree set
        X.list<-seq(1,length(tree_sets),1)
        Y.list<-seq(1,length(tree_sets),1)

        #Creating the combinations list
        xl<-rep(X.list[1:length(tree_sets)], each=length(tree_sets))
        yl<-rep(c(Y.list[1:length(tree_sets)]),length(tree_sets))
        Cb_list<-data.frame(x=xl, y=yl)
        Cb_list<-Cb_list[-which(Cb_list[,1] == Cb_list[,2]),]

        #Calculating the matrix (depending on output)

        #output == FALSE
        if (output == FALSE) {
        
            Comp_list<-NULL

            for (i in 1:nrow(Cb_list)) {
                Comp_list[[i]]<-FUN.set.tree.diff(tree_sets[Cb_list[i,1]], tree_sets[Cb_list[i,2]], metrics, sample_size)

                #Renaming
                names(Comp_list)[i]<-paste(names(tree_sets[Cb_list[i,1]]),names(tree_sets[Cb_list[i,2]]), sep='_')

                #Print when done (conditional)
                if(verbose == TRUE) {
                    cat(format(Sys.time(), '%H:%M:%S'), '-', names(tree_sets[Cb_list[i,1]]), 'compared to', names(tree_sets[Cb_list[i,2]]), ';', sample_size, 'random draws', '\n')
                }
            }

            return(Comp_list)

        } else {

            #output == TRUE
            Comp_element<-NULL

            for (i in 1:nrow(Cb_list)) {
                Comp_element<-FUN.set.tree.diff(tree_sets[Cb_list[i,1]], tree_sets[Cb_list[i,2]], metrics, sample_size)

                #Saving each element
                write.csv(Comp_element, file=paste(paste(names(tree_sets[Cb_list[i,1]]), names(tree_sets[Cb_list[i,2]]), opt_output, sep='_'), 'csv', sep='.'))

                #Print when done (conditional)
                if(verbose == TRUE) {
                    cat(format(Sys.time(), '%H:%M:%S'), '-', names(tree_sets[Cb_list[i,1]]), 'compared to', names(tree_sets[Cb_list[i,2]]), ';', sample_size, 'random draws', '\n')
                }
            }

            cat('All comparisons are saved in the', opt_output, 'chain', '\n')
        }

        cat('Analyse completed', '\n')
        cat('Use sumBTC to summarise the comparisons', '\n')

    }


    #All to ref method
    if(method=='all_to_ref') {

        #Renaming the reference tree set
        ref.tree<-list(opt_reference) #is x
        names(ref.tree)<-'Ref'
        Y.list<-seq(1,length(tree_sets),1)

        #Calculating the matrix 

        #output == FALSE
        if (output == FALSE) {
        
            Comp_list<-NULL

            for (i in Y.list) {
                Comp_list[[i]]<-FUN.set.tree.diff(ref.tree, tree_sets[i], metrics, sample_size)

                #Renaming
                names(Comp_list)[i]<-paste('ref.tree',names(tree_sets[i]), sep='_')

                #Print when done (conditional)
                if(verbose == TRUE) {
                    cat(format(Sys.time(), '%H:%M:%S'), '-', names(tree_sets[i]), 'compared to', names(ref.tree), ';', sample_size, 'random draws', '\n')
                }
            }

            return(Comp_list)

        } else {

            #output == TRUE
            Comp_element<-NULL

            for (i in Y.list) {
                Comp_element<-FUN.set.tree.diff(ref.tree, tree_sets[i], metrics, sample_size)
                
                #Saving each element
                write.csv(Comp_element, file=paste(paste(names(tree_sets[i]), names(ref.tree), opt_output, sep='_'), 'csv', sep='.'))

                #Print when done (conditional)
                if(verbose == TRUE) {
                    cat(format(Sys.time(), '%H:%M:%S'), '-', names(tree_sets[i]), 'compared to', names(ref.tree), ';', sample_size, 'random draws', '\n')
                }
           }
        
            cat('All comparisons are saved in the', opt_output, 'chain', '\n')

        }

        cat('Analyse completed', '\n')
        cat('Use sumBTC to summarise the comparisons', '\n')

    }

#end

}