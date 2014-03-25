##########################
#Prune fossils on one or multiple trees
##########################
#Remove the fossils from one or multiple trees. The trees can be in R environment or present out of R environment under a pattern name.
#v.0.3
#Update: allow to read trees in newick or nexus file using option nexus=TRUE/FALSE
#Update: write options allows to indicate whether to write and return the results or whether the function just returns the results
#Update: fixed the nexus=FALSE option, the prune phylogeny is now returned
##########################
#SYNTAX :
#<phy> a Phylo object or a multiPhylo or the name of a pattern of multiPhylo files. If phy is a pattern, the output will be saved as *.pruned and the output will be verbose, also, no return will be given.
#<fossil_name> a vector containing the list of fossils to remove or a string of character containing the fossil's name pattern.
#<write> whether to write the results or not (default is write=FALSE). N.B. if phy is a pattern, write option is automatically set to TRUE
#<nexus> whether the trees from the pattern chain are in newick or nexus format (default is nexus=TRUE)
##########################
#----
#guillert(at)tcd.ie -  25/03/2014
##########################
#Requirements:
#-R 3
#-R package "ape"
##-To install the R packages:
##install.packages('ape')
##########################

prune.fossil<-function(phy, fossil_names, write=FALSE, nexus=TRUE){

#HEADER

#Loading the libraries
    require(ape)

#DATA INPUT

#fossil_name
    if(class(fossil_names) !='character') {
    stop('Fossil_name must be a vector containing the list or a pattern of the fossils to prune')
    } else {

        #Is fossil_name a pattern or a list?
        if(length(fossil_names) > 1 ) {
            fossil.list=FALSE
        } else {
            fossil.list=TRUE
            fossil.pattern=fossil_names
        }
    }

#phy
    if(class(phy) =='phylo') {
        phy.list=FALSE
    } else {
        if(class(phy) == 'multiPhylo') {
            phy.list=FALSE
        } else {
            if(class(phy) =='character') {
                phy.list=TRUE
                write=TRUE
                phy.pattern=phy
            } else {    
                stop('Phy must be a phylo or multiPhylo object or a name pattern')
            }
        }
    } 

#nexus
    if(class(write) != 'logical') {
        stop('write option must be logical')
    }

#nexus
    if(class(nexus) != 'logical') {
        stop('nexus option must be logical')
    }

#FUNCTIONS

    FUN.prune.fossil<-function(phy, fossil_names, tree.set, fossil.list) {
        if(tree.set==FALSE){

            #Only one tree as input
            if(fossil.list==FALSE){

                #Fossils names is provided
                dropped<-drop.tip(phy, fossil_names)
                return(dropped)
        
            } else {

                #Fossil names is not provided
                prune<-phy$tip.label[grep(fossil.pattern,phy$tip.label)]
                dropped<-drop.tip(phy, prune)
                return(dropped)
            } 
    
        } else {

            #Phy is a multiPhylo object
            if(fossil.list==FALSE){

                #Fossil names is provided
                dropped.trees<-lapply(phy, function (trees) drop.tip(trees, fossil_names))
                class(dropped.trees)<-'multiPhylo'
                return(dropped.trees)

            } else {

                #Fossil names is not provide
                prune<-phy[[1]]$tip.label[grep(fossil.pattern,phy[[1]]$tip.label)]
                dropped.trees<-lapply(phy, function (trees) drop.tip(trees, prune))
                class(dropped.trees)<-'multiPhylo'
                return(dropped.trees)
            }
        }
    }

#PRUNING THE TREES

    if(phy.list==FALSE){

        #Given phy is a phylo or multiPhylo object (i.e. phy is already in the R environment)
        if(class(phy) == 'phylo') {

            #Run the FUN.prune.fossil function on a single tree
            tree.set=FALSE
            drop<-FUN.prune.fossil(phy, fossil_names, tree.set, fossil.list)

            if(write == TRUE) {
                if(nexus == TRUE) {
                    write.nexus(drop, file="pruned_tree.tre")
                } else {
                    write.tree(drop, file="pruned_tree.tre")
                }
            }

            return(drop)
        }

        if(class(phy) == 'multiPhylo') {

            #Run the FUN.prune.fossil function on multiple trees
            tree.set=FALSE
            drop<-NULL
            for (i in 1:length(phy)) {
                drop[[i]]<-FUN.prune.fossil(phy[[i]], fossil_names, tree.set, fossil.list)
            }
            class(drop)<-'multiPhylo'
            names(drop)<-names(phy)

            if(write == TRUE) {
                if(nexus == TRUE) {
                    write.nexus(drop, file="pruned_trees.tre")
                } else {
                    write.tree(drop, file="pruned_trees.tre")
                }
            }

            return(drop)
        }

    } else {

        #Given phy is a pattern and has to be imported (i.e. phy is not in the R environment)

        #Creating the file list
        phy.list<-list.files(pattern=phy.pattern)

        for (n in 1:length(phy.list)){

            #Name of the tree
            tree.name.list<-strsplit(phy.list[n], phy.pattern)[[1]]
            tree.name.pos<-grep('[:alpha:]', tree.name.list)
            tree.name<-tree.name.list[tree.name.pos]

            #Reading the tree
            if (nexus==TRUE) {
                phy<-read.nexus(phy.list[n])
            } else {
                phy<-read.tree(phy.list[n])
            }

            #is phy a set of trees?
            if(class(phy) == 'multiPhylo') {
                tree.set=TRUE
            } else {
                tree.set=FALSE
            }

            #Pruning the tree
            drop<-FUN.prune.fossil(phy, fossil_names, tree.set, fossil.list)
            cat(format(Sys.time(), "%H:%M:%S"), "-",tree.name,"pruned","\n")
            if(nexus == TRUE) {
                write.nexus(drop, file=paste("pruned-", tree.name, ".", phy.pattern, sep=""))
            } else {
                write.tree(drop, file=paste("pruned-", tree.name, ".", phy.pattern, sep=""))
            }
        }
    }
#end
}