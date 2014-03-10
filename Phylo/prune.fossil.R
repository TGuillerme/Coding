##########################
#Prune fossils on one or multiple trees
##########################
#Remove the fossils from one or multiple trees
#v.0.1
##########################
#SYNTAX :
#<phy> a Phylo object or a multiPhylo or the name of a pattern of multiPhylo files. If phy is a pattern, the output will be saved as *.pruned and the output will be verbose.
#<fossil_name> a vector containing the list of fossils to remove or a string of character containing the fossil's name pattern.
##########################
#----
#guillert(at)tcd.ie -  06/03/2014
##########################
#Requirements:
#-R 3
#-R package "ape"
##-To install the R packages:
##install.packages('ape')
##########################

prune.fossil<-function(phy, fossil_names){

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
            fossil.list=TRUE
        } else {
            fossil.list=FALSE
            fossil.pattern=fossil_names
        }
    }

#phy
    if(class(phy) =='phylo') {
        tree.set=FALSE
        phy.list=FALSE
    } else {
        if(class(phy) == 'multiPhylo') {
            tree.set=TRUE
            phy.list=FALSE
        } else {
            if(class(phy) =='character') {
                tree.set=TRUE
                phy.list=TRUE
                phy.pattern=phy
            } else {    
                stop('Phy must be a phylo or multiPhylo object or a name pattern')
            }
        }
    } 

#FUNCTIONS

    FUN.prune.fossil<-function(phy, fossil_names, tree.set, fossil.list) {
        if(tree.set==FALSE){

            #Only one tree as input
            if(fossil.list==TRUE){

                #Fossils names is provided
                dropped<-drop.tip(phy, fossil_names)
            removeeturn(dropped)
        
            } else {

                #Fossil names is not provided
                prune<-phy[[4]][grep(fossil.pattern,phy[[4]])]
                dropped<-drop.tip(phy, prune)
                return(dropped)
            } 
    
        } else {

            #Phy is a multiPhylo object
            if(fossil.list==TRUE){

                #Fossil names is provided
                dropped.trees<-lapply(phy, function (trees) drop.tip(trees, fossil_names))
                class(dropped.trees)<-'multiPhylo'
                return(dropped.trees)

            } else {

                #Fossil names is not provide
                prune<-phy[[1]][[4]][grep(fossil.pattern,phy[[1]][[4]])]
                dropped.trees<-lapply(phy, function (trees) drop.tip(trees, prune))
                class(dropped.trees)<-'multiPhylo'
                return(dropped.trees)
            }
        }
    }

#PRUNING THE TREES

    if(phy.list==FALSE){

        #Run the simple FUN.prune.fossil function
        drop<-FUN.prune.fossil(phy, fossil_names, tree.set, fossil.list)
        return(drop)
    } else {

        #Creating the file list
        phy.list<-list.files(pattern=phy.pattern)

        for (n in 1:length(phy.list)){

            #Name of the tree
            tree.name<-strsplit(phy.list[n], phy.pattern)[[1]][1]

            #Reading the tree
            phy<-read.nexus(phy.list[n])

            #Pruning the tree
            drop<-FUN.prune.fossil(phy, fossil_names, tree.set, fossil.list)
            cat(format(Sys.time(), "%H:%M:%S"), "-",tree.name,"pruned","\n")
            write.nexus(drop, file=paste(tree.name, ".pruned", phy.pattern, sep=""))
        }
    }
#end
}