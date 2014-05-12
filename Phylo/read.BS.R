##########################
#Bootstrap read
##########################
#Read (and plots) the distribution of bootstraps of multiple sets of trees
#v.1
#To do: make the function run for trees in R environment
##########################
#SYNTAX :
#<phy> can be either a phylo or multiPhylo object or a chain name of nexus or newick tree out of R environment.
#<plot> whether to plot the results or not (default = TRUE).
#<save> whether to save the plot or not (default = FALSE).
#<rm.na> whether to remove the NA (splits with no support value, e.g. the tree root) (default=TRUE)
##########################
#----
#guillert(at)tcd.ie - 12/05/2014
##########################
#Requirements:
#-R 3
#-R package 'ape'
#########################


#IN DEVELOPEMENT USE FIRST COMMIT



read.BS<-function(phy, plot=TRUE, save=FALSE, rm.na=TRUE){

	#DATA INPUT

	#phy
	if(class(phy) ==  "phylo"){
		single.tree <- TRUE
		pattern.tree <- FALSE
	} else {
		if(class(phy) == "multiPhylo"){
			single.tree <- FALSE
			pattern.tree <- FALSE
		} else {
			if(class(phy) != "character"){
				stop("Input must be a 'phylo' or 'multiPhylo' object or a pattern of filename(s) containing 'phylo' or 'multiPhylo' object(s).")
			} else {
				pattern.tree <- TRUE
			}
		}
	}
	if(pattern.tree== TRUE) {
    	pattern.list<-list.files(pattern=phy)
	}

	#Import the trees if pattern.tree==TRUE
    if(pattern.tree == TRUE) {

        #Set the number of trees
        if(length(pattern.list) == 1){
            single.tree <- TRUE
        } else {
            if(length(pattern.list) == 0) {
                stop('Tree pattern not found')
            } else {
                single.tree <- FALSE
            }
        }

        #Select the tree format
        tree.format<-readChar(pattern.list[1],10)
        tree.format<-grep('#NEXUS', tree.format, ignore.case=TRUE)
        if(length(tree.format) == 0){
            tree.nexus<-FALSE
        } else {
            tree.nexus<-TRUE
        }
    }


    #plot
    if(class(plot) != 'logical'){
        stop('Plot should be TRUE or FALSE')
    }

    #save
    if(class(save) != 'logical'){
        stop('Save should be TRUE or FALSE')
    }

    #rm.na
    if(class(rm.na) != 'logical'){
        stop('Save should be TRUE or FALSE')
    }

	#FUNCTION

    #Reading phylogenies from a given pattern
    FUN.read.pat.tree<-function(pattern.list, single.tree, tree.nexus) {
        tree.names<-pattern.list
        pat.tree<-NULL

        if(single.tree == TRUE) {

            #phylo object
            if(tree.nexus == FALSE) {
                pat.tree<-read.tree(pattern.list, tree.names=tree.names)
            } else {
                pat.tree<-read.nexus(pattern.list, tree.names=tree.names)
            }

        } else {

            #multiPhylo object
            pat.tree<-as.list(tree.names)

            for (f in 1:length(tree.names)) {
                if(tree.nexus == FALSE) {
                    pat.tree[[f]]<-read.tree(pattern.list[f])
                } else {
                    pat.tree[[f]]<-read.nexus(pattern.list[f])
                }
                names(pat.tree)<-tree.names
                class(pat.tree)<-'multiPhylo'
            }
        }
        return(pat.tree)
    }

	FUN.read.BS<-function(phy, single.tree, rm.na) {
		#phy=pattern

        #add names to multiPhylo object if none provided
        if(single.tree == FALSE) {
            if(length(names(phy)) == 0) {
                tree.names<-seq(1:length(phy))
                names(phy)<-tree.names
            }
        
            #Create the BS.list
            BS.list<-as.list(names(phy))
            for (i in 1:length(BS.list)) {
                BS.list[[i]]<-NA
            }
            names(BS.list)<-names(phy)
        
            #Extract the Bootstraps
            for(i in 1:length(phy)) {

                if(length(phy[[i]]) == 4) {
                    warning('Support values must in in fifth position of the phylo/multiPhylo object (default)')
                    stop('No support values found in the given tree')
                }

                if(rm.na == FALSE) {
                    BS.list[[i]]<-as.numeric(phy[[i]][[5]])
                } else {
                    BS.list[[i]]<-as.numeric(phy[[i]][[5]])[-c(which(is.na(as.numeric(phy[[i]][[5]]))))]
                }
            }

        } else {
            #only one tree
            if(length(phy) == 4) {
                warning('Support values must in in fifth position of the phylo/multiPhylo object (default)')
                stop('No support values found in the given tree')
            }

            if(rm.na == FALSE) {
                BS.list<-as.numeric(phy[[5]])
            } else {
                BS.list<-as.numeric(phy[[5]])[-c(which(is.na(as.numeric(phy[[5]]))))]
            }
        }

        return(BS.list)
    }

	#BOOSTRAP READ

    #Load the trees if necessary
    if(pattern.tree == TRUE) {
        phy<-FUN.read.pat.tree(pattern.list, single.tree, tree.nexus)
    }

    #Read the bootstrap values from the tree(s)
	BS.list<-FUN.read.BS(phy, single.tree, rm.na)

    #plot
	if(plot==TRUE){
		boxplot(BS.list,las=2,ylab="Boostrap")}

    #save
	if(save==TRUE){
		pdf("Bootstraps.pdf")
		boxplot(BS.list,las=2,ylab="Boostrap")
		dev.off()
	}

	#Output
	BSlist<-list(Boostraps=unlist(BS.list, use.names=FALSE), Details=BS.list)	
	return(BSlist)

}
#End
