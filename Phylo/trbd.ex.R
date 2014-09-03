##########################
#Builds a Birth-Death tree with a given number of extant and extinct species
##########################
#Building a Birth-Death tree by rejection sampling process using tree.bd code from diversitree package (FitzJohn, 2013) with a maximum living taxa stop and allows a fixed number of associated fossils.
#v.2.1.1
#Update: redesigned the rejection sampling algorithm
#Update: added outgroup option
#Update: added a verbose option
#Update: user can now specify birth-death parameters
#Update: fixed the birth-death parametrization
#Update: fixed birth.death=NULL
##########################
#SYNTAX:
#<extant> the number of extant species
#<extinct> the number of extinct species
#<error> an error or tolerance for variation in the number of fossil (default=0).
#<birth.death> the birth death parameters. If not specified, birth death parameters are sampled from a random distribution with birth > death (greatly speeds up the function)
#<outgroup> whether to add an extra outgroup or not. The outgroup will be an extant taxa (default=FALSE)
#<rename> changing the taxa names to contain the same number of digit (default=FALSE)
#<verbose> whether to be verbose or not (default=FALSE). Verbose can also be 'system' to print the messages at a system level
##########################
#----
#guillert(at)tcd.ie - 03/09/2014
##########################
#Requirements:
#-R 3
#-R package 'ape'
#-R package 'diversitree'
#########################


trbd.ex<-function(extant, extinct, birth.death=NULL, error=0, outgroup=FALSE, rename=FALSE, verbose=FALSE)
{
#REQUIREMENTS
    require(ape)
    require(diversitree)

#DATA INPUT
    #Extant
    if(missing(extant)) {
        stop('Extant argument should be a numerical value')
    } else {
        if(class(extant) !='numeric') {
            stop('Extant argument should be a numerical value')
        } else {
            if(length(extant) !=1) {
                stop('Extant argument should be a numerical value')
            }
        }
    }
        
    #Extinct
    if(missing(extinct)) {
        stop('Extinct argument should be a numerical value')
    } else {
        if(class(extinct) !='numeric') {
            stop('Extinct argument should be a numerical value')
        } else {
            if(length(extinct) !=1) {
                stop('Extinct argument should be a numerical value')
            }
        }
    }

    #Error
    if(class(error) !='numeric') {
        stop('Error argument should be a value within the interval 0-1')
    } else {
        if(error > 1) {
            stop('Error argument should be a value within the interval 0-1')
        } else {
            if(error < 0) {
                stop('Error argument should be a value within the interval 0-1')
            }
        }
    }

    #Birth-death
    if(is.null(birth.death)) {
        generate.birth.death=TRUE
    } else {
        generate.birth.death=FALSE
        message('Birth death parameters are fixed by the user:\nThis configuration may slow down the tree generation')
        if(class(birth.death) != 'numeric') {
            stop('Birth-death must be a vector of two values')
        } else {
            if(length(birth.death) != 2) {
                stop('Birth-death must be a vector of two values')
            } else {
                if(birth.death[1] < birth.death[2]) {
                    message('Birth parameter < Death parameter:\nThis configuration may slow down the tree generation')
                }
            }
        }
    }

    #Outgroup
    if(class(outgroup) != 'logical') {
        stop('Outgroup argument should be logical')
    }

    #Rename
    if(class(rename) != 'logical') {
        stop('Rename argument should be logical')
    }

    #verbose
    if(verbose != 'system') {
        if(class(verbose) != 'logical') {
            stop('Verbose argument should be logical')
        }
    }

#FUNCTIONS

    #Generating lambda and mu parameter
    fun.parameters<-function(){
        lambda<-runif(1)
        mu<-runif(1,0,lambda)
        return(cbind(lambda, mu))
    }

    #Conditions for the rejection-sampling algorithm

    #Condition 1: the birth death process generates a tree (output = 'phylo')   
    fun.tbdcon1<-function(extant, birth.death, verbose) {
        if(is.null(birth.death)) {
            trbd.tmp<-tree.bd(fun.parameters(), max.taxa=extant, include.extinct=TRUE)
        } else {
            trbd.tmp<-tree.bd(birth.death, max.taxa=extant, include.extinct=TRUE)
        }      
        if(verbose==TRUE){
            cat('.')
        }
        if(verbose=='system'){
            message('.', appendLF=FALSE)
        }


        #Conditional loop (trbd.tmp = 'phylo')
        while (class(trbd.tmp) !='phylo') {
            if(is.null(birth.death)) {
                trbd.tmp<-tree.bd(fun.parameters(), max.taxa=extant, include.extinct=TRUE)
            } else {
                trbd.tmp<-tree.bd(birth.death, max.taxa=extant, include.extinct=TRUE)
            } 
            if(verbose==TRUE){
                cat('.')
            }
            if(verbose=='system'){
                message('.', appendLF=FALSE)
            }

        }

        return(trbd.tmp)
    }
            
    #Condition 2: the birth death process generates the right number of extant taxa given in the input
    fun.tbdcon2<-function(extant, birth.death, verbose) {
        if(is.null(birth.death)) {
            trbd.tmp<-fun.tbdcon1(extant, fun.parameters(), verbose)
        } else {
            trbd.tmp<-fun.tbdcon1(extant, birth.death, verbose)
        }
        if(verbose==TRUE){
            cat('.')
        }
        if(verbose=='system'){
            message('.', appendLF=FALSE)
        }


        #Conditional loop (length living taxa in trbd.tmp == extant)
        while (length(grep('sp', trbd.tmp$tip.label)) != extant) {
            if(is.null(birth.death)) {
                trbd.tmp<-fun.tbdcon1(extant, fun.parameters(), verbose)
            } else {
                trbd.tmp<-fun.tbdcon1(extant, birth.death, verbose)
            }
            if(verbose==TRUE){
                cat('.')
            }
            if(verbose=='system'){
                message('.', appendLF=FALSE)
            }

        }

        return(trbd.tmp)
    }   

    #Condition 3: that the birth death process generates at least the number of extinct species given in the input (- error)
    fun.tbdcon3<-function(extant, extinct, error, birth.death, verbose) {
        if(is.null(birth.death)) {
            trbd.tmp<-fun.tbdcon2(extant, fun.parameters(), verbose)
        } else {
            trbd.tmp<-fun.tbdcon2(extant, birth.death, verbose)
        }
            if(verbose==TRUE){
            cat('.')
        }
        if(verbose=='system'){
            message('.', appendLF=FALSE)
        }

        #Conditional loop (length fossil taxa in trbd.tmp < fossil - error)
        while((length(trbd.tmp$tip.label)-extant) < extinct-extinct*error) {
            if(is.null(birth.death)) {
                trbd.tmp<-fun.tbdcon2(extant, fun.parameters(), verbose)
            } else {
                trbd.tmp<-fun.tbdcon2(extant, birth.death, verbose)
            }
            if(verbose==TRUE){
                cat('.')
            }
            if(verbose=='system'){
                message('.', appendLF=FALSE)
            }

        }

        return(trbd.tmp)
    }
            

    #Applying the three conditions and when all the conditions are encountered, check if extinct=extant (+/- error), else randomly prune extinct species until extinct=extant (+/- error)
    fun.trbd.ex<-function(extant, extinct, error, birth.death, verbose){
        if(verbose==TRUE){
            cat('Generating a conditional birth death tree\n')
        }
        if(verbose=='system'){
            message('Generating a conditional birth death tree', appendLF=TRUE)
        }


        #Running through the three conditions
        if(is.null(birth.death)) {
            trbd.tmp<-fun.tbdcon3(extant, extinct, error, fun.parameters(), verbose)
        } else {
            trbd.tmp<-fun.tbdcon3(extant, extinct, error, birth.death, verbose)
        }
        #Remove extra fossil taxa
        if (extinct-extinct*error <= (length(trbd.tmp$tip.label)-extant) & (length(trbd.tmp$tip.label)-extant) <= extinct+extinct*error) {
            trbd<-trbd.tmp
        } else {
            trbd<-drop.tip(trbd.tmp, c(trbd.tmp$tip.label[c(sample(grep('ex', trbd.tmp$tip.label), (length(grep('ex', trbd.tmp$tip.label))-extinct)))]))
        }
        
        #Returning the tree
        if(verbose==TRUE){
            cat('Done\n')
        }
        if(verbose=='system'){
            message('Done', appendLF=TRUE)
        }

        return(trbd)
    }

    #Outgroup function
    fun.outgroup<-function(tree){
        #Setting the root edge as the mean tree edge length
        tree$root.edge<-mean(tree$edge.length)

        #Generating the outgroup
        outgroup<-list(edge=matrix(c(2,1),1,2),
            tip.label='sp0',
            edge.length=mean(tree$edge.length)+max(node.depth.edgelength(tree)),
            Nnode=1)
        class(outgroup)<-'phylo'

        #Combining the outgroup and the tree
        tree<-bind.tree(tree, outgroup ,position=mean(tree$edge.length))
        return(tree)
    }

    #Renaming function
    fun.rename<-function(tree){
        for (i in 1:length(tree[[3]]))
        if(nchar(tree[[3]][i])==3) {
            tree[[3]][i]<-paste(substr(tree[[3]][i],1,2), "000", substr(tree[[3]][i],3,nchar(tree[[3]][i])), sep="")
        } else {
            if(nchar(tree[[3]][i])==4) {
                tree[[3]][i]<-paste(substr(tree[[3]][i],1,2), "00", substr(tree[[3]][i],3,nchar(tree[[3]][i])), sep="")
            } else {
                if(nchar(tree[[3]][i])==5) {
                    tree[[3]][i]<-paste(substr(tree[[3]][i],1,2), "0", substr(tree[[3]][i],3,nchar(tree[[3]][i])), sep="")
                }
            }
        }   
        return(tree)
    }

#GENERATING A BIRTH DEATH TREE WITH A GIVEN NUMBER OF LIVING AND FOSSIL TAXA
    tree<-fun.trbd.ex(extant, extinct, error, birth.death, verbose)

    if(outgroup==TRUE){
        tree<-fun.outgroup(tree)
    }

    if(rename==TRUE){
        tree<-fun.rename(tree)
    }

    #Output
    return(tree)
}