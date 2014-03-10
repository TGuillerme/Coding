##########################
#Binary+Rooted+Ultrametric Tree
##########################
#Make sure that the tree is dichotomous, rooted and ultrametric
#----
#guillerme.thomas(at)gmail.com - 8/11/2012
##########################

bruTree<-function(phylo,outgroup){
	{if((-is(phylo, "phylo"))==0) {stop("Input tree is not a 'phylo' object")} else	
		{if((-is(outgroup, "character"))==0) {stop("Input outgroup is not a 'character' object")} else
			ifelse(is.binary.tree(phylo), NA, phylo<-multi2di(phylo, random=F))
			ifelse(is.rooted(phylo), NA, phylo<-root(phylo, outgroup=outgroup, resolve.root=T))
			ifelse(is.ultrametric(phylo), NA, phylo<-chronoMPL(phylo))
			cat("Is the tree binary? \n")
			print(is.binary.tree(phylo))
			cat("Is the tree rooted? \n")
			print(is.rooted(phylo))
			cat("Is the tree ultrametric? \n")
			print(is.ultrametric(phylo))
			return(phylo)
			}
		}
		}
