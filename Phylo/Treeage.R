##########################
#Tree ages
##########################
#Extract the node and the tips ages of a chronogram
#v1.1
#Update: Syntax and typos
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<scale> the scale of the tree (i.e. the age of the root)
#<type> either 'past' if the tree must be a dated tree (units = time to past; tips=0), or 'present' if the tree is an absolute dated tree (units = time to present; root=0)
##########################
#----
#guillert(at)tcd.ie 30/06/2014
#Modified from [R-sig-phylo] nodes and taxa depth II - 21/06/2011 - Paolo Piras - ppiras(at)uniroma3.it
##########################


Treeage<-function(tree, scale=1, type='past'){

#INPUT

	#tree
	if(class(tree) !='phylo') {
		stop('Tree argument should be a phylo object')
	}

	#ultrametric
	ultrametric<-is.ultrametric(tree)

	#scale
	if(class(scale) !='numeric') {
		stop('Scaling argument should be a numerical value')
	} else {
		if(length(scale) !=1) {
			stop('Scaling argument should be a numerical value')
		}
	}

	#type
	if(class(type) !='character') {
		stop('Type argument should be past or present')
	} else {
			if(type !='past') {
				if(type !='present') {
				stop('Type argument should be past or present')
			}
		}
	}

#FUNCTIONS

	#Calculating the tips and the edges age
	FUN.ages.table<-function(tree){
		ages<-dist.nodes(tree)[length(tree$tip.label)+1,]
		edges<-c(rbind(tree$tip.label),'root',c((length(tree$tip.label)+2):length(dist.nodes(tree)[,1])))
		ages.table<-data.frame(ages=ages,edges=edges)
		return(ages.table)
	}	

	#Scaling the ages (default=1)
	FUN.ages.scale<-function(ages.table, scale){
		ages.table$ages<-ages.table$ages/max(ages.table$ages)
		ages.table$ages<-ages.table$ages*scale
		return(ages.table)
	}

#CALCULATE THE EDGES AGE

	if(scale == 0) {
		ages.table<-FUN.ages.table(tree)
	} else {
		ages.table<-FUN.ages.scale(FUN.ages.table(tree), scale)
	}

	#Type
	if(type == 'past'){
		tree.height<-max(ages.table$ages)
		ages.table$ages<-round(abs(ages.table$ages-tree.height), digit=7)
	} else {
		ages.table$ages<-round(ages.table$ages, digit=7)
	}

	#Output
	return(ages.table)

	#Example
	example=FALSE
	if(example == TRUE){
		#Examples
		##Generate a birth-death tree with fossil and living species
		library(diversitree)
		tree<-tree.bd(c(1,0.3), max.taxa=20, include.extinct=TRUE)

		##Calculate the edges age by setting the root a 65 Mya
		Treeage(tree, scale=65)

		##Ploting the distribution of the node ages
		hist(ages[-c(1:length(tree$tip.label)),1], xlab="Time", main="Divergence frequency per time")

		##Calculate when the fossil went extinct
		ages.fossil<-Treeage(tree, scale=65, type='past')
		tree.fossil<-tree
		tree.fossil$tip.label[grep("ex",ages.fossil[,2])]<-paste(tree.fossil$tip.label[grep("ex",ages.fossil[,2])],round(ages.fossil[grep("ex",ages.fossil[,2]),1], digit=2), "Mya", sep=" ")
		plot(tree.fossil)

		##Ploting the node age from root
		ages.nodes<-Treeage(tree, scale=65, type='present') #change 'present' into 'past' to plot a classical chronogram
		tree.nodes<-tree
		tree.nodes$node.label<-round(ages.nodes[-c(1:length(tree.nodes$tip.label)),1], digit=2)
		plot(tree.nodes)
		nodelabels(tree.nodes$node.label, cex=0.6)
	}
}
