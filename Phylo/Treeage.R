##########################
#Tree ages
##########################
#Extract the node and the tips ages of a chronogram
#v1.0
#----
#guillert(at)tcd.ie 21/09/2013
#Modified from [R-sig-phylo] nodes and taxa depth II - 21/06/2011 - Paolo Piras - ppiras(at)uniroma3.it
##########################


Treeage<-function(tree,scale=1, type='past'){

#Data checking
if(class(tree) !='phylo') {
	stop('Tree argument should be a phylo object')}
			
if(class(scale) !='numeric') {
	stop('Scaling argument should be a numerical value')} else {
		if(length(scale) !=1) {
			stop('Scaling argument should be a numerical value')}
	}

if(class(type) !='character') {
	stop('Type argument should be past or present')} else {
		if(type !='past') {
			if(type !='present') {
			stop('Type argument should be past or present')}
	}}

		
#Calculating the tips and the edges age (from Paolo Piras - 21/06/2011 )
ages.table<-data.frame(ages=as.matrix(dist.nodes(tree))[,length(tree$tip)+1],tips=c(rbind(tree$tip.label),c((length(tree$tip)+1):max(length(dist.nodes(tree)[,1])))))		

#Scaling the ages (default=1)
ages.table$ages<-ages.table$ages/max(ages.table$ages)
ages.table$ages<-ages.table$ages*scale

#Type
if(type == 'past'){
	ages.table$ages<-abs(ages.table$ages-max(ages.table$ages))}

#Output
return(ages.table)
}
