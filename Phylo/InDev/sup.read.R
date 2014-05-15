###################
#Support reading
###################

sup.read<-function(pattern,plot=TRUE){
	
warning('In developement')

	names<-list.files(pattern=pattern)
	if (length(names)==1){
		
		#Posterior probabilities from a mbtrConv.sh nexus file
		
		#Loading the function BS.read
		PP.read<-function(pattern,plot=TRUE){
		trees<-read.nexus(pattern)
		for (i in 1:length(trees)) {trees[[i]][[4]]<-as.numeric(trees[[i]][[4]])}
		names<-names(trees)
		PP.list<-as.list(names)
		
		for (i in 1:length(trees)){PP.list[[i]]<-trees[[i]][[4]]}
		names(PP.list)<-names
		
			if(plot==TRUE){
			boxplot(PP.list,las=2,ylim=c(0,1), ylab="Posterior probabilities values",main="Trees Posterior probabilities comparison")
			for (i in 1:length(PP.list)) {points(i,length(PP.list[[i]])/length(trees[[i]][[2]])-1,pch=20, col="red")}
			}
			else{
				pdf("PosteriorProbs.pdf")
				boxplot(PP.list,las=2,ylim=c(0,1),ylab="Posterior probabilities values",main="Trees Posterior probabilities comparison")
				for (i in 1:length(PP.list)) {points(i,length(PP.list[[i]])/length(trees[[i]][[2]])-1,pch=20, col="red")}
				dev.off()
			}
		
		return(PP.list)}

#Adding the resolution (number of nodes)
bla<-PP.read("Sum_trees.tre")
vect<-rep(NA,length(bla))
for (i in 1:length(bla)){vect[i]<-((length(bla[[i]])/50))}
for(i in 1:length(vect)) {points(i,vect[i], pch=20, col="red")}


		
	#Running the function BS.read		
	sup.list<-PP.read(pattern,plot)

		
	} else {
		
		#Boostrap values from multiple tree files
		
		#Loading the function BS.read
		BS.read<-function(pattern,plot=TRUE){
		names<-list.files(pattern=pattern)
		BS.list<-as.list(names)
		tree<-read.tree(list.files(pattern=pattern)[1])
		
		for (i in 1:length(BS.list)){BS.list[[i]]<-rep(NA,(tree[[2]]))}
		names(BS.list)<-names

		for (f in 1:length(list.files(pattern=pattern))){
			trees<-read.tree(list.files(pattern=pattern)[f])
				for (j in 1: tree[[2]]){
					BS.list[[f]][j]<-as.numeric(trees[[5]])[j]}
			}
	
		if(plot==TRUE){
			boxplot(BS.list,las=2,ylab="Boostrap values",main="Trees bootstrap comparison")}
			else{
				pdf("Bootstraps.pdf")
				boxplot(BS.list,las=2,ylab="Boostrap values",main="Trees bootstrap comparison")
				dev.off()
			}
		
		return(BS.list)}
		
		#Running the function BS.read
		sup.list<-BS.read(pattern,plot)
		
		}
	
	return(sup.list)

}
