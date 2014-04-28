##########################
#Bootstrap read
##########################
#Read (and plots) the distribution of bootstraps of multiple sets of trees
#v.0.1
#To do: make the function run for trees in R environment
##########################
#SYNTAX :
#<pattern> the chain name of the trees to read.
#<plot> whether to plot the results or not (default = TRUE).
##########################
#----
#guillert(at)tcd.ie - 28/04/2014
##########################
#Requirements:
#-R 3
#########################



read.BS<-function(pattern,plot=TRUE){
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
	
	BSlist<-list(Boostraps=unlist(BS.list, use.names=FALSE), Details=BS.list)	
	return(BSlist)
	
	}