##########################
#Bind randomly trees together
##########################
#Binding randomly trees together with a provided number of trees and a root age.
#----
#guillerme.thomas(at)gmail.com - 02/01/2012
##########################

RandTreeBind<-function(x,y,sample,root.age,print.progress=T)
	{
	require(ape)
	#data checking
	##trees
	if(((-is(x, "multiPhylo"))==0)) {stop("'x' is not a multiPhylo object")} else
		{if (((-is(y, "multiPhylo"))==0)) {stop("'y' is not a multiPhylo object")} else ok<-"ok"}

	##sample
	if (((-is(sample, "numeric"))==0)) {stop("sample size is not numeric")} else
		{if (sample>(length(x)+1)) {stop("sample size is > than the trees provided in 'x'")} else
			{if (sample>(length(y)+1)) {stop("sample size is > than the trees provided in 'y'")} else ok<-"ok"}}

	##root.age
	if (((-is(sample, "numeric"))==0)) {stop("root age is not numeric")} else {ok<-"ok"}

	#sampling the trees
	randX<-sample(1:length(x), sample, replace=F)
	randY<-sample(1:length(y), sample, replace=F)
	#multiPhylo
	w<-rmtree(sample,2)

	#binding trees
	for (i in 1:sample){
		##root.age
		x[[randX[i]]]$root.edge<-root.age-diag(vcv.phylo(x[[randX[i]]]))[[i]]
		y[[randY[i]]]$root.edge<-root.age-diag(vcv.phylo(y[[randY[i]]]))[[i]]
		##binding
		w[[i]]<-x[[randX[i]]]+y[[randY[i]]]
		if (print.progress) {cat("tree", i, format(Sys.time(), "%H:%M:%S"), "\n")}
		}
	return(w)
	}
