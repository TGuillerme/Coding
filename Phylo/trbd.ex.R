##########################
#Builds a Birth-Death tree with the same number of extant and extinct species
##########################
#Building a Birth-Death tree by rejection sampling process using tree.bd algorithm from diversitree package (FitzJohn, 2013) with a maximum living taxa stop and allows a fixed number of associated fossils. The lambda (speciation rate) and mu (extinction rate) parameters are sampeled from a uniform distribution with lambda>mu.
#Optional arguments are:
#-error (default=0) which is a possible error margin for the generated number of extinct species to match with the given extinct number
#-rename (default=TRUE) which is the possibility of renaming the tips with three digits (from "sp0000" and "ex0000" to "sp9999" and "ex9999")
#v1.0
#----
#Thomas Guillerme - guillert(at)tcd.ie - 12/09/2013
##########################


trbd.ex<-function(extant, extinct, error=0, rename=TRUE)
{
	#Loading the packages
	require(ape)
	require(diversitree)

	#input check
	##Number of extant species should be numerical
	if(missing(extant)) {
		stop('Extant argument should be a numerical value')} else {
			if(class(extant) !='numeric') {
				stop('Extant argument should be a numerical value')} else {
					if(length(extant) !=1) {
						stop('Extant argument should be a numerical value')}
				}
		}
		
	##Number of extinct species should be numerical
	if(missing(extinct)) {
		stop('Extinct argument should be a numerical value')} else {
			if(class(extinct) !='numeric') {
				stop('Extinct argument should be a numerical value')} else {
					if(length(extinct) !=1) {
						stop('Extinct argument should be a numerical value')}
				}
		}

	##error should be a probability
	if(class(error) !='numeric') {
		stop('Error argument should be a value within the interval 0-1')} else {
			if(error >= 1.0000001) {
				stop('Error argument should be a value within the interval 0-1')} else {
					if(error <= -0.000001) {
						stop('Error argument should be a value within the interval 0-1')}
				}
		}

	##Rename should be logical
	if(class(rename) != 'logical') {
		stop('Rename argument should be logical')
	}

	#Algorithm functions: generating a birth death tree with the 3 following conditions
	##Setting the parameters lambda and mu
	fun.pars<-function(x=1){
		lambda<-runif(x)
		mu<-runif(x,0,lambda)
		return(cbind(lambda, mu))}

	##Building the birth death tree.
	##Condition 1 is that the birth death process don't fail (output = 'phylo')		
	fun1<-function(extant)
		{pars<<-fun.pars() ; trbd.tmp<<-tree.bd(pars, max.taxa= extant, include.extinct=TRUE)
			while (class(trbd.tmp) !='phylo')
			{pars<<-fun.pars() ; trbd.tmp<<-tree.bd(pars, max.taxa= extant, include.extinct=TRUE)}}
			
	###Condition 2 is that the birth death process generates the number of extant species given in the input
	fun2<-function(extant)
		{fun1(extant)
			while (length(grep('sp', trbd.tmp$tip.label)) != extant)
			{fun1(extant)}}	

	###Condition 3 is that the birth death process generates at least the number of extinct species given in the input (- error)
	fun3<-function(extant,extinct,error)
		{fun2(extant)
			while ((length(trbd.tmp$tip.label)-extant) < extinct-extinct*error)
			{fun2(extant)}}

	###When all the conditions are encountered, check if extinct=extant (+/- error), else randomly prune extinct species until extinct=extant (+/- error)
	fun4<-function(extant,extinct,error){fun3(extant, extinct, error)
		if (extinct-extinct*error <= (length(trbd.tmp$tip.label)-extant) & (length(trbd.tmp$tip.label)-extant) <= extinct+extinct*error) {
			trbd<<-trbd.tmp} else {
			trbd<<-drop.tip(trbd.tmp, c(trbd.tmp$tip.label[c(sample(grep('ex', trbd.tmp$tip.label), (length(grep('ex', trbd.tmp$tip.label))-extinct)))]))}}

	#Runing the tree
	trbd.tmp<-NULL
	trbd<-NULL
	fun4(extant, extinct, error)

	#Renaming the labels (optional, default=TRUE)
	if(rename==TRUE){

		for (i in 1:length(trbd[[3]]))
		if(nchar(trbd[[3]][i])==3){
			trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), "000", substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep="")} else {
			if(nchar(trbd[[3]][i])==4){
				trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), "00", substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep="")} else {
				if(nchar(trbd[[3]][i])==5){
					trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), "0", substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep="")}
				}
			}		
	
	}

	#Output
	return(trbd)
}