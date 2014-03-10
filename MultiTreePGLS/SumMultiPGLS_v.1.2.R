##########################
#Summary MultiTreePGLS output
##########################
#Extract summary from MultiTreePGLS output (Bayesian or ML framework)
#v1.2
#----
#guillert(at)tcd.ie & healyke(at)tcd.ie - 13/02/2013
##########################


SumMultiPGLS<-function(chain.name){
#Data checking
	if(((-is(chain.name, "character"))==0)) {stop("Chain name is not correct")} else {ok<-"ok"}
	if(length(list.files(pattern=chain.name))==0) {stop("No models corresponding to chain.name")} else {ok<-"ok"}
	load(list.files(pattern=chain.name)[[1]])
	if(((-is(model, "pgls"))==0)) {if(((-is(model, "MCMCglmm"))==0)){stop("Model format should be 'pgls' or 'MCMCglmm'")}else{Mod.typ<-"MCMCglmm"}} else {Mod.typ<-"pgls"}
	
#Summary
	if(Mod.typ=="pgls"){
	##Creating the vectors
	coef.names <- names(model$model$coef)
	Coefficient.est <- matrix(0,nrow = (length(list.files(pattern=chain.name))), ncol = (length(coef.names)), dimnames = list(c(),c(coef.names)))
	p.values <- matrix(0,nrow = (length(list.files(pattern=chain.name))), ncol = (length(coef.names)), dimnames = list(c(),c(coef.names)))
	Lambda.est <- matrix(0,nrow = (length(list.files(pattern=chain.name))), ncol = 1, dimnames = list(c(),c("Lambda.est")))
	##Summary
	for (i in 1:length(list.files(pattern=chain.name))){
		load(list.files(pattern=chain.name)[[i]])
		cat(format(Sys.time(), "%H:%M:%S"), "-", chain.name, "from tree", as.character(i), "loaded", "\n")
		for(t in 1:(length(model$model$coef[]))){
			Coefficient.est[i,t]<-rnorm(1,model$model$coef[t],model$sterr[t])
			p.values[i,t]<-summary(model)$coefficients[t,4]
			Lambda.est[i]<-model$param[2]}}
	##output				
	return<-list(Lambda.est=Lambda.est, Coefficient.est=Coefficient.est, p.values=p.values)}
	else{
	##Creating the vectors
	Gcovariances<-vector("list", length(list.files(pattern=chain.name)))
	Rcovariances <-vector("list", length(list.files(pattern=chain.name)))
	Sol<-vector("list", length(list.files(pattern=chain.name)))		
	Residuals<-vector("list", length(list.files(pattern=chain.name)))
	Phylo.term<-vector("list", length(list.files(pattern=chain.name)))

	##Summary
	for (i in 1:length(list.files(pattern=chain.name))){
	load(list.files(pattern=chain.name)[[i]])
	cat(format(Sys.time(), "%H:%M:%S"), "-", chain.name, "from tree", as.character(i), "loaded", "\n")
		Gcovariances[[i]]<-summary(model)$Gcovariances
		Rcovariances[[i]]<-summary(model)$Rcovariances		
		Sol[[i]]<-model$Sol
		Residuals[[i]]<-model$VCV[,"units"]
		Phylo.term[[i]]<-model$VCV[,"animal"]}
	##Output
	return<-list(Gcovariances=Gcovariances, Rcovariances=Rcovariances, Sol=Sol, Residuals=Residuals, Phylo.term=Phylo.term)}
return(return)}