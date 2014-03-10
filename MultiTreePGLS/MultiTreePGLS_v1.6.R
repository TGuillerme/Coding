##########################
#Run PGLS on multiPhylo object
##########################
#Running a given PGLS or MCMCglmm (Baysian) model on a set of given trees. Results are saved in a list of model files. outputs based on modeling framework with baysian approach allowing posterior distribution while PGLS only allows subsampling of mmodel output (currently used fro MCMCglmm).
#v1.6
##########################
#SYNTAX :
#<data.frame> an object of class data.frame
#<species.col> any string of characters that is the name of the column containing the species names in the data.frame object
#<formula> an object of class formula
#<trees> an object of class multiPhylo
#<output> any string of characters that will be used as chain name for the models output
#<framework> should be "Bayesian" for a Bayesian analysis or "ML" for a Maximum Likelihood analysis
#<opt_ngen> an optional numerical value that is the maximum number of generations for the MCMC - this option is ignored if framework="ML"
#<opt_converge> an optional numerical value that is the convergence value for the MCMC - this option is ignored if framework="ML"
#<opt_thin> an optional numerical value that is the thinning interval for the MCMC - this option is ignored if framework="ML"
##########################
#----
#guillert(at)tcd.ie & healyke(at)tcd.ie - 20/01/2014
##########################
#Requirements:
#-R 3
#-R package "ape"
#-R package "caper"
#-R package "mvtnorm"
#-R package "MCMCglmm"
#-R package "modeest"
#-R package "coda"
##-To install the R packages:
##install.packages(c('ape','caper','mvtnorm', 'MCMCglmm', 'modeest'))
##########################


MultiTreePGLS<-function(data.frame, species.col, formula, trees, output="Model1",framework="Bayesian", opt_ngen, opt_converge, opt_thin)
{
#Step 1: HEADER

	#Loading the libraries
	require(ape)
	require(caper)
	require(mvtnorm)
	require(coda)

#Step 2: DATA INPUT

	#Data input checking

		#framework

			#checks if framework input is chain of characters; either "ML" for PGLS or "Bayesian" for MCMCglmm	
			if(class(framework) !='character') {
				stop("framework should be 'ML' or 'Bayesian'")
			}

			#Is the framework "ML" or "Bayesian"
			if(framework=="ML") {
				ok<-"ok"} else {
					if(framework=="Bayesian") {
						require(MCMCglmm)
						require(modeest) } else {
							stop("framework should be 'ML' or 'Bayesian'")
					}
			}
		#species.col

			#Is species column a character string?
			if(class(species.col) !='character') {
				stop("species.col does not match in data.frame") 
			}

		#data.frame
	
			#Is the data.frame a data.frame object?
			if(class(data.frame) !='data.frame') {
				stop("'data.frame' is not a data.frame object")
			}

			#Does species column match with data.frame	
			if(any((names(data.frame)) == species.col) == FALSE) {
				stop("species.col does not match in data.frame")
			}

			#Modifying the data.frame object to match with the next steps

				#Changing the species column name
				names(data.frame)<-sub(species.col,"sp.col",names(data.frame))
				data.frame["animal"]<-NA
				data.frame$animal<-data.frame$sp.col

		#formula

			#Is formula of class formula?
			if(class(formula) !='formula') {
				stop("'formula' is not a formula object")
			}

		#trees

			#Are the trees of class multiPhylo
			if(class(trees) !='multiPhylo') {
				stop("trees are not a multiPhylo object")
			}

			#Are the species in the trees matching with the species in the data.frame?

				#species names in trees and in the data.frame...
				sp.tree<-lapply(trees, function (x) x$tip.label)
				Sp.tree<-lapply(sp.tree, function (x) sort(x))
				sp.data<-sort(data.frame$sp.col)

				#...are they matching?
				trees<-lapply(trees, function (x) if(length(which(is.na(match(x$tip.label, sp.data))))==0) {x<-x} else {
				x<-drop.tip(x, which(is.na(match(x$tip.label, as.character(data.frame$sp.col)))))
				cat("Trees tip label are now corresponding to species column")})
				class(trees)<-"multiPhylo"	

		#output

			#Is output a character string?
			if(class(output) !='character') {
				warning("Output name is not a string of characters")
			}

			#creating the file names were to save each tree [i]
			file.names<-vector("list", length(trees))
			for (i in (1:length(trees))){
				file.names[[i]]<-paste(output, as.character("tree"),as.character(i), sep="_")
				file.names[[i]]<-paste(file.names[[i]], as.character("R"), sep=".")
			}
		
		#optionals
		
			#opt_ngen
			if(missing(opt_ngen)) {
				ngen<-500000} else {
					if(is(opt_ngen, "numeric")) {
						ngen<-opt_ngen} else {
							if(framework=="Bayesian") {
								stop("opt_ngen is not numeric")
							}
					}
			}
	
			#opt_converge
			if(missing(opt_converge)) {
				converge<-1.1} else {
					if(is(opt_converge, "numeric")) {
						converge<-opt_converge } else {
							if(framework=="Bayesian") {
								stop("opt_converge is not numeric")
							}
					}
			}
	
			#opt_thin
			if(missing(opt_thin)) {
				opt_thin<-100 } else {
					if(is(opt_thin, "numeric")) {
						thinn<-opt_thin} else {
							if(framework=="Bayesian") {
								stop("opt_thin is not numeric")
							}
					}
			}

#Step 3: RUNNING THE MODELS
		
	#running the models (overwriting)

		#Maximum likelihood framework.
		if(framework=='ML') {
			for (i in (1:length(trees))) {

				#Model running using pgls function (caper) on each tree [i]
				model<-pgls(formula, comparative.data(data=data.frame, phy=trees[[i]], names.col="sp.col", vcv=TRUE, vcv.dim=3), lambda='ML', bounds=list(lambda=c(1e-03,1),kappa=c(1e-06,3),delta=c(1e-06,3)))
				
				#Saving the model ran on each tree [i]
				save(model, file=as.character(file.names[[i]]))

				#Printing the time
				cat(format(Sys.time(), "%H:%M:%S"), "-", output, "on tree", as.character(i), "done", "\n")
			}
		}
	
		#Bayesian framework
			#set up uniformative prior. see Jarrod Hadfield's notes on prior parametrers used. 
		else { prior<-list(R = list(V = 1/2, nu=0.002), G = list(G1=list(V = 1/2, nu=0.002)))
			for (i in (1:length(trees))){

				#Model running using MCMCglmm function (MCMCglmm) on each tree [i] on two independent chains

					#Chain 1
					model<-MCMCglmm(formula, random=~animal,pedigree=trees[[i]],prior=prior,data=comparative.data(data=data.frame, phy=trees[[i]], names.col="sp.col", vcv=FALSE)$data,verbose=FALSE,family=c("gaussian"),nitt=ngen,burnin=as.integer(ngen/6),thin=thinn)
				
					#Chain 2
					model.1<-MCMCglmm(formula, random=~animal,pedigree=trees[[i]],prior=prior,data=comparative.data(data=data.frame, phy=trees[[i]], names.col="sp.col", vcv=FALSE)$data,verbose=FALSE,family=c("gaussian"),nitt=ngen/50,burnin=as.integer(ngen/300),thin=thinn/50)
					
					#Convergence check using Gelman and Rubins diagnoses set to return true or false based on level of scale reduction set (default == 1.1)
					convergence<-gelman.diag(mcmc.list(as.mcmc(model.1$Sol[1:(length(model.1$Sol[,1])),]),as.mcmc(model$Sol[1:(length(model.1$Sol[,1])),])))

				#Saving the model ran on each tree [i]
				save(model, file=as.character(file.names[[i]]))

				#Printing the time and the ESS + convergence diagnosis
				cat(format(Sys.time(), "%H:%M:%S"), "-", output, "on tree", as.character(i), "done", "\n")
				cat(format(Sys.time(), "%H:%M:%S"), "-", "Effective sample size is >1000:",all(effectiveSize(model$Sol[])>1000), "\n")
				cat(format(Sys.time(), "%H:%M:%S"), "-", "All levels converged:",all(convergence$psrf[,1]<converge), "\n")
			}
		}

#Step 4: OUTPUT						
	
	#Output
	cat(format(Sys.time(), "%H:%M:%S"), "-", output, "run and saved on all trees","\n")
	cat("Use SumMultiPGLS function to summarise the output", "\n")

#End
}