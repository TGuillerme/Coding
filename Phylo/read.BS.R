##########################
#Bootstrap read
##########################
#Read (and plots) the distribution of bootstraps of multiple sets of trees
#v.0.1
#To do: make the function run for trees in R environment
##########################
#SYNTAX :
#<phy> can be either a phylo or multiPhylo object or a chain name of nexus or newick tree out of R environment.
#<plot> whether to plot the results or not (default = TRUE).
#<save> whether to save the plot or not (default = FALSE).
##########################
#----
#guillert(at)tcd.ie - 28/04/2014
##########################
#Requirements:
#-R 3
#########################


#IN DEVELOPEMENT USE FIRST COMMIT



read.BS<-function(phy,plot=TRUE,save=FALSE){

	#IN DEVELOPEMENT USE FIRST COMMIT
	stop('IN DEVELOPEMENT USE FIRST COMMIT')

	#DATA INPUT

	#phy
	if(class(phy) ==  "phylo"){
		single.tree <- TRUE
		pattern.tree <- FALSE
	} else {
		if(class(phy) == "multiPhylo"){
			single.tree <- FALSE
			pattern.tree <- FALSE
		} else {
			if(class(phy) != "character"){
				stop("phy bad")
			} else {
				pattern.tree <- TRUE
			}
		}
	}
	if(pattern.tree== TRUE) {
    	pattern.list<-list.files(pattern=phy)
	}

	#Read pattern and check if tree is nexus/newick and phylo/multiphylo



    #plot
    if(class(plot) != 'logical'){
        stop('Plot should be TRUE or FALSE')
    }

    #save
    if(class(save) != 'logical'){
        stop('Save should be TRUE or FALSE')
    }

	#FUNCTION

	FUN.read.BS<-function(phy) {
		#phy=pattern
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
	}

	#BOOSTRAP READ

	BS.list<-FUN.read.BS(phy)

	if(plot==TRUE){
		boxplot(BS.list,las=2,ylab="Boostrap")}

	if(save==TRUE){
		pdf("Bootstraps.pdf")
		boxplot(BS.list,las=2,ylab="Boostrap")
		dev.off()
	}

	#Output
	BSlist<-list(Boostraps=unlist(BS.list, use.names=FALSE), Details=BS.list)	
	return(BSlist)

#Draft


	if(class(data) == 'character') {
        data.is.chain<-TRUE
    } else {
        if(class(data) != 'list') {
            stop('The data input is not a list')
        } else {
            data.is.chain<-FALSE
        }
    }

    #If data is a chain name, load the data from the .csv files
    if(data.is.chain == TRUE) {
        chain.list<-list.files(pattern=data)
        data.tmp<-NULL
        
        #Load the csv files
        for (i in 1:length(list.files(pattern=data))) {
            data.tmp[[i]]<-read.csv(list.files(pattern=data)[i], row.names=1)
        }

        #Rename the csv files list
        names(data.tmp)<-unlist(strsplit(chain.list, split=paste('_',data,'.csv',sep='')))
        data<-data.tmp
        data.tmp<-NULL

    }


}
#End












	    if(class(data) == 'character') {
        data.is.chain<-TRUE
    } else {
        if(class(data) != 'list') {
            stop('The data input is not a list')
        } else {
            data.is.chain<-FALSE
        }
    }

    #If data is a chain name, load the data from the .csv files
    if(data.is.chain == TRUE) {
        chain.list<-list.files(pattern=data)
        data.tmp<-NULL
        
        #Load the csv files
        for (i in 1:length(list.files(pattern=data))) {
            data.tmp[[i]]<-read.csv(list.files(pattern=data)[i], row.names=1)
        }

        #Rename the csv files list
        names(data.tmp)<-unlist(strsplit(chain.list, split=paste('_',data,'.csv',sep='')))
        data<-data.tmp
        data.tmp<-NULL

    }