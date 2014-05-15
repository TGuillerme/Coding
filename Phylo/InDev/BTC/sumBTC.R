##########################
#Summarise BTC function output
##########################
#Calculate the CI and the mode of the difference between two set of trees.
#v.0.2
#Update: the three core functions have been isolated
##########################
#SYNTAX :
#<data> the output list of BTC. Can be a R object (if the list is stored in R) or a chain name (if the output was not stored in R).
#<metric.number> the number of the column in data with the metric of interest.
#<probs> a vector of three probabilities levels (default = c(95,75,50)).
#<plot> whether to plot the results or not (default = TRUE).
#<save.details> whether to save the details of each comparison (default=FALSE). If TRUE, saves a density plot for each comparison. The chain name will be the one given in 'data'.
##########################
#----
#guillert(at)tcd.ie - 19/02/2014
##########################
#Requirements:
#-R 3
#-R package 'hdrcde'
##########################

sumBTC<-function(data, metric.number, probs=c(95, 75, 50), plot=TRUE, save.details=FALSE){

warning('In developement')

#HEADER

#Loading the libraries
    require(hdrcde)

#DATA INPUT

#Data
    #Is the data a list or a chain of characters
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

#Metric
    #Is the metric.number a number?
    if(class(metric.number) == 'numeric'){
        if(length(metric.number) > 1) {
            stop('Only one metric number should be given')
        }
    } else {
        stop('Metric number should be numeric')
    }

#Probabilities
    if(class(probs) == 'numeric'){
        if(length(probs) != 3){
            stop('Probs should contain three probabilities levels')
        }
    } else {
        stop('Probs are not numeric')
    }

#Plot
    #Is plot logical
    if(class(plot) != 'logical'){
        stop('Plot should be TRUE or FALSE')
    }

#save.details
    if(class(save.details) != 'logical'){
        stop('Plot should be TRUE or FALSE')
    }
    if(save.details == TRUE){
        cat('Save.details options is enabled and will generate a pdf file for each comparison', '\n')
        cat('This options might increase the computational time','\n')
    }

#FUNCTIONS

    #Calculates the hdr for each comparison set
    FUN.hdr<-function(data, metric.number, probs) {
        hdr.results<-NULL

        for (i in 1:length(data)) {
            hdr.results[[i]]<-hdr(data[[i]][,metric.number], probs, h = bw.nrd0(data[[i]][,metric.number])) 
        }

        names(hdr.results)<-names(data)

        return(hdr.results)
    }

    #Density Plot function (from densityplot.R by Andrew Jackson - a.jackson@tcd.ie)
    FUN.densityplot <- function (data, metric.number, probs, hdr.results) {

        #Transform dat into a column format
        dat<-data.frame(seq(1:length(data[[1]][,metric.number])))
        for (i in 1:length(data)) {
            dat[,i]<-data[[i]][metric.number]
        }

        #Renaming dat
        names(dat)<-names(data)

        #Set the column numbers
        n<-ncol(dat)

        #Set the y axis as the label name
        ylabels<-colnames(data[[1]][metric.number])

        #Set up the plot
        ylims<-c(min(dat) - 0.1*min(dat), max(dat) + 0.1*(max(dat)))
        xspc<-0.5
        plot(1,1, xlab='', ylab=ylabels, main=paste('','', sep=''), xlim= c(1 - xspc, n + xspc), ylim=ylims, type='n', xaxt='n')
        axis(side = 1, at = 1:n, labels = (as.character(names(dat))), las=2, cex=0.75)

        #Set the colors (grayscale)
        clr = gray((9:1)/10)
        clrs <- rep(clr, 5)

        #Set the scale
        scl=1

        #Plotting the data distribution
        for (j in 1:n) {
            temp <- hdr.results[[j]]
            line_widths <- seq(2, 20, by = 4) * scl
            bwd <- c(0.1, 0.15, 0.2, 0.25, 0.3) * scl

            for (k in 1:length(probs)) {
                temp2 <- temp$hdr[k, ]
                polygon(c(j - bwd[k], j - bwd[k], j + bwd[k], j + bwd[k]), c(min(temp2[!is.na(temp2)]), max(temp2[!is.na(temp2)]), max(temp2[!is.na(temp2)]), min(temp2[!is.na(temp2)])), col = clrs[k])
                points(j,temp$mode,pch=19)
            }
        }
    }  

    #Saves the details for each comparison (graph and values)
    FUN.details.save<-function(data, metric.number, probs) {
        #Compute hdr for one metric
        hdr.saving<-NULL

        for (i in 1:length(data)) {
            pdf(paste(names(data[i]), colnames(data[[1]][metric.number]), 'pdf', sep='.'))
            hdr.saving[[i]]<-hdr.den(data[[i]][,metric.number], prob=probs, xlab=colnames(data[[1]][metric.number]), main=paste(names(data[i]), 'comparison',sep=' '))
            dev.off()
            cat(paste(names(data[i]), colnames(data[[1]][metric.number]), 'pdf', sep='.'), 'density plot saved', '\n')
        }

        #rename the list
        names(hdr.saving)<-paste(names(data),colnames(data[[1]][metric.number]), sep='.')
    }

#SUMMARIZE THE LIST

    #Calculate the hdr for each comparion
    hdr.results<-FUN.hdr(data, metric.number, probs)

    #Optional plot
    if (plot == TRUE) {
        FUN.densityplot(data, metric.number, probs, hdr.results)
    }

    #Optional saving
    if (save.details == TRUE) {
        FUN.details.save(data, metric.number, probs)
    }

    #Output
    return(hdr.results)

#End
}