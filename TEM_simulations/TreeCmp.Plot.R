##########################
#Plot TreeCmp output
##########################
#Calculate the CI and the mode of a given metric using TreeCmp class objects
#Function modified from sumBTC (guillert(at)tcd.ie - 19/02/2014)
#v.0.1
#Update: the three core functions have been isolated
##########################
#SYNTAX :
#<data> an object of the class TreeCmp
#<metric.number> the name of the metric of interest
#<probs> a vector probabilities levels (default = c(95,75,50)).
#<plot> whether to plot the results or not (default = TRUE).
#<save.details> whether to save the details of each comparison (default=FALSE). If TRUE, saves a density plot for each comparison. The chain name will be the one given in 'data'.
##########################
#----
#guillert(at)tcd.ie - 21/03/2014
##########################
#Requirements:
#-R 3
#-R package 'hdrcde'
#-R TreeCmp objects
##########################

TreeCmp.Plot<-function(data, metric, probs=c(95, 75, 50), plot=TRUE, save.details=FALSE) {

#HEADER

#Loading the libraries
    require(hdrcde)

#DATA INPUT

    #data
    if(class(data) != 'TreeCmp') {
        stop('The data must be a TreeCmp object')
    } else {
        data.length<-length(data)
        if(data.length == 0) {
            stop('The provided TreeCmp object is empty')
        } else {
            data.columns<-length(data[[1]])
            if(data.columns == 0) {
                stop('The provided TreeCmp object is empty')
            }
        }
    }

    #metric
    if(class(metric) != 'character') {
        stop('Provided metric name is not found in the data')
    } else {
        if(length(metric) != 1) {
            stop('Only one metric name can be provided')
        } else {
            metric.column<-grep(metric, colnames(data[[1]]))
            if(length(metric.column) == 0) {
                stop('Provided metric name is not found in the data')
            } else {
                data.rows<-NULL
                for (i in 1:data.length) {
                    data.rows[i]<-length(data[[i]][,metric.column])
                }
            }
        }
    }

    #probs
    if(class(probs) == 'numeric'){
        if(length(probs) != 3){
            stop('Probs should contain three probabilities levels')
        }
    } else {
        stop('Probs are not numeric')
    }

    #plot
    if(class(plot) != 'logical'){
        stop('Plot is not logical')
    }

    #save.details
    if(class(save.details) != 'logical'){
        stop('Save.details is not logical')
    }
    if(save.details == TRUE){
        cat('Save.details options is enabled and will generate a pdf file for each comparison', '\n')
        cat('This options might increase the computational time','\n')
    }

#FUNCTIONS

    #Calculates the hdr for each comparison set
    FUN.hdr<-function(data, metric.column, probs) {
        hdr.results<-NULL

        for (i in 1:length(data)) {
            hdr.results[[i]]<-hdr(data[[i]][,metric.column], probs, h = bw.nrd0(data[[i]][,metric.column])) 
        }

        names(hdr.results)<-names(data)

        return(hdr.results)
    }

    #Density Plot function (from densityplot.R by Andrew Jackson - a.jackson@tcd.ie)
    FUN.densityplot<-function (data, metric.column, data.rows, probs, hdr.results) {

        #Transform dat into a column format
        dat<-matrix(NA, nrow=max(data.rows), ncol=data.length)
        dat<-as.data.frame(dat)
        for (i in 1:length(data)) {
            #The number of rows in the data frame is equal to the maximum number of rows in the list. If the elements in the list don't have the same number of rows, NAs are added.
            dat[,i]<-c(data[[i]][,metric.column], rep(NA,(max(data.rows)-length(data[[i]][,metric.column]))))
        }

        #Renaming dat
        names(dat)<-names(data)

        #Set the column numbers
        n<-ncol(dat)

        #Set the y axis as the label name
        ylabels<-colnames(data[[1]][metric.column])

        #Set up the plot
        ylims<-c(min(dat, na.rm=TRUE) - 0.1*min(dat, na.rm=TRUE), max(dat, na.rm=TRUE) + 0.1*(max(dat, na.rm=TRUE)))
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
    FUN.details.save<-function(data, metric.column, probs) {
        #Compute hdr for one metric
        hdr.saving<-NULL

        for (i in 1:length(data)) {
            pdf(paste(names(data[i]), colnames(data[[1]][metric.column]), 'pdf', sep='.'))
            hdr.saving[[i]]<-hdr.den(data[[i]][,metric.column], prob=probs, xlab=colnames(data[[1]][metric.column]), main=paste(names(data[i]), 'comparison',sep=' '))
            dev.off()
            cat(paste(names(data[i]), colnames(data[[1]][metric.column]), 'pdf', sep='.'), 'density plot saved', '\n')
        }

        #rename the list
        names(hdr.saving)<-paste(names(data),colnames(data[[1]][metric.column]), sep='.')
    }

#SUMMARIZE THE LIST

    #Calculate the hdr for each comparion
    hdr.results<-FUN.hdr(data, metric.column, probs)

    #Optional plot
    if (plot == TRUE) {
        FUN.densityplot(data, metric.column, data.rows, probs, hdr.results)
    }

    #Optional saving
    if (save.details == TRUE) {
        FUN.details.save(data, metric.column, probs)
    }

    #Output
    return(hdr.results)

#End
}