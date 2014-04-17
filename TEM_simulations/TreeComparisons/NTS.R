##########################
#Calculate the Normalised Tree Similarity
##########################
#Calculate the NTS as defined by Bogdanowicz et al. 2012 : For any tree with n taxa compared using a tree distance metric m, NTSm represents the similarity score between the two trees given the expected distance between two random Yule trees with n taxa..
#v.0.2
#Update: now supports TreeCmp.Read entries
#Update: scaling option has been removed
##########################
#SYNTAX :
#<diff> A value of difference/similarity between two trees for any metric (multiple metrics can be providen)
#<rand> A value of the expected difference/similarity between two random trees for any metric (has to be the same metric(s) used for diff argument). If more than one value is provided for the expected difference/similarity, a mean is calculated per metric.
##########################
#----
#guillert(at)tcd.ie - 18/03/2014
##########################
#Requirements:
#-R 3
##########################

NTS<-function(diff, rand) {

#DATA INPUT

    #Is diff a TreeCmp object
    if(class(diff) == 'TreeCmp') {
        diff.length<-length(diff)
        metric.length<-length(diff[[1]])-2
        metric.names<-colnames(diff[[1]])[3:(metric.length+2)]
    } else {
    #Is diff numeric
        if(class(diff) != 'numeric') {
            stop('Diff has to be a numeric or TreeCmp class object')
        
        } else {
            metric.length<-1
            if(length(diff) >= 2) {
                diff.vector<-TRUE
            } else {
                diff.vector<-FALSE
            }
        }
    }

    #Is rand a data.frame
    if(class(rand) == 'data.frame') {
        rand.length<-length(rand)-2
        rand.names<-colnames(rand)[3:(rand.length+2)]
        
        if (all(rand.names == metric.names)) {

            #Calculating the means for the randoms
            rands<-NULL
            for (i in 1:(rand.length)) {
                rands[i]<-mean(rand[,(2+i)])
            }

        } else {
            stop('Metrics in rand and diff do not match')
        }

    } else {

        #Is rand numeric
        if(class(rand) != 'numeric') {
            stop('The mean value of the random difference/similarity is not numeric')
        } else {
            if(length(rand) != 2) {
                if(metric.length == 1) {
                    rand<-mean(rand)
                    warning('Random difference/similarity given is now the mean of the given vector')
                } else {
                    if(length(rand) == metric.length) {
                        rands<-rand
                    } else {
                        stop('Provided number of random values does not match the number of metrics')
                    }
                }
            }
        }
    }



#FUNCTION

#Normalised Tree Similarity
    FUN.NTS<-function(diff, rand) {
        nts<-( (rand-diff) / rand )
    }

#CALCULATE THE NTS
    #NTS for a single metric
    if(metric.length == 1) {
        nts<-FUN.NTS(diff,rand)
        return(nts)
    
    #NTS for a TreeCmp object
    } else {
        nts<-diff
        for (i in 1:diff.length) {
            for (j in 1:metric.length){
                nts[[i]][,j+2]<-FUN.NTS(diff[[i]][,j+2],rands[j])
            }
        }
    }
    #Output
    
    return(nts)
#end
}