##########################
#Calculate the Normalised Tree Similarity
##########################
#Calculate the NTS as defined by Bogdanowicz et al. 2012 : For any tree with n taxa compared using a tree distance metric m, NTSm represents the similarity score between the two trees given the expected distance between two random Yule trees with n taxa..
#v.0.1
##########################
#SYNTAX :
#<diff> A value of difference/similarity between two trees for the metric m (e.g. RF)
#<rand> The mean value of the difference/similarity between two random trees for the same metric m as used in 'diff'
#<scale> whether or not to scale the NTS (default=FALSE), if the scaling option is used, scale must be equal to n (the number of tips of the compared phylogenies)
##########################
#----
#guillert(at)tcd.ie - 18/03/2014
##########################
#Requirements:
#-R 3
##########################

NTS<-function(diff, rand, scale=FALSE){

#HEADER

#Loading the libraries

#DATA INPUT

#Data
    #Is diff numeric
    if(class(diff) != 'numeric') {
        stop('The difference/similarity value is/are not numeric')
    } else {
        if(length(diff) >= 2) {
            diff.vector<-TRUE
        } else {
            diff.vector<-FALSE
        }
    }

    #Is rand numeric
    if(class(rand) != 'numeric') {
        stop('The mean value of the random difference/similarity is not numeric')
    } else {
        if(length(rand) >= 2) {
            rand<-mean(rand)
            warning('Random difference/similarity given is now the mean of the given vector')
        }
    }

    #Scale
    if(scale == FALSE) {
        scale<-1
    } else {
        if(class(scale) != 'numeric') {
           stop('Scale should be either FALSE or a numerical value')
        }
    }

#FUNCTION

#NTS
    FUN.NTS<-function(diff, rand, scale) {
        nts<-( (rand-diff) / rand ) / scale
        return(nts)
    }

#CALCULATE THE NTS
    nts<-FUN.NTS(diff,rand,scale)

    #Output
    return(nts)
#end
}