##########################
#Calculate the mode of a distribution
##########################
#v.1
##########################
#SYNTAX :
#<x> a vector of numerical values.
##########################
#----
#guillert(at)tcd.ie - 26/05/2014
##########################
#Requirements:
#-R 3
##########################

Mode<-function(x){
        as.numeric(names(sort(-table(x))[1]))
        }