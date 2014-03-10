##########################
#Measurement error checking
##########################
#Checking the measurement error following Cooper & Purvis 2009 method
#v1.1
#----
#sfinlay(at)tcd.ie & guillert(at)tcd.ie - 15/04/2013
##########################


MesErr<-function(data, ID=1, measurement=2, value=3, coeff.var=5, spread=25, filename="filename")
{
#Input checking
	if (class(data) != "data.frame") {
	stop("'data' is not a data.frame object")
	}
	if (class(ID) != "numeric") {
	stop("ID argument is not numeric")
	}
	if (class(measurement) != "numeric") {
	stop("measurement argument is not numeric")
	}
	if (class(value) != "numeric") {
	stop("value argument is not numeric")
	}
	if (class(coeff.var) != "numeric") {
	stop("maximum variation coefficient argument is not numeric")
	}
	if (class(value) != "numeric") {
	stop("maximum spread percentage argument is not numeric")
	}
	if (class(data[[value]]) != "numeric") {
	stop("Values are not numeric")
	}
	if (class(filename) != "character") {
	stop("Output name is not character")
	}

#Selecting all the replicates of one measurement on one specimen
for (i in 1:length(levels(data[[ID]]))){
	for (j in 1:length(levels(data[[measurement]]))){

		one.measure<-data[[value]][data[[ID]] == levels(data[[ID]])[2] & data[[measurement]] == levels(data[[measurement]])[1]]

		#median
		med<-median(one.measure, na.rm=T)

		#coef.var
		coef.var<-(sd(one.measure,na.rm=T)/mean(one.measure,na.rm=T))*100

	##Spread
		###Removing NA
		if(!is.na(one.measure[1])){
		
		###3 replicats
		if(length(one.measure)==3){
		
			####Spread
			if((median(one.measure, na.rm=T)-min(one.measure, na.rm=T)) <= (max(one.measure, na.rm=T)-median(one.measure, na.rm=T))){
				per.spread<-((median(one.measure, na.rm=T)-min(one.measure, na.rm=T))/(max(one.measure, na.rm=T)-min(one.measure, na.rm=T)))*100}
    			else{
    			per.spread<-((max(one.measure, na.rm=T)-median(one.measure, na.rm=T))/(max(one.measure, na.rm=T)-min(one.measure, na.rm=T)))*100}
		
		} else {
		###5 replicats
			####Spread
			sorted<-sort(one.measure)		
			diff<-rep(NA,(length(levels(specimen[[measurement]]))-1))
			for(k in 1:(length(levels(specimen[[measurement]])))-1){
			diff[k]<-sorted[k+1]-sorted[k]}
			sortdiff<-sort(diff)
			e<-(max(one.measure)-min(one.measure))
			per.spread<-(sortdiff[1]/e + sortdiff[2]/e+ sortdiff[2]/e)*100
			}
		}

	output<-data.frame(ID=levels(data[[ID]])[i] ,measurement=levels(specimen[[measurement]])[j],median=Median, spread=per.spread, var=coef.var)
	write.table(file=filename,output,col.names=F, row.names=F,sep="\t",quote=F,append=T)	
	
}
}	
	
	##Error checking
		###data frame
		output<-read.table(filename, col.names=c("ID", "measurement", "median", "spread", "var"))
	
		###Error checking
		Error<-rep(FALSE, nrow(output))
		for (i in 1:nrow(output)){
			if(is.na(output$var[i])){
				Error[i]<-NA
				} else {
					if(output$var[i]>=coeff.var){
						if(output$spread[i]>=spread){
							Error[i]<-TRUE
							} else {
								Error[i]<-FALSE}
						} else {
							Error[i]<-FALSE}
				}
		}
		errors<-which(Error,TRUE)


#Output
	ret<-list(errors=output[c(errors),c(1,2)], details=output)
	cat(length(errors),"measurement errors encountered","\n")
	print(output[c(errors),c(1,2)])	
	return(ret)
	}