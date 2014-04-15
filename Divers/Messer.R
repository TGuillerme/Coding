##########################
#Measurement error checking
##########################
#Checking the measurement error following Cooper & Purvis 2009 method (spread) or Yezerinac et al 1992 (variance)
#v.2
#Update: Cooper & Purvis and Yezerinac method allowed
##########################
#SYNTAX :
#<data> a data.frame object containing the data
#<variable> the name of the variable to analyse. Must match with at least one element in data.
#<format> the format of the data: whether the measurements are displayed in a row or in a column: can be either 'row' or 'column'
#	the row format must be:
#	ID	Variables	Measures
#	a 	b	x
#
#	the column format must be:
#	ID	Variable 	Variable
#	a 	x	y
#
#<error> the amount of error tolerated for not rejecting a measure (default=5). For the spread method, the spread tolerance is calculated as 25+error.
#	If the error is set to >25, then the spread tolerance is not correct following Cooper & Purvis 2009.
#<method> which method to use: can be 'spread' (Cooper & Purvis 2009) or 'variance' (Yezerinac et al 1992) or 'both' (default='both')
##########################
#----
#guillert(at)tcd.ie - 14/04/2014
##########################
#Requirements:
#-R 3

Messer<-function(data, format, variable, error=5, method='both')
{

#INPUT

	#data
	if (class(data) != "data.frame") {
		stop("'data' is not a data.frame object")
	}

	#format
	if(class(format) != "character") {
		stop("Format must be 'row' or 'column'")
	} else {
		if(format != 'row') {
			if(format != 'column') {
				stop("Format must be 'row' or 'column'")
			}
		}
	}

	#variable
	if (format == "row") {
		var.test<-grep(variable, colnames(data))
		if(length(var.test) == 0) {
			stop("Variable not found")
		} else {
			var=var.test
		}
	} else {
		var.test<-grep(variable, data[,2])
		if(length(var.test) == 0) {
			stop("Variable not found")
		} else {
			var=var.test
		}
	}

	#error
	if (class(error) != "numeric") {
		stop("Tolerated error is not numeric")
		if(error > 25) {
			warning("Tolerated error is >25: spread tolerance is now >50 and no more following Cooper and Purvis (2009) recommendations")
		}
	}

	#method
	if(class(method) != "character") {
		stop("Method must be 'spread', 'variance' or 'both'")
	} else {
		if(method != 'both') {
			if(method != 'spread') {
				if(method != 'variance') {
					stop("Method must be 'spread', 'variance' or 'both'")
				}
			}			
		}
	}


#FUNCTIONS

	#transform the table in column/row format for spread/variance method
	FUN.transform.table<-function(data, format, method, variable, var) {

		#Column format required for method spread
		if(method=='spread') {
			if(format=='column') {
				transformed.data<-data[c(var),]
			} else {

				#transform to column format
				transformed.data<-data.frame('ID'=data[,1], 'Variable'=rep(variable, nrow(data)), 'Measurements'=data[,var])
			}
		
		#method == 'variance'
		} else {
			if(format=='row') {
				transformed.data<-data[,c(1, var)]
			} else {

				#transform to row format
				transformed.data<-data.frame('ID'=data[c(var),1], variable=data[c(var),3])
				names(transformed.data)[2]<-variable
			}
		}

		#output
		return(transformed.data)
	}

	#Spread measurement error method (Cooper & Purvis 2009)
	FUN.method.spread<-function(transformed.data) {

		#Input
		specimens<-levels(transformed.data[,1])
		n.specimens<-length(specimens)

		#Measure the median for each specimen
		spec.median<-NULL
		for(n in 1:n.specimens){
			spec.median[n]<-median(transformed.data[c(grep(specimens[n], transformed.data[,1])),3], na.rm=TRUE)
		}

		#Measure the coefficient of variation for each specimen
		spec.coeff<-NULL
		for(n in 1:n.specimens){
			spec.coeff[n]<-( sd(transformed.data[c(grep(specimens[n], transformed.data[,1])),3], na.rm=TRUE) / mean(transformed.data[c(grep(specimens[n], transformed.data[,1])),3], na.rm=TRUE) )*100
		}

		#Measure the spread for each specimen
		spec.spread<-NULL
		for(n in 1:n.specimens){
			measurements<-transformed.data[c(grep(specimens[n], transformed.data[,1])),3]

			#remove eventual NAs
			if(any(is.na(measurements))) {
				measurements<-measurements[-which(is.na(measurements))]
			}

			#sorted list of observations (i(1)<i(2)<...<i(n))
			measurements<-sort(measurements)
			i.measurements<-length(measurements)

			#Calculate the pairwise differences as diff(1)=|i(1)-i(2)|, diff(2)=|i(2)-i(3)|, diff(n-1)=|i(n-1)-i(n)|
			meas.diff<-NULL
			for (i in 1:(i.measurements-1)){
				meas.diff[i]<-abs(measurements[1+(i-1)]-measurements[1+i])
			}

			#calculate the spread as ((diff(1) + diff(2) + .. + diff(n-2)) / diff(n-1)) * 100
			spec.spread[n]<-( sum(meas.diff[-length(meas.diff)]) / abs(measurements[1]-measurements[i.measurements]))*100
		}

		#Return the results
		results<-data.frame('ID'=specimens, 'variable'=rep(variable, n.specimens), 'median'=spec.median, 'coeff.var'=spec.coeff, 'spread'=spec.spread, 'error'=rep(NA, n.specimens))
		return(results)
	}

	#Variance measurement error method (Yezerinac et al 1992)

	FUN.method.variance<-function(transformed.data){
		
		#Input
		specimens<-levels(transformed.data[,1])
		n.specimens<-length(specimens)
		observations<-length(transformed.data[,1])

		#measurements per specimens
		n.measurements<-NULL
		for (n in 1:n.specimens) {
			n.measurements[n]<-length(grep(specimens[n], transformed.data[,1]))
		}

		if(length(levels(as.factor(n.measurements))) != 1){
			warning('The measurements per specimens varies. \nThe Sum squared among is approximated by using:\nmean of (mean squared within - mean squared among) / number of measures per specimen.')
		}

		#aov model
		model<-summary(aov(transformed.data[,2]~transformed.data[,1]))

		#Mean squared
		MS.within<-model[[1]][2,3]
		MS.among<-model[[1]][1,3]

		#Sum squared
		SS.within<-model[[1]][2,2]
		SS.among<-mean((MS.among-MS.within)/n.measurements) #Sokal and Rohlf 1981

		#Calculating the variable error
		var.error<-SS.within/(SS.within + SS.among)*100 #Yezerinac et al 1992

		#return the result
		return(var.error)
	}

#MEASUREMENT ERROR CHECKING

	#Calculating the spread only
	if (method == 'spread') {
		results<-FUN.method.spread(FUN.transform.table(data, format, 'spread', variable, var))

		#Error checking
		for(n in 1:nrow(results)){
			if(results[n,4]>=error){
				if(results[n,5]>=(25+error)){ #Acceptance spread is set a 25+error as a rule of thumb (acceptance should be < 50)
					results[n,6]<-TRUE
				} else {
					results[n,6]<-FALSE
				}
			} else {
				results[n,6]<-FALSE
			}
		}

		return(results)
	}

	#Calculating the variance only
	if (method == 'variance') {
		results<-FUN.method.variance(FUN.transform.table(data, format, 'variance', variable, var))
		cat('Measurement error of the variable is:', results, '%', '\n')

		return(results)
	}

	#Calculating both

	if (method == 'both') {

		#Calculating the spread
		results<-FUN.method.spread(FUN.transform.table(data, format, 'spread', variable, var))

		#Error checking
		for(n in 1:nrow(results)){
			if(results[n,4]>=error){
				if(results[n,5]>=(25+error)){ #Acceptance spread is set a 25+error as a rule of thumb (acceptance should be < 50)
					results[n,6]<-TRUE
				} else {
					results[n,6]<-FALSE
				}
			} else {
				results[n,6]<-FALSE
			}
		}

		#Calculating the variance
		results.variance<-FUN.method.variance(FUN.transform.table(data, format, 'variance', variable, var))

		#returning both
		result.list<-list(variable.error=results.variance, specimen.error=results)
		return(result.list)
	}

#End
}