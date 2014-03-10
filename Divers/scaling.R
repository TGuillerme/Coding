##########################
#Transforming points in to length vectors
##########################
#Transforms vector coordinates into vector length and scale it with provided scale bar vectors coordinates
#----
#guillerme.thomas(at)gmail.com - 8/06/2012
##########################


scaling<-function(x1,x2,x3)
	{
	#Input synthax
	if(missing(x1)) {stop("No input file name provided")} else
		{if(missing(x2)) {stop("No output file name provided")} else
			{{if(missing(x3)) {stop("No scale provided")} else {bli<-"bli"}}
			data<-read.table(x1, header=T, dec=".")
			#Input checking
			if((-is(data,"data.frame"))==0) {stop("Input file provided is not a data frame")} else
			{if(length(data)==9) {bla<-"bla"} else {stop("Input file need to contain only 9 columns : X1=name, X2-3=Xcoordinates, X4-5=Ycoordinates, X6-7=Xscale, X8-9=Yscale")}}
			{if(is(x3, "numeric")) {blu<-"blu"} else {stop("Provided scale is not numeric")}}
			#Scaling
			n<-nrow(data)
			Leng<-rep(NA,n)
			Length<-rep(NA,n)
			Scale<-rep(NA,n)
			for (i in 1:n){
				Leng[i]<-sqrt((data[i,4]-data[i,2])^2+(data[i,5]-data[i,3])^2)
				Scale[i]<-sqrt((data[i,8]-data[i,6])^2+(data[i,9]-data[i,7])^2)
				Length[i]<-Leng[i]*x3/Scale[i]
				}
			Z<-data.frame(names=data[,1], length=Length)
			write.table(Z, file=x2)
			cat("New data written in file", "\n",x2, "\n")
			return(Z)
			}
		}
	}