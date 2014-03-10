##########################
#REMOVING TIPS
##########################
#Remove tips from a tree to match with a data set
#----
#guillerme.thomas(at)gmail.com - 8/11/2012
##########################


rmTip<-function(x1, x2, x3, x4)
	{
	#Library checking
	require(ape)
	
	#Input synthax
	if(missing(x1)) {stop("rmTip(Arg1, Arg2, Arg3, Arg4) ; Arg1=phylo object ; Arg2=a data frame with names of species in first column ; Arg3=the first column of interest name ; Arg4=the second column of interest name")} else
	
	#Input checking
		{if((-is(x1, "phylo"))==0) {stop("Arg1 is not a 'phylo' object")} else
			{if((-is(x2, "data.frame"))==0) {stop("Arg2 is not a 'data.frame' object")} else

				#Tree checking
				{if((-is.rooted(x1))==0)  {stop("The tree must be rooted use function : root(<arg1>, <outgroup>, resolve.root=TRUE)")} else
					{if((-is.binary.tree(x1))==0) {warning("Tree is now set to dichotomous") ; x1<-multi2di(x1)}}
						
						{
						#Selecting NA lines in Arg2$Arg3
						NAx3<-which(is.na(x2[,x3]))
						NAx4<-which(is.na(x2[,x4]))
						X2<-x2[-c(NAx3,NAx4),]				
					
						#Removing tips from the tree (taxa with <NA>)
						tiplab<-x1$tip.label
						datnames<-X2[,1]
						a<-datnames[match(tiplab, datnames)]
						droptip<-which(is.na(a))
						X1<-drop.tip(x1,droptip)
						
						#Removing species that are not in the tree
						tipl<-X1$tip.label
						b<-tipl[match(datnames,tipl)]
						bb<-which(is.na(b))
						Xx2<-X2[-c(bb),]
						
						#Tips checking
						nsp1<-("Number of species ?")
						FinalSpNumber<-(nrow(Xx2))
						finc1<-("Same number of tips and of species ?")
						finc2<-(nrow(Xx2) == length(tipl))
						
						#Output data
						outputData<<-Xx2
						outputTree<<-X1
							cat("New data written in", "\n","outputData", "\n")
							cat("New tree written in", "\n", "outputTree", "\n")
							cat("Same number of species in the new data and the new tree ?", "\n", finc2, "\n")
							cat("Final number of species", "\n", FinalSpNumber, "\n")
							}

				 }
							
							
				}				
			}
		 }