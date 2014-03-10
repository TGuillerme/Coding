### Name: divComp
### Title: Compute the median of 4 chronograms with median on divergence times + median of the 95% credibility intervals bounds (as written by PhyloBayes).
### Author: Emmanuel J. P. DOUZERY ; October 2011
### What's next? Ultrametricity

#`divComp` <- function	(
#						chronog1 = 61t35g_exon-gb_cat_gtr_gam-ln_bd_rp8318_calib01_min_sb5.dat01.chronogram,
#						chronog2 = 61t35g_exon-gb_cat_gtr_gam-ln_bd_rp8318_calib01_min_sb5.dat02.chronogram,
#						chronog3 = 61t35g_exon-gb_cat_gtr_gam-ln_bd_rp8318_calib01_min_sb5.dat01.chronogram,
#						chronog4 = 61t35g_exon-gb_cat_gtr_gam-ln_bd_rp8318_calib01_min_sb5.dat02.chronogram,
#						) 

library(ape)
library(picante)
library(vegan)
library(nlme)
						
chronog1<-read.tree("61t35g_exonR-gb_cat_gtr_gam-ln_bd_rp8318_calib10_lar01_sb5_mat3_dsom_mtmam.coe01.postmeandates.tre")
chronog2<-read.tree("61t35g_exonR-gb_cat_gtr_gam-ln_bd_rp8318_calib10_lar01_sb5_mat3_dsom_mtmam.coe02.postmeandates.tre")
chronog3<-read.tree("61t35g_exonR-gb_cat_gtr_gam-ln_bd_rp8318_calib10_lar01_sb5_mat3_dsom_mtmam.coe01.postmeandates.tre")
chronog4<-read.tree("61t35g_exonR-gb_cat_gtr_gam-ln_bd_rp8318_calib10_lar01_sb5_mat3_dsom_mtmam.coe02.postmeandates.tre")


{

	require(ape)
	require(picante)


	# Create the median chronogram
	
	chronog_median <- chronog1

	num_nod <- length(chronog_median $ node.label)
	num_tip <- length(chronog_median $ tip.label)
	num_branch <- length(chronog_median $ edge.length)
	
	chronog_median $ node.label <- NULL
	chronog_median $ edge.length <- NULL


	# Extract the time intervals
	
	chronog1 $ edge.length -> delta_time1
	chronog2 $ edge.length -> delta_time2
	chronog3 $ edge.length -> delta_time3
	chronog4 $ edge.length -> delta_time4


	# Extract 95% CredI on node ages

	chronog1 $ node.label -> upp_low_age1
	chronog2 $ node.label -> upp_low_age2
	chronog3 $ node.label -> upp_low_age3
	chronog4 $ node.label -> upp_low_age4

	sapply( strsplit (upp_low_age1, "_"), function (x) as.numeric (x) ) -> upp_low_age1
	sapply( strsplit (upp_low_age2, "_"), function (x) as.numeric (x) ) -> upp_low_age2
	sapply( strsplit (upp_low_age3, "_"), function (x) as.numeric (x) ) -> upp_low_age3
	sapply( strsplit (upp_low_age4, "_"), function (x) as.numeric (x) ) -> upp_low_age4

#	c() -> upp1 ; for (i in 1:num_nod) { upp_low_age1[[i]][1] -> upp1[i] }
#	c() -> upp2 ; for (i in 1:num_nod) { upp_low_age2[[i]][1] -> upp2[i] }
#	c() -> upp3 ; for (i in 1:num_nod) { upp_low_age3[[i]][1] -> upp3[i] }
#	c() -> upp4 ; for (i in 1:num_nod) { upp_low_age4[[i]][1] -> upp4[i] }

#	c() -> low1 ; for (i in 1:num_nod) { upp_low_age1[[i]][2] -> low1[i] }
#	c() -> low2 ; for (i in 1:num_nod) { upp_low_age1[[i]][2] -> low2[i] }
#	c() -> low3 ; for (i in 1:num_nod) { upp_low_age1[[i]][2] -> low3[i] }
#	c() -> low4 ; for (i in 1:num_nod) { upp_low_age1[[i]][2] -> low4[i] }

	c() -> upp1 ; for (i in 1:num_nod) { upp_low_age1[1,i] -> upp1[i] }
	c() -> upp2 ; for (i in 1:num_nod) { upp_low_age2[1,i] -> upp2[i] }
	c() -> upp3 ; for (i in 1:num_nod) { upp_low_age3[1,i] -> upp3[i] }
	c() -> upp4 ; for (i in 1:num_nod) { upp_low_age4[1,i] -> upp4[i] }

	c() -> low1 ; for (i in 1:num_nod) { upp_low_age1[2,i] -> low1[i] }
	c() -> low2 ; for (i in 1:num_nod) { upp_low_age2[2,i] -> low2[i] }
	c() -> low3 ; for (i in 1:num_nod) { upp_low_age3[2,i] -> low3[i] }
	c() -> low4 ; for (i in 1:num_nod) { upp_low_age4[2,i] -> low4[i] }



	# Correct for the root credibility interval (avoid 'NA' at root)
	
	root_age1 <- dist.nodes(chronog1)[1,num_tip+1]
	root_age2 <- dist.nodes(chronog2)[1,num_tip+1]
	root_age3 <- dist.nodes(chronog3)[1,num_tip+1]
	root_age4 <- dist.nodes(chronog4)[1,num_tip+1]

	upp1[1] <- root_age1
	upp2[1] <- root_age2
	upp3[1] <- root_age3
	upp4[1] <- root_age4

	low1[1] <- root_age1
	low2[1] <- root_age2
	low3[1] <- root_age3
	low4[1] <- root_age4


	# Calculate the median of time intervals => to be used for branch lengths of the median chronogram

	cbind(delta_time1,delta_time2,delta_time3,delta_time4) -> time_x4

	c() -> time_median ; for (i in 1:num_branch) { median(time_x4[i]) -> time_median[i] ; print(i) ; print(time_median[i]) }
	print("Median of delta time intervals (= branch lengths): written")


	# Calculate the median of upper and lower bounds of 95% credibility intervals

	cbind(upp1,upp2,upp3,upp4) -> upp_x4
	cbind(low1,low2,low3,low4) -> low_x4

	c() -> upp_median ; for (i in 1:num_nod) { median(upp_x4[i,]) -> upp_median[i] }
	print("Median of upper bounds of 95% credibility intervals: written")

	c() -> low_median ; for (i in 1:num_nod) { median(low_x4[i,]) -> low_median[i] }
	print("Median of lower bounds of 95% credibility intervals: written")

	paste(upp_median, low_median, sep = "_") -> upp_low_median


	# Write the median chronogram

	chronog_median $ edge.length <- time_median
	chronog_median $ node.label <- upp_low_median
	
	write.tree(chronog_median, file = "61t35g_exonR-gb_cat_gtr_gam-ln_bd_rp8318_calib10_lar01_sb5_mat3_dsom_mtmam.postmeandates.median.tre")

	print("Median of the 4 chronograms: written in file 'chronogram.median'")

	
	# Compute patristic distances

	dist.nodes(chronog1) -> pat1
	dist.nodes(chronog2) -> pat2
	dist.nodes(chronog3) -> pat3
	dist.nodes(chronog4) -> pat4


	# Extract node ages

	c() -> node_age1 ; for (i in (num_tip+1) : (num_tip+num_nod) ) { pat1[i, internal2tips(chronog1,i)[1] ] -> node_age1[i] }
	c() -> node_age2 ; for (i in (num_tip+1) : (num_tip+num_nod) ) { pat2[i, internal2tips(chronog2,i)[1] ] -> node_age2[i] }
	c() -> node_age3 ; for (i in (num_tip+1) : (num_tip+num_nod) ) { pat3[i, internal2tips(chronog3,i)[1] ] -> node_age3[i] }
	c() -> node_age4 ; for (i in (num_tip+1) : (num_tip+num_nod) ) { pat4[i, internal2tips(chronog4,i)[1] ] -> node_age4[i] }


	# Restriction to internal nodes

	node_age1[ (num_tip+1) : (num_tip+num_nod) ] -> node_age1
	node_age2[ (num_tip+1) : (num_tip+num_nod) ] -> node_age2
	node_age3[ (num_tip+1) : (num_tip+num_nod) ] -> node_age3
	node_age4[ (num_tip+1) : (num_tip+num_nod) ] -> node_age4


	# Comparison of nodes ages

	cbind(node_age1,node_age2,node_age3,node_age4) -> ages_x4
	
	c() -> max_age_x4 ; for (i in 1:num_nod) { max(ages_x4[i,]) -> max_age_x4[i] }
	c() -> min_age_x4 ; for (i in 1:num_nod) { min(ages_x4[i,]) -> min_age_x4[i] }
	
	diff_age <- max_age_x4 - min_age_x4


	# divComp Summary: among-chronogram divergence times comparison => checking for convergence

	print("Maximum divergence time difference (max_diff)") ;	print(max(diff_age))
	print("Median divergence time difference") ; 	print(median(diff_age))
	print("Minimum divergence time difference") ;	print(min(diff_age))


	# Details on max_diff

	match( max( diff_age ), diff_age ) -> max_diff

	print("max_diff involves the following divergence times:")
	print( ages_x4 [ max_diff, ] )
	
	print("max_diff involves the MRCA of the following taxa:")
	print( internal2tips( chronog_median, max_diff + num_tip, return.names = T ) )

	("REMINDER: Check the ultrametricity of 'chronogram.median'")

}
