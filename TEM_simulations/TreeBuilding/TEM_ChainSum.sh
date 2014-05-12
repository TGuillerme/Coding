##########################
#Building the tree sets for analysing the TEM simulation results
##########################
#SYNTAX :
#sh TEM_ChainSum.sh <Number of species> <Burnin> <Chain name>
#<Number of species> the number of species per tree
#<Burnin> the proportion of first trees to ignore as a burnin
#<Chain Name> creates a folder named after the chain name to store the results
##########################
#Create .treeset files with all the trees from two mb chains.
#version: 0.4
#Update: Allows a burnin proportion value
#Update: Error in the burnin is now fixed
#Update: Name change and cleaning
#----
#guillert(at)tcd.ie - 17/04/2014
##########################


#Set the variables

nsp=$1
burnin=$2
name=$3
    
header=$nsp
let "header += 5"

FirstTree=$nsp
let "FirstTree += 6"

#Create the folder

mkdir ${name}_treesets

#Create the .treeset file
for f in *.run1.t
    do prefix=$(basename $f .run1.t)

    echo $prefix

    #print the header 
    head -$header ${prefix}.run1.t > ${name}_treesets/${prefix}.treeset

    #Burnin
    BurnTree=$'FirstTree'
    length=$(grep '[&U]' $f | wc -l)
    let "length *= $burnin"
    let "length /= 100"
    let "BurnTree += $length"

    #Add the two list of trees
    sed -n ''"$BurnTree"',$p' ${prefix}.run1.t | sed '$d' >> ${name}_treesets/${prefix}.treeset

    sed -n ''"$BurnTree"',$p' ${prefix}.run2.t >> ${name}_treesets/${prefix}.treeset

done

#end