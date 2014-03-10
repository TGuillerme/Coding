##########################
#Uses TreeCmp java script to compare the distributions of Bayesian trees
##########################
#SYNTAX: TEM_TreeCmp_v <reference.treeset> <input.treeset> <number of species> <number of random draws> <output>
#with:
#<reference.treeset> a nexus file with the reference trees
#<input.treeset> a nexus file with the trees to compare to the reference tree
#<number of species> number of species per trees
#<number of random draws> number of random draws for the comparison
#<output> a chain name for the output
##########################
#version: 0.2
#Update: improved output format
#----
#guillert(at)tcd.ie - 26/02/2014
##########################
#Requirements:
#-R 3.x
#-TreeCmp java script
#-http://www.la-press.com/treecmp-comparison-of-trees-in-polynomial-time-article-a3300
#-TreeCmp folder to be installed at the level of the analysis
##########################

#INPUT

#Input values
REFtreeset=$1
INPtreeset=$2
species=$3
draws=$4
output=$5

    #testers
    #REFtreeset=$(echo 'Tree1.tre')
    #INPtreeset=$(echo 'Tree2.tre')
    #species=$(echo '51')
    #draws=$(echo '10')
    #output=$(echo 'test_output')

#RANDOM DRAWS LIST

#Make the nexus file header
header=$species
let "header += 5"
head -$header $REFtreeset > HEADER_${output}.Tree.tmp

#Make the random draws in both list of trees
REFntrees=$(grep '\[\&U\]' $REFtreeset | wc -l)
INPntrees=$(grep '\[\&U\]' $INPtreeset | wc -l)

#Create the list of trees to sample. If the provided number of random draws is higher than the number of trees, the sample is done with replacement.
echo "if($REFntrees < $draws) {
        REFrep=TRUE } else {
        REFrep=FALSE}
    if($INPntrees < $draws) {
        INPrep=TRUE } else {
        INPrep=FALSE }    
    write(sample(seq(1:$REFntrees), $draws, replace=REFrep), file='REFtreeset.sample', ncolumns=1)
    write(sample(seq(1:$INPntrees), $draws, replace=INPrep), file='INPtreeset.sample', ncolumns=1) " | R --no-save

#TREE COMPARISONS

#Creates the files of single trees and compare them one to one
for n in $(seq 1 $draws)
do
    #Creates the ref tree file for one draw
    cp HEADER_${output}.Tree.tmp REFtreeset_tree${n}.Tree.tmp
    REFrand=$(sed -n ''"${n}"'p' REFtreeset.sample)
    let "REFrand += $header"
    sed -n ''"$REFrand"'p' $REFtreeset >> REFtreeset_tree${n}.Tree.tmp
    echo 'end;' >> REFtreeset_tree${n}.Tree.tmp

    #Creates the input tree file for one draw
    cp HEADER_${output}.Tree.tmp INPtreeset_tree${n}.Tree.tmp
    INPrand=$(sed -n ''"${n}"'p' INPtreeset.sample)
    let "INPrand += $header"
    sed -n ''"$INPrand"'p' $INPtreeset >> INPtreeset_tree${n}.Tree.tmp
    echo 'end;' >> INPtreeset_tree${n}.Tree.tmp

    #Make the comparison using the TreeCmp java script on all rooted metrics
    java -jar TreeCmp/bin/TreeCmp.jar -r REFtreeset_tree${n}.Tree.tmp -d mc rc ns tt -i  INPtreeset_tree${n}.Tree.tmp -o ${output}_draw${n}.Cmp.tmp
done

#Removing the extra trees
rm *.Tree.tmp

#SUMMARIZING THE COMPARISONS

#Comparison header
sed -n '1p' ${output}_draw1.Cmp.tmp | sed $'s/Tree/Ref.trees\t\Input.trees/g' > ${output}.Cmp

#Add the values from each comparison
for n in $(seq 1 $draws)
do
    REFrand=$(sed -n ''"${n}"'p' REFtreeset.sample)
    INPrand=$(sed -n ''"${n}"'p' INPtreeset.sample)
    sed -n '2p' ${output}_draw${n}.Cmp.tmp | sed 's/^./'"$REFrand"'@'"$INPrand"'/' | sed $'s/@/\t/' >> ${output}.Cmp
done

#Removing extra comparisons
rm *.Cmp.tmp

#Removing sample files
rm *.sample
#end