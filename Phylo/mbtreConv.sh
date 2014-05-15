##########################
#Concatenate and convert MrBayes consensus trees (*.con.tre) into a unique <Chain>_trees.tre
##########################
#SYNTAX :
#sh mbtreConv.sh <Chain>
#with
#<Chain> optional chain name
##########################
#Concatenate and convert MrBayes consensus trees (*.con.tre) into a unique <CHAIN>_trees.tre file with species names, branch length and posterior probabilities as node labels
#version: 1.0
#----
#guillert(at)tcd.ie - 14/11/2013
##########################



#Input
Chain=${1:-Sum}



#Creating the temporary folder 
mkdir .tmp_mbtrConv
cp *.con.tre .tmp_mbtrConv/
cd .tmp_mbtrConv/
ls *.con.tre > ls.tmp
Tree=$(sed -n '1p' ls.tmp)



#Put all the trees in a same file, keep only the posterior probability node labels and rename them



##Selecting the number of species
ntax=$(grep 'dimensions ntax=' ${Tree} | sed 's/dimensions ntax=//g' | sed 's/;//g' | sed 's/ //g')


##Selecting the header part of the nexus file (taxa list + translate)
nexushead=$ntax
let "nexushead += $ntax"
let "nexushead +=11"
sed -n '1,'"$nexushead"'p' ${Tree} > trees.tmp



##Copy all the trees with the nexus header
for f in *.con.tre ;
do echo $f ;
echo $f > tree.tre ;
grep "tree con" $f >> tree.tre ;



##Removing the annotations and keeping only the posterior probabilities as node labels
sed $'s/\[\&prob=/\\\n\[\&prob=/g' tree.tre |sed $'s/]:/]\\\n:/g' | sed 's/,prob_stddev=.*]/]/g' | sed 's/\[&length_mean=.*]//g' > tree.tmp ;
tr '\n' '@' < tree.tmp | sed 's/@//g' | sed 's/\([0-9]\)\[&prob=.\....................\]/\1/g' | sed 's/\[&prob=//g' | sed 's/]:/:/g' >> trees.tmp ;
rm tree.tre ;
rm tree.tmp ;
done ;



##Renaming the trees
echo "end;" >> trees.tmp

sed 's/.con.tre   tree con_50_majrule//g' trees.tmp | sed 's/M[0-9]*_/tree /g' | sed 's/\[ID: [0-9]*\]/[Nexus file generated using mbtreConv.sh: '"$(date)"']/g' > ${Chain}_trees.tre
rm trees.tmp



#Saving
cp ${Chain}_trees.tre ../
cd ..
rm -R .tmp_mbtrConv



#End