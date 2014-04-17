##########################
#TreeCmp-cluster.tasker
##########################
#SYNTAX: sh TreeCmp-cluster.tasker <from> <to> <treestype>
#with:
#<from> a nexus file with the reference trees
#<to> a nexus file with the trees to compare to the reference tree
#<treestype> number of species per trees
##########################
#version: 0.1
#----
#guillert(at)tcd.ie - 17/04/2014
##########################
#Requirements:
#-TEM_TreeCmp.sh script
#-R 3.x
#-TreeCmp java script
#-http://www.la-press.com/treecmp-comparison-of-trees-in-polynomial-time-article-a3300
#-TreeCmp folder to be installed at the level of the analysis
##########################

#!/bin/sh
#INPUT

#Input values
from=$1
to=$2
treestype=$3


#CREATING THE TASK FILES
for n in $(seq $from  $to)
do
    #remove the '\' in first line
    echo "#@/bin/sh
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH -p compute
#SBATCH -J TCmp${n}
source /etc/profile.d/modules.sh
export http_proxy=http://proxy.tchpc.tcd.ie:8080
module load cports6 openmpi/1.6.5-gnu4.8.2 
module load cports gsl/1.16-gnu

##########################
#TASK FILE - Chain${n}
##########################

sh TreeCmp_Chain${n}.task ;" > TreeCmp_Chain${n}.tmp

sed 's/@/!/g' TreeCmp_Chain${n}.tmp >  TreeCmp_Chain${n}.job

rm  TreeCmp_Chain${n}.tmp

#done

#Modified for pruned trees. Use prefix folder '_treesets_living/fossil' or '_treesets'.
#Modified for pruned trees. Use prefix file '..treesets' or '.treesets'
#Modified for pruned trees. Use 'pruned-M${prefix_folder}_L00F00C00..treeset' or 'M${prefix_folder}_L00F00C00.treeset'

#for n in $(seq 20  29)
#do

    echo "for folder in Chain${n}_treesets_living
do
    prefix_folder=\$(basename \$folder _treesets_living)
    echo \$prefix_folder

    #Get the right folder
    for file in \${folder}/*.living
    do
        #In the right folder get the right trees
        prefix_file=\$(basename \$file .living)
        sh TEM_TreeCmp.sh \${prefix_folder}_treesets_living/pruned-M\${prefix_folder}_L00F00C00.living \${file} 51 1000 \${prefix_file} #or \${prefix_folder}.True_tree.tre.nex M\${prefix_folder}_L00F00C00.treeset
    done
done " > TreeCmp_Chain${n}.task

done