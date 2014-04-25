##########################
#TreeCmp-cluster.tasker
##########################
#SYNTAX: sh TreeCmp-cluster.tasker <from> <to> <treestype>
#with:
#<from> start chain
#<to> end chain
#<treestype> treesets type (fossil or living)
##########################
#version: 0.2
#----
#guillert(at)tcd.ie - 22/04/2014
#Update: treestype allows to specify fossil or living treesets.
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

#Modified for pruned trees. Use prefix folder '_treesets_${treestype}/fossil' or '_treesets'.
#Modified for pruned trees. Use prefix file '..treesets' or '.treesets'
#Modified for pruned trees. Use 'pruned-M${prefix_folder}_L00F00C00..treeset' or 'M${prefix_folder}_L00F00C00.treeset'

#for n in $(seq 20  29)
#do

    echo "for folder in Chain${n}_treesets
do
    prefix_folder=\$(basename \$folder _treesets)
    echo \$prefix_folder

    #Get the right folder
    for file in \${folder}/*.treeset
    do
        #In the right folder get the right trees
        prefix_file=\$(basename \$file .treeset)
        sh TEM_TreeCmp.sh \${prefix_folder}_treesets/\${prefix_folder}.True_tree.tre.nex \${file} 51 1000 \${prefix_file} #or \${prefix_folder}.True_tree.tre.nex M\${prefix_folder}_L00F00C00.treeset
    done
done " > TreeCmp_Chain${n}.task

done