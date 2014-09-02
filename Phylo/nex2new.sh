#!/bin/sh
##########################
#Transforms a chain of nexus trees into newick trees (or the other way around)
##########################
#SYNTAX :
#sh TEM_NexToNew.sh <input chain> <output format>
#with
#<input chain> the chain of trees to convert
#<output format> the new output format (either nexus or newick)
##########################
#version: 0.1
TEM_NexToNew_version="TEM_NexToNew.sh v0.1"
#----
#guillert(at)tcd.ie - 15/08/2014
##########################
#Requirements:
#-R 3
#-R package "ape"
##########################

#INPUT
chain=$1
format=$2

#R SETTINGS
echo "library(ape) ; trees<-list.files(pattern='$chain')
for (tree in 1:length(trees)) {" > TEM_NexToNew.R

if echo $format | grep 'newick' >/dev/null
then
    echo "write.tree(read.nexus(trees[tree]), file=paste(trees[tree], '.tre', sep='')) ; message('.', appendLF=FALSE) }" >> TEM_NexToNew.R
else
    echo "write.nexus(read.tree(trees[tree]), file=paste(trees[tree], '.nex', sep='')) ; message('.', appendLF=FALSE) }" >> TEM_NexToNew.R
fi

#CONVERTING THE FILES
length=$(ls *$chain* | wc -l | sed 's/[[:space:]]//g')
echo "Converting $length trees into $format."
R --no-save < TEM_NexToNew.R >/dev/null
rm TEM_NexToNew.R
echo "Done."

#End