#!/usr/bin/sh
##########################
#Swap the citation types between sysbiol or vancouver for a tex file
##########################
#SYNTAX:
#sh CiteChange.sh <tex>
#with:
#<tex> being the tex file
#########################
#version 0.1
#----
#guillert(at)tcd.ie - 18/05/2015
##########################

#Input value
tex=$1

#Check which is the actual cite format
if grep 'bibliographystyle' ${tex} | grep 'vancouver' > /dev/null
then
    #Change from vancouver to sysbio style
    echo "Changing bibliography style from 'vancouver' to 'sysbio'."
    sed 's/\\bibliographystyle{vancouver}/\\bibliographystyle{sysbio}/g' $tex > ${tex}.tmp

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #enable natbib
    header=$(grep -n "documentclass" $tex | sed 's/\:\\documentclass.*//g')
    sed ''"$header"'s/$/\'$'\n\\\usepackage\{natbib\}/' ${tex} > ${tex}.tmp

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #change all the \cite into \citep
    sed 's/\\cite/\\citep/g' ${tex} > ${tex}.tmp 

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #correct the e.g. (1)
    sed 's/(e.g. \\citep/\\citep[e.g.][]/g' ${tex} | sed 's/})/}/g' > ${tex}.tmp 

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #Done!
    echo "Remember to add the 'sysbio.bst' file to the current directory."

else
    #Change from sysbio to vancouver style
    echo "Changing bibliography style from 'sysbio' to 'vancouver':"
    sed 's/\\bibliographystyle{sysbio}/\\bibliographystyle{vancouver}/g' $tex > ${tex}.tmp

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #disable natbib
    sed 's/\\usepackage{natbib}//g' ${tex} > ${tex}.tmp

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #change all the \citep into \cite
    sed 's/\\citep/\\cite/g' ${tex} > ${tex}.tmp 

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #correct the e.g. (1)
    sed 's/\\cite\[e.g.\]\[\]/(e.g. \cite/g' ${tex} > ${tex}.tmp 

    #rename the tex file
    mv ${tex}.tmp ${tex}

    #Done!
    echo "Brackets not closed for examples."
    echo "Remember to add the 'vancouver.bst' file to the current directory."
fi