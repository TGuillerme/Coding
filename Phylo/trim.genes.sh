#!/bin/sh

##########################
# Script for using trimming gene alignements
##########################
#SYNTAX: TreeCmp -i <input> -o <output> -s <start> -e <end>
#with:
#-i <input> the name of the input file
#-o <output> the name of the output file
#-s <start> the position from which to start
#-e <end> the position from which to end
##########################
#guillert(at)tcd.ie - 2018/01/16
##########################

#INPUT
## Input values
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input)
        INPUT="$2"
        shift
        ;;
    -o|--output)
        OUTPUT="$2"
        shift
        ;;
    -s|--start)
        START="$2"
        shift
        ;;
    -e|--end)
        END="$2"
        ;;
        *)

        ;;
esac
shift
done

## Set default arguments

## Number of taxa
NTAXA=$(grep '>' $INPUT | wc -l)

## Making the fasta into a single line
sed '/>/s/>/@>/g' $INPUT | sed '/>/s/$/@/g' | sed 's/\n//g' > output.trim.gene.tmp
tr -cd "[:print:]" < output.trim.gene.tmp > output2.trim.gene.tmp
tr '@' '\n'  < output2.trim.gene.tmp > output.trim.gene.tmp
sed '1d'  output.trim.gene.tmp > output2.trim.gene.tmp

## Splitting the taxa and data
grep '>' output2.trim.gene.tmp > taxa.trim.gene.tmp
grep -v '>' output2.trim.gene.tmp > data.trim.gene.tmp

## Removing the end of the sequence
cut -c$START-$END data.trim.gene.tmp > data2.trim.gene.tmp

## Recombining the files
paste -d"@" taxa.trim.gene.tmp data2.trim.gene.tmp > output.trim.gene.tmp
tr '@' '\n'  < output.trim.gene.tmp > $OUTPUT

## Remove temporaries
rm -f *.trim.gene.tmp

## end