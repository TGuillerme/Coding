##########################
#Total Evidence Method simulations
##########################
#SYNTAX:
#sh TEM_cluster-tasker.sh <Living Species> <Molecular characters> <Morphological characters> <Chain name>
#with:
#<Living Species> being any entire number of living species to put into the matrices
#<Molecular characters> being any entire number of molecular characters to put into the matrices
#<Morphological characters> being any entire number of morphological characters to put into the matrices
#<Chain name>
#Notes on TEM_treesim.sh options:
#<Evolutionary model> as been fixed as 'HKY'
#<Method> as been fixed as 'Bayesian'
#<Generations millions> as been fixed as '50'
#<Number of simulations> as been fixed as '1'
#<Input> option as been disabled
#<Generations millions> being any entire number of mega generations (10e6)
#<Number of simulations> being any entire number of repetitions of the simulations
#<Chain name> being any string of characters used as the chain name
##########################
#Simulate a Total Evidence Method with variation of three parameters:
#-the number of morphologicaly coded "living" taxa (NL)
#-the number of missing data for the "fossil" taxa (NF)
#-the number of morphological data (NC)
#
#GENERATOR v1.3
#Generates the matrices and the MrBayes command lines
#########################
#Update: Now creates a series of 25 jobs to submit
#Update: Updated to use mpirun on the new version of the cluster
#
#----
#guillert(at)tcd.ie - 17/04/2014
##########################
#Requirements:
#-Shell script TEM_matsim.sh
#-R 3.0.1
#-R package "ape"
#-R package "phyclust"
#-R package "diversitree"
#-R pachage "MCMCpack"
#-seqConverter.pl
#-raxmlHPC-PTHREADS-SSE3
#-MrBayes 3.2.2
##-To install the R packages:
##echo "install.packages(c('ape','phyclust','diversitree','MCMCpack'))" | R --no-save
##########################



#Step 0 - INPUT

#Input values
LivingSp=$1
MolecularChar=$2
MorphologChar=$3
Model=$'HKY' #Fixed
Method=$'Bayesian' #Fixed
NGen=$'50' #Fixed
Simulations=$'1' #Fixed
Chain=$4

#Creating secondary inputs
FossilSp=$LivingSp
TotalSp=$LivingSp

let "TotalSp += $FossilSp"
let "TotalSp +=1"

TotalChar=$MolecularChar
let "TotalChar += $MorphologChar"
Nburn=$NGen
let "Nburn *=10"
let "Nburn *=25"
let "Nburn /=100"

#Step 1 - INPUT

#Initializing the loop
for n in $(seq 1 $Simulations)
do
    mkdir ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}
    ed 's/MM/M'"${Chain}"'/g' TEM_matsim.sh | sed 's/Simulation_@.log/Simulation_'"${Chain}"'.log/g' > ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}/Matsim_${Chain}.sh
    cp $Input ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}/
    cd ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}

#Step 2 - PREPARING THE MATRICES FOR MRBAYES

#Step 2.1 - Using the Matsim script to create the 125 submatrices

    sh Matsim_${Chain}.sh $LivingSp $MolecularChar $MorphologChar $Model $Input
    rm M${Chain}F.phylip
    rm M${Chain}L.phylip
    rm MatBas.phylip
    FossilSp=$(sed -n '8p' Simulation_${Chain}.log | sed 's/Number of fossil species = //g')
    TotalSp=$FossilSp
    let "TotalSp += $LivingSp"

    if grep 'Matrices generated using Matsim_v4.3.1' Simulation_${Chain}.log > /dev/null

    then
        let "TotalSp +=1"
    else
        echo "nothing" > /dev/null
    fi

    #Setting the evolutionary model for MrBayes
    if grep 'Chosen model = HKY'  Simulation_${Chain}.log > /dev/null
    then nst=$'2'
    else nst=$'6'
    fi

    #Step 2.2 - Transforming the phylip submatrices in nexus ones

    #Transforming the phylip in nexus
    for f in *.phylip
    do
        echo $f
        seqConverter.pl -d${f} -ip -on -ri
    done

    #Adding header in the nexus file
    MolChar1=$MolecularChar
    let "MolChar1 += 1"
    MorphoChar00=$MorphologChar
    let "MorphoChar00 += $MolecularChar"
    MorphoChar10=$MorphologChar
    let "MorphoChar10 *= 90"
    let "MorphoChar10 /= 100"
    let "MorphoChar10 += $MolecularChar"
    MorphoChar25=$MorphologChar
    let "MorphoChar25 *= 75"
    let "MorphoChar25 /= 100"
    let "MorphoChar25 += $MolecularChar"
    MorphoChar50=$MorphologChar
    let "MorphoChar50 *= 50"
    let "MorphoChar50 /= 100"
    let "MorphoChar50 += $MolecularChar"
    MorphoChar75=$MorphologChar
    let "MorphoChar75 *= 25"
    let "MorphoChar75 /= 100"
    let "MorphoChar75 += $MolecularChar"

    #IF STATEMENT = 'datatype = protein' - IF-BLOCK a START
    if grep 'format datatype = protein' M${Chain}_L00F00C00.nex > /dev/null
    then echo 'datatype = protein' > /dev/null

        for f in *C00.nex ; do sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C10.nex ; do sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C25.nex ; do sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C50.nex ; do sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C75.nex ; do sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;

    else echo 'datatype = protein' > /dev/null
    fi

    #IF STATEMENT = 'datatype = nucleotide' - IF-BLOCK b START
    if grep 'format datatype = nucleotide' M${Chain}_L00F00C00.nex > /dev/null
    then echo 'datatype = nucleotide' > /dev/null

        for f in *C00.nex ; do sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C10.nex ; do sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C25.nex ; do sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C50.nex ; do sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C75.nex ; do sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;

    else echo 'datatype = nucleotide' > /dev/null
    fi

    #IF STATEMENT = 'datatype = DNA' - IF-BLOCK c START
    if grep 'format datatype = DNA' M${Chain}_L00F00C00.nex > /dev/null
    then echo 'datatype = DNA' > /dev/null

        for f in *C00.nex ; do sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C10.nex ; do sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C25.nex ; do sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C50.nex ; do sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;
        for f in *C75.nex ; do sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp ; done ;

    else echo 'datatype = DNA' > /dev/null
    fi

    rm *.nex

    #Checking the seqConverter.pl names bug
    for f in *.nex.tmp
    do
        prefix=$(basename $f .nex.tmp)
        echo converting ${prefix}
        sed 's/\([a-z]\)\([a-z]\)\([0-9]\)\([0-9]\)\([0-9]\)[[:space:]]\([0-9]\)/\1\2\3\4\5\6	/g' $f > ${prefix}.nex
    done
    rm *.nex.tmp

    #Step 3 - CREATING THE MOLECULAR BACKBONE TREE

    #Step 3.1 - Creating the molecular backbone tree with RaXML

    #Setting the outgroup
    #sed -n '2p' M${Chain}_L00F00C00.phylip > out.tmp
    #sed -n 's/\([A-z]*\)[[:space:]].*/\1/p' out.tmp > outgroup.txt
    #rm *.tmp
    #outgroup=$(sed -n '1p' outgroup.txt)
    outgroup=$'sp0000'

    echo "##########################

    TREE BUILDING

    Chosen method = $Method
    Outgroup = $outgroup
    ==========================" >> Simulation_${Chain}.log

    #IF STATEMENT = 'Chosen method: ML' - IF-BLOCK 1 START
    if grep 'Chosen method = ML'  Simulation_${Chain}.log > /dev/null
    then echo 'ML method'

    #[ML method, fully available in TEM_treesim]

    else echo 'Bayesian method'

        #Step 3b - PREPARING THE MCMC STARTING TREE FROM THE TRUE TREE - Removing the branch length from the true tree
        echo "library(ape) ; Start.tree<-read.tree('True_tree.tre') ; Start.tree[[4]]<-rep(1, length(Start.tree[[4]])) ; write.tree(Start.tree, 'Start_tree.tre')" | R --no-save
        StartTree=$(sed -n '1p' Start_tree.tre | sed 's/:1):0;/:1);/g' )

        for f in *.nex
        do echo "
Begin trees;

tree Start_tree = $StartTree

End;" >> ${f}
        done

        #Step 4 - PREPARING THE MRBAYES CODE

        #Step 4.1 - Create the default command - <CHAIN> and <NUMBER> arguments to fil in

        #mrbayes.cmd template
        echo "begin mrbayes;
[Data input]
set autoclose=yes nowarn=yes;
log start filename=<CHAIN>.log;
execute <CHAIN>.nex;
charset DNA = 1-$MolecularChar;
charset morphology = $MolChar1-<NUMBER>;
partition favored = 2: DNA, morphology;
set partition = favored;

[Model settings]
outgroup $outgroup ;
prset applyto=(1) Shapepr=Exponential(0.5) Tratiopr = beta(80,40);
prset applyto=(2) Shapepr=Exponential(0.5);
lset applyto=(1) nst=${nst} rates=gamma Ngammacat=4;
lset applyto=(2) nst=1 rates=gamma Ngammacat=4;

[MCMC settings]
startvals tau=Start_tree V=Start_tree ;
mcmc nruns=2 Nchains=4 ngen=${NGen}000000 samplefreq=10000 printfreq=50000 diagnfreq=500000 Stoprule=YES stopval=0.01 mcmcdiagn=YES file=<CHAIN>;
sump Filename=<CHAIN> Relburnin=YES Burninfrac=0.25;
sumt Filename=<CHAIN> Relburnin=YES Burninfrac=0.25;
end;" > base-cmd.tmp

        #Step 5 - RUN MR BAYES

        #Step 5.1 - Preparing the commands

        #Giving the <CHAIN> and <NUMBER> arguments function of the input nexus file
        for f in *C00.nex ; do prefix=$(basename $f .nex) ; sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar00"'/g' > ${prefix}.cmd ; done ;
        for f in *C10.nex ; do prefix=$(basename $f .nex) ; sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar10"'/g' > ${prefix}.cmd ; done ;
        for f in *C25.nex ; do prefix=$(basename $f .nex) ; sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar25"'/g' > ${prefix}.cmd ; done ;
        for f in *C50.nex ; do prefix=$(basename $f .nex) ; sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar50"'/g' > ${prefix}.cmd ; done ;
        for f in *C75.nex ; do prefix=$(basename $f .nex) ; sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar75"'/g' > ${prefix}.cmd ; done ;

        rm base-cmd.tmp
    fi

    #Step 5.2 - Run MrBayes

    #Generate the task files (8 cores)
    echo "/#!/bin/sh
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -p compute
#SBATCH -J ${Chain}
source /etc/profile.d/modules.sh
export http_proxy=http://proxy.tchpc.tcd.ie:8080
module load cports6 openmpi/1.6.5-gnu4.8.2 
module load cports gsl/1.16-gnu

##########################
#TASK FILE <N> - <CHAIN>
##########################

mpirun -np 8 mb M${Chain}_L00<CHAIN>.cmd ;
mpirun -np 8 mb M${Chain}_L10<CHAIN>.cmd ;
mpirun -np 8 mb M${Chain}_L25<CHAIN>.cmd ;
mpirun -np 8 mb M${Chain}_L50<CHAIN>.cmd ;
mpirun -np 8 mb M${Chain}_L75<CHAIN>.cmd ;" | sed 's/\/\#!/\#!/g' > job.template


    sed 's/<N>/1/g' job.template | sed 's/<CHAIN>/F00C00/g' > ${Chain}_1.sh.job
    sed 's/<N>/2/g' job.template | sed 's/<CHAIN>/F10C00/g' > ${Chain}_2.sh.job
    sed 's/<N>/3/g' job.template | sed 's/<CHAIN>/F25C00/g' > ${Chain}_3.sh.job
    sed 's/<N>/4/g' job.template | sed 's/<CHAIN>/F50C00/g' > ${Chain}_4.sh.job
    sed 's/<N>/5/g' job.template | sed 's/<CHAIN>/F75C00/g' > ${Chain}_5.sh.job
    sed 's/<N>/6/g' job.template | sed 's/<CHAIN>/F00C10/g' > ${Chain}_6.sh.job
    sed 's/<N>/7/g' job.template | sed 's/<CHAIN>/F10C10/g' > ${Chain}_7.sh.job
    sed 's/<N>/8/g' job.template | sed 's/<CHAIN>/F25C10/g' > ${Chain}_8.sh.job
    sed 's/<N>/9/g' job.template | sed 's/<CHAIN>/F50C10/g' > ${Chain}_9.sh.job
    sed 's/<N>/10/g' job.template | sed 's/<CHAIN>/F75C10/g' > ${Chain}_10.sh.job
    sed 's/<N>/11/g' job.template | sed 's/<CHAIN>/F00C25/g' > ${Chain}_11.sh.job
    sed 's/<N>/12/g' job.template | sed 's/<CHAIN>/F10C25/g' > ${Chain}_12.sh.job
    sed 's/<N>/13/g' job.template | sed 's/<CHAIN>/F25C25/g' > ${Chain}_13.sh.job
    sed 's/<N>/14/g' job.template | sed 's/<CHAIN>/F50C25/g' > ${Chain}_14.sh.job
    sed 's/<N>/15/g' job.template | sed 's/<CHAIN>/F75C25/g' > ${Chain}_15.sh.job
    sed 's/<N>/16/g' job.template | sed 's/<CHAIN>/F00C50/g' > ${Chain}_16.sh.job
    sed 's/<N>/17/g' job.template | sed 's/<CHAIN>/F10C50/g' > ${Chain}_17.sh.job
    sed 's/<N>/18/g' job.template | sed 's/<CHAIN>/F25C50/g' > ${Chain}_18.sh.job
    sed 's/<N>/19/g' job.template | sed 's/<CHAIN>/F50C50/g' > ${Chain}_19.sh.job
    sed 's/<N>/20/g' job.template | sed 's/<CHAIN>/F75C50/g' > ${Chain}_20.sh.job
    sed 's/<N>/21/g' job.template | sed 's/<CHAIN>/F00C75/g' > ${Chain}_21.sh.job
    sed 's/<N>/22/g' job.template | sed 's/<CHAIN>/F10C75/g' > ${Chain}_22.sh.job
    sed 's/<N>/23/g' job.template | sed 's/<CHAIN>/F25C75/g' > ${Chain}_23.sh.job
    sed 's/<N>/24/g' job.template | sed 's/<CHAIN>/F50C75/g' > ${Chain}_24.sh.job
    sed 's/<N>/25/g' job.template | sed 's/<CHAIN>/F75C75/g' > ${Chain}_25.sh.job

    rm job.template

    echo 'Next step:'
    echo "for f in *.sh.job ; do echo $f ; sbatch $f ; done ;"

    cd ..
done

#end