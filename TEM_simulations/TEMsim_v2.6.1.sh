##########################
#Total Evidence Method simulations
##########################
#SYNTAX:
#sh TEMsim [<Living Species> <Input Matrix>] <Molecular characters> <Morphological characters> <Evolutionary model> <Method> <Replicates> <Number of simulations> <Chain name> <Input>
#with:
#<Living Species> being any entire number of living species to put into the matrices.
#<Input matrix> an optional matrix in phylip format to split. The first species should be the outgroup.
#<Molecular characters> being any entire number of molecular characters to put into the matrices.
#<Morphological characters> being any entire number of morphological characters to put into the matrices. Is ignored if an input matrix is given.
#<Evolutionary model> can be chosen between HKY or GTR as an evolutionary model to build the matrices.
#<Method> can be chosen between ML or Bayesian.
#<Replicates> being either a number of bootstraps (if method is ML) or any entire number of mega generations (10e6) (if method is Bayesian).
#<Number of simulations> being any entire number of repetitions of the simulations.
#<Chain name> being any string of characters used as the chain name.
##########################
#Simulate a Total Evidence Method with variation of three parameters:
#-the number of morphologicaly coded "living" taxa (NL)
#-the number of missing data for the "fossil" taxa (NF)
#-the number of morphological data (NC)
#version: 2.6.1
#Update: Use the full random Matsim_v4.4.sh script
#Update: Only one outgroup is given in input
#Update: Uses MrBayes
#Update: Build a backbone for MrBayes input
#Update: Evolutionary model for the matrices construction and bayesian method for the tree construction can be precised
#Update: Possibility to chose the method between ML, CON (Bayesian with backbone) or TEM (Bayesian Total Evidence)
#Update: Optional Input file is accepted
#Update: The same evolutionary model is used for building the matrix through Marsim and to run the trees
#Update: Methods are now restricted to ML or Bayesian only
#Update: On the multicore version, the number of runs is fixed a 2 and the number of chains is fixed at 4
#Update: When using bayesian method, rates are fixed to Gamma distribution with a given alpha prior (normal distribution centered on the given mean with 5% sd) and 4 distinct gamma categories
#Update: When using bayesian method, the True tree topology simulated by Matsim si given as a starting tree for the mcmc
#Update: Code clean and tidied up
#----
#guillert(at)tcd.ie - 04/03/2014
##########################
#Requirements:
#-Shell script Matsim_v4.4.sh
#-R =< 3.0.1
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



#Step 1 - INPUT



#Input values
LivingSp=$1
MolecularChar=$2
MorphologChar=$3
Model=$4
Method=$5
NGen=$6
Simulations=$7
Chain=$8

#LivingSp=RonquistMat.phylip
#MolecularChar=354
#MorphologChar=5097
#Model=HKY
#Method=ML
#NGen=100
#Simulations=1
#Chain=test


#Testing if input matrix is provided
input=$LivingSp
if [ "$input" -eq "$input" ] 2>/dev/null; then
  Input='NOINPUT'
else
  Input=$LivingSp
fi



echo $Input > input.test
if grep "NOINPUT" input.test > /dev/null
then
    FossilSp=$LivingSp
    TotalSp=$LivingSp
    let "TotalSp += $FossilSp"
    let "TotalSp +=1"
else 
    TotalSp=$LivingSp
fi
rm input.test



#Creating secondary inputs
TotalChar=$MolecularChar
let "TotalChar += $MorphologChar"
Nburn=$NGen
let "Nburn *=10"
let "Nburn *=25"
let "Nburn /=100"



#Initializing the loop
for n in $(seq 1 $Simulations)
do
    #Creates the output folder
    mkdir ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}_run${n}
    #Copy the Matsim script
    sed 's/MM/C'"${n}"'/g' Matsim_v4.4.sh | sed 's/Simulation_@.log/Simulation_'"${n}"'.log/g' > ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}_run${n}/Matsim_${n}.sh
    cp $Input ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}_run${n}/
    cd ${TotalSp}t_${TotalChar}c_${Model}_${Method}_${Chain}_run${n}



    #Step 2 - CREATING THE MATRICES



    #Using the Matsim script to create the matrices
    sh Matsim_${n}.sh $LivingSp $MolecularChar $MorphologChar $Model $Input
    rm C${n}F.phylip
    rm C${n}L.phylip
    rm MatBas.phylip
    FossilSp=$(sed -n '8p' Simulation_${n}.log | sed 's/Number of fossil species = //g')
    TotalSp=$FossilSp
    let "TotalSp += $LivingSp"

    

    #Setting the evolutionary model for MrBayes
    if grep 'Chosen model = HKY'  Simulation_${n}.log > /dev/null
    then
        nst=$'2'
    else
        nst=$'6'
    fi



    #Transforming into nexus format
    for f in *.phylip
    do
        echo $f
        seqConverter.pl -d${f} -ip -on -ri
    done



    #Adding header in the nexus file
    if grep "Input matrix:"  Simulation_${n}.log > /dev/null
    then
        MolChar1=$MolecularChar
        MorphoChar00=$MorphologChar
        let "MorphoChar00 += $MolecularChar"
        let "MorphoChar00 += 1"
        MorphoChar10=$MorphologChar
        let "MorphoChar10 *= 90"
        let "MorphoChar10 /= 100"
        let "MorphoChar10 += $MolecularChar"
        let "MorphoChar10 += 1"
        MorphoChar25=$MorphologChar
        let "MorphoChar25 *= 75"
        let "MorphoChar25 /= 100"
        let "MorphoChar25 += $MolecularChar"
        let "MorphoChar25 += 1"
        MorphoChar50=$MorphologChar
        let "MorphoChar50 *= 50"
        let "MorphoChar50 /= 100"
        let "MorphoChar50 += $MolecularChar"
        let "MorphoChar50 += 1"
        MorphoChar75=$MorphologChar
        let "MorphoChar75 *= 25"
        let "MorphoChar75 /= 100"
        let "MorphoChar75 += $MolecularChar"
        let "MorphoChar75 += 1"
    else
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
    fi


    #Convert to nexus
    if grep 'Chosen method = Bayesian' Simulation_${n}.log > /dev/null
    then

        #format datatype = protein
        if grep 'format datatype = protein' C${n}_L00F00C00.nex > /dev/null
            then
            echo 'datatype = protein' > /dev/null
            for f in *C00.nex
            do
                sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C10.nex
            do
                sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C25.nex
            do
                sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C50.nex
            do
                sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C75.nex
            do
                sed 's/format datatype = protein gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
        else 
            echo 'datatype = protein' > /dev/null
        fi



        #format datatype = nculeotide
        if grep 'format datatype = nucleotide' C${n}_L00F00C00.nex > /dev/null
        then
            echo 'datatype = nucleotide' > /dev/null
            for f in *C00.nex
            do
                sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C10.nex
            do 
                sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C25.nex
            do
                sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C50.nex
            do
                sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C75.nex
            do
                sed 's/format datatype = nucleotide gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
        else
            echo 'datatype = nucleotide' > /dev/null
        fi



        #format datatype = DNA
        if grep 'format datatype = DNA' C${n}_L00F00C00.nex > /dev/null
        then
            echo 'datatype = DNA' > /dev/null
            for f in *C00.nex
            do
                sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar00"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C10.nex
            do 
                sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar10"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C25.nex
            do
                sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar25"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C50.nex
            do
                sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar50"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
            for f in *C75.nex
            do
                sed 's/format datatype = DNA gap = - missing = ?;/format datatype=mixed(DNA:1-'"$MolecularChar"',standard:'"$MolChar1"'-'"$MorphoChar75"') interleave=yes gap=- missing=?;/g' $f > ${f}.tmp
            done
        else
            echo 'datatype = DNA' > /dev/null
        fi

    else
        echo "nothing" > /dev/null
    fi

    #Checking the seqConverter.pl names bug
    rm *.nex
    for f in *.nex.tmp
    do
        prefix=$(basename $f .nex.tmp)
        echo converting ${prefix}
        sed 's/\([a-z]\)\([a-z]\)\([0-9]\)\([0-9]\)\([0-9]\)[[:space:]]\([0-9]\)/\1\2\3\4\5\6	/g' $f > ${prefix}.nex
    done
    rm *.nex.tmp



    #Step 3 - RUNNING THE TREES



    #Setting the outgroup
    sed -n '2p' C${n}_L00F00C00.phylip > out.tmp
    sed -n 's/\([A-z]*\)[[:space:]].*/\1/p' out.tmp > outgroup.txt
    rm *.tmp
    outgroup=$(sed -n '1p' outgroup.txt)



    echo "##########################" >> Simulation_${n}.log
    echo "TREE BUILDING" >> Simulation_${n}.log
    echo "Chosen method = $Method" >> Simulation_${n}.log
    echo "Outgroup = $outgroup" >> Simulation_${n}.log
    echo "==========================" >> Simulation_${n}.log



    #Which method?
    if grep 'Chosen method = ML' Simulation_${n}.log > /dev/null
    then
        #Running ML trees with fast bootstraps
        echo 'ML method'



        #Creating partition set
        echo "DNA, set1 = 1-$MolecularChar" > part00.set
        echo "MULTI, set2 = $MolChar1 - $MorphoChar00" >> part00.set
        echo "DNA, set1 = 1-$MolecularChar" > part10.set
        echo "MULTI, set2 = $MolChar1 - $MorphoChar10" >> part10.set
        echo "DNA, set1 = 1-$MolecularChar" > part25.set
        echo "MULTI, set2 = $MolChar1 - $MorphoChar25" >> part25.set
        echo "DNA, set1 = 1-$MolecularChar" > part50.set
        echo "MULTI, set2 = $MolChar1 - $MorphoChar50" >> part50.set
        echo "DNA, set1 = 1-$MolecularChar" > part75.set
        echo "MULTI, set2 = $MolChar1 - $MorphoChar75" >> part75.set



        #Running the trees
        echo "TIMER:" >> Simulation_${n}.log
        echo "ML trees building: START" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C00.phylip
        do
            prefix=$(basename $f .phylip)
            echo ${prefix}
            raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $f -n ${prefix} -m GTRGAMMA -q part00.set -o $outgroup -x 12345 -# $NGen
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done

        echo "ML trees building: 20%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C10.phylip
        do
            prefix=$(basename $f .phylip)
            echo ${prefix}
            raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $f -n ${prefix} -m GTRGAMMA -q part10.set -o $outgroup -x 12345 -# $NGen
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done

        echo "ML trees building: 40%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C25.phylip
        do
            prefix=$(basename $f .phylip)
            echo ${prefix}
            raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $f -n ${prefix} -m GTRGAMMA -q part25.set -o $outgroup -x 12345 -# $NGen
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log ; done ;

        echo "ML trees building: 60%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C50.phylip
        do
            prefix=$(basename $f .phylip)
            echo ${prefix}
            raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $f -n ${prefix} -m GTRGAMMA -q part50.set -o $outgroup -x 12345 -# $NGen
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done

        echo "ML trees building: 80%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C75.phylip
        do
            prefix=$(basename $f .phylip)
            echo ${prefix}
            raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $f -n ${prefix} -m GTRGAMMA -q part75.set -o $outgroup -x 12345 -# $NGen
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done

        echo "ML trees building: 100%" >> Simulation_${n}.log
        date >> Simulation_${n}.log


        
    else
        #Running Bayesian trees with start tree
        echo 'Bayesian method'



        #Starting tree
        echo "library(ape)
            Start.tree<-read.tree('True_tree.tre')
            Start.tree[[4]]<-rep(1, length(Start.tree[[4]]))
            write.tree(Start.tree, 'Start_tree.tre')" | R --no-save

        StartTree=$(sed -n '1p' Start_tree.tre | sed 's/:1):0;/:1);/g')

        for f in *.nex
        do
            echo "Begin trees;
            ">> ${f}
            echo "tree Start_tree = $StartTree
            ">> ${f}
            echo "End;">> ${f}
        done



        #Preparing the MrBayes priors
        if grep 'Chosen model = HKY'  Simulation_${n}.log > /dev/null
        then
            #Using HKY model
            DGamma=$(grep "Molecular rates distribution" Simulation_${n}.log | sed 's/Molecular rates distribution (gamma) alpha = //g')
            DGammaInf=$(echo "$DGamma*0.9" | bc -l)
            DGammaSup=$(echo "$DGamma*1.1" | bc -l)

            MGamma=$(grep "Morphological rates distribution" Simulation_${n}.log | sed 's/Morphological rates distribution (gamma) alpha = //g')
            MGammaInf=$(echo "$MGamma*0.9" | bc -l)
            MGammaSup=$(echo "$MGamma*1.1" | bc -l)
        else
            #Using GTR model
            DGamma=$(grep "Molecular rates distribution" Simulation_${n}.log | sed 's/Molecular rates distribution (gamma) alpha = //g')
            DGammaInf=$(echo "$DGamma*0.9" | bc -l)
            DGammaSup=$(echo "$DGamma*1.1" | bc -l)

            MGamma=$(grep "Morphological rates distribution" Simulation_${n}.log | sed 's/Morphological rates distribution (gamma) alpha = //g')
            MGammaInf=$(echo "$MGamma*0.9" | bc -l)
            MGammaSup=$(echo "$MGamma*1.1" | bc -l)
        fi



        #Creating the mrbayes.cmd template
        echo "begin mrbayes;" > base-cmd.tmp
        echo "[Data input]" >> base-cmd.tmp
        echo "set autoclose=yes nowarn=yes;" >> base-cmd.tmp
        echo "log start filename=<CHAIN>.log;" >> base-cmd.tmp
        echo "execute <CHAIN>.nex;" >> base-cmd.tmp
        echo "charset DNA = 1-$MolecularChar;" >> base-cmd.tmp
        echo "charset morphology = $MolChar1-<NUMBER>;" >> base-cmd.tmp
        echo "partition favored = 2: DNA, morphology;" >> base-cmd.tmp
        echo "set partition = favored;" >> base-cmd.tmp
        echo "" >> base-cmd.tmp
        echo "[Model settings]" >> base-cmd.tmp
        echo "outgroup $outgroup ;" >> base-cmd.tmp
        echo "prset applyto=(1) Shapepr=Exponential($DGamma) Tratiopr = beta(80,40);" >> base-cmd.tmp
        echo "prset applyto=(2) Shapepr=Exponential($MGamma);" >> base-cmd.tmp
        echo "lset applyto=(1) nst=${nst} rates=gamma Ngammacat=4;" >> base-cmd.tmp
        echo "lset applyto=(2) nst=1 rates=gamma Ngammacat=4;" >> base-cmd.tmp
        echo "" >> base-cmd.tmp
        echo "[MCMC settings]" >> base-cmd.tmp
        echo "startvals tau=Start_tree V=Start_tree ;" >> base-cmd.tmp
        echo "mcmc nruns=2 Nchains=4 ngen=${NGen}000000 samplefreq=10000 printfreq=50000 diagnfreq=500000 Stoprule=YES stopval=0.01 mcmcdiagn=YES file=<CHAIN>;" >> base-cmd.tmp
        echo "sump Filename=<CHAIN>;" >> base-cmd.tmp
        echo "sumt Filename=<CHAIN> burnin=${Nburn}00000;" >> base-cmd.tmp
        echo "end;" >> base-cmd.tmp



        #Giving the <CHAIN> and <NUMBER> arguments function of the input nexus file
        for f in *C00.nex
        do
            prefix=$(basename $f .nex)
            sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar00"'/g' > ${prefix}.cmd
        done
        for f in *C10.nex
        do
            prefix=$(basename $f .nex)
            sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar10"'/g' > ${prefix}.cmd
        done
        for f in *C25.nex
        do
            prefix=$(basename $f .nex)
            sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar25"'/g' > ${prefix}.cmd
        done
        for f in *C50.nex
        do
            prefix=$(basename $f .nex)
            sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar50"'/g' > ${prefix}.cmd
        done
        for f in *C75.nex
        do
            prefix=$(basename $f .nex)
            sed 's/<CHAIN>/'"${prefix}"'/g' base-cmd.tmp | sed 's/<NUMBER>/'"$MorphoChar75"'/g' > ${prefix}.cmd
        done



        #launch the MCMC
        echo "Bayesian trees building: START" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C00.cmd
        do
            mb $f
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done
        echo "Bayesian trees building: 20%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C10.cmd
        do
            mb $f
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done
        echo "Bayesian trees building: 40%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C25.cmd
        do
            mb $f
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done
        echo "Bayesian trees building: 60%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C50.cmd
        do
            mb $f
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done
        echo "Bayesian trees building: 80%" >> Simulation_${n}.log
        date >> Simulation_${n}.log

        for f in *C75.cmd
        do
            mb $f
            echo "$f completed" >> Simulation_${n}.log
            date >> Simulation_${n}.log
        done
        cho "Bayesian trees building: 100%" >> Simulation_${n}.log
        date >> Simulation_${n}.log



    fi



    #Step 4 - CLEANING



    #Matsim
    tar cf Matrices_archives.tar *.phylip *.nex
    rm *.phylip ; rm *.nex
    cp Simulation_$n.log Simulation_$n.tmp
    rm Simulation_$n.log
    tar cf Log_archives.tar *.log
    rm *.log
    cp Simulation_$n.tmp Simulation_$n.log
    rm Simulation_$n.tmp

    #MrBayes
    cp bayes_constraint.nex bayes_constraint.tmp
    rm bayes_constraint.nex
    cp bayes_constraint.tmp bayes_constraint.nex
    rm bayes_constraint.tmp
    tar cf MrBayesCmd_archives.tar *.cmd
    rm *.cmd
    tar cf ConsensusTrees_archives.tar *.con.tre
    rm *.con.tre
    tar cf Trees_archives.tar *.run1.t *.run2.t *.tstat
    rm *.run1.t ; rm *.run2.t ; rm *.tstat
    tar cf Param_archives.tar *.run1.p *.run2.p *.pstat
    rm *.run1.p ; rm *.run2.p ; rm *.pstat
    tar cf MrBayes_archives.tar *.ckp *.ckp~ *.lstat *.mcmc *.parts *.trprobs *.vstat
    rm *.ckp ; rm *.ckp~ ; rm *.lstat ; rm *.mcmc ; rm *.parts ; rm *.trprobs ; rm *.vstat
 
    #RAxML
    tar cf Runnings_archives.tar *_run*
    rm *_run*
    tar cf RAxML.tar RAxML_*
    rm RAxML_*
    tar ML_trees.tar RAxML_bipartitions.*
    tar cf RAxML_archives.tar RAxML_*
    rm RAxML_*


    #Closing the loop
    cd ..
done



#end