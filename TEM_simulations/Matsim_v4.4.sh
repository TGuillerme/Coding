##########################
#Simulating 125 phylip files with variation in molecular and morphological data
##########################
#SYNTAX :
#sh Matsim [<Living Species> <Input matrix>] <Molecular characters> <Morphological characters> <Evolutionary model>
#with
#<Living Species> being any entire number of living species to put into the matrices.
#<Input matrix> an optional matrix in phylip format to split. The first species should be the outgroup.
#<Molecular characters> being any entire number of molecular characters to put into the matrices.
#<Morphological characters> being any entire number of morphological characters to put into the matrices.
#<Evolutionary model> can be chosen between HKY or GTR as an evolutionary model to build the matrices. Is ignored if an input matrix is given.
##########################
#Simulating 125 phylip matrices with combined "fossil" and "living" taxa and insert gaps in the morphological part according to 3 parameters with 5 states each.
#version: 4.4
#Update: The tree shape and the fossil/living distribution is fully randomized.
#Update: Both parts of the complete matrix (morphological and molecular) are created independently from the same given "true" tree.
#Update: The HKY model is used for the molecular matrix and a ARD for the morphological one.
#Update: All variables are given in input
#Update: Possibility of using a GTR or HKY model for building the matrix.
#Update: Possibility to put an input matrix and skip the matrix building step.
#Update: Minimal amount of species deleted in NL parameter is always 1.
#Update: Minimal data per tip is set to 5% on the fossils. If the matrix is generated, the random gapping process is reran for all the species if any fossil is <5% data, if the matrix is given in input, fossils with <5% data are pruned.
#Update: A true outgroup is now generate with the True tree
#Update: Characters (molecular and morphological) rates distribution is now fixed
#Update: Morphological characters states frequencies is now fixed to 2 and 3 states only
#Update: A fixed gamma distribution is now given to the morphological characters rates
#Update: Morphological states changing rates are now all equal per character (changing model from ARD to ER)
#Update: Morphological characters are now saved in a MorphologicalRates.txt file
#Update: The NC gapping parameters states as been changed to 0, 10, 25, 50, 75 of deleted characters (instead of previous 0, 20, 40, 60, 80)
#Update: Does not allow invariant morphological characters in the morphological matrix
#Update: Code tidied up and cleaned
#Update: Improvement in the input matrix reading 
#----
#guillert(at)tcd.ie - 05/03/2014
##########################
#Requirements:
#-seqConverter.pl
#-R =< 3.0.1
#-R package "ape"
#-R package "phyclust"
#-R package "diversitree"
#-R pachage "MCMCpack"
##-To install the R packages:
##echo "install.packages(c('ape','phyclust','diversitree','MCMCpack'))" | R --no-save
##########################



#Step 1 - INPUT



#Input values
LivingSp=$1
MolecularChar=$2
MorphologChar=$3
Model=$4



#Testing if input matrix is provided
test=$LivingSp

if [ "$test" -eq "$test" ] 2>/dev/null; then
  Input='NOINPUT'
else
  Input=$LivingSp
fi



#Writing the header
echo $Input > test.tmp
if grep "NOINPUT" test.tmp > /dev/null
then
    echo "No input" > /dev/null
else 
    dat=$(date)
    echo "##########################" > Simulation_@.log
    echo "Input matrix: $Input" >> Simulation_@.log
    echo "$dat" >> Simulation_@.log
fi



#Transforming the input matrix
if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    echo 'Using the input matrix'



    #Creating the input attributes
    seqConverter.pl -d${Input} -ope -ri > seqConverter.log
    sed -n '2p' seqConverter.log | sed 's/Converting file //g' | sed 's/.phylip \.\.\.//g' > matrix.name
    Inprefix=$(sed -n '1p' matrix.name)
    rm seqConverter.log
    rm matrix.name
    rm ${Inprefix}_new.phylip



    #Reasigning the values of characters and species
    TotalSp=$(sed -n '1p' $Inprefix.phylip | sed 's/\(.\)[[:space:]]\(.\)*/\1/g') 
    TotalChar=$(sed -n '1p' $Inprefix.phylip | sed 's/'"${TotalSp}"'[[:space:]]//g') 
    LivingSp=$(grep "[[:space:]].*[ACGT]" $Inprefix.phylip | wc -l | sed 's/[[:space:]]//g')
    FossilSp=$TotalSp
    let "FossilSp -= $LivingSp"



    #Printing the parameters
    echo "PARAMETERS:" >> Simulation_@.log
    echo "Number of living species = $LivingSp" >> Simulation_@.log
    echo "Number of fossil species = $FossilSp" >> Simulation_@.log
    echo "Number of molecular characters = $MolecularChar" >> Simulation_@.log
    echo "Number of morphological characters = $MorphologChar" >> Simulation_@.log



    #Separating the species and the data


#TESTER separate the name/data sets


    #Living (select the species with DNA)
    grep "[[:space:]].*[ACGT]" $Inprefix.phylip | sed -n '2,$p' > living.subset
    sed 's/[[:space:]].*//g' living.subset > living.name
    sed 's/.*[[:space:]]//g' living.subset > living.data



    #Fossil (select the species without DNA (-v))
    grep -v "[[:space:]].*[ACGT]" $Inprefix.phylip | sed -n '2,$p' > fossil.subset
    sed 's/[[:space:]].*//g' fossil.subset > fossil.name
    sed 's/.*[[:space:]]//g' fossil.subset > fossil.data



    #Outgroup (the first group of the list)
    grep "[[:space:]].*[ACGT]" $Inprefix.phylip | sed -n '1p' > outgroup.subset
    sed 's/[[:space:]].*//g' outgroup.subset > outgroup.name
    sed 's/.*[[:space:]]//g' outgroup.subset > outgroup.data



    #Names list
    echo "taxa" > input.name
    cat outgroup.name living.name fossil.name >> input.name



    #Test if the morphological matrix is before or after the molecular one
    cut -c1-$MorphologChar living.data > morpho.test
    if grep "[0-9]" morpho.test > /dev/null
    then #Morphological matrix comes first



        #Living
        let "MorphologChar -= 1"
        cut -c1-$MorphologChar living.data > living_morpholog.data
        let "MorphologChar += 1"
        let "MolecularChar += $MorphologChar"
        cut -c$MorphologChar-$MolecularChar living.data > living_molecular.data
        let "MolecularChar -= $MorphologChar"


        #Fossil
        let "MorphologChar -= 1"
        cut -c1-$MorphologChar fossil.data > fossil_morpholog.data
        let "MorphologChar += 1"
        let "MolecularChar += $MorphologChar"
        cut -c$MorphologChar-$MolecularChar fossil.data > fossil_molecular.data
        let "MolecularChar -= $MorphologChar"



        #Outgroup
        let "MorphologChar -= 1"
        cut -c1-$MorphologChar outgroup.data > outgroup_morpholog.data
        let "MorphologChar += 1"
        let "MolecularChar += $MorphologChar"
        cut -c$MorphologChar-$MolecularChar outgroup.data > outgroup_molecular.data
        let "MolecularChar -= $MorphologChar"


    rm morpho.test
    else #Molecular matrix comes first



        #Living
        let "MolecularChar -= 1"
        cut -c1-$MolecularChar living.data > living_molecular.data
        let "MolecularChar += 1"
        let "MorphologChar += $MolecularChar"
        cut -c$MolecularChar-$MorphologChar living.data > living_morpholog.data
        let "MorphologChar -= $MolecularChar"



        #Fossil
        let "MolecularChar -= 1"
        cut -c1-$MolecularChar fossil.data > fossil_molecular.data
        let "MolecularChar += 1"
        let "MorphologChar += $MolecularChar"
        cut -c$MolecularChar-$MorphologChar fossil.data > fossil_morpholog.data
        let "MorphologChar -= $MolecularChar"




        #Outgroup
        let "MolecularChar -= 1"
        cut -c1-$MolecularChar outgroup.data > outgroup_molecular.data
        let "MolecularChar += 1"
        let "MorphologChar += $MolecularChar"
        cut -c$MolecularChar-$MorphologChar outgroup.data > outgroup_morpholog.data
        let "MorphologChar -= $MolecularChar"


    rm morpho.test
    fi



    #Building the final files to be used in the removing data part. Names are removed from the files, replaced by '@@@@@@' and then replaced by the original names at the end of the removing data part.
    sed 's/^/@@@@@@ /' living_molecular.data > living_molecular.phylip.tmp
    sed 's/^/@@@@@@ /' living_morpholog.data > living_morpholog.phylip.tmp
    sed 's/^/@@@@@@ /' fossil_molecular.data | sed '$ d' > fossil_molecular.phylip.tmp
    sed 's/^/@@@@@@ /' fossil_morpholog.data | sed '$ d' > fossil_morpholog.phylip.tmp
    sed 's/^/@@@@@@ /' outgroup_molecular.data  > outgroup_molecular.phylip.tmp
    sed 's/^/@@@@@@ /' outgroup_morpholog.data  > outgroup_morpholog.phylip.tmp
    cp fossil_morpholog.phylip.tmp fossil_morpholog.phylip
    cp living_morpholog.phylip.tmp living_morpholog.phylip



#Building a simulated matrix
else
    echo 'Building simulated matrix'



    #Combining the input values
    FossilSp=$LivingSp
    TotalSp=$LivingSp
    let "TotalSp += $FossilSp"
    TotalChar=$MolecularChar
    let "TotalChar += $MorphologChar"



    #Generating the parameters
    echo "
        library(MCMCpack) ;
        Freq<-runif(4) ;
        write(Freq/sum(Freq), file='Nfreq.rand') ;
        write(runif(1,0.1,0.5), file='Gamma.rand') ;
        write(as.vector(rdirichlet(1, c(2,2,2,2,2,2))), file='DirRates.rand', ncolumns=6) ;
        write(as.vector(rdirichlet(1, c(10,10,10,10))), file='DirFreq.rand')" | R --no-save
    Nfreq=$(sed -n '1p' Nfreq.rand | sed 's/ /,/g')
    #Gamma=$(sed -n '1p' Gamma.rand) #Random Gamma is desactivated and fixed to 0.5
    DGamma=$'0.5'
    MGamma=$'0.5'
    DirRates=$(sed -n '1p' DirRates.rand | sed 's/ /,/g')
    DirFreq=$(sed -n '1p' DirFreq.rand | sed 's/ /,/g')
    dat=$(date)
    


    #Saving the parameters
    echo "##########################" > Simulation_@.log
    echo "Matrices generated using Matsim_v4.4" >> Simulation_@.log
    echo "$dat" >> Simulation_@.log
    echo "##########################" >> Simulation_@.log

    echo "PARAMETERS:" >> Simulation_@.log
    echo "Number of living species = $LivingSp" >> Simulation_@.log
    echo "Number of fossil species = $FossilSp" >> Simulation_@.log
    echo "Number of molecular characters = $MolecularChar" >> Simulation_@.log
    echo "Number of morphological characters = $MorphologChar" >> Simulation_@.log



    #Generate the original (full) molecular and morphological matrices through R
    echo "
        #R code
        require(ape)
        require(diversitree)
        require(phyclust)

        #R# Loading the R code
        ##########################
        #Builds a Birth-Death tree with the same number of extant and extinct species
        ##########################
        #Building a Birth-Death tree by rejection sampling process using tree.bd algorithm from diversitree package (FitzJohn, 2013) with a maximum living taxa stop and allows a fixed number of associated fossils.
        #The lambda (speciation rate) and mu (extinction rate) parameters are sampeled from a uniform distribution with lambda>mu.
        #v1.0
        #----
        #Thomas Guillerme - guillert(at)tcd.ie - 12/09/2013
        ##########################
        trbd.ex<-function(extant, extinct, error=0, rename=TRUE)
        {
            #Loading the packages
            require(ape)
            require(diversitree)

            #Algorithm functions: generating a birth death tree with the 4 following conditions
            #Setting the parameters lambda and mu
            FUN.parameters<-function(x=1)
            {
                lambda<-runif(x)
                mu<-runif(x,0,lambda)
                return(cbind(lambda, mu))
            }

            #Building the birth death tree.
            #Condition 1 is that the birth death process don't fail (output = 'phylo')
            FUN.condition1<-function(extant)
            {
                pars<<-FUN.parameters()
                trbd.tmp<<-tree.bd(pars, max.taxa= extant, include.extinct=TRUE)
                while (class(trbd.tmp) !='phylo') {
                    pars<<-FUN.parameters()
                    trbd.tmp<<-tree.bd(pars, max.taxa= extant, include.extinct=TRUE)
                }
            }

            #Condition 2 is that the birth death process generates the number of extant species given in the input
            FUN.condition2<-function(extant)
            {
                FUN.condition1(extant)
                while (length(grep('sp', trbd.tmp[[3]])) != extant) {
                    FUN.condition1(extant)
                }
            }

            #Condition 3 is that the birth death process generates at least the number of extinct species given in the input (- error)
            FUN.condition3<-function(extant,extinct,error)
            {
                FUN.condition2(extant)
                while ((length(trbd.tmp[[3]])-extant) < extinct-extinct*error) {
                    FUN.condition2(extant)
                }
            }

            #When all the conditions are encountered, check if extinct=extant (+/- error), else randomly prune extinct species until extinct=extant (+/- error)
            FUN.condition4<-function(extant,extinct,error)
            {
                FUN.condition3(extant, extinct, error)
                if (extinct-extinct*error <= (length(trbd.tmp[[3]])-extant) & (length(trbd.tmp[[3]])-extant) <= extinct+extinct*error) {
                    trbd<<-trbd.tmp
                } else {
                    trbd<<-drop.tip(trbd.tmp, c(trbd.tmp[[3]][c(sample(grep('ex', trbd.tmp[[3]]), (length(grep('ex', trbd.tmp[[3]]))-extinct)))]))
                }
            }

            #Runing the tree
            trbd.tmp<-NULL
            trbd<-NULL
            FUN.condition4(extant, extinct, error)

            #Renaming the labels (optional, default=TRUE)
            if(rename==TRUE)
            {
                for (i in 1:length(trbd[[3]])) {
                    if(nchar(trbd[[3]][i])==3) {
                        trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), '000', substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep='')
                    } else {
                        if(nchar(trbd[[3]][i])==4) {
                            trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), '00', substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep='')
                        } else {
                            if(nchar(trbd[[3]][i])==5) {
                                trbd[[3]][i]<-paste(substr(trbd[[3]][i],1,2), '0', substr(trbd[[3]][i],3,nchar(trbd[[3]][i])), sep='')
                            }
                        }
                    }
                }
            }

            #Output
            return(trbd)
        }

        #R# Creating the 'True' input tree
        tree<-trbd.ex($LivingSp, $FossilSp)
        trbd<-rtree($TotalSp)
        trbd[[1]]<-tree[[1]]
        trbd[[2]]<-tree[[3]]
        trbd[[3]]<-tree[[5]]
        trbd[[4]]<-tree[[2]]

        #R# Creating the outgroup
        trbd\$root.edge<-mean(trbd\$edge.length)
        tip<-list(edge=matrix(c(2,1),1,2),
        tip.label='sp0000',
        edge.length=mean(trbd\$edge.length)+max(node.depth.edgelength(trbd)), Nnode=1)
        class(tip)<-'phylo'
        trbd<-bind.tree(trbd,tip,position=mean(trbd\$edge.length))

        write.tree(trbd, 'True_tree.tre')
        write(pars, 'tree.param', sep=',')


        #R# Creating the molecular characters matrix
        require(phyclust) # -m=model KHY ; -f=Nucleotide frequencies ; -t=transversion/translation ratio (2) ; -a=Gamma distribution shape ; -u=Number of taxa ; -l=Number of molecular characters
        Model=\"${Model}\"
        if(Model == 'GTR') {
            MatrixDNA<-seqgen(opts='-mGTR -r$DirRates -f$DirFreq -s1 -a$DGamma -u$TotalSp -l$MolecularChar', rooted.tree=trbd)
        } else {
            MatrixDNA<-seqgen(opts='-mHKY -f$Nfreq -t2 -s1 -a$DGamma -u$TotalSp -l$MolecularChar', rooted.tree=trbd)
        }
        write(sort(MatrixDNA), file='matrix_molecular.phy')

        #R# Creating the morphological characters matrix
        require(ape)
        m<-$MorphologChar
        MatrixMorpho<-as.matrix(data.frame(matrix(nrow=length(trbd[[2]]),ncol=m), row.names=trbd[[2]]))

        MGam<-rgamma(m,shape=$MGamma)
        write(MGam, file='MorphoGammaRates.txt', ncolumns=1)

        for(i in 1:m) {
            MatrixMorpho[,i]<- rTraitDisc(trbd, model='ER', k=sample(c(2:3), 1, replace=TRUE, prob=c(0.85,0.15)), rate=MGam[i])
            while (length(levels(as.factor(MatrixMorpho[,i]))) == 1) {
                MatrixMorpho[,i] <-rTraitDisc(trbd, model='ER', k=sample(c(2:3), 1, replace=TRUE, prob=c(0.85,0.15)), rate=MGam[i])
            }
        }

        #Save the table
        write.table(as.data.frame(MatrixMorpho), 'matrix_morpholog.tmp')
        " > matgen.R.tmp



    #Running
    R --no-save < matgen.R.tmp



    #Printing the parameters
    Tpar=$(sed -n '1p' tree.param)
    rm tree.param
    echo "==========================" >> Simulation_@.log
    echo "Birth-Death parameters (lambda,mu)=$Tpar" >> Simulation_@.log
    echo "$Model" > model.tmp
    if grep "HKY" model.tmp > /dev/null
    then
        echo "Molecular data:" >> Simulation_@.log
        echo "Chosen model = $Model" >> Simulation_@.log
        echo "Nucleotide frequencies = $Nfreq" >> Simulation_@.log
        echo "Molecular rates distribution (gamma) alpha = $DGamma" >> Simulation_@.log
    else
        echo "Molecular data:" >> Simulation_@.log
        echo "Chosen model = $Model" >> Simulation_@.log
        echo "Nucleotide frequencies (Dirichlet distribution) = $DirFreq" >> Simulation_@.log
        echo "Rates (Dirichlet distribution) = $DirRates" >> Simulation_@.log
        echo "Molecular rates distribution (gamma) alpha = $DGamma" >> Simulation_@.log
    fi
    rm model.tmp



    echo "Morphological data:" >> Simulation_@.log
    echo "Chosen model = Equal Rates" >> Simulation_@.log
    echo "States frequencies = 2:0.85, 3:0.15" >> Simulation_@.log
    echo "Morphological rates distribution (gamma) alpha = $MGamma" >> Simulation_@.log
    echo "All rates are saved in the MorphoGammaRates.txt file" >> Simulation_@.log
    echo "##########################" >> Simulation_@.log

    echo "TIMER:" >> Simulation_@.log
    echo "Matrices generations start:" >> Simulation_@.log
    date >> Simulation_@.log



    #Adding the outgroup
    let "TotalSp += 1"



    #Merge the two molecular and morphological matrices
    #Modifying the morphological matrix (transforming the characters from 1->0 sorted and cleaned)
    tail -n +2 matrix_morpholog.tmp | sort -r | sed 's/"......"//g' | sed 's/[[:space:]]//g' | sed 's/1/0/g' | sed 's/2/1/g' | sed 's/3/2/g' | sed 's/4/3/g' | sed 's/5/4/g' | sed 's/6/5/g' | sed 's/7/6/g' | sed 's/8/7/g' | sed 's/9/8/g' > matrix_morpholog.phy.tmp
    echo "head" > head.tmp
    cat head.tmp matrix_morpholog.phy.tmp > matrix_morpholog.phy

    #Merging the matrices
    paste -d\\0 matrix_molecular.phy matrix_morpholog.phy > matrix_combined.tmp
    sed '1d' matrix_combined.tmp > matrix_combined.phy.tmp
    echo "$TotalSp $TotalChar" > matrix_combined.phylip  ; cat matrix_combined.phy.tmp >> matrix_combined.phylip



    #Splitting the complete matrix



    #Selecting the fossils
    let "FossilSp += 1"
    sed -n '2,'"${FossilSp}"'p' matrix_combined.phylip > fossil_combined.phy
    let "FossilSp -= 1"



    #Selecting the outgroup
    let "FossilSp += 2"
    sed -n ''"${FossilSp}"'p' matrix_combined.phylip > outgroup_combined.phy
    let "FossilSp -= 2"



    #Selecting the livings
    let "LivingSp += 3"
    let "TotalSp += 1"
    sed -n ''"${LivingSp}"','"${TotalSp}"'p' matrix_combined.phylip > living_combined.phy
    let "LivingSp -= 3"
    let "TotalSp -= 1"



    #Outgroup
    #DNA
    let "MolecularChar += 10"
    cut -c1-7,8-$MolecularChar outgroup_combined.phy > outgroup_molecular.phylip.tmp
    #Morphological
    let "MolecularChar += 1"
    let "MorphologChar += $MolecularChar"
    cut -c1-7,$MolecularChar-$MorphologChar outgroup_combined.phy > outgroup_morpholog.phylip.tmp
    let "MorphologChar -= $MolecularChar"
    let "MolecularChar -= 11"



    #Living species
    #DNA
    let "MolecularChar += 10"
    cut -c1-7,8-$MolecularChar living_combined.phy > living_molecular.phylip.tmp
    #Morphological
    let "MolecularChar += 1"
    let "MorphologChar += $MolecularChar"
    cut -c1-7,$MolecularChar-$MorphologChar living_combined.phy > living_morpholog.phylip.tmp
    let "MorphologChar -= $MolecularChar"
    let "MolecularChar -= 11"



    #Fossil species
    #DNA (including missing data '?')
    let "MolecularChar += 10"
    cut -c1-7,8-$MolecularChar fossil_combined.phy > fossil_molecular.phyl.tmp
    sed 's/[ACGT]/?/g' fossil_molecular.phyl.tmp > fossil_molecular.phylip.tmp
    #Morphological
    let "MolecularChar += 1"
    let "MorphologChar += $MolecularChar"
    cut -c1-7,$MolecularChar-$MorphologChar fossil_combined.phy > fossil_morpholog.phylip.tmp
    let "MorphologChar -= $MolecularChar"
    let "MolecularChar -= 11"



    #Renaming
    cp fossil_morpholog.phylip.tmp fossil_morpholog.phylip
    cp living_morpholog.phylip.tmp living_morpholog.phylip
fi


#Step 2 - INCLUDING GAPS IN THE COMPLETE MATRIX



#Building the 5 NL parameters on the living_morpholog matrix (NL= number of coded living species for morphological characters)



#NL00 (remove nothing)
cp living_morpholog.phylip living_morpholog_NL00.phy
LivSp=$LivingSp



#NL10 (remove 10% of the species)
LivingSpNL=$LivingSp



#Preparing the temporary removing data script (NL*.sh)
echo "LivingSpNL=\$'$LivingSpNL'" > NL10.s

let 'LivingSpNL *=10'
let 'LivingSpNL /=100'

echo "
    let 'LivingSpNL *=10'
    let 'LivingSpNL /=100'

    #Condition: at least one species must be selected
    if ((LivingSpNL == 0))
        then
            let 'LivingSpNL +=1'
        else
            echo "nothing" > /dev/null
    fi



    #Generating a list of random numbers without replacement in R
    echo 'write(sample(1:$LivSp, $LivingSpNL), \"NL10.rand\", ncolumns=1)' | R --no-save



    #Extracting the random numbers from R and saving them as a shell argument
    for n in \$(seq 1 \$LivingSpNL)
    do
        echo \"NL10\"
        let \"NL10_\${n}=\$(sed -n ''\"\${n}\"'p' NL10.rand)\"
        echo \"NL10_\${n}=\\$'\$NL10_\${n}'\"
    done" >> NL10.s



#Passing the arguments into 'sed' and replace the randomly selected character scores by '?\?\?\?\?'
for n in $(seq 1 $LivingSpNL)
do
    echo "sed ''\"\$NL10_${n}\"' s/[0-9][0-9][0-9][0-9][0-9]/\?\?\?\?\?/g' |" >> NL10.t
done
sed '1 s/|/ living_morpholog.phylip |/' NL10.t | sed '$ s/|/ > living_morpholog_NL10.phy/' > NL10.tmp




if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    #Modify the temporary removing data script by replacing 5 nucleotides/? by only 1
    echo 'Using the input matrix' > /dev/null
    cat NL10.s NL10.tmp | sed 's/\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]/\[0-9\]/g' | sed 's/\\?\\?\\?\\?\\?/\\?/g' > NL10.sh
    sh NL10.sh
    sed 's/-/\?/g' living_morpholog_NL10.phy > living_morpholog_NL10.phy.tmp
    rm living_morpholog_NL10.phy
    cp living_morpholog_NL10.phy.tmp living_morpholog_NL10.phy
    rm living_morpholog_NL10.phy.tmp
else
    echo 'Simulated matrix' > /dev/null
    cat NL10.s NL10.tmp > NL10.sh
    sh NL10.sh
fi



rm NL10.t ; rm NL10.s ; rm NL10.tmp ; rm NL10.sh



#NL25 (remove 25% of the species) (Script follow the indications above for the temporary removing data script (NL*.sh))
LivingSpNL=$LivingSp

echo "LivingSpNL=\$'$LivingSpNL'" > NL25.s

let 'LivingSpNL *=25'
let 'LivingSpNL /=100'

echo "
    let 'LivingSpNL *=25'
    let 'LivingSpNL /=100'

    #Condition: at least one species must be selected
    if ((LivingSpNL == 0))
        then
            let 'LivingSpNL +=1'
        else
            echo "nothing" > /dev/null
    fi

    #Generating a list of random numbers without replacement in R
    echo 'write(sample(1:$LivSp, $LivingSpNL), \"NL25.rand\", ncolumns=1)' | R --no-save

    #Extracting the random numbers from R and saving them as a shell argument
    for n in \$(seq 1 \$LivingSpNL)
    do
        echo \"NL25\"
        let \"NL25_\${n}=\$(sed -n ''\"\${n}\"'p' NL25.rand)\"
        echo \"NL25_\${n}=\\$'\$NL25_\${n}'\"
    done" >> NL25.s

#Passing the arguments into 'sed' and replace the randomly selected character scores by '?\?\?\?\?'
for n in $(seq 1 $LivingSpNL)
do
    echo "sed ''\"\$NL25_${n}\"' s/[0-9][0-9][0-9][0-9][0-9]/\?\?\?\?\?/g' |" >> NL25.t
done
sed '1 s/|/ living_morpholog.phylip |/' NL25.t | sed '$ s/|/ > living_morpholog_NL25.phy/' > NL25.tmp



if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    #Modify the temporary removing data script by replacing 5 nucleotides/? by only 1
    echo 'Using the input matrix' > /dev/null
    cat NL25.s NL25.tmp | sed 's/\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]/\[0-9\]/g' | sed 's/\\?\\?\\?\\?\\?/\\?/g' > NL25.sh
    sh NL25.sh
    sed 's/-/\?/g' living_morpholog_NL25.phy > living_morpholog_NL25.phy.tmp
    rm living_morpholog_NL25.phy
    cp living_morpholog_NL25.phy.tmp living_morpholog_NL25.phy
    rm living_morpholog_NL25.phy.tmp
else
    echo 'Simulated matrix' > /dev/null
    cat NL25.s NL25.tmp > NL25.sh
    sh NL25.sh
fi

rm NL25.t ; rm NL25.s ; rm NL25.tmp ; rm NL25.sh



##NL50 (remove 50% of the species)
LivingSpNL=$LivingSp

echo "LivingSpNL=\$'$LivingSpNL'" > NL50.s

let 'LivingSpNL *=50'
let 'LivingSpNL /=100'

echo "
    let 'LivingSpNL *=50'
    let 'LivingSpNL /=100'

    #Condition: at least one species must be selected
    if ((LivingSpNL == 0))
        then
            let 'LivingSpNL +=1'
        else
            echo "nothing" > /dev/null
    fi

    #Generating a list of random numbers without replacement in R
    echo 'write(sample(1:$LivSp, $LivingSpNL), \"NL50.rand\", ncolumns=1)' | R --no-save

    #Extracting the random numbers from R and saving them as a shell argument
    for n in \$(seq 1 \$LivingSpNL)
    do
        echo \"NL50\"
        let \"NL50_\${n}=\$(sed -n ''\"\${n}\"'p' NL50.rand)\"
        echo \"NL50_\${n}=\\$'\$NL50_\${n}'\"
    done" >> NL50.s

#Passing the arguments into 'sed' and replace the randomly selected character scores by '?\?\?\?\?'
for n in $(seq 1 $LivingSpNL)
do
    echo "sed ''\"\$NL50_${n}\"' s/[0-9][0-9][0-9][0-9][0-9]/\?\?\?\?\?/g' |" >> NL50.t
done
sed '1 s/|/ living_morpholog.phylip |/' NL50.t | sed '$ s/|/ > living_morpholog_NL50.phy/' > NL50.tmp



if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    #Modify the temporary removing data script by replacing 5 nucleotides/? by only 1
    echo 'Using the input matrix' > /dev/null
    cat NL50.s NL50.tmp | sed 's/\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]/\[0-9\]/g' | sed 's/\\?\\?\\?\\?\\?/\\?/g' > NL50.sh
    sh NL50.sh
    sed 's/-/\?/g' living_morpholog_NL50.phy > living_morpholog_NL50.phy.tmp
    rm living_morpholog_NL50.phy
    cp living_morpholog_NL50.phy.tmp living_morpholog_NL50.phy
    rm living_morpholog_NL50.phy.tmp
else
    echo 'Simulated matrix' > /dev/null
    cat NL50.s NL50.tmp > NL50.sh
    sh NL50.sh
fi

rm NL50.t ; rm NL50.s ; rm NL50.tmp ; rm NL50.sh



##NL75 (remove 75% of the species)
LivingSpNL=$LivingSp

echo "LivingSpNL=\$'$LivingSpNL'" > NL75.s

let 'LivingSpNL *=75'
let 'LivingSpNL /=100'

echo "
    let 'LivingSpNL *=75'
    let 'LivingSpNL /=100'

    if ((LivingSpNL == 0))
        then
            let 'LivingSpNL +=1'
        else
            echo "nothing" > /dev/null
    fi

    #Generating a list of random numbers without replacement in R
    echo 'write(sample(1:$LivSp, $LivingSpNL), \"NL75.rand\", ncolumns=1)' | R --no-save

    #Extracting the random numbers from R and saving them as a shell argument
    for n in \$(seq 1 \$LivingSpNL)
    do
        echo \"NL75\"
        let \"NL75_\${n}=\$(sed -n ''\"\${n}\"'p' NL75.rand)\"
        echo \"NL75_\${n}=\\$'\$NL75_\${n}'\"
    done" >> NL75.s

#Passing the arguments into 'sed' and replace the randomly selected character scores by '?\?\?\?\?'
for n in $(seq 1 $LivingSpNL)
do
    echo "sed ''\"\$NL75_${n}\"' s/[0-9][0-9][0-9][0-9][0-9]/\?\?\?\?\?/g' |" >> NL75.t
done
sed '1 s/|/ living_morpholog.phylip |/' NL75.t | sed '$ s/|/ > living_morpholog_NL75.phy/' > NL75.tmp




if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    #Modify the temporary removing data script by replacing 5 nucleotides/? by only 1
    echo 'Using the input matrix' > /dev/null
    cat NL75.s NL75.tmp | sed 's/\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]/\[0-9\]/g' | sed 's/\\?\\?\\?\\?\\?/\\?/g' > NL75.sh
    sh NL75.sh
    sed 's/-/\?/g' living_morpholog_NL75.phy > living_morpholog_NL75.phy.tmp
    rm living_morpholog_NL75.phy
    cp living_morpholog_NL75.phy.tmp living_morpholog_NL75.phy
    rm living_morpholog_NL75.phy.tmp
else
    echo 'Simulated matrix' > /dev/null
    cat NL75.s NL75.tmp > NL75.sh
    sh NL75.sh
fi

rm NL75.t ; rm NL75.s ; rm NL75.tmp ; rm NL75.sh



#Building the 5 NF parameters on the fossil_morpholog matrix (NF= number of missing data in the fossil matrix)



#Preparing the R code
echo "
    #R#Phylip file input
    Inp<-read.delim('fossil_morpholog.phylip',header=F)
    
    #Rransforming the phylip file in a matrix
    m<-matrix(NA,nrow(Inp),length(noquote(strsplit(as.character(Inp[1,]), NULL)[[1]])))
    for (i in 1:nrow(Inp)){
        m[i,]<-noquote(strsplit(as.character(Inp[i,]), NULL)[[1]])
    }
    
    #Introducing missing data
    M<-as.vector(m[,-c(1:7)])
    M[sample(1:length(M), @@@)]<-'?'
    Ma<-matrix(M,ncol=ncol(m[,-c(1:7)]),nrow=nrow(m))
    MA<-noquote(Ma)
        
    #Transforming the matrix into a phylip file
    sp<-rep(NA, nrow(m))
    for (i in 1:nrow(m)){
        sp[i]<-paste(m[i,1],m[i,2],m[i,3],m[i,4],m[i,5],m[i,6], sep='')
    }
    cha<-rep(NA, nrow(m))
    for (i in 1:nrow(m)){
        cha[i]<-paste(as.vector(Ma[i,]),sep='',collapse='')
    }
    
    #Phylip output
    write.table(data.frame(X25_100=cha,row.names=sp), file='fossil_morpholog_NF@.phy.tmp')" > NF@.R.tmp



#NF00 (remove nothing)
cp fossil_morpholog.phylip fossil_morpholog_NF00.phy



#NF10 (remove 10% of the cells)
#Preparing the R code
FossilChar=$FossilSp
let "FossilChar *=$MorphologChar"
let "FossilChar *=10"
let "FossilChar /=100"
sed 's/@@@/'"$FossilChar"'/g' NF@.R.tmp | sed 's/NF@/NF10/g' > NF10.R.tmp



if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    #Modifying the row names to match with the input
    echo 'Using the input matrix' > /dev/null
    sed 's/,row.names=sp/,row.names=sprintf\("%06d",1:'"$FossilSp"')/g' NF10.R.tmp > NF10.R.tmpBIS
    rm NF10.R.tmp
    cp NF10.R.tmpBIS NF10.R.tmp
    rm NF10.R.tmpBIS
else 
    echo 'Simulated matrix' > /dev/null
fi



#Running the R code
R --no-save < NF10.R.tmp




#Correcting the file
sed 's/"//g' fossil_morpholog_NF10.phy.tmp | sed '1d' > fossil_morpholog_NF10.phy




#Accept/Reject algorithm: every fossil should have a least 5% data
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MorphologChar
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c8-$MorphologChar fossil_morpholog_NF10.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
        
    while grep "[0-9]" NF_rejection.test > /dev/null
    do
        R --no-save < NF10.R.tmp
        sed 's/"//g' fossil_morpholog_NF10.phy.tmp | sed '1d' > fossil_morpholog_NF10.phy
        cut -c8-$MorphologChar fossil_morpholog_NF10.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
    done
    rm NF_rejection.test

else
    echo "nothing" > /dev/null
fi



#NF25
#Preparing the R code
FossilChar=$FossilSp
let "FossilChar *=$MorphologChar"
let "FossilChar *=25"
let "FossilChar /=100"
sed 's/@@@/'"$FossilChar"'/g' NF@.R.tmp | sed 's/NF@/NF25/g' > NF25.R.tmp

if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    echo 'Using the input matrix' > /dev/null
    sed 's/,row.names=sp/,row.names=sprintf\("%06d",1:'"$FossilSp"')/g' NF25.R.tmp > NF25.R.tmpBIS
    rm NF25.R.tmp
    cp NF25.R.tmpBIS NF25.R.tmp
    rm NF25.R.tmpBIS
else
    echo 'Simulated matrix' > /dev/null
fi

#Runing
R --no-save < NF25.R.tmp

#Correcting the file
sed 's/"//g' fossil_morpholog_NF25.phy.tmp | sed '1d' > fossil_morpholog_NF25.phy

#Accept/Rejet
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MorphologChar
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c8-$MorphologChar fossil_morpholog_NF25.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test

    while grep "[0-9]" NF_rejection.test > /dev/null
    do
        R --no-save < NF25.R.tmp
        sed 's/"//g' fossil_morpholog_NF25.phy.tmp | sed '1d' > fossil_morpholog_NF25.phy
        cut -c8-$MorphologChar fossil_morpholog_NF25.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
    done
    rm NF_rejection.test
else
    echo "nothing" > /dev/null
fi



#NF50
#Preparing the R code
FossilChar=$FossilSp
let "FossilChar *=$MorphologChar"
let "FossilChar *=50"
let "FossilChar /=100"
sed 's/@@@/'"$FossilChar"'/g' NF@.R.tmp | sed 's/NF@/NF50/g' > NF50.R.tmp

if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    echo 'Using the input matrix' > /dev/null
    sed 's/,row.names=sp/,row.names=sprintf\("%06d",1:'"$FossilSp"')/g' NF50.R.tmp > NF50.R.tmpBIS
    rm NF50.R.tmp
    cp NF50.R.tmpBIS NF50.R.tmp
    rm NF50.R.tmpBIS
else
    echo 'Simulated matrix' > /dev/null
fi

#Runing
R --no-save < NF50.R.tmp

#Correcting the file
sed 's/"//g' fossil_morpholog_NF50.phy.tmp | sed '1d' > fossil_morpholog_NF50.phy

#Accept/Rejet
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MorphologChar
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c8-$MorphologChar fossil_morpholog_NF50.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test

    while grep "[0-9]" NF_rejection.test > /dev/null
    do
        R --no-save < NF50.R.tmp
        sed 's/"//g' fossil_morpholog_NF50.phy.tmp | sed '1d' > fossil_morpholog_NF50.phy
        cut -c8-$MorphologChar fossil_morpholog_NF50.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
    done
    rm NF_rejection.test
else
    echo "nothing" > /dev/null
fi



#NF75
#Preparing the R code
FossilChar=$FossilSp
let "FossilChar *=$MorphologChar"
let "FossilChar *=75"
let "FossilChar /=100"
sed 's/@@@/'"$FossilChar"'/g' NF@.R.tmp | sed 's/NF@/NF75/g' > NF75.R.tmp

if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then
    echo 'Using the input matrix' > /dev/null
    sed 's/,row.names=sp/,row.names=sprintf\("%06d",1:'"$FossilSp"')/g' NF75.R.tmp > NF75.R.tmpBIS
    rm NF75.R.tmp
    cp NF75.R.tmpBIS NF75.R.tmp
    rm NF75.R.tmpBIS
else
   echo 'Simulated matrix' > /dev/null
fi

#Runing
R --no-save < NF75.R.tmp

#Correcting the file
sed 's/"//g' fossil_morpholog_NF75.phy.tmp | sed '1d' > fossil_morpholog_NF75.phy

#Accept/Rejet
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MorphologChar
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c8-$MorphologChar fossil_morpholog_NF75.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
    
    while grep "[0-9]" NF_rejection.test > /dev/null
    do
        R --no-save < NF75.R.tmp
        sed 's/"//g' fossil_morpholog_NF75.phy.tmp | sed '1d' > fossil_morpholog_NF75.phy
        cut -c8-$MorphologChar fossil_morpholog_NF75.phy | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NF_rejection.test
    done
    rm NF_rejection.test
else
    echo "nothing" > /dev/null
fi



#Pairwise combination of the sub NL (step 2.2) and NF (step 2.3) matrices



#Combine outgroup + NL00 + NF00/10/25/50/75
echo "$TotalSp $MorphologChar" > MM_L00F00.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L00F00.phyl ; cat living_morpholog_NL00.phy >> MM_L00F00.phyl ; cat fossil_morpholog_NF00.phy >> MM_L00F00.phyl
echo "$TotalSp $MorphologChar" > MM_L00F10.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L00F10.phyl ; cat living_morpholog_NL00.phy >> MM_L00F10.phyl ; cat fossil_morpholog_NF10.phy >> MM_L00F10.phyl
echo "$TotalSp $MorphologChar" > MM_L00F25.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L00F25.phyl ; cat living_morpholog_NL00.phy >> MM_L00F25.phyl ; cat fossil_morpholog_NF25.phy >> MM_L00F25.phyl
echo "$TotalSp $MorphologChar" > MM_L00F50.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L00F50.phyl ; cat living_morpholog_NL00.phy >> MM_L00F50.phyl ; cat fossil_morpholog_NF50.phy >> MM_L00F50.phyl
echo "$TotalSp $MorphologChar" > MM_L00F75.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L00F75.phyl ; cat living_morpholog_NL00.phy >> MM_L00F75.phyl ; cat fossil_morpholog_NF75.phy >> MM_L00F75.phyl



#Combine outgroup + NL10 + NF00/10/25/50/75
echo "$TotalSp $MorphologChar" > MM_L10F00.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L10F00.phyl ; cat living_morpholog_NL10.phy >> MM_L10F00.phyl ; cat fossil_morpholog_NF00.phy >> MM_L10F00.phyl
echo "$TotalSp $MorphologChar" > MM_L10F10.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L10F10.phyl ; cat living_morpholog_NL10.phy >> MM_L10F10.phyl ; cat fossil_morpholog_NF10.phy >> MM_L10F10.phyl
echo "$TotalSp $MorphologChar" > MM_L10F25.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L10F25.phyl ; cat living_morpholog_NL10.phy >> MM_L10F25.phyl ; cat fossil_morpholog_NF25.phy >> MM_L10F25.phyl
echo "$TotalSp $MorphologChar" > MM_L10F50.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L10F50.phyl ; cat living_morpholog_NL10.phy >> MM_L10F50.phyl ; cat fossil_morpholog_NF50.phy >> MM_L10F50.phyl
echo "$TotalSp $MorphologChar" > MM_L10F75.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L10F75.phyl ; cat living_morpholog_NL10.phy >> MM_L10F75.phyl ; cat fossil_morpholog_NF75.phy >> MM_L10F75.phyl



#Combine outgroup + NL25 + NF00/10/25/50/75
echo "$TotalSp $MorphologChar" > MM_L25F00.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L25F00.phyl ; cat living_morpholog_NL25.phy >> MM_L25F00.phyl ; cat fossil_morpholog_NF00.phy >> MM_L25F00.phyl
echo "$TotalSp $MorphologChar" > MM_L25F10.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L25F10.phyl ; cat living_morpholog_NL25.phy >> MM_L25F10.phyl ; cat fossil_morpholog_NF10.phy >> MM_L25F10.phyl
echo "$TotalSp $MorphologChar" > MM_L25F25.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L25F25.phyl ; cat living_morpholog_NL25.phy >> MM_L25F25.phyl ; cat fossil_morpholog_NF25.phy >> MM_L25F25.phyl
echo "$TotalSp $MorphologChar" > MM_L25F50.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L25F50.phyl ; cat living_morpholog_NL25.phy >> MM_L25F50.phyl ; cat fossil_morpholog_NF50.phy >> MM_L25F50.phyl
echo "$TotalSp $MorphologChar" > MM_L25F75.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L25F75.phyl ; cat living_morpholog_NL25.phy >> MM_L25F75.phyl ; cat fossil_morpholog_NF75.phy >> MM_L25F75.phyl



#Combine outgroup + NL50 + NF00/10/25/50/75
echo "$TotalSp $MorphologChar" > MM_L50F00.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L50F00.phyl ; cat living_morpholog_NL50.phy >> MM_L50F00.phyl ; cat fossil_morpholog_NF00.phy >> MM_L50F00.phyl
echo "$TotalSp $MorphologChar" > MM_L50F10.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L50F10.phyl ; cat living_morpholog_NL50.phy >> MM_L50F10.phyl ; cat fossil_morpholog_NF10.phy >> MM_L50F10.phyl
echo "$TotalSp $MorphologChar" > MM_L50F25.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L50F25.phyl ; cat living_morpholog_NL50.phy >> MM_L50F25.phyl ; cat fossil_morpholog_NF25.phy >> MM_L50F25.phyl
echo "$TotalSp $MorphologChar" > MM_L50F50.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L50F50.phyl ; cat living_morpholog_NL50.phy >> MM_L50F50.phyl ; cat fossil_morpholog_NF50.phy >> MM_L50F50.phyl
echo "$TotalSp $MorphologChar" > MM_L50F75.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L50F75.phyl ; cat living_morpholog_NL50.phy >> MM_L50F75.phyl ; cat fossil_morpholog_NF75.phy >> MM_L50F75.phyl



#Combine outgroup + NL75 + NF00/10/25/50/75
echo "$TotalSp $MorphologChar" > MM_L75F00.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L75F00.phyl ; cat living_morpholog_NL75.phy >> MM_L75F00.phyl ; cat fossil_morpholog_NF00.phy >> MM_L75F00.phyl
echo "$TotalSp $MorphologChar" > MM_L75F10.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L75F10.phyl ; cat living_morpholog_NL75.phy >> MM_L75F10.phyl ; cat fossil_morpholog_NF10.phy >> MM_L75F10.phyl
echo "$TotalSp $MorphologChar" > MM_L75F25.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L75F25.phyl ; cat living_morpholog_NL75.phy >> MM_L75F25.phyl ; cat fossil_morpholog_NF25.phy >> MM_L75F25.phyl
echo "$TotalSp $MorphologChar" > MM_L75F50.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L75F50.phyl ; cat living_morpholog_NL75.phy >> MM_L75F50.phyl ; cat fossil_morpholog_NF50.phy >> MM_L75F50.phyl
echo "$TotalSp $MorphologChar" > MM_L75F75.phyl ; cat outgroup_morpholog.phylip.tmp >> MM_L75F75.phyl ; cat living_morpholog_NL75.phy >> MM_L75F75.phyl ; cat fossil_morpholog_NF75.phy >> MM_L75F75.phyl



#Including the NC parameters to the combined NL+NF matrices (NC= number of overall morphological characters (number of columns in the matrix))



MorphoChar=$MorphologChar
let "MorphoChar += 7"



#NC00 = 10% of morphological characters relative to molecular ones (keep all the characters)
for f in *.phyl
do
    prefix=$(basename $f .phyl)
    echo ${prefix}C00
    cp $f ${prefix}C00.phyli
done



#NC10
MChar10=$MorphologChar
let "MChar10 *=90"
let "MChar10 /=100"



#Select the columns to remove
echo "write(sample(8:$MorphoChar, $MChar10), 'NC10.rand', ncolumns=$MorphologChar)" | R --no-save
NC10=$(sed -n 'p' NC10.rand | sed 's/ /,/g')



#Remove the columns with accept/reject
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MChar10
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c1-7,$NC10 MM_L00F75.phyl > C10_rejection.test
    let "MChar10 +=8"
    let "LivingSp +=3"
    cut -c8-$MChar10 C10_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
    let "LivingSp -=3"
    let "MChar10 -=8"

    while grep "[0-9]" NC_rejection.test > /dev/null
    do
        echo "write(sample(8:$MorphoChar, $MChar10), 'NC10.rand', ncolumns=$MorphologChar)" | R --no-save
        NC10=$(sed -n 'p' NC10.rand | sed 's/ /,/g')
        cut -c1-7,$NC10 MM_L00F75.phyl > C10_rejection.test
        let "MChar10 +=8"
        let "LivingSp +=3"
        cut -c8-$MChar10 C10_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
        let "LivingSp -=3"
        let "MChar10 -=8"
    done

    rm NC_rejection.test
    rm C10_rejection.test

else
    echo "nothing" > /dev/null
fi


#Apply to all
for f in *.phyl
do
    prefix=$(basename $f .phyl)
    echo ${prefix}C10
    cut -c1-7,$NC10 $f > ${prefix}C10.phyli
done



#NC25
MChar25=$MorphologChar
let "MChar25 *=75"
let "MChar25 /=100"



#Select the columns to remove
echo "write(sample(8:$MorphoChar, $MChar25), 'NC25.rand', ncolumns=$MorphologChar)" | R --no-save
NC25=$(sed -n 'p' NC25.rand | sed 's/ /,/g')



#Remove the columns with accept/reject
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MChar25
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c1-7,$NC25 MM_L00F75.phyl > C25_rejection.test
    let "MChar25 +=8"
    let "LivingSp +=3"
    cut -c8-$MChar25 C25_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
    let "LivingSp -=3"
    let "MChar25 -=8"

    while grep "[0-9]" NC_rejection.test > /dev/null
    do
        echo "write(sample(8:$MorphoChar, $MChar25), 'NC25.rand', ncolumns=$MorphologChar)" | R --no-save
        NC25=$(sed -n 'p' NC25.rand | sed 's/ /,/g')
        cut -c1-7,$NC25 MM_L00F75.phyl > C25_rejection.test
        let "MChar25 +=8"
        let "LivingSp +=3"
        cut -c8-$MChar25 C25_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
        let "LivingSp -=3"
        let "MChar25 -=8"
    done

    rm NC_rejection.test
    rm C25_rejection.test

else
    echo "nothing" > /dev/null
fi

##Apply to all
for f in *.phyl
do
    prefix=$(basename $f .phyl)
    echo ${prefix}C25
    cut -c1-7,$NC25 $f > ${prefix}C25.phyli
done



#NC50
MChar50=$MorphologChar
let "MChar50 *=50"
let "MChar50 /=100"



#Select the columns to remove
echo "write(sample(8:$MorphoChar, $MChar50), 'NC50.rand', ncolumns=$MorphologChar)" | R --no-save
NC50=$(sed -n 'p' NC50.rand | sed 's/ /,/g')




#Remove the columns with accept/reject
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MChar50
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c1-7,$NC50 MM_L00F75.phyl > C50_rejection.test
    let "MChar50 +=8"
    let "LivingSp +=3"
    cut -c8-$MChar50 C50_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
    let "LivingSp -=3"
    let "MChar50 -=8"

    while grep "[0-9]" NC_rejection.test > /dev/null
    do
        echo "write(sample(8:$MorphoChar, $MChar50), 'NC50.rand', ncolumns=$MorphologChar)" | R --no-save
        NC50=$(sed -n 'p' NC50.rand | sed 's/ /,/g')
        cut -c1-7,$NC50 MM_L00F75.phyl > C50_rejection.test
        let "MChar50 +=8"
        let "LivingSp +=3"
        cut -c8-$MChar50 C50_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
        let "LivingSp -=3"
        let "MChar50 -=8"
    done

    rm NC_rejection.test
    rm C50_rejection.test

else
    echo "nothing" > /dev/null
fi

#Apply to all
for f in *.phyl
do
    prefix=$(basename $f .phyl)
    echo ${prefix}C50
    cut -c1-7,$NC50 $f > ${prefix}C50.phyli
done



#NC75
MChar75=$MorphologChar
let "MChar75 *=25"
let "MChar75 /=100"



#Select the columns to remove
echo "write(sample(8:$MorphoChar, $MChar75), 'NC75.rand', ncolumns=$MorphologChar)" | R --no-save
NC75=$(sed -n 'p' NC75.rand | sed 's/ /,/g')




#Remove the columns with accept/reject
if grep "Matrices generated using "  Simulation_@.log > /dev/null
then
    Rejection=$MChar75
    let "Rejection *=5"
    let "Rejection /=100"
    cut -c1-7,$NC75 MM_L00F75.phyl > C75_rejection.test
    let "MChar75 +=8"
    let "LivingSp +=3"
    cut -c8-$MChar75 C75_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
    let "LivingSp -=3"
    let "MChar75 -=8"

    while grep "[0-9]" NC_rejection.test > /dev/null
    do
        echo "write(sample(8:$MorphoChar, $MChar75), 'NC75.rand', ncolumns=$MorphologChar)" | R --no-save
        NC75=$(sed -n 'p' NC75.rand | sed 's/ /,/g')
        cut -c1-7,$NC75 MM_L00F75.phyl > C75_rejection.test
        let "MChar75 +=8"
        let "LivingSp +=3"
        cut -c8-$MChar75 C75_rejection.test | sed -n ''"$LivingSp"',$p' | awk -F "[0-9]" '{print NF-1}' | awk '$1 < '"$Rejection"'' > NC_rejection.test
        let "LivingSp -=3"
        let "MChar75 -=8"
    done

    rm NC_rejection.test
    rm C75_rejection.test

else
    echo "nothing" > /dev/null
fi

#Apply to all
for f in *.phyl
do
    prefix=$(basename $f .phyl)
    echo ${prefix}C75
    cut -c1-7,$NC75 $f > ${prefix}C75.phyli
done



#Step 3 - CREATING THE 125 FINAL MATRICES



#Combining the Morphological and the Molecular matrices



#Combining the molecular data
echo "$TotalSp $MolecularChar" > MD.phy ; cat outgroup_molecular.phylip.tmp >> MD.phy ; cat living_molecular.phylip.tmp >> MD.phy ; cat fossil_molecular.phylip.tmp >> MD.phy



#Combining DNA data to the morphological ones
for f in *.phyli
do
    prefix=$(basename $f .phyli)
    echo loading ${prefix}
    paste MD.phy $f > ${prefix}.phylip_tmp
done

for f in *.phylip_tmp
do
    prefix=$(basename $f .phylip_tmp)
    echo combining ${prefix}
    sed 's/[[:space:]]......[[:space:]]//g' $f > ${prefix}.phylip.tmp
done



#Removing the files not to convert
rm fossil_molecular.phylip.tmp
rm fossil_morpholog.phylip.tmp
rm living_molecular.phylip.tmp
rm living_morpholog.phylip.tmp
rm outgroup_molecular.phylip.tmp
rm outgroup_morpholog.phylip.tmp



#Removing extra characters
for f in *.phylip.tmp
do
    prefix=$(basename $f .phylip.tmp)
    echo saving ${prefix}

    if grep "Input matrix: $Input"  Simulation_@.log > /dev/null

    #Renaming the taxa and adding a dummy morphological character (sed 's/$/0/g')
    then
        paste input.name $f | sed 's/@@@@@@ //g' | sed 's/[[:space:]][0-9][0-9][[:space:]]//g' | sed 's/[[:space:]]??[[:space:]]//g' | sed '1d' | sed '$d' | sed 's/$/0/g' > ${prefix}.phylip2.tmp
    else
        sed 's/[[:space:]][0-9][0-9][[:space:]]//g' $f | sed 's/[[:space:]]??[[:space:]]//g' | sed '1d' > ${prefix}.phylip2.tmp
    fi
done




#Proper Phylip format
if grep "Input matrix: $Input"  Simulation_@.log > /dev/null
then #Add the dummy character
    let "TotalChar +=1"
    echo "$TotalSp $TotalChar" > headC00.tmp
    let "MChar10 += $MolecularChar"
    let "Mchar10 +=1"
    echo "$TotalSp $MChar10" > headC10.tmp
    let "MChar25 += $MolecularChar"
    let "Mchar25 +=1"
    echo "$TotalSp $MChar25" > headC25.tmp
    let "MChar50 += $MolecularChar"
    let "Mchar50 +=1"
    echo "$TotalSp $MChar50" > headC50.tmp
    let "MChar75 += $MolecularChar"
    let "Mchar75 +=1"
    echo "$TotalSp $MChar75" > headC75.tmp    
else
    echo "$TotalSp $TotalChar" > headC00.tmp
    let "MChar10 += $MolecularChar"
    echo "$TotalSp $MChar10" > headC10.tmp
    let "MChar25 += $MolecularChar"
    echo "$TotalSp $MChar25" > headC25.tmp
    let "MChar50 += $MolecularChar"
    echo "$TotalSp $MChar50" > headC50.tmp
    let "MChar75 += $MolecularChar"
    echo "$TotalSp $MChar75" > headC75.tmp
fi

for f in *C00.phylip2.tmp
do
    prefix=$(basename $f .phylip2.tmp)
    cat headC00.tmp $f > ${prefix}.phylip
done
for f in *C10.phylip2.tmp
do
    prefix=$(basename $f .phylip2.tmp)
    cat headC10.tmp $f > ${prefix}.phylip
done
for f in *C25.phylip2.tmp
do
    prefix=$(basename $f .phylip2.tmp)
    cat headC25.tmp $f > ${prefix}.phylip
done
for f in *C50.phylip2.tmp
do
    prefix=$(basename $f .phylip2.tmp)
    cat headC50.tmp $f > ${prefix}.phylip
done
for f in *C75.phylip2.tmp
do
    prefix=$(basename $f .phylip2.tmp)
    cat headC75.tmp $f > ${prefix}.phylip
done



#Cleaning the folder (keep the randoms)
echo "Saving random seeds"
echo "Cleaning temporary files"
rm matrix_combined.phylip
rm living_morpholog.phylip
rm fossil_morpholog.phylip
rm *.name
rm *.tmp
rm *.phy
rm *.test
rm *.phyl
rm *.phyli
rm *.phylip_tmp
rm *.data
rm *.subset
tar cf randoms.tar *.rand MorphoGammaRates.txt
rm MorphoGammaRates.txt
rm *.rand
echo "OUTPUT FILE NAMES:"
echo "Matrice Number _ NL+State Number NF+State number NC+State Number"
echo "example:"
echo "M2_L10N25C50"
echo "M2 = Matrice number 2"
echo "L10 = NL 10%"
echo "F25 = NF 25%"
echo "C50 = NC 50%"



echo "Matrices generations end:" >> Simulation_@.log
date >> Simulation_@.log


echo "##########################" >> Simulation_@.log

echo "OUTPUT FILE NAMES:
Matrice Number _ NL+State Number NF+State number NC+State Number
==========================
Example:
M2_L10F25C50
L10 = NL 10%
F25 = NF 25%
C50 = NC 50%
==========================
All the random numbers are stored in random.rar file" >> Simulation_@.log



#end