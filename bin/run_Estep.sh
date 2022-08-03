#!/usr/bin/bash

wkdir=$1
line=$2
N=$3
pv=$4
burnin=$5
Nmcmc=$6

# summary stat directories
LDdir=$7
Scoredir=$8

# Input files
hfile=$9
annoDir=${10}
#annoCode=${11}
AnnoNumber=${11}
initype=${12}

# echo Curerent working directory
# pwd

echo Run Estep with $line $N $pv

/home/jchen/bfGWAS/BFGWAS_QUAN_ver2/bin/Estep_mcmc -inputSS \
-score ${Scoredir}/${line}.score.txt.gz \
-LDcorr ${LDdir}/${line}.LDcorr.txt.gz \
-a ${annoDir}/Anno_${line}.gz \
-hfile ${hfile} \
-n ${N} -pv ${pv} -maf 0.01 -r2 0.001 -bvsrm -smin 0 -smax 10 -win 100 \
-o ${line} -w ${burnin} -s ${Nmcmc} -initype ${initype} -AnnoNumber ${AnnoNumber}\
> ${wkdir}/OUT/${line}.output.txt

exit
