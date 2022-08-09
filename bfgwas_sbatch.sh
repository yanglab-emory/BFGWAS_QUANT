#!/bin/bash
#SBATCH --job-name=BFGWAS
#SBATCH --partition=yanglab,day-long-cpu,week-long-cpu,month-long-cpu,preemptable
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=./SlurmOut/BFGWAS.%A_%a.out.txt
#SBATCH --error=./SlurmErr/SBFGWAS.%A_%a.err.txt
#SBATCH --mail-user=jingjing.yang@emory.edu
#SBATCH --mail-type=END,FAIL

#  wet working directory
wkdir=$1
echo Working directory: $wkdir
cd $wkdir

mkfile=$2
echo Makefile: $mkfile

j=$3
echo Number of jobs in parallel: ${j}

###### Run BFGWAS_QUANT
BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT"
srun ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}



