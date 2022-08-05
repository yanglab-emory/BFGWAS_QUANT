BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT"

############# Test Simulation Data with 4 blocks
##### Set directories for filehead, GWAS summary data, Annotation Data, and Reference LD files
wkdir=${BFGWAS_dir}/Example/Test_wkdir
cd $wkdir
filehead=${BFGWAS_dir}/Example/ExData/filehead_4block.txt
Zscore_dir=${BFGWAS_dir}/Example/ExData/Zscore
LDdir=${BFGWAS_dir}/Example/ExData/RefLD
anno_dir=${BFGWAS_dir}/Example/ExData/Anno

############# Test Simulation Data with  20 blocks
 LDdir=/home/jyang51/YangLabData/SharedData/AMP-AD/BFGWAS_LDdir
 filehead=/projects/YangLabData/jyang/BFGWAS_Test/sim_20blocks_filehead.txt
 Zscore_dir=/home/jyang51/YangLabData/jyang/BFGWAS_Test/test_sim/Sim_Zscore
 anno_dir=/home/jyang51/YangLabData/jyang/BFGWAS_Test/test_sim/Anno

############################################################
############ Set input argument values
Nsample=1893 # sample size
Anum=4 # 4 annotations
hfile=${BFGWAS_dir}/Example/ExData/hypval_4anno.txt #  Initial prior parameter values
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number

wkdir=/projects/YangLabData/jyang/BFGWAS_Test/test_sim; cd $wkdir
mkfile=${wkdir}/simu_BFGWAS.mk ## Make file directory

########### Generate make file with all jobs
${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf 0.01 --Nsample ${Nsample} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

############ Step 3 ###############
## Clean all jobs in the make file when you need rerun everything
# make -f ${wkdir}/simu_BFGWAS.mk clean; rm make.output make.err Rout.txt hypval.current output/** OUT/** Eoutput/**
# rm ${wkdir}/SlurmErr/** ${wkdir}/SlurmOut/**

######### Submit the job for running the makefile
# j=4 # Number of cores to request
j=20 # Number of cores to request
# ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

sbatch ${BFGWAS_dir}/bfgwas_sbatch.sh ${wkdir} ${mkfile} ${j}

squeue
sinfo

###### Test E-step of a single block
/projects/YangLabData/Software/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/Zscore/WGS_1898_samples_CHR_19_29790947_30727954.Zscore.txt.gz -LDcorr /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/RefLD/WGS_1898_samples_CHR_19_29790947_30727954.LDcorr.txt.gz -a /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/Anno/Anno_WGS_1898_samples_CHR_19_29790947_30727954.txt.gz -hfile /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/hypval.current -maf 0.01 -n 1893 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_29790947_30727954 -w 10 -s 10 -AnnoNumber 4 -seed 2022

###### Test M-step
Rscript --vanilla /projects/YangLabData/Software/BFGWAS_QUANT/bin/Mstep.r /projects/YangLabData/jyang/BFGWAS_Test/test_sim/Eoutput/hyptemp1.txt 1 2 1 1893 /projects/YangLabData/jyang/BFGWAS_Test/test_sim/hypval.current /projects/YangLabData/jyang/BFGWAS_Test/test_sim/Eoutput/paramtemp1.txt.gz /projects/YangLabData/jyang/BFGWAS_Test/test_sim/Eoutput/EM_result.txt >> /projects/YangLabData/jyang/BFGWAS_Test/test_sim/Rout.txt


############ Generate LDcorr matrix and GWAS Zscore for simulation data (20 blocks)
pheno=${BFGWAS_dir}/Example/ExData/sim_pheno.txt
geno_dir=/home/jyang51/YangLabData/SharedData/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF

# filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
#Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore
#LDdir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD
# LDwindow=1000000

## Run all blocks sequencially
## >32GB memory might be needed
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

qsub -q b.q -j y -pe smp 4 -wd ${wkdir} -N GetRefLD -t 1-20 -tc 14 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/run_GetRefLD_array.sh