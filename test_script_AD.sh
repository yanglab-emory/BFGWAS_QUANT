#########################################################
############## Test on GWAS AD data #####################
#########################################################

BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT_08_03"

############ Step 1 ###############
LDdir=/home/jyang51/YangLabData/SharedData/AMP-AD/BFGWAS_LDdir
anno_dir=/home/jyang51/YangLabData/SharedData/AMP-AD/Segmented_Quant_Anno

#filehead=/home/jyang51/YangLabData/jyang/BFGWAS_Test/sim_20blocks_filehead.txt
filehead=/home/jyang51/YangLabData/SharedData/AMP-AD/FileHead_1703.txt

Nsample=1893 # Set as reference data (LD) sample size
Anum=10 # 10 annotations
maf=0.001
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number

Zscore_dir=/home/jyang51/YangLabData/SharedData/AMP-AD/GWAS_AD_Data; wkdir="/home/jyang51/YangLabData/jyang/BFGWAS_Test/test_AD"
# Zscore_dir=/projects/YangLabData/jyang/BFGWAS_Test/test_rosmap/Zscore; wkdir="/home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap"
cd $wkdir
mkfile=${wkdir}/BFGWAS_AD.mk
hfile=/home/jyang51/YangLabData/jyang/BFGWAS_Test/hypval_10anno.txt


############ Step 2 ###############
########### Generate make file with all jobs
${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf ${maf} --Nsample ${Nsample} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

############ Step 3 ###############
## Clean all jobs in the make file when you need rerun everything
make -f ${mkfile} clean
rm make.output make.err Rout.txt hypval.current output/** OUT/** Eoutput/**

######### Submit the job for running the makefile
j=100 # Number of cores to request
# srun ${BFGWAS_SS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

sbatch ${BFGWAS_dir}/bfgwas_sbatch.sh ${wkdir} ${mkfile} ${j}

squeue
sinfo

################ Test a single block ################


################ Test M-step ########################
Rscript --vanilla /projects/YangLabData/Software/BFGWAS_QUANT_08_03/bin/Mstep.r /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/hyptemp0.txt 0 2 1 1893 /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/hypval.current /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/paramtemp0.txt.gz /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/EM_result.txt





