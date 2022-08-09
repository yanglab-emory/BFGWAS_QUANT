BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT"

############# Test Simulation Data with 4 blocks
##### Set directories for filehead, GWAS summary data, Annotation Data, and Reference LD files
filehead=${BFGWAS_dir}/Example/ExData/filehead_4block.txt
Zscore_dir=${BFGWAS_dir}/Example/ExData/Zscore
LDdir=${BFGWAS_dir}/Example/ExData/RefLD

############ Step 1 ###############
####  Generate LDcorr matrix and GWAS Zscore for simulation data (4 blocks)
pheno=${BFGWAS_dir}/Example/ExData/sim_pheno.txt
geno_dir=/home/jyang51/YangLabData/SharedData/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF
LDwindow=1000000

#### Run all blocks sequencially; >32GB memory might be needed
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

###### Run one blocks
line="WGS_1898_samples_CHR_19_13471127_14486347"
${BFGWAS_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} \
                    -maf ${maf} -GTfield ${GTfield} \
                    -o ${line} -LDwindow ${LDwindow} -saveSS -zipSS

############ Step 2 ###############
#### Set input argument values
Nsample=1893 # sample size
Anum=4 # 4 annotations
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number
anno_dir=${BFGWAS_dir}/Example/ExData/Anno
hfile=${BFGWAS_dir}/Example/ExData/hypval_4anno.txt #  Initial prior parameter values

wkdir=${BFGWAS_dir}/Example/Test_wkdir
cd $wkdir
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
## Clean all jobs in the make file when you need re-run everything
make -f ${wkdir}/simu_BFGWAS.mk clean;
rm make.output make.err Rout.txt hypval.current

######### Submit the job for running the makefile
j=4 # Number of cores to request

## Run makefile
# ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

## Submit job to run makefile by sbatch
sbatch ${BFGWAS_dir}/bfgwas_sbatch.sh ${wkdir} ${mkfile} ${j}

#################################################################
###### Test E-step of a single block
/projects/YangLabData/Software/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/Zscore/WGS_1898_samples_CHR_19_28557893_29790947.Zscore.txt.gz -LDcorr /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/RefLD/WGS_1898_samples_CHR_19_28557893_29790947.LDcorr.txt.gz -a /projects/YangLabData/Software/BFGWAS_QUANT/Example/ExData/Anno/Anno_WGS_1898_samples_CHR_19_28557893_29790947.txt.gz -hfile /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/hypval.current -maf 0.01 -n 1893 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_28557893_29790947 -w 10000 -s 10000 -AnnoNumber 4 -seed 2022


###### Test M-step for one EM iteration
Rscript --vanilla /projects/YangLabData/Software/BFGWAS_QUANT/bin/Mstep.r /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/hyptemp0.txt 0 2 1 1893 /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/hypval.current /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp0.txt.gz /projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/EM_result.txt



