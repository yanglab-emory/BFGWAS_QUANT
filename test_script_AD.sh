#########################################################
############## Test on GWAS AD data #####################
#########################################################

BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT"

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
j=30 # Number of cores to request
# srun ${BFGWAS_SS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

sbatch ${BFGWAS_dir}/bfgwas_sbatch.sh ${wkdir} ${mkfile} ${j}

squeue
sinfo

################ Test a single block ################


################ Test M-step ########################
Rscript --vanilla /projects/YangLabData/Software/BFGWAS_QUANT/bin/Mstep.r /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/hyptemp0.txt 0 2 1 1893 /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/hypval.current /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/paramtemp0.txt.gz /home/jyang51/YangLabData/jyang/BFGWAS_Test/test_rosmap/Eoutput/EM_result.txt



##### ROC Plot with Mayo Clinic Data

## Extract Genotype data
paramFile="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/Eoutput/paramtemp3.txt"
head -n1 $paramFile | awk '{print $1,$2,$4,$5,$7,$8,$9,$11}' > /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/AD_sig.txt
tail -n+2 ${paramFile} | awk '{if($7>0.001){print $1,$2,$4,$5,$7,$8,$9,$11}}' | sort -nk1 -nk2 >> /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/AD_sig.txt

for chr in {1..22} ; do
# zcat MayoLOADGWAS_CHR${chr}_filtered.dose.geno.gz | bgzip > MayoLOADGWAS_CHR${chr}_filtered.dose.geno.bgzip &
# mv MayoLOADGWAS_CHR${chr}_filtered.dose.geno.bgzip MayoLOADGWAS_CHR${chr}_filtered.dose.geno.txt.gz
tabix -p vcf -f MayoLOADGWAS_CHR${chr}_filtered.dose.geno.txt.gz &
done

tabix MayoLOADGWAS_CHR22_filtered.dose.geno.txt.gz 22:17066020-17066020 | cut -f 1-10

############# Test One block
# WGS_1898_samples_CHR_19_19877471_20905757
# WGS_1898_samples_CHR_19_16374416_18409862

/home/jyang/GIT/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/GWAS_AD_Data/WGS_1898_samples_CHR_19_19877471_20905757.Zscore.txt.gz -LDcorr /mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/RefLD/WGS_1898_samples_CHR_19_19877471_20905757.LDcorr.txt.gz -a /mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/Segmented_Quant_Anno/Anno_WGS_1898_samples_CHR_19_19877471_20905757.txt.gz -hfile /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/hypval.current -n 10000 -maf 0.001 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_19877471_20905757 -w 10000 -s 10000 -AnnoNumber 10 -seed 2022


Rscript --vanilla /home/jyang/GIT/BFGWAS_QUANT/bin/Mstep.r /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/Eoutput/hyptemp0.txt 0 1.0001 1 10000 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/hypval.current /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/Eoutput/paramtemp0.txt /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/Eoutput/EM_result.txt


############ Generate LDcorr matrix and GWAS Zscore for simulation data
pheno=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF/fake_pheno.txt
#pheno=/mnt/YangFSS/data/jchen/BFGWAS_simulation/01-SimulationOld/simulation_10casual_h2_50/wkdir_1/pheno.txt
filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Temp_Zscore
# LDwindow=1000000

## Run all blocks sequencially
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

### Submit Array jobs
qsub -q b.q -j y -pe smp 6 -wd ${wkdir} -N GetRefLD -t 1-20 -tc 14 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/run_GetRefLD_array.sh

###### Extract AD GWAS data per genome block
cd /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/GWAS_AD_Data
cat /mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF/lddetect_1KG_EUR.bed | while read line ; do
echo $line
chr=`echo ${line} | awk '{print $1}'`
start_pos=`echo ${line} |awk '{print $2}'`
end_pos=`echo ${line} | awk '{print $3}'`
file_prefix=`echo ${line} | awk '{print $4}'`
# echo ${chr}:${start_pos}-${end_pos}
tabix -h /mnt/YangFSS/data2/Sumstats/AD/Douglas_NG_2021/PGCALZ2sumstatsExcluding23andMe.txt.gz ${chr}:${start_pos}-${end_pos} > temp.ZScore.txt
## Format Zscore statistics file
echo -e "#CHROM\tPOS\tID\tREF\tALT\tN\tMAF\tZscore" > temp.txt
awk 'NR>1{print $1"\t"$2"\t"$1":"$2":"$3":"$4"\t"$3"\t"$4"\t"$7"\tNA""\t"$5}' temp.ZScore.txt  >> temp.txt
mv temp.txt ${file_prefix}.Zscore.txt
bgzip -f ${file_prefix}.Zscore.txt
tabix -b2 -e2 -S1 -s1 -f ${file_prefix}.Zscore.txt.gz
done

#### Get median sample size
zcat /mnt/YangFSS/data2/Sumstats/AD/Douglas_NG_2021/PGCALZ2sumstatsExcluding23andMe.txt.gz | cut -f7 | tail -n+2 | sort -n  | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

