######## Please use our data disk (/mnt/YangFSS/data) for input and output of all BFGWAS jobs
######## Example scripts are under /home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts
######### Tool directory


for chr in {100..100};
do
  ######## Specify file directories and variable values for BFGWAS
  N=1893  # sample size
  number=$(($chr + 1))
  pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/simulation_10casual_h2_25/simulation_info.txt | awk '{print $5}')
  # genome-block names
  filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19.bed
  #filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19_new.bed

  # phenotype directory
  pheno=/mnt/YangFSS/data/jchen/simulation_10casual_h2_25/wkdir_${chr}/pheno.txt
  # genotype directory
  geno_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF

  ######### Please first generate summary statistics (LD and score statistics files) files for best computation speed
  # LD coefficient directory
  LD_dir=/mnt/YangFSS/data/jchen/simulation_20blocks/LD_statistic
  mkdir ${LD_dir}
  # Score statistics directory
  Score_dir=/mnt/YangFSS/data/jchen/simulation_10casual_h2_25/wkdir_${chr}/Score
  if [[ ! -e $Score_dir ]]; then
    mkdir $Score_dir
  elif [[ ! -d $Score_dir ]]; then
    echo "$Score_dir already exists" 1>&2
  fi
  #LDwindow=1000000
  LDwindow=1

  # Submit Array jobs with GetRefLD.sh
  qsub -j y -wd ${Score_dir} -N Score_${chr} -t 1-20 -tc 100 -pe smp 5 /home/jchen/bfGWAS/GetRefLD_Jingjing.sh ${geno_dir} ${pheno} ${filehead} ${LD_dir} ${Score_dir} ${LDwindow}
done


####### Command to run BFGWAS for genome-wide blocks
##  run_BFGWAS.sh will call the perl script gen_mkf.pl to generate a makefile with all BFGWAS_SS jobs run by run_make.sh
## The executible file /home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc was run by run_Estep.sh, which is called by make
## Then clean all previous outputs and re-run the generated Makefile
## Please test with ~10 genome-blocks instead of running the whole genome-wide data

for chr in {1..1};
do
  bfGWAS_SS_dir="/home/jchen/bfGWAS/BFGWAS_QUAN_ver2"
  ######## Specify file directories and variable values for BFGWAS
  N=1893  # sample size
  number=$(($chr + 1))
  pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_10casual_h2_25/simulation_info.txt | awk '{print $5}')
  # genome-block names
  filehead=/mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_20blocks/FileHead_1703_anno_chr19.bed
  #filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19_new.bed

  # phenotype directory
  pheno=/mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_10casual_h2_25/wkdir_${chr}/pheno.txt
  # genotype directory
  geno_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF

  ######### Please first generate summary statistics (LD and score statistics files) files for best computation speed
  # LD coefficient directory
  LD_dir=/mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_20blocks/LD_statistic
  # Score statistics directory
  Score_dir=/mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_10casual_h2_25/wkdir_${chr}/Score
  LDwindow=1
  Anum=4
  hfile=/home/jchen/bfGWAS/BFGWAS_QUAN_ver2/hypval.current

  annoDir=/mnt/YangFSS/data/jchen/BFGWAS_simulation/simulation_20blocks/Anno

  # Set up Working directory
  wkdir=/mnt/YangFSS/data/jchen/BFGWAS_simulation/BFGWAS_QUAN_Test
  #wkdir=/mnt/YangFSS/data/jchen/simulation_10casual_h2_25/wkdir_${chr}
  if [[ ! -e $wkdir ]]; then
    mkdir $wkdir
  elif [[ ! -d $wkdir ]]; then
    echo "$wkdir already exists" 1>&2
  fi
  cd $wkdir

  ###### Call gen_mkf.pl to generate Make file
  # directory for run_Estep.sh and run_Mstep.sh
  EMdir=/home/jchen/bfGWAS/BFGWAS_QUAN_ver2/bin
  # Specify computation specs
  em=3 # EM steps
  burnin=10000  # Burn-in iterations in MCMC
  Nmcmc=10000  # MCMC iteration number
  mkfile=${wkdir}/${study}_BFGWAS.mk

  /home/jchen/bfGWAS/BFGWAS_QUAN_ver2/bin/gen_mkf.pl \
  --EMdir ${EMdir} --hyp ${hfile} \
  -n ${N} --pv ${pheno_var} \
  -w ${wkdir} --geno sumstat \
  -f ${filehead} -l local \
  --ad ${annoDir} \
  --LDdir ${LD_dir} --Scoredir ${Score_dir} \
  -j BFGWAS_${chr} --em ${em} -b ${burnin} -N ${Nmcmc} \
  --mf ${mkfile} --AnnoNumber ${Anum}

done

work_dir=/mnt/YangFSS/data/jchen/BFGWAS_simulation/BFGWAS_QUAN_Test
j=12 # Number of cores to run the job in parallel
qsub -q b.q -j y -pe smp ${j} -t 1-1 -tc 4 -wd ${work_dir} -N BFGWAS /mnt/YangFSS/data/jchen/BFGWAS_simulation/BFGWAS_QUAN_Test/run_make.sh ${j}
