#!/usr/bin/bash

#######################################################################
### Input Arguments for Generate GWAS Summary Stat and LD Files
#######################################################################
# --wkdir: Working directory with writing access
# --toolE: Emcmc executible file directory
# --gene_dir: Directory for all genotype (VCF) files
# --pheno : Path for phenotype file
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --GTfield: Genotype data format: "GT" or "ES"
# --filehead: Path for the text file of all file head names of genome block files
# --LDdir: Directory to save all LD files
# --Zscore_dir: Directory to save all GWAS summary score statistic files
# --LDwindow: Window size to generate the LD correlation values
# --maf: MAF threshold to exclude rare variants
#######################################################################

VARS=`getopt -o "" -a -l \
wkdir:,toolE:,geno_dir:,pheno:,genofile_type:,GTfield:,filehead:,LDdir:,Zscore_dir:,LDwindow:,maf: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Terminating....." >&2
    exit 1
fi

eval set -- "$VARS"

while true
do
    case "$1" in
	    --wkdir|-wkdir) wkdir=$2; shift 2;;
        --toolE|-toolE) toolE=$2; shift 2;;
        --geno_dir|-geno_dir) geno_dir=$2; shift 2;;
        --pheno|-pheno) pheno=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
		--GTfield|-GTfield) GTfield=$2; shift 2;;
        --filehead|-filehead) filehead=$2; shift 2;;
        --LDdir|-LDdir) LDdir=$2; shift 2;;
        --Zscore_dir|-Zscore_dir) Zscore_dir=$2; shift 2;;
        --LDwindow|-LDwindow) LDwindow=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
LDwindow=${LDwindow:-1000000} # 1000KB
LDwindow_i=${LDwindow}
maf=${maf:-0.001}
genofile_type=${genofile_type:-vcf}
GTfield=${GTfield:-GT}

#### Create output directory for LD and Score summary stat files if not existed
mkdir -p $wkdir
cd $wkdir
mkdir -p ${LDdir}
mkdir -p ${Score_dir}

####################

##### run ./BFGWAS/bin/Estep_mcmc to generate GWAS summary score statistic and LD files
echo Generate GWAS summary Zscore statistic and LD correlation files with genome blocks under: ${geno_dir}
echo Genotype file type: $genofile_type
echo Phenotype: ${pheno}
echo Genetic variants with MAF greater than ${maf}
echo Under working directory: $wkdir
echo Input LDwidnow: ${LDwindow}

if [ ! -s ${pheno} ] ; then
	echo Phenotype file $pheno dose not exist or is empty. Please provide a phenotype file.
	exit 1
fi

cat ${filehead} | while read line ; do
	echo Genome block: $line
	LDwindow_i=${LDwindow}

	# Set LDwindow = 1 if LD file exit; Exit this script if Score statistic file exist
	if [ -s ${LDdir}/${line}.LDcorr.txt.gz ] ; then
		echo ${LDdir}/${line}.LDcorr.txt.gz exists!
		LDwindow_i=1
	fi

	if [ "${genofile_type}" == "vcf" ] ; then
		if [ -s ${geno_dir}/${line}.vcf.gz ] ; then
			echo Genotype field: $GTfield
			echo LDwidnow: ${LDwindow_i}
			### With input genotype file in VCF format
			${toolE} -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} -maf ${maf} -GTfield ${GTfield} \
					-o ${line} -LDwindow ${LDwindow_i} -saveSS -zipSS
		else
			echo Genotype block file ${geno_dir}/${line}.vcf.gz dose not exist or is empty. Please check.
			exit 1
		fi
	elif [ "${genofile_type}" == "geno" ] ; then
			if [ -s ${geno_dir}/${line}.geno.gz ] ; then
				echo LDwidnow: ${LDwindow_i}
				### With input genotype file in dosage format
				${toolE} -g ${geno_dir}/${line}.geno.gz -p ${pheno} -maf ${maf} \
						-o ${line} -LDwindow ${LDwindow_i} -saveSS -zipSS
			else
				echo Genotype block file ${geno_dir}/${line}.geno.gz dose not exist or is empty. Please check.
				exit 1
			fi
	else
		echo Please specify \'--genofile_type vcf\' ...
		exit 1
	fi

	echo Successfully generate GWAS summary score statistic and LD files for genome block $line.

##### Copy summary stat back


	if [ ! -s ${LDdir}/${line}.LDcorr.txt.gz ] ; then
		rsync  ./output/${line}.LDcorr.txt.gz* ${LDdir}/
	fi

	if [ ! -s ${Zscore_dir}/${line}.Zscore.txt.gz ] ; then
		rsync  ./output/${line}.Zscore.txt.gz* ${Zscore_dir}/
	else
		echo ${Zscore_dir}/${line}.Zscore.txt.gz exists.
	fi

	## Remove temperary output directory
	rm -f ./output/${line}*

done

exit 0



