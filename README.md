Bayesian Functional Genome-wide Association Study using Multivariate Quantitative Annotations (BFGWAS_QUANT)
=======
**BFGWAS_QUANT** is developed based on our previous tool [**BFGWAS**](https://github.com/yjingj/BFGWAS_SS) tool as proposed  by [Yang J. et. al, AJHG, 2017](https://doi.org/10.1016/j.ajhg.2017.08.002). **BFGWAS_QUANT** is based on a multivariable Bayesian variable selection model, accounting for multivariate quantitative functional annotations and LD to fine-map and prioritize GWAS hits.

- With **individual-level GWAS data**, **BFGWAS_QUANT** first computes single variant Z-score statistics (GWAS summary statistics) and LD correlation files from individual-level GWAS data, and then loads these two summary statistics files and known GWAS sample size to run MCMC.
- With **summary-level GWAS data**,  **BFGWAS_QUANT** uses Z-score statistics file and LD correlation files generated from a reference panel (*by BFGWAS_QUANT*) to run MCMC.
- **BFGWAS_QUANT** runs MCMC in parallel for multiple genome-blocks, which is handled through a `Makefile` generated by a PERL script `./bin/gen_mkf.pl`. EM-MCMC algorithm is handled based on dependency jobs in the generated `Makefile`. Users will only need to submit the generated `Makefile` as a job for the analysis results.


---

- [Software Installation](#software-installation)
- [Input Files](#input-files)
	- [Individual-level GWAS Data Input Files](#individual-level-gwas-data-input-files)
		- [1. Genotype VCF Files](#1-genotype-vcf-files)
		- [2. Phenotype File](#2-phenotype-file)
	- [Summary-level GWAS Data Input Files](#summary-level-gwas-data-input-files)
		- [1. GWAS Zscore File](#1-gwas-zscore-file)
		- [2. Reference LD File](#2-reference-LD-file)
		- [3. Annotation File](#3-annotation-file)
	- [Other Input Files](#other-input-files)
		- [1. Genome Block Prefix File](#1-genome-block-prefix-file)
		- [2. Prior Parameter File](#2-prior-parameter-file)
- [Example Usage](#example-usage)
	- [Step 1. Obtain GWAS Zscore and LD Files](#step-1-obtain-gwas-zscore-and-ld-files)
	- [Step 2. Generate Makefile](#step-2-generate-makefile)
	- [Step 3. Run Makefile](#step-3-run-makefile)
- [Output Files](#output-files)
- [Analyse BFGWAS_QUANT Results](#analyse-bfgwas_quant-results)
- [Remarks](#remarks)

---


## Software Installation

### 1. Compile **Estep_mcmc**
* Install required C++ libraries C++ libraries **zlib**, **gsl**, **eigen3**, **lapack**, **atlas**, **blas**. Please install these libraries to your system and include the library path `-I[path to libraries]` accordingly in the C++ compilation command line in the `Makefile`.

* Compile C++ library `./libStatGen/libStatGen.a` under your system by using the following commands:

```
cd BFGWAS_QUANT/libStatGen/;
make clean;
make
```

* Compile C++ source code for the executable file `./bin/Estep_mcmc` that will be used to run the Estep MCMC algorithm, by using the following commands under `BFGWAS_QUANT/` directory:

```
cd BFGWAS_QUANT/;
make clean ;
make
```

* Even though a compiled executable file `./bin/Estep_mcmc` from our cluster is provided on GITHUB, please still compile one for your own system. The `BFGWAS_QUANT/Makefile` might need to be adapted for one's own system.

### 2. Additional Requirements
* Tool [**TABIX**](https://www.htslib.org/doc/tabix.html) : used to segment VCF, GWAS Zscore, Annotation files.
* [**R**](https://www.r-project.org/) library [**data.table**](https://github.com/Rdatatable/data.table/wiki/Installation), [**tidyverse**](https://www.tidyverse.org/), [**optimx**](https://cran.r-project.org/web/packages/optimx/index.html) : used in the M-step R script `./bin/Mstep.r`.

## Input Files
**Example data files are provided under `./Example/ExData/`.**

### Individual-level GWAS Data Input Files

#### 1. Genotype VCF Files
- **[VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf)** are required for using individual-level GWAS data, or for generating reference LD correlation files. The VCF genotype files should be segmented into one VCF file per genome block (variants of the same chromosome should be in the same block), sorted by position, and then zipped by `bgzip` with file names `[filehead].vcf.gz`.
- All VCF genotype files should be stored under the same parent directory.
- Genotype files are supposed to be segmented based on LD. Genome blocks are expected to be approximately independent with ~5K -- 10K SNPs per block.
	- Segmentation information of **hg19** derived by [LDetect](https://bitbucket.org/nygcresearch/ldetect/src/master/) are provided in `./Example/ExData/Segmentations/lddetect_1KG_*_hg19.bed ` for EUR, AFR, and ASN populations.
	- Segmentation information of **GRCH38** are available from [LDblocks_GRCh38](https://github.com/jmacdon/LDblocks_GRCh38).

#### 2. Phenotype File
- **Phenotype file** is a two column text file, with sample IDs (the same IDs as in the genotype VCF files) in the first column and phenotype values in the second column
- No column header is needed
- Example file: `/Example/ExData/sim_pheno.txt`

	| ...    | ...  |
	|-------------|-------------|
	| 1000    | 0.64 |
	| 1005   | -0.55 |
	|...	| ... |

### Summary-level GWAS Data Input Files

#### 1. GWAS Zscore File
- **GWAS Zscore files** are also segmented based on LD, sorted by position, zipped by `bgzip` with `[filehead].Zscore.txt.gz`, which can be either generated using individual-level GWAS data, or constructed from GWAS summary data.
- The same segmentation information used to segment individual-level VCF files can also be used here to construct segmented GWAS Zscore files.
- **GWAS Zscore files** are tab separated text files with first 8 columns: `#CHROM POS ID REF ALT N MAF Z_SCORE`, denoting chromosome number, base pair coordinate, variant ID, reference allele, alternative allele, sample size, minor allele frequency, Z-score statistic values by single variant tests. Unknown sample size and minor allele frequency values can be denoted as `NA`.
	- `MAF` will be used to filter out rare variant if provided.
- Example files: `/Example/ExData/Zscore/`

	| #CHROM   | POS    | ID    | REF    | ALT    | N   | MAF    | Z_SCORE   |
	|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|
	| 19    | 29791987 | 19:29791987:C:A | C |A | 1893 | 0.021 | 0.505 |
	| 19   | 29792094 |  19:29792094:G:A | G | A | 1893 | 0.400 | 0.139 |
	|...	| ... |  ... |  ... |  ... |  ... |  ... |  ... |

#### 2. Reference LD File
- **Reference LD files** are segmented using the same segmentation information for segmenting **GWAS Zscore files** or **VCF files**, with names `[filehead].LDcorr.txt.gz`.
- **Reference LD files** have to be generated by `./bin/Estep_mcmc` from either individual-level GWAS genotype VCF files of test samples or a reference cohort such as [**1000 Genome**](https://www.internationalgenome.org/category/vcf/) of the same ancestry. Reference LD files generated using a reference cohort with different ancestry from the test GWAS cohort will cause errors.
-   **Reference LD files** are tab separated text files with 9 columns: `#ORDER CHROM POS ID REF ALT N MAF CORR`, denoting the order of variants, chromosome number, base pair coordinate, variant ID, reference allele, alternative allele, sample size, minor allele frequency, correlation between variant in the current row and variants within a specified right hand size `$LDwindow` of the current variant (separated by `,`, started with `1.0` as standardized correlation).
-  Example files: `/Example/ExData/RefLD/`. For example, the correlation between `19:29791987:C:A` and `19:29792094:G:A` is `-0.123`.

	| #ORDER | CHROM   | POS    | ID    | REF    | ALT    | N   | MAF    | CORR   |
	|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|
	|0| 19    | 29791987 | 19:29791987:C:A | C |A | 1893 | 0.021 | 1.0,-0.123,-0.202,... |
	|1| 19   | 29792094 |  19:29792094:G:A | G | A | 1893 | 0.400 | 1.0,0.640,0.117,... |
	|...	| ... |  ... |  ... |  ... |  ... |  ... |  ... |

#### 3. Annotation File
- **Annotation files** are also segmented using the same segmentation information for segmenting **GWAS Zscore files** or **VCF files**, with names `Anno_[filehead].txt.gz`.
- **Annotation files** are tab separated text files with 6 columns: `#CHROM POS ID REF ALT Annotation`, denoting chromosome number, base pair coordinate, variant ID, reference allele, alternative allele, multivariate quantitative annotations (separated by `,`).

	| #CHROM   | POS    | ID    | REF    | ALT    | Annotation   |
	|-------------|-------------|-------------|-------------|-------------|-------------|	| 19    | 29791987 | 19:29791987:C:A | C |A | 0,0,0,-1.5646 |
	| 19   | 29792094 |  19:29792094:G:A | G | A | 0,0,0,-2.2281  |
	|...	| ... |  ... |  ... |  ... |  ... |  ... |  ... |

### Other Input Files 
#### 1. Genome Block Prefix File

* **Genome block prefix file** of VCF genotype files as in `./Example/ExData/filehead_4block.txt` is required. Each row of the list file is the file head of the VCF file of one genome block as in `[filehead].vcf.gz`. Note that the VCF file extension suffix `.vcf.gz` should not be included.
* The same set of file heads as in the **genome block prefix file** should be used to name **VCF, Zscore, LD, Annotation files**.

#### 2. Prior Parameter File
- **Prior parameter file** includes the values of initial enrichment coefficients **`a`** (separated by `,`) and a fixed value for `tau_beta`. See example file as in ``./Example/ExData/hypval_4anno.txt`. 
- The first value in the **`a`** row is for `alpha_0` with a recomended value in `[-13.8, -9]`, which controls for the benchmark of model sparsity when enrichment parameters are `0` and will not be updated during M-step. 
- Default `tau_beta` value is recommended to be set as `0.1` if individual-level GWAS data are used, and `1` if summary-level GWAS data are used.


	| #hyper_parameter   |   |
	|-------------|-------------|
	| **a**   | -13.8,0,0,0,0 |
	|tau_beta   | 0.1 |



## Example Usage

### **Example commands: `./test_script.sh`**

#### Step 1. Obtain GWAS Zscore and LD Files
- Shell script `./bin/GetRefLD.sh` can be used to generate GWAS Zscore and LD files sequentially
- The following commands can be used to generate GWAS Zscore and LD files for one genome block:

	```
	BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT"
	pheno=${BFGWAS_dir}/Example/ExData/sim_pheno.txt
	geno_dir=/home/jyang51/YangLabData/SharedData/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF
	LDwindow=1000000
	maf=0.001
	GTfield=GT
	line="WGS_1898_samples_CHR_19_13471127_14486347"
	${BFGWAS_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} \
	                    -maf ${maf} -GTfield ${GTfield} \
	                    -o ${line} -LDwindow ${LDwindow} -saveSS -zipSS
	```

#### Step 2. Generate Makefile
- A Makefile will be generated by `./bin/gen_mkf.pl` to rap MCMC jobs (E-steps) and updating enrichment coefficient estimates (M-steps), and enable parallel computation.
- In the following example commands, `Nsample` sets GWAS sample size, `Anum` sets the number of annotations, `em` sets number of EM iterations, `Nburnin` sets the number of burn-in MCMC iterations, `Nmcmc` sets the number of MCMC iterations, `anno_dir` specifies the directory of all annotation files,

```
BFGWAS_dir="/projects/YangLabData/Software/BFGWAS_QUANT" # Tool directory

filehead=${BFGWAS_dir}/Example/ExData/filehead_4block.txt # Genome block prefix file
Zscore_dir=${BFGWAS_dir}/Example/ExData/Zscore # Zscore file directory
LDdir=${BFGWAS_dir}/Example/ExData/RefLD # LD file directory

Nsample=1893 # GWAS sample size
Anum=4 # Number of annotations
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number

anno_dir=${BFGWAS_dir}/Example/ExData/Anno # Annotation file directory
hfile=${BFGWAS_dir}/Example/ExData/hypval_4anno.txt #  Initial prior parameter values

wkdir=${BFGWAS_dir}/Example/Test_wkdir # Working directory
cd $wkdir
mkfile=${wkdir}/simu_BFGWAS.mk ## Makefile directory

########### Generate make file with all jobs
${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf 0.01 --Nsample ${Nsample} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}
```


#### Step 3. Run Makefile
- Run on computation node:

```
j=4
make -f ${wkdir}/simu_BFGWAS.mk clean; # Clean target files
${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}
```
- Run by `sbatch`:

```
j=4
make -f ${wkdir}/simu_BFGWAS.mk clean; # Clean target files
sbatch ${BFGWAS_dir}/bfgwas_sbatch.sh ${wkdir} ${mkfile} ${j}
```

## Output Files
- **Variant specific estimates** of effect sizes (column `Beta`) and causal posterior probabilities (CPP, column `Pi`) are provided under `${wkdir}/Eoutput/paramtemp${em}.txt.gz`.
	- Column `mBeta` : marginal genetic effect size estimate by single variant regression model
	- Column `ChisqTest` : chisquare test statistic value by single variant test
	- Column `Pval_svt` : p-value by single variant test
	- Column 	`Rank` : rank of p-values from the smallest to the largest per genome block
	- Columns `Anno_*` : annotation values per variant
- **Enrichment coefficient estimates** (column `avec`, separated by `,`) are provided in `${wkdir}/Eoutput/EM_result.txt`


## Analyse **BFGWAS_QUANT** Results
- Example R script for analyzing the **BFGWAS_QUANT** results is provided in `./Example/AnalyzeResults/Analysis.r`
- Example plots are provided under `./Example/AnalyzeResults/`
- Example Rscript for making Manhattan plot for millions of SNPs: `./bin/mp_plot.r`


## Remarks
- **BFGWAS_QUANT** tool is derived by assuming all test GWAS samples are of the same ancestry. It is critical to use reference LD of the same population as the GWAS test data.
- Around `32GB` memory might be needed to generate reference LD files.
- GWAS sample size `Nsample` is recommended to be set as the sample size that are used to generate the `LDcorr` files.
- GWAS `Zscore` statistics can be easily derived by using the GWAS `p-values` and effect size signs.
- If **BFGWAS_QUANT** tool fail to complete all EM iterations, one could use the output files generated by the last success EM iteration. The failed scenarios are likely due to the lack of associations.
