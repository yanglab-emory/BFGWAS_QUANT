# rm(list=ls(all=TRUE))
library(data.table)
library(tidyverse)

####### Source Util Functions and set data directories
source("/projects/YangLabData/Software/BFGWAS_QUANT/bin/R_funcs.r")
setwd("/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults")

######## Compare results
#paramdata_bfgwas = LoadEMdata(filename="/projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp3.txt.gz", header = TRUE)
paramdata_bfgwas = LoadEMdata(filename="/projects/YangLabData/jyang/BFGWAS_Test/test_sim/Eoutput/paramtemp3.txt.gz", header = TRUE)
dim(paramdata_bfgwas)
sum(paramdata_bfgwas$Pi)

paramdata_bfgwas$Zscore = sqrt(paramdata_bfgwas$Chisq) * sign(paramdata_bfgwas$mBeta)
var(paramdata_bfgwas$Zscore)


## Manhantton plot
paramdata_bfgwas_sig <- filter(paramdata_bfgwas, Pi > 0.1)
dim(paramdata_bfgwas_sig)
paramdata_bfgwas_sig

###### Manhantton Plot
ggplot(paramdata_bfgwas, aes(x=POS, y = -log10(Pval))) +
	geom_point(shape = 21, fill = "blue", color = "blue") +
	geom_hline(yintercept=-log10(5e-8)) +
 	geom_point(data = paramdata_bfgwas_sig, aes(x=POS, y = -log10(Pval), color = Pi))  +
 	scale_color_gradient(low="yellow", high="red")  +
	facet_grid(cols = vars(CHR), scales = "free")
ggsave("/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/mp.pdf")

###### Effect size plot: Bayesian estimates vs. marginal estimates
ggplot(paramdata_bfgwas[paramdata_bfgwas$Beta>0, ], aes(x = mBeta, y = Beta, col = Pi)) +
	geom_point() + geom_abline(intercept=0, slope = 1) +
	scale_color_gradient(low="blue", high="red")
ggsave("/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/beta.pdf")

######## Results of hyper parameter estimates #######
test_hyp <- fread("/projects/YangLabData/jyang/BFGWAS_Test/test_sim/Eoutput/EM_result.txt", header = TRUE)
# test_hyp <- fread("/projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/EM_result.txt", header = TRUE)
print(test_hyp)

######## END ################










