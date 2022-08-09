library(data.table)
library(tidyverse)

######## Read BFGWAS Results #############
paramdata_bfgwas = fread("/projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp3.txt.gz", header = TRUE)
dim(paramdata_bfgwas)
sum(paramdata_bfgwas$Pi)

##### Subsetting for significant SNPs with CPP > 0.1 #############
paramdata_bfgwas_sig <- filter(paramdata_bfgwas, Pi > 0.1)
dim(paramdata_bfgwas_sig)
sum(paramdata_bfgwas_sig$Pi) # Estimating total number of causal SNPs

###### Manhantton Plot for thousands of SNPs #############
ggplot(paramdata_bfgwas, aes(x=POS, y = -log10(Pval_svt))) +
	geom_point(shape = 21, fill = "blue", color = "blue") +
	geom_hline(yintercept=-log10(5e-8)) +
 	geom_point(data = paramdata_bfgwas_sig, aes(x=POS, y = -log10(Pval_svt), color = Pi))  +
 	scale_color_gradient(low="yellow", high="red")  +
	facet_grid(cols = vars(`#CHR`), scales = "free")
ggsave("/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/mp.pdf")

###### Effect size plot: Bayesian estimates vs. marginal estimates ############
ggplot(paramdata_bfgwas, aes(x = mBeta, y = Beta)) +
	geom_point(shape = 21, fill = "blue", color = "blue") +
	geom_abline(intercept=0, slope = 1) +
	geom_point(data = paramdata_bfgwas_sig, aes(x = mBeta, y = Beta, color = Pi))  +
	scale_color_gradient(low="yellow", high="red")
ggsave("/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/beta.pdf")

######## Results of hyper parameter estimates #######
test_hyp <- fread("/projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/EM_result.txt", header = TRUE)
print(test_hyp)

Enrich_df <- test_hyp$avec %>% str_split(",") %>% unlist() %>% matrix(nrow = nrow(test_hyp)) %>% data.frame()
colnames(Enrich_df) <- paste0("Anno", 1:ncol(Enrich_df))
rownames(Enrich_df) <- paste0("EM", 1:nrow(Enrich_df))
print(Enrich_df)

########  Manhantton plot for millions of SNPs ########
sig_level="5e-8"
paramFile="/projects/YangLabData/Software/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp3.txt.gz"
out_prefix="/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/sim"
log_file="/projects/YangLabData/Software/BFGWAS_QUANT/Example/AnalyzeResults/sim_mp.log.txt"

system(paste(
"R CMD BATCH --no-save --no-restore ",
"'--args",
" paramFile=\"",paramFile,"\"",
" sig_level=\"", sig_level,"\"",
" out_prefix=\"", out_prefix,"\"",
"' ",
" /projects/YangLabData/Software/BFGWAS_QUANT/bin/mp_plot.r ",
log_file, sep=""), wait=FALSE)

######## END ################










