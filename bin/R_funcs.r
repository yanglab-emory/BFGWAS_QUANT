####### Functions for logistic model
options(stringsAsFactors=F)
library(tidyverse)
library(data.table)

##########
calcMAF <- function(x) {
    n = sum(!is.na(x))
    if(n >0){
        maf  = sum(x, na.rm = TRUE) / (2 * sum(!is.na(x))) # unfold maf
    }else{
        maf = NA # all samples has NA genotypes
    }
	return(maf)
}

require(scales)
cubic_trans = function() trans_new("cubic", function(x) x^(1/3), function(x) x^(3))

## Read Eoutput/hyptemp*.txt
LoadEMdata <- function(filename, header = FALSE){
    paramdata = fread(filename, sep = "\t", header = header)
    setnames(paramdata, c("CHR", "POS", "ID", "REF", "ALT", "MAF", "Pi", "Beta", "mBeta", "Chisq", "Pval", "Rank", "Anno"))
    setkey(paramdata, "ID")
    return(paramdata)
}

# Read Eoutput/EM_result.txt
LoadEMhyp <- function(filename, header=FALSE){
    EMhyp_data <- read.table(filename, sep = "\t", header = header)
    return(EMhyp_data)
}

######## Revise the following functions
getCI <- function(est, est_se, alpha){
    z_alpha <- -qnorm((1-alpha)/2)
    pi_low <- est - z_alpha * est_se
    if(pi_low < 0) pi_low = 0
    pi_up <- est + z_alpha * est_se
    return(c(pi_low, pi_up))
}


### Load Annotation file
LoadAnnodata <- function(filename, header = FALSE){
    Annodata = fread(filename, sep = "\t", header = header)
    setnames(Annodata, c("chr", "bp", "ID", "ref", "alt", "Annotation"))
    return(Annodata)
}

getPi <- function(a, A){
	a <- matrix(a, ncol = 1)
	return(1/(1 + exp(-A %*% a)))
}

loglike_lambda <- function(lambda, lambda_bar, r, A, S){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	lambda <- lambda - lambda_bar
	l = sum(r * log(pi_pre) + (1-r) * log(1 - pi_pre)) - 0.5 * t(lambda) %*% S %*% (lambda)
	return(l)
}

loglike_lambda_g <- function(lambda, lambda_bar, r, A, S){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	lambda <- lambda - lambda_bar
	g = t(A) %*% (r - pi_pre) - S %*% lambda
	return(g)
}

### Hessian matrix, i.e., Fisher's information
FisherInfo_lambda <- function(lambda, A){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	h=diag(0, dim(A)[2])
	for(i in 1:dim(A)[1]){
		h = h + pi_pre[i] * (1 - pi_pre[i]) * (A[i, ] %*% t(A[i, ]))
	}
	return(h)
}

###########
weight <- function(x){
	return(x/sum(x))
}

loglike_v <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	l = sum(r * (0.5 * log(tau / s) - 0.5 * tau * beta^2/s ) ) + sum(-(a+1) * log(v) - b / v)
	return(l)
}

loglike_v_g <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	g = t(W) %*% (r * (-0.5/s + 0.5*tau*(beta/s)^2)) + (-(a+1) / v + b/v^2)
	return(g)
}

FisherInfo_v <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	h = diag(0, length(v))

	for(i in which(r>0)){
		s2 = s[i]^2
		s3 = s2 * s[i]
		h = h + r[i] * (0.5 / s2 - tau * beta[i]^2/s3) * ( W[i, ] %*% t(W[i, ]) )
	}
	# h = t(W) %*% diag(as.vector(r * (0.5 / s^2 - tau * beta^2/s^3))) %*% W 
		#+ diag(as.vector((a+1)/v^2 - 2*b/v^3))
	return(-h)
}

######################## log prior functions
logprior_sigma <- function(a, b, x){ return(-(1+a) * log(x) - b/x) }
logprior_pi <- function(a, b, x){ return((a-1) * log(x) + (b-1) * log(1 - x)) }
logprior_lambda <- function(s, x){ return(- 0.5 * x^2 / s) }
