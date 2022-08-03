/*
	Bayesian Functional GWAS with Summary Statistics --- MCMC (BFGWAS_SS:MCMC)
    Copyright (C) 2018  Jingjing Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "calcSS.h"


void CALCSS::CopyFromParam (PARAM &cPar)
{

    zipSS=cPar.zipSS;
    file_out=cPar.file_out;

    ni_total=cPar.ni_total;
    ns_total=cPar.ns_total;
    ni_test=cPar.ni_test;
    ns_test=cPar.ns_test;
    n_type = cPar.n_type;

    LDwindow=cPar.LDwindow;

    indicator_idv=cPar.indicator_idv;
    indicator_snp=cPar.indicator_snp;

    SNPmean = cPar.SNPmean;
    pheno_mean = cPar.pheno_mean;
    pheno_var = cPar.pheno_var;
    snp_pos = cPar.snp_pos;

    return;
}


//calculat summary statistics of score statistics and LD matrix
// Assume genotype X has SNPs in rows and samples in columns
void CALCSS::GetSS(gsl_matrix *X, gsl_vector *y, vector< vector<double> > &LD, vector<double> &beta, vector<double> &Z_SCORE, vector<double> &pval, vector<pair<size_t, double> > &pos_ChisqTest)
{
    cout << "\nStart calculating summary statistics ... \n";
    // X and y are standardized
    cout << "ni_test = " << ni_test << endl;
    cout << "ns_test = " << ns_test << endl;
    double yty;
    gsl_blas_ddot(y, y, &yty);
    pheno_var = yty / ((double)(ni_test-1)) ;
    cout << "Standardized pheno_var = " << pheno_var << "\n";

    // define used variables 
    gsl_vector *xvec_i = gsl_vector_alloc(ni_test);
    gsl_vector *xvec_j = gsl_vector_alloc(ni_test);

    double xtx_ij, xty,  beta_i, chisq_i, r2, z;

    // cout << "calculate beta, score statistics by the order of chr/bp ... \n";
    beta.clear();
    LD.clear();
    pval.clear();
    pos_ChisqTest.clear();

    if(isnan(pheno_var)==1 || pheno_var == 0)
    {
        cout << "Phenotype variance = 0. Only save LDcorr file...\n";
        for (size_t i=0; i<ns_test; ++i) {
        //calculate xtx_i
        //Lei's change
            gsl_matrix_get_row(xvec_i, X, snp_pos[i].pos);
            // saving X'X to LD
            LD.push_back(vector<double>()); // save correlation
            LD[i].push_back(1.0); // save correlation 1.0 for the diagnal values
            if(i < (ns_test-1) ){
                //calculate xtx_ij
                for(size_t j=(i+1); j < ns_test; ++j){
                    if( (snp_pos[j].chr == snp_pos[i].chr) && (snp_pos[j].bp <= snp_pos[i].bp + LDwindow) )
                    {
                        //Lei's change
                        gsl_matrix_get_row(xvec_j, X, snp_pos[j].pos);
                        gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
                        r2 = xtx_ij / ((double)ni_test) ;
                        LD[i].push_back( r2 ); // Correlation between x_i and x_j
                    }
                    else{break;}
                }
            }
        }
    }
    else{
        for (size_t i=0; i<ns_test; ++i) {
            //Lei's change
            gsl_matrix_get_row(xvec_i, X, snp_pos[i].pos);
            //calculate effect-size
            if(xvec_i->size != y->size){cerr << "Genotype length dose not equal to phenotype length!\n Some samples in the genotype file may not have genotype data!\n Please check your phenotype and genotype input files!\n"; exit(-1);}
            gsl_blas_ddot(xvec_i, y, &xty);
            beta_i = xty / ((double)ni_test);
            beta.push_back(beta_i); // effect size

            //z-score
            z = beta_i * sqrt((double)ni_test);
            Z_SCORE.push_back(z);

            // chisq_i = ((double)ni_test)*(log(yty)-log(yty-xty*xty/xtx_i)); // LRT statistic
            chisq_i = z * z; // Score test statistic
            pval.push_back( gsl_cdf_chisq_Q (chisq_i, 1.0) ); // pvalue needed for BVSRM
            pos_ChisqTest.push_back( make_pair(i, chisq_i) ) ; // pos_ChisqTest needed for BVSRM

            // saving X'X to LD
            LD.push_back(vector<double>()); // save correlation
            LD[i].push_back(1.0); // diagnal correlation

            if(i < (ns_test-1) ){
                //calculate xtx_ij
                for(size_t j=(i+1); j < ns_test; ++j){
                    if( (snp_pos[j].chr == snp_pos[i].chr) && (snp_pos[j].bp <= snp_pos[i].bp + LDwindow) )
                    {
                        //Lei's change
                        gsl_matrix_get_row(xvec_j, X, snp_pos[j].pos);
                        gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
                        r2 = xtx_ij / ((double)ni_test) ;
                        LD[i].push_back( r2 ); // Correlation between x_i and x_j
                    }
                    else{break;}
                }
            }
        }
    }
    gsl_vector_free(xvec_i);
    gsl_vector_free(xvec_j);
    return;
}

double Conv_xtx2_r2(const double &xtx_ij, const vector<double> &xtx_vec, const size_t &i, const size_t &j){
    // get correlation as r2
    double r2 = 0.0;
    if(xtx_ij != 0.0){
        if( (xtx_vec[i] > 0.0) && (xtx_vec[j] > 0.0) ){
            r2 = xtx_ij / sqrt( xtx_vec[i] * xtx_vec[j] );
        }
    }
    return r2;
}

// JY updated 06/15/2022
void CALCSS::WriteSS(const vector< vector<double> > &LD, const vector<double> &beta, const vector<double> &Z_SCORE,const vector<double> &pval)
{
    cout << "\nStart writing summary statistics ... \n";
    String fout = file_out.c_str();

    // output files matches RareMetalWorker outputs
    String cov_file_str = "./output/" + fout;
    String score_file_str = "./output/" + fout;

    IFILE cov_out=NULL;
    IFILE score_out=NULL;
    if( (isnan(pheno_var) != 1) && (pheno_var > 0)) {
        cout << "\nThe pheno variance > 0, writing z_score and LDcorr files. \n";
        if(zipSS){
            cov_file_str +=".LDcorr.txt.gz";
            cov_out = ifopen(cov_file_str, "w", InputFile::BGZF);
            score_file_str += ".Zscore.txt.gz";
            score_out = ifopen(score_file_str, "w", InputFile::BGZF);
            if(cov_out == NULL || score_out == NULL){
                perror("Fail to open LD or Zscore file!!! \n");
            }
        }
        else{
            cov_file_str +=".LDcorr.txt";
            cov_out = ifopen(cov_file_str, "w", InputFile::UNCOMPRESSED);
            score_file_str += ".Zscore.txt";
            score_out = ifopen(score_file_str, "w", InputFile::UNCOMPRESSED);
            if(cov_out == NULL || score_out == NULL){
                perror("Fail to open LD or Zscore file!!! \n");
            }
        }

        ifprintf(score_out, "#CHROM\tPOS\tID\tREF\tALT\tN\tMAF\tZ_SCORE\tmBeta\tPVALUE\n");
        // assuming variants have unique CHR:POS
        ifprintf(cov_out, "#ORDER\tCHROM\tPOS\tID\tREF\tALT\tN\tMAF\tCORR\n");
                          //order
        cout << "\nStart writting LDcorr and Zscore files ... \n";
        //Write files by the order of chr/bp in snp_pos
        for(size_t i=0; i<ns_test; i++){
            ifprintf(score_out, "%s\t%ld\t%s\t%s\t%s\t%u\t%.3e\t%.3e\t%.3e\t%.3e\n", snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf, Z_SCORE[i], beta[i], pval[i]);
            ifprintf(cov_out, "%u\t%s\t%ld\t%s\t%s\t%s\t%u\t%.3e\t", i,snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf);
            for(size_t j=0; j<LD[i].size(); j++){
                ifprintf(cov_out, "%.3e,", LD[i][j]);
            }
            ifprintf(cov_out, "\n");
        }
        ifclose(cov_out);
        ifclose(score_out);

        // tabix zipped files
        // printf("Lei's note: Cannot tabix");
        String cmd;
        int sys_status=1;
        if(zipSS){
            printf("Tabixing .LDcorr.txt.gz files ... \n");
            cmd = String("tabix -c \"#\" -s 2 -b 3 -e 3 -f ") + cov_file_str;
            sys_status = system(cmd.c_str());
            if ( sys_status == 0 ) {
                printf( "LD correlation file %s has been tabixed\n", cov_file_str.c_str() );
            }
            else {
                printf("Unable to tabix %s\n", cov_file_str.c_str());
            }
            printf("Tabixing .Zscore.txt.gz files ... \n");
            cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + score_file_str;
            sys_status = system(cmd.c_str());
            if ( sys_status == 0 ) {
                printf( "ZScore statistic file %s has been tabixed\n", score_file_str.c_str() );
            }
            else {
                printf("Unable to tabix %s\n", score_file_str.c_str());
            }
        }
    }
    else{
        cout << "\nPhenotype variance = 0. Only write LDcorr file. \n";
        if(zipSS){
            cov_file_str +=".LDcorr.txt.gz";
            cov_out = ifopen(cov_file_str, "w", InputFile::BGZF);

            if(cov_out == NULL){
                perror("Fail to open LD file!!! \n");
        }
        }
        else{
            cov_file_str +=".LDcorr.txt";
            cov_out = ifopen(cov_file_str, "w", InputFile::UNCOMPRESSED);

            if(cov_out == NULL ){
                perror("Fail to open LD file!!! \n");
            }
        }
        // assuming variants have unique CHR:POS
        ifprintf(cov_out, "#ORDER\tCHROM\tPOS\tID\tREF\tALT\tN\tMAF\tCORR\n");
        //Write files by the order of chr/bp in snp_pos
        for(size_t i=0; i<ns_test; i++){
            ifprintf(cov_out, "%u\t%s\t%ld\t%s\t%s\t%s\t%u\t%.3e\t", i,snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf);
            for(size_t j=0; j<LD[i].size(); j++){
                ifprintf(cov_out, "%.3e,", LD[i][j]);
            }
            ifprintf(cov_out, "\n");
        }
        ifclose(cov_out);

        // tabix zipped files
        String cmd;
        int sys_status=1;

        if(zipSS){
            printf("Tabixing .LDcorr.txt.gz files ... \n");
            cmd = String("tabix -c \"#\" -s 2 -b 3 -e 3 -f ") + cov_file_str;
            sys_status = system(cmd.c_str());
            if ( sys_status == 0 ) {
                printf( "LD correlation output %s has been tabixed\n", cov_file_str.c_str() );
            }
            else {
                printf("Unable to tabix %s\n", cov_file_str.c_str());
            }
        }
    }

    return;
}

void getXty(const vector<double> &beta, const vector<double> &xtx, vector <double> &Xty)
{
	//n is the sample size
    cout << "Calculate Xty ... \n";
    Xty.clear();
    for(size_t i=0; i<beta.size(); i++){
        Xty.push_back( beta[i] * xtx[i] );
    }
    return;
}

void getPval(const vector<double> &beta, const vector<double> &beta_sd, vector <double> &pval, vector<pair<size_t, double> > &pos_ChisqTest)
{
    cout << "Calculate pval ... \n";
    pval.clear();
    pos_ChisqTest.clear();
    double pval_i, chisq_i;

    for(size_t i=0; i<beta.size(); i++){
        chisq_i = pow(beta[i] / beta_sd[i], 2);
        pos_ChisqTest.push_back( make_pair(i, chisq_i) );

        pval_i = gsl_cdf_chisq_Q (chisq_i, 1);
        pval.push_back(pval_i);
    }
    return;
}

// (Updated 06/15/2022 JY)
// get correlation between pos_i and pos_j from reference
double getXtX(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j)
{
    double xtx_ij = 0.0;

    if(pos_i == pos_j){
        xtx_ij = 1.0;
    }
    else
    {
        if( (pos_j - pos_i) > 0.0 && (pos_j - pos_i) < LD[pos_i].size()  )
            {
                xtx_ij = LD[pos_i][pos_j - pos_i];
            }     
        else if( (pos_i - pos_j) > 0.0 && (pos_i - pos_j) < LD[pos_j].size() )
            {
                xtx_ij = LD[pos_j][pos_i - pos_j];
            }
    }

    return xtx_ij;
}


// (Updated 06/15/2022 JY)
double CalcResVar(const gsl_matrix *D_cond, const gsl_vector * beta_cond)
{
    size_t s_size = D_cond->size1;
    gsl_vector *D_beta = gsl_vector_alloc (s_size);
    gsl_blas_dgemv(CblasNoTrans, 1, D_cond, beta_cond, 0, D_beta);
    double R2;
    gsl_blas_ddot (D_beta, beta_cond, &R2);
    if( R2 > 1.0){
        R2 = 1.0 ;
    }else if(R2 < 0.0){
        R2 = 0.0 ;
    }
    return (1.0 - R2);
}

// (Updated 06/15/2022 JY)
// Calculate beta_hat based on multivariable model from mbeta
void CalcBeta(const gsl_matrix *D_cond, const gsl_vector * mbeta_cond, gsl_vector * beta_cond)
{
    size_t s_size = mbeta_cond->size;
    gsl_matrix *Dinv = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(Dinv, D_cond);
    if(LapackSolve(Dinv, mbeta_cond, beta_cond)!=0)
       EigenSolve(Dinv, mbeta_cond, beta_cond);
    gsl_matrix_free(Dinv);
    return ;
}












