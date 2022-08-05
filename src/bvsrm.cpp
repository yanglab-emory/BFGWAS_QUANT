/*
    Bayesian Functional GWAS --- MCMC (BFGWAS_QUANT:MCMC)
    Copyright (C) 2022  Jingjing Yang

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
#include "bvsrm.h"

// define to_string function : convert to string
// In C++, a template is a way to introduce a concept known as "generics".
// With generics, you no longer need to concern your self with creating a function for each type,
// the type will be deduced from the template function signature.
template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

void BVSRM::CopyFromParam (PARAM &cPar)
{
    VcfSampleID_test = cPar.VcfSampleID_test;
    Anum = cPar.Anum;
    e = cPar.e;
    vscale = cPar.vscale;
    FIXHYP = cPar.FIXHYP;
    saveSS = cPar.saveSS;
    LDwindow = cPar.LDwindow;

    iniType = cPar.iniType;
    iniSNPfile = cPar.iniSNPfile;
    hypfile = cPar.hypfile;
    SNPmean = cPar.SNPmean;

    win = cPar.win;
    nadd_accept = cPar.nadd_accept;
    ndel_accept= cPar.ndel_accept;
    nswitch_accept= cPar.nswitch_accept;
    nother_accept= cPar.nother_accept;
    nadd= cPar.nadd;
    ndel= cPar.ndel;
    nswitch= cPar.nswitch;
    nother= cPar.nother;

      a_mode=cPar.a_mode;
      d_pace=cPar.d_pace;

    file_bfile=cPar.file_bfile;
    file_geno=cPar.file_geno;
    file_vcf = cPar.file_vcf;
    file_out=cPar.file_out;
    file_anno=cPar.file_anno;

      l_min=cPar.h_min;
      l_max=cPar.h_max;
      pheno_mean=cPar.pheno_mean;
    pheno_var=cPar.pheno_var;

      time_UtZ=0.0;
      time_Omega=0.0;
      n_accept=0;
    region_pip = 0.0;
    Switch_Flag = 0;

      h_min=cPar.h_min;
      h_max=cPar.h_max;
      h_scale=cPar.h_scale;
      rho_min=cPar.rho_min;
      rho_max=cPar.rho_max;
      rho_scale=cPar.rho_scale;
      logp_min=cPar.logp_min;
      logp_max=cPar.logp_max;
      logp_scale=cPar.logp_scale;

      s_min=cPar.s_min;
      s_max=cPar.s_max;
      w_step=cPar.w_step;
      s_step=cPar.s_step;
      n_mh=cPar.n_mh;
      randseed=cPar.randseed;
      trace_G=cPar.trace_G;

      ni_total=cPar.ni_total;
      ns_total=cPar.ns_total;
      ni_test=cPar.ni_test;
      ns_test=cPar.ns_test;

      indicator_idv=cPar.indicator_idv;
      indicator_snp=cPar.indicator_snp;
      snpInfo=cPar.snpInfo;
    snp_pos = cPar.snp_pos;

    // summary stat related
    refLD = cPar.refLD;
    mbeta = cPar.mbeta;
    //mbeta_SE = cPar.mbeta_SE;
    pval_vec = cPar.pval_vec;
    pos_ChisqTest = cPar.pos_ChisqTest;
	//xtx_vec = cPar.xtx_vec;
    //snp_var_vec = cPar.snp_var_vec;
    ni_effect_vec = cPar.ni_effect_vec;
    return;
}


void BVSRM::CopyToParam (PARAM &cPar)
{
    cPar.time_Omega=time_Omega;
    cPar.time_Proposal=time_Proposal;
    cPar.cHyp_initial=cHyp_initial;
    cPar.n_accept=n_accept;
    cPar.pheno_mean=pheno_mean;
    cPar.pheno_var = pheno_var;
    cPar.randseed=randseed;

    cPar.nadd_accept=nadd_accept;
    cPar.ndel_accept=ndel_accept;
    cPar.nswitch_accept=nswitch_accept;
    cPar.nother_accept=nother_accept;
    cPar.nadd=nadd;
    cPar.ndel=ndel;
    cPar.nswitch=nswitch;
    cPar.nother=nother;
    cPar.region_pip = region_pip;

    return;
}

bool comp_vec (size_t a, size_t b)
{
    return (a < b);
}

bool comp_pi (pair<string, double> a, pair<string, double> b)
{
    return (a.second > b.second);
}

//JY add to write initial significant SNP id out
void BVSRM::WriteIniRank (const vector<string> &iniRank)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".inirank.txt";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    for (size_t i=0; i<iniRank.size(); ++i) {

            outfile<<iniRank[i]<<endl;

    }

    outfile.clear();
    outfile.close();
    return;
}

void BVSRM::WriteIniSNP (const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".iniSNP";
    //cout << "write iniSNP at " << file_str << endl;

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    for (size_t i=0; i<rank.size(); ++i) {
        outfile<< snp_pos[SNPrank_vec[rank[i]].second].rs <<endl;
    }

    outfile.clear();
    outfile.close();
    return;
}

void BVSRM::WriteMCMC(const vector<string> &snps_mcmc){
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".mcmc";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    for (size_t i=0; i<snps_mcmc.size(); ++i) {
        outfile<< snps_mcmc[i] <<endl;
    }

    outfile.clear();
    outfile.close();
    return;
}


void BVSRM::WriteParam(vector<pair<double, double> > &beta_g, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr, const vector<double> &Z_scores, const vector<double> pval_lrt)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".paramtemp";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    //outfile<<"chr"<<"\t" <<"bp"<<"\t" << "markerID"<<"\t"<<"REF"<<"\t" <<"ALT"<<"\t" << "maf" << "\t" <<"gamma" << "\t" <<"beta_bfgwas"<<"\t"<< "mbeta_SE" << "\t" << "LRT" << "\t" << "pval_lrt"  << "\t" << "rank" << endl;

    size_t pos, r;
    vector< pair<string , double> > pi_vec;
    string rs;
    double pi_temp;

    for (size_t i=0; i<snp_pos.size(); ++i) {

        // save the data along the order of all variants, snp_pos is sorted by order
        rs = snp_pos[i].rs;
        outfile<<snp_pos[i].chr<<"\t"<<snp_pos[i].bp<<"\t"<<rs<<"\t"<< snp_pos[i].a_major<<"\t"<<snp_pos[i].a_minor<<"\t" ;
        outfile << scientific << setprecision(3)  << snp_pos[i].maf << "\t";
        

        //Lei's note (7/10)
        //for (size_t j=0; j < n_type; j++) {
            //if (snp_pos[i].indicator_func[j]) {
                //outfile << j << "\t";
                //break;
            //}
            //else if(j == (n_type - 1)) outfile << "NA" << "\t";
        //}

        //Or just do not use?
        
        pos = snp_pos[i].pos;
        //beta_g is saved by position
        if (beta_g[pos].second!=0) {
            pi_temp = beta_g[pos].second/(double)s_step;
            outfile << pi_temp <<"\t" << beta_g[pos].first/beta_g[pos].second<< "\t" << mbeta_SE[pos] << "\t" ;
        }
        else {
            pi_temp = 0.0;
            outfile << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t";
        }
        pi_vec.push_back(make_pair(rs, pi_temp));

        if ( SNPorder_vec[i].first  != pos) {
            cerr << "ERROR: SNPorder_vec[i].first not equal to pos... \n";
            exit(-1);
        }
        r = mapOrder2Rank[i]; //map to rank
        outfile << scientific << setprecision(3) << pos_loglr[r].second << "\t"<< pval_lrt[pos] << "\t" ;
        outfile << r << endl;
    }
    outfile.clear();
    outfile.close();

    return;
}


void BVSRM::WriteParam_SS(vector<pair<double, double> > &beta_g, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_ChisqTest, const vector<double> pval, const vector< vector<double> > &LD)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".paramtemp";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    outfile<<"#CHR\tPOS\tID\tREF\tALT\tMAF\tPi\tBeta\tmBeta\tChisqTest\tPval_svt\tRank\t";
    for(size_t j = 0; j < Anum; ++j)
    {
        if(j < (Anum - 1)){
            outfile<<"Anno_" << j << "\t";
        }
        else{
            outfile<<"Anno_" << j << endl;
        }
    }

    size_t r, n_causal;
    vector<size_t> rank_vec; // Save positions of potential causal SNPs
    map<size_t, size_t> mapRank2Vec; // map Rank to position in rank_vec

    double pi_temp;
    beta_mcmc.assign(ns_test, 0.0); // Bayesian estimates of beta
    em_gamma = 0.0; //save sum of gamma;
    sumbeta2 = 0.0; // save sum of gamma * beta2

    n_causal = 0;
    for (size_t i=0; i<ns_test; ++i) {
        if ( (beta_g[i].second / (double)s_step) >= 0.01)
        {
            r = mapPos2Rank[i];
            rank_vec.push_back(r) ;
            mapRank2Vec[r]=n_causal;
            n_causal++;
        }
    }
    size_t r_size = rank_vec.size();
    cout << "WriteParam_SS: Selected # SNPs with PIP > 0.01 : " << r_size << endl;

    size_t pos_i, pos_j;
    double xtx_ij;
    if(r_size > 1) {
        gsl_matrix *D_gamma = gsl_matrix_alloc(r_size, r_size);
        gsl_vector *mbeta_gamma = gsl_vector_alloc(r_size);
        gsl_vector *beta_hat = gsl_vector_alloc(r_size); // calculate posterior beta
        for( size_t ii=0; ii < r_size; ++ii){
            pos_i = mapRank2pos[rank_vec[ii]];
            gsl_vector_set(mbeta_gamma, ii, mbeta[pos_i]);
            gsl_matrix_set(D_gamma, ii, ii, 1);
            for(size_t jj=(ii+1); jj < (r_size); ++jj ){
                pos_j = mapRank2pos[rank_vec[jj]];
                xtx_ij = getXtX(LD, pos_i, pos_j);
                gsl_matrix_set(D_gamma, ii, jj, xtx_ij);
                gsl_matrix_set(D_gamma, jj, ii, xtx_ij);
            }
        }
        gsl_vector_view D_diag = gsl_matrix_diagonal(D_gamma);
        gsl_vector_add_constant(&D_diag.vector, tau_beta);
        // posterior Bayesian estimates of beta
        if(LapackSolve(D_gamma, mbeta_gamma, beta_hat)!=0)
            EigenSolve(D_gamma, mbeta_gamma, beta_hat);

        for (size_t i=0; i<ns_test; ++i) {
            r = mapPos2Rank[i]; //map to rank
            outfile<<snp_pos[i].chr<<"\t"<<snp_pos[i].bp<<"\t"<<snp_pos[i].rs<<"\t"<< snp_pos[i].a_major<<"\t"<<snp_pos[i].a_minor<<"\t" ;
            if(snp_pos[i].maf < 0.0 || snp_pos[i].maf > 1.0){
                outfile << scientific << "NA" << "\t";
            }else{
                outfile << scientific << setprecision(3)  << snp_pos[i].maf << "\t";
            }

            if ( beta_g[i].second > 0.0 ) {
                pi_temp = beta_g[i].second/(double)s_step;
                if(mapRank2Vec[r]){
                    beta_mcmc[i] = gsl_vector_get(beta_hat, mapRank2Vec[r]);
                }else{
                    beta_mcmc[i] = beta_g[i].first/beta_g[i].second;
                }
                outfile << pi_temp << "\t" << beta_mcmc[i] << "\t" << mbeta[i]  << "\t" ;
                em_gamma = em_gamma + pi_temp ;
                sumbeta2 = sumbeta2 + pi_temp * beta_mcmc[i] * beta_mcmc[i] ;
            }
            else{
                pi_temp = 0.0;
                outfile << 0.0 << "\t" << 0.0 << "\t" << mbeta[i] << "\t";
            }
            outfile << scientific << setprecision(3) << pos_ChisqTest[r].second << "\t"<< pval[r] << "\t" ;
            outfile << r << "\t";

            for(size_t j = 0; j < snp_pos[i].annoscore.size(); ++j)
            {
                if(j < (snp_pos[i].annoscore.size() -1))
                {
                    outfile << snp_pos[i].annoscore[j] << "\t" ;
                }
                else{
                    outfile << snp_pos[i].annoscore[j] << endl ;
                }
            }
        }
        gsl_matrix_free(D_gamma);
        gsl_vector_free(mbeta_gamma);
        gsl_vector_free(beta_hat);
    }
    else{
        for (size_t i=0; i<ns_test; ++i) {
            r = mapPos2Rank[i]; //map to rank
            outfile<<snp_pos[i].chr<<"\t"<<snp_pos[i].bp<<"\t"<<snp_pos[i].rs<<"\t"<< snp_pos[i].a_major<<"\t"<<snp_pos[i].a_minor<<"\t" ;
            if(snp_pos[i].maf < 0.0 || snp_pos[i].maf > 1.0){
                outfile << scientific << "NA" << "\t";
            }else{
                outfile << scientific << setprecision(3)  << snp_pos[i].maf << "\t";
            }

            if ( beta_g[i].second > 0.0 ) {
                pi_temp = beta_g[i].second/(double)s_step;
                beta_mcmc[i] = beta_g[i].first/beta_g[i].second;
                outfile << pi_temp << "\t" << beta_mcmc[i] << "\t" << mbeta[i]  << "\t" ;
                em_gamma = em_gamma + pi_temp ;
                sumbeta2 = sumbeta2 + pi_temp * beta_mcmc[i] * beta_mcmc[i] ;
            }else{
                pi_temp = 0.0;
                outfile << 0.0 << "\t" << 0.0 << "\t" << mbeta[i] << "\t";
            }
            outfile << scientific << setprecision(3) << pos_ChisqTest[r].second << "\t"<< pval[r] << "\t" ;
            outfile << r << "\t";

            for(size_t j = 0; j < snp_pos[i].annoscore.size(); ++j)
            {
                if(j < (snp_pos[i].annoscore.size() -1))
                {
                    outfile << snp_pos[i].annoscore[j] << "\t" ;
                }
                else{
                    outfile << snp_pos[i].annoscore[j] << endl ;
                }
            }
        }
    }
    outfile.clear();
    outfile.close();

    return;
}


//writing the annotation score file for M steps
void BVSRM::WriteAnnoScore(vector<pair<double, double> > &beta_g,  const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".AnnoScore.txt";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    size_t pos;
    string rs;
    double pi_temp;

    for (size_t i=0; i<snp_pos.size(); ++i) {

        rs = snp_pos[i].rs;
        outfile<<rs<<"\t";

        pos = snp_pos[i].pos;
        //beta_g is saved by position
        if (beta_g[pos].second!=0) {
            pi_temp = beta_g[pos].second/(double)s_step;
            outfile << pi_temp <<"\t" << beta_g[pos].first/beta_g[pos].second ;
        }
        else {
            pi_temp = 0.0;
            outfile << 0.0 << "\t" << 0.0 ;
        }

        for(size_t j=0; j < Anum; j++){
              outfile << "\t" << snp_pos[i].annoscore[j];
        }

        outfile<<endl;

    }
    outfile.clear();
    outfile.close();

    return;
}


void BVSRM::WriteFGWAS_InputFile(const vector<SNPPOS> &snp_pos, const vector<double> &Z_scores)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".ifgwas";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    //header of the space-delimited text file
    outfile<<"SNPID"<<" "<<"CHR"<<" " <<"POS"<<" " << "F" << " " << "Z"<< " " << "SE" << " " << "N" << " " << "nonsynonymous" << endl;
    
    size_t pos;
    
    for (size_t i=0; i<ns_test; ++i) {
        // save the data along the order of all variants, snp_pos is sorted by order
        pos = snp_pos[i].pos;

        outfile<< snp_pos[i].rs <<" "<< snp_pos[i].chr<<" " <<snp_pos[i].bp << " ";
        
        outfile << scientific << setprecision(6)  << snp_pos[i].maf << " " << Z_scores[pos] << " " << mbeta_SE[pos] << " ";

        outfile << ni_test << " ";

        //if(snp_pos[i].indicator_func[0]) outfile << 1 ;
        //else outfile << 0;
        
        outfile << endl;
    }

    outfile.clear();
    outfile.close();
}

// lei's change
void BVSRM::WriteGenotypeFile(gsl_matrix *X, const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".geno";

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    //write header with VcfSampleID_test
    outfile<<"ID"<<"\t"<<"CHROM"<<"\t" <<"POS"<<"\t" << "REF"<< "\t" << "ALT"  << "\t";
    for (size_t i=0; i<ni_test; i++) {
        if (i ==(ni_test-1)) {
            outfile << VcfSampleID_test[i] << endl;
        }
        else outfile << VcfSampleID_test[i] << "\t";
    }

    size_t pos;
    string rs;
    double geno_j;
    
    for (size_t i=0; i<ns_test; ++i) {
        
        pos = snp_pos[i].pos;
        // save the data along the order of all variants, snp_pos is sorted by order
        rs = snp_pos[i].rs;
        outfile<< rs <<"\t"<< snp_pos[i].chr<<"\t" <<snp_pos[i].bp << "\t";
        outfile << snp_pos[i].a_major << "\t" <<snp_pos[i].a_minor << "\t";

        /*for (size_t j=0; j < n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                outfile << j << "\t";
                break;
            }
            else if(j == (n_type - 1)) outfile << "NA" << "\t";
        }*/
        //outfile << scientific << setprecision(6)  << snp_pos[i].maf << "\t";

        for (size_t j=0; j < ni_test; j++) {
            //c = X[pos][j];
            //geno_j = UcharTable[(int)c].second;
            geno_j = gsl_matrix_get(X,pos,j);
            if(geno_j < 0.0 || geno_j > 2.0){
                cout << "ERROR: genotype = " << geno_j << endl;
                exit(-1);
            }else{
                if (geno_j == 0.0) geno_j = 0;
                else if (geno_j == 2.0) geno_j = 2;
                else if (geno_j == 1.0) geno_j = 1;

                if (j == (ni_test-1))
                    outfile << fixed << setprecision(2)  << geno_j << endl;
                else
                    outfile << fixed << setprecision(2) << geno_j << "\t";
            }
        }
    }

    outfile.clear();
    outfile.close();
}


void BVSRM::SetPgamma (size_t p_gamma_top)
{
    double p, q;

    if(p_gamma_top < 50){
        p_gamma_top=100;
    } else if(p_gamma_top > 300){
        p_gamma_top=300;
    }

  if((ns_test-p_gamma_top) < 0){
        p = 1.0 / double(ns_test);
        for (size_t i=0; i<ns_test; ++i) {
            p_gamma[i] = p;
        }
  } else{
    p = 0.9 / double(p_gamma_top);
    q = 0.1 / ((double)(ns_test-p_gamma_top));
    for (size_t i=0; i<ns_test; ++i) {
        if(i < p_gamma_top) p_gamma[i] = p;
            else p_gamma[i] = q;
    }
  }
    return;
}




//currently used
//lei's change
void BVSRM::SetXgamma (gsl_matrix *Xgamma, gsl_matrix *X, vector<size_t> &rank)
{
    size_t pos;
    for (size_t i=0; i<rank.size(); ++i) {
        pos=SNPrank_vec[rank[i]].first;
        gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
        //copy the column
        gsl_matrix_get_col(&Xgamma_col.vector,X,pos);

        //getGTgslVec(X, &Xgamma_col.vector, pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
    }
    return;
}

//currently used
//lei's change
void BVSRM::SetXgamma (gsl_matrix *X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, const vector<size_t> &rank_new, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    double d;
    // cout << "X_add set start" << endl;

    //rank_old and rank_new are sorted already inside PorposeGamma
    //calculate vectors rank_remove and rank_add
    //  size_t v_size=max(rank_old.size(), rank_new.size());
    //make sure that v_size is larger than repeat
    size_t v_size=20;
    vector<size_t> rank_remove(v_size), rank_add(v_size), rank_union(s_max+v_size);
    vector<size_t>::iterator it;

    it=set_difference (rank_old.begin(), rank_old.end(), rank_new.begin(), rank_new.end(), rank_remove.begin());
    rank_remove.resize(it-rank_remove.begin());
    stable_sort(rank_remove.begin(), rank_remove.end(), comp_vec);
   // cout << "rank_remove: "; PrintVector(rank_remove);

    it=set_difference (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_add.begin());
    rank_add.resize(it-rank_add.begin());
    stable_sort(rank_add.begin(), rank_add.end(), comp_vec);
  //  cout << "rank_add: "; PrintVector(rank_add);

    it=set_union (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_union.begin());
    rank_union.resize(it-rank_union.begin());
    stable_sort(rank_union.begin(), rank_union.end(), comp_vec);
   // cout << "rank_union: "; PrintVector(rank_union);

    //map rank_remove and rank_add
    map<size_t, int> mapRank2in_remove, mapRank2in_add;
    for (size_t i=0; i<rank_remove.size(); i++) {
        mapRank2in_remove[rank_remove[i]]=1;
    }
    for (size_t i=0; i<rank_add.size(); i++) {
        mapRank2in_add[rank_add[i]]=1;
    }

    //obtain the subset of matrix/vector
    gsl_matrix_const_view Xold_sub=gsl_matrix_const_submatrix(X_old, 0, 0, X_old->size1, rank_old.size());
    gsl_matrix_const_view XtXold_sub=gsl_matrix_const_submatrix(XtX_old, 0, 0, rank_old.size(), rank_old.size());
    gsl_vector_const_view Xtyold_sub=gsl_vector_const_subvector(Xty_old, 0, rank_old.size());

    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(X_new, 0, 0, X_new->size1, rank_new.size());
    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
    gsl_vector_view Xtynew_sub=gsl_vector_subvector(Xty_new, 0, rank_new.size());

    if (rank_remove.size()==0 && rank_add.size()==0) {
        gsl_matrix_memcpy(&Xnew_sub.matrix, &Xold_sub.matrix);
        gsl_matrix_memcpy(&XtXnew_sub.matrix, &XtXold_sub.matrix);
        gsl_vector_memcpy(&Xtynew_sub.vector, &Xtyold_sub.vector);
        //cout << "rank_old = rank_new; " << "Xgamma_new set success" << endl;
    } else {
        size_t i_old, j_old, i_new, j_new, i_add, j_add, i_flag, j_flag;
        if (rank_add.size()==0) {
            //only delete a snp
            i_old=0; i_new=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}

                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);

                d=gsl_vector_get (Xty_old, i_old);
                gsl_vector_set (Xty_new, i_new, d);

                j_old=i_old; j_new=i_new;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}

                    d=gsl_matrix_get(XtX_old, i_old, j_old);

                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}

                    j_old++; j_new++;
                }
                i_old++; i_new++;
            }
            //cout << "X_add = NULL; " << "Xgamma_new set success" << endl;
        } else {
            //rank_add has length > 0
            gsl_matrix *X_add=gsl_matrix_alloc(X_old->size1, rank_add.size() );
            gsl_matrix *XtX_aa=gsl_matrix_alloc(X_add->size2, X_add->size2);
            gsl_matrix *XtX_ao=gsl_matrix_alloc(X_add->size2, rank_old.size());
            gsl_vector *Xty_add=gsl_vector_alloc(X_add->size2);

            //get X_add
            SetXgamma (X_add, X, rank_add);
            //cout << "X_add set success" << endl;

            //get t(X_add)X_add and t(X_add)X_temp
            clock_t time_start=clock();

            //somehow the lapack_dgemm does not work here
            //#ifdef WITH_LAPACK
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_add, 0.0, XtX_aa);
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_old, 0.0, XtX_ao);

            //#else
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, X_add, 0.0, XtX_aa);
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, &Xold_sub.matrix, 0.0, XtX_ao);
            //#endif
            gsl_blas_dgemv(CblasTrans, 1.0, X_add, y, 0.0, Xty_add);

            time_Omega += (clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

            //save to X_new, XtX_new and Xty_new
            i_old=0; i_new=0; i_add=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}
                if (mapRank2in_add.count(rank_new[i_new])!=0) {i_flag=1; //within x_add
                } else {i_flag=0; //within x_common
                }

                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                if (i_flag==1) {
                    gsl_vector_view Xcopy_col=gsl_matrix_column(X_add, i_add);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                } else {
                    gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                }
                //  cout << "Xgamma_new set success" << endl;

                if (i_flag==1) {
                    d=gsl_vector_get (Xty_add, i_add);
                } else {
                    d=gsl_vector_get (Xty_old, i_old);
                }
                gsl_vector_set (Xty_new, i_new, d);
                // cout << "Xty_new set success" << endl;

                j_old=i_old; j_new=i_new; j_add=i_add;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}
                    if (mapRank2in_add.count(rank_new[j_new])!=0) {j_flag=1;} else {j_flag=0;}

                    if (i_flag==1 && j_flag==1) {
                        d=gsl_matrix_get(XtX_aa, i_add, j_add);
                    } else if (i_flag==1 && j_flag==0) {
                        d=gsl_matrix_get(XtX_ao, i_add, j_old);
                    } else if (i_flag==0 && j_flag==1) {
                        d=gsl_matrix_get(XtX_ao, j_add, i_old);
                    } else {
                        d=gsl_matrix_get(XtX_old, i_old, j_old);
                    }

                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}

                    j_new++;
                    if (j_flag==1) {j_add++;} else {j_old++;}
                }
                //cout << "XtX_new success" << endl;
                i_new++; if (i_flag==1) {i_add++;} else {i_old++;}
            }
            // cout << "X_gamma set success" << endl;

            gsl_matrix_free(X_add);
            gsl_matrix_free(XtX_aa);
            gsl_matrix_free(XtX_ao);
            gsl_vector_free(Xty_add);
        }

    }

    rank_remove.clear();
    rank_add.clear();
    rank_union.clear();
    mapRank2in_remove.clear();
    mapRank2in_add.clear();

    return;
}


//currently used in propose gamma
void BVSRM::SetXgammaDel (const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const vector<size_t> &rank_old, size_t col_id, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    size_t s_size = rank_old.size();
    size_t s2;

    if (col_id==0) {
        s2 = s_size-1;
        gsl_matrix_const_view X2old = gsl_matrix_const_submatrix(X_old, 0, 1, ni_test, s2);
        gsl_matrix_const_view XtX22_sub = gsl_matrix_const_submatrix(XtX_old, 1, 1, s2, s2);
        gsl_vector_const_view Xty2_sub = gsl_vector_const_subvector(Xty_old, 1, s2);

        gsl_matrix_view Xnew2_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, s2);
        gsl_matrix_view XtXnew22_sub=gsl_matrix_submatrix(XtX_new, 0, 0, s2, s2);
        gsl_vector_view Xtynew2_sub=gsl_vector_subvector(Xty_new, 0, s2);

        gsl_matrix_memcpy(&Xnew2_sub.matrix, &X2old.matrix);
        gsl_matrix_memcpy(&XtXnew22_sub.matrix, &XtX22_sub.matrix);
        gsl_vector_memcpy(&Xtynew2_sub.vector, &Xty2_sub.vector);
    }
    else if(col_id == (s_size-1)){

        gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, col_id);
        gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, col_id, col_id);
        gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, col_id);

        gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, col_id);
        gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, col_id, col_id);
        gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, col_id);

        gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
        gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
        gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);
    }
    else{
        s2 = s_size - col_id - 1;
        gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, col_id);
        gsl_matrix_const_view X2old = gsl_matrix_const_submatrix(X_old, 0, col_id+1, ni_test, s2);

        gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, col_id, col_id);
        gsl_matrix_const_view XtX12_sub = gsl_matrix_const_submatrix(XtX_old, 0, col_id+1, col_id, s2);
        gsl_matrix_const_view XtX21_sub = gsl_matrix_const_submatrix(XtX_old, col_id+1, 0, s2, col_id);
        gsl_matrix_const_view XtX22_sub = gsl_matrix_const_submatrix(XtX_old, col_id+1, col_id+1, s2, s2);

        gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, col_id);
        gsl_vector_const_view Xty2_sub = gsl_vector_const_subvector(Xty_old, col_id+1, s2);

        gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, col_id);
        gsl_matrix_view Xnew2_sub=gsl_matrix_submatrix(X_new, 0, col_id, ni_test, s2);

        gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, col_id, col_id);
        gsl_matrix_view XtXnew12_sub=gsl_matrix_submatrix(XtX_new, 0, col_id, col_id, s2);
        gsl_matrix_view XtXnew21_sub=gsl_matrix_submatrix(XtX_new, col_id, 0, s2, col_id);
        gsl_matrix_view XtXnew22_sub=gsl_matrix_submatrix(XtX_new, col_id, col_id, s2, s2);

        gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, col_id);
        gsl_vector_view Xtynew2_sub=gsl_vector_subvector(Xty_new, col_id, s2);

        gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
        gsl_matrix_memcpy(&Xnew2_sub.matrix, &X2old.matrix);

        gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
        gsl_matrix_memcpy(&XtXnew12_sub.matrix, &XtX12_sub.matrix);
        gsl_matrix_memcpy(&XtXnew21_sub.matrix, &XtX21_sub.matrix);
        gsl_matrix_memcpy(&XtXnew22_sub.matrix, &XtX22_sub.matrix);

        gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);
        gsl_vector_memcpy(&Xtynew2_sub.vector, &Xty2_sub.vector);
    }

}

void BVSRM::SetXgammaAdd (gsl_matrix *X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, size_t ranki, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    double xty;
    size_t s_size = rank_old.size();
    size_t pos = SNPrank_vec[ranki].first;

    if (s_size==0) {
        cerr << "setXgammaAdd rank_old has size 0\n";
        exit(-1);
    }
    //copy rank_old
    gsl_matrix_const_view X1old = gsl_matrix_const_submatrix(X_old, 0, 0, ni_test, s_size);
    gsl_matrix_const_view XtX11_sub = gsl_matrix_const_submatrix(XtX_old, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty1_sub = gsl_vector_const_subvector(Xty_old, 0, s_size);

    gsl_matrix_view Xnew1_sub=gsl_matrix_submatrix(X_new, 0, 0, ni_test, s_size);
    gsl_matrix_view XtXnew11_sub=gsl_matrix_submatrix(XtX_new, 0, 0, s_size, s_size);
    gsl_vector_view Xtynew1_sub=gsl_vector_subvector(Xty_new, 0, s_size);

    gsl_matrix_memcpy(&Xnew1_sub.matrix, &X1old.matrix);
    gsl_matrix_memcpy(&XtXnew11_sub.matrix, &XtX11_sub.matrix);
    gsl_vector_memcpy(&Xtynew1_sub.vector, &Xty1_sub.vector);

    //create ranki
    gsl_vector_view xvec = gsl_matrix_column(X_new, s_size);
    //getGTgslVec(X, &xvec.vector, pos, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
    //lei's change
    //copy the column
    gsl_matrix_get_col(&xvec.vector,X,pos);


    gsl_vector_view Xtx_col = gsl_matrix_subcolumn(XtX_new, s_size, 0, s_size);
    gsl_vector_view Xtx_row = gsl_matrix_subrow(XtX_new, s_size, 0, s_size);
    gsl_blas_dgemv(CblasTrans, 1.0, &X1old.matrix, &xvec.vector, 0.0, &Xtx_col.vector);
    gsl_vector_memcpy(&Xtx_row.vector, &Xtx_col.vector);
    gsl_matrix_set(XtX_new, s_size, s_size, xtx_vec[pos]);

    gsl_blas_ddot(&xvec.vector, y, &xty);
    gsl_vector_set(Xty_new, s_size, xty);

}


//end of currently used

double BVSRM::CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2)
{
    double pve, var_y;

    gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
    gsl_vector *Xty=gsl_vector_alloc (UtXgamma->size2);
    gsl_vector *OiXty=gsl_vector_alloc (UtXgamma->size2);

    gsl_matrix_set_identity (Omega);
    gsl_matrix_scale (Omega, 1.0/sigma_a2);

#ifdef WITH_LAPACK
    lapack_dgemm ((char *)"T", (char *)"N", 1.0, UtXgamma, UtXgamma, 1.0, Omega);
#else
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, UtXgamma, UtXgamma, 1.0, Omega);
#endif
    gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma, Uty, 0.0, Xty);

    CholeskySolve(Omega, Xty, OiXty);

    gsl_blas_ddot (Xty, OiXty, &pve);
    gsl_blas_ddot (Uty, Uty, &var_y);

    pve/=var_y;

    gsl_matrix_free (Omega);
    gsl_vector_free (Xty);
    gsl_vector_free (OiXty);

    return pve;
}


void BVSRM::setHyp(double tau_beta_temp)
{

    // Default initial values
    //cout << "Anum = " << Anum <<endl;
    tau = 1.0 ;
    avector.assign((Anum+1), 0.0);
    avector[0] = -12.0;
    tau_beta = tau_beta_temp;

    //cout << "load fixed parameter a & b vector from : " << hypfile << endl;
    string line;
    char *pch, *nch;
    size_t anno_idx;
    double theta_value = 0.0;

    if (hypfile.empty()) {
        cout << "Did not specify hypefile, use default values ...\n";
    }
    else{
        ifstream infile(hypfile.c_str(), ifstream::in);
        if(!infile) {
            cout << "Error opening file " << hypfile << endl;
            exit(-1);
        }
        //cout << "load hyp from hypfile... " << hypfile << endl;
        while (!safeGetline(infile, line).eof()) {
            if (line[0] == '#') {
                continue;
            }
            else if (line[0] == 'a'){
                pch = (char *)line.c_str();
                nch = strchr(pch, '\t');
                pch = (nch == NULL) ? NULL : nch+1;
                avector[0] = strtod(pch, NULL);
                anno_idx = 1;
                while(pch != NULL){
                    nch = strchr(pch, ',');
                    pch = (nch == NULL) ? NULL : nch+1;
                    if (anno_idx <= Anum ){
                        avector[anno_idx] = strtod(pch, NULL);
                        anno_idx++;
                    }
                }
            }
            else if (line[0] == 't'){
                pch = (char *)line.c_str();
                nch = strchr(pch, '\t');
                pch = (nch == NULL) ? NULL : nch+1;
                tau_beta = strtod(pch, NULL);
            }
            else{ break; }
        }
        infile.clear();
        infile.close();
    }
    cout << "Prior avector : ";
    PrintVector(avector) ;
    cout << "Prior tau_beta = " << tau_beta << endl;

    theta.clear();
    // setup theta vector from enrichment coefficient avector
    //cout << "snp_pos.size() = " << snp_pos.size() << endl;
    //cout << "ns_test = " << ns_test << endl;
    for (size_t i=0; i < snp_pos.size(); ++i){
        theta_value = avector[0];
        if(snp_pos[i].annoscore.size() > 0){
            for (size_t j=0; j < snp_pos[i].annoscore.size(); ++j){
              theta_value = theta_value + avector[j+1] * snp_pos[i].annoscore[j];
            }
        }else{
            cout << snp_pos[i].key << endl;
            cout << "snp_pos[i].annoscore is empty ... \n" ;
        }
        theta_value = exp(theta_value) / (1.0 + exp(theta_value));
        theta.push_back(theta_value);
        /*if(i < 10){
            cout << snp_pos[i].key << endl;
            cout << "snp_pos[i].annoscore : " ;
            PrintVector(snp_pos[i].annoscore);
            cout << "Prior CPP theta = " << theta[i] << endl;
        }*/
    }

    sum_theta = 0.0 ;
    // cout << "\n theta.size() = " << theta.size() << endl;
    for (size_t i=0; i < theta.size(); ++i){
        sum_theta = sum_theta + theta[i];
    }
    cout << "Prior sum_theta: " << sum_theta << endl;

    // scale prior causal probability
   // if ( ((double)(snp_pos.size()) / 1000000.0) < sum_theta) {
     //  cout << "Scale prior CPP to 1e-6 per SNP ... \n";
    //  for (size_t i=0; i < snp_pos.size(); ++i){
    //    theta[i] = (theta[i] / sum_theta) * ((double)(snp_pos.size()) / 1000000.0);
    //  }
   // }

    log_theta.clear();
    log_qtheta.clear();
    sum_logqtheta = 0.0;
    for (size_t i=0; i < snp_pos.size(); ++i){
        log_theta.push_back(log(theta[i]));
        log_qtheta.push_back(log(1.0 - theta[i]));
        sum_logqtheta += log_qtheta[i] ;
    }
    cout << "sum of log( 1 - CPP ) of all SNPs : " << sum_logqtheta << endl;
}


//InitialMCMC currently used
//Lei's change
/*
void BVSRM::InitialMCMC (gsl_matrix *X, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos)

{
    //double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
    double q_genome=gsl_cdf_chisq_Qinv(5e-8, 1);
    //cout << "significant chisquare value : " << q_genome << endl;
    cHyp.n_gamma=0;
    for (size_t i=0; i<pos_loglr.size(); ++i) {
        if (pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
    }
    //cout << "number of snps before adjust = " << cHyp.n_gamma << endl;
    if (cHyp.n_gamma<30) {cHyp.n_gamma=30;}
    if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
    if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}


    if (!iniSNPfile.empty() && iniType == 0) {

        ifstream infile(iniSNPfile.c_str(), ifstream::in);
        if(!infile) {cout << "Error opening file " << iniSNPfile << endl; exit(-1);}
        string lineid;
        rank.clear();
        size_t orderj, rankj;

        cout << "Start loading initial snp IDs from " << iniSNPfile << "\n";
        while (!safeGetline(infile, lineid).eof()) {
            orderj = 0;
            for (size_t i=0; i < snp_pos.size(); i++) {
                if (snp_pos[i].rs.compare(lineid) == 0) {
                    orderj=i;
                    rankj = SNPorder_vec[orderj].second;
                    rank.push_back(rankj);
                    //cout << lineid << " with rank = " << rankj;
                    //snp_pos[orderj].printMarker();
                    break;
                }
            }
        }
        infile.close();
        infile.clear();

        if (rank.size() == 0) {
            for (size_t i=0; i<cHyp.n_gamma; ++i) {
                rank.push_back(i);
            }
        } //take rank 0 if tracked no SNPs from the iniSNPfile

        cHyp.n_gamma = rank.size();
    }
    else if(iniType == 0) {iniType = 1;}
    else if(iniType == 1) {
        cout << "Start with top variants.\n";
        rank.clear();
        for (size_t i=0; i<cHyp.n_gamma; ++i) {
            rank.push_back(i);
        }
    } // Start with most significant variants from SVT
    else if(iniType == 2) {
        cout << "Start with top SVT SNPs/ColinearTest.\n";
        size_t posr, j=0, i=0;
        double xtx;

        rank.clear();
        rank.push_back(0);
        posr = SNPrank_vec[0].first;
        xtx = xtx_vec[posr];
        // cout << "rank added: " << 0 << ", ";

        gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
        gsl_vector * xvec = gsl_vector_alloc(ni_test);
        //lei's change
        //copy the column
        gsl_matrix_get_col(xvec,X,posr);
        //getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); //get geno column
        
        gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
        gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);

        do{
       // for (size_t i=1; i < cHyp.n_gamma; ++i){
            gsl_matrix_set_col(Xr, i, xvec);
            gsl_matrix_set(XtXr, i, i, xtx);

            if (i>0) {
                gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
                gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
                gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);

                gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, i, 0, i);
                gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, i, 0, i);
                gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
                gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);

            }

            if ((rank.size() < cHyp.n_gamma)) {

                do{
                    j++; // Consider rank j
                    //cout << "consider rank j" << j << endl;
                } while( (j < snp_pos.size()) && ColinearTest(X, Xr, XtXr, j, rank.size()));

                rank.push_back(j);
                posr = SNPrank_vec[j].first;
                xtx = xtx_vec[posr];
                //lei's change
                //copy the column
                gsl_matrix_get_col(xvec,X,posr);
                //getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); // get geno column

            }
            i++;

        }while(i < (cHyp.n_gamma));

        //PrintMatrix(XtXr, cHyp.n_gamma, cHyp.n_gamma);

        gsl_matrix_free(Xr);
        gsl_matrix_free(XtXr);
        gsl_vector_free(xvec);
        gsl_vector_free(Xtxvec);

    } // Start with most significant variants from SVT
    else if(iniType == 3){
        cout << "\nStart with Step-wise selected variants. \n";
    vector<pair<size_t, double> > rank_loglr;
    size_t posr, radd;

    size_t topMarkers=500;
    if(ns_test<500){topMarkers = ns_test;}
        
    for (size_t i=1; i<topMarkers; ++i) {
        rank_loglr.push_back(make_pair(i, pos_loglr[i].second));
    }
    cout << endl;

    rank.clear();
    rank.push_back(0);
    posr = SNPrank_vec[0].first;
    //cout << "rank added: " << 0 << " with LRT "<< pos_loglr[0].second << "," ;

    gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
    gsl_vector * xvec = gsl_vector_alloc(ni_test);
    //lei's change
    //copy the column
    gsl_matrix_get_col(xvec,X,posr);

    //getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); //get geno column
    
    gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
    gsl_vector * Xtyr = gsl_vector_alloc(cHyp.n_gamma);

    gsl_vector * yres = gsl_vector_alloc(ni_test);
    gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);

    double xty, yty;
    gsl_blas_ddot(Uty, Uty, &yty);

    for (size_t i=1; i < s_max; ++i){
        gsl_matrix_set_col(Xr, (i-1), xvec);
        gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);

        gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
        gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);

        gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, (i-1), 0, i);
        gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, (i-1), 0, i);
        gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
        gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);

        gsl_blas_ddot(xvec, Uty, &xty);
        gsl_vector_set(Xtyr, (i-1), xty);

        // calculate conditional yres
        CalcRes(Xr, Uty, XtXr, Xtyr, yres, i, yty);
        for (size_t j=0; j<rank_loglr.size(); ++j) {
            posr = SNPrank_vec[rank_loglr[j].first].first;
            //lei's change
            //copy the column
            gsl_matrix_get_col(xvec,X,posr);
            //getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable); // get geno column

            rank_loglr[j].second = CalcLR(yres, xvec, posr);
        }
        stable_sort (rank_loglr.begin(), rank_loglr.end(), comp_lr); //sort the initial rank.

        if (rank_loglr[0].second > q_genome) {
            radd = rank_loglr[0].first;
           // cout << "XtXr : "; PrintMatrix(XtXr, rank.size(), rank.size());
            if (ColinearTest(X, Xr, XtXr, radd, rank.size())) {
                continue;
            }
            else {
                posr = SNPrank_vec[radd].first;
                //lei's change
                //copy the column
                gsl_matrix_get_col(xvec,X,posr);
                //getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);

                rank.push_back(radd);
                rank_loglr.erase(rank_loglr.begin());
                //cout << "rank added: " << radd << " with LRT "<< rank_loglr[0].second << "," ;
            }
        }
        else break;
    }
        cHyp.n_gamma = rank.size();
        //cout <<"Initial XtX: \n"; PrintMatrix(XtXr, cHyp.n_gamma, cHyp.n_gamma);
        gsl_matrix_free(Xr);
        gsl_matrix_free(XtXr);
        gsl_vector_free(Xtyr);
        gsl_vector_free(xvec);
        gsl_vector_free(yres);
        gsl_vector_free(Xtxvec);
    }
    //cout << "number of snps = " << cHyp.n_gamma << endl;
    //stable_sort (rank.begin(), rank.end(), comp_vec); //sort the initial rank.
    cout << "Starting model has variants with ranks: \n"; PrintVector(rank);
    
    cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
    cHyp.h=pve_null;

    if (cHyp.logp==0) {cHyp.logp=-0.000001;}
    if (cHyp.h<1e-8) {cHyp.h=0.1;}

    gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
    SetXgamma (UtXgamma, X, rank);

    cout << "trace_G = " << trace_G << endl;
    //double sigma_a2;
    //if (trace_G!=0) {
    //    sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp));
    //} else {
    //    sigma_a2=cHyp.h*1.0/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
    //}
    //if (sigma_a2==0) {sigma_a2=0.025;}
    // cout << "initial sigma_a2 = " << sigma_a2 << endl;

    //cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
    gsl_matrix_free (UtXgamma);

    //if (cHyp.rho>1.0) {cHyp.rho=1.0;}
    if (cHyp.h<h_min) {cHyp.h=h_min;}
    if (cHyp.h>h_max) {cHyp.h=h_max;}
    //if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
    //if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
    if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
    if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}

    //cout << "start setHyp... \n";
    double tau_beta_temp = 0.1;
    setHyp(tau_beta_temp);

    cHyp.theta = theta;
    cHyp.log_theta = log_theta;

    // cout<<"initial value of theta_vec = "; PrintVector(theta);
    cout<<"Initially selected number of variants in the model = "<<cHyp.n_gamma<<endl;
    return;
}
*/


//used in EM_BLock
double BVSRM::CalcLikegamma(const vector<size_t> &rank)
{
    double loglikegamma = sum_logqtheta;
    size_t pos_i;
    if(rank.size() > 0)
    {
        for (size_t i=0; i < rank.size(); i++) {
            pos_i = mapRank2pos[rank[i]] ;
            loglikegamma =  loglikegamma + (log_theta[pos_i] - log_qtheta[pos_i]);
        }
    }
    return loglikegamma;
}

/*
// for the new model
double BVSRM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag, double &loglike)
{
    //conditioning on hyper parameters: subvar, log_theta
    double logpost=0.0;
    double d;
    double logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;

    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);

    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);

    //calculate Omega
    gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
   // cout << "Omega : "; PrintMatrix(Omega, 5, 5);
    CalcXVbeta(Omega, &sigma_sub.vector);
    gsl_vector_view Omega_diag = gsl_matrix_diagonal(Omega);
    gsl_vector_add_constant(&Omega_diag.vector, 1.0);

    if(LapackSolve(Omega, &Xty_sub.vector, beta_hat)!=0)
       EigenSolve(Omega, &Xty_sub.vector, beta_hat);
    logdet_O=LapackLogDet(Omega);

    //cout << "beta_hat from solve : "; PrintVector(beta_hat);
    gsl_vector_mul(beta_hat, &sigma_sub.vector);
    //cout << "beta_hat: "; PrintVector(beta_hat);
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);

    double bxy;
    gsl_blas_ddot (&Xty_sub.vector, beta_hat, &bxy);
    double R2 = bxy / yty;

    if (R2 > 1.0 || R2 < -0.0) {

        //cout << "R2 in CalcPosterior = " << R2 << endl;
        Error_Flag=1;
    }
    else{
        Error_Flag=0;
        gsl_vector_memcpy(&beta_sub.vector, beta_hat);
        gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);
        gsl_blas_ddot (Xb, Xb, &d);
        if (a_mode==11) {
            cHyp.pve=d/(double)ni_test;
            //cHyp.pve/=cHyp.pve+1.0/tau;
            // cHyp.pge=1.0;
        }
    }

    logpost = tau * bxy;
    loglike = -0.5 * ((double)cHyp.n_gamma * logrv + (double)cHyp.m_gamma[0] * log_subvar[0] + (double)cHyp.m_gamma[1] * log_subvar[1] - logpost);

    logpost = -0.5 * (logdet_O - logpost);

    gsl_matrix_free (Omega);
    gsl_vector_free (beta_hat);

    return logpost;
}
*/

/*
// for the new model
//calculate likelihood P(Y | gamma, subvar, theta)
double BVSRM::CalcLikelihood (const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag)
{
    //double sigma_a2=cHyp.h/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test * trace_G);

    double loglike=0.0;
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;

  if (s_size == 0) {
        //calculate likelihood if ngamma=0
        if (a_mode==11) {loglike-=0.5*(double)ni_test*log(yty);}
        else {loglike-=0.5*yty;}
  }

  else{
   // gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);

      gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
      gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
      gsl_vector *beta_hat=gsl_vector_alloc (s_size);
      gsl_vector *Xty_temp=gsl_vector_alloc (s_size);

      //calculate Omega
      gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
      CalcXVbeta(Omega, &sigma_sub.vector);
      gsl_matrix_set_identity (M_temp);
      gsl_matrix_add (Omega, M_temp);

      //calculate beta_hat
      gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
      logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);    //solve Omega * beta_hat = Xty for beta_hat
      // Omega was inverted here
      // logdet_0 = det(Omega)
      gsl_vector_mul(beta_hat, &sigma_sub.vector);
      gsl_blas_ddot (Xty_temp, beta_hat, &d);
      P_yy-=d;
      if (P_yy <= 0) {
          Error_Flag = 1;
          cout << "Error in calclikelihood: P_yy = " << P_yy << endl;
          cout << "h = "<<setprecision(6) << cHyp.h << "; rho: " << cHyp.rho_vec[0] << ", " << cHyp.rho_vec[1];
          cout << "theta: "<<setprecision(6) << exp(cHyp.log_theta[0]) << ", " << exp(cHyp.log_theta[1]);
          cout << "; subvar: "<<setprecision(6) << cHyp.subvar[0] << ", " << cHyp.subvar[1];
          cout << "beta_hat: "; PrintVector(beta_hat);

          cout << "set beta_hat to 0\n";
          gsl_vector_set_zero(beta_hat);
          P_yy = yty;
         // exit(-1);
      }
      else {Error_Flag = 0;}

    loglike = -0.5 * logdet_O;
    if (a_mode==11) {loglike -= 0.5 * (double)ni_test * log(P_yy);}
    else {loglike -= 0.5*P_yy;}

    gsl_matrix_free (Omega);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (Xty_temp);
  }

    return loglike;
}
*/


//currently used
double BVSRM::ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, gsl_matrix *X, const gsl_vector *z, const gsl_matrix *Xgamma_old, const gsl_matrix *XtX_old, const gsl_vector *Xtz_old, const double &ztz, gsl_matrix *Xgamma_new, gsl_matrix *XtX_new, gsl_vector *Xtz_new)
{
    map<size_t, int> mapRank2in;
    double unif, logp = 0.0;
    size_t r_add, r_remove, col_id, r;

    rank_new.clear();
    if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"size wrong"<<endl;}
    if (cHyp_old.n_gamma!=0) {
        for (size_t i=0; i<rank_old.size(); ++i) {
            r=rank_old[i];
            rank_new.push_back(r);
            mapRank2in[r]=1;
        }
    }
    cHyp_new.n_gamma=cHyp_old.n_gamma;

    //for (size_t i=0; i<repeat; ++i) {

        unif=gsl_rng_uniform(gsl_r);

        if (unif < 0.33 && cHyp_new.n_gamma<s_max) {flag_gamma=1;}
        else if (unif>=0.33 && unif < 0.67 && cHyp_new.n_gamma>s_min) {flag_gamma=2;}
        else if (unif>=0.67 && cHyp_new.n_gamma>0 && cHyp_new.n_gamma<ns_test) {flag_gamma=3;}
        else {flag_gamma=0;}

        if(flag_gamma==1)  {//add a snp;
            // cout << "add a snp" << endl;

          /*  if (rank_old.size() > 0) {
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0) \
                    || (ColinearTest(X, Xgamma_old, XtX_old, r_add, rank_old.size()))
                    );
            }
            else {*/
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0));
           // }

            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }

            mapRank2in[r_add]=1;
            rank_new.push_back(r_add);
            cHyp_new.n_gamma++;
            logp += -log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);

            if (rank_old.size()>0) {
                SetXgammaAdd(X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, r_add, Xgamma_new, XtX_new, Xtz_new);
            }
            else{
                SetXgamma (Xgamma_new, X, rank_new);
                CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
            }
           // cout << "XtX from set calcXtX: \n"; PrintMatrix(XtX_new, rank_new.size(), rank_new.size());
           // cout << "succesfully added a snp" << endl;

        }
        else if (flag_gamma==2) {//delete a snp;
            //  cout << "delete a snp" << endl;

            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];

            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            prob_total+=p_gamma[r_remove];

            mapRank2in.erase(r_remove);
            rank_new.erase(rank_new.begin()+col_id);
            logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
            cHyp_new.n_gamma--;

            if (rank_new.size() > 0) {
                SetXgammaDel(Xgamma_old, XtX_old, Xtz_old, rank_old, col_id, Xgamma_new, XtX_new, Xtz_new);
            }
            // cout << "succesfully deleted a snp" << endl;
        }
        else if (flag_gamma==3) {//switch a snp;
            // cout << "switch a snp" << endl;
            long int o_add, o_remove;
            long int o_rj, o_aj;
            size_t j_add, j_remove, o;

            gsl_ran_discrete_t *gsl_s, *gsl_a; //JY added dynamic gsl_s
            double *p_BFr = new double[ns_neib];
            double *p_BFa = new double[ns_neib];

            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];//careful with the proposal
            if(mapRank2in.count(r_remove) == 0) {cout << "wrong proposal of r_remove;" << endl; exit(-1);}
            o_remove = SNPrank_vec[r_remove].second;
            rank_new.erase(rank_new.begin()+col_id);
            size_t s_size = rank_new.size();
            mapRank2in.erase(r_remove);

            //cout << "s_size = "<< s_size << "; s_size+1 = " << s_size+1 << endl;
            gsl_matrix *Xgamma_temp=gsl_matrix_alloc (ni_test, s_size+1);
            gsl_matrix *XtX_gamma=gsl_matrix_alloc (s_size+1, s_size+1);
            gsl_vector *Xtz_gamma=gsl_vector_alloc (s_size+1);
            gsl_vector *z_res = gsl_vector_alloc(ni_test);

            //cout << "Switch step rank_old:"; PrintVector(rank_old);
            //cout <<"XtX_old: "; PrintMatrix(XtX_old, rank_old.size(), rank_old.size());
            //cout << "temp rank_new:"; PrintVector(rank_new);
            if (s_size > 0) {
                //SetXgamma (Xgamma_temp, X, rank_new);
                //CalcXtX (Xgamma_temp, z, s_size, XtX_gamma, Xtz_gamma);
                //cout << "XtX from set calcXtX: \n"; PrintMatrix(XtX_gamma, s_size, s_size);

                SetXgammaDel(Xgamma_old, XtX_old, Xtz_old, rank_old, col_id, Xgamma_temp, XtX_gamma, Xtz_gamma);
                //cout << "XtX from set XgammaDel success: \n";
                //PrintMatrix(XtX_gamma, s_size, s_size);

                CalcRes(Xgamma_temp, z, XtX_gamma, Xtz_gamma, z_res, s_size, ztz);
                gsl_s = MakeProposal(o_remove, p_BFr, X, z_res, mapRank2in);
            }
            else {
                gsl_s = MakeProposal(o_remove, p_BFr, X, z, mapRank2in);
            }

            j_add = gsl_ran_discrete(gsl_r, gsl_s);
            o_add = (o_remove - win) + j_add;
            if((o_add < 0) || (o_add >= (long int)ns_test) || (o_add == (long int)o_remove))
              //  cout << "ERROR proposing switch snp"; //new snp != removed snp
            r_add = SNPorder_vec[(size_t)o_add].second;

            //cout << "XtX from set Xgamma: \n"; PrintMatrix(XtX_gamma, s_size, s_size);
            if (s_size>0) {

                /*if (ColinearTest(X, Xgamma_temp, XtX_gamma, r_add, s_size)) {
                    flag_gamma=-1;
                    //cout << "Failed colinear test in switch" << endl;
                }
                else{*/
                    gsl_a = MakeProposal(o_add, p_BFa, X, z_res, mapRank2in);

                    double prob_total_remove=1.0;
                    double prob_total_add=1.0;

                    for (size_t ii=0; ii<rank_new.size(); ++ii) {
                        r = rank_new[ii];
                        o = SNPrank_vec[r].second;
                        o_rj = ((long int)o - o_remove) + win;
                        o_aj = ((long int)o - o_add) + win;
                        if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
                        if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
                    }

                    j_remove = o_remove - o_add + win;
                    logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                    logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)

                    SetXgammaAdd(X, Xgamma_temp, XtX_gamma, Xtz_gamma, z, rank_new, r_add, Xgamma_new, XtX_new, Xtz_new);
                   // cout << "XtX from setXgammaAdd success: \n";
                    //PrintMatrix(XtX_new, s_size+1, s_size+1);

                    mapRank2in[r_add]=1;
                    rank_new.push_back(r_add);

                    gsl_ran_discrete_free(gsl_a);
                //}
            }
            else{
                //construct gsl_s, JY
                //cout << "o_add = " << o_add <<  "; r_add = "<<r_add << endl;
                gsl_a = MakeProposal(o_add, p_BFa, X, z, mapRank2in);

                double prob_total_remove=1.0;
                double prob_total_add=1.0;

                for (size_t ii=0; ii<rank_new.size(); ++ii) {
                    r = rank_new[ii];
                    o = SNPrank_vec[r].second;
                    o_rj = ((long int)o - o_remove) + win;
                    o_aj = ((long int)o - o_add) + win;
                    if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
                    if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
                }

                j_remove = o_remove - o_add + win;
                logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)

                mapRank2in[r_add]=1;
                rank_new.push_back(r_add);
                SetXgamma (Xgamma_new, X, rank_new);
                CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);

                gsl_ran_discrete_free(gsl_a);
            }

            gsl_matrix_free(Xgamma_temp);
            gsl_matrix_free(XtX_gamma);
            gsl_vector_free(Xtz_gamma);
            gsl_vector_free(z_res);

            gsl_ran_discrete_free(gsl_s);

            delete[] p_BFr;
            delete[] p_BFa;
            // cout << "successfully switched a snp" << endl;
        }

        else {logp+=0.0;}//do not change
    //}
    mapRank2in.clear();
    return logp;
}

// updated (06/15/2022) by JY
void BVSRM::WriteHyptemp(gsl_vector *LnPost){
    
    double em_logpost = 0.0, logpost_max =  gsl_vector_max(LnPost);
    // cout << "logpost_max = " << logpost_max << endl;
    for (size_t i=0; i < s_step; i++) {
        em_logpost += exp(gsl_vector_get(LnPost, i) - logpost_max);
    }
    em_logpost /= double(s_step);
    em_logpost = log(em_logpost) + logpost_max;
    
    //save E(file_out, lnpost, GV, rv, n[i], Gvec[i], m[i], sigma2[i])
    string file_hyp;
    file_hyp = "./output/" + file_out;
    file_hyp += ".hyptemp";

    ofstream outfile_hyp;
    // GV: phenotype variation explained by genotypes, r2

    // write *.hyptemp
    outfile_hyp.open (file_hyp.c_str(), ofstream::out);
    if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
    cout << "Start writing hyptemp ... \n" ;
    cout << "Sum CPP = " << em_gamma << endl;

    // Add a file header starting with # (06/15/2022)
    outfile_hyp << "#genome_block_prefix\tlog_post_likelihood\tr2\tsum_gamma\tsum_beta2\n";
    outfile_hyp << file_out << "\t"; //genome_block_prefix
    // log_pos_likelihood ; r2 ;
    outfile_hyp << scientific << setprecision(3) << em_logpost << "\t" << (GV / (double)s_step) << "\t" << em_gamma << "\t"  << sumbeta2 << endl;

    outfile_hyp.clear();
    outfile_hyp.close();
}



/*
//MCMC function (old version; not used)
// lei's change
void BVSRM::MCMC (gsl_matrix *X, const gsl_vector *y, bool original_method) {

    if (original_method) {
        cout << "Run MCMC...\n";
    }
    // cout << "# of unique function types = " << n_type << endl;

    //new model related
    gsl_vector *sigma_subvec_old = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_old);
    gsl_vector *sigma_subvec_new = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_new);
    gsl_vector *LnPost = gsl_vector_alloc(s_step); //save logPost...

    em_gamma.assign(n_type, 0.0);
        //save sum of m_q, beta2_q;
    GV = 0.0; 
    sumbeta2.assign(n_type, 0.0);
    
    gsl_vector *Xb_new=gsl_vector_alloc (ni_test);
    gsl_vector *Xb_old=gsl_vector_alloc (ni_test);
    gsl_vector *z_hat=gsl_vector_alloc (ni_test);
    gsl_vector *z=gsl_vector_alloc (ni_test);

    gsl_matrix *Xgamma_old=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_old=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_old=gsl_vector_alloc (s_max);
    gsl_vector *beta_old=gsl_vector_alloc (s_max);

    gsl_matrix *Xgamma_new=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_new=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_new=gsl_vector_alloc (s_max);
    gsl_vector *beta_new=gsl_vector_alloc (s_max);

    double ztz=0.0, mean_z=0.0;
    gsl_vector_memcpy (z, y);

    gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    pheno_var = ztz / ((double)(ni_test-1)) ;

    //cout << "ztz = " << ztz << "; phenotype variance = " << ztz / pheno_var <<endl;
    gsl_vector_scale(z, 1.0 / sqrt(pheno_var)); // standardize phenotype z
    gsl_blas_ddot(z, z, &ztz); // calculate ztz for one more time after standardization

    //Initialize variables for MH
    double logPost_new, logPost_old, loglike_new, loglike_old, loglikegamma;
    double logMHratio;
    vector<size_t> rank_new, rank_old;
    class HYPBSLMM cHyp_old, cHyp_new;
    bool Error_Flag=0;

    if (a_mode==13) {
        pheno_mean=0.0;
    }
    vector<pair<double, double> > beta_g; //save beta estimates
    for (size_t i=0; i<ns_test; i++) {
        beta_g.push_back(make_pair(0.0, 0.0));
    }

    // UcharTable and vector<SNPPOS> snp_pos were created in CALCSS SS

    vector<pair<size_t, double> > pos_loglr;
    vector<double> pval_lrt;
    vector<double> Z_scores;

    cout << "Calculating Z_scores, standard errors of effect-sizes, LRT statistics, pvals ... \n";
    MatrixCalcLmLR (indicator_snp,X, z, pos_loglr, ns_test, ni_test, SNPmean, Gvec, xtx_vec, Z_scores, mbeta, mbeta_SE, pval_lrt, snp_pos); //calculate trace_G or Gvec, Z_scores, beta_SE
 //calculate trace_G or Gvec
    trace_G = VectorSum(Gvec) / double(ns_test);
    cout << "Trace of Genotype Matrix = " << trace_G << endl;

    for(size_t i=0; i < n_type; i++){
        Gvec[i] /= mFunc[i];
    } 
    //cout << "Avg trace of genotypes per category : " ; PrintVector(Gvec); 
    
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // order snp_pos by chr/bp
    stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio

    //Jingjing's edit, create maps between rank and order
    size_t pos;
    for (size_t i=0; i<ns_test; ++i) {
        mapRank2pos[i]=pos_loglr[i].first;
        mapPos2Rank[pos_loglr[i].first] = i;

        mapOrder2pos[i] = snp_pos[i].pos;
        mapPos2Order[snp_pos[i].pos] = i;
    }
    
    for (size_t i=0; i<ns_test; ++i) {
        pos = mapRank2pos[i];
        mapRank2Order[i] = mapPos2Order[pos];

        pos = mapOrder2pos[i];
        mapOrder2Rank[i] = mapPos2Rank[pos];
    }

    SNPorder_vec.clear();
    SNPrank_vec.clear();
    for (size_t i=0; i<ns_test; i++) {
        SNPorder_vec.push_back(make_pair(snp_pos[i].pos, mapOrder2Rank[i]));
        SNPrank_vec.push_back(make_pair(pos_loglr[i].first, mapRank2Order[i]));
    }


    //end of Jingjing's edit

    //Calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
    gsl_rng_env_setup();
    const gsl_rng_type * gslType;
    gslType = gsl_rng_default;
    if (randseed<0)
    {
        time_t rawtime;
        time (&rawtime);
        tm * ptm = gmtime (&rawtime);

        randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
    }
    gsl_r = gsl_rng_alloc(gslType);
    gsl_rng_set(gsl_r, randseed);
    p_gamma = new double[ns_test]; //defined in bvsrm.h

    size_t p_gamma_top=0;
    for(size_t i=0; i < pval_lrt.size(); i++){
        if (pval_lrt[i] < 5e-8) p_gamma_top++;
    }
    cout << "Number of variants with p-value < 5e-8 : " << p_gamma_top << endl;

    // calculate discrete distribution for gamma
    SetPgamma (p_gamma_top);

    gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma); // set up proposal function for gamma

    pheno_var *= ztz / ((double)(ni_test-1));
    rv = pheno_var;
    //cout<< "Fix Residual Variance = " << rv << endl;
    //cout << "tau = " << tau << "; log(2pi*rv) = " <<logrv << endl;

    //Initial parameters

    cout << "\nStart initializing MCMC ... \n";
    InitialMCMC (X, z, rank_old, cHyp_old, pos_loglr, snp_pos); // Initialize rank and cHyp
       
    inv_subvar.assign(n_type, 0.0), log_subvar.assign(n_type, 0.0);
    for(size_t i=0; i < n_type; i++){
        inv_subvar[i] = (1.0 / subvar[i]); 
        log_subvar[i] = (log(subvar[i])); 
    }
    cout << "inv_subvar = "; PrintVector(inv_subvar);
    //cout << "log_subvar = "; PrintVector(log_subvar);

    if (cHyp_old.n_gamma > 0) {
        SetXgamma (Xgamma_old, X, rank_old);
        CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
    }

    //cout << "Set m_gamma... \n";
    set_mgamma(cHyp_old, rank_old, snp_pos);
    cout << "Initial number of selected variants per category : ";
    PrintVector(cHyp_old.m_gamma); 
    //cout << "Set sigma_subvec... \n";
    // getSubVec(sigma_subvec_old, rank_old, snp_pos);
    //PrintVector(sigma_subvec_old, rank_old.size());

    cHyp_initial=cHyp_old;
    gsl_vector_memcpy(sigma_subvec_new, sigma_subvec_old);

    //Calculate first loglikelihood
    //cout << "first calculating logpost ... \n";
    if (cHyp_old.n_gamma==0) {
        loglikegamma = CalcLikegamma(rank_old);
        logPost_old = CalcPosterior (ztz, cHyp_old) + loglikegamma;
        loglike_old = loglikegamma;
    }
    else {
        loglikegamma = CalcLikegamma(rank_old);
        logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag, loglike_old) + loglikegamma;
        loglike_old += loglikegamma;
    }
    if (!Error_Flag) {
       // cout <<  "logPost_old = " << logPost_old << endl;
    }
    else {
        cerr << "Failed at initialMCMC...\n";
        exit(-1);
    }
    //cout <<  "Initial logPost_old = " << logPost_old << endl;

    //calculate centered z_hat, and pve
    if (a_mode==13) {
        if (cHyp_old.n_gamma==0) {
            CalcCC_PVEnZ (z_hat, cHyp_old);
        }
        else {
            CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
        }
    }

    //Start MCMC
    size_t k_save_sample=0;
    w_pace=3000;
    int accept; // accept_theta; naccept_theta=0,
    size_t total_step=w_step+s_step;
    size_t repeat=1;
    int flag_gamma=0;
    double accept_percent, betai; // accept_theta_percent;


    cHyp_new = cHyp_old;
    rank_new = rank_old;

    vector <string> snps_mcmc; // save locations of included snps per iteration
    string snps_mcmc_temp;
    size_t order_i;

    for (size_t t=0; t<total_step; ++t) {

       if (t%d_pace==0 || t==total_step-1) {ProgressBar ("\nRunning MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));
                cout << endl;
            }
        //      if (t>10) {break;}

        if (a_mode==13) {
            SampleZ (y, z_hat, z); //sample z
            mean_z=CenterVector (z);
            gsl_blas_ddot(z,z,&ztz);

            //First proposal need to be revised
            if (cHyp_old.n_gamma==0) {
                loglikegamma = CalcLikegamma(rank_old);
                logPost_old = CalcPosterior (ztz, cHyp_old) + loglikegamma;
                loglike_old = loglikegamma;
            } else {
                gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
                gsl_vector_view Xtz_sub=gsl_vector_subvector(Xtz_old, 0, rank_old.size());
                gsl_blas_dgemv (CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0, &Xtz_sub.vector);
                loglikegamma = CalcLikegamma(rank_old);
                logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag, loglike_old) + loglikegamma;
                loglike_old += loglikegamma;
            }
        }

        //////// Set repeat number
        //if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
        //else {repeat=1;}
        //cout << "n_mh = " << n_mh << endl;

        for (size_t i=0; i<n_mh; ++i) {

            //cout << "\n \n propose gamam...\n";
            //cout << "old rank: "; PrintVector(rank_old);
            //repeat = 1;
            logMHratio = ProposeGamma (rank_old, rank_new, cHyp_old, cHyp_new, repeat, X, z, Xgamma_old, XtX_old, Xtz_old, ztz, Xgamma_new, XtX_new, Xtz_new); //JY
           // rank_new.clear(); cHyp_new.n_gamma=0;
            //cout << "propose new rank: "; PrintVector(rank_new);
            //cout << "flag_gamma = " << flag_gamma << endl;
            //cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
            //cout << "propose gamma logMHratio = "<<logMHratio << "; MHratio = " << exp(logMHratio) << endl;

        if (flag_gamma > 0) {

            if(flag_gamma==1) nadd++;
            else if(flag_gamma==2) ndel++;
            else  nswitch++;

            if (rank_new.size() > 0) {
                set_mgamma(cHyp_new, rank_new, snp_pos);
                getSubVec(sigma_subvec_new, rank_new, snp_pos);
                loglikegamma = CalcLikegamma(rank_new);
                logPost_new = CalcPosterior (Xgamma_new, XtX_new, Xtz_new, ztz, Xb_new, beta_new, cHyp_new, sigma_subvec_new, Error_Flag, loglike_new) + loglikegamma;
                loglike_new += loglikegamma;
            }
            else {
                cHyp_new.m_gamma.assign(n_type, 0);
                loglikegamma = CalcLikegamma(rank_new);
                logPost_new = CalcPosterior (ztz, cHyp_new) + loglikegamma;
                loglike_new = loglikegamma;
            }
           // cout << "new m_gamma: " << cHyp_new.m_gamma[0] << ", "<< cHyp_new.m_gamma[1]<< endl;

             // cout << "Calcposterior success." << endl;
            if (!Error_Flag) {
                logMHratio += logPost_new-logPost_old;
                //cout <<"logPost_old = " << logPost_old<< "; logPost_new = "<< logPost_new<< "\n logMHratio = " << logMHratio<< "; MHratio = " << exp(logMHratio) << endl;
                if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio)
                    { accept=1; if (flag_gamma < 4) n_accept++;}
                else {accept=0;}
            }
            else{
                accept=0;
            }
        }
        else{
            nother++;
            accept = 0;
        }

            //cout << "accept = " << accept << endl;

            if (accept==1) {
                    if(flag_gamma==1) nadd_accept++;
                    else if(flag_gamma==2) ndel_accept++;
                    else if(flag_gamma==3) nswitch_accept++;
                    else nother_accept++;

                    logPost_old=logPost_new;
                    loglike_old = loglike_new;
                    cHyp_old.n_gamma = cHyp_new.n_gamma;
               // cout << "cHyp_old.m_gamma = "; PrintVector(cHyp_old.m_gamma);
                    cHyp_old.m_gamma = cHyp_new.m_gamma;
                //cout << "cHyp_new.m_gamma = "; PrintVector(cHyp_new.m_gamma);
                //cout << "cHyp_old.m_gamma = "; PrintVector(cHyp_old.m_gamma);
                    cHyp_old.pve = cHyp_new.pve;
                    //cHyp_old.rv = cHyp_new.rv;
                    //cHyp_old.pge = cHyp_new.pge;
                    gsl_vector_memcpy (Xb_old, Xb_new);
                rank_old.clear();
                for (size_t i=0; i<rank_new.size(); i++) {
                    rank_old.push_back(rank_new[i]);
                }
                if (rank_old.size() != rank_new.size()) {
                    cerr << "Error: rank_old size != rank_new size\n";
                    exit(-1);
                }
                //cout << "Accept proposal: "; PrintVector(rank_old);

                if(rank_old.size()>0){
                    gsl_vector_view sigma_oldsub=gsl_vector_subvector(sigma_subvec_old, 0, rank_old.size());
                    gsl_vector_view sigma_newsub=gsl_vector_subvector(sigma_subvec_new, 0, rank_old.size());
                    gsl_vector_memcpy(&sigma_oldsub.vector, &sigma_newsub.vector);

                    gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXold_sub=gsl_matrix_submatrix(XtX_old, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtzold_sub=gsl_vector_subvector(Xtz_old, 0, rank_new.size());
                    gsl_vector_view betaold_sub=gsl_vector_subvector(beta_old, 0, rank_new.size());

                    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(Xgamma_new, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtznew_sub=gsl_vector_subvector(Xtz_new, 0, rank_new.size());
                    gsl_vector_view betanew_sub=gsl_vector_subvector(beta_new, 0, rank_new.size());

                    gsl_matrix_memcpy(&Xold_sub.matrix, &Xnew_sub.matrix);
                    gsl_matrix_memcpy(&XtXold_sub.matrix, &XtXnew_sub.matrix);
                    gsl_vector_memcpy(&Xtzold_sub.vector, &Xtznew_sub.vector);
                    gsl_vector_memcpy(&betaold_sub.vector, &betanew_sub.vector);

                    gsl_vector_memcpy(Xb_old, Xb_new);
                }
                //else{
                  //  gsl_vector_set_zero(Xb_old); // set Xb = 0
                //}
            } else {
                cHyp_new.n_gamma = cHyp_old.n_gamma;
                cHyp_new.m_gamma = cHyp_old.m_gamma;
                rank_new.clear();
            }
          //  cout << "copy data from new propose -> old " << endl;
        } //end of n_mh

        //calculate z_hat, and pve
        if (a_mode==13) {
            if (cHyp_old.n_gamma==0) {
                CalcCC_PVEnZ (z_hat, cHyp_old);
            }
            else {
                CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
            }
            //sample mu and update z hat
            gsl_vector_sub (z, z_hat);
            mean_z+=CenterVector(z);
            mean_z+=gsl_ran_gaussian(gsl_r, sqrt(1.0/(double) ni_test) );
            gsl_vector_add_constant (z_hat, mean_z);
        }

         //if (t % 10 == 0 && t > w_step) {
         if ( (t>0) && (t % w_pace == 0) && (t > w_step) ) {
             accept_percent = (double)n_accept/(double)((t+1) * n_mh);
             //cout << "cHyp_old.n_gamma= " << cHyp_old.n_gamma << endl;
             cout << "acceptance percentage = " << setprecision(6) << accept_percent << endl ;
             cout << "# of selected variants per category: " << endl; PrintVector(cHyp_old.m_gamma);
             //cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); cout << endl;
             cout << "loglike: " << loglike_old << endl;
        }

        //Save data
        if (t<w_step) {continue;}
        else {
            //save loglikelihood
            gsl_vector_set (LnPost, k_save_sample, loglike_old);
            GV += cHyp_old.pve;

            if (cHyp_old.n_gamma > 0){

                region_pip += 1.0 ; //count increase if the model has >0 SNPs
                snps_mcmc_temp="";

                for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
                    // beta_g saved by position
                    pos=SNPrank_vec[rank_old[i]].first;
                    order_i = SNPrank_vec[rank_old[i]].second;

                    betai = gsl_vector_get(beta_old, i);
                    beta_g[pos].first += betai;
                    beta_g[pos].second += 1.0;

                    snps_mcmc_temp += string(snp_pos[order_i].rs) + string(":") + string(snp_pos[order_i].chr) + string(":") + to_string(snp_pos[order_i].bp) + string(":") + to_string(snp_pos[order_i].a_major) + string(":") + to_string(snp_pos[order_i].a_minor)+ string(";");
                }
                snps_mcmc.push_back(snps_mcmc_temp);

                //if(cHyp_old.m_gamma[0]>0)
                   // sample_sigma0.push_back(make_pair(sumbeta2[0] /(double)cHyp_old.m_gamma[0], cHyp_old.m_gamma[0] ));
                //if(cHyp_old.m_gamma[1]>0)
                   // sample_sigma1.push_back(make_pair(sumbeta2[1] /(double)cHyp_old.m_gamma[1], cHyp_old.m_gamma[1] ));

                for(size_t j=0; j < n_type; j++){
                    if(cHyp_old.m_gamma[j] > 0) 
                        em_gamma[j] += (double)cHyp_old.m_gamma[j];
                    }
                k_save_sample++;
                }
            }
    }

    cout << "MCMC completed ... " << endl << endl;
    region_pip = region_pip / double(s_step);
    cout << "region_pip = " << setprecision(5) << region_pip << endl;

    accept_percent = (double)n_accept/(double)(total_step * n_mh);
    cout << "gamma acceptance percentage = " << accept_percent << endl ;
    cout << "# of selected variants per category: "; PrintVector(cHyp_old.m_gamma);
    cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); 
    cout << "loglike: " << loglike_old << endl;
    cout << "k_save_sample = " << k_save_sample << endl;


    //save all marker information
    // WriteGenotypeFile(X, snp_pos);
    //WriteIniSNP(rank_old, snp_pos);
    //WriteParam(beta_g, snp_pos, pos_loglr, Z_scores, pval_lrt);
    //WriteFGWAS_InputFile(snp_pos, Z_scores);


    //Save temp EM results
    WriteHyptemp(LnPost);
    WriteParam(beta_g, snp_pos, pos_loglr, Z_scores, pval_lrt);
    WriteAnnoScore(beta_g, snp_pos);
    //# include a function to calculate and save fitted phenotype values


    WriteMCMC(snps_mcmc); // save all active SNPs from MCMC

   // gsl_matrix_free(Result_hyp);
   // gsl_matrix_free(Result_gamma);
    // gsl_matrix_free(Sample_m);

    gsl_vector_free(sigma_subvec_old);
    gsl_vector_free(sigma_subvec_new);
    gsl_vector_free(LnPost);

    gsl_vector_free(z_hat);
    gsl_vector_free(z);
    gsl_vector_free(Xb_new);
    gsl_vector_free(Xb_old);

    //gsl_vector_free(Xb_mcmc);

    gsl_matrix_free(Xgamma_old);
    gsl_matrix_free(XtX_old);
    gsl_vector_free(Xtz_old);
    gsl_vector_free(beta_old);

    gsl_matrix_free(Xgamma_new);
    gsl_matrix_free(XtX_new);
    gsl_vector_free(Xtz_new);
    gsl_vector_free(beta_new);

    delete [] p_gamma;
    beta_g.clear();

    return;
}
*/
// end of current version


void BVSRM::SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z)
{
    double d1, d2, z_rand=0.0;
    for (size_t i=0; i<z->size; ++i) {
        d1=gsl_vector_get (y, i);
        d2=gsl_vector_get (z_hat, i);
        //y is centerred for case control studies
        if (d1<=0.0) {
            //control, right truncated
            do {
                z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
            } while (z_rand>0.0);
        }
        else {
            do {
                z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
            } while (z_rand<0.0);
        }

        gsl_vector_set (z, i, z_rand);
    }

    return;
}


//JY edit start
void BVSRM::CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX, const gsl_vector *Xtz, gsl_vector *z_res, const size_t &s_size, const double &ztz){

    gsl_matrix_const_view X_gsub=gsl_matrix_const_submatrix(Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_gsub = gsl_matrix_const_submatrix(XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xtz_gsub = gsl_vector_const_subvector(Xtz, 0, s_size);

    gsl_vector *beta_gamma_hat = gsl_vector_alloc(s_size);
    gsl_matrix *XtXtemp = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(XtXtemp, &XtX_gsub.matrix);

    double SSR, R2 ;

    if(LapackSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat) != 0)
        EigenSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat);
    gsl_blas_ddot(&Xtz_gsub.vector, beta_gamma_hat, &SSR);
    R2 = (SSR / ztz);

    /* double lambda = 0.0;
    int k=0;

   if ((R2 < -0.0) || (R2 > 1.0)) {
        cout << "R2 in Calcres = " << R2 << endl;
        PrintMatrix(XtXtemp, s_size, s_size);
        WriteMatrix(&X_gsub.matrix, "_Xres");
        WriteMatrix(XtXtemp, "_XtXres");
        WriteVector(beta_gamma_hat, "_bres");
        WriteVector(z_res, "_zres");
        //exit(-1);
    }*/

if ((R2 < -0.1) || (R2 > 1.1)) {
    /*
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(XtX, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.01;
    // cout << "labmda = " << lambda << endl;
    gsl_vector_view XtXtemp_diag = gsl_matrix_diagonal(XtXtemp);

    while ((R2 < -0.1) || (R2 > 1.1)){

        cout << "add lambda: negative R2 in calcresidual ...  "<< R2 << "; k = " << k << endl;

        gsl_vector_add_constant(&XtXtemp_diag.vector, lambda);
        if(LapackSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat) != 0)
            EigenSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat);
        gsl_blas_ddot(&Xtz_gsub.vector, beta_gamma_hat, &SSR);
        R2 = (SSR / ztz);
        k++;
        cout << "now R2 = " << R2 << "; k = " << k << endl;
        if (k > 9) {
            break;
        }
    }

    if ((R2 < 0.0) || (R2 > 1.0)) {
        gsl_vector_memcpy(z_res, z);
    }*/
    gsl_vector_memcpy(z_res, z);

}
else if ( (R2 < 0.0) || (R2 > 1.0) ){
    gsl_vector_memcpy(z_res, z);
}
else{
    gsl_blas_dgemv(CblasNoTrans, 1.0, &X_gsub.matrix, beta_gamma_hat, 0.0, z_res);
    gsl_vector_scale(z_res, -1.0);
    gsl_vector_add(z_res, z);
}

    gsl_matrix_free(XtXtemp);
    gsl_vector_free(beta_gamma_hat);
    return;
}


//calculate likelihood ratio statistic
double BVSRM::CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec, size_t posj){
    double LR;
    double xtz_res, ztz_res, xtx = xtx_vec[posj];

    gsl_blas_ddot(z_res, z_res, &ztz_res);
    gsl_blas_ddot(x_vec, z_res, &xtz_res);
    //cout << "ztz_res = " << ztz_res << "; xtx = " << xtx << "; xtz_res = " << xtz_res << endl;
    //double ixtx = 1.0 / xtx;
    //double bhat = ixtx * xtz_res;
    //double V = (ixtx * ztz_res - bhat * bhat) / (ni_test);
    //double Z2 = bhat * bhat / V;
    //double VpW = V + Wvar;
    LR = (ni_test)*(log(ztz_res)-log(ztz_res-xtz_res*xtz_res/xtx));
    //sqrt(VpW / V) * exp(-0.5 * Z2 * Wvar / VpW);
    //cout << "log LR = " << BF << ", ";
    return (LR);
}
//lei's change
gsl_ran_discrete_t * BVSRM::MakeProposal(const size_t &o, double *p_BF, gsl_matrix *X, const gsl_vector *z_res, const map<size_t, int> &mapRank2in)
{
    gsl_vector *xvec = gsl_vector_alloc(ni_test);

    long int orderj;
    size_t posj, rank_j;
    vector<int> j_ind;
    double pmax, psum=0.0, countj = 0.0;

    for (size_t j=0; j < ns_neib; ++j){
        orderj = (o - win) + j;
        if((orderj >= 0) && (orderj < (long int)ns_test))
            rank_j = SNPorder_vec[(size_t)orderj].second;
        if((orderj >= 0) && (orderj < (long int)ns_test) && (j != win) && (mapRank2in.count(rank_j) == 0)){
            posj = SNPorder_vec[(size_t)orderj].first;
            //lei's change
            //copy the column
            gsl_matrix_get_col(xvec,X,posj);

            //getGTgslVec(X, xvec, posj, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);

            p_BF[j]=CalcLR(z_res, xvec, posj); //calc loglr
            j_ind.push_back(1);
            countj += 1.0;
        }
        else{
            p_BF[j] = -std::numeric_limits<double>::max();
            j_ind.push_back(0);
        }
    }
    pmax = *std::max_element(p_BF, p_BF+(ns_neib));

    for (size_t j=0; j < ns_neib; ++j){
        if(j_ind[j]==1){
            p_BF[j]=exp(p_BF[j]- pmax);
            psum += p_BF[j];
        }
        else{p_BF[j] = 0.0;}
    }
    //
    psum = 1.0/psum;
    for(size_t j=0; j < ns_neib; ++j){
        p_BF[j] *= psum;
    }

    gsl_vector_free(xvec);

    return (gsl_ran_discrete_preproc(ns_neib, p_BF));
}

//lei's change
bool BVSRM::ColinearTest(gsl_matrix *X, const gsl_matrix * Xtemp, const gsl_matrix * XtX_temp, size_t r_add, size_t s_size)
{
    bool colinear = 0;
    double vreg;
    size_t pos = SNPrank_vec[r_add].first;
    double xtx = xtx_vec[pos];


    gsl_vector *beta_temp = gsl_vector_alloc(s_size);
    gsl_vector *Xtx_temp = gsl_vector_alloc(s_size);
    gsl_vector *xvec_temp = gsl_vector_alloc(ni_test);
    //lei's change
    //copy the column
    gsl_matrix_get_col(xvec_temp,X,pos);
    //getGTgslVec(X, xvec_temp, pos, ni_test, ns_test, SNPmean,CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);

    //gsl_blas_ddot(xvec_temp, xvec_temp, &xtx);
    //cout << "xtx_vec[pos] = "<< xtx_vec[pos] << "; xtx = " << xtx << endl;

    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xtemp, 0, 0, Xtemp->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX_temp, 0, 0, s_size, s_size);
    //cout << "XtX in colinear test: \n "; PrintMatrix(&XtX_sub.matrix, s_size, s_size);
    gsl_blas_dgemv(CblasTrans, 1.0, &Xgamma_sub.matrix, xvec_temp, 0.0, Xtx_temp);

    if (LapackSolve(&XtX_sub.matrix, Xtx_temp, beta_temp) !=0)
        EigenSolve(&XtX_sub.matrix, Xtx_temp, beta_temp);

    gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
    //cout << "vreg = " << vreg << endl;

    double tR2 = 0.95;
    double R2 = (vreg / xtx);
//    int k=0;
//    double lambda = 0.0;
//    gsl_matrix *XtXlu = gsl_matrix_alloc(s_size, s_size);

    if ( (R2 >= tR2) && (R2 <= 1.1) ) {
        colinear = 1;
       // cout << "R2 in ColinearTest = " << R2 << endl;
    }
    else if ((R2 < -0.0) || (R2 > 1.0)){
       colinear = 1;
        //cout << "R2 in ColinearTest = " << R2 << "; k = " << k << endl;
        /*
        //PrintMatrix(&XtX_sub.matrix, s_size, s_size);
        WriteMatrix(&Xgamma_sub.matrix, "_Xct");
        WriteMatrix(&XtX_sub.matrix, "_XtXct");
        WriteVector(beta_temp, "_bct");
        WriteVector(xvec_temp, "_xvec");

        gsl_matrix_memcpy(XtXlu, &XtX_sub.matrix);
        for (size_t i=0; i<s_size; ++i) {
            lambda += gsl_matrix_get(XtXlu, i, i);
        }
        lambda /= (double)s_size;
        lambda *= 0.01;
        gsl_vector_view XtXlu_diag = gsl_matrix_diagonal(XtXlu);

        while ((R2 < -0.1) || (R2 > 1.1)) {

            gsl_vector_add_constant(&XtXlu_diag.vector, lambda);
            if (LapackSolve(XtXlu, Xtx_temp, beta_temp) !=0)
                    EigenSolve(XtXlu, Xtx_temp, beta_temp);
            gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
            R2 = (vreg / xtx);
            k++;
            cout << "now R2 = " << R2 << "; k = " << k << endl;

            if (k>9) {
                cout << "reached k = " << k << endl;
                WriteMatrix(&XtX_sub.matrix, "_XtX_ct");
                WriteVector(Xtx_temp, "_Xtx_ct");
                WriteVector(beta_temp, "_bct");
                break;
            }
        }

        if ((R2 >= tR2) && (R2 <= 1.1)) {
            colinear = 1;
        } else {colinear = 0;} */
    }
    else {
        colinear = 0;
    }

 //   gsl_matrix_free(XtXlu);
    gsl_vector_free(xvec_temp);
    gsl_vector_free(beta_temp);
    gsl_vector_free(Xtx_temp);

    return colinear;
}


bool BVSRM::ColinearTest_SS(const gsl_matrix *D_temp, const gsl_vector * Xtx_temp, gsl_vector * beta_temp)
{
    bool colinear = 0;
    double R2;

    if (LapackSolve(D_temp, Xtx_temp, beta_temp) != 0)
        EigenSolve(D_temp, Xtx_temp, beta_temp);
    gsl_blas_ddot(Xtx_temp, beta_temp, &R2);
    if(R2 < 0.0 || R2 > 1.0) {
        //cout << "ColinearTest_SS: Conditional R2 = " << R2 << endl;
        //perror("Conditional R2 either negative or greater than 1." );
    }
    //
    if ( R2 >= 0.95 ) {
        colinear = 1;
        //cout << "R2 in ColinearTest = " << R2 << endl;
    }

    return colinear;
}


//below fits MCMC for rho=1
void BVSRM::CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty)
{
  time_t time_start=clock();
  gsl_matrix_const_view X_sub=gsl_matrix_const_submatrix(X, 0, 0, X->size1, s_size);
  gsl_matrix_view XtX_sub=gsl_matrix_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_view Xty_sub=gsl_vector_subvector(Xty, 0, s_size);

  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &X_sub.matrix, &X_sub.matrix, 0.0, &XtX_sub.matrix);
  gsl_blas_dgemv(CblasTrans, 1.0, &X_sub.matrix, y, 0.0, &Xty_sub.vector);

  time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

  return;
}


//calculate pve and pge, and calculate z_hat for case-control data
void BVSRM::CalcCC_PVEnZ (gsl_vector *z_hat, class HYPBSLMM &cHyp)
{
  gsl_vector_set_zero(z_hat);
  cHyp.pve=0.0;
  //cHyp.pge=1.0;
  return;
}


//calculate pve and pge, and calculate z_hat for case-control data
void BVSRM::CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp)
{
    double d;

    gsl_blas_ddot (Xb, Xb, &d);
    cHyp.pve=d/(double)ni_test;
    //cHyp.pve/=cHyp.pve+1.0;
    //cHyp.pge=1.0;

    gsl_vector_memcpy (z_hat, Xb);

    return;
}



void BVSRM::SetSSgamma(const vector< vector<double> > &LD, const vector<double> &mbeta, const vector <size_t> &rank, gsl_matrix *D_gamma, gsl_vector *mbeta_gamma)
{
    size_t pos_i, pos_j;
    size_t r_size = rank.size();
    double xtx_ij;
    //cout << "r_size = " << r_size << endl;

    if(r_size > 1) {
        for( size_t i=0; i < r_size; ++i){
            pos_i = mapRank2pos[rank[i]];
            //cout << "pos_i: " << pos_i ;
            gsl_vector_set(mbeta_gamma, i, mbeta[pos_i]);
            gsl_matrix_set(D_gamma, i, i, 1);

            //cout << "; pos_j : " << endl;
            for(size_t j=(i+1); j < (r_size); ++j ){
                pos_j = mapRank2pos[rank[j]];
                //cout << pos_j << "," ;
                xtx_ij = getXtX(LD, pos_i, pos_j);
                gsl_matrix_set(D_gamma, i, j, xtx_ij);
                gsl_matrix_set(D_gamma, j, i, xtx_ij);
            }
            // cout << endl;
        }
    }else{
        pos_i = mapRank2pos[rank[0]];
        //cout << " and position : " << pos_i << endl;
        gsl_vector_set(mbeta_gamma, 0, mbeta[pos_i]);
        gsl_matrix_set(D_gamma, 0, 0, 1.0);
    }
    return;
}

void BVSRM::SetSSgammaAdd (const vector< vector<double> > &LD, const vector<double> &mbeta, const gsl_matrix *D_old, const gsl_vector *mbeta_old, const vector<size_t> &rank_old, size_t ranki, gsl_matrix *D_new, gsl_vector *mbeta_new)
{
    double xtx_i;
    size_t s_size = rank_old.size();
    size_t pos = mapRank2pos[ranki];
    size_t pos_i; // position of the newly proposed SNP rank

    if (s_size==0) {
        cerr << "SetSSgammaAdd rank_old has size 0\n";
        exit(-1);
    }
    //copy rank_old
    gsl_matrix_const_view D11_sub = gsl_matrix_const_submatrix(D_old, 0, 0, s_size, s_size);
    gsl_vector_const_view mbeta1_sub = gsl_vector_const_subvector(mbeta_old, 0, s_size);
    
    gsl_matrix_view D_new11_sub=gsl_matrix_submatrix(D_new, 0, 0, s_size, s_size);
    gsl_vector_view mbeta_new1_sub=gsl_vector_subvector(mbeta_new, 0, s_size);
    
    gsl_matrix_memcpy(&D_new11_sub.matrix, &D11_sub.matrix);
    gsl_vector_memcpy(&mbeta_new1_sub.vector, &mbeta1_sub.vector);
    
    //Set SS for ranki
    gsl_matrix_set(D_new, s_size, s_size, 1.0);
    for(size_t i=0; i < s_size; i++){
        pos_i = mapRank2pos[rank_old[i]];
        xtx_i = getXtX(LD, pos_i, pos); // from calcSS.cpp
        gsl_matrix_set(D_new, s_size, i, xtx_i);
        gsl_matrix_set(D_new, i, s_size, xtx_i);
    }
    gsl_vector_set(mbeta_new, s_size, mbeta[pos]);
    return;
}

// set Xgamma for MCMC delete step
void BVSRM::SetSSgammaDel (const gsl_matrix *D_old, const gsl_vector *mbeta_old, const vector<size_t> &rank_old, size_t col_id, gsl_matrix *D_new, gsl_vector *mbeta_new)
{
    size_t s_size = rank_old.size();
    size_t s2;

    if (col_id==0) {
        s2 = s_size-1;
        gsl_matrix_const_view D22_sub = gsl_matrix_const_submatrix(D_old, 1, 1, s2, s2);
        gsl_vector_const_view mbeta2_sub = gsl_vector_const_subvector(mbeta_old, 1, s2);
        
        gsl_matrix_view D_new22_sub=gsl_matrix_submatrix(D_new, 0, 0, s2, s2);
        gsl_vector_view mbeta_new2_sub=gsl_vector_subvector(mbeta_new, 0, s2);
        
        gsl_matrix_memcpy(&D_new22_sub.matrix, &D22_sub.matrix);
        gsl_vector_memcpy(&mbeta_new2_sub.vector, &mbeta2_sub.vector);
    }
    else if(col_id == (s_size-1)){

        gsl_matrix_const_view D11_sub = gsl_matrix_const_submatrix(D_old, 0, 0, col_id, col_id);
        gsl_vector_const_view mbeta1_sub = gsl_vector_const_subvector(mbeta_old, 0, col_id);
        
        gsl_matrix_view D_new11_sub=gsl_matrix_submatrix(D_new, 0, 0, col_id, col_id);
        gsl_vector_view mbeta_new1_sub=gsl_vector_subvector(mbeta_new, 0, col_id);
        
        gsl_matrix_memcpy(&D_new11_sub.matrix, &D11_sub.matrix);
        gsl_vector_memcpy(&mbeta_new1_sub.vector, &mbeta1_sub.vector);
    }
    else{
        s2 = s_size - col_id - 1;
        
        gsl_matrix_const_view D11_sub = gsl_matrix_const_submatrix(D_old, 0, 0, col_id, col_id);
        gsl_matrix_const_view D12_sub = gsl_matrix_const_submatrix(D_old, 0, col_id+1, col_id, s2);
        gsl_matrix_const_view D21_sub = gsl_matrix_const_submatrix(D_old, col_id+1, 0, s2, col_id);
        gsl_matrix_const_view D22_sub = gsl_matrix_const_submatrix(D_old, col_id+1, col_id+1, s2, s2);
        
        gsl_vector_const_view mbeta1_sub = gsl_vector_const_subvector(mbeta_old, 0, col_id);
        gsl_vector_const_view mbeta2_sub = gsl_vector_const_subvector(mbeta_old, col_id+1, s2);

        gsl_matrix_view Dnew11_sub=gsl_matrix_submatrix(D_new, 0, 0, col_id, col_id);
        gsl_matrix_view Dnew12_sub=gsl_matrix_submatrix(D_new, 0, col_id, col_id, s2);
        gsl_matrix_view Dnew21_sub=gsl_matrix_submatrix(D_new, col_id, 0, s2, col_id);
        gsl_matrix_view Dnew22_sub=gsl_matrix_submatrix(D_new, col_id, col_id, s2, s2);
        
        gsl_vector_view mbeta_new1_sub=gsl_vector_subvector(mbeta_new, 0, col_id);
        gsl_vector_view mbeta_new2_sub=gsl_vector_subvector(mbeta_new, col_id, s2);
        
        gsl_matrix_memcpy(&Dnew11_sub.matrix, &D11_sub.matrix);
        gsl_matrix_memcpy(&Dnew12_sub.matrix, &D12_sub.matrix);
        gsl_matrix_memcpy(&Dnew21_sub.matrix, &D21_sub.matrix);
        gsl_matrix_memcpy(&Dnew22_sub.matrix, &D22_sub.matrix);
        
        gsl_vector_memcpy(&mbeta_new1_sub.vector, &mbeta1_sub.vector);
        gsl_vector_memcpy(&mbeta_new2_sub.vector, &mbeta2_sub.vector);
    }
    return;

}

// checked 06/15/2022
double BVSRM::CalcLR_cond_SS(const double &rtr, const size_t pos_j, const vector< vector<double> > &LD, const vector<double> &mbeta, const vector <size_t> &rank_cond, const gsl_vector *beta_cond, gsl_vector * Xtx_j)
{
    size_t pos_i;
    double lrt=0.0, xtx_ij, Xtxb_j, xtr_j;

    for(size_t i=0; i<rank_cond.size(); i++)
    {
        pos_i = mapRank2pos[rank_cond[i]];
        xtx_ij = getXtX(LD, pos_i, pos_j); // correlation
        gsl_vector_set(Xtx_j, i, xtx_ij);
    }

    gsl_blas_ddot(Xtx_j, beta_cond, &Xtxb_j);
    xtr_j = mbeta[pos_j] - Xtxb_j; // conditional genetic effect
    lrt = rtr - xtr_j * xtr_j ; // residual r2 in the conditional model
    if( lrt <= 0)
    {
        // cout << "rtr = " << rtr << "; xtr_j = " << xtr_j << "; lrt = " << lrt << endl;
        lrt = 0.0000001;
    }
    lrt = (double)(ni_test) * (log(rtr) - log(lrt) ); // Likelihood Ratio Test Statistic
    return lrt;
}



// Propose rank for switch step (switching pos) when X_cond not empty
// Based on conditional log likelihood test statistics
gsl_ran_discrete_t * BVSRM::MakeProposalSS(const vector< vector<double> > &LD, const vector <double> &mbeta, const size_t &pos, double *p_cond, const map<size_t, int> &mapRank2in, const gsl_vector * beta_cond, const double &rtr, const vector<size_t> rank_cond)
{

    size_t pos_j, rank_j;
    vector<bool> j_ind;
    double p_max, p_sum=0.0;

    gsl_vector * Xtx_j = gsl_vector_alloc(beta_cond->size); // save (X_cond)'x_j

    for (size_t j=0; j < ns_neib; ++j){

        pos_j = (pos - win) + j;
        if((pos_j >= 0) && (pos_j < ns_test))
        {
            rank_j = mapPos2Rank[pos_j];

            if( (j != win) && (mapRank2in.count(rank_j) == 0) )
            {
                p_cond[j]=CalcLR_cond_SS(rtr, pos_j, LD, mbeta, rank_cond, beta_cond, Xtx_j); //calc loglr
                j_ind.push_back(1); //propose candidates
            }
            else{
                p_cond[j] = -std::numeric_limits<double>::max();
                j_ind.push_back(0);
            }
        }else{
            p_cond[j] = -std::numeric_limits<double>::max();
            j_ind.push_back(0);
        }
    }
    p_max = *std::max_element(p_cond, p_cond+(ns_neib));

    for (size_t j=0; j < ns_neib; ++j){
        if(j_ind[j]){
            p_cond[j]=exp(p_cond[j]- p_max);
            p_sum += p_cond[j];
        }
        else{p_cond[j] = 0.0;}
    }
    //
    p_sum = 1.0/p_sum;
    for(size_t j=0; j < ns_neib; ++j){
        p_cond[j] *= p_sum;
    }

    gsl_vector_free(Xtx_j);

    return (gsl_ran_discrete_preproc(ns_neib, p_cond));
}

// Propose rank for switch step (switching pos) when X_cond is empty
gsl_ran_discrete_t * BVSRM::MakeProposalSS(const size_t &pos, double *p_cond, const map<size_t, int> &mapRank2in)
{

    size_t pos_j, rank_j;
    vector<bool> j_ind;
    double p_max, p_sum=0.0;

    for (size_t j=0; j < ns_neib; ++j){

        pos_j = (pos - win) + j;
        if((pos_j >= 0) && (pos_j < ns_test))
        {
            rank_j = mapPos2Rank[pos_j];

            if( (j != win) && (mapRank2in.count(rank_j) == 0) )
            {
                p_cond[j]= pos_ChisqTest[rank_j].second; // obtain loglr
                j_ind.push_back(1); //propose candidates
            }
            else{
                p_cond[j] = -std::numeric_limits<double>::max();
                j_ind.push_back(0);
            }
        }
        else{
            p_cond[j] = -std::numeric_limits<double>::max();
            j_ind.push_back(0);
        }
    }
    p_max = *std::max_element(p_cond, p_cond+(ns_neib));

    for (size_t j=0; j < ns_neib; ++j){
        if(j_ind[j]){
            p_cond[j]=exp(p_cond[j]- p_max);
            p_sum += p_cond[j];
        }
        else{p_cond[j] = 0.0;}
    }
    //
    p_sum = 1.0/p_sum;
    for(size_t j=0; j < ns_neib; ++j){
        p_cond[j] *= p_sum;
    }

    return (gsl_ran_discrete_preproc(ns_neib, p_cond));
}


// MCMC_SS function (Updated 06/15/2022 JY)
void BVSRM::MCMC_SS (const vector< vector<double> > &LD, const vector<double> &mbeta) {
    
    cout << "\nRunning MCMC with Summary Statistics: LD correlation matrix and marginal effect sizes. \n";
    cout << "Assuming the marginal effect sizes were obtained with standardized genotype and phenotype vectors." << endl;
    //cout << "Unique function types = " << n_type << endl;
    //cout << "ni_test = " << ni_test << endl;
    //cout << "ns_test = " << ns_test << endl;
    //cout << "snp_pos size = " << snp_pos.size() << endl;
    //cout << "ni_effect_vec size = " << ni_effect_vec.size() << endl;
    //cout << snp_pos[0].key <<"; " << snp_pos[2].key <<"; " << snp_pos[3].key <<"; " << endl;
    //cout << snp_pos[0].maf <<"; " << snp_pos[2].maf <<"; " << snp_pos[3].maf <<"; " << endl;

    // standardized phenotype with variance 1
    rv = 1.0;

    clock_t time_start;
    time_Proposal = 0.0;
    time_Omega = 0.0;

    //new model related
    gsl_vector *sigma_subvec_old = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_old);
    gsl_vector *sigma_subvec_new = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_new);
    gsl_vector *LnPost = gsl_vector_alloc(s_step); //save logPost...
    GV = 0.0;  // save regression R2
    
    // let D represent correlation matrix; mbeta
    gsl_matrix *D_old=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *mbeta_old=gsl_vector_alloc (s_max);
    gsl_vector *beta_old=gsl_vector_alloc (s_max);
    
    gsl_matrix *D_new=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *mbeta_new=gsl_vector_alloc (s_max);
    gsl_vector *beta_new=gsl_vector_alloc (s_max);

    //Initialize variables for MH
    double logPost_new=0.0, logPost_old=0.0, loglikegamma;
    double logMHratio;
    vector<size_t> rank_new, rank_old;
    class HYPBSLMM cHyp_old, cHyp_new;
    bool Error_Flag=0;

    vector<pair<double, double> > beta_g; //ysave beta estimates
    for (size_t i=0; i<ns_test; i++) {
        beta_g.push_back(make_pair(0.0, 0.0));
    }

    // UcharTable, vector<SNPPOS> snp_pos were, pos_ChisqTest created in CALCSS SS
    //cout << "Sort pos_ChisqTest, pval_vec by association evidence ... \n ";
    stable_sort (pos_ChisqTest.begin(), pos_ChisqTest.end(), comp_lr); // sort ChisqStat

    //cout << "Sort pval_vec by association evidence ... \n ";
    stable_sort (pval_vec.begin(), pval_vec.end()); // sort pval
    // PrintVector(pval_vec, 10);

    //cout << "Generate maps ... \n";
    size_t pos; // rank based on chisq test statistic is more stable than based on pvalue
    for (size_t i=0; i<ns_test; ++i) {
        mapRank2pos[i]=pos_ChisqTest[i].first;
        mapPos2Rank[pos_ChisqTest[i].first] = i;

        mapRank2Order[i]=pos_ChisqTest[i].first;
        mapOrder2Rank[pos_ChisqTest[i].first] = i;
    }

    //Calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
    //cout << "Calculate proposal distribution for gamma \n";
    gsl_rng_env_setup();
    const gsl_rng_type * gslType;
    gslType = gsl_rng_default;
    if (randseed<0)
    {
        time_t rawtime;
        time (&rawtime);
        tm * ptm = gmtime (&rawtime);
        randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
    }
    gsl_r = gsl_rng_alloc(gslType);
    gsl_rng_set(gsl_r, randseed);
    p_gamma = new double[ns_test]; // defined in bvsrm.h

    size_t p_gamma_top=0;
    for(size_t i=0; i < pval_vec.size(); i++){
        if (pval_vec[i] < 5e-8) { p_gamma_top++; }
        else{ break; }
    }
    cout << "Number of variants with p-value < 5e-8 : " << p_gamma_top << endl;

    // calculate discrete distribution for gamma
    // cout << "ns_test = " << ns_test;
    SetPgamma (p_gamma_top);
    gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma); // set up proposal function for gamma

    //Initial parameters
    cout << "\nStart initializing MCMC ... \n";
    InitialMCMC_SS (LD, rank_old, cHyp_old, pval_vec); // Initialize rank and cHyp
    if (cHyp_old.n_gamma > 0) SetSSgamma(LD, mbeta, rank_old, D_old, mbeta_old);
    cHyp_initial=cHyp_old;
    // Calculate likelihood with the top significant SNP
    loglikegamma = CalcLikegamma(rank_old);
    logPost_old = CalcPosterior_SS(D_old, mbeta_old, beta_old, cHyp_old, Error_Flag) + loglikegamma;
    if (Error_Flag) {
        cerr << "CalcPosterior_SS Failed at initialMCMC...\n";
        exit(-1);
    }
    cout <<  "Initial logPost_old = " << logPost_old << endl;

    //Start MCMC
    region_pip = 0.0;
    size_t k_save_sample=0;
    w_pace=2000;
    int accept; // accept_theta; naccept_theta=0,
    size_t total_step=w_step+s_step;
    cout << "total mcmc iterations = " << total_step << endl;
    size_t repeat=1;
    flag_gamma=0; // defined in bvsrm.h
    double accept_percent, betai; // accept_theta_percent;

    vector <string> snps_mcmc; // save locations of included snps per iteration
    string snps_mcmc_temp;

    for (size_t t=0; t<total_step; ++t) {
        
       if (t>0 && (t%d_pace==0 || t==total_step-1) ) {
                ProgressBar ("\nRunning MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));
                cout << endl;
            }

        for (size_t i=0; i<n_mh; ++i) {
            //cout << "\n \n propose gamam...\n";
            //cout << "old rank: "; PrintVector(rank_old);
            //cout <<"beta_old: "; PrintVector(beta_old, rank_old.size());
            //repeat = 1;
            cHyp_new = cHyp_old;
            rank_new = rank_old;
            time_start = clock();
            logMHratio = ProposeGamma_SS (rank_old, rank_new, cHyp_old, cHyp_new, repeat, LD, mbeta, D_old, mbeta_old, D_new, mbeta_new); //JY
           // cout << "After ProposeGamma_SS : propose gamma logMHratio = "<<logMHratio << "; MHratio = " << exp(logMHratio) << endl ;
            time_Proposal += (clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            //cout << "propose new rank: "; PrintVector(rank_new);
            //cout << "flag_gamma = " << flag_gamma << endl;
            //cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
            //cout << "propose gamma logMHratio = "<<logMHratio << "; MHratio = " << exp(logMHratio) << endl;

            accept = 0;
            if (flag_gamma > 0) {

                if(flag_gamma==1) nadd++;
                else if(flag_gamma==2) ndel++;
                else  nswitch++;

                loglikegamma = CalcLikegamma(rank_new);
                //   cout << "loglikegamma = " << loglikegamma << " in the non-Null model \n";
                if (rank_new.size() > 0) {
                    logPost_new = CalcPosterior_SS (D_new, mbeta_new, beta_new, cHyp_new, Error_Flag) + loglikegamma;
                 //   cout << "Logpos of the newly proposed non-Null model is " << logPost_new << endl;
                }
                else{
                    logPost_new = loglikegamma;
                    cHyp_new.pve=0.0;
                    //cout << "Logpos of the null model is " << logPost_new << endl;
                }
               //cout << "cHyp_old.m_gamma = "; PrintVector(cHyp_old.m_gamma);
               //cout << "cHyp_new.m_gamma = "; PrintVector(cHyp_new.m_gamma);
               //cout <<"beta_new: "; PrintVector(beta_new, rank_new.size());

                //  cout << "Calcposterior success." << endl;
                if (!Error_Flag) {
                    logMHratio += logPost_new-logPost_old;
                   // cout <<"logPost_old = " << logPost_old<< "; logPost_new = "<< logPost_new<< "\n logMHratio = " << logMHratio<< "; MHratio = " << exp(logMHratio) << endl;
                    if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio)
                        { accept=1;  n_accept++; }
                }
            }
            else if (flag_gamma == 0) { nother++; }
            time_Omega += (clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

            //cout << "accept = " << accept << endl;

            if (accept==1) {
                    if(flag_gamma==1) nadd_accept++;
                    else if(flag_gamma==2) ndel_accept++;
                    else if(flag_gamma==3) nswitch_accept++;
                    else nother_accept++;

                    logPost_old = logPost_new;
                    cHyp_old.pve = cHyp_new.pve;
                    cHyp_old.n_gamma = cHyp_new.n_gamma;
                    rank_old.clear();
                    for (size_t i=0; i<rank_new.size(); i++) {
                        rank_old.push_back(rank_new[i]); // copy rank_new to rank_old
                    }
                if (rank_old.size() != rank_new.size()) {
                    cerr << "Error: rank_old size != rank_new size\n";
                    exit(-1);
                }
                //cout << "Accept proposal: "; PrintVector(rank_old);

                if(rank_old.size()>0){
                    gsl_matrix_view Dold_sub=gsl_matrix_submatrix(D_old, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view mbeta_old_sub=gsl_vector_subvector(mbeta_old, 0, rank_new.size());
                    gsl_vector_view betaold_sub=gsl_vector_subvector(beta_old, 0, rank_new.size());
                    
                    gsl_matrix_view Dnew_sub=gsl_matrix_submatrix(D_new, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view mbeta_new_sub=gsl_vector_subvector(mbeta_new, 0, rank_new.size());
                    gsl_vector_view betanew_sub=gsl_vector_subvector(beta_new, 0, rank_new.size());
                    
                    gsl_matrix_memcpy(&Dold_sub.matrix, &Dnew_sub.matrix);
                    gsl_vector_memcpy(&mbeta_old_sub.vector, &mbeta_new_sub.vector);
                    gsl_vector_memcpy(&betaold_sub.vector, &betanew_sub.vector);
                }
            } //cout << "copy data from new propose -> old " << endl;
        } //end of n_mh

         //if (t % 10 == 0 && t > w_step) {
         if (t>0 && t % w_pace == 0 && t > w_step) {
             cout << "MCMC iteration after burn-ins = " << t - w_step << endl;
             cout << "From last mcmc: cHyp_old.n_gamma = " << cHyp_old.n_gamma << endl;
             cout << "From last mcmc: beta_hat: "; PrintVector(beta_old, rank_old.size());
             cout << "From last mcmc: posterior loglike = " << logPost_old << endl;
             accept_percent = (double)n_accept/(double)((t+1) * n_mh);
             cout << "acceptance percentage = " << setprecision(6) << accept_percent << endl ;
        }

        //Save data
        if (t<w_step) {continue;}
        else {
            //save loglikelihood to LnPost
            gsl_vector_set (LnPost, k_save_sample, logPost_old);
            GV += cHyp_old.pve;

            if (cHyp_old.n_gamma > 0){
                region_pip += 1.0; //count increase if the model has >0 SNPs
                snps_mcmc_temp="";
                for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
                    // beta_g saved by position
                    pos=mapRank2pos[rank_old[i]];
                    betai = gsl_vector_get(beta_old, i);
                    beta_g[pos].first += betai;
                    beta_g[pos].second += 1.0;
                    snps_mcmc_temp += string(snp_pos[pos].rs) + string(":") + string(snp_pos[pos].chr) + string(":") + to_string(snp_pos[pos].bp) + string(":") + to_string(snp_pos[pos].a_major) + string(":") + to_string(snp_pos[pos].a_minor)+ string(";");
                }
                snps_mcmc.push_back(snps_mcmc_temp);
                k_save_sample++;
              }
            }
    }

    cout<< "MCMC completed ... " << endl << endl;
    // cout << "Total Accepted mcmc steps = " << region_pip << "; s_step = " << s_step << endl;
    region_pip = region_pip / (double)s_step;
    cout << "region_pip = " << setprecision(5) << region_pip << endl;

    accept_percent = (double)n_accept/(double)(total_step * n_mh);
    cout << "Gamma Proposal Acceptance Percentage = " << accept_percent << endl ;
    //cout << "Last iteration beta_hat: "; PrintVector(beta_old, rank_old.size());
    //cout << "Last iteration Posterior loglike: " << logPost_old << endl;
    cout << "Accepted proposals: k_save_sample = " << k_save_sample << endl;
    // cout << "Write snps_mcmc ... \n";
    // cout << "snps_mcmc length: " << snps_mcmc.size() << endl;
    // WriteMCMC(snps_mcmc); // save all active SNPs from MCMC

    // Write paramtemp
    cout << "Write paramtemp ... \n";
    // calculate em_gamma; beta_mamc; sumbeta2;
    WriteParam_SS(beta_g, snp_pos, pos_ChisqTest, pval_vec, LD);

    //Save temp EM results
    cout << "Write hyptemp ... \n";
    WriteHyptemp(LnPost);

    gsl_vector_free(sigma_subvec_old);
    gsl_vector_free(sigma_subvec_new);
    gsl_vector_free(LnPost);
        
    gsl_matrix_free(D_old);
    gsl_vector_free(mbeta_old);
    gsl_vector_free(beta_old);
    
    gsl_matrix_free(D_new);
    gsl_vector_free(mbeta_new);
    gsl_vector_free(beta_new);

    delete [] p_gamma;
    // beta_g.clear();

    return;
}
// end of MCMC_SS version

// JY edit 06/15/2022
//InitialMCMC_SS with Summary Statistics
void BVSRM::InitialMCMC_SS (const vector< vector<double> > &LD, vector<size_t> &rank, class HYPBSLMM &cHyp, const vector<double> &pval)

{
    cout << "Start with Step-wise selected variants.\n";
        vector< pair<size_t, double> > rank_loglr;
        size_t pos_r, pos_j, radd, s_size;
        double rtr;

        double sig_lr = gsl_cdf_chisq_Qinv(5e-8, 1) ;
        // cout << "Genome-wide significant LRT is " << sig_lr << endl;

        size_t topMarkers = 500;
        //cout << "ns_test = " << ns_test << endl;
        if(topMarkers > ns_test) topMarkers = ns_test;
        //cout << "pos_ChisqTest size = " << pos_ChisqTest.size() << endl;

        rank.clear();
        rank.push_back(0);
        //cout << "Initial rank size " << rank.size() << endl;

        //excluded top SSNP
        for(size_t i=1; i < topMarkers; i++){
            rank_loglr.push_back( make_pair(i, pos_ChisqTest[i].second) );
        } // rank_loglr: pair of rank and ChisqTest

        //cout << "s_max = " << s_max << endl;
        gsl_matrix *D_cond=gsl_matrix_alloc (s_max, s_max);
        gsl_vector *mbeta_cond=gsl_vector_alloc (s_max);
        gsl_vector *Xtx_cond=gsl_vector_alloc (s_max);
        gsl_vector *beta_cond = gsl_vector_alloc (s_max);

        for(size_t i = 1; i < s_max; i++)
        {
            s_size = rank.size();
            SetSSgamma(LD, mbeta, rank, D_cond, mbeta_cond);
           // cout << "SetSSgamma with rank  " << endl; PrintVector(rank);
            //for(size_t l = 0; l < rank.size(); l++){
            //    cout << "bp = " << snp_pos[ mapRank2pos[rank[l]] ].bp << ";" ;
            //}
           // cout << "\n D_cond : \n" ; PrintMatrix(D_cond, s_size, s_size);
           // cout << "\n mbeta_cond : \n"; PrintVector(mbeta_cond, s_size);

            gsl_matrix_const_view D_cond_temp = gsl_matrix_const_submatrix(D_cond, 0, 0, s_size, s_size);
            gsl_vector_const_view mbeta_cond_temp = gsl_vector_const_subvector(mbeta_cond, 0, s_size);

            // calculate beta-hat
            gsl_vector_view beta_cond_temp = gsl_vector_subvector(beta_cond, 0, s_size);
            CalcBeta(&D_cond_temp.matrix, &mbeta_cond_temp.vector, &beta_cond_temp.vector);

            // calculate conditioned residual variance
            gsl_vector_const_view beta_cond_const = gsl_vector_const_subvector(beta_cond, 0, s_size);
            rtr = CalcResVar(&D_cond_temp.matrix, &beta_cond_const.vector);
          //  cout << "rtr = " << rtr << endl;

            gsl_vector_view Xtx_cond_temp = gsl_vector_subvector(Xtx_cond, 0, s_size);
            for(size_t j=0; j < rank_loglr.size(); j++){
                pos_j = mapRank2pos[ rank_loglr[j].first ];
                rank_loglr[j].second = CalcLR_cond_SS(rtr, pos_j, LD, mbeta, rank, &beta_cond_const.vector, &Xtx_cond_temp.vector);
            }
            stable_sort(rank_loglr.begin(), rank_loglr.end(), comp_lr); //sort conditional LRT statistics
            // cout << "Next top conditioned LRT is " << rank_loglr[0].second << endl;

            if(rank_loglr[0].second > sig_lr )
            {
                radd = rank_loglr[0].first;
                pos_r = mapRank2pos[radd];
                SetXtx(LD, rank, pos_r, &Xtx_cond_temp.vector); // BVSRM function

                gsl_vector_const_view Xtx_cond_const = gsl_vector_const_subvector(Xtx_cond, 0, s_size);
                if ( ColinearTest_SS(&D_cond_temp.matrix, &Xtx_cond_const.vector, &beta_cond_temp.vector) ) continue ;
                else{
                    rank.push_back(radd); // include rank r into initial model
                    rank_loglr.erase(rank_loglr.begin());
                }
            }
            else  {break;}

        }

   // cout << "\n D_cond : \n" ; PrintMatrix(D_cond, s_size, s_size);
   // cout << "\n mbeta_cond : \n"; PrintVector(mbeta_cond, s_size);
    gsl_matrix_free(D_cond);
    gsl_vector_free(mbeta_cond);
    gsl_vector_free(Xtx_cond);
    gsl_vector_free(beta_cond);

    cHyp.n_gamma = rank.size();
   // cout << "number of snps = " << cHyp.n_gamma << endl;
    cout << "Initial model with ranks: \n"; PrintVector(rank);
    cout << endl << endl;

    cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
    if (cHyp.logp==0) {cHyp.logp=-0.000001;}
    cHyp.h = 0.1;
    //if (cHyp.h<h_min) {cHyp.h=h_min;}
    //if (cHyp.h>h_max) {cHyp.h=h_max;}
    if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
    if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}

    //cout << "start setHyp... \n";
    double tau_beta_temp = 1.0;
    setHyp(tau_beta_temp);

    // cout<<"Prior causal prabobility per SNP = "; PrintVector(theta, 10);
    // cout<<"Initially selected number of variants in the model = "<<cHyp.n_gamma<<endl;
    return;
}

// Calculate posterior likelihood with summary statistics
double BVSRM::CalcPosterior_SS (const gsl_matrix *D, const gsl_vector *mbeta, gsl_vector *beta, class HYPBSLMM &cHyp, bool &Error_Flag)
{
    //conditioning on hyper parameters: subvar, log_theta
    double loglike=0.0;
    size_t s_size = cHyp.n_gamma;
    Error_Flag=0;
    
    gsl_matrix_const_view D_sub=gsl_matrix_const_submatrix (D, 0, 0, s_size, s_size);
    gsl_vector_const_view mbeta_sub=gsl_vector_const_subvector (mbeta, 0, s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);

    // Calculate Sigma_beta = (Omega)^(-1)
    gsl_matrix_memcpy(Omega, &D_sub.matrix);
    gsl_vector_view Omega_diag = gsl_matrix_diagonal(Omega);
    gsl_vector_add_constant(&Omega_diag.vector, tau_beta);

    // calculate logdet(Omega)
    double logdet_O = 0.0;
    logdet_O = LapackLogDet(Omega);
   // cout << "logdet_O = " << logdet_O << endl;

    // posterior estimates of beta_hat
    if(LapackSolve(Omega, &mbeta_sub.vector, beta_hat)!=0)
       EigenSolve(Omega, &mbeta_sub.vector, beta_hat);
    //cout << "beta_hat: "; PrintVector(beta_hat);
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
    gsl_vector_memcpy(&beta_sub.vector, beta_hat);
    double bSb;
    gsl_blas_ddot (&mbeta_sub.vector, beta_hat, &bSb);

    // Start here 06/04/2022
    gsl_vector *D_beta_hat = gsl_vector_alloc (s_size);
    gsl_blas_dgemv(CblasNoTrans, 1, &D_sub.matrix, beta_hat, 0, D_beta_hat);
    double R2;
    gsl_blas_ddot (D_beta_hat, beta_hat, &R2);
   // cout << "Regression R2 in CalcPosterior = " << R2 << endl;
     
    if (R2 > 1.0 ) {
        R2 = 1.0;
    }else if(R2 < 0.0){
        R2 = 0.0;
    }
    Error_Flag=0;
    cHyp.pve = R2; // Calculate pve

    loglike = 0.5 * ((double)s_size * log(tau_beta) +  (double)ni_test * bSb - logdet_O ); // log posterior likelihood
   // cout << "Posterior Loglike  = " << loglike << endl;

    gsl_matrix_free (Omega);
    gsl_vector_free (beta_hat);
    gsl_vector_free (D_beta_hat);

    return loglike;
}

bool BVSRM::CheckR2 (const vector<size_t> &rank, size_t ranki, const vector< vector<double> > &LD, const vector<double> &Xty, gsl_matrix *XtX_cond, gsl_vector *Xtx_cond, gsl_vector *Xty_cond, gsl_vector *beta_cond)
{
    bool colinear = 0;
    size_t s_size = rank.size();
    size_t pos_r;
    SetSSgamma(LD, Xty, rank, XtX_cond, Xty_cond);
    gsl_matrix_const_view XtX_cond_temp = gsl_matrix_const_submatrix(XtX_cond, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_cond_temp = gsl_vector_const_subvector(Xty_cond, 0, s_size);
    // calculate beta-hat
    gsl_vector_view beta_cond_temp = gsl_vector_subvector(beta_cond, 0, s_size);
    CalcBeta(&XtX_cond_temp.matrix, &Xty_cond_temp.vector, &beta_cond_temp.vector);

    pos_r = mapRank2pos[ranki];
    gsl_vector_view Xtx_cond_temp = gsl_vector_subvector(Xtx_cond, 0, s_size);
    SetXtx(LD, rank, pos_r, &Xtx_cond_temp.vector);
    gsl_vector_const_view Xtx_cond_const = gsl_vector_const_subvector(Xtx_cond, 0, s_size);
    if ( ColinearTest_SS(&XtX_cond_temp.matrix, &Xtx_cond_const.vector, &beta_cond_temp.vector) ) {
        colinear = 1;
    }
    return colinear;
}

// Propose new indicator ranks with summary statistics
double BVSRM::ProposeGamma_SS (const vector<size_t> &rank_old, vector<size_t> &rank_new, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, const vector< vector<double> > &LD, const vector<double> &mbeta, const gsl_matrix *D_old, const gsl_vector *mbeta_old, gsl_matrix *D_new, gsl_vector *mbeta_new)
{
    map<size_t, int> mapRank2in;
    double unif, logp = 0.0;
    size_t r_add, r_remove, col_id, r;
    if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"rank_old size wrong"<<endl;}
    rank_new.clear();
    if (rank_old.size() > 0) {
        for (size_t i=0; i<rank_old.size(); ++i) {
            r=rank_old[i];
            rank_new.push_back(r);
            mapRank2in[r]=1;
        }
    }
    cHyp_new.n_gamma=cHyp_old.n_gamma;
    //for (size_t i=0; i<repeat; ++i) {
        if(rank_old.size() == 0){
            flag_gamma=1; // add a SNP
        }
        else if(rank_old.size() == s_max){
            unif=gsl_rng_uniform(gsl_r);
            if (unif<=0.5) {flag_gamma=2;} // delete a SNP
            else if (unif>0.5) {flag_gamma=3;} // switch a SNP
        }
        else{
            unif=gsl_rng_uniform(gsl_r);
            if (unif < 0.33) {flag_gamma=1;}
            else if (unif>=0.33 && unif < 0.67) {flag_gamma=2;}
            else {flag_gamma=3;}
        }

        if(flag_gamma==1)  {//add a snp;
            //cout << "adding a snp ... \n" ;
            do {
                r_add=gsl_ran_discrete (gsl_r, gsl_t);
            } while ((mapRank2in.count(r_add)!=0));
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            mapRank2in[r_add]=1;
            rank_new.push_back(r_add);
            cHyp_new.n_gamma++;
            logp += -log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);
            if (rank_old.size()>0) {
                SetSSgammaAdd(LD, mbeta, D_old, mbeta_old, rank_old, r_add, D_new, mbeta_new);
            }
            else{ SetSSgamma (LD, mbeta, rank_new, D_new, mbeta_new); }
           //cout << "succesfully added a snp" << endl;
        }
        else if (flag_gamma==2) {//delete a snp;
            ///cout << "delete a snp" << endl;
            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            prob_total+=p_gamma[r_remove];
            mapRank2in.erase(r_remove);
            rank_new.erase(rank_new.begin()+col_id);
            logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
            cHyp_new.n_gamma--;
            if (rank_new.size() > 0) {
                SetSSgammaDel(D_old, mbeta_old, rank_old, col_id, D_new, mbeta_new);
            }
        }
        else if (flag_gamma==3) {//switch a snp;
            long int pos_add, pos_remove, pos_rj, pos_aj;
            size_t j_add, j_remove, o;
            double rtr;
            
            gsl_ran_discrete_t *gsl_s, *gsl_a; //JY added dynamic gsl_s
            double *p_cond_remove = new double[ns_neib];
            double *p_cond_add = new double[ns_neib];

            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma); //switch candidate
            r_remove=rank_new[col_id];//careful with the proposal
            pos_remove = mapRank2pos[r_remove];
            rank_new.erase(rank_new.begin()+col_id); // delete switch candidate
            size_t s_size = rank_new.size(); // conditioned SNPs
            mapRank2in.erase(r_remove);

            // conditional SS
            //cout << "Switch step rank_old:"; PrintVector(rank_old);
            //cout << "temp rank_new:"; PrintVector(rank_new);
            if (s_size > 0) {
                gsl_matrix *D_cond=gsl_matrix_alloc (s_size, s_size);
                gsl_vector *mbeta_cond=gsl_vector_alloc (s_size);
                gsl_vector *beta_cond = gsl_vector_alloc (s_size);
                SetSSgammaDel(D_old, mbeta_old, rank_old, col_id, D_cond, mbeta_cond);
                CalcBeta(D_cond, mbeta_cond, beta_cond); // from calcSS.cpp

                rtr = CalcResVar(D_cond, beta_cond); // residual variance from calcSS.cpp
                gsl_s = MakeProposalSS(LD, mbeta, pos_remove, p_cond_remove, mapRank2in, beta_cond, rtr, rank_new);

                j_add = gsl_ran_discrete(gsl_r, gsl_s);
                pos_add = (pos_remove - win) + j_add;
                while((pos_add < 0) || (pos_add >= (long int)ns_test) || (pos_add == pos_remove)){
                        cout << "j_add = " << j_add << "; pos_add = " << pos_add << endl;
                    j_add = gsl_ran_discrete(gsl_r, gsl_s);
                    pos_add = (pos_remove - win) + j_add;
                }
                r_add = mapPos2Rank[pos_add];
                gsl_a = MakeProposalSS(LD, mbeta, pos_add, p_cond_add, mapRank2in, beta_cond, rtr, rank_new);

                double prob_total_remove=1.0;
                double prob_total_add=1.0;
                for (size_t ii=0; ii<rank_new.size(); ++ii) {
                    r = rank_new[ii];
                    o = mapRank2pos[r];
                    pos_rj = ((long int)o - pos_remove) + win;
                    pos_aj = ((long int)o - pos_add) + win;
                    if(pos_aj >= 0 && pos_aj < (long int)ns_neib) prob_total_add -= p_cond_add [pos_aj];
                    if(pos_rj >= 0 && pos_rj < (long int)ns_neib) prob_total_remove -= p_cond_remove[pos_rj];
                }

                j_remove = pos_remove - pos_add + win;
                logp += log( p_cond_add[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                logp -= log( p_cond_remove[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)
                
                SetSSgammaAdd(LD, mbeta, D_cond, mbeta_cond, rank_new, r_add, D_new, mbeta_new);
                //cout << "XtX from setSSgammaAdd success: \n";
                // PrintMatrix(XtX_new, s_size+1, s_size+1);

                mapRank2in[r_add]=1;
                rank_new.push_back(r_add);
                
                gsl_matrix_free(D_cond);
                gsl_vector_free(mbeta_cond);
                gsl_vector_free(beta_cond);

            }
            else {
                //cout << "switch a snp 2" << endl;
                gsl_s = MakeProposalSS(pos_remove, p_cond_remove, mapRank2in);
                j_add = gsl_ran_discrete(gsl_r, gsl_s);
                pos_add = (pos_remove - win) + j_add;
                while((pos_add < 0) || (pos_add >= (long int)ns_test) || (pos_add == (long int)pos_remove)){
                    j_add = gsl_ran_discrete(gsl_r, gsl_s);
                    pos_add = (pos_remove - win) + j_add;
                }
                r_add = mapPos2Rank[pos_add];
                //cout << "o_add = " << pos_add <<  "; r_add = "<<r_add << endl;
                gsl_a = MakeProposalSS(pos_add, p_cond_add, mapRank2in);

                double prob_total_remove=1.0;
                double prob_total_add=1.0;
                for (size_t ii=0; ii<rank_new.size(); ++ii) {
                    r = rank_new[ii];
                    o = mapRank2pos[r];
                    pos_rj = ((long int)o - pos_remove) + win;
                    pos_aj = ((long int)o - pos_add) + win;
                    if(pos_aj >= 0 && pos_aj < (long int)ns_neib) prob_total_add -= p_cond_add[pos_aj];
                    if(pos_rj >= 0 && pos_rj < (long int)ns_neib) prob_total_remove -= p_cond_remove[pos_rj];
                }

                j_remove = pos_remove - pos_add + win;
                logp += log( p_cond_add[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
                logp -= log( p_cond_remove[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)

                mapRank2in[r_add]=1;
                rank_new.push_back(r_add);
                SetSSgamma (LD, mbeta, rank_new, D_new, mbeta_new);
            }

            gsl_ran_discrete_free(gsl_a);
            gsl_ran_discrete_free(gsl_s);

            delete[] p_cond_remove;
            delete[] p_cond_add;
            //cout << "successfully switched a snp" << endl;
        }

        else {logp+=0.0;}//do not change
    //}
    mapRank2in.clear();
    return logp;
}


void BVSRM::SetXtX(const vector< vector<double> > &LD, const vector<size_t> rank, gsl_matrix *XtX){
    double xtx_ij = 0.0;
    size_t pos_i, pos_j;
    for(size_t i=0; i<rank.size(); i++){
        pos_i = mapRank2pos[ rank[i] ] ;
        gsl_matrix_set(XtX, i, i, xtx_vec[pos_i]);
        for(size_t j= (i+1); j<rank.size(); j++){
            pos_j = mapRank2pos[ rank[j] ] ;
            xtx_ij = getXtX(LD, pos_i, pos_j);
            gsl_matrix_set(XtX, i, j, xtx_ij);
            gsl_matrix_set(XtX, j, i, xtx_ij);
        }
    }
    return ;
}

double BVSRM::Findmaxyty(const vector<size_t> &rank, const size_t s_size){
    double yty_max = yty;
    size_t pos_i;
    if(s_size > 0){
        for(size_t i = 0; i < s_size; i++){
            pos_i = mapRank2pos[ rank[i] ] ;
            if(yty_vec[pos_i] > yty_max)
                { yty_max = yty_vec[pos_i]; }
        }
    }
    return yty_max;
}


void BVSRM::SetXtx(const vector< vector<double> > &LD, const vector<size_t> rank, const size_t &pos_j, gsl_vector *Xtx_temp){
    double xtx_ij = 0.0;
    size_t pos_i;
    for(size_t i=0; i<rank.size(); i++){
        pos_i = mapRank2pos[ rank[i] ] ;
        xtx_ij = getXtX(LD, pos_i, pos_j);
        gsl_vector_set(Xtx_temp, i, xtx_ij);
    }
    return ;
}