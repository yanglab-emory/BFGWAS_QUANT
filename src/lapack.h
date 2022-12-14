/*
    Bayesian Functional GWAS --- MCMC (bfGWAS:MCMC)
    Copyright (C) 2016  Jingjing Yang

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

#ifndef __LAPACK_H__                
#define __LAPACK_H__



using namespace std;


void lapack_float_cholesky_decomp (gsl_matrix_float *A);
void lapack_cholesky_decomp (gsl_matrix *A);
void lapack_float_cholesky_solve (gsl_matrix_float *A, const gsl_vector_float *b, gsl_vector_float *x);
void lapack_cholesky_solve (gsl_matrix *A, const gsl_vector *b, gsl_vector *x);
void lapack_sgemm (char *TransA, char *TransB, float alpha, const gsl_matrix_float *A, const gsl_matrix_float *B, float beta, gsl_matrix_float *C);
void lapack_dgemm (char *TransA, char *TransB, double alpha, const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C);
void lapack_float_eigen_symmv (gsl_matrix_float *A, gsl_vector_float *eval, gsl_matrix_float *evec, const size_t flag_largematrix);
void lapack_eigen_symmv (gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec, const size_t flag_largematrix);

double EigenDecomp (gsl_matrix *G, gsl_matrix *U, gsl_vector *eval, const size_t flag_largematrix);
double EigenDecomp (gsl_matrix_float *G, gsl_matrix_float *U, gsl_vector_float *eval, const size_t flag_largematrix);

double CholeskySolve(gsl_matrix *Omega, const gsl_vector *Xty, gsl_vector *OiXty);
double CholeskySolve(gsl_matrix_float *Omega, gsl_vector_float *Xty, gsl_vector_float *OiXty);

void LUDecomp (gsl_matrix *LU, gsl_permutation *p, int *signum);
void LUDecomp (gsl_matrix_float *LU, gsl_permutation *p, int *signum);
void LUInvert (const gsl_matrix *LU, const gsl_permutation *p, gsl_matrix *inverse);
void LUInvert (const gsl_matrix_float *LU, const gsl_permutation *p, gsl_matrix_float *inverse);
double LULndet (gsl_matrix *LU);
double LULndet (gsl_matrix_float *LU);
void LUSolve (const gsl_matrix *LU, const gsl_permutation *p, const gsl_vector *b, gsl_vector *x);
void LUSolve (const gsl_matrix_float *LU, const gsl_permutation *p, const gsl_vector_float *b, gsl_vector_float *x);
void topdm(gsl_matrix *Omega);
void Ginv(gsl_matrix *XtX_gtemp);
void Ginv_logdet(gsl_matrix *Omega, double &logdet); // jy added 06/2022

void EigenSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta, const double &lambda);
void EigenInverse(gsl_matrix *XtX);
void CholeskyInverse(gsl_matrix *XtX);
void EigenSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta);
double CalcLogdet(const gsl_matrix *Omega);
void GSLSolve(const gsl_matrix *XtX, const gsl_vector *Xty, gsl_vector *beta_hat, const double lambda);
double LapackCholSolve(gsl_matrix *Omega, const gsl_vector *Xty, gsl_vector *OiXty);
double LapackLogDet(const gsl_matrix *Omega);
int LapackSolve(const gsl_matrix *A, const gsl_vector *b, gsl_vector *x);
//int LapackSolve(gsl_matrix *A, gsl_vector *b, gsl_vector *x);

#endif



