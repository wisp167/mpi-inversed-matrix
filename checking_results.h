#ifndef CHECKING_RESULTS_H_INCLUDED
#define CHECKING_RESULTS_H_INCLUDED
#pragma once
using namespace std;
void matr_mult_norm(int n, int m, double* a, double* b, double* temp1, double* pc);
double find_r(double* a, double* b, int n,int m, double* temp1, double* temp2);
double find_r_par(double* a, double* res, int kk, int p, int n, int m, double* block1, double* block2, double* block3, double* block4, double* temp1, double* temp2);
double norm1(double* a, int n);
void matrix_mult_vector(double* a, double* b, double* c);
void matrix_mult_vector1(double* a, double*b, double*c);
#endif // CHECKING_RESULTS_INCLUDED
