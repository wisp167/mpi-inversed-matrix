#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED
#pragma once
#define OMPI_SKIP_MPICXX 1
#include "mpi.h"
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/sysinfo.h>
using namespace std;
#include "print_matr.h"
#include "init_formula.h"
#include "init_matr.h"
#include "checking_results.h"
const double eps = 0.12*1e-15;
struct major{
	double mj;
	int mj_i, mj_j;
};
int MPI_f(MPI_Comm com, int n, int m, double &t1, double &t2, double &r1, double &r2, int s, int r, int kk, int p, string &filename, double *a, double*res, double *temp1, double *temp2, double* block1, double* block2, double* block3, double *block4, int* perm1, int* perm2, int* pperm1, int* pperm2);
void matr_tr(double* a, int n, int m, int p, int k, double*block, double* block1, MPI_Comm com);
double get_cpu_time();
double get_full_time();
int print_array(double* a, int n, int m, int printed_rows, int max_print);
int read_matrix(double* a, int n, int m, int p, int k, string &filename, double* buf, MPI_Comm com);
void print_matr(double* a, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com);
void print_matr_t(double* a, int n, int m, int p, int k, double* buf, double*block1, int max_print, MPI_Comm com);
void err_output(char* str);
void get_block(double* a, double* block, int n, int m, int i, int j);
void set_block(double *a, double* block, int n, int m, int i, int j);
void apply_perm(int n, int m, double *res, int* perm1, int* perm2);
void apply_perm_par(MPI_Comm com, int n, int m, int kk, int p, double *res, double* temp1, double* block, int* perm1, int* perm2);
void swap_blocks(int n, double* corner1, double* corner2, int n1, int m1);
void mult_by_minus_one(double* block, int n, int m);
double* get_corner(double* a, int n, int m, int i, int j);
void matr_swap_columns(int n, int m, double*a, int i, int j);
void print_matr_t(double* a, int n, int m, int p, int k, double *temp1, double*block, double &r1, double &r2, int r, MPI_Comm com);
void get_block_tr(double* a, double* block, int n, int m, int i, int j);
void matr_swap_rows(int n, int m, double* a, int i, int j);
bool major_element(int n, int m, int step, double &global_mx, double* a, int* perm1, int* perm2, int* pperm1, int* pperm2, double* block1);
void matr_sub(double *corner, double *b, int n, int n1, int m);
void matr_mult(double* c, double* a, double *b, int n1, int k, int n2);
void matr_swap_columns_par(int n, int m, int kk, int p, double*a, int i, int j);
void matr_swap_rows_par(MPI_Comm com, int n, int m, int kk, int p, double* a, double* temp1, double* block, int i, int j);
bool inversed_matr(int n, int m, double *a, double* block1, double* block2, double* block3, int* perm1, int* perm2, int* pperm1, int* pperm2);
bool inversed_matr_par(MPI_Comm com, int n, int m, int kk, int p, double *a, double* temp1, double *temp2, double* block1, double* block2, double* block3, int* perm1, int* perm2, int* pperm1, int* pperm2);
void matr_swap_rows_par(int n, int m, double* a, int i, int j, int k, int p);
double reduce_sum_mx(int p, double val, int i, int j, int &ans_i, int &ans_j, double* block, double* a, int n, int m, double& s);
double norm1_par(double* a, int n, int m, int kk, int p, double *temp1, double *temp2, MPI_Comm com);
#endif // HEADER_H_INCLUDED
