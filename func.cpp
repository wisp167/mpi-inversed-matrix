#include "header.h"

void err_output(char* str){
	printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", str, 20, (double)-1, (double)-1, (double)0, (double)0, 0, 0, 0, 0);
}
double get_cpu_time(){
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec+buf.ru_utime.tv_usec/1e6;
}
double get_full_time(){
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec+t.tv_usec/1e6;
}
inline int l2g(int m, int p, int k, int i_loc){
	int i_loc_m = i_loc/m;
	int i_glob_m = i_loc_m*p+k;
	return i_glob_m*m+i_loc%m;
}
inline int g2l(int m, int p, int i_glob){
	int i_glob_m= i_glob/m;
	int i_loc_m = i_glob_m/p;
	return i_loc_m*m+i_glob%m;
}
inline int get_max_rows(int n, int m, int p){
	int b = (n+m-1)/m;
	return (b+p-1)/p;
}
inline int get_rows(int n, int m, int p, int k){
	int b = (n+m-1)/m;
	return (b%p < k+1) ? b/p : b/p+1;
}
inline int get_rows_l(int n, int m, int p, int kk){
	int k = n/m;
	int l = n-m*k;
	int rows = get_rows(n, m, p, kk);
	return rows*m+(l && (k%p == kk))*(l-m);
}
inline int get_k(int m, int p, int i_glob){
	int i_glob_m= i_glob/m;
	return i_glob_m%p;
}
int print_array(double* a, int n, int m, int printed_rows, int max_print){
	if(printed_rows >= max_print) return 0;
	int p_n= (n > max_print ? max_print : n);
	int p_m = (printed_rows+m < max_print ? m : max_print-printed_rows);
	int i, j;
	for(i = 0; i < p_m; ++i){
		for(j = 0; j < p_n; ++j){
			printf(" %10.3e", a[i*n+j]);
		}
		printf("\n");
	}
	return p_m;
}
void print_matr(double* a, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com){
	int main_k = 0;
	int b, max_b = (n+m-1)/m;
	int printed_rows = 0;
	for(b = 0; b < max_b; ++b){
		int owner = b%p;
		int rows = min(m, n-b*m);
		int b_loc = b/p;
		if(k == main_k){
			if(owner == main_k){
				printed_rows+= print_array(a+b_loc*m*n, n, rows, printed_rows, max_print);
			}
			else{
				MPI_Status st;
				MPI_Recv(buf, n*rows, MPI_DOUBLE, owner, 0, com, &st);
				printed_rows += print_array(buf, n, rows, printed_rows, max_print);
			}
		}
		else if(k == owner){
			MPI_Send(a+b_loc*m*n, n*rows, MPI_DOUBLE, main_k, 0, com);
		}
	}
}
void print_matr_t(double* a, int n, int m, int p, int k, double* buf, double*block1, int max_print, MPI_Comm com){
	int main_k = 0;
	int b, max_b = (n+m-1)/m;
	int printed_rows = 0;
	int h;
	int l = n-(n/m)*m;
	for(b = 0; b < max_b; ++b){
		int row = (b*m+m <= n ? m : n-b*m);
		if(k == main_k){
			MPI_Status st;
			int counter =0;
			for(int i = 0; i < max_b; ++i){
				h = (l && (i == max_b-1) ? l : m);
				if(i%p == main_k){
					for(int ii = 0; ii < row; ++ii){
						for(int jj = 0; jj < h; ++jj){
							buf[i*m + ii*n + jj] = a[counter*m*n + b*m + ii + jj*n];
						}
					}
					++counter;
				}
				else{
					MPI_Recv(block1, h*row, MPI_DOUBLE, i%p, 0, com, &st);
					for(int ii = 0; ii < row; ++ii){
						for(int jj = 0; jj < h; ++jj){
							buf[i*m + ii*n + jj] = block1[ii+jj*row];
						} 
					}
				}
			}
			printed_rows += print_array(buf, n, row, printed_rows, max_print);
		}
		else{
			int te = get_rows(n, m, p, k);
			for(int i = 0; i < te; ++i){
				h = (l && (i == te-1) && (k == (max_b-1)%p) ? l : m);
				for(int ii = 0; ii < row; ++ii){
					for(int jj = 0; jj < h; ++jj){
						block1[ii+jj*row] = a[i*m*n + b*m + ii + jj*n];
					} 
				}
				MPI_Send(block1, h*row, MPI_DOUBLE, main_k, 0, com);
			}
		}
	}
}
int read_array(FILE *fp, double *a, int len){
	for(int i = 0; i < len; ++i){
		if(fscanf(fp, "%lf", a+i) != 1){
			return -2;
		}
	}
	return 0;
}
int read_array_t(FILE *fp, double *a, int n, int m, int row){
	for(int z = 0; z < n*row; ++z){
		int i = z/n;
		int j = z%n;
		if(fscanf(fp, "%lf", a+m*row*(j/m) + (j%m)*row + i) != 1){
			return -2;
		}
	}
	return 0;
}
int read_matrix(double* a, int n, int m, int p, int k, string &filename, double* buf, MPI_Comm com){
	int main_k = 0;
	FILE *fp = nullptr;
	int err = 0;
	if(k == main_k){
		fp = fopen(filename.c_str(), "r");
		if(fp == nullptr) err = 1;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if(err) return err;
	memset(buf, 0, n*m*sizeof(double));
	int max_b = (n+m-1)/m;
	int b;
	for(b = 0; b < max_b; ++b){
		int owner = b%p;
		int row = (b*m+m <= n ? m : n-b*m);
		int b_loc = b/p;
		if(k == main_k){
			err+=read_array(fp, buf, n*row);
			if(owner == main_k){
				memcpy(a+b_loc*m*n, buf, n*row*sizeof(double));
			}
			else{
				MPI_Send(buf, n*row, MPI_DOUBLE, owner, 0, com);
			}
		}
		else if(k == owner){
			MPI_Status st;
			MPI_Recv(a+b_loc*n*m, n*row, MPI_DOUBLE, main_k, 0, com, &st);
		}
	}
	if(k == main_k){
		fclose(fp);
		fp = nullptr;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if(err) return err;
	return 0;
}
int read_matrix_t(double* a, int n, int m, int p, int k, string &filename, double* buf, double* block1, MPI_Comm com){
	int main_k = 0;
	FILE *fp = nullptr;
	int err = 0;
	if(k == main_k){
		fp = fopen(filename.c_str(), "r");
		if(fp == nullptr) err = 1;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if(err) return err;
	memset(buf, 0, n*m*sizeof(double));
	int max_b = (n+m-1)/m;
	int b;
	int l = n-(n/m)*m;
	int h; int counter = 0;
	for(b = 0; b < max_b; ++b){
		//int owner = b%p;
		int row = (b*m+m <= n ? m : n-b*m);
		//int b_loc = b/p;
		if(k == main_k){
			err+=read_array_t(fp, buf, n, m, row); // {block, block, block...}
			counter = 0;
			for(int i = 0; i < max_b; ++i){
				h = (l && (i == max_b-1) ? l : m);
				if(i%p == main_k){
					for(int ii = 0; ii < row; ++ii){
						for(int jj = 0; jj < h; ++jj){
							a[counter*m*n + b*m + ii + jj*n] = buf[i*m*row + ii+jj*row];
						} 
					}
					++counter;
				}
				else{
					MPI_Send(buf + row*m*i, row*h, MPI_DOUBLE, i%p, 0, com);
				}
			}
		}
		else{
			MPI_Status st;
			int te = get_rows(n, m, p, k);
			for(int i = 0; i < te; ++i){
				h = (l && (i == te-1) && (k == (max_b-1)%p) ? l : m);
				MPI_Recv(block1, h*row, MPI_DOUBLE, main_k, 0, com, &st);
				for(int ii = 0; ii < row; ++ii){
					for(int jj = 0; jj < h; ++jj){
						a[i*m*n + b*m + ii + jj*n] = block1[ii+jj*row];
					} 
				}
			}
		}
	}
	if(k == main_k){
		fclose(fp);
		fp = nullptr;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if(err) return err;
	return 0;
}
void init_matrix(double* a, int n, int m, int p, int kk, int s, double (*f) (int, int, int, int)){
	int i_loc, j_loc, i_glob, j_glob,rows;
	rows = get_rows(n, m, p, kk);
	int k = n/m;
	int l = n-m*k;
	for(i_loc = 0; i_loc < rows*m+(l && (k%p == kk))*(l-m); i_loc++){
		i_glob = l2g(m, p, kk, i_loc);
		for(j_loc = 0; j_loc < n; ++j_loc){
			j_glob = j_loc;
			a[i_loc*n+j_loc] = (*f)(i_glob+1, j_glob+1, n, s);
		}
	}
}
int MPI_f(MPI_Comm com, int n, int m, double &t1, double &t2, double &r1, double &r2, int s, int r, int kk, int p, string& filename, double *a, double* res, double* temp1, double *temp2, double* block1, double* block2, double* block3, double *block4, int* perm1, int* perm2, int* pperm1, int* pperm2){
	if(!s){
		if(read_matrix_t(a, n, m, p, kk, filename, temp1, block1, com)){
			return -1;
        }
        if(read_matrix(res, n, m, p, kk, filename, temp1, com)){
			return -1;
        }
	}
	else{
		init_matrix(a, n, m, p, kk, s, *init_formula);
		init_matrix(res, n, m, p, kk, s, *init_formula);
	}
	print_matr_t(a, n, m, p, kk, temp1, block1, r, com);
	if(!kk){
		printf("\n");
	}
	swap(a, res);
	double t;
	MPI_Barrier(com);
	if(p == 1){
		t1 = MPI_Wtime();
		if(inversed_matr(n, m, res, block1, block2, block3, perm1, perm2, pperm1, pperm2)){
			return -1;
		}
		t1 = MPI_Wtime()-t1;
		for(int i  = 0; i < n; i++){
			for(int j = 0; j < i; ++j){
				swap(res[i*n+j], res[j*n+i]);
			}
		}
		t2 = MPI_Wtime();
		r1 = find_r(a, res, n, m, temp1, block1);
		r2 = find_r(res, a, n, m, temp1, block1);
		t2 = MPI_Wtime()-t2;
		for(int i  = 0; i < n; i++){
			for(int j = 0; j < i; ++j){
				swap(res[i*n+j], res[j*n+i]);
			}
		}
	}
	else{
		t1 = MPI_Wtime();
		if(inversed_matr_par(com, n, m, kk, p, res, temp1, temp2, block1, block2, block3, perm1, perm2, pperm1, pperm2)){
			return -1;
		}
		t1 = MPI_Wtime()-t1;
		t = MPI_Wtime();
		r1 = find_r_par(a, res, kk, p, n, m, block1, block2, block3, block4, temp1, temp2);
		t = MPI_Wtime()-t;
		matr_tr(a, n, m, p, kk, block1, block2, com);
		matr_tr(res, n, m, p, kk, block1, block2, com);
		t2 = MPI_Wtime();
		r2 = find_r_par(res, a, kk, p, n, m, block1, block2, block3, block4, temp1, temp2);
		t2 = MPI_Wtime()-t2+t;
		matr_tr(a, n, m, p, kk, block1, block2, com);
		matr_tr(res, n, m, p, kk, block1, block2, com);
	}
	return 0;
}

inline void get_block(double* a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((i < k) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				block[i*w+j] = corner[i*n+j];
		}
	}
}
inline void get_block_l(double* a, double* block, int n, int m, int i, int j){
        double* corner = get_corner(a, n, m, i, j);
        int k = n/m;
        int l = n-k*m;
        int h,w;
        w = ((j < k) ? m : l);
        h = ((l == 0) ? m : l);
        for(int i = 0; i < h; ++i){
                for(int j = 0; j < w; ++j){
                                block[i*w+j] = corner[i*n+j];
                }
        }
}
inline void get_block_tr(double* a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((i < k) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				block[j*h+i] = corner[i*n+j];
		}
	}
}
inline void get_block_tr_l(double* a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((l == 0) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				block[j*h+i] = corner[i*n+j];
		}
	}
}
inline void mult_by_minus_one(double* block, int n, int m){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			block[i*m+j] = -block[i*m+j];
		}
	}
}
inline double* get_corner(double* a, int n, int m, int i, int j){
	return a + i*n*m+j*m;
}
inline void set_block(double *a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((i < k) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				corner[i*n+j] = block[i*w+j];
		}
	}
}
inline void set_block_l(double *a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((l == 0) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				corner[i*n+j] = block[i*w+j];
		}
	}
}
inline void set_block_tr(double *a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((i < k) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				corner[i*n+j] = block[j*h+i];
		}
	}
}
inline void set_block_tr_l(double *a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((l == 0) ? m : l);
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < w; ++j){
				corner[i*n+j] = block[j*h+i];
		}
	}
}
inline void matr_sub(double *corner, double *b, int n, int n1, int m){
	for(int i = 0; i < n1; ++i){
		for(int j = 0; j < m; ++j){
			corner[i*n+j]-=b[i*m+j];
		}
	}
}
inline void matr_mult(double* pc, double* pa, double *pb, int n1, int k, int n2){ //to add rotations
	double ss;
	double* pt, *bt;
	for(int i = 0; i < n1; ++i){
		pt = pc+i*n2;
		for(int j = 0; j < n2; ++j){
			pt[j] = 0;
		}
		for(int j = 0; j < k; ++j){
			bt = pb + j*n2;
			ss = pa[i*k+j];
			for(int z= 0 ; z < n2; ++z){
				pt[z]+=ss*bt[z];
			}
		}
	}
	return;
	int r, t, q, v = n1, h = n2,v3, h3;
	double s, s00, s01, s02, s10, s20, s11, s22, s12, s21;
	for(r = 0; r < v; ++r){
		for(t = 0; t < h; ++t){
			pc[r*h+t] = 0;
		}
	}
	v3 = v%3; h3 = h%3;
	for(r = 0; r < v3; ++r){
		for(t = 0; t < h3; ++t){
			s = 0;
			for(q= 0; q < k; ++q){
				s+=pa[r*k+q]*pb[q*n2+t];
			}
			pc[r*h+t]+=s;
		}
		for(; t < h; t+=3){
			s00 = 0; s01 = 0; s02 = 0;
			for(q = 0; q < k; ++q){
				s00+=pa[r*k+q]*pb[q*n2+t];
				s01+=pa[r*k+q]*pb[q*n2+t+1];
				s02+=pa[r*k+q]*pb[q*n2+t+2];
			}
			pc[r*h+t]+=s00;
			pc[r*h+t+1]+=s01;
			pc[r*h+t+2]+=s02;
		}
	}
	for(;r<v;r+=3){
		for(t = 0; t<h3; ++t){
			s00 = s10 = s20 = 0;
			for(q = 0; q < k; q++){
				s00+=pa[r*k+q]*pb[q*n2+t];
				s10+=pa[(r+1)*k+q]*pb[q*n2+t];
				s20+=pa[(r+2)*k+q]*pb[q*n2+t];
			}
			pc[r*h+t]+=s00;pc[(r+1)*h+t]+=s10;pc[(r+2)*h+t]+=s20;
		}
		for(;t<h;t+=3){
			s01=s02=s10=s20=s11=s12=s21=s22=s00=0;
			for(q=0;q<k;q++){
				s00+=pa[r*k+q]*pb[q*n2+t];
				s01+=pa[r*k+q]*pb[q*n2+t+1];
				s02+=pa[r*k+q]*pb[q*n2+t+2];
				s10+=pa[(r+1)*k+q]*pb[q*n2+t];
				s11+=pa[(r+1)*k+q]*pb[q*n2+t+1];
				s12+=pa[(r+1)*k+q]*pb[q*n2+t+2];
				s20+=pa[(r+2)*k+q]*pb[q*n2+t];
				s21+=pa[(r+2)*k+q]*pb[q*n2+t+1];
				s22+=pa[(r+2)*k+q]*pb[q*n2+t+2];
			}
			pc[r*h+t]+=s00;
			pc[r*h+t+1]+=s01;
			pc[r*h+t+2]+=s02;
			pc[(r+1)*h+t]+=s10;
			pc[(r+1)*h+t+1]+=s11;
			pc[(r+1)*h+t+2]+=s12;
			pc[(r+2)*h+t]+=s20;
			pc[(r+2)*h+t+1]+=s21;
			pc[(r+2)*h+t+2]+=s22;
		}
	}
}
inline void swap_blocks(int n, double* corner1, double* corner2, int n1, int m1){
	double temp;
	for(int i = 0; i < n1; ++i){
		for(int j = 0; j < m1; ++j){
			temp = corner2[i*n+j];
			corner2[i*n+j] = corner1[i*n+j];
			corner1[i*n+j] = temp;
		}
	}
}
void print_matr_t(double* a, int n, int m, int p, int k, double *temp1, double*block, double &r1, double &r2, int r, MPI_Comm com){
	print_matr_t(a, n, m, p, k, temp1, block, r, com);
	double l = (double)9/10;
	double t= r1;r1*=1-min(max(r1, r2),(double)1)*l;r2*=1-min(max(t, r2),(double)1)*l;
	return;
}
void matr_swap_rows_par(MPI_Comm com, int n, int m, int kk, int p, double* a, double*temp1, double* block, int i, int j){
	if((i == j) || ((i%p != kk) && (j%p != kk))){
		return;
	}
	int t;
	int ii = g2l(1, p, i), jj = g2l(1, p, j);
	//int rows = get_rows(n, m, p, kk);
	//int rrows = (get_rows_l(n, m ,p, kk)%m != 0);
	if(m == 1){
		if(i % p == j%p){
			for(int z = 0; z < n; z++){
				swap(a[ii*n+z], a[jj*n+z]);
			}
		}
		else{
			MPI_Status st;
			if(kk == i%p){
				MPI_Send(a+n*ii, n, MPI_DOUBLE, j%p, 0, com);
				MPI_Recv(temp1, n, MPI_DOUBLE, j%p, 0, com, &st);
				t = ii;
				for(int z= 0; z < n; ++z){
					a[n*t+z] = temp1[z];
				}
			}
			else{
				MPI_Recv(temp1, n, MPI_DOUBLE, i%p, 0, com, &st);
				MPI_Send(a+n*jj, n, MPI_DOUBLE, i%p, 0, com);
				t = jj;
				for(int z= 0; z < n; ++z){
					a[n*t+z] = temp1[z];
				}
			}
		}
	}
	else{
		int k = n/m;
		int l = n-k*m;
		double *corner1, *corner2;
		if(i % p == j%p){
			for(int z = 0; z < k+(l!=0); z++){
				corner1 = get_corner(a, n, m, jj, z);
				corner2 = get_corner(a, n, m, ii, z);
				swap_blocks(n, corner1, corner2, m, ((z == k) ? l : m));
			}
		}
		else{
			MPI_Status st;
			if(kk == i%p){
				MPI_Send(a+n*m*ii, n*m, MPI_DOUBLE, j%p, 0, com);
				MPI_Recv(temp1, n*m, MPI_DOUBLE, j%p, 0, com, &st);
				t = ii;
				for(int z= 0; z < k+(l!=0); ++z){
					get_block(temp1, block, n, m, 0, z);
					set_block(a, block, n, m, t, z);
				}
			}
			else{
				MPI_Recv(temp1, n*m, MPI_DOUBLE, i%p, 0, com, &st);
				MPI_Send(a+n*m*jj, n*m, MPI_DOUBLE, i%p, 0, com);
				t = jj;
				for(int z= 0; z < k+(l!=0); ++z){
					get_block(temp1, block, n, m, 0, z);
					set_block(a, block, n, m, t, z);
				}
			}
		}
	}
	return;
}
void matr_swap_rows(int n, int m, double* a, int i, int j){
	double t;
	if(m == 1){
		for(int z = 0; z < n; z++){
			t = a[j*n+z];
			a[j*n+z] = a[i*n+z];
			a[i*n+z] = t;
		}
	}
	else{
		int k = n/m;
		int l = n-k*m;
		double *corner1, *corner2;
		for(int z = 0; z < k+(l!=0); z++){
			corner1 = get_corner(a, n, m, j, z);
			corner2 = get_corner(a, n, m, i, z);
			swap_blocks(n, corner1, corner2, m, ((z == k) ? l : m));
		}
	}
	return;
}
void matr_swap_columns_par(int n, int m, int kk, int p, double*a, int i, int j){
	if(i == j){
		return;
	}
	int rows = get_rows(n, m, p, kk);
	int rrows = (get_rows_l(n, m, p, kk)%m != 0);
	if(m == 1){
		for(int z = 0; z < rows; z++){
			swap(a[z*n+i], a[z*n+j]);
		}
	}
	else{
		int k = n/m;
		int l =n-k*m;
		double *corner1, *corner2;
		for(int z = 0; z < rows; z++){
			corner1 = get_corner(a, n, m, z, i);
			corner2 = get_corner(a, n, m, z, j);
			swap_blocks(n, corner1, corner2, ((rrows && (z == rows-1)) ? l : m), m);
		}
	}
	return;
}
void matr_swap_columns(int n, int m, double*a, int i, int j){
	if(m == 1){
		double t;
		for(int z = 0; z < n; ++z){
			t = a[z*n+j];
			a[z*n+j] = a[z*n+i];
			a[z*n+i] = t;
		}
	}
	else{
		int k = n/m;
		int l =n-k*m;
		double *corner1, *corner2;
		for(int z = 0; z < k+(l!=0); ++z){
			corner1 = get_corner(a, n, m, z, i);
			corner2 = get_corner(a, n, m, z, j);
			swap_blocks(n, corner1, corner2, ((z == k) ? l : m), m);
		}
	}
	return;
}
void apply_perm(int n, int m, double *res, int* perm1, int* perm2){
	if(m == 1){
		for(int i = n-1; i >= 0; --i){
			matr_swap_rows(n, 1, res, i, perm2[i]);
			matr_swap_columns(n, 1, res, i, perm1[i]);
		}
	}
	else{
		int k = n/m;
		for(int i = k-1; i >= 0; --i){
			matr_swap_rows(n, m, res, i, perm2[i]);
			matr_swap_columns(n, m, res, i, perm1[i]);
		}
	}
	return;
}
void apply_perm_par(MPI_Comm com, int n, int m, int kk, int p, double *res, double* temp1, double* block, int* perm1, int* perm2){
	if(m == 1){
		for(int i = n-1; i >= 0; --i){
			matr_swap_rows_par(com, n, 1, kk, p, res, temp1, block, i, perm2[i]);
			matr_swap_columns_par(n, 1, kk, p, res, i, perm1[i]);
		}
	}
	else{
		int k = n/m;
		for(int i = k-1; i >= 0; --i){
			matr_swap_rows_par(com, n, m, kk, p, res, temp1, block, i, perm2[i]);
			matr_swap_columns_par(n, m, kk, p, res, i, perm1[i]);
		}
	}
	return;
}
inline int start(int q, int k, int p, int n){
	if(q%p > k){
		return min(n, q-q%p+k+p);
	}
	return min(n, q-q%p+k);
}
void matr_tr(double* a, int n, int m, int p, int kk, double*block, double *block1, MPI_Comm com){
	int k = n/m;
	int l = n-k*m;
	int rows = get_rows(n, m, p, kk);
	int rrows = (get_rows_l(n, m ,p, kk)%m !=0);
	int ri, ii;
	int sz = k+(l!=0);
	MPI_Status st;
	//MPI_Request rq;
	//p = min(p, sz);
	int ans = 0;
	int ans1 = 0;
	for(int i = 0; i < sz-1; ++i){
		if(i%p == kk){
			ii = g2l(1, p, i);
			for(int j = i; j < sz; ++j){
				ri = (l && (j == sz-1));
				get_block(a, block, n, m, ii, j);
				if(j%p!=kk){
					MPI_Send(block, m*m, MPI_DOUBLE, j%p, 0, com);
					MPI_Recv(block, m*m, MPI_DOUBLE, j%p, 0, com, &st);
					set_block_tr(a, block, n, m, ii, j);
				}
				else{
					int jj = g2l(1, p, j);
					ri = (rrows && (jj == rows-1));
					if(!ri){
						get_block(a, block1, n, m, jj, i);
						set_block_tr(a, block, n, m, jj, i);
					}
					else{
						get_block_l(a, block1, n, m, jj, i);
						set_block_tr_l(a, block, n, m, jj, i);
					}
					set_block_tr(a, block1, n, m, ii, j);
				}
			}
		}
		else if(start(i, kk, p, sz) < sz){
			for(int j = g2l(1, p, start(i, kk, p, sz)); j < rows; ++j){
				ri = (rrows && (j == rows-1));
				if(!ri){
					get_block(a, block, n, m, j, i);
				}
				else{
					get_block_l(a, block, n, m, j, i);
				}
				if(i%p!=kk){
					MPI_Recv(block1, m*m, MPI_DOUBLE, i%p, 0, com, &st);
					MPI_Send(block, m*m, MPI_DOUBLE, i%p, 0, com);
				}
				if(!ri){
					set_block_tr(a, block1, n, m, j, i);
				}
				else{
					set_block_tr_l(a, block1, n, m, j, i);
				}
			}
		}
	}
	if(kk == ((sz-1)%p)){
		if(l){
			get_block_l(a, block, n, m, rows-1, sz-1);
			set_block_tr_l(a, block, n, m, rows-1, sz-1);
		}
		else{
			get_block(a, block, n, m, rows-1, sz-1);
			set_block_tr(a, block, n, m, rows-1, sz-1);
		}
	}
	MPI_Allreduce(&ans, &ans1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}
void Max(void *in, void *out, int* len, MPI_Datatype *type){
	int *l;
	MPI_Datatype *m = nullptr;
	swap(len, l); swap(len, l);
	swap(type, m); swap(type, m);
	major *inn = (major*) in;
	major *outt = (major*) out;
	if(inn->mj < outt->mj){
		return;
	}
	if((inn->mj > outt->mj) || (inn->mj_j < outt->mj_j)){
		*outt = *inn;
	}
}










bool major_element(int n, int m, int step, double &global_mx, double* a, int* perm1, int* perm2, int* pperm1, int* pperm2, double* block1){
	double mj;
	int mj_i = step;
	int mj_j = step;
	if(m == 1){
		mj = fabs(a[step*n+step]);
		for(int i = step; i < n; ++i){
			for(int j = step; j < n; ++j){
				if(fabs(a[i*n+j]) > mj){
					mj = fabs(a[i*n+j]);
					mj_i = i;
					mj_j = j;
				}
			}
		}
		if(step == 0){
			if(!(mj > 0)){return false;}
		}
		if(fabs(mj) < eps*global_mx){
			return false;
		}
	}
	else{
		int k = n/m;
		//int l = n-k*m;
		double s;
		int i = step, j = step;
		bool stat = true;
		mj = -1;
		for(int i1=step; i1 < k; ++i1){
			for(j =step; j < k; ++j){
				i = step;
				get_block(a, block1, n, m, i, j);
				s = norm1(block1, m);
				if(!(s < global_mx*eps*2) && !inversed_matr(m, 1, block1, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
					stat = false;
					s = norm1(block1, m);
					if(mj > 0){
						if(mj >= s){
							mj = s;
							mj_i = i;
							mj_j = j;
						}
					}
					else{
						mj = s;
						mj_i = i;
						mj_j = j;
					}
				}
			}
		}
		if(stat){
			return false;
		}
	}
	perm1[step] = mj_i;
	perm2[step] = mj_j;
	if(mj_i != step){
		matr_swap_rows(n, m, a, mj_i, step);
	}
	if(mj_j != step){
		matr_swap_columns(n, m, a, mj_j, step);
	}
	return true;
}


bool inversed_matr(int n, int m, double *a, double* block1, double* block2, double* block3, int* perm1, int* perm2, int* pperm1, int* pperm2){
	if((m < 1) || (m > n)){
		return true;
	}
	int k = n/m;
	int l = n-k*m;
	if(perm1 != nullptr){
		for(int i = 0; i < k+(l != 0); ++i){
			perm1[i] = i;
			perm2[i] = i;
		}
	}
	double* corner;
	int i, j;
	double global_mx = norm1(a, n);
	if(m == 1){
		double s;
		for(int q = 0; q < n; ++q){
			if(!major_element(n, 1, q, global_mx, a, perm1, perm2, nullptr, nullptr, nullptr)){
				return true;
			}
			s = 1/a[q*n+q];
			a[q*n+q]=s;
			for(i = 0; i < q; ++i){
				a[i*n+q]*=s;
			}
			for(i = q+1; i < n; ++i){
				a[i*n+q]*=s;
			}	
			for(i = 0; i < q; ++i){
				for(j = 0; j < q; ++j){
					a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < n; ++j){
					a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
			}
			for(i = q+1; i < n; ++i){
				for(j = 0; j < q; ++j){
					a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < n; ++j){
					a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
			}
			for(i = 0; i < q; ++i){
				a[q*n+i]*=-s;
			}
			for(int i = q+1; i < n; ++i){
				a[q*n+i]*=-s;
			}
		}
		apply_perm(n, 1, a, perm1, perm2);
		return false;
	}
	if(l == 0){
		for(int q= 0;  q < k; ++q){
			if(!major_element(n, m, q, global_mx, a, perm1, perm2, pperm1, pperm2, block1)){
				return true;
			}
			get_block(a, block2, n, m, q, q);
			if(inversed_matr(m, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
				return true;
			}
			set_block(a, block2, n, m, q, q);
			for(i = 0; i < q; ++i){
				get_block(a, block1, n, m, i, q);
				matr_mult(block3, block1, block2, m, m, m);
				set_block(a, block3, n, m, i, q);
				//a[q*n+i]*=s;
			}
			for(i = q+1; i < k; ++i){
				get_block(a, block1, n, m, i, q);
				matr_mult(block3, block1, block2, m, m, m);
				set_block(a, block3, n, m, i, q);
				//a[q*n+i]*=s;
			}
			for(i = 0; i < q; ++i){
				for(int j = 0; j < q; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < k; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
			}
			for(i = q+1; i < k; ++i){
				for(j = 0; j < q; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < k; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
			}
			get_block(a, block2, n, m, q, q);
			mult_by_minus_one(block2, m, m);
			for(i = 0; i < q; ++i){
				get_block(a, block1, n, m, q, i);
				matr_mult(block3, block2, block1, m, m, m);
				set_block(a, block3, n, m, q, i);
				//a[i*n+q]*=-a[q*n+q];
			}
			for(i = q+1; i < k; ++i){
				get_block(a, block1, n, m, q, i);
				matr_mult(block3, block2, block1, m, m, m);
				set_block(a, block3, n, m, q, i);
				//a[i*n+q]*=-a[q*n+q];
			}
		}
		apply_perm(n, m, a, perm1, perm2);
		return false;
	}
	else{
		for(int q= 0;  q < k; ++q){
			if(!major_element(n, m, q, global_mx, a, perm1, perm2, pperm1, pperm2, block1)){
				return true;
			}
			get_block(a, block2, n, m, q, q);
			if(inversed_matr(m, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
				return true;
			}
			set_block(a, block2, n, m, q, q);
			for(i = 0; i < q; ++i){
				get_block(a, block1, n, m, i, q);
				matr_mult(block3, block1, block2, m, m, m);
				set_block(a, block3, n, m, i, q);
				//a[q*n+i]*=s;
			}
			for(i = q+1; i < k; ++i){
				get_block(a, block1, n, m, i, q);
				matr_mult(block3, block1, block2, m, m, m);
				set_block(a, block3, n, m, i, q);
				//a[q*n+i]*=s;
			}
			get_block(a, block1, n, m, k, q);
			matr_mult(block3, block1, block2, l, m, m);
			set_block(a, block3, n, m, k, q);
			for(i = 0; i < q; ++i){
				for(int j = 0; j < q; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < k; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				get_block(a, block1, n, m, i, q);
				get_block(a, block2, n, m, q, k);
				matr_mult(block3, block1, block2, m, m, l);
				corner = get_corner(a, n, m, i, k);
				matr_sub(corner, block3, n, m, l);
			}
			for(i = q+1; i < k; ++i){
				for(j = 0; j < q; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				for(j = q+1; j < k; ++j){
					get_block(a, block1, n, m, i, q);
					get_block(a, block2, n, m, q, j);
					matr_mult(block3, block1, block2, m, m, m);
					corner = get_corner(a, n, m, i, j);
					matr_sub(corner, block3, n, m, m);
					//a[i*n+j]-=a[i*n+q]*a[q*n+j];
				}
				get_block(a, block1, n, m, i, q);
				get_block(a, block2, n, m, q, k);
				matr_mult(block3, block1, block2, m, m, l);
				corner = get_corner(a, n, m, i, k);
				matr_sub(corner, block3, n, m, l);
			}
			for(j = 0; j < q; ++j){
				get_block(a, block1, n, m, k, q);
				get_block(a, block2, n, m, q, j);
				matr_mult(block3, block1, block2, l, m, m);
				corner = get_corner(a, n, m, k, j);
				matr_sub(corner, block3, n, l, m);
				//a[i*n+j]-=a[i*n+q]*a[q*n+j];
			}
			for(j = q+1; j < k; ++j){
				get_block(a, block1, n, m, k, q);
				get_block(a, block2, n, m, q, j);
				matr_mult(block3, block1, block2, l, m, m);
				corner = get_corner(a, n, m, k, j);
				matr_sub(corner, block3, n, l, m);
				//a[i*n+j]-=a[i*n+q]*a[q*n+j];
			}
			get_block(a, block1, n, m, k, q);
			get_block(a, block2, n, m, q, k);
			matr_mult(block3, block1, block2, l, m, l);
			corner = get_corner(a, n, m, k, k);
			matr_sub(corner, block3, n, l, l);
			//a[i*n+j]-=a[i*n+q]*a[q*n+j];
			get_block(a, block2, n, m, q, q);
			for(i = 0; i < q; ++i){
				get_block(a, block1, n, m, q, i);
				mult_by_minus_one(block1, m, m);
				matr_mult(block3, block2, block1, m, m, m);
				set_block(a, block3, n, m, q, i);
				//a[i*n+q]*=-a[q*n+q];
			}
			for(i = q+1; i < k; ++i){
				get_block(a, block1, n, m, q, i);
				mult_by_minus_one(block1, m, m);
				matr_mult(block3, block2, block1, m, m, m);
				set_block(a, block3, n, m, q, i);
				//a[i*n+q]*=-a[q*n+q];
			}
			get_block(a, block1, n, m, q, k);
			mult_by_minus_one(block1, m, l);
			matr_mult(block3, block2, block1, m, m, l);
			set_block(a, block3, n, m, q, k);
			//a[i*n+q]*=-a[q*n+q];
		}
		//       q = k
		int q = k;
		get_block(a, block2, n, m, q, q);
		if((l == 1) && (norm1(block2, l) < eps*global_mx)){
			return true;
		}
		if(inversed_matr(l, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
			return true;
		}
		set_block(a, block2, n, m, q, q);
		for(i = 0; i < k; ++i){
			get_block(a, block1, n, m, i, q);
			matr_mult(block3, block1, block2, m, l, l);
			set_block(a, block3, n, m, i, q);
			//a[q*n+i]*=s;
		}
		for(i = 0; i < k; ++i){
			for(int j = 0; j < k; ++j){
				get_block(a, block1, n, m, i, q);
				get_block(a, block2, n, m, q, j);
				matr_mult(block3, block1, block2, m, l, m);
				corner = get_corner(a, n, m, i, j);
				matr_sub(corner, block3, n, m, m);
				//a[i*n+j]-=a[i*n+q]*a[q*n+j];
			}
		}
		get_block(a, block2, n, m, q, q);
		for(i = 0; i < k; ++i){
			get_block(a, block1, n, m, q, i);
			mult_by_minus_one(block1, l, m);
			matr_mult(block3, block2, block1, l, l, m);
			set_block(a, block3, n, m, q, i);
			//a[i*n+q]*=-a[q*n+q];
		}
		apply_perm(n, m, a, perm1, perm2);
		return false;
	}
}




















bool major_element_par(MPI_Comm com, int n, int m, int kk, int p, int step, double &global_mx, double* a, double* temp1, int* perm1, int* perm2, int* pperm1, int* pperm2, double* block1){
	double mj;
	int mj_i = step;
	int mj_j = step;
	int rows = get_rows(n, m, p, kk);
	MPI_Datatype ctype;
	MPI_Datatype typesig[2] = {MPI_DOUBLE, MPI_INT};
	int lens[2] = {1, 2};
	MPI_Aint displacements[2] = {0, 8};
	MPI_Type_create_struct(2, lens, displacements, typesig, &ctype);
    MPI_Type_commit(&ctype);
	major maj, ans;
	MPI_Op my_MAX;
    MPI_Op_create(Max, 1, &my_MAX);
    int st;
	if(m == 1){
		mj = 0;
		st = start(step, kk, p, n);
		if(st < n){
			for(int i = g2l(1, p, st); i < rows; ++i){
				for(int j = step; j < n; ++j){
					if(fabs(a[i*n+j]) > mj){
						mj = fabs(a[i*n+j]);
						mj_i = l2g(m, p, kk, i);
						mj_j = j;
					}
				}
			}
		}
		maj.mj = mj; maj.mj_i = mj_i; maj.mj_j = mj_j;
		MPI_Allreduce(&maj, &ans, 1, ctype, my_MAX, com);
		MPI_Op_free(&my_MAX);
		swap(ans, maj);
		mj_i = maj.mj_i; mj_j = maj.mj_j; mj = maj.mj;
		if(step == 0){
			if(!(mj > 0)){return false;}
		}
		if(fabs(mj) < eps*global_mx){
			return false;
		}
		/*
		MPI_Status st;
		if(!(mj > ans)){
			stat = kk;
			for(int i = 0; i < p; ++i){
				if(i != kk){
					MPI_Send(&stat, 1, MPI_INT, i, 0, com);
				}
			}
		}
		else{
			MPI_Recv(&stat, 1, MPI_INT, MPI_ANY_SOURCE, 0, com, &st);
		}
		MPI_Bcast(&mj_i, 1, MPI_INT, stat, com);
		MPI_Bcast(&mj_j, 1, MPI_INT, stat, com);
		*/
	}
	else{
		int k = n/m;
		//int l = n-k*m;
		double s;
		int i = step, j = step, rrows = (get_rows_l(n ,m ,p , kk)%m!=0);
		bool stat = true;
		mj = 0;
		st = start(step, kk, p, k);
		if(st == step){
			for(int i1= g2l(1, p, st); i1 < rows-rrows; i1++){
				for(j = step; j < k; j++){
					i = g2l(1, p, st);
					get_block(a, block1, n, m, i, j);
					s = norm1(block1, m);
					if(!(s < global_mx*eps*2) && !inversed_matr(m, 1, block1, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
						stat = false;
						s = norm1(block1, m);
						if(mj > 0){
							if(mj >= s){
								mj = s;
								mj_i = l2g(1, p, kk, i);
								mj_j = j;
							}
						}
						else{
							mj = s;
							mj_i = l2g(1, p, kk, i);
							mj_j = j;
						}
					}
				}
			}
		}
		if(!stat){
			mj = 1/mj;
		}
		maj.mj = mj; maj.mj_i = mj_i; maj.mj_j = mj_j;
		MPI_Allreduce(&maj, &ans, 1, ctype, my_MAX, com);
		MPI_Op_free(&my_MAX);
		swap(ans, maj);
		mj_i = maj.mj_i; mj_j = maj.mj_j; mj = maj.mj;
		if(mj < eps*eps*global_mx){
			return false;
		}
	}
	perm1[step] = mj_i;
	perm2[step] = mj_j;
	if(mj_i != step){
		matr_swap_rows_par(com, n, m, kk, p, a, temp1, block1, mj_i, step);
	}
	if(mj_j != step){
		matr_swap_columns_par(n, m, kk, p, a, mj_j, step);
	}
	return true;
}
bool inversed_matr_par(MPI_Comm com, int n, int m, int kk, int p, double *a, double* temp1, double *temp2, double* block1, double* block2, double* block3, int* perm1, int* perm2, int* pperm1, int* pperm2){
	if((m < 1) || (m > n)){
		return true;
	}
	int k = n/m;
	int l = n-k*m;
	if(perm1 != nullptr){
		for(int i = 0; i < k; ++i){
			perm1[i] = i;
			perm2[i] = i;
		}
	}
	double* corner;
	int i, j, rcounter = 0;
	MPI_Status st;
	int rows = get_rows(n, m, p, kk);
	int rrows = (get_rows_l(n, m, p, kk)%m != 0);
	double s = 0;
	double global_mx = norm1_par(a, n, m, kk, p, temp1, temp2, com);
	if(m == 1){
		for(int q = 0; q < n; ++q){
			if(!major_element_par(com, n, 1, kk, p, q, global_mx, a, temp1, perm1, perm2, nullptr, nullptr, nullptr)){
				return true;
			}
			if(kk == q%p){
				for(i = 0; i < n; ++i){
					temp1[i] = a[rcounter*n+i];
				}
				for(int i = 0; i < kk; ++i){
					MPI_Send(temp1, n, MPI_DOUBLE, i, 0, com);
				}
				for(int i = kk+1; i < p; ++i){
					MPI_Send(temp1, n, MPI_DOUBLE, i, 0, com);
				}
			}
			else{
				MPI_Recv(temp1, n, MPI_DOUBLE, q%p, 0, com, &st);
			}
			//MPI_Bcast(temp1, n, MPI_DOUBLE, q%p, com);
			s = 1/temp1[q];
			if(kk == q%p){
				a[rcounter*n+q]=s;
				for(i = 0; i < rcounter; i++){
					a[i*n+q]*=s;
				}
				for(i = rcounter+1; i < rows; i++){
					a[i*n+q]*=s;
				}
				for(i = 0; i < rcounter; i++){
					for(j = 0; j < q; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
					for(j = q+1; j < n; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
				}
				for(i = rcounter+1; i < rows; i++){
					for(j = 0; j < q; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
					for(j = q+1; j < n; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
				}
				for(i = 0; i < q; ++i){
					a[rcounter*n+i]*=-s;
				}
				for(int i = q+1; i < n; ++i){
					a[rcounter*n+i]*=-s;
				}
				rcounter++;
			}
			else{
				for(i = 0; i < rows; i++){
					a[i*n+q]*=s;
				}
				for(i = 0; i < rows; i++){
					for(j = 0; j < q; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
					for(j = q+1; j < n; ++j){
						a[i*n+j]-=a[i*n+q]*temp1[j];
					}
				}
			}
		}
		apply_perm_par(com, n, 1, kk, p, a, temp1, nullptr, perm1, perm2);
		return false;
	}
    if(l == 0){
		for(int q= 0;  q < k; ++q){
			if(!major_element_par(com, n, m, kk, p, q, global_mx, a, temp1, perm1, perm2, pperm1, pperm2, block2)){
					return true;
			}
			if(kk == q%p){
				for(i = 0; i < k; ++i){
					get_block(a, block2, n, m, rcounter, i);
					set_block(temp1, block2, n, m, 0, i);
					//temp1[i] = a[rcounter*n*m+i];
				}
				for(int i = 0; i < kk; ++i){
					MPI_Send(temp1, n*m, MPI_DOUBLE, i, 0, com);
				}
				for(int i = kk+1; i < p; ++i){
					MPI_Send(temp1, n*m, MPI_DOUBLE, i, 0, com);
				}
			}
			else{
				MPI_Recv(temp1, n*m, MPI_DOUBLE, q%p, 0, com, &st);
			}
			//MPI_Bcast(temp1, n*m, MPI_DOUBLE, q%p, com);
			get_block(temp1, block2, n, m, 0, q);
			if(inversed_matr(m, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
				return true;
			}
			if(kk == q%p){
				set_block(a, block2, n, m, rcounter, q);
				for(i = 0; i < rcounter; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				for(i = rcounter+1; i < rows; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				for(i = 0; i < rcounter; i++){
					for(int j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
				}
				for(i = rcounter+1; i < rows; i++){
					for(j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
				}
				get_block(a, block2, n, m, rcounter, q);
				mult_by_minus_one(block2, m, m);
				for(i = 0; i < q; i++){
					get_block(a, block1, n, m, rcounter, i);
					matr_mult(block3, block2, block1, m, m, m);
					set_block(a, block3, n, m, rcounter, i);
					//a[i*n+q]*=-a[q*n+q];
				}
				for(i = q+1; i < k; i++){
					get_block(a, block1, n, m, rcounter, i);
					matr_mult(block3, block2, block1, m, m, m);
					set_block(a, block3, n, m, rcounter, i);
					//a[i*n+q]*=-a[q*n+q];
				}
				++rcounter;
			}
			else{
				for(i = 0; i < rows; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				for(i = 0; i < rows; i++){
					for(int j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
				}
			}
		}
		apply_perm_par(com, n, m, kk, p, a, temp1, block1, perm1, perm2);
		return false;
	}
	else{
		for(int q= 0;  q < k; ++q){
			if(!major_element_par(com, n, m, kk, p, q, global_mx, a, temp1, perm1, perm2, pperm1, pperm2, block2)){
				return true;
			}
			if(kk == q%p){
				for(i = 0; i <= k; ++i){
					get_block(a, block2, n, m, rcounter, i);
					set_block(temp1, block2, n, m, 0, i);
					//temp1[i] = a[rcounter*n*m+i];
				}
				for(int i = 0; i < kk; ++i){
					MPI_Send(temp1, n*m, MPI_DOUBLE, i, 0, com);
				}
				for(int i = kk+1; i < p; ++i){
					MPI_Send(temp1, n*m, MPI_DOUBLE, i, 0, com);
				}
			}
			else{
				MPI_Recv(temp1, n*m, MPI_DOUBLE, q%p, 0, com, &st);
			}
			//MPI_Bcast(temp1, n*m, MPI_DOUBLE, q%p, com);
			get_block(temp1, block2, n, m, 0, q);
			if(inversed_matr(m, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
				return true;
			}
			if(kk == q%p){
                set_block(a, block2, n, m, rcounter, q);
				for(i = 0; i < rcounter; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				for(i = rcounter+1; i < rows-rrows; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				if(rrows){
					get_block(a, block1, n, m, rows-1, q);
					matr_mult(block3, block1, block2, l, m, m);
					set_block(a, block3, n, m, rows-1, q);
				}
				for(i = 0; i < rcounter; i++){
					for(int j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					get_block(a, block1, n, m, i, q);
					get_block(temp1, block2, n, m, 0, k);
					matr_mult(block3, block1, block2, m, m, l);
					corner = get_corner(a, n, m, i, k);
					matr_sub(corner, block3, n, m, l);
				}
				for(i = rcounter+1; i < rows-rrows; i++){
					for(j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					get_block(a, block1, n, m, i, q);
					get_block(temp1, block2, n, m, 0, k);
					matr_mult(block3, block1, block2, m, m, l);
					corner = get_corner(a, n, m, i, k);
					matr_sub(corner, block3, n, m, l);
				}
				if(rrows){
					for(j = 0; j < q; ++j){
                        get_block(a, block1, n, m, rows-1, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, l, m, m);
						corner = get_corner(a, n, m, rows-1, j);
						matr_sub(corner, block3, n, l, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
                        get_block(a, block1, n, m, rows-1, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, l, m, m);
						corner = get_corner(a, n, m, rows-1, j);
						matr_sub(corner, block3, n, l, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
                    get_block(a, block1, n, m, rows-1, q);
					get_block(temp1, block2, n, m, 0, k);
					matr_mult(block3, block1, block2, l, m, l);
					corner = get_corner(a, n, m, rows-1, k);
					matr_sub(corner, block3, n, l, l);
				}
				get_block(a, block2, n, m, rcounter, q);
                for(i = 0; i < q; i++){
					get_block(a, block1, n, m, rcounter, i);
					mult_by_minus_one(block1, m, m);
					matr_mult(block3, block2, block1, m, m, m);
					set_block(a, block3, n, m, rcounter, i);
					//a[i*n+q]*=-a[q*n+q];
				}
                for(i = q+1; i < k; i++){
					get_block(a, block1, n, m, rcounter, i);
					mult_by_minus_one(block1, m, m);
					matr_mult(block3, block2, block1, m, m, m);
					set_block(a, block3, n, m, rcounter, i);
					//a[i*n+q]*=-a[q*n+q];
				}
				get_block(a, block1, n, m, rcounter, k);
				mult_by_minus_one(block1, m, l);
				matr_mult(block3, block2, block1, m, m, l);
				set_block(a, block3, n, m, rcounter, k);
				++rcounter;
			}
			
			
			
			
			
			
			else{
				for(i = 0; i < rows-rrows; i++){
					get_block(a, block1, n, m, i, q);
					matr_mult(block3, block1, block2, m, m, m);
					set_block(a, block3, n, m, i, q);
					//a[q*n+i]*=s;
				}
				if(rrows){
					get_block(a, block1, n, m, rows-1, q);
					matr_mult(block3, block1, block2, l, m, m);
					set_block(a, block3, n, m, rows-1, q);
				}
				for(i = 0; i < rows-rrows; i++){
					for(int j = 0; j < q; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
						get_block(a, block1, n, m, i, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, m, m, m);
						corner = get_corner(a, n, m, i, j);
						matr_sub(corner, block3, n, m, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					get_block(a, block1, n, m, i, q);
					get_block(temp1, block2, n, m, 0, k);
					matr_mult(block3, block1, block2, m, m, l);
					corner = get_corner(a, n, m, i, k);
					matr_sub(corner, block3, n, m, l);
				}
				if(rrows){
					for(j = 0; j < q; ++j){
                        get_block(a, block1, n, m, rows-1, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, l, m, m);
						corner = get_corner(a, n, m, rows-1, j);
						matr_sub(corner, block3, n, l, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
					for(j = q+1; j < k; ++j){
                        get_block(a, block1, n, m, rows-1, q);
						get_block(temp1, block2, n, m, 0, j);
						matr_mult(block3, block1, block2, l, m, m);
						corner = get_corner(a, n, m, rows-1, j);
						matr_sub(corner, block3, n, l, m);
						//a[i*n+j]-=a[i*n+q]*a[q*n+j];
					}
                    get_block(a, block1, n, m, rows-1, q);
					get_block(temp1, block2, n, m, 0, k);
					matr_mult(block3, block1, block2, l, m, l);
					corner = get_corner(a, n, m, rows-1, k);
					matr_sub(corner, block3, n, l, l);
                }
			}
        }
		//       q = k
		int q = k;
		if(kk == q%p){
			for(i = 0; i <= k; ++i){
				get_block(a, block2, n, m, rcounter, i);
				set_block_l(temp1, block2, n, m, 0, i);
				//temp1[i] = a[rcounter*n*m+i];
			}
		}
		MPI_Bcast(temp1, n*l, MPI_DOUBLE, q%p, com);
		get_block_l(temp1, block2, n, m, 0, q);
		if((l == 1) && (norm1(block2, l) < eps*global_mx)){
				return true;
		}
		if(inversed_matr(l, 1, block2, nullptr, nullptr, nullptr, pperm1, pperm2, nullptr, nullptr)){
			return true;
		}
		if(kk == q%p){
			set_block(a, block2, n, m, rcounter, q);
		}
		for(i = 0; i < rows-rrows; i++){
			get_block(a, block1, n, m, i, q);
			matr_mult(block3, block1, block2, m, l, l);
			set_block(a, block3, n, m, i, q);
			//a[q*n+i]*=s;
		}
		for(i = 0; i < rows-rrows; i++){
			for(int j = 0; j < k; ++j){
				get_block(a, block1, n, m, i, q);
				get_block_l(temp1, block2, n, m, 0, j);
				matr_mult(block3, block1, block2, m, l, m);
				corner = get_corner(a, n, m, i, j);
				matr_sub(corner, block3, n, m, m);
				//a[i*n+j]-=a[i*n+q]*a[q*n+j];
			}
		}
		if(kk == q%p){
			get_block_l(a, block2, n, m, rcounter, q);
			for(i = 0; i < k; i++){
				get_block_l(a, block1, n, m, rcounter, i);
				mult_by_minus_one(block1, l, m);
				matr_mult(block3, block2, block1, l, l, m);
				set_block_l(a, block3, n, m, rcounter, i);
				//a[i*n+q]*=-a[q*n+q];
			}
		}
		apply_perm_par(com, n, m, kk, p, a, temp1, block1, perm1, perm2);
		return false;
	}
}











