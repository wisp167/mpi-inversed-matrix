#include "checking_results.h"
#include "header.h"
#include "print_matr.h"
double norm1(double* a, int n){
	double norm;
	double mx = 0;
	for(int j = 0; j < n; ++j){
		norm = 0;
		for(int i = 0; i < n; ++i){
			norm+=fabs(a[i*n+j]);
		}
		mx = max(norm, mx);
	}
	return mx;
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
double norm1_par(double* a, int n, int m, int kk, int p, double *temp1, double *temp2, MPI_Comm com){
	double norm;
    double mx = 0;
    int rows = get_rows_l(n ,m, p, kk);
	for(int j = 0; j < n; ++j){
		norm = 0;
		for(int i = 0; i < rows; ++i){
			norm+=fabs(a[i*n+j]);
		}
		temp1[j] = norm;
	}
	MPI_Allreduce(temp1, temp2, n, MPI_DOUBLE, MPI_SUM, com);
	for(int i = 0; i < n; ++i){
		mx = max(mx, temp2[i]);
	}
	return mx;
}
/*									// useless
void matrix_mult_vector(double* a, double* b, double* c, int n, int m, int p){ //rows x vec
	int j, i2;
	double s;
	int k = n/m;
	static pthread_mutex_t mu= PTHREAD_MUTEX_INITIALIZER;
	for(int i = k*m; i < n; i+=p*m){
		int h = (i+m < n ? m : n-i);
		for(i2 = i; i2 < i+h; i2++){
			s = 0;
			for(j = 0; j < n; ++j){
				s+=a[i2*n+j]*b[j];
			}
			pthread_mutex_lock(&mu);
			c[i2] += s;
			pthread_mutex_unlock(&mu);
		}
	}
}
void matrix_mult_vector1(double* a, double*b, double*c, int n, int m, int p){ //columns x vec
	int j, i2;
	double s, bb;
	int k = n/m;
	static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
	for(int i = k*m; i < n; i+=p*m){
		int h = (i+m < n ? m : n-i);
		for(i2 = i; i2 < i+h; i2++){
			s = 0;
			bb = b[i2];
			for(j = 0; j < n; ++j){
				s+=a[i2*n+j]*bb;
			}
			pthread_mutex_lock(&mu);
			c[j] += s;
			pthread_mutex_unlock(&mu);
		}
	}
}
*/
inline void upd(double* block, int m1, int m2, double* temp2, int i){
	double s=  0;
	for(int j= 0; j < m2; ++j){
		s = 0;
		for(int z = 0; z < m1; ++z){
			s+=fabs(block[z*m2+j]);
		}
		temp2[i+j] += s;
	}
	return;
}
inline void upd1(double* block, int m1, int m2, double* temp2, int i){
	double s=  0;
	for(int z= 0; z < m1; ++z){
		s = 0;
		for(int j = 0; j < m2; ++j){
			s+=fabs(block[z*m2+j]);
		}
		temp2[i+z] += s;
	}
	return;
}
inline double* get_corner(double* a, int n, int m, int i, int j){
	return a + i*n*m+j*m;
}
inline void get_block(double* a, double* block, int n, int m, int i, int j){
	double* corner = get_corner(a, n, m, i, j);
	int k = n/m;
	int l = n-k*m;
	int h,w;
	w = ((j < k) ? m : l);
	h = ((i < k) ? m : l);
	for(int q1 = 0; q1 < h; ++q1){
		for(int q2 = 0; q2 < w; ++q2){
			block[q1*w+q2] = corner[q1*n+q2];
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
double find_r_par(double* a, double* res, int kk, int p, int n, int m, double* block1, double* block2, double* block3, double* block4, double* temp1, double* temp2){
	MPI_Comm com = MPI_COMM_WORLD;
	int k = n/m, l = n-k*m;
	int sz = k + (l != 0);
	int rows = get_rows(n, m, p, kk);
	int rrows = (get_rows_l(n, m, p, kk)%m != 0);
	if(n > 11000){
		return 0;
	}
	for(int i = 0; i < n; i++){
		temp2[i] = 0;
		temp1[i] = 0;
	}
	int h, lf, r;
	int curr = kk;
	for(int i = 0; i < m*m; ++i){
		block4[i] = 0;
	}
	MPI_Status st;
	//MPI_Request rq;
	bool stat = false;
	for(int ii = 0; ii < p; ++ii){
		int k1 = (kk+ii)%p;
		int rows1 = get_rows(n, m, p, k1);
		int rrows1 = (get_rows_l(n, m, p, k1) % m != 0);
		for(int t = 0; t < rows1; ++t){
			curr = k1+t*p;
			stat = (rrows1 && (t == rows1-1));
			for(int q = 0; q < rows; q++){
				r = ((curr != k) ? m : l);
				lf = ((rrows && (q == rows-1)) ? l : m);
				if(lf == m){
					for(int j = 0; j < sz; ++j){
						get_block(a, block1, n, m, q, j);
						if(!stat){
							get_block_tr(res, block2, n, m, t, j);
						}
						else{
							get_block_tr_l(res, block2, n, m, t, j);
						}
						h = ((j != k) ? m : l);
						matr_mult(block3, block1, block2, lf, h, r);
						for(int z = 0; z < lf*r; ++z){
							block4[z] += block3[z];
						}
					}
				}
				else{
					for(int j = 0; j < sz; ++j){
						get_block_l(a, block1, n, m, q, j);
						if(!stat){
							get_block_tr(res, block2, n, m, t, j);
						}
						else{
							get_block_tr_l(res, block2, n, m, t, j);
						}
						h = ((j != k) ? m : l);
						matr_mult(block3, block1, block2, lf, h, r);
						for(int z = 0; z < lf*r; ++z){
							block4[z] += block3[z];
						}
					}
				}
				
				if(curr == (kk + q*p)){
					for(int w = 0; w < r;++w){
						block4[w*r+w]--;
					}
				}
				for(int w = 0; w < lf*r; ++w){
					block3[w] = block4[w];
					block4[w]  = 0;
				}
				upd(block3, lf,r, temp2, curr*m);
				for(int z = curr*m; z < curr*m+r; ++z){
					temp1[z]+=fabs(temp2[z]);
					temp2[z] = 0;
				}
			}
		}
		
		MPI_Sendrecv_replace(res, n*(n/p + m), MPI_DOUBLE, (kk+p-1)%p, 0, (kk+1)%p, 0, com, &st);
	}
	double ans = 0;
	MPI_Allreduce(temp1, temp2, n, MPI_DOUBLE, MPI_SUM, com);
	for(int i = 0; i < n; ++i){
		ans = max(ans, temp2[i]);
	}
	return ans;
}
/*
double find_r2(double* res, double* a, int kk, int p, int n, int m, double* block1, double* block2, double* block3, double* block4, double* temp1, double* temp2){
	int k = n/m, l = n-k*m;
	int sz = k+(l!=0);
	if(n > 11000){
		return 0;
	}
	for(int i = kk*m; i < n; i+=p*m){
		int h = (i+m < n ? m : n-i);
		for(int j = i; j < i+h; ++j){
			temp2[j] = 0;
			temp1[j] = 0;
		}
	}
	int h, lf, r;
	reduce_sum(p);
	int curr = kk;
	for(int i = 0; i < m*m; ++i){
		block4[i] = 0;
	}
	for(int i = 0; i < sz; ++i){
		curr = (kk+i)%sz;
		for(int q = kk; q < sz; q+=p){
			r = ((curr != k) ? m : l);
			lf = ((q != k) ? m : l);
			for(int j = 0; j < sz; ++j){
				get_block_tr(res, block1, n, m, j, q);
				get_block(a, block2, n, m, j, curr);
				h = ((j != k) ? m : l);
				matr_mult(block3, block1, block2, lf, h, r);
				for(int i = 0; i < lf*r; ++i){
					block4[i] += block3[i];
				}
			}
			if(curr==q){
				for(int w = 0; w < r;++w){
					block4[w*r+w]--;
				}
			}
			for(int w = 0; w < lf*r; ++w){
				block3[w] = block4[w];
				block4[w]  = 0;
			}
			upd(block3, lf,r, temp2, curr*m);
			for(int z = curr*m; z <curr*m+r; ++z){
				temp1[z]+=fabs(temp2[z]);
				temp2[z] = 0;
			}
		}
		reduce_sum(p);
	}
	double ans = 0;
	for(int i = 0; i < n; ++i){
		ans = max(ans, temp1[i]);
	}
	reduce_sum(p);
	return ans;
}
void fast_matr_mult_test(int n, int m, double* a, double* b, double* pc){
	int k = n/m, bl;
	int l = n-k*m, v, h, v3, h3, i, j, s, r, t, q, ah;
	double *pa, *pb;
	double s00, s01, s02, s10, s20, s11, s22, s12, s21;
	bl = ((l!=0) ? k+1 : k);
	for(i = 0; i < bl; ++i){
		for(j = 0; j < bl; ++j){
			v = (i < k ? m : l); h = (j < k ? m : l);
			for(r = 0; r < v; ++r){
				for(t = 0; t < h; ++t){
					pc[r*h+t] = 0;
				}
			}
			for(s = 0; s < bl; ++s){
				ah = (s < k ? m : l);
				pa = a+i*n*m+s*m;
				pb = b+s*n*m+j*m;
				v3 = v%3; h3 = h%3;
				for(r = 0; r < v3; ++r){
					for(t = 0; t < h3; ++t){
						s = 0;
						for(q= 0; q < h; ++q){
							s+=pa[r*n+q]*pb[q*n+t];
						}
						pc[r*h+t]+=s;
					}
					for(; t < h; t+=3){
						s00 = 0; s01 = 0; s02 = 0;
						for(q = 0; q < bl; ++q){
							s00+=pa[r*n+q]*pb[q*n+t];
							s01+=pa[r*n+q]*pb[q*n+t+1];
							s02+=pa[r*n+q]*pb[q*n+t+2];
						}
						pc[r*h+t]+=s00;
						pc[r*h+t+1]+=s01;
						pc[r*h+t+2]+=s02;
					}
				}
				for(;r<v;r+=3){
					for(t = 0; t<h3; ++t){
						s00 = s10 = s20 = 0;
						for(q = 0; q < h; q++){
							s00+=pa[r*n+q]*pb[q*n+t];
							s10+=pa[(r+1)*n+q]*pb[q*n+t];
							s20+=pa[(r+2)*n+q]*pb[q*n+t];
						}
						pc[r*h+t]+=s00;pc[(r+1)*h+t]+=s10;pc[(r+2)*h+t]+=s20;
					}
					for(;t<h;t+=3){
						s01=s02=s10=s20=s11=s12=s21=s22=s00=0;
						for(q=0;q<ah;q++){
							s00+=pa[r*n+q]*pb[q*n+t];
							s01+=pa[r*n+q]*pb[q*n+t+1];
							s02+=pa[r*n+q]*pb[q*n+t+2];
							s10+=pa[(r+1)*n+q]*pb[q*n+t];
							s11+=pa[(r+1)*n+q]*pb[q*n+t+1];
							s12+=pa[(r+1)*n+q]*pb[q*n+t+2];
							s20+=pa[(r+2)*n+q]*pb[q*n+t];
							s21+=pa[(r+2)*n+q]*pb[q*n+t+1];
							s22+=pa[(r+2)*n+q]*pb[q*n+t+2];
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
		}
	}
}
*/
void matr_mult_norm(int n, int m, double* a, double* b, double* temp1, double* pc){
	int k = n/m, bl, t1, t2, t3;
	int l = n-k*m, v, h, v3, h3, i, j,s, r, t, q, ah;
	double *pa, *pb, ss;
	double s00, s01, s02, s10, s20, s11, s22, s12, s21;
	bl = ((l!=0) ? k+1 : k);
	for(i = 0; i < n; ++i){
		temp1[i] = 0;
	}
	for(i = 0; i < bl; ++i){
		for(j = 0; j < bl; ++j){
			v = (i < k ? m : l); h = (j < k ? m : l);
			for(s = 0; s < v*h; ++s){
				pc[s]=0;
			}
			for(s = 0; s < bl; ++s){
				ah = (s < k ? m : l);
				pa = a+i*n*m+s*m;
				pb = b+s*n*m+j*m;
				v3 = v%3; h3 = h%3;
				for(r = 0; r < v3; ++r){
					for(t = 0; t < h3; ++t){
						ss = 0;
						for(q= 0; q < ah; ++q){
							ss+=pa[r*n+q]*pb[q*n+t];
						}
						pc[r*h+t]+=ss;
					}
					for(; t < h; t+=3){
						s00 = 0; s01 = 0; s02 = 0;
						for(q = 0; q < ah; ++q){
							s00+=pa[r*n+q]*pb[q*n+t];
							s01+=pa[r*n+q]*pb[q*n+t+1];
							s02+=pa[r*n+q]*pb[q*n+t+2];
						}
						pc[r*h+t]+=s00;
						pc[r*h+t+1]+=s01;
						pc[r*h+t+2]+=s02;
					}
				}
				for(;r<v;r+=3){
					for(t = 0; t<h3; ++t){
						s00 = s10 = s20 = 0;
						for(q = 0; q < ah; q++){
							s00+=pa[r*n+q]*pb[q*n+t];
							s10+=pa[(r+1)*n+q]*pb[q*n+t];
							s20+=pa[(r+2)*n+q]*pb[q*n+t];
						}
						pc[r*h+t]+=s00;pc[(r+1)*h+t]+=s10;pc[(r+2)*h+t]+=s20;
					}
					for(;t<h;t+=3){
						s01=s02=s10=s20=s11=s12=s21=s22=s00=0;
						for(q=0;q<ah;q++){
							s00+=pa[r*n+q]*pb[q*n+t];
							s01+=pa[r*n+q]*pb[q*n+t+1];
							s02+=pa[r*n+q]*pb[q*n+t+2];
							s10+=pa[(r+1)*n+q]*pb[q*n+t];
							s11+=pa[(r+1)*n+q]*pb[q*n+t+1];
							s12+=pa[(r+1)*n+q]*pb[q*n+t+2];
							s20+=pa[(r+2)*n+q]*pb[q*n+t];
							s21+=pa[(r+2)*n+q]*pb[q*n+t+1];
							s22+=pa[(r+2)*n+q]*pb[q*n+t+2];
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
			t3 = j*m;
			if(i!=j){
				for(t1 = 0; t1 < v; ++t1){
					for(t2 =0; t2 < h; ++t2){
						temp1[t3+t2]+=fabs(pc[t1*h+t2]);
					}
				}
			}
			else{
				for(t1 = 0; t1 < v; ++t1){
					for(t2 = 0; t2 < t1; ++t2){
						temp1[t3+t2]+=fabs(pc[t1*h+t2]);
					}
					temp1[t3+t1]+=fabs(pc[t1*h+t1]-1);
					for(t2 = t1+1; t2 < h; ++t2){
						temp1[t3+t2]+=fabs(pc[t1*h+t2]);
					}
				}
			}
		}
	}
}

double find_r(double* a, double* b, int n, int m, double* temp1, double* temp2){
	if(n > 11000){
		return 0;
	}
	matr_mult_norm(n, m, a, b, temp1, temp2);
	double mx = 0;
	for(int i = 0; i < n; ++i){
		mx = max(mx, temp1[i]);
	}
	return mx;
}
