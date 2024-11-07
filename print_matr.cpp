#include "print_matr.h"
void print_matr1(double* a, int n, int m, int r){
	for(int i = 0; i < min(r, n); ++i){
		for(int j = 0; j < min(r, m)-1; ++j){
			printf("%10.3e ", a[i*n+j]);
		}
		printf("%10.3e\n", a[i*n+min(r, m)-1]);
	}
}
