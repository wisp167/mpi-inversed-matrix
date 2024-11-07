#include "init_matr.h"
bool init_by_form(double *a, int n, int s){
	if(s < 0 || s > 4){
		return false;
	}
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			a[i*n+j] = init_formula(i+1, j+1, n, s);
		}
	}
	return true;
}

bool init_matr(double* a, int n, int s, ifstream &fin){
	if(s){
		return init_by_form(a, n, s);
	}
	for(int i = 0; i < n*n; ++i){
		if(!(fin >> a[i])){
			return false;
		}
	}
	return true;
}
void init_E(double *a, int n){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			a[i*n+j] = ((i == j) ? 1 : 0);
		}
	}
}
