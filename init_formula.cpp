#include "init_formula.h"
double init_formula(int i, int j, int n, int s){
	if(s == 1){
		return n-max(j,i)+1;
	}
	else if(s == 2){
		return max(i, j);
	}
	else if(s == 3){
		return abs(i-j);
	}
	return (double)1/(i+j-1);
}
