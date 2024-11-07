#include "header.h"


int main(int argc, char* argv[]){
	int n, m, p, r, s, k;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &k);
	MPI_Comm com = MPI_COMM_WORLD;
	int* temp = (int*)malloc (sizeof (int) * 4);
	if(argc < 5){
		free(temp);
		if(!k){
			err_output(argv[0]);
		}
		MPI_Finalize();
		return -1;
	}
	string ss;
	for(int i = 0; i < 4; ++i){
		ss = argv[i+1];
		if((sscanf(argv[i+1], "%d", &temp[i]) != 1) || (to_string(temp[i]) != ss)){
			free(temp);
			if(!k){
				err_output(argv[0]);
			}
			MPI_Finalize();
			return -1;
		}
	}
	n = temp[0]; m = temp[1]; r = temp[2]; s = temp[3];
	free(temp);
	string filename;
	if((s < 0) || (s > 4) || ((s == 0) && (argc < 6)) || (n < 1)){
		if(!k){
			err_output(argv[0]);
		}
		MPI_Finalize();
		return -1;
	}
	if(s == 0){
		filename = argv[5];
	}
	double *a = nullptr, *res = nullptr;
	double *temp1 = nullptr, *temp2 = nullptr;
	int *perm1 = nullptr, *perm2 = nullptr, *pperm1 = nullptr, *pperm2 = nullptr;
	double *block1 = nullptr, *block2 = nullptr, *block3 = nullptr, *block4 = nullptr;
	try{
		a = (double*)malloc (sizeof (double) * n*(max(3, n/p) + m));
		res = (double*)malloc(sizeof(double)*n*(max(3, n/p) + m));
		perm1 = (int*) malloc(sizeof(int)*n);
		perm2 = (int*) malloc(sizeof(int)*n);
		temp1 = (double*) malloc(sizeof(double)*n*m);
		temp2 = (double*) malloc(sizeof(double)*n);
		pperm1 = (int*) malloc(sizeof(int)*n);
		pperm2 = (int*) malloc(sizeof(int)*n);
		block1 = (double*) malloc(sizeof(double)*m*m);
		block2 = (double*) malloc(sizeof(double)*m*m);
		block3 = (double*) malloc(sizeof(double)*m*m);
		block4 = (double*) malloc(sizeof(double)*m*m);
	}
	catch(std::bad_alloc&){
		if(!k){
			err_output(argv[0]);
		}
		MPI_Finalize();
		return -3;
	}
	double t1 = 0, t2 = 0, r1 = 0, r2 = 0;
	if(MPI_f(com, n, m, t1, t2, r1, r2, s, r, k, p, filename, a, res, temp1, temp2, block1, block2, block3, block4, perm1, perm2, pperm1, pperm2)){
		free(res); free(a); free(perm1); free(perm2); free(temp1); free(temp2);
		free(pperm1); free(pperm2); free(block1); free(block2); free(block3); free(block4);
		if(!k){
			err_output(argv[0]);
		}
		MPI_Finalize();
		return -5;
	}
	MPI_Barrier(com);
    print_matr_t(a, n, m, p, k, temp1, block1, r1, r2, r, com);
	free(res); free(a); free(perm1); free(perm2); free(temp1); free(temp2);
	free(pperm1); free(pperm2); free(block1); free(block2); free(block3); free(block4);
	if(!k){
		printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 20, r1, r2, t1, t2, s, n, m, p);
	}
	MPI_Finalize();
	return 0;
}
