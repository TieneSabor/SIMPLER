#include "basic_linear.h"

TYPE **dmatrix(int nr, int nc){
    int size = nr*nc;
	TYPE **m = (TYPE **)malloc(nr*sizeof(TYPE *));
	//if(m == NULL) {Serial.print("pointer malloc error \r\n");}
	TYPE *data = malloc(size*sizeof(TYPE));
	//if(m == NULL) {Serial.print("data malloc error \r\n");}
	for(int i=0;i<nr;i++){
		m[i] = data;
		data += nc;
	}
	return m;
}

void freedmat(TYPE** A, int m, int n){
	free(A[0]);
	free(A);
}

TYPE **deyes(int n){
	TYPE **result = dmatrix(n,n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i==j){
				result[i][j] = 1;
			}
			else{
				result[i][j] = 0;
			}
		}
	}
	return result;
}

/* return a m*n zero matrix
*/
TYPE **dzero(int m, int n){
	TYPE **result = dmatrix(m,n);
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			result[i][j] = 0;
		}
	}
	return result;
}

void makeeyes(TYPE** A, int m){
	for(int i=0;i<m;i++){
		for(int j=0;j<m;j++){
			if(i==j){
				A[i][j] = 1;
			}
			else{
				A[i][j] = 0;
			}
		}
	}
}

void makezero(TYPE** A, int m, int n){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			A[i][j] = 0;
		}
	}
}

void printdmat(TYPE **A, int nr, int nc){
	printf("The matrix %d * %d: \r\n", nr, nc);
	for(int i=0; i<nr; i++){
		for(int j=0; j<nc; j++){
			printf("%f, ", A[i][j]);
		}
		printf("\r\n");
	}
}

void matadd(TYPE **A, TYPE **B, TYPE **C, int p, int q){
	for(int i=0;i<p;i++){
		for(int j=0;j<q;j++){
			C[i][j] = A[i][j]+B[i][j];
		}
	}
}

void matsub(TYPE **A, TYPE **B, TYPE **C, int p, int q){
	for(int i=0;i<p;i++){
		for(int j=0;j<q;j++){
			C[i][j] = A[i][j]-B[i][j];
		}
	}
}

void matmult(TYPE **A, TYPE **B, TYPE **C, int p, int q, int r){
	makezero(C,p,r);
	for(int i = 0; i<p; i++){
		for(int k=0; k<q; k++){
			for(int j=0; j<r; j++){
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void matscal(TYPE **A, TYPE **B, TYPE s, int m, int n){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			B[i][j] = s*A[i][j];
		}
	}
}

void transpose(TYPE **A, TYPE** B, int nr, int nc){
    for(int i=0; i<nr; i++){
        for(int j=0; j<nc; j++){
            B[j][i] = A[i][j];
        }
    }
}

void Cholesky_getL(TYPE **A, TYPE** L, int n){
	for(int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			if(i==j){
				TYPE sum=0;
				for(int k=0;k<j;k++){
					sum+=pow(L[j][k],2);
				}
				L[i][j] = sqrt(A[i][j]-sum);
			}
			else{
				TYPE sum=0;
				for(int k=0;k<j;k++){
					sum+=L[i][k]*L[j][k];
				}
				L[i][j] = (A[i][j]-sum)/(L[j][j]+MIN);
			}
		}
	}
}

void Solve_LowTri_Sys(TYPE **L, TYPE **x, TYPE **b, int n, int p){
	//TYPE **x = dmatrix(n,1);
	for(int i=0;i<n;i++){
		TYPE sum = 0;
		for(int j=0;j<i;j++){
			sum += L[i][j]*x[j][p];
		}
		x[i][p] = (b[i][p] - sum)/(L[i][i]+MIN);
	}
}
void Solve_UpTri_Sys(TYPE **U, TYPE **x, TYPE **b, int n, int p){
	for(int i=n-1;i>=0;i--){
		TYPE sum = 0;
		for(int j=n-1;j>i;j--){
			sum += U[i][j]*x[j][p];
		}
		x[i][p] = (b[i][p] - sum)/(U[i][i]+MIN);
	}
}

void Cholesky_Inv(TYPE **A, TYPE **CX, TYPE **Ainv, TYPE **L, TYPE **LT, TYPE **I, int n){
	Cholesky_getL(A, L, n);
	transpose(L, LT, n, n);
	for(int i=0;i<n;i++){
		Solve_LowTri_Sys(L, CX, I, n, i);
	}
	for(int i=0;i<n;i++){
		Solve_UpTri_Sys(LT, Ainv, CX, n, i);
	}
}
