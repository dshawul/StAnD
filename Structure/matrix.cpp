#include "matrix.h"

/*
matrix alloc/free
*/
void mat_alloc(MATRIX& A,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	register VECTOR M;

	A = new DOUBLE* [n];
	M = new DOUBLE[n * m];
	for(i = 0;i < n;i++) {
		A[i] = M;
	    M += m;
	}
	for(i = 0;i < n;i++)
		for(j = 0;j < m;j++)
			A[i][j] = 0;
}
void vec_alloc(VECTOR& A,INTEGER n) {
	register INTEGER i;
	A = new DOUBLE[n];
	for(i = 0;i < n;i++)
		A[i] = 0;
}
void sec_alloc(MATRIX& A,const MATRIX K,INTEGER x,INTEGER y,INTEGER n,INTEGER m) {
	register INTEGER i;
	A = new DOUBLE* [n];
	for(i = 0;i < n;i++)
		A[i] = K[i + y] + x;
}
void vec_free(MATRIX& A,int is_section) {
	if(!is_section) {
		delete[] A[0];
	}
	delete A;
	A = 0;
}
void vec_free(VECTOR& A) {
	delete[] A;
	A = 0;
}
/*
matrix operations
*/
void multiply(const MATRIX A,const MATRIX B,MATRIX C,INTEGER n,INTEGER m,INTEGER l) {
	register INTEGER i,j,k;
	DOUBLE sum;
   	for(i = 0; i < n;i++) {
		for(j = 0; j < l;j++) {
			sum = 0;
			for(k = 0;k < m;k++)
				sum = sum + A[i][k] * B[k][j];
			C[i][j] = sum;
		}
	}
}
void multiply(const MATRIX A,const VECTOR B,VECTOR C,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	DOUBLE sum;
   	for(i = 0; i < n;i++) {
		sum = 0;
		for(j = 0; j < m;j++)
			sum = sum + A[i][j] * B[j];
		C[i] = sum;
	}
}
void multiply(const VECTOR A,const VECTOR B,DOUBLE& c,INTEGER n) {
	c = 0;
   	for(int i = 0; i < n;i++) {
		c += A[i] * B[i];
	}
}
void multiply(const VECTOR A,const VECTOR B,MATRIX C,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < m;j++) {
			C[i][j] = A[i] * B[j];
		}
	}
}
void multiply(const MATRIX A,const DOUBLE f,MATRIX C,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < m;j++) {
			C[i][j] = f * A[i][j];
		}
	}
}
void multiply(const VECTOR A,const DOUBLE f,VECTOR C,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++) {
		C[i] = f * A[i];
	}
}
void multiply(const MATRIX A,const DOUBLE f,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < m;j++) {
			A[i][j] *= f;
		}
	}
}
void multiply(const VECTOR A,const DOUBLE f,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++) {
		A[i] *= f;
	}
}
void add(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			A[i][j] += B[i][j];
}
void add(const VECTOR A,const VECTOR B,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++)
		A[i] += B[i];
}
void sub(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			A[i][j] -= B[i][j];
}
void sub(const VECTOR A,const VECTOR B,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++)
		A[i] -= B[i];
}
void equ(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			A[i][j] = B[i][j];
}
void equ(const VECTOR A,const VECTOR B,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++)
		A[i] = B[i];
}
void equsec(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m,INTEGER x1,INTEGER y1,INTEGER x2,INTEGER y2) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			A[i + y1][j + x1] = B[i + y2][j + x2];
}
void clear(const MATRIX A,INTEGER n,INTEGER m) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			A[i][j] = 0;
}
void clear(const VECTOR A,INTEGER n) {
	register INTEGER i;
   	for(i = 0; i < n;i++)
		A[i] = 0;
}
void clear_raw(INTEGER k,const MATRIX A,INTEGER n) {
	for(INTEGER i = 0;i < n;i++)
		A[k][i] = 0;
}
void clear_col(INTEGER k,const MATRIX A,INTEGER n) {
	for(INTEGER i = 0;i < n;i++)
		A[i][k] = 0;
}
void transpose(const MATRIX A,MATRIX B,INTEGER n,INTEGER m) {
	register INTEGER i,j;
   	for(i = 0; i < n;i++)
		for(j = 0; j < m;j++)
			B[j][i] = A[i][j];
}
void transpose(MATRIX A,INTEGER n,INTEGER m) {
	register INTEGER i,j;
	DOUBLE temp;
   	for(i = 0; i < n;i++)
		for(j = i + 1; j < m;j++) {
			temp = A[j][i];
			A[j][i] = A[i][j];
            A[i][j] = temp;
		}
}
void switch_rows(MATRIX A,INTEGER m,INTEGER i,INTEGER j) {
	register int k;
	register DOUBLE temp;
    for(k = 0; k < m;k++) {
		temp = A[i][k];
		A[i][k] = A[j][k];
		A[j][k] = temp;
	}
}
void switch_cols(MATRIX A,INTEGER n,INTEGER i,INTEGER j) {
	register int k;
	register DOUBLE temp;
    for(k = 0; k < n;k++) {
		temp = A[k][i];
		A[k][i] = A[k][j];
		A[k][j] = temp;
	}
}
/*
static condensation
*/
void condense(MATRIX K,INTEGER N,INTEGER x) {
	register INTEGER i,j;
	for(i = 0;i < N;i++) {
		if(i == x) continue;
		for(j = 0;j < N;j++) {
			if(j == x) continue;
			K[i][j] -= (K[i][x] * K[x][j]) / K[x][x];
		}
	}
	CLEAR(K,N,x);
}
void condense(MATRIX K,VECTOR F,INTEGER N,INTEGER x) {
	register INTEGER i,j;
	for(i = 0;i < N;i++) {
        if(i == x) continue;
		for(j = 0;j < N;j++) {
			if(j == x) continue;
			K[i][j] -= (K[i][x] * K[x][j]) / K[x][x];
		}
		F[i] -= (F[x] * K[i][x]) / K[x][x];
	}
	CLEAR(K,N,x);
	F[x] = 0;
}
bool condense(MATRIX K,MATRIX k,MATRIX T,INTEGER N,INTEGER DUN) {

	register int i,code;
	bool ret = true;
	MATRIX Koo,Kot,Kto;
	
	sec_alloc(Koo,K,DUN,DUN,N - DUN,N - DUN);
	sec_alloc(Kot,K,0,DUN,N - DUN,DUN);
	sec_alloc(Kto,K,DUN,0,DUN,N - DUN);

	code = invertP(Koo,N - DUN,N);
	if(code == 1) {
		ret = false;
		goto END;
	}
	multiply(Koo,Kot,&T[DUN],N - DUN,N - DUN,DUN);
	multiply(&T[DUN],-1,N - DUN,DUN);
    for(i = 0;i < DUN;i++) T[i][i] = 1;

	multiply(Kto,&T[DUN],k,DUN,N - DUN,DUN);
	add(k,K,DUN,DUN);
END:
	vec_free(Koo,1);
	vec_free(Kot,1);
	vec_free(Kto,1);
	
	return ret;
}