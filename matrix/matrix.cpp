#include "matrix.h"

#define IMSLCALL extern "C" void __stdcall
#define IMSLINT extern "C" INTEGER __stdcall

IMSLCALL DLFTRB(const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,const INTEGER&,VECTOR,const INTEGER&,INTEGER*);
IMSLCALL DLFTQS(const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,VECTOR,const INTEGER&);

IMSLCALL DLSLSF(const INTEGER&,VECTOR,const INTEGER&,const VECTOR,VECTOR);
IMSLCALL DLFSQS(const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,const VECTOR,VECTOR);
IMSLCALL DLFSRB(const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,const INTEGER&,INTEGER*,const VECTOR,const INTEGER&,VECTOR);

IMSLCALL DLINRG(const INTEGER&,VECTOR,const INTEGER&,VECTOR,const INTEGER&);
IMSLCALL DLINDS(const INTEGER&,VECTOR,const INTEGER&,VECTOR,const INTEGER&);

IMSLCALL DEVCSF(const INTEGER&,VECTOR,const INTEGER&,VECTOR,VECTOR,const INTEGER&);

IMSLCALL DSBMV(const char&,const INTEGER&,const INTEGER&,const INTEGER&,const DOUBLE&,const VECTOR,const INTEGER&,
			   const VECTOR,const INTEGER&,const DOUBLE&,VECTOR,const INTEGER&);

IMSLCALL DCSBRB(const INTEGER&,const VECTOR,const INTEGER&,const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,const INTEGER&);
IMSLCALL DSBAND(const INTEGER&,const char&,const INTEGER&,const INTEGER*,VECTOR,VECTOR,const INTEGER&,const DOUBLE&,
				const INTEGER&,VECTOR,VECTOR,const INTEGER&,VECTOR,const INTEGER&,const INTEGER&,const char*,const INTEGER&,const char&,const INTEGER&,const INTEGER&,
				const DOUBLE&,VECTOR,const INTEGER&,VECTOR,const INTEGER&,const INTEGER*,VECTOR,VECTOR,
				const INTEGER&,const INTEGER*,const INTEGER&);

IMSLCALL ERSET(const INTEGER&,const INTEGER&,const INTEGER&);
IMSLINT IERCD();

/*
vector 
*/
void vec_alloc(VECTOR& A,INTEGER n) {
	register INTEGER i;
	A = new DOUBLE[n];
	for(i = 0;i < n;i++)
		A[i] = 0;
}

void vec_free(VECTOR& A) {
	delete[] A;
	A = 0;
}

/*
using imsl
*/
int invert(const MATRIX A,MATRIX B,INTEGER N,INTEGER LDA,INTEGER LDB) {
	if(!LDA) LDA = N;
	if(!LDB) LDB = N;
    ERSET(0,0,0);
	DLINRG(N,A[0],LDA,B[0],LDB);
	return IERCD();
}
int invertP(MATRIX A,INTEGER N,INTEGER LDA) {
    ERSET(0,0,0);
	DLINDS(N,A[0],LDA,A[0],LDA);
	return IERCD();
}
int solve(const MATRIX A,INTEGER N,const VECTOR B,VECTOR X) {
	ERSET(0,0,0);
	DLSLSF(N,A[0],N,B,X);
	return IERCD();
}
int factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT) {
    ERSET(0,0,0);
	DLFTRB(N,A,3 * NB - 2,NB - 1,NB - 1,A,3 * NB - 2,IPVT);
	return IERCD();
}
int solve_from_factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT,const VECTOR B,VECTOR X) {
	ERSET(0,0,0);
	DLFSRB(N,A,3 * NB - 2,NB - 1,NB - 1,IPVT,B,1,X);
	return IERCD();
}
int factor_banded(const MATRAN A,INTEGER N,INTEGER NB) {
    ERSET(0,0,0);
	DLFTQS(N,A,NB,NB - 1,A,NB);
	return IERCD();
}
int solve_from_factor_banded(const MATRAN A,INTEGER N,INTEGER NB,const VECTOR B,VECTOR X) {
	ERSET(0,0,0);
	DLFSQS(N,A,NB,NB - 1,B,X);
	return IERCD();
}
int multiply_banded(const MATRAN A,const VECTOR B,VECTOR C,INTEGER N,INTEGER NB) {
	ERSET(0,0,0);
	DSBMV('U',1,N, NB - 1,1.0,A,NB,B,1,0.0,C,1);
	return IERCD();
}
int convert_symm_full(const MATRAN A,MATRAN B,INTEGER N,INTEGER NB,INTEGER LDB) {
	ERSET(0,0,0);
	DCSBRB(N,A,NB,NB - 1,B,LDB,NB - 1,NB - 1);
	return IERCD();
}
int find_eigen(MATRIX A,INTEGER N,VECTOR eigval,VECTOR eigvec) {
	ERSET(0,0,0);
	DEVCSF(N,A[0],N,eigval,eigvec,N);
	return IERCD();
}
/*
find eigen values using arpack library
*/
int __find_eigen(MATRAN k,MATRAN m,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR eigval,VECTOR eigvec) {

	INTEGER lda = 3 * NB - 3 + 1;
	VECTOR rfac,resid,workl,workd,V;
	INTEGER iparam[11],*iwork,lworkl = ncv * ncv + 8 * ncv,info = 0;
    INTEGER* select;
	
	/*call arpack dsband*/
	vec_alloc(V,ncv * N);
	vec_alloc(rfac,lda * N);
	vec_alloc(resid,N);
    vec_alloc(workl,lworkl);
	vec_alloc(workd,3 * N);
	iwork = new INTEGER[N];
	select = new INTEGER[ncv];

	iparam[3 - 1] = 300;
	iparam[7 - 1] = 3;

    DSBAND( 1, 'A',1, select, eigval, eigvec, N, sigma, N, k, m, 
		lda, rfac, NB - 1, NB - 1, "LM",2, 'G',1, nev, tolerance, 
		resid, ncv, V, N, iparam, workd, workl, lworkl, 
		iwork, info);
	
	/*free memory*/
	delete[] select;
	delete[] iwork;
	vec_free(V);
	vec_free(workl);
	vec_free(workd);
	vec_free(resid);
	vec_free(rfac);
	return info;
}

int find_eigen(MATRAN K,VECTOR M,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR eigval,VECTOR eigvec) {

    MATRAN k,m;
	INTEGER lda = 3 * NB - 3 + 1;
	INTEGER idiag,j;
	INTEGER info = 0;

	/*construct k and m*/
	vec_alloc(k,lda * N);
	vec_alloc(m,lda * N);

	convert_symm_full(K,&k[NB - 1],N,NB,lda);

	idiag = 2 * NB - 2;
	for(j = 0;j < N;j++) {
		m[idiag + j * lda] = M[j];
	}
	
	info = __find_eigen(k,m,N,NB,nev,ncv,sigma,tolerance,eigval,eigvec);

	/*free*/
	vec_free(k);
	vec_free(m);

	return info;
}

int find_buckling(MATRAN K,MATRAN G,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR eigval,VECTOR eigvec) {

    MATRAN k,g;
	INTEGER lda = 3 * NB - 3 + 1;
	INTEGER info = 0;

	/*construct k and m*/
	vec_alloc(k,lda * N);
	vec_alloc(g,lda * N);

	convert_symm_full(K,&k[NB - 1],N,NB,lda);
	convert_symm_full(G,&g[NB - 1],N,NB,lda);

	info = __find_eigen(k,g,N,NB,nev,ncv,sigma,tolerance,eigval,eigvec);

	/*free*/
	vec_free(k);
	vec_free(g);

	return info;
}
