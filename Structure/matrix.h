#ifndef __MATRIX_H
#define __MATRIX_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

#undef DOUBLE
#undef INTEGER

typedef long INTEGER;
typedef double DOUBLE;
typedef DOUBLE* VECTOR;
typedef DOUBLE* MATRAN;
typedef DOUBLE** MATRIX;

/*macros*/
#define CLEAR(A,N,i) {\
	clear_raw(i,A,N);\
    clear_col(i,A,N);\
};

#define BANDU(x,y) ((NB - 1) + (y - x) + x * NB)
#define BANDI(x,y) ((x >= y) ? BANDU(x,y) : BANDU(y,x))

/*local*/
void vec_alloc(VECTOR& A,INTEGER n);
void mat_alloc(MATRIX& A,INTEGER n,INTEGER m);
void sec_alloc(MATRIX& A,const MATRIX K,INTEGER x,INTEGER y,INTEGER n,INTEGER m);
void vec_free(MATRIX& A,int is_section = 0);
void vec_free(VECTOR& A);
void multiply(const MATRIX A,const MATRIX B,MATRIX C,INTEGER n,INTEGER m,INTEGER l);
void multiply(const MATRIX A,const VECTOR B,VECTOR C,INTEGER n,INTEGER m);
void multiply(const VECTOR A,const VECTOR B,DOUBLE& c,INTEGER n); 
void multiply(const MATRIX A,const DOUBLE f,MATRIX C,INTEGER n,INTEGER m);
void multiply(const VECTOR A,const DOUBLE f,VECTOR C,INTEGER n);
void multiply(const MATRIX A,const DOUBLE f,INTEGER n,INTEGER m);
void multiply(const VECTOR A,const DOUBLE f,INTEGER n);
void multiply(const VECTOR A,const VECTOR B,MATRIX C,INTEGER n,INTEGER m);
void add(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m);
void add(const VECTOR A,const VECTOR B,INTEGER n);
void sub(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m);
void sub(const VECTOR A,const VECTOR B,INTEGER n);
void equ(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m);
void equ(const VECTOR A,const VECTOR B,INTEGER n);
void equsec(const MATRIX A,const MATRIX B,INTEGER n,INTEGER m,INTEGER x1,INTEGER y1,INTEGER x2,INTEGER y2);
void clear(const MATRIX A,INTEGER n,INTEGER m);
void clear(const VECTOR A,INTEGER n);
void clear_raw(INTEGER k,const MATRIX A,INTEGER n);
void clear_col(INTEGER k,const MATRIX A,INTEGER n);
void transpose(const MATRIX A,MATRIX B,INTEGER n,INTEGER m);
void transpose(MATRIX A,INTEGER n,INTEGER m);
void switch_rows(MATRIX A,INTEGER m,INTEGER i,INTEGER j);
void switch_cols(MATRIX A,INTEGER n,INTEGER i,INTEGER j); 
void condense(MATRIX K,INTEGER N,INTEGER x);
void condense(MATRIX K,VECTOR F,INTEGER N,INTEGER x);

/*imported*/
int invert(const MATRIX A,MATRIX B,INTEGER N,INTEGER LDA = 0,INTEGER LDB = 0);
int invertP(MATRIX A,INTEGER N,INTEGER LDA);
int solve(const MATRIX A,INTEGER N,const VECTOR B,VECTOR X);
int solve_from_factor_banded(const MATRAN A,INTEGER N,INTEGER NB,const VECTOR B,VECTOR X);
int solve_from_factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT,const VECTOR B,VECTOR X);
int factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT);
int factor_banded(const MATRAN A,INTEGER N,INTEGER NB);
int multiply_banded(const MATRAN A,const VECTOR B,VECTOR C,INTEGER N,INTEGER NB);
bool condense(MATRIX K,MATRIX k,MATRIX T,INTEGER N,INTEGER DUN);
int convert_symm_full(const MATRAN A,MATRAN B,INTEGER N,INTEGER NB,INTEGER LDB);
int find_eigen(MATRAN K,VECTOR M,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR d,VECTOR v);
int find_buckling(MATRAN K,MATRAN G,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR d,VECTOR v);
int find_eigen(MATRIX K,INTEGER N,VECTOR eigval,VECTOR eigvec);


#endif

