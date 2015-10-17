#ifndef __MATRIX_H
#define __MATRIX_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

#define DLLExport __declspec(dllexport)

#undef DOUBLE
#undef INTEGER

typedef long INTEGER;
typedef double DOUBLE;
typedef DOUBLE* VECTOR;
typedef DOUBLE* MATRAN;
typedef DOUBLE** MATRIX;

DLLExport int invert(const MATRIX A,MATRIX B,INTEGER N,INTEGER LDA = 0,INTEGER LDB = 0);
DLLExport int invertP(MATRIX A,INTEGER N,INTEGER LDA);
DLLExport int solve(const MATRIX A,INTEGER N,const VECTOR B,VECTOR X);
DLLExport int solve_from_factor_banded(const MATRAN A,INTEGER N,INTEGER NB,const VECTOR B,VECTOR X);
DLLExport int solve_from_factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT,const VECTOR B,VECTOR X);
DLLExport int factor_banded_full(const MATRAN A,INTEGER N,INTEGER NB,INTEGER* IPVT);
DLLExport int factor_banded(const MATRAN A,INTEGER N,INTEGER NB);
DLLExport int multiply_banded(const MATRAN A,const VECTOR B,VECTOR C,INTEGER N,INTEGER NB);
DLLExport bool condense(MATRIX K,MATRIX k,MATRIX T,INTEGER N,INTEGER DUN);
DLLExport int convert_symm_full(const MATRAN A,MATRAN B,INTEGER N,INTEGER NB,INTEGER LDB);
DLLExport int find_eigen(MATRAN K,VECTOR M,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR d,VECTOR v);
DLLExport int find_buckling(MATRAN K,MATRAN G,INTEGER N,INTEGER NB,INTEGER nev,INTEGER ncv,DOUBLE sigma,DOUBLE tolerance,VECTOR d,VECTOR v);
DLLExport int find_eigen(MATRIX K,INTEGER N,VECTOR eigval,VECTOR eigvec);


#endif

