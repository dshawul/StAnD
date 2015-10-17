#ifndef __UTILS_H
#define __UTILS_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

/*
print
*/
void open_file(const char* name,bool append,bool log = false);
void close_file(bool log = false);
void print(char* format,...);
void print_matrix(const MATRIX A,int n,int m);
void print_vector(const VECTOR A,int n);
void print_vector(const VECTOR A,int n,int m);
/*
rotation
*/
RPoint GetGlobal(RPoint& p,MATRIX Rt);
void GetPointRotMatrix(RPoint rp,MATRIX R);
void GetLineRotMatrix(RPoint rp,DOUBLE alpha,MATRIX R);
void GetPlaneRotMatrix(RPoint rp,DOUBLE alpha,MATRIX R);
/*
end
*/
#endif