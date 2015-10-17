#include <math.h>
#include <stdarg.h>
#include <stdio.h>

#include "matrix.h"
#include "rpoint.h"
#include "utils.h"

/*
print
*/
FILE* logfile = NULL;
FILE* specialfile = NULL;

void open_file(const char* name,bool append,bool log) {
	if(log) logfile = fopen(name,append ? "a" : "w");
	else specialfile = fopen(name,append ? "a" : "w");
}
void close_file(bool log) {
	if(log) {
		fclose(logfile);
	    logfile = NULL;
	} else {
		fclose(specialfile);
	    specialfile = NULL;
	}
}
void print(char* format,...) {

	FILE* myfile = (specialfile ? specialfile : logfile);

	va_list ap;
	va_start(ap, format);
	vfprintf(myfile, format, ap);
	fflush(myfile);
	va_end(ap);
}

void print_matrix(const MATRIX A,int n,int m) {
	register int i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < m;j++)
			print("%15.6f ",A[i][j]);
		print("\n");
	}
    print("\n");
}
void print_vector(const VECTOR A,int n) {
	register int i;
	for(i = 0; i < n;i++) {
		print("%15.6f\n",A[i]);
	}
    print("\n");
}
void print_vector(const VECTOR A,int n,int m) {
#define __SQUARE
#ifndef __SQUARE
	register int i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < m;j++)
			print("%15.6g ",A[i * m + j]);
		print("\n");
	}
    print("\n");
#else
	int NB = m;
	register int i,j;
	for(i = 0; i < n;i++) {
		for(j = 0; j < n;j++)
			print("%15.6g ",A[BANDI(i,j)]);
		print("\n");
	}
    print("\n");
#endif
#undef __SQUARE
}
/*
Transformations
*/
RPoint GetGlobal(RPoint& p,MATRIX Rt) {
	RPoint rp;
	DOUBLE vec[3],vec1[3];
	vec[0] = p.x;
	vec[1] = p.y;
	vec[2] = p.z;
    multiply(Rt,vec,vec1,3,3);
    rp.x = vec1[0];
    rp.y = vec1[1];
	rp.z = vec1[2];
	return rp;
}
static void FormRotMatrix(RPoint vx,RPoint vy,RPoint vz,MATRIX R) {
	vx = vx.unitvector();
	vy = vy.unitvector();
	vz = vz.unitvector();

	R[0][0] = vx.x;
	R[0][1] = vx.y;
	R[0][2] = vx.z;
	R[1][0] = vy.x;
	R[1][1] = vy.y;
	R[1][2] = vy.z;
	R[2][0] = vz.x;
	R[2][1] = vz.y;
	R[2][2] = vz.z;
}
void GetPointRotMatrix(RPoint rp,MATRIX R) {
	MATRIX Rx,Ry,Rz,temp;

	mat_alloc(Rx,3,3);
	mat_alloc(Ry,3,3);
	mat_alloc(Rz,3,3);
	mat_alloc(temp,3,3);

	Rx[0][0] = 1;
	Rx[1][1] = Rx[2][2] = cos(rp.x);
	Rx[1][2] = sin(rp.x);
	Rx[2][1] = -sin(rp.x);

	Ry[1][1] = 1;
	Ry[0][0] = Ry[2][2] = cos(rp.y);
	Ry[2][0] = sin(rp.y);
	Ry[0][2] = -sin(rp.y);

	Rz[2][2] = 1;
	Rz[0][0] = Rz[1][1] = cos(rp.z);
	Rz[0][1] = sin(rp.z);
	Rz[1][0] = -sin(rp.z);

    multiply(Rx,Ry,temp,3,3,3);
	multiply(temp,Rz,R,3,3,3);

	vec_free(Rx);
	vec_free(Ry);
	vec_free(Rz);
	vec_free(temp);
}
void GetLineRotMatrix(RPoint vx,DOUBLE alpha,MATRIX R) {
	RPoint vy,vz;
	//vy
	vz.set(0,0,1);
	vy = cross(vz,vx);
	if(EQUAL(vy.magnitude(),0)) {
		vz.set(-1,0,0);
		vy = cross(vz,vx);
	}
	vy = rotate(vy,vx.unitvector(),alpha);
	
    //vz
	vz = cross(vx,vy);

	FormRotMatrix(vx,vy,vz,R);
}
void GetPlaneRotMatrix(RPoint vz,DOUBLE alpha,MATRIX R) {
	RPoint vx, vy;
	//vy
	vx.set(1,0,0);
	vy = cross(vz,vx);
	if(EQUAL(vy.magnitude(),0)) {
		vx.set(0,1,0);
	}
	vx = rotate(vx,vz.unitvector(),alpha);
	
	vy = cross(vz,vx);

	//vx
	vx = cross(vy,vz);

	FormRotMatrix(vx,vy,vz,R);
}
