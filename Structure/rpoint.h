#ifndef __RPOINT_H
#define __RPOINT_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

#undef DOUBLE
typedef double DOUBLE;

#undef pow
#define pow(x,y) pow((double)x,(double)y)
#undef log
#define log(x) log((double)x)
#undef sqrt
#define sqrt(x) sqrt((double)x)


#define INVALID                          -1
#define PI                    3.14159265359
#define SMALLNUMBER                   1e-10
#define VERYSMALLNUMBER               1e-14
#define EXTREMESMALLNUMBER            1e-18
#define EQUAL(x,y)         (fabs(DOUBLE((x) - (y))) <= SMALLNUMBER)
#define VEQUAL(x,y)        (fabs(DOUBLE((x) - (y))) <= VERYSMALLNUMBER)
#define EEQUAL(x,y)        (fabs(DOUBLE((x) - (y))) <= EXTREMESMALLNUMBER)

/*
RPoint: a point in 3d
*/
struct RPoint {
	DOUBLE x,y,z;
	RPoint();
	RPoint(DOUBLE x1,DOUBLE y1);
	RPoint(DOUBLE x1,DOUBLE y1, DOUBLE z1);
	void set(DOUBLE x1,DOUBLE y1, DOUBLE z1);
	DOUBLE magnitude();
    RPoint unitvector();
	friend RPoint ptprojection(RPoint& p,RPoint& p1,RPoint& p2);
	friend int ptinbetween(RPoint& p1,RPoint& p2,RPoint& p);
	friend DOUBLE distance(RPoint& p,RPoint& q);
	friend DOUBLE dot(RPoint& p,RPoint& q);
	friend RPoint cross(RPoint& p,RPoint& q);
	friend int OnSameLine(RPoint& p1,RPoint& p2,RPoint& p3);
	friend int Perpendicular(RPoint& p1,RPoint& p2,RPoint& p3);
	friend RPoint rotate(RPoint p, RPoint unitv,DOUBLE q);
	friend RPoint perspective(RPoint p, DOUBLE f);
	friend RPoint aperspective(RPoint p,DOUBLE f);
	friend DOUBLE triangle_area(RPoint&,RPoint&,RPoint&);
    friend DOUBLE quad_area(RPoint&,RPoint&,RPoint&,RPoint&);
	friend int clockwise(const RPoint,const RPoint,const RPoint);
	friend int operator == (const RPoint& left,const RPoint& right);
	friend int operator > (const RPoint& left,const RPoint& right);
	friend int operator < (const RPoint& left,const RPoint& right);
	friend RPoint operator + (const RPoint& left,const RPoint& right);
	friend RPoint operator - (const RPoint& left,const RPoint& right);
	friend RPoint operator * (const RPoint& left,const DOUBLE& right);
	friend bool Intersect(RPoint p1,RPoint p2,RPoint p3,RPoint p4,RPoint& p);
	friend RPoint RDtoCT(RPoint& rp);
    friend RPoint CTtoRD(RPoint& rp);
};

double myatan(double y,double x);

#endif