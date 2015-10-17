#include <math.h>
#include "rpoint.h"

RPoint::RPoint() {
	x = 0;y = 0;z = 0;
}
RPoint::RPoint(DOUBLE x1,DOUBLE y1) {
	x = x1;y = y1;
}
RPoint::RPoint(DOUBLE x1,DOUBLE y1, DOUBLE z1) {
	x = x1;y = y1;z = z1;
}
void RPoint::set(DOUBLE x1,DOUBLE y1, DOUBLE z1) {
	x = x1;y = y1;z = z1;
}
DOUBLE RPoint::magnitude() {
	return sqrt(x * x + y * y + z * z);
}
RPoint RPoint::unitvector() {
	DOUBLE d = magnitude();
	return (*this * (1.0 / d));
}
int ptinbetween(RPoint& p1,RPoint& p2,RPoint& p) {
	if(((p.x - p1.x ) * ( p.x - p2.x) <= 0 )
		&& ((p.y - p1.y ) * ( p.y - p2.y) <= 0 )
		&& ((p.z - p1.z ) * ( p.z - p2.z) <= 0 ))
		return 1;
	return 0;
}
DOUBLE distance(RPoint& p,RPoint& q) {
	return sqrt(pow(p.x - q.x, 2) + pow(p.y - q.y, 2) + pow(p.z - q.z, 2));
}
DOUBLE dot(RPoint& p,RPoint& q) {
	return (p.x * q.x + p.y * q.y + p.z * q.z);
}
RPoint cross(RPoint& p,RPoint& q) {
	RPoint rp;
	rp.x = p.y * q.z - p.z * q.y;
	rp.y = p.z * q.x - p.x * q.z;
	rp.z = p.x * q.y - p.y * q.x;
	return rp;
}
int OnSameLine(RPoint& p1,RPoint& p2,RPoint& p3) {
	RPoint rp = cross(p2 - p1,p3 - p1);
	if(EQUAL(rp.x,0) && EQUAL(rp.y,0) && EQUAL(rp.z,0)) return true;
	return false;
}
int Perpendicular(RPoint& p1,RPoint& p2,RPoint& p3) {
	if(dot(p2 - p1,p3 - p1) == 0) return true;
	return false;
}
int operator == (const RPoint& left,const RPoint& right) {
	return (EQUAL(left.x,right.x) && EQUAL(left.y,right.y) && EQUAL(left.z,right.z));
}
int operator > (const RPoint& left,const RPoint& right) {
	if(EQUAL(left.x,right.x)) {
        if(EQUAL(left.y,right.y)) {
			if(EQUAL(left.z,right.z)) {
			} else  if((left.z - right.z) > SMALLNUMBER) {
				return true;
			}
		} else  if((left.y - right.y) > SMALLNUMBER) {
			return true;
		}
	} else if((left.x - right.x) > SMALLNUMBER) {
		return true;
	}
	return false;
}
int operator < (const RPoint& left,const RPoint& right) {
	if(EQUAL(left.x,right.x)) {
        if(EQUAL(left.y,right.y)) {
			if(EQUAL(left.z,right.z)) {
			} else  if((left.z - right.z) < SMALLNUMBER) {
				return true;
			}
		} else  if((left.y - right.y) < SMALLNUMBER) {
			return true;
		}
	} else if((left.x - right.x) < SMALLNUMBER) {
		return true;
	}
	return false;
}
RPoint operator + (const RPoint& left,const RPoint& right) {
	RPoint result;
	result.x = left.x + right.x;
	result.y = left.y + right.y;
	result.z = left.z + right.z;
	return result;
}
RPoint operator - (const RPoint& left,const RPoint& right) {
	RPoint result;
	result.x = left.x - right.x;
	result.y = left.y - right.y;
	result.z = left.z - right.z;
	return result;
}
RPoint operator * (const RPoint& left,const DOUBLE& right) {
	RPoint result;
	result.x = left.x * right;
	result.y = left.y * right;
	result.z = left.z * right;
	return result;
}
RPoint rotate(RPoint p, RPoint v,DOUBLE q) {
    RPoint p1;
	p1.x = p.x * (cos(q) + (1 - cos(q)) * v.x * v.x) + 
		   p.y * ((1 - cos(q)) * v.x * v.y - sin(q) * v.z) + 
		   p.z * ((1 - cos(q)) * v.x * v.z + sin(q) * v.y);

    p1.y = p.x * ((1 - cos(q)) * v.x * v.y + sin(q) * v.z) + 
		   p.y * (cos(q) + (1 - cos(q)) * v.y * v.y) + 
		   p.z * ((1 - cos(q)) * v.y * v.z - sin(q) * v.x);

	p1.z = p.x * ((1 - cos(q)) * v.z * v.x - sin(q) * v.y) + 
		   p.y * ((1 - cos(q)) * v.z * v.y + sin(q) * v.x) + 
		   p.z * (cos(q) + (1 - cos(q)) * v.z * v.z);
	return p1;
}
RPoint perspective(RPoint p, DOUBLE f) {
	RPoint p1;
	DOUBLE w = (p.z / f + 1);
	p1.x = p.x / w;
	p1.y = p.y / w;
	p1.z = 0;
	return p1;
}
RPoint aperspective(RPoint p,DOUBLE f) {
    RPoint p1;
	DOUBLE w = (p.z / f + 1);
	p1.x = w * p.x;
	p1.y = w * p.y;
	p1.z = p.z;
	return p1;
}
RPoint ptprojection(RPoint& p,RPoint& p1,RPoint& p2) {
	DOUBLE L = distance(p1,p2),L1;
	RPoint e = p2 - p1,v,r;
	e = e.unitvector();
	v = p - p1;
    L1 = dot(v,e);
    r = p1 + (p2 - p1) * (L1 / L);
	return r;
}
DOUBLE triangle_area(RPoint& p1,RPoint& p2,RPoint& p3) {
	RPoint A = cross(p3 - p1,p2 - p1);
	return fabs(0.5 * A.magnitude());
}
DOUBLE quad_area(RPoint& p1,RPoint& p2,RPoint& p3,RPoint& p4) {
	return (triangle_area(p1,p2,p4) + triangle_area(p2,p3,p4));
}
int clockwise(const RPoint p1,const RPoint p2,const RPoint p3) {
	RPoint v = cross(p2 - p1,p3 - p1);
	return ((v.x + v.y + v.z) >= 0);
}
bool Intersect(RPoint p1,RPoint p2,RPoint p3,RPoint p4,RPoint& p) {
	DOUBLE m1,m2,x,y;
	if(!EQUAL(p1.x, p2.x)) m1 = DOUBLE(p2.y - p1.y) / (p2.x - p1.x);
	if(!EQUAL(p3.x, p4.x)) m2 = DOUBLE(p4.y - p3.y) / (p4.x - p3.x);
	
	if(EQUAL(p1.x, p2.x))  {
		x = p1.x;
		y = m2 * x + (p3.y - m2 * p3.x);
	} else if(EQUAL(p3.x, p4.x))  {
		x = p3.x;
		y = m1 * x + (p1.y - m1 * p1.x);
	} else {
		if(EQUAL(m1,m2)) 
			return false;
		x = ((p3.y - m2 * p3.x) - (p1.y - m1 * p1.x)) / (m1 - m2);
		y = m1 * x + (p1.y - m1 * p1.x);
	}

    if( ((x >= p1.x && x <= p2.x) || (x >= p2.x && x <= p1.x)) &&
        ((y >= p1.y && y <= p2.y) || (y >= p2.y && y <= p1.y)) &&
	    ((x >= p3.x && x <= p4.x) || (x >= p4.x && x <= p3.x)) &&
        ((y >= p3.y && y <= p4.y) || (y >= p4.y && y <= p3.y)) )  {
        p.x = x;
		p.y = y;
	    return true;
	}

	return false;
}
double myatan(double y,double x) {
    double d = atan(fabs(y / x));
	if(y >= 0) {
		if(x >= 0) return d;
		else return PI - d;
	} else {
		if(x >= 0) return 2 * PI - d;
		else return d + PI;
	}
}
RPoint RDtoCT(RPoint& rp) {
	RPoint myp;
	myp.x = rp.x * cos(rp.y * (PI / 180)); 
	myp.y = rp.x * sin(rp.y * (PI / 180));
	myp.z = rp.z; 
	return myp;
}
RPoint CTtoRD(RPoint& rp) {
	RPoint myp;
	myp.x = sqrt(rp.x * rp.x + rp.y * rp.y); 
	myp.y = 180 * myatan(rp.y , rp.x) / PI;
	myp.z = rp.z; 
    return myp;
}