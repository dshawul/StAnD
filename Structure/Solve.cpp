#include "common.h"

UBMP32 ENTITY::TotalLoadCases;
UBMP32 JOINT::TotalNumber;
UBMP32 MEMBER::TotalNumber;
UBMP32 SLAB::TotalNumber;
BOOL MEMBER::add_geometric_stiffness = FALSE;
BOOL MEMBER::only_geometric_stiffness = FALSE;
int MEMBER::p_delta;
const UINT SLAB::NSTRESS = 23;
PLOADCOMBOLIST BASE_LOAD::cloadcombolist;
DOUBLE BASE_LOAD::factor;
DOUBLE BASE_LOAD::load_step_factor = 1;
DETAILING* SECTION::pDocDetailing;

/*
Allocate memory for analysis
*/
void JOINT::AllocVectors() {
	forces_all = new VECTOR[TotalLoadCases];
    disps_all = new VECTOR[TotalLoadCases];
}
void JOINT::FreeVectors() {
	delete[] forces_all;
    delete[] disps_all;
}
void MEMBER::AllocVectors() {
	forces_all = new MATRIX[TotalLoadCases];
	disps_all = new MATRIX[TotalLoadCases];
	for(UINT i = 0;i < TotalLoadCases;i++) {
		forces_all[i] = 0;
		disps_all[i] = 0;
	}
	vec_alloc(station,nDiv);
	design_data = new DESIGNDATA[nDiv];
}
void MEMBER::FreeVectors() {
	delete[] forces_all;
	delete[] disps_all;
	vec_free(station);
	delete[] design_data;
}
void SLAB::AllocVectors() {
	forces_all = new VECTOR[TotalLoadCases];
    disps_all = new VECTOR[TotalLoadCases];
}
void SLAB::FreeVectors() {
	delete[] forces_all;
	delete[] disps_all;
}
void JOINT::AllocAnalysis(int index) {
	vec_alloc(forces_all[index],6); 
	vec_alloc(disps_all[index],6); 
}
void JOINT::FreeAnalysis(int index) {
	vec_free(forces_all[index]); 
	vec_free(disps_all[index]); 
}
void MEMBER::AllocAnalysis(int index) {
	mat_alloc(forces_all[index],nDiv,6); 
    mat_alloc(disps_all[index],nDiv,6); 
}
void MEMBER::FreeAnalysis(int index) {
	vec_free(forces_all[index]); 
    vec_free(disps_all[index]); 
}
void SLAB::AllocAnalysis(int index) {
	vec_alloc(forces_all[index],NSTRESS * NJ); 
    vec_alloc(disps_all[index],6 * NJ); 
}
void SLAB::FreeAnalysis(int index) {
    vec_free(forces_all[index]); 
    vec_free(disps_all[index]); 
}
/*
SECTION
*/
static DOUBLE GetI(int type,int index,DOUBLE w,DOUBLE h,DOUBLE r) {
	if(type == RECTANGULAR) {
		DOUBLE b;
		if(index == IUX) {
			if(w <= h) {
				b = 1 / 3.0 - (0.21f * w / h) * (1 - pow(w,4) / (12 * pow(h,4)));
				return (b * h * pow(w,3));
			} else {
				b = 1 / 3.0 - (0.21f * h / w) * (1 - pow(h,4) / (12 * pow(w,4)));
				return (b * w * pow(h,3));
			}
		} else if(index == IUY) {
			return ((w * pow(h,3)) / 12);
		} else if(index == IUZ) {
			return ((h * pow(w,3)) / 12);
		} else if(index == IRY) {
			return ((w * pow(h,2)) / 4);
		} else if(index == IRZ) {
			return ((h * pow(w,2)) / 4);
		}
	} else {
		if(index == IUX) {return (PI * pow(r,4) / 2);}
		else if(index == IUY) {return (PI * pow(r,4) / 4);}
		else if(index == IUZ) {return (PI * pow(r,4) / 4);}
		else if(index == IRY) {return (4 * pow(r,3) / 3);}
		else if(index == IRZ) {return (4 * pow(r,3) / 3);}
	}
	return 0;
}
static DOUBLE GetS(DOUBLE w,DOUBLE h) {
	return ((w * pow(h,2)) / 2);
}

void SECTION::GetData() {

	if(type == RECTANGULAR) {
	    A = w * h;

		Ay = A;
		Az = A;

        Ix = GetI(RECTANGULAR,IUX,w,h,r);
		Iy = GetI(RECTANGULAR,IUY,w,h,r);
		Iz = GetI(RECTANGULAR,IUZ,w,h,r);
		Zy = Iy / (h / 2);
		Zz = Iz / (w / 2);
		Sy = GetI(RECTANGULAR,IRY,w,h,r);
		Sz = GetI(RECTANGULAR,IRZ,w,h,r);
	} else if(type == CIRCULAR) {
		A = PI * pow(r,2);

		Ay = A;
		Az = A;

		Ix = GetI(CIRCULAR,IUX,w,h,r);
		Iy = GetI(CIRCULAR,IUY,w,h,r);
		Iz = GetI(CIRCULAR,IUZ,w,h,r);
		Zy = Iy / (r);
		Zz = Iz / (r);
		Sy = GetI(CIRCULAR,IRY,w,h,r);
		Sz = GetI(CIRCULAR,IRZ,w,h,r);
	} else if(type == WIDEFLANGE) {
		A = w * (h + tf) - (w - tw) * (h - tf);

		if(built_up) {
			Ay = 2 * tf * w;
			Az = A - Ay;
		} else {
            Ay = 2 * tf * w;
			Az = A - 2 * tf * w + (tw) * tf;
		}

		Ix = 2 * GetI(RECTANGULAR,IUX,w,tf,r) + GetI(RECTANGULAR,IUX,tw,h - tf,r);
		Iy = GetI(RECTANGULAR,IUY,w,h + tf,r) - GetI(RECTANGULAR,IUY,w - tw,h - tf,r);
		Iz = GetI(RECTANGULAR,IUZ,w,h + tf,r) - (GetI(RECTANGULAR,IUZ,w,h - tf,r) - GetI(RECTANGULAR,IUZ,tw,h - tf,r));
		Zy = Iy / ((h + tf) / 2);
		Zz = Iz / (w / 2);
		Sy = GetI(RECTANGULAR,IRY,w,h + tf,r) - GetI(RECTANGULAR,IRY,w - tw,h - tf,r);
		Sz = GetI(RECTANGULAR,IRZ,w,h + tf,r) - (GetI(RECTANGULAR,IRZ,w,h - tf,r) - GetI(RECTANGULAR,IRZ,tw,h - tf,r));
	} else if(type == CHANNEL) {
		DOUBLE x,hm,bm;
		hm = h + tf;
		bm = w + tw / 2;
		A = bm * hm - (bm - tw) * (hm - 2 * tf);

		if(built_up) {
			Ay = 2 * tf * w;
			Az = A - Ay;
		} else {
            Ay = 2 * tf * w;
			Az = A - 2 * tf * w + (tw) * tf;
		}

		x = (2 * bm * bm * tf + (hm - 2 * tf) * tw * tw) / (2 * A);
		Ix = 2 * GetI(RECTANGULAR,IUX,bm,tf,r) + GetI(RECTANGULAR,IUX,tw,hm - 2 * tf,r);
		Iy = GetI(RECTANGULAR,IUY,bm,hm,r) - GetI(RECTANGULAR,IUY,bm - tw,hm - 2 * tf,r);
		Iz = (2 * tf * pow(bm,3) + (hm - 2 * tf) * pow(tw,3)) / 3 - (A * pow(x,2));
		Zy = Iy / (hm / 2);
		Zz = Iz / (bm - x);
		Sy = GetI(RECTANGULAR,IRY,bm,hm,r) - GetI(RECTANGULAR,IRY,bm - tw,hm - 2 * tf,r);
		x = A / (4 * tf);
		if(x > bm - tw) x = bm - tw + (A / 2 - 2 * (bm - tw) * tf) / hm;
		Sz = 2 * GetS(tf,x) + GetS(hm,bm - x) - GetS(hm - 2 * tf,bm - x - tw);
    } else if(type == DOUBLE_CHANNEL) {
		DOUBLE hm,bm;
		hm = h + tf;
		bm = w + tw / 2;
		A = 2 * (bm * hm - (bm - tw) * (hm - 2 * tf));

		if(built_up) {
			Ay = 4 * tf * w;
			Az = A - Ay;
		} else {
            Ay = 4 * tf * w;
			Az = A - 4 * tf * w + (2 * tw) * tf;
		}

		Ix = 2 * (2 * GetI(RECTANGULAR,IUX,bm,tf,r) + GetI(RECTANGULAR,IUX,tw,hm - 2 * tf,r));
		Iy = 2 * (GetI(RECTANGULAR,IUY,bm,hm,r) - GetI(RECTANGULAR,IUY,bm - tw,hm - 2 * tf,r));
		Iz = GetI(RECTANGULAR,IUZ,2 * bm + bd,2 * tf,r) - GetI(RECTANGULAR,IUZ,bd,2 * tf,r) +
			 GetI(RECTANGULAR,IUZ,2 * tw + bd,hm - tf,r) - GetI(RECTANGULAR,IUZ,bd,hm - tf,r);
		Zy = Iy / (hm / 2);
		Zz = Iz / (bm);
		Sy = 2 * (GetI(RECTANGULAR,IRY,bm,hm,r) - GetI(RECTANGULAR,IRY,bm - tw,hm - 2 * tf,r));
		Sz = GetI(RECTANGULAR,IRZ,2 * bm + bd,2 * tf,r) - GetI(RECTANGULAR,IRZ,bd,2 * tf,r) +
			 GetI(RECTANGULAR,IRZ,2 * tw + bd,hm - tf,r) - GetI(RECTANGULAR,IRZ,bd,hm - tf,r);
	} else if(type == TEE) {
		DOUBLE y,hm,bm;
		hm = h + tf / 2;
		bm = w;
		A = bm * hm - (bm - tw) * (hm - tf);

		Ay = tf * bm;
        Az = tw * hm;

		y = hm - (pow(hm, 2) * tw + pow(tf,2) * (bm - tw)) / (2 * A);
		Ix = GetI(RECTANGULAR,IUX,bm,tf,r) + GetI(RECTANGULAR,IUX,tw,hm - tf,r);
		Iy = (tw * pow(y,3) + bm * pow(hm - y, 3) - (bm - tw) * pow(hm - tf - y, 3)) / 3;
		Iz = GetI(RECTANGULAR,IUZ,bm,hm,r) - (GetI(RECTANGULAR,IUZ,bm,hm - tf,r) - GetI(RECTANGULAR,IUZ,tw,hm - tf,r));
		Zy = Iy / (y);
		Zz = Iz / (bm / 2);
		y = A / (2 * tw);
		if(y > hm - tf) y = hm - tf + (A / 2 - (hm - tf) * tw) / bm;
		Sy = GetS(tw,y) + GetS(bm,hm - y) - GetS(bm - tw,hm - y - tf);
		Sz = GetI(RECTANGULAR,IRZ,bm,hm,r) - (GetI(RECTANGULAR,IRZ,bm,hm - tf,r) - GetI(RECTANGULAR,IRZ,tw,hm - tf,r));
	} else if(type == ANGLE) {
		DOUBLE x,y,hm,bm;
		hm = h + tf / 2;
		bm = w + tw / 2;
		A = bm * hm - (bm - tw) * (hm - tf);

		Ay = tf * bm;
        Az = tw * hm;

		x = bm - (tw * (2 * (bm - tw) + hm) + pow(bm - tw,2)) / (2 * (bm - tw + hm));
        y = hm - (tf * (2 * (hm - tf) + bm) + pow(hm - tf,2)) / (2 * (hm - tf + bm)); 
		Ix = GetI(RECTANGULAR,IUX,bm,tf,r) + GetI(RECTANGULAR,IUX,tw,hm - tf / 2,r);
		Iy = (tw * pow(y,3) + bm * pow(hm - y, 3) - (bm - tw) * pow(hm - tf - y, 3)) / 3;
		Iz = (tf * pow(x,3) + hm * pow(bm - x, 3) - (hm - tf) * pow(bm - tw - x, 3)) / 3;
		Zy = Iy / (y);
		Zz = Iz / (x);
		x = A / (2 * tf);
		if(x > bm - tw) x = bm - tw + (A/2 - (bm - tw) * tf) / (hm);
		Sz = GetS(tf,x) + GetS(hm,bm - x) - GetS(hm - tf,bm - x - tw);
		y = A / (2 * tw);
		if(y > hm - tf) y = hm - tf + (A/2 - (hm - tf) * tw) / (bm);
		Sy = GetS(tw,y) + GetS(bm,hm - y) - GetS(bm - tw,hm - y - tf);
	} else if(type == DOUBLE_ANGLE) {
		DOUBLE y,hm,bm;
		hm = h + tf / 2;
		bm = w;
		A = 2 * (bm * hm - (bm - tw) * (hm - tf));

		Ay = 2 * tf * bm;
        Az = 2 * tw * hm;

        y = hm - (tf * (2 * (hm - tf) + bm) + pow(hm - tf,2)) / (2 * (hm - tf + bm)); 
		Ix = 2 * (GetI(RECTANGULAR,IUX,bm,tf,r) + GetI(RECTANGULAR,IUX,tw,hm - tf / 2,r));
		Iy = 2 * (tw * pow(y,3) + bm * pow(hm - y, 3) - (bm - tw) * pow(hm - tf - y, 3)) / 3;
		Iz = GetI(RECTANGULAR,IUZ,2 * bm + bd,tf,r) - GetI(RECTANGULAR,IUZ,bd,tf,r) +
			 GetI(RECTANGULAR,IUZ,2 * tw + bd,hm - tf,r) - GetI(RECTANGULAR,IUZ,bd,hm - tf,r);
		Zy = Iy / (y);
		Zz = Iz / (bm);

		y = A / (4 * tw);
		if(y > hm - tf) y = hm - tf + (A/2 - 2 * (hm - tf) * tw) / (2 * bm);
		Sy = GetS(2 * tw,y) + GetS(2 * bm,hm - y) - GetS(2 * bm - 2 * tw,hm - y - tf);
		Sz = GetI(RECTANGULAR,IRZ,2 * bm + bd,tf,r) - GetI(RECTANGULAR,IRZ,bd,tf,r) +
			 GetI(RECTANGULAR,IRZ,2 * tw + bd,hm - tf,r) - GetI(RECTANGULAR,IRZ,bd,hm - tf,r);
	} else if(type == BOX) {
		DOUBLE hm,bm;
		hm = h + tf;
		bm = w + tw;
		A = hm * bm - (hm - 2 * tf) * (bm - 2 * tw);

		Ay = A * bm / (hm + bm);
		Az = A * hm / (hm + bm);

		Ix = 2 * pow(h,2) * pow(w,2) * tf * tw /(w * tw + h * tf);
		Iy = GetI(RECTANGULAR,IUY,w + tw,h + tf,r) - GetI(RECTANGULAR,IUY,w - tw,h - tf,r );
		Iz = GetI(RECTANGULAR,IUZ,w + tw,h + tf,r) - GetI(RECTANGULAR,IUZ,w - tw,h - tf,r );
		Zy = Iy / ((h + tf) / 2);
		Zz = Iz / ((w + tw) / 2);
		Sy = GetI(RECTANGULAR,IRY,w + tw,h + tf,r) - GetI(RECTANGULAR,IRY,w - tw,h - tf,r );
		Sz = GetI(RECTANGULAR,IRZ,w + tw,h + tf,r) - GetI(RECTANGULAR,IRZ,w - tw,h - tf,r );
	} else if(type == PIPE) {
		DOUBLE rm;
		rm = r + tw / 2;
		A = PI * pow(rm,2) - PI * pow(rm - tw,2);

		Ay = 2 * A / PI;
		Az = 2 * A / PI;

		Ix = GetI(CIRCULAR,IUX,w,h,r + tw / 2) - GetI(CIRCULAR,IUX,w,h,r - tw / 2);
		Iy = GetI(CIRCULAR,IUY,w,h,r + tw / 2) - GetI(CIRCULAR,IUY,w,h,r - tw / 2);
		Iz = GetI(CIRCULAR,IUZ,w,h,r + tw / 2) - GetI(CIRCULAR,IUZ,w,h,r - tw / 2);
		Zy = Iy / (r + tw/2);
		Zz = Iz / (r + tw/2);
		Sy = GetI(CIRCULAR,IRY,w,h,r + tw / 2) - GetI(CIRCULAR,IRY,w,h,r - tw / 2);
		Sz = GetI(CIRCULAR,IRZ,w,h,r + tw / 2) - GetI(CIRCULAR,IRZ,w,h,r - tw / 2);
	}
	ry = sqrt(Iy / A);
    rz = sqrt(Iz / A);
}
/*
Loads
*/
void JLOAD::ApplyFactor() {
    for(int i = IUX;i <= IRZ;i++) 
		Q[i] *= factor;
}
void JLOAD::RemoveFactor() {
	for(int i = IUX;i <= IRZ;i++) 
		Q[i] /= factor;
}
void LOAD::ApplyFactor() {
	P *= factor;
	P1 *= factor;
}
void LOAD::RemoveFactor() {
	P /= factor;
	P1 /= factor;
}
BOOL BASE_LOAD::is_considered() {
	if(!cloadcombolist)
		return FALSE;

	LOADCOMBO* pldcombo;
	POSITION pos = cloadcombolist->GetHeadPosition();
	while(pos) {
		pldcombo = &cloadcombolist->GetNext(pos);
		if(EQUAL(pldcombo->FS,0))
			continue;
		if(pldcombo->loadcase == loadcase) {
			factor = (pldcombo->FS * BASE_LOAD::load_step_factor);
			return TRUE;
		}
	}
	return FALSE;
}
/*
JOINT class
*/

void JOINT::ApplyLoad(VECTOR Qg,VECTOR Dg) {
	JLOAD* pload;
	POSITION pos;
	VECTOR vec;

	MATRIX R,R1,R2;

	mat_alloc(R,3,3);
	mat_alloc(R1,3,3);
	mat_alloc(R2,3,3);
	GetRotMatrix(R);

	pos = load.GetHeadPosition();
	while(pos) {
        pload = &load.GetNext(pos);
		if(pload->type == MASS) continue;

		/*joint load*/
		if(!pload->is_joint_load())
			continue;

		/*apply factor*/
		if(!pload->is_considered())
			continue;

		pload->ApplyFactor();
		/*rotate to local axis*/
		if(pload->system) {
			GetPointRotMatrix(pload->system->rotation,R1);
			multiply(R,R1,R2,3,3,3);
			equ(R,R2,3,3);
		}

		/*add to Q/D*/
		vec = (pload->type == FORCE) ? Qg : Dg;
		if(!pload->system) {
			for(int i = IUX; i <= IRZ;i++) {
				if(pload->type == DISPLACEMENT && !(restraint & (1 << i))) continue;
				if(number[i] != INVALID) vec[number[i]] += pload->Q[i];
			}
		} else {
			RPoint p(pload->Q[IUX],pload->Q[IUY],pload->Q[IUZ]);
			p = GetGlobal(p,R);
			if(pload->type == DISPLACEMENT) {
			   if(!(restraint & (1 << IUX))) p.x = 0;
               if(!(restraint & (1 << IUY))) p.y = 0;
			   if(!(restraint & (1 << IUZ))) p.z = 0;
			}
			if(number[IUX] != INVALID) vec[number[IUX]] += p.x;
			if(number[IUY] != INVALID) vec[number[IUY]] += p.y;
			if(number[IUZ] != INVALID) vec[number[IUZ]] += p.z;
			p.set(pload->Q[IRX],pload->Q[IRY],pload->Q[IRZ]);
			p = GetGlobal(p,R);
			if(pload->type == DISPLACEMENT) {
			   if(!(restraint & (1 << IRX))) p.x = 0;
               if(!(restraint & (1 << IRY))) p.y = 0;
			   if(!(restraint & (1 << IRZ))) p.z = 0;
			}
			if(number[IRX] != INVALID) vec[number[IRX]] += p.x;
			if(number[IRY] != INVALID) vec[number[IRY]] += p.y;
			if(number[IRZ] != INVALID) vec[number[IRZ]] += p.z;
		}
		/*remove factor*/
		pload->RemoveFactor();
	}
	
	vec_free(R);
	vec_free(R1);
	vec_free(R2);
}
void JOINT::GetRotMatrix(MATRIX T) {
	GetPointRotMatrix(rotation,T);
}
void ELEMENT::AssembleJoints(JOINT** jt,int NJ,MATRAN K,MATRIX Kl,int NB) {
	UINT i,j,x,y,TOTAL = 6 * NJ;
	for(i = 0;i < TOTAL;i++) {
		x = jt[i / 6]->number[i % 6]; 
		if(x == INVALID) continue;

		for(j = i;j < TOTAL;j++) {
			y = jt[j / 6]->number[j % 6];
			if(y == INVALID) continue;
			
			K[BANDI(x,y)] += Kl[i][j];
		}
	}
}
/*
RIGID BODY
*/
DOUBLE CONSTRAINT::PenaltyWeight;

bool CONSTRAINT::DetermineAxis(RPoint& rp,RPoint centroid) {
	MATRIX I;
	RPoint p;
    POSITION pos;
	JOINT* joint;

	mat_alloc(I,3,3);

	/*second moment tensor*/
	pos = jointlist.GetHeadPosition();
	while(pos) {
		joint = jointlist.GetNext(pos);
		p = joint->p - centroid;
		I[0][0] += pow(p.y,2) + pow(p.z,2);
        I[1][1] += pow(p.x,2) + pow(p.z,2);
		I[2][2] += pow(p.x,2) + pow(p.y,2);
		I[0][1] += -(p.x) * (p.y);	
		I[0][2] += -(p.x) * (p.z);
		I[1][2] += -(p.y) * (p.z);
	}
	I[1][0] = I[0][1];
	I[2][0] = I[0][2];
	I[2][1] = I[1][2];
	
	/*find moment and directions*/
	DOUBLE sigma[3];
	DOUBLE dir[9];
	find_eigen(I,3,sigma,dir);

	/*largest eigen value*/
	int i,maxi,mini,index;
	DOUBLE maxv = -1e16,minv = 1e16;

	for(i = 0;i < 3;i++) {
		if(fabs(sigma[i]) > maxv) {
			maxv = fabs(sigma[i]);
			maxi = i;
		}
		if(fabs(sigma[i]) < minv) {
			minv = fabs(sigma[i]);
		    mini = i;
		}
	}

	if(type == DIAPHRAGM || type == PLATE) index = maxi;
	else index = mini;

	/*not unique*/
	int count = 0;
	for(i = 0;i < 3;i++) {
		if(EQUAL(sigma[i],sigma[index]))
			count++;	
	}
	if(count > 1) {
		vec_free(I);
        return false;
	}

	/*set rotation*/
	rp.x = dir[index * 3 + 0];
	rp.y = dir[index * 3 + 1];
    rp.z = dir[index * 3 + 2];

    vec_free(I);
    return true;
}

void CONSTRAINT::AddToK(MATRAN K,UINT NB) {
	MATRIX R,Rm1,Rm2,Rs,T,Tl,Cl,Cg;
	VECTOR v,va,vb;
	JOINT *joint,*slave,*master1,*master2;
	POSITION pos,pos1;
	UINT i,NF,count;
	RPoint rp,delta;
	DOUBLE f1,f2;

	count = jointlist.GetCount();
	if(count <= 1)
		return;

	/*find close joints and weld them*/
	if(type == WELD) {
		DOUBLE dist;

		pos = jointlist.GetHeadPosition();
		while(pos) {
			master1 = jointlist.GetNext(pos);
			
			pos1 = jointlist.GetHeadPosition();
			while(pos1) {
				slave = jointlist.GetNext(pos1);
				
				dist = distance(master1->p,slave->p);
				if(dist <= tolerance) {
					CONSTRAINT cons;
					cons.type = BODY;
					cons.constraint = constraint;
					cons.axis = axis;
					cons.jointlist.AddTail(master1);
					cons.jointlist.AddTail(slave);
					AddToK(K,NB);
				}
			}
		}
		return;
	}

	/*alloc*/
	if(type == LINE) NF = 18;
	else NF = 12;

	vec_alloc(v,NF);
	vec_alloc(va,NF);
	vec_alloc(vb,NF);
	mat_alloc(Cl,NF,NF);
	mat_alloc(Cg,NF,NF);
	mat_alloc(T,NF,NF);
	mat_alloc(Tl,NF,NF);
	mat_alloc(R,3,3);
	mat_alloc(Rm1,3,3);
	mat_alloc(Rm2,3,3);
	mat_alloc(Rs,3,3);

	/*centroid*/
	RPoint centroid;
	pos = jointlist.GetHeadPosition();
	while(pos) {
		joint = jointlist.GetNext(pos);
		centroid = centroid + joint->p;
	}
	centroid = centroid * (1.0 / jointlist.GetCount());

	/*rotation matrix*/
	if(isauto && DetermineAxis(rp,centroid)) {
		if(type == DIAPHRAGM || type == PLATE) GetPlaneRotMatrix(rp,0,R);
		else if(type == BEAM || type == ROD) GetLineRotMatrix(rp,0,R);
	} else {
		GetPointRotMatrix(axis->rotation,R);
		if(type == EQUAL || type == LOCAL || type == BODY || type == LINE);
		else {
			if(constraint & UX) rp.set(1,0,0);
			else if(constraint & UY) rp.set(0,1,0);
			else if(constraint & UZ) rp.set(0,0,1);

			if(type == DIAPHRAGM || type == PLATE) GetPlaneRotMatrix(rp,0,Rm1);
			else if(type == BEAM || type == ROD) GetLineRotMatrix(rp,0,Rm1);
			multiply(R,Rm1,Rm2,3,3,3);
			equ(R,Rm2,3,3);
		}
	}
	/*rotation to global system*/
	if(type == LOCAL) {
		for(i = 0;i < NF;i++)
			Tl[i][i] = 1;
	} else {
		equsec(Tl,R,3,3,0,0,0,0);
		equsec(Tl,R,3,3,3,3,0,0);			
		equsec(Tl,R,3,3,6,6,0,0);
		equsec(Tl,R,3,3,9,9,0,0);
		if(type == LINE) {
			equsec(Tl,R,12,12,0,0,0,0);
			equsec(Tl,R,15,15,3,3,0,0);	
		}
		transpose(Tl,NF,NF);
	}

	/*determine masters*/
	master1 = master2 = 0;
	if(type == LINE) {
		DOUBLE dist,dist1 = 0,dist2 = 0;
		pos = jointlist.GetHeadPosition();
		while(pos) {
			joint = jointlist.GetNext(pos);
			dist = distance(joint->p,centroid);
			if(dist > dist1) {
				master2 = master1;
				master1 = joint;
				dist2 = dist1;
				dist1 = dist2;
			} else if(dist > dist2) {
				master2 = joint;
				dist2 = dist;
			}
		}
		master1->GetRotMatrix(Rm1);
		master2->GetRotMatrix(Rm2);
	} else {
		master1 = jointlist.GetHead();
		master1->GetRotMatrix(Rm1);
	}

	/*constrained dofs*/
	UINT myconstraint;
	if(type == BEAM) myconstraint = UZ | RY | RZ;
	else if(type == ROD) myconstraint = UX;
	else if(type == DIAPHRAGM) myconstraint = UX | UY | RZ;
	else if(type == PLATE) myconstraint = RX | RY | UZ;
	else myconstraint = constraint;

	/*fore each slave*/
	pos = jointlist.GetHeadPosition();
	while(pos) {
		slave = jointlist.GetNext(pos);
		if(master1 && master1 == slave) continue;
        if(master2 && master2 == slave) continue;

		slave->GetRotMatrix(Rs);
		delta.x = slave->p.x - master1->p.x;
		delta.y = slave->p.y - master1->p.y;
		delta.z = slave->p.z - master1->p.z;

		delta = GetGlobal(delta,Rs); //check Rs?

		clear(Cg,NF,NF);
		clear(T,NF,NF);

		if(type == LOCAL) {
			for(i = 0;i < NF;i++)
				T[i][i] = 1;
		} else {
			equsec(T,Rm1,3,3,0,0,0,0);
			equsec(T,Rm1,3,3,3,3,0,0);			
			equsec(T,Rs,3,3,6,6,0,0);
			equsec(T,Rs,3,3,9,9,0,0);
			if(type == LINE) {
				equsec(T,Rm1,12,12,0,0,0,0);
			    equsec(T,Rm1,15,15,3,3,0,0);	
			}
		}

		/*form equations for different constraints*/
		for(i = IUX;i <= IRZ;i++) {
			if(myconstraint & (1 << i)) {
				clear(v,NF);
				if(type == EQUAL) {
					v[i] = 1;
					v[i + 6] = -1;
				} else if(type == LOCAL) {
					v[i] = 1;
					v[i + 6] = -1;
				} else if(type == ROD) {
					v[i] = 1;
					v[i + 6] = -1;
				} else if(type == LINE) {
					f1 = distance(slave->p,master2->p) / distance(master1->p,master2->p);
					f2 = distance(slave->p,master1->p) / distance(master1->p,master2->p);
					v[i] = -f1;
					v[i + 6] = 1;
					v[i + 12] = -f2;
				} else if(type == BODY) {
					v[i] = 1;
					v[i + 6] = -1;
					if(i <= IUZ) {	
						if(i == IUX) {
							v[IRY] = delta.z;				
							v[IRZ] = -delta.y;
						} else if(i == IUY) {
							v[IRZ] = delta.x;
							v[IRX] = -delta.z;
						} else if(i == IUZ) {
							v[IRX] = delta.y;
						    v[IRY] = -delta.x;	
						}
					}
				} else if(type == DIAPHRAGM) {
					v[i] = 1;
					v[i + 6] = -1;
					if(i == IUX) v[IRZ] = -delta.y;
					else if(i == IUY) v[IRZ] = delta.x;
				} else if(type == PLATE) {
					v[i] = 1;
					v[i + 6] = -1;
					if(i == IUZ) {							
						v[IRX] = delta.y;
						v[IRY] = -delta.x;
					}
				} else if(type == BEAM) {
					v[i] = 1;
					v[i + 6] = -1;
					if(i == IUY) v[IRZ] = delta.x;
					else if(i == IUZ) v[IRY] = -delta.x;
				}
			}
						
			multiply(Tl,v,va,NF,NF);
			multiply(T,va,vb,NF,NF);
			multiply(vb,vb,Cl,NF,NF);
			add(Cg,Cl,NF,NF);
		}
		/*assemble values*/
		multiply(Cg,PenaltyWeight,NF,NF);

		JOINT* jt[3];
		jt[0] = master1;
		jt[1] = slave;
		count = 2;
		if(type == LINE) {
			jt[2] = master2;
			count = 3;
		}
		ELEMENT::AssembleJoints(jt,count,K,Cg,NB);
		/*end*/
	};

	/*free*/
	vec_free(v);
	vec_free(va);
	vec_free(vb);
	vec_free(Cl);
	vec_free(Cg);
	vec_free(T);
	vec_free(Tl);
	vec_free(R);
	vec_free(Rm1);
	vec_free(Rm2);
	vec_free(Rs);
}


void SLAB::AddEdgeConstraint(PJOINTLIST joints,SYSTEM* global,MATRAN K,UINT NB) {
	POSITION pos;
	JOINT* joint;
	UINT i,next;
	RPoint p1,p2;
	pos = joints->GetHeadPosition();
	while(pos) {
	    joint = &joints->GetNext(pos);
		
		for(i = 0;i < NJ;i++) {
			next = (i == NJ - 1) ? 0  : (i + 1);
			p1 = jt[i]->p;
			p2 = jt[next]->p;
			if((p1 == joint->p) || (p2 == joint->p))
				break;
			
			CONSTRAINT cons;
			cons.type = CONSTRAINT::LINE;
			cons.constraint = ALLUR;
			cons.axis = global;
			
			if(ptinbetween(p1,p2,joint->p) && OnSameLine(p1,p2,joint->p)) {
				cons.jointlist.AddTail(jt[i]);
				cons.jointlist.AddTail(jt[next]);
				cons.jointlist.AddTail(joint);
				cons.AddToK(K,NB);
			}
		}
	}
}
