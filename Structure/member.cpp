#include "common.h"

/*
MEMBER
*/
void MEMBER::GetRestraint(VECTOR d) {
	UBMP8 i;
	MATRIX T;
	VECTOR D;
	mat_alloc(T,12,12);
	vec_alloc(D,12);
	for(i = 0;i < 12;i++) {
		if(i < 6) D[i] = (j1->restraint & (1 << i));
		else D[i] = (j2->restraint & (1 << (i - 6)));
	}
	
	GetTransformationMatrix(T);
	multiply(T,D,d,12,12);
	for(i = 0;i < 12;i++) {
		if(EQUAL(d[i],0)) d[i] = 0;
	}

	vec_free(T);
	vec_free(D);
}

DOUBLE MEMBER::GetLength() {
	if(is_curved) {
		DOUBLE L = 0;
		RPoint p1,p;
		POSITION pos = curved_list.GetHeadPosition();
        p1 = curved_list.GetNext(pos);
		while(pos) {
			p = curved_list.GetNext(pos);
			L += distance(p1,p);
			p1 = p;
		}
		return L;
	} else {
		return distance(j1->p,j2->p);
	}
}

int MEMBER::FindCurvedPoint(DOUBLE x,DOUBLE* len) {
	int i,count;
    DOUBLE L = 0;
	RPoint p1,p;
	POSITION pos;

	count = curved_list.GetCount();
	pos = curved_list.GetHeadPosition();
	p1 = curved_list.GetNext(pos);
	for(i = 1;i < count;i++) {
		p = curved_list.GetNext(pos);
		if(len) *len = L;
		L += distance(p1,p);
		if(EQUAL(L,x)) {
			if(!pos || (p1 = curved_list.GetNext(pos)) == p);
			else i++;
			break;
		} else if(L > x) {
			break;
		}
		p1 = p;
	}
	return (i - 1);
}

RPoint MEMBER::FindPoint(DOUBLE x) {
	RPoint rp;
	DOUBLE L,pL;
	if(!is_curved) {
		L = GetLength();
		rp =  j1->p + (j2->p - j1->p) * (x / L);
		return rp;
	}

	RPoint p1,p;
	POSITION pos;
	L = 0;
	pL = 0;
	pos = curved_list.GetHeadPosition();
	p1 = curved_list.GetNext(pos);
	while(pos) {
		p = curved_list.GetNext(pos);
		L += distance(p1,p);
		if(EQUAL(L,x)) {
			rp = p + j1->p;
			return rp;
		} else if(L > x) {
			rp = p1 + (p - p1) * ((x - pL) / (L - pL)) + j1->p;
			return rp;
		}
		p1 = p;
		pL = L;
	}
	return j2->p;
}

void MEMBER::SwapEnds() {
	POSITION pos;

	JOINT* temp = j1;
	j1 = j2;
	j2 = temp;
	
	int trel = nrelease;
	nrelease = frelease;
	frelease = trel;
	
	if(e_section) {

		SECTION* tsec = section;
		section = e_section;
		e_section = tsec;
		
		DOUBLE tf = start_offset;
		start_offset = end_offset;
		end_offset = tf;
	}

	if(is_curved) {
		RPoint* rp;
		pos = curved_list.GetHeadPosition();
		while(pos) {
			rp = &curved_list.GetNext(pos);
			*rp = *rp + (j1->p - j2->p);
		}
	}


	LOAD* pload;
	DOUBLE L = GetLength();
	pos = load.GetHeadPosition();
	while(pos) {
		pload = &load.GetNext(pos);
		if(pload->type == CONCENTRATED) {
			pload->x = L - pload->x;
		} else if(pload->type == TRAPEZOIDAL) {
			pload->x = L - pload->x;
			pload->x1 = L - pload->x1;
			DOUBLE temp = pload->P;
			pload->P = pload->P1;
			pload->P1 = temp;
		}
	}
}

void MEMBER::ApplyReleases(VECTOR Aml) {
	MATRIX Ktemp;
	mat_alloc(Ktemp,12,12);
	equ(Ktemp,K_fem,12,12);
	int i;
	for(i = 0;i < 6;i++) {
		if(nrelease & (1 << i))  condense(Ktemp,Aml,12,i);
	}
	for(i = 0;i < 6;i++) {
		if(frelease & (1 << i))  condense(Ktemp,Aml,12,i + 6);
	}
	vec_free(Ktemp);
}

void MEMBER::CalculateLocalConcentratedFEM(VECTOR Aml,DOUBLE P,DOUBLE a,DOUBLE L,int dir) {
	DOUBLE b = L - a;
	clear(Aml,12);

	int i,j;
	MATRIX t,t1,S,Tjk,Tjl,TlkT,Fall,Tml;
	VECTOR v;
	
	vec_alloc(v,6);
	mat_alloc(t,6,6);
	mat_alloc(t1,6,6);
	mat_alloc(S,6,6);
	mat_alloc(Tjk,6,6);
	mat_alloc(Tjl,6,6);
	mat_alloc(TlkT,6,6);
	mat_alloc(Fall,6,6);
	mat_alloc(Tml,12,6);
	
	for(i = 0;i < 6;i++) {
		Tjk[i][i] = 1;
		Tjl[i][i] = 1;
		TlkT[i][i] = 1;
	}
	Tjk[IRZ][IUY] = L;
	Tjk[IRY][IUZ] = -L;
	Tjl[IRZ][IUY] = a;
	Tjl[IRY][IUZ] = -a;
	TlkT[IUY][IRZ] = b;
	TlkT[IUZ][IRY] = -b;
	
	v[dir] = P;
	
	/*Tmlk*/
	DetermineFlexibility(t,L,FALSE);
	invert(t,S,6);
	
	DetermineFlexibility(Fall,a,FALSE);
	
	multiply(S,TlkT,t,6,6,6);
	multiply(t,Fall,t1,6,6,6);
	multiply(t1,-1,t,6,6);
	for(i = 0;i < 6;i++) {
		for(j = 0;j < 6;j++) {
			Tml[i + 6][j] = t[i][j];
		}
	}
	
	/*Tmlj*/
	multiply(Tjk,t,t1,6,6,6);
	add(t1,Tjl,6,6);
	multiply(t1,-1,t,6,6);
	for(i = 0;i < 6;i++) {
		for(j = 0;j < 6;j++) {
			Tml[i][j] = t[i][j];
		}
	}
	
	multiply(Tml,v,Aml,12,6);
	
	vec_free(v);
	vec_free(t);
	vec_free(t1);
	vec_free(S);
	vec_free(Tjk);
	vec_free(Tjl);
	vec_free(TlkT);
	vec_free(Fall);
	vec_free(Tml);

	ApplyReleases(Aml);
}

void MEMBER::CalculateLocalFEM(LOAD* pload,VECTOR Aml) {
    DOUBLE L = GetLength();
	clear(Aml,12);

	VECTOR temp;
	vec_alloc(temp,12);

    if(pload->type == CONCENTRATED) {
		CalculateLocalConcentratedFEM(Aml,pload->P,pload->x,L,pload->dir);
	} else {
		DOUBLE delta,x,p,P,P1,X,X1,slope;

		if(pload->type == UNIFORM) {
			P = pload->P;
			P1 = P;
			X = 0;
			X1 = L;
			slope = 0;
		} else {
			P = pload->P;
			P1 = pload->P1;
			X = pload->x;
			X1 = pload->x1;
			slope = (P1 - P) / (X1 - X);
		}

		/*simpson integraion*/
		const int DIV = 15 * 4 - 1;
		CalculateLocalConcentratedFEM(temp,P,X,L,pload->dir);
		add(Aml,temp,12);
		CalculateLocalConcentratedFEM(temp,P1,X1,L,pload->dir);
		add(Aml,temp,12);

		delta = (X1 - X) / (DIV - 1);
		for(int i = (0 + 1);i < (DIV - 1);i++) {
			x = X + i * delta;
			p = P + slope * (x - X);
			CalculateLocalConcentratedFEM(temp,p,x,L,pload->dir);
			if((i % 2) == 0) multiply(temp,2,temp,12);
			else multiply(temp,4,temp,12);
			add(Aml,temp,12);
		}
		multiply(Aml,delta / 3,Aml,12);
		/*end*/
	} 

	vec_free(temp);
}
void MEMBER::AssignSections() {
	DOUBLE L = GetLength();
	UBMP32 i,j;
	LOAD* pload;
	POSITION pos;
	DOUBLE v;

	/*sections at points of loading change*/
	UINT count = 0;
	pos = load.GetHeadPosition();
	while(pos) {
		pload = &load.GetNext(pos);
        if(pload->type == CONCENTRATED) {
			if(!EQUAL(pload->x,0) && !EQUAL(pload->x,L)) {
				station[count++] = pload->x - 1e-5;
				station[count++] = pload->x + 1e-5; 
			}
		} else if(pload->type == TRAPEZOIDAL) {
			if(!EQUAL(pload->x,0) && !EQUAL(pload->x,L)) {
				station[count++] = pload->x - 1e-5;
				station[count++] = pload->x + 1e-5; 
			}
			if(!EQUAL(pload->x1,0) && !EQUAL(pload->x1,L)) {
				station[count++] = pload->x1 - 1e-5;
				station[count++] = pload->x1 + 1e-5; 
			}
		}
	}

	UINT total = nDiv - count,start = count;
	for(i = 0; i < total;i++) {
		station[count++] = i * L / (total - 1);
	}

	station[count - 1] -= 1e-5;
    station[start] += 1e-5;

	for(i = 0; i < nDiv;i++) {
		for(j = i; j < nDiv;j++) {
            if(station[j] < station[i]) {
				v = station[i];
				station[i] = station[j];
				station[j] = v;
			}
		}
	}
}
/*
Internal force calculation
*/
void MEMBER::CalculateLocalConcentratedIF(VECTOR F,DOUBLE P,DOUBLE a,DOUBLE x,int dir) {
	clear(F,6);

	switch(dir) {
	case IUX:
		if(x >= a) F[IUX] = -P;
		break;
	case IUY:
		if(x >= a) {
			F[IUY] = -P;
			F[IRZ] = P * (x - a);
		}
		break;
	case IUZ:
		if(x >= a) {
			F[IUZ] = -P;
			F[IRY] = -P * (x - a);
		}
		break;
	case IRX:
		if(x >= a) F[IRX] = -P;
		break;
	case IRY:
		if(x >= a) F[IRY] = -P;
		break;
	case IRZ:
		if(x >= a) F[IRZ] = -P;
		break;
	};
}
void MEMBER::CalculateLocalIF(LOAD* pload,DOUBLE x, VECTOR F) {
	DOUBLE L = GetLength();
	clear(F,6);

    if(pload->type == CONCENTRATED) {
		CalculateLocalConcentratedIF(F,pload->P,pload->x,x,pload->dir);
	} else {
		DOUBLE P,P1,X,X1,slope;

		if(pload->type == UNIFORM) {
			P = pload->P;
			P1 = P;
			X = 0;
			X1 = L;
			slope = 0;
		} else {
			P = pload->P;
			P1 = pload->P1;
			X = pload->x;
			X1 = pload->x1;
			slope = (P1 - P) / (X1 - X);
		}
		
		if(X < x) {
			X1 = min(X1,x);
			P1 = P + slope * (X1 - X);
			CalculateLocalConcentratedIF(F,((P + P1) / 2) * (X1 - X),(X + X1) / 2,x,pload->dir);
		}
	} 
}

void MEMBER::CalculateInternalForceAndDisplacement(VECTOR q,VECTOR d,BOOL add_to_prev) {

	/*variables*/
	MATRIX F = forces;
	MATRIX D = disps;
	UINT i,j;
	DOUBLE L,x,delta;
	LOAD* pload;
	POSITION pos;

	/*non-linear analysis*/
	if(add_to_prev) {
		for(i = IUX;i <= IRZ;i++) {
			F[0][i] += -q[i];
            F[nDiv - 1][i] += q[i + 6];
			D[0][i] += d[i];
			D[nDiv - 1][i] += d[i + 6];
		}
		return;
	}
   	/*load rotation matrices*/
	VECTOR Fml;
	MATRIX T;
	MATRIX Tt;

	vec_alloc(Fml,6);
	mat_alloc(T,12,12);
	mat_alloc(Tt,12,12);
	
	L = GetLength();
	delta = L / (nMinDiv - 1);

	for(i = 0; i < nDiv;i++) {
		x = station[i];
		
		/*init*/
		F[i][IUX] = -q[IUX];
		F[i][IUY] = -q[IUY];
		F[i][IUZ] = -q[IUZ];
		F[i][IRX] = -q[IRX];
		F[i][IRY] = -q[IRY] - q[IUZ] * x;
		F[i][IRZ] = -q[IRZ] + q[IUY] * x;

	    /*local member forces*/
		pos = load.GetHeadPosition();
		while(pos) {
			pload = &load.GetNext(pos);

			/*joint load*/
		    if(!pload->is_member_load())
			    continue;

			/*apply factor*/
			if(!pload->is_considered())
				continue;

			pload->ApplyFactor();
		
			if(pload->system) {
				MATRIX T;
                VECTOR Fm,Fm11,Fm12;
				vec_alloc(Fm,12);
				vec_alloc(Fm11,12);
				vec_alloc(Fm12,12);
				mat_alloc(T,12,12);

			    GetLoadTransMatrix(T,pload->system);

				Fm[pload->dir] = pload->P;
				multiply(T,Fm,Fm11,12,12);
				Fm[pload->dir] = pload->P1;
				multiply(T,Fm,Fm12,12,12);

				for(int j = 0;j < 12;j++) {
					if(Fm11[j]) {
						LOAD mload = *pload;
						mload.P = Fm11[j];
                        mload.P1 = Fm12[j];
						mload.dir = IUX + (j % 6);
						mload.system = 0;
						CalculateLocalIF(&mload,x,Fml);
						add(F[i],Fml,6);
					}
				}

				vec_free(Fm12);
				vec_free(Fm11);
				vec_free(Fm); 
				vec_free(T);
			} else {
				CalculateLocalIF(pload,x,Fml);
				add(F[i],Fml,6);
			}
			/*remove factor*/
			pload->RemoveFactor();
		}
	}

	/*adjust forces*/
	for(j = IUX;j <= IRZ;j++) {
		delta = F[nDiv - 1][j] - q[j + 6];
		for(i = 0; i < nDiv;i++) {
			F[i][j] -= (station[i] / L) * delta;
			if(VEQUAL(F[i][j],0)) F[i][j] = 0;
		}
	}

    /*end section properties*/
    DOUBLE E1,E2,G1,G2,Ix1,Ix2,Iy1,Iy2,Iz1,Iz2;
	DOUBLE EI[2][3];
	SECTION* me_section = e_section ? e_section : section;
	
	E1 = section->material->E;
	G1 = section->material->G;
	Ix1 = section->fIx * section->Ix;
	Iy1 = section->fIy * section->Iy;
	Iz1 = section->fIz * section->Iz;
	
    E2 = me_section->material->E;
    G2 = me_section->material->G;
	Ix2 = me_section->fIx * me_section->Ix;
	Iy2 = me_section->fIy * me_section->Iy;
	Iz2 = me_section->fIz * me_section->Iz;
	
	EI[0][IUX] = G1 * Ix1 + (G2 * Ix2 - G1 * Ix1) * start_offset;
	EI[1][IUX] = G2 * Ix2 + (G1 * Ix1 - G2 * Ix2) * end_offset;
	if(e_section) {
		if(EIyy == LINEAR) {
			EI[0][IUY] = E1 * Iy1 + (E2 * Iy2 - E1 * Iy1) * start_offset;
			EI[1][IUY] = E2 * Iy2 + (E1 * Iy1 - E2 * Iy2) * end_offset;
		} else if(EIyy == PARABOLIC) { 
			EI[0][IUY] = pow(pow(E1 * Iy1,0.5) + (pow(E2 * Iy2,0.5) - pow(E1 * Iy1,0.5)) * start_offset, 2);
			EI[1][IUY] = pow(pow(E2 * Iy2,0.5) + (pow(E1 * Iy1,0.5) - pow(E2 * Iy2,0.5)) * end_offset, 2);
		} else if(EIyy == CUBIC) { 
			EI[0][IUY] = pow(pow(E1 * Iy1,1 / 3.0) + (pow(E2 * Iy2,1 / 3.0) - pow(E1 * Iy1,1 / 3.0)) * start_offset, 3);
			EI[1][IUY] = pow(pow(E2 * Iy2,1 / 3.0) + (pow(E1 * Iy1,1 / 3.0) - pow(E2 * Iy2,1 / 3.0)) * end_offset, 3);
		}
		
		if(EIzz == LINEAR) {
			EI[0][IUZ] = E1 * Iz1 + (E2 * Iz2 - E1 * Iz1) * start_offset;
			EI[1][IUZ] = E2 * Iz2 + (E1 * Iz1 - E2 * Iz2) * end_offset;
		} else if(EIzz == PARABOLIC) { 
			EI[0][IUZ] = pow(pow(E1 * Iz1,0.5) + (pow(E2 * Iz2,0.5) - pow(E1 * Iz1,0.5)) * start_offset, 2);
			EI[1][IUZ] = pow(pow(E2 * Iz2,0.5) + (pow(E1 * Iz1,0.5) - pow(E2 * Iz2,0.5)) * end_offset, 2);
		} else if(EIzz == CUBIC) { 
			EI[0][IUZ] = pow(pow(E1 * Iz1,1 / 3.0) + (pow(E2 * Iz2,1 / 3.0) - pow(E1 * Iz1,1 / 3.0)) * start_offset, 3);
			EI[1][IUZ] = pow(pow(E2 * Iz2,1 / 3.0) + (pow(E1 * Iz1,1 / 3.0) - pow(E2 * Iz2,1 / 3.0)) * end_offset, 3);
		}
	} else {
		EI[0][IUY] = E1 * Iy1 + (E2 * Iy2 - E1 * Iy1) * start_offset;
		EI[1][IUY] = E2 * Iy2 + (E1 * Iy1 - E2 * Iy2) * end_offset;
		EI[0][IUZ] = E1 * Iz1 + (E2 * Iz2 - E1 * Iz1) * start_offset;
		EI[1][IUZ] = E2 * Iz2 + (E1 * Iz1 - E2 * Iz2) * end_offset;
	}

	/*rotations*/
    DOUBLE sum,EIaverage,Faverage,x1,x2,y1,y2;
	UINT dir;
	for(j = IRX;j <= IRZ;j++) {
		if(j == IRZ) dir = IUY;
		else if(j == IRY) dir = IUZ;
		else dir = IUX;

		sum = 0;
		for(i = 1; i < nDiv;i++) {
			x1 = station[i - 1];
			x2 = station[i];
			y1 = EI[0][j - 3];
			y2 = EI[1][j - 3];
			EIaverage = ((y1 + y2) + (y1 - y2) * (x1 + x2)) / 2;
			Faverage = F[i - 1][j] + (3 * F[i - 1][dir] + F[i][dir]) * (x2 - x1) / 8;
			Faverage = (F[i - 1][j] + 4 * Faverage + F[i][j]) / 6;
			sum += Faverage * (x2 - x1) / EIaverage;
			
			D[i][j] = sum;
		}
	}

	/*adjust rotations*/
	for(j = IRX;j <= IRZ;j++) {
		delta = D[nDiv - 1][j];
		for(i = 0; i < nDiv;i++) {
			D[i][j] -= (station[i] / L) * delta;
			if(VEQUAL(D[i][j],0)) D[i][j] = 0;
		}
	}
	for(i = IRX;i <= IRZ;i++) {
		D[0][i] = d[i];
		D[nDiv - 1][i] = d[i + 6];
	}
	for(i = 1; i < nDiv - 1;i++) {
		for(j = IRX;j <= IRZ;j++) {
			D[i][j] += D[0][j] + (D[nDiv - 1][j] - D[0][j]) * station[i] / L;
		}
	}

	/*displacements*/
	DOUBLE Daverage;
	for(j = IUX;j <= IUZ;j++) {
		if(j == IUY) dir = IRZ;
		else if(j == IUZ) dir = IRY;
		else dir = IRX;
		sum = 0;
		for(i = 1; i < nDiv;i++) {
			x1 = station[i - 1];
			x2 = station[i];
		
			y1 = EI[0][dir - 3];
			y2 = EI[1][dir - 3];
			EIaverage = ((y1 + y2) + (y1 - y2) * (x1 + x2)) / 2;
			Daverage = D[i - 1][dir] + (3 * F[i - 1][dir] + F[i][dir]) * (x2 - x1) / (8 * EIaverage);
			Daverage = (D[i - 1][dir] + 4 * Daverage + D[i][dir]) / 6;
            
			if(j == IUZ) sum += -Daverage * (x2 - x1);
			else sum += Daverage * (x2 - x1);
			D[i][j] = sum;
		}
	}

	/*adjust displacements*/
	for(j = IUX;j <= IUZ;j++) {
		delta = D[nDiv - 1][j];
		for(i = 0; i < nDiv;i++) {
			D[i][j] -= (station[i] / L) * delta;
			if(VEQUAL(D[i][j],0)) D[i][j] = 0;
		}
	}
	for(i = IUX;i <= IUZ;i++) {
		D[0][i] = d[i];
		D[nDiv - 1][i] = d[i + 6];
	}
	for(i = 1; i < nDiv - 1;i++) {
		for(j = IUX;j <= IUZ;j++) {
			D[i][j] += D[0][j] + (D[nDiv - 1][j] - D[0][j]) * station[i] / L;
		}
	}
	/*free*/
	vec_free(Fml);
	vec_free(T);
	vec_free(Tt);
}
/*
calcuate internal force / deflection at x
*/
DOUBLE MEMBER::_CalculateInternal(BOOL bforce,DOUBLE x, int dir1,int myindex) {
	MATRIX myforces;
	UINT i;

	if(myindex == -1) {
		if(bforce) myforces = forces;
		else myforces = disps;
	} else {
		if(bforce) myforces = forces_all[myindex];
		else myforces = disps_all[myindex];
	}

	RPoint v;
	for(i = 0; i < nDiv;i++) {
		if(station[i] >= x) {
			break;
		}
	}
	if(i == 0 || EQUAL(x,station[i])) {
		return myforces[i][dir1];
	} else if(i == nDiv) {
		return myforces[i - 1][dir1];
	} else {
		return myforces[i - 1][dir1] + (myforces[i][dir1] - myforces[i - 1][dir1]) * (x - station[i - 1]) / (station[i] - station[i - 1]);
	}
}
DOUBLE MEMBER::CalculateForce(DOUBLE x, int dir1,int myindex) {
	return _CalculateInternal(TRUE,x,dir1,myindex);
}
DOUBLE MEMBER::CalculateDeflection(DOUBLE x,int dir1,int myindex) {
	return _CalculateInternal(FALSE,x,dir1,myindex);
}
void MEMBER::CalculateDeflection(DOUBLE x,RPoint& v,int index) {
	v.x = CalculateDeflection(x,IUX,index);
    v.y = CalculateDeflection(x,IUY,index);
    v.z = CalculateDeflection(x,IUZ,index);
}
/*
Calculate fixed end moment
*/
void MEMBER::CalculateFEM(VECTOR FEM) {
	LOAD* pload;
	DOUBLE L,value;

	L = GetLength();

	/*stiffness matrix*/
	mat_alloc(K_fem,12,12);
	GetLocalStiffnessMatrix(K_fem,1);

	/*load rotation matrices*/
	register int i;
	VECTOR Aml;
	MATRIX T;
	MATRIX Tt;
	MATRIX Kl;

	vec_alloc(Aml,12);
	mat_alloc(T,12,12);
	mat_alloc(Tt,12,12);
	mat_alloc(Kl,12,12);

	equ(Kl,K_fem,12,12);

	POSITION pos = load.GetHeadPosition();
	while(pos) {
		pload = &load.GetNext(pos);
		
		/*apply factor*/
		if(!pload->is_considered())
			continue;
		
		/*temperature and fabrication errors*/
		if(pload->is_strain_load()) {
			if(pload->type == TEMPERATURE) 
				value = section->material->alphac * pload->P;
			else
				value = pload->P;

			int dir;
			dir = IUX;
			FEM[dir] += Kl[dir][dir] * L * value;
			FEM[dir + 6] += -Kl[dir + 6][dir + 6] * L * value;
		}
		/*member load*/
		if(!pload->is_member_load())
			continue;
	
		pload->ApplyFactor();

		/*transform global load in to local load*/
		if(pload->system) {
			MATRIX T;
			VECTOR Am,Am11,Am12,Am2;
			vec_alloc(Am,12);
			vec_alloc(Am11,12);
			vec_alloc(Am12,12);
			vec_alloc(Am2,12);
			mat_alloc(T,12,12);

			GetLoadTransMatrix(T,pload->system);
			
			Am[pload->dir] = pload->P;
			multiply(T,Am,Am11,12,12);
			Am[pload->dir] = pload->P1;
			multiply(T,Am,Am12,12,12);
			
			for(i = 0;i < 12;i++) {
				if(!EQUAL(Am11[i],0)) {
					LOAD mload = *pload;
					mload.P = Am11[i];
					mload.P1 = Am12[i];
					mload.dir = IUX + (i % 6);
					mload.system = 0;
				    CalculateLocalFEM(&mload,Aml);
					add(FEM,Aml,12);
				}
			}

			vec_free(Am2);
			vec_free(Am11);
			vec_free(Am12);
			vec_free(Am);
			vec_free(T);
		} else {
			VECTOR Am;
			vec_alloc(Am,12);
 		    CalculateLocalFEM(pload,Aml);
			add(FEM,Aml,12);
			vec_free(Am);
		}
		/*remove factor*/
		pload->RemoveFactor();
		/*end*/
	}

	/*free*/
	vec_free(K_fem);
	vec_free(Aml);
	vec_free(T);
	vec_free(Tt);
	vec_free(Kl);
}
void MEMBER::ApplyLoad(VECTOR Q) {
	VECTOR FEM,temp;
	MATRIX T,Tt;

    vec_alloc(FEM,12);
	vec_alloc(temp,12);
	mat_alloc(T,12,12);
	mat_alloc(Tt,12,12);

	CalculateFEM(temp);
	GetTransformationMatrix(T);
	transpose(T,Tt,12,12);
	multiply(Tt,temp,FEM,12,12);

	for(int i = IUX;i <= IRZ;i++) {
		if(j1->number[i] != INVALID) Q[j1->number[i]] -= FEM[i];
        if(j2->number[i] != INVALID) Q[j2->number[i]] -= FEM[i + 6];
	}

	vec_free(FEM);
	vec_free(temp);
	vec_free(T);
	vec_free(Tt);
}

void MEMBER::CalcMF(VECTOR Qg,VECTOR Dg,VECTOR q,VECTOR d) {
	
	RPoint p;
	MATRIX Kl;
	MATRIX T;
	VECTOR q0;
	VECTOR D;
	
	mat_alloc(Kl,12,12);
	mat_alloc(T,12,12);
	vec_alloc(D,12);
	vec_alloc(q0,12);

	for(int i = IUX;i <= IRZ;i++) {
		if(j1->number[i] != INVALID) D[i] = Dg[j1->number[i]];
		if(j2->number[i] != INVALID) D[i + 6] = Dg[j2->number[i]];
	}
	
	GetLocalStiffnessMatrix(Kl);
	GetTransformationMatrix(T);
	
	multiply(T,D,d,12,12);
	multiply(Kl,d,q,12,12);
	CalculateFEM(q0);
	add(q,q0,12);

	vec_free(q0);
	vec_free(D);
	vec_free(T);
	vec_free(Kl);
}
/*
Member for a 3D frame
*/
void MEMBER::GetRotMatrix(MATRIX T,int which,int start,int end) {
	RPoint vx,rp1,rp2;
	if(start != -1) {
		rp1 = curved_list.GetAt(curved_list.FindIndex(start));
        rp2 = curved_list.GetAt(curved_list.FindIndex(end));
	} else {
		rp1 = j1->p;
		rp2 = j2->p;
	}
	
	vx = rp2 - rp1;
	if(which) {
		MATRIX R;
		mat_alloc(R,3,3);
		if(which == 1) j1->GetRotMatrix(R);
		else j2->GetRotMatrix(R);
		vx = GetGlobal(vx,R);
		vec_free(R);
	}
	GetLineRotMatrix(vx,alpha,T);
}
void MEMBER::GetLoadTransMatrix(MATRIX T,SYSTEM* system) {
	MATRIX R;
	mat_alloc(R,3,3);
	GetRotMatrix(R,0);
	if(system) {
		MATRIX R1,R2;
		mat_alloc(R1,3,3);
		mat_alloc(R2,3,3);
		GetPointRotMatrix(system->rotation,R1);
		multiply(R,R1,R2,3,3,3);
		equ(R,R2,3,3);
		vec_free(R1);
		vec_free(R2);
	}
	
	clear(T,12,12);
	equsec(T,R,3,3,0,0,0,0);
	equsec(T,R,3,3,3,3,0,0);
	equsec(T,R,3,3,6,6,0,0);
	equsec(T,R,3,3,9,9,0,0);

	vec_free(R);
}
void MEMBER::GetTransformationMatrix(MATRIX T) {
	MATRIX R;
	mat_alloc(R,3,3);
	GetRotMatrix(R,1);
	equsec(T,R,3,3,0,0,0,0);
	equsec(T,R,3,3,3,3,0,0);
	GetRotMatrix(R,2);
	equsec(T,R,3,3,6,6,0,0);
	equsec(T,R,3,3,9,9,0,0);
	vec_free(R);
}

void MEMBER::DetermineFlexibility(MATRIX t,DOUBLE len,BOOL invert) {
	DOUBLE a,b,L,x1,x2,
		   A1,A2,Ay1,Ay2,Az1,Az2,
		   E1,E2,G1,G2,EA1,EA2,
		   GIx1,GIx2,EIy1,EIy2,EIz1,EIz2,
		   Ix1,Ix2,Iy1,Iy2,Iz1,Iz2;
	SECTION* me_section = e_section ? e_section : section;

	x2 = len;
	x1 = 0;

	L = GetLength();

	/*end section properties*/
	E1 = section->material->E;
	G1 = section->material->G;
	A1 = section->fA * section->A;
	Ay1 = section->fAy * section->Ay;
    Az1 = section->fAz * section->Az;
	Ix1 = section->fIx * section->Ix;
	Iy1 = section->fIy * section->Iy;
	Iz1 = section->fIz * section->Iz;

    E2 = me_section->material->E;
    G2 = me_section->material->G;
	A2 = me_section->fA * me_section->A;
	Ay2 = section->fAy * me_section->Ay;
    Az2 = section->fAz * me_section->Az;
	Ix2 = me_section->fIx * me_section->Ix;
	Iy2 = me_section->fIy * me_section->Iy;
	Iz2 = me_section->fIz * me_section->Iz;

	GIx1 = G1 * Ix1 + (G2 * Ix2 - G1 * Ix1) * start_offset;
	GIx2 = G2 * Ix2 + (G1 * Ix1 - G2 * Ix2) * end_offset;
	EA1  = E1 * A1  + (E2 * A2  - E1 * A1)  * start_offset;
	EA2  = E2 * A2  + (E1 * A1  - E2 * A2)  * end_offset;
	if(e_section) {
		if(EIyy == LINEAR) {
			EIy1 = E1 * Iy1 + (E2 * Iy2 - E1 * Iy1) * start_offset;
			EIy2 = E2 * Iy2 + (E1 * Iy1 - E2 * Iy2) * end_offset;
		} else if(EIyy == PARABOLIC) { 
			EIy1 = pow(pow(E1 * Iy1,0.5) + (pow(E2 * Iy2,0.5) - pow(E1 * Iy1,0.5)) * start_offset, 2);
			EIy2 = pow(pow(E2 * Iy2,0.5) + (pow(E1 * Iy1,0.5) - pow(E2 * Iy2,0.5)) * end_offset, 2);
		} else if(EIyy == CUBIC) { 
			EIy1 = pow(pow(E1 * Iy1,1 / 3.0) + (pow(E2 * Iy2,1 / 3.0) - pow(E1 * Iy1,1 / 3.0)) * start_offset, 3);
			EIy2 = pow(pow(E2 * Iy2,1 / 3.0) + (pow(E1 * Iy1,1 / 3.0) - pow(E2 * Iy2,1 / 3.0)) * end_offset, 3);
		}

		if(EIzz == LINEAR) {
			EIz1 = E1 * Iz1 + (E2 * Iz2 - E1 * Iz1) * start_offset;
			EIz2 = E2 * Iz2 + (E1 * Iz1 - E2 * Iz2) * end_offset;
		} else if(EIzz == PARABOLIC) { 
			EIz1 = pow(pow(E1 * Iz1,0.5) + (pow(E2 * Iz2,0.5) - pow(E1 * Iz1,0.5)) * start_offset, 2);
			EIz2 = pow(pow(E2 * Iz2,0.5) + (pow(E1 * Iz1,0.5) - pow(E2 * Iz2,0.5)) * end_offset, 2);
		} else if(EIzz == CUBIC) { 
			EIz1 = pow(pow(E1 * Iz1,1 / 3.0) + (pow(E2 * Iz2,1 / 3.0) - pow(E1 * Iz1,1 / 3.0)) * start_offset, 3);
			EIz2 = pow(pow(E2 * Iz2,1 / 3.0) + (pow(E1 * Iz1,1 / 3.0) - pow(E2 * Iz2,1 / 3.0)) * end_offset, 3);
		}
	} else {
		EIy1 = E1 * Iy1 + (E2 * Iy2 - E1 * Iy1) * start_offset;
		EIy2 = E2 * Iy2 + (E1 * Iy1 - E2 * Iy2) * end_offset;
		EIz1 = E1 * Iz1 + (E2 * Iz2 - E1 * Iz1) * start_offset;
		EIz2 = E2 * Iz2 + (E1 * Iz1 - E2 * Iz2) * end_offset;
	}
	/*end*/

	if(!e_section || EQUAL(EIy1 , EIy2)) {
		DOUBLE q = (12 * EIy1 * len) / (G1 * Az1) + 4 * pow(len,3);
		t[IUZ][IUZ] = q / (12 * EIy1);
		t[IUZ][IRY] = t[IRY][IUZ] = pow(len,2) / (2 * EIy1);
		t[IRY][IRY] = len / (EIy1);
	} else {
		if(EIyy == LINEAR) {
			if(!invert) {
				a = EIy1 + (EIy2 - EIy1) * len / L;
				b = (EIy1 - EIy2) / L;
			} else {
				a = EIy2 + (EIy1 - EIy2) * len / L;
				b = (EIy2 - EIy1) / L;
			}
			t[IUZ][IUZ] = (x2 * x2 / (2 * b) - x2 * a / (b * b) + a * a * log(a + b * x2) / pow(b,3)) - 
				(x1 * x1 / (2 * b) - x1 * a / (b * b) + a * a * log(a + b * x1) / pow(b,3));
			t[IUZ][IRY] = t[IRY][IUZ] = (x2 / b - a * log(a + b * x2) / (b * b)) -
				(x1 / b - a * log(a + b * x1) / (b * b));
			t[IRY][IRY] = log(a + b * x2) / b - 
				log(a + b * x1) / b;
		} else if(EIyy == PARABOLIC) {
			if(!invert) {
				a = pow(EIy1,0.5) + (pow(EIy2,0.5) - pow(EIy1,0.5)) * len / L;
				b = (pow(EIy1,0.5) - pow(EIy2,0.5)) / L;
			} else {
				a = pow(EIy2,0.5) + (pow(EIy1,0.5) - pow(EIy2,0.5)) * len / L;
				b = (pow(EIy2,0.5) - pow(EIy1,0.5)) / L;
			}
			t[IUZ][IUZ] = (x2 / (b * b) - (a * a) / (pow(b,3) * (a + b * x2)) - 2 * a * log(a + b * x2) / pow(b,3)) -
				(x1 / (b * b) - (a * a) / (pow(b,3) * (a + b * x1)) - 2 * a * log(a + b * x1) / pow(b,3));
			t[IUZ][IRY] = t[IRY][IUZ] = (a / ((a + b * x2) * b * b) + log(a + b * x2)/ (b * b)) - 
				(a / ((a + b * x1) * b * b) + log(a + b * x1)/ (b * b));
			t[IRY][IRY] = -1 / ((a + b * x2) * b) -
				-1 / ((a + b * x1) * b);
		} else {
			if(!invert) {
				a = pow(EIy1,1/3.0) + (pow(EIy2,1/3.0) - pow(EIy1,1/3.0)) * len / L;
				b = (pow(EIy1,1/3.0) - pow(EIy2,1/3.0)) / L;
			} else {
				a = pow(EIy2,1/3.0) + (pow(EIy1,1/3.0) - pow(EIy2,1/3.0)) * len / L;
				b = (pow(EIy2,1/3.0) - pow(EIy1,1/3.0)) / L;
			}
			
			t[IUZ][IUZ] = (2 * a / (pow(b,3) * (a + b * x2)) - a * a / (2 * pow(b,3) * pow(a + b * x2,2)) + log(a + b * x2) / pow(b,3)) -
				(2 * a / (pow(b,3) * (a + b * x1)) - a * a / (2 * pow(b,3) * pow(a + b * x1,2)) + log(a + b * x1) / pow(b,3));
			t[IUZ][IRY] = t[IRY][IUZ] = (-1 / ((a + b * x2) * b * b) + a / (2 * pow(a + b * x2,2) * b * b)) -  
				(-1 / ((a + b * x1) * b * b) + a / (2 * pow(a + b * x1,2) * b * b));
			t[IRY][IRY] = -1 / (2 * pow(a + b * x2,2) * b) -
				-1 / (2 * pow(a + b * x1,2) * b);
		}
	}

	if(!e_section || EQUAL(EIz1 , EIz2)) {
		DOUBLE q = (12 * EIz1 * len) / (G1 * Ay1) + 4 * pow(len,3);
		t[IUY][IUY] = q / (12 * EIz1);
		t[IUY][IRZ] = t[IRZ][IUY] = pow(len,2) / (2 * EIz1);
		t[IRZ][IRZ] = len / (EIz1);
	} else {
		if(EIzz == LINEAR) {
			if(!invert) {
				a = EIz1 + (EIz2 - EIz1) * len / L;
				b = (EIz1 - EIz2) / L;
			} else {
				a = EIz2 + (EIz1 - EIz2) * len / L;
				b = (EIz2 - EIz1) / L;
			}
			t[IUY][IUY] = (x2 * x2 / (2 * b) - x2 * a / (b * b) + a * a * log(a + b * x2) / pow(b,3)) - 
				(x1 * x1 / (2 * b) - x1 * a / (b * b) + a * a * log(a + b * x1) / pow(b,3));
			t[IUY][IRZ] = t[IRZ][IUY] = (x2 / b - a * log(a + b * x2) / (b * b)) -
				(x1 / b - a * log(a + b * x1) / (b * b));
			t[IRZ][IRZ] = log(a + b * x2) / b - 
				log(a + b * x1) / b;
		} else if(EIzz == PARABOLIC) {
			if(!invert) {
				a = pow(EIz1,0.5) + (pow(EIz2,0.5) - pow(EIz1,0.5)) * len / L;
				b = (pow(EIz1,0.5) - pow(EIz2,0.5)) / L;
			} else {
				a = pow(EIz2,0.5) + (pow(EIz1,0.5) - pow(EIz2,0.5)) * len / L;
				b = (pow(EIz2,0.5) - pow(EIz1,0.5)) / L;
			}
			
			t[IUY][IUY] = (x2 / (b * b) - (a * a) / (pow(b,3) * (a + b * x2)) - 2 * a * log(a + b * x2) / pow(b,3)) -
				(x1 / (b * b) - (a * a) / (pow(b,3) * (a + b * x1)) - 2 * a * log(a + b * x1) / pow(b,3));
			t[IUY][IRZ] = t[IRZ][IUY] = (a / ((a + b * x2) * b * b) + log(a + b * x2)/ (b * b)) - 
				(a / ((a + b * x1) * b * b) + log(a + b * x1)/ (b * b));
			t[IRZ][IRZ] = -1 / ((a + b * x2) * b) -
				-1 / ((a + b * x1) * b);
		} else {
			if(!invert) {
				a = pow(EIz1,1/3.0) + (pow(EIz2,1/3.0) - pow(EIz1,1/3.0)) * len / L;
				b = (pow(EIz1,1/3.0) - pow(EIz2,1/3.0)) / L;
			} else {
				a = pow(EIz2,1/3.0) + (pow(EIz1,1/3.0) - pow(EIz2,1/3.0)) * len / L;
				b = (pow(EIz2,1/3.0) - pow(EIz1,1/3.0)) / L;
			}
			
			t[IUY][IUY] = (2 * a / (pow(b,3) * (a + b * x2)) - a * a / (2 * pow(b,3) * pow(a + b * x2,2)) + log(a + b * x2) / pow(b,3)) -
				(2 * a / (pow(b,3) * (a + b * x1)) - a * a / (2 * pow(b,3) * pow(a + b * x1,2)) + log(a + b * x1) / pow(b,3));
			t[IUY][IRZ] = t[IRZ][IUY] = (-1 / ((a + b * x2) * b * b) + a / (2 * pow(a + b * x2,2) * b * b)) -  
				(-1 / ((a + b * x1) * b * b) + a / (2 * pow(a + b * x1,2) * b * b));
			
			t[IRZ][IRZ] = -1 / (2 * pow(a + b * x2,2) * b) -
				-1 / (2 * pow(a + b * x1,2) * b);
		}
	}
	if(invert) {
		t[IUY][IRZ] *= -1;
		t[IRZ][IUY] *= -1;
	} else {
		t[IUZ][IRY] *= -1;
		t[IRY][IUZ] *= -1;
	}

	if(EQUAL(EA1 , EA2)) {
		t[IUX][IUX] = len / (EA1);
	} else {
		if(!invert) {
			a = EA1 + (EA2 - EA1) * len / L;
			b = (EA1 - EA2) / L;
		} else {
			a = EA2 + (EA1 - EA2) * len / L;
			b = (EA2 - EA1) / L;
		}
		t[IUX][IUX] = log(a + b * x2) / b - 
		          log(a + b * x1) / b;
	}

	if(EQUAL(GIx1 , GIx2)) {
		t[IRX][IRX] = len / (GIx1);
	} else {
		if(invert) {
			a = GIx1 + (GIx2 - GIx1) * len / L;
			b = (GIx1 - GIx2) / L;
		} else {
			a = GIx2 + (GIx1 - GIx2) * len / L;
			b = (GIx2 - GIx1) / L;
		}
		t[IRX][IRX] = log(a + b * x2) / b - 
		          log(a + b * x1) / b;
	}
}
void MEMBER::ApplyReleases(MATRIX K) {
	int i;
	for(i = 0;i < 6;i++) {
		if(nrelease & (1 << i)) condense(K,12,i);
	}
	for(i = 0;i < 6;i++) {
		if(frelease & (1 << i)) condense(K,12,i + 6);
	}
}
void MEMBER::ApplyGeometricStiffness(MATRIX K) {
	DOUBLE N,L,a;

	L = GetLength();
	N = (forces[0][IUX] + forces[nDiv - 1][IUX]) / 2;

	if(MEMBER::p_delta == LINEARIZED) {
		a = N / L;
		
		K[IUY][IUY] += a;
		K[6 + IUY][6 + IUY] += a;
		K[IUY][6 + IUY] += -a;
		K[6 + IUY][IUY] += -a;
		
		K[IUZ][IUZ] += a;
		K[6 + IUZ][6 + IUZ] += a;
		K[IUZ][6 + IUZ] += -a;
		K[6 + IUZ][IUZ] += -a;
	} else {
		DOUBLE fac1,fac2,fac3,fac4,fac;
		fac = N / (30.0 * L);
		fac1 = 36 * fac;
		fac2 = 3 * L * fac;
		fac3 = 4 * L * L * fac;
		fac4 = -L * L * fac;
		
		MATRIX Kl;
		mat_alloc(Kl,12,12);
				
		Kl[1][1] = fac1;
		Kl[2][2] = fac1;
		Kl[4][4] = fac3;
		Kl[5][5] = fac3;
		Kl[7][7] = fac1;
		Kl[8][8] = fac1;
		Kl[10][10] = fac3;
		Kl[11][11] = fac3;
		
		Kl[2][4] = Kl[4][2] = -fac2;
		Kl[1][5] = Kl[5][1] = fac2;
		Kl[1][7] = Kl[7][1] = -fac1;
		Kl[2][8] = Kl[8][2] = -fac1;
		Kl[5][7] = Kl[7][5] = -fac2;
		Kl[4][8] = Kl[8][4] = fac2;
		Kl[2][10] = Kl[10][2] = -fac2;
		Kl[1][11] = Kl[11][1] = fac2;
		Kl[4][10] = Kl[10][4] = fac4;
		Kl[5][11] = Kl[11][5] = fac4;
		Kl[7][11] = Kl[11][7] = -fac2;
		Kl[8][10] = Kl[10][8] = fac2;
		
		add(K,Kl,12,12);
		vec_free(Kl);
	}
}
void MEMBER::GetLocalStiffnessMatrix(MATRIX K,int skip) {
	clear(K,12,12);

	if(!MEMBER::only_geometric_stiffness) {
		int i,j;
		
		MATRIX t,t1,skk,T,Tt;
		mat_alloc(T,6,6);
		mat_alloc(Tt,6,6);
		mat_alloc(t,6,6);
		mat_alloc(t1,6,6);
		mat_alloc(skk,6,6);
		
		
		DOUBLE L = GetLength();
		for(i = 0;i < 6;i++) T[i][i] = 1;
		T[IRZ][IUY] = L;
		T[IRY][IUZ] = -L;
		transpose(T,Tt,6,6);
		
		DetermineFlexibility(t,L,FALSE);
		invert(t,skk,6);
		
		for(i = 0;i < 6;i++) {
			for(j = 0;j < 6;j++) {
				K[6 + i][6 + j] = skk[i][j];
			}
		}
		
		multiply(T,skk,t,6,6,6);
		for(i = 0;i < 6;i++) {
			for(j = 0;j < 6;j++) {
				K[i][6 + j] = K[6 + j][i] = -t[i][j];
			}
		}
		
		
		multiply(T,skk,t1,6,6,6);
		multiply(t1,Tt,t,6,6,6);
		for(i = 0;i < 6;i++) {
			for(j = 0;j < 6;j++) {
				K[i][j] = t[i][j];
			}
		}
		
		vec_free(Tt);
		vec_free(T);
		vec_free(t);
		vec_free(t1);
		vec_free(skk);
    } 

	if(MEMBER::add_geometric_stiffness)
		ApplyGeometricStiffness(K);

	if(!skip)
	    ApplyReleases(K);
}
/*
Assemble joint stiffness matrix
*/
void MEMBER::Assemble(MATRAN K,MATRIX Kl,int NB) {
	JOINT* jt[2];
	jt[0] = j1;
	jt[1] = j2;
	ELEMENT::AssembleJoints(jt,2,K,Kl,NB);
}
void MEMBER::AddToK(MATRAN K,int NB) {
	MATRIX Kl;
	MATRIX Kg;
	MATRIX T;
	MATRIX Tt;
	MATRIX temp;

	mat_alloc(Kl,12,12);
    mat_alloc(Kg,12,12);
	mat_alloc(T,12,12);
	mat_alloc(Tt,12,12);
	mat_alloc(temp,12,12);

	GetLocalStiffnessMatrix(Kl);
	GetTransformationMatrix(T);
	transpose(T,Tt,12,12);

	multiply(Tt,Kl,temp,12,12,12);
	multiply(temp,T,Kg,12,12,12);

    Assemble(K,Kg,NB);
	 
	vec_free(temp);
	vec_free(Tt);
	vec_free(T);
	vec_free(Kg);
	vec_free(Kl);
}
/*
Join members
*/
void MEMBER::Join(MEMBERPLIST* list,PMEMBERLIST members) {
	DOUBLE d1,d2,s1,s2;
	LOAD* pload,*pload1,myload;
	BOOL found;
    POSITION pos1,pos2,pos3;
	MEMBER member,*pmember,*head,*tail;

	head = list->GetHead();
	tail = list->GetTail();
	member = *head;
	member.nrelease = head->nrelease;
    member.frelease = tail->frelease;
	member.section = head->section;
    member.e_section = tail->e_section;
	member.start_offset = head->start_offset;
	member.end_offset = tail->end_offset;

	pos1 = list->GetHeadPosition();
	list->GetNext(pos1);
	while(pos1) {
		pmember = list->GetNext(pos1);

		d1 = member.GetLength();
		d2 = d1 + pmember->GetLength();

        pos2 = pmember->load.GetHeadPosition();
		while(pos2) {
			pload = &pmember->load.GetNext(pos2);

			if(pload->type == CONCENTRATED) {
				pload->x += d1;

				/*add*/
				found = FALSE;
				pos3 = member.load.GetHeadPosition();
				while(pos3) {
					pload1 = &member.load.GetNext(pos3);
					if(pload1->loadcase == pload->loadcase 
						&& pload1->type == pload->type 
						&& pload1->system == pload->system
						&& pload1->dir == pload->dir
						&& EQUAL(pload1->x,pload->x)) {
						pload1->P += pload->P;
						found = TRUE;
						break;
					}
				}
				if(!found)
					member.load.AddTail(*pload);
				/*end*/
			} else if(pload->type == UNIFORM || pload->type == TRAPEZOIDAL) {
				if(pload->type == TRAPEZOIDAL) {
					pload->x += d1;
					pload->x1 += d1;
				} else {
					pload->type = TRAPEZOIDAL;
					pload->P1 = pload->P;
					pload->x = d1;
					pload->x1 = d2;
				}

				s1 = (pload->P1 - pload->P) / (pload->x1 - pload->x);
				/*add*/
				found = FALSE;
				pos3 = member.load.GetHeadPosition();
				while(pos3) {
					pload1 = &member.load.GetNext(pos3);
					if(pload1->loadcase == pload->loadcase 
						&& pload1->system == pload->system
						&& pload1->dir == pload->dir) {
						
						if(pload1->type == UNIFORM) s2 = 0;
						else s2 = (pload1->P1 - pload1->P) / (pload1->x1 - pload1->x);
						if(EQUAL(s1,s2)) {
							if(EQUAL(pload1->P1,pload->P)) {
								pload1->P1 = pload->P1;
								pload1->x1 = pload->x1;
								if(EQUAL(pload1->P,pload1->P1))
									pload1->type = UNIFORM;
								found = TRUE;
								break;
							}
							if(EQUAL(pload1->P,pload->P1)) {
								pload1->P = pload->P;
								pload1->x = pload->x;
								if(EQUAL(pload1->P,pload1->P1))
									pload1->type = UNIFORM;
								found = TRUE;
								break;
							}
						}
					}
				}
				if(!found)
					member.load.AddTail(*pload);
				/*end*/
			}
		}
		
		//will it be curved
		if(!member.is_curved) {
			if(pmember->is_curved || !OnSameLine(member.j1->p,member.j2->p,pmember->j2->p)) {
				member.is_curved = TRUE;
				member.curved_list.AddTail(member.j1->p - member.j1->p);
				member.curved_list.AddTail(member.j2->p - member.j1->p);
			}
		} 

		//add new member
		if(member.is_curved) {
			if(pmember->is_curved) {
				RPoint rp;
				pos3 = pmember->curved_list.GetHeadPosition();
				pmember->curved_list.GetNext(pos3);
				while(pos3) {
					rp = pmember->curved_list.GetNext(pos3) + pmember->j1->p;
					member.curved_list.AddTail(rp - member.j1->p);
				}
			} else {
				member.curved_list.AddTail(pmember->j2->p - member.j1->p);
			}
		}
		//connect
		member.j2 = pmember->j2;
	}
	member.j1->mconnect++;
    member.j2->mconnect++;
	members->AddTail(member);
}
/*
Break member
*/
void MEMBER::Break(JOINTPLIST* list,PMEMBERLIST members,MEMBERPLIST_PLIST* mesh_list) {
	
	POSITION pos1,pos2;
	LOAD* pload,myload;
	JOINT* myj1,*myj2;
	MEMBER member;
	DOUBLE L = GetLength(),len;
	DOUBLE d1,d2;

	member.CopyNormalData(this);

	MEMBERPLIST* mlist;
	if(mesh_list) {
		mlist = new MEMBERPLIST;
		mesh_list->AddTail(mlist);
	}

	len = 0;
	pos1 = list->GetHeadPosition();
	myj1 = list->GetNext(pos1);
	while(pos1) {
        myj2 = list->GetNext(pos1); 

		member.j1 = myj1;
		member.j2 = myj2;
		member.start_offset =  start_offset + (distance(member.j1->p,j1->p) / L) * (1 - start_offset - end_offset);
		member.end_offset =  end_offset + (distance(member.j2->p,j2->p) / L) * (1 - start_offset - end_offset);;
		if(member.j1->p == j1->p) member.nrelease = nrelease;
		else member.nrelease = 0;
		if(member.j2->p == j2->p) member.frelease = frelease;
		else member.frelease = 0;

		/*load*/
		member.load.RemoveAll();
		pos2 = load.GetHeadPosition();
		while(pos2) {
			pload = &load.GetNext(pos2);

			if(pload->type == CONCENTRATED) {
				d1 = len;
				d2 = len + distance(member.j2->p,member.j1->p);
				
				if((EQUAL(pload->x,d1) || pload->x > d1) 
					&& (EQUAL(pload->x,d2) || pload->x <= d2)) {
					myload = *pload;
					myload.x -= d1;
					member.load.AddTail(myload);
				}
			} else if(pload->type == UNIFORM) {
				member.load.AddTail(*pload);
			} else if(pload->type == TRAPEZOIDAL) {
				d1 = len;
				d2 = len + distance(member.j2->p,member.j1->p);
				
				if(pload->x < d2 && pload->x1 > d1) {
					myload = *pload;
					if(myload.x < d1) myload.x = d1;
					if(myload.x1 > d2) myload.x1 = d2;
					myload.P = pload->P + (pload->P1 - pload->P) * ((myload.x - pload->x) / (pload->x1 - pload->x));
                    myload.P1 = pload->P1 + (pload->P - pload->P1) * (1 - (myload.x1 - pload->x) / (pload->x1 - pload->x));
					myload.x -= d1;
					myload.x1 -= d1;
					member.load.AddTail(myload);
				}
			}
		}

		len += member.GetLength();

		member.name.Format("%d",++MEMBER::TotalNumber);
		member.j1->mconnect++;
        member.j2->mconnect++;
		members->AddTail(member);
		if(mesh_list) {
			MEMBER* pmem = &members->GetTail();
			mlist->AddTail(pmem);
		}
		myj1 = myj2;
	}
}
/*
Mesh member
*/
BOOL MEMBER::Mesh(PMEMBERLIST members,PJOINTLIST joints,MEMBERPLIST_PLIST* mesh_list) {
	JOINTPLIST list;
	JOINT *joint, *pjoint;
	POSITION pos1,pos2;

	DOUBLE L = GetLength();

	/*joint list*/
	list.AddTail(j1);
	
	/*Intermediate joints*/
	if(DivAtInterim) {
		pos1 = joints->GetHeadPosition();
		while(pos1) {
			joint = &joints->GetNext(pos1);
			if(j1->p == joint->p ||  j2->p == joint->p)
				continue;
			if(!(joint->mconnect + joint->sconnect))
				continue;
			if(ptinbetween(j1->p,j2->p,joint->p) && 
				OnSameLine(j1->p,j2->p,joint->p)) {
				list.AddTail(joint);
			}
		}
	}

	/*divisions*/
	JOINT jn;
	for(UINT i = 1;i < nMinFrameDiv;i++) {
		jn.p = j1->p + (j2->p - j1->p) * (i / DOUBLE(nMinFrameDiv));
		
		if(!DivAtInterim) {
			pos2 = joints->GetHeadPosition();
			while(pos2) {
				pjoint = &joints->GetNext(pos2);
				if(pjoint->p == jn.p) {
					list.AddTail(pjoint);
					goto SKIP;
				}
			}
		}
		
		pos1 = list.GetHeadPosition();
		while(pos1) {
			joint = list.GetNext(pos1);
			if(joint->p == jn.p) goto SKIP;
		}
		joints->AddTail(jn);
		joint = &joints->GetTail();
		joint->name.Format("%d",++JOINT::TotalNumber);
		list.AddTail(joint);
SKIP:;
	}

	list.AddTail(j2);

	/*sort*/
	DOUBLE d1,d2;
	JOINT* tj,**mj1,**mj2;
	pos1 = list.GetHeadPosition();
	while(pos1) {
		mj1 = &list.GetNext(pos1);
		
		pos2 = pos1;
		while(pos2) {
			mj2 = &list.GetNext(pos2);
			
			d1 = distance((*mj1)->p,j1->p);
			d2 = distance((*mj2)->p,j1->p); 
			if(d2 < d1) {
				tj = *mj1;
				*mj1 = *mj2;
				*mj2 = tj;
			}
		}
	}
	/*break member*/
	if(list.GetCount() > 2)	{
		Break(&list,members,mesh_list);
		return TRUE;
	}
	return FALSE;
}
/*
Mesh curved member
*/
BOOL MEMBER::MeshCurved(PMEMBERLIST members,PJOINTLIST joints,MEMBERPLIST_PLIST* mesh_list) {
	JOINTPLIST list;
	JOINT *joint;
	POSITION pos1,pos2;

	/*joint list*/
	list.AddTail(j1);
	
	/*curved member*/
	JOINT jn;
	pos1 = curved_list.GetHeadPosition();
	curved_list.GetNext(pos1);
	while(pos1) {
		jn.p = curved_list.GetNext(pos1) + j1->p;
		if(!pos1) break;
		
		BOOL found = FALSE;
		pos2 = joints->GetHeadPosition();
		while(pos2) {
			joint = &joints->GetNext(pos2);
			if(joint->p == jn.p) {
				found = TRUE;
				break;
			}
		}

        if(!found) {
			joints->AddTail(jn);
			joint = &joints->GetTail();
			joint->name.Format("%d",++JOINT::TotalNumber);
		}
		list.AddTail(joint);
	}
	
	list.AddTail(j2);
	
	/*break member*/
	if(list.GetCount() > 2)	{
		Break(&list,members,mesh_list);
		return TRUE;
	}
	return FALSE;
}
