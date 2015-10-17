#include "common.h"
/*
slab area
*/
DOUBLE SLAB::GetArea() {
	return quad_area(jt[0]->p,jt[1]->p,jt[2]->p,jt[3]->p);
}
/*
Mesh slab
*/
JOINT* AddJoint(RPoint& target,PJOINTLIST joints) {
	JOINT* joint;
	POSITION pos;
	pos = joints->GetHeadPosition();
	while(pos) {
		joint = &joints->GetNext(pos);
		if(joint->p == target)
			return joint;
	}

	JOINT jn;
	jn.p = target;
	jn.name.Format("%d",++JOINT::TotalNumber);
	joints->AddTail(jn);
	joint = &joints->GetTail();
	return joint;
}
void SLAB::MeshPoly(RPoint* rp,UINT np,UINT n,PJOINTLIST joints, PSLABLIST slabs, SLABPLIST* mlist) {
	UINT i,k;
	RPoint p0,midp[NPOLY];
	
	for(i = 0;i < np;i++) 
		p0 = p0 + rp[i];
	p0 = p0 * (1.0 / np);
	
	for(i = 0;i < np - 1;i++) 
		midp[i] = (rp[i] + rp[i + 1]) * 0.5;
	midp[i] = (rp[i] + rp[0]) * 0.5;
	
	if(n <= 0) {
		SLAB slab = *this;
		slab.nDivx = 1;
		slab.nDivy = 1;

		slab.NJ = section->type;
		if(section->type == ASECTION::PLATE4 || 
			section->type == ASECTION::SHELL4 || 
			section->type == ASECTION::Q4D) 
			slab.NJ = 4;

		for(i = 0;i < 4;i++) 
			slab.jt[i] = AddJoint(rp[i],joints);
		if(slab.NJ >= ASECTION::Q5) slab.jt[4] = AddJoint(midp[0],joints);
		if(slab.NJ >= ASECTION::Q6) slab.jt[5] = AddJoint(midp[1],joints);
		if(slab.NJ >= ASECTION::Q7) slab.jt[6] = AddJoint(midp[2],joints);
		if(slab.NJ >= ASECTION::Q8) slab.jt[7] = AddJoint(midp[3],joints);
		if(slab.NJ >= ASECTION::Q9) slab.jt[8] = AddJoint(p0,joints);
		
		for(k = 0;k < slab.NJ;k++)
			slab.jt[k]->sconnect++;
		slab.name.Format("%d",++SLAB::TotalNumber);
		slabs->AddTail(slab);
		
		if(mlist) {
			SLAB* pslab = &slabs->GetTail();
			mlist->AddTail(pslab);
		}
		return;
	}
	
	RPoint myrp[4];
	for(i = 0;i < np;i++) {
		myrp[0] = rp[i];
		myrp[1] = midp[i];
		myrp[2] = p0;
		myrp[3] = midp[(i == 0) ? np - 1 : i - 1];
		MeshPoly(myrp,4,n - 1,joints,slabs,mlist);
	}
}
void SLAB::MeshRect(PSLABLIST slabs,PJOINTLIST joints,SLABPLIST* mlist) {
    RPoint rp[4],sp[4];
	UINT i,j;
	
	for(i = 0;i < nDivy;i++) {
		sp[0] = jt[0]->p + (jt[3]->p - jt[0]->p) * (DOUBLE(i) / nDivy);
		sp[1] = jt[1]->p + (jt[2]->p - jt[1]->p) * (DOUBLE(i) / nDivy);
        sp[2] = jt[0]->p + (jt[3]->p - jt[0]->p) * (DOUBLE(i + 1) / nDivy);
		sp[3] = jt[1]->p + (jt[2]->p - jt[1]->p) * (DOUBLE(i + 1) / nDivy);
		for(j = 0;j < nDivx;j++) {
			rp[0] = sp[0] + (sp[1] - sp[0]) * (DOUBLE(j) / nDivx);
			rp[1] = sp[0] + (sp[1] - sp[0]) * (DOUBLE(j + 1) / nDivx);
			rp[2] = sp[2] + (sp[3] - sp[2]) * (DOUBLE(j + 1) / nDivx);
			rp[3] = sp[2] + (sp[3] - sp[2]) * (DOUBLE(j) / nDivx);
            			
			MeshPoly(rp,4,0,joints,slabs,mlist);
		}
	}
}
BOOL SLAB::Mesh(PSLABLIST slabs,PJOINTLIST joints,SLABPLIST_PLIST* mesh_list) {
	SLABPLIST* mlist;
	if(mesh_list) {
		mlist = new SLABPLIST;
		mesh_list->AddTail(mlist);
	}

	if(NJ == 4) {
		MeshRect(slabs,joints,mlist);
	} else {
		UINT i,depth = UINT(log(max(nDivx,nDivy)) / log(2));
		RPoint rp[NPOLY];
		for(i = 0;i < NJ;i++)
			rp[i] = jt[i]->p;
		MeshPoly(rp,NJ,depth,joints,slabs,mlist);
	}

	return TRUE;
}
/*
rotation matrix for slab
*/
void SLAB::GetRotMatrix(MATRIX R,int which) {
	UINT i;
	RPoint vx,vy,vz,p0,px,py;
	for(i = 0;i < NJ;i++) p0 = p0 + jt[i]->p;
	p0 = p0 * (1.0 / NJ);
	px = (jt[0]->p + jt[1]->p) * 0.5;
	py = (jt[1]->p + jt[2]->p) * 0.5;

	//vz
	vx = px - p0;
	vy = py - p0;
	vz = cross(vx,vy);
	if(which) {
		MATRIX R;
		mat_alloc(R,3,3);
		jt[which - 1]->GetRotMatrix(R);
		vz = GetGlobal(vz,R);
		vec_free(R);
	}

	GetPlaneRotMatrix(vz,alpha,R);
}
void SLAB::GetLoadTransMatrix(MATRIX T, SYSTEM* system) {
	UINT i,TOTAL = 6 * NJ;
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
	
	clear(T,TOTAL,TOTAL);
	for(i = 0;i < NJ;i++) {
		equsec(T,R,3,3,6 * i,6 * i,0,0);
		equsec(T,R,3,3,6 * i + 3,6 * i + 3,0,0);
	}

	vec_free(R);
}
void SLAB::GetTransformationMatrix(MATRIX T) {
	UINT i,TOTAL = 6 * NJ;
	MATRIX R;
	mat_alloc(R,3,3);
	
	clear(T,TOTAL,TOTAL);
	for(i = 0;i < NJ;i++) {
		GetRotMatrix(R,i + 1);
		equsec(T,R,3,3,6 * i,6 * i,0,0);
		equsec(T,R,3,3,6 * i + 3,6 * i + 3,0,0);
	}

	vec_free(R);
}
/*
Material property
*/
void SLAB::GetE(MATRIX E) {
	DOUBLE nu = section->material->nu;
	DOUBLE Em = section->material->E;
	DOUBLE h = section->h;
	DOUBLE factor;

	clear(E,3,3);

	if(section->state_type == ASECTION::PLANE_STRESS) {
		factor = Em / (1 - nu * nu);
		E[0][0] = 1 * factor;
		E[1][1] = 1 * factor;
		E[2][2] = (1 - nu) * factor / 2;
		E[0][1] = E[1][0] = nu * factor;
	} else if(section->state_type == ASECTION::PLANE_STRAIN) {
		factor = Em / ((1 + nu) * (1 - 2 * nu));
		E[0][0] = (1 - nu) * factor;
		E[1][1] = (1 - nu) * factor;
		E[2][2] = (1 - 2 * nu) * factor / 2;
		E[0][1] = E[1][0] = nu * factor;
	} else if(section->state_type == ASECTION::PLATE_ACTION) {
		factor = Em * pow(h,3) / (12 * (1 - nu * nu));
        E[0][0] = 1 * factor;
		E[1][1] = 1 * factor;
		E[2][2] = (1 - nu) * factor / 2;
		E[0][1] = E[1][0] = nu * factor;
	}
}
void SLAB::GetShapeFunction(VECTOR Nf,DOUBLE e,DOUBLE n,UINT ctype) {
	
	DOUBLE N[9];
	N[0] = (1 - e) * (1 - n) / 4;
	N[1] = (1 + e) * (1 - n) / 4;
	N[2] = (1 + e) * (1 + n) / 4;
	N[3] = (1 - e) * (1 + n) / 4;
	N[4] = (1 - e * e) * (1 - n) / 2;
	N[5] = (1 + e) * (1 - n * n) / 2;
	N[6] = (1 - e * e) * (1 + n) / 2;
	N[7] = (1 - e) * (1 - n * n) / 2;
	N[8] = (1 - e * e) * (1 - n * n);
	
	for(UINT i = 0;i < 9;i++) {
		Nf[i] = N[i];
	}
	
	if(ctype >= ASECTION::Q5) {
		Nf[0] += -N[4] / 2;
		Nf[1] += -N[4] / 2;
	}
	if(ctype >= ASECTION::Q6) {
		Nf[1] += -N[5] / 2;
		Nf[2] += -N[5] / 2;
	}
	if(ctype >= ASECTION::Q7) {
		Nf[2] += -N[6] / 2;
		Nf[3] += -N[6] / 2;
	}
	if(ctype >= ASECTION::Q8) {
		Nf[0] += -N[7] / 2;
		Nf[3] += -N[7] / 2;
	}
	if(ctype >= ASECTION::Q9) {
		Nf[0] += -N[8] / 4;
		Nf[1] += -N[8] / 4;
		Nf[2] += -N[8] / 4;
		Nf[3] += -N[8] / 4;
		Nf[4] += -N[8] / 2;
		Nf[5] += -N[8] / 2;
		Nf[6] += -N[8] / 2;
		Nf[7] += -N[8] / 2;
	}
}
void SLAB::GetShapeFunctionDerivative(VECTOR dN,bool et,DOUBLE e,DOUBLE n,UINT ctype) {
	DOUBLE N[9];
	if(et) {
		N[0] = -(1 - n) / 4;
		N[1] = (1 - n) / 4;
		N[2] = (1 + n) / 4;
		N[3] = -(1 + n) / 4;
		N[4] = (-2 * e) * (1 - n) / 2;
		N[5] = (1 - n * n) / 2;
		N[6] = (-2 * e) * (1 + n) / 2;
		N[7] = -(1 - n * n) / 2;
		N[8] = (-2 * e) * (1 - n * n);
	} else {
		N[0] = -(1 - e) / 4;
		N[1] = -(1 + e) / 4;
		N[2] = (1 + e) / 4;
		N[3] = (1 - e) / 4;
		N[4] = -(1 - e * e) / 2;
		N[5] = (1 + e) * (-2 * n) / 2;
		N[6] = (1 - e * e) / 2;
		N[7] = (1 - e) * (-2 * n) / 2;
		N[8] = (1 - e * e) * (-2 * n);
	}
	
	for(UINT i = 0;i < 9;i++) {
		dN[i] = N[i];
	}
	
	if(ctype >= ASECTION::Q5) {
		dN[0] += -N[4] / 2;
		dN[1] += -N[4] / 2;
	}
	if(ctype >= ASECTION::Q6) {
		dN[1] += -N[5] / 2;
		dN[2] += -N[5] / 2;
	}
	if(ctype >= ASECTION::Q7) {
		dN[2] += -N[6] / 2;
		dN[3] += -N[6] / 2;
	}
	if(ctype >= ASECTION::Q8) {
		dN[0] += -N[7] / 2;
		dN[3] += -N[7] / 2;
	}
	if(ctype >= ASECTION::Q9) {
		dN[0] += -N[8] / 4;
		dN[1] += -N[8] / 4;
		dN[2] += -N[8] / 4;
		dN[3] += -N[8] / 4;
		dN[4] += -N[8] / 2;
		dN[5] += -N[8] / 2;
		dN[6] += -N[8] / 2;
		dN[7] += -N[8] / 2;
	}
}
/*
Jacobian
*/
void SLAB::GetLocalJointCoordinates(RPoint* p,UINT ctype,RPoint* selrp) {
	UINT i;

	/*Get local joint coordinates*/
	RPoint p0;
	MATRIX R;
	mat_alloc(R,3,3);
	GetRotMatrix(R,0);

	for(i = 0;i < NJ;i++) {
        p[i] = GetGlobal(jt[i]->p,R);
		p0 = p0 + p[i];
	}
	p0 = p0 * (1.0 / NJ);

	/*DKT plates and shells*/
	if(ctype > NJ) {
		p[4] = (p[0] + p[1]) * 0.5;
		p[5] = (p[1] + p[2]) * 0.5;
		p[6] = (p[2] + p[3]) * 0.5;
		p[7] = (p[3] + p[0]) * 0.5;
	}

	/*shift*/
	for(i = 0;i < ctype ;i++) {
		p[i] = p[i] - p0;
	}

	/*selected point*/
	if(selrp) {
		*selrp = GetGlobal(*selrp,R) - p0;
	}

	vec_free(R);
}
void SLAB::__CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn,RPoint* p) {
	UINT i,ctype;
	ctype = section->type;
	if(ctype == ASECTION::PLATE4 ||
		ctype == ASECTION::SHELL4 ||
		ctype == ASECTION::Q4D)
		ctype = ASECTION::Q8;

	GetShapeFunctionDerivative(dNe,true,e,n,ctype);
	GetShapeFunctionDerivative(dNn,false,e,n,ctype);
	GetLocalJointCoordinates(p,ctype);

	/*jacobian*/
	clear(J,4);
	for(i = 0;i < ctype;i++) {
		J[0] += dNe[i] * p[i].x;
		J[1] += dNe[i] * p[i].y;
		J[2] += dNn[i] * p[i].x;
		J[3] += dNn[i] * p[i].y;
	}
	Jdet = J[0] * J[3] - J[2] * J[1];
}
void SLAB::CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn,RPoint* p) {
	__CalculateJacobian(J,Jdet,e,n,dNe,dNn,p);
}
void SLAB::CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn) {
	RPoint p[9];
	__CalculateJacobian(J,Jdet,e,n,dNe,dNn,p);
}
void SLAB::CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n) {
	DOUBLE dNe[9],dNn[9];
	RPoint p[9];
	__CalculateJacobian(J,Jdet,e,n,dNe,dNn,p);
}
/*
Quad DKT element
*/
void SLAB::GetDKTQStrainDisplacementMatrix(MATRIX BD,DOUBLE& Jdet,DOUBLE e,DOUBLE n) {
	DOUBLE dNe[9],dNn[9],J[4];
	RPoint p[9],p0;
	UINT i,j;

	CalculateJacobian(J,Jdet,e,n,dNe,dNn,p);

	/*variables*/
	DOUBLE xij,yij,lij;
	DOUBLE A[4],B[4],C[4],D[4],E[4];
	for(i = 0;i < 4;i++) {
		j = (i == 3) ? 0 : i + 1;
		p0 = p[i] - p[j];
		xij = p0.x;
		yij = p0.y;
		lij = p0.magnitude();
		A[i] = -pow(xij,1) / pow(lij,2);
		B[i] = (0.25 * pow(xij,2) - 0.5 * pow(yij,2)) / pow(lij,2);
		C[i] = -pow(yij,1) / pow(lij,2);
		D[i] = (0.25 * pow(yij,2) - 0.5 * pow(xij,2)) / pow(lij,2);
		E[i] = 0.75 * xij * yij / pow(lij,2);
	}

	/*Hxe,Hxn*/
	DOUBLE* dN;
	DOUBLE Hx[2][12],Hy[2][12];

	for(i = 0;i < 2;i++) {
		if(i == 0) dN = dNe;
		else dN = dNn;
		Hx[i][0] = 1.5 * (A[0] * dN[4] - A[3] * dN[7]);
		Hx[i][3] = 1.5 * (A[1] * dN[5] - A[0] * dN[4]);
		Hx[i][6] = 1.5 * (A[2] * dN[6] - A[1] * dN[5]);
		Hx[i][9] = 1.5 * (A[3] * dN[7] - A[2] * dN[6]);
		Hx[i][2] = (dN[0] - B[0] * dN[4] - B[3] * dN[7]);
		Hx[i][5] = (dN[1] - B[1] * dN[5] - B[0] * dN[4]);
		Hx[i][8] = (dN[2] - B[2] * dN[6] - B[1] * dN[5]);
        Hx[i][11] = (dN[3] - B[3] * dN[7] - B[2] * dN[6]);
		Hx[i][1] = (E[0] * dN[4] + E[3] * dN[7]);
		Hx[i][4] = (E[1] * dN[5] + E[0] * dN[4]);
		Hx[i][7] = (E[2] * dN[6] + E[1] * dN[5]);
		Hx[i][10] = (E[3] * dN[7] + E[2] * dN[6]);
	}

	for(i = 0;i < 2;i++) {
		if(i == 0) dN = dNe;
		else dN = dNn;
		Hy[i][0] = 1.5 * (C[0] * dN[4] - C[3] * dN[7]);
		Hy[i][3] = 1.5 * (C[1] * dN[5] - C[0] * dN[4]);
		Hy[i][6] = 1.5 * (C[2] * dN[6] - C[1] * dN[5]);
		Hy[i][9] = 1.5 * (C[3] * dN[7] - C[2] * dN[6]);
		Hy[i][1] = -(dN[0] - D[0] * dN[4] - D[3] * dN[7]);
		Hy[i][4] = -(dN[1] - D[1] * dN[5] - D[0] * dN[4]);
		Hy[i][7] = -(dN[2] - D[2] * dN[6] - D[1] * dN[5]);
		Hy[i][10] = -(dN[3] - D[3] * dN[7] - D[2] * dN[6]);
		Hy[i][2] = -(E[0] * dN[4] + E[3] * dN[7]);
		Hy[i][5] = -(E[1] * dN[5] + E[0] * dN[4]);
		Hy[i][8] = -(E[2] * dN[6] + E[1] * dN[5]);
        Hy[i][11] = -(E[3] * dN[7] + E[2] * dN[6]);
	}

	DOUBLE j11,j12,j21,j22;
	j11 = J[3] / Jdet;
	j12 = -J[1] / Jdet;
	j21 = -J[2] / Jdet;
	j22 = J[0] / Jdet;

	for(i = 0;i < 12;i++) {
		BD[0][i] = j11 * Hx[0][i] + j12 * Hx[1][i];
        BD[1][i] = j21 * Hy[0][i] + j22 * Hy[1][i]; 
		BD[2][i] = j11 * Hy[0][i] + j12 * Hy[1][i] + 
			       j21 * Hx[0][i] + j22 * Hx[1][i]; 
	}
}
void SLAB::GetQXStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n) {
	UINT i,r,s,SIZE = 2 * NJ;
	DOUBLE dNe[9],dNn[9],J[4],dNx[9],dNy[9];

	CalculateJacobian(J,Jdet,e,n,dNe,dNn);

	for(i = 0;i < NJ;i++) {
		dNx[i] = (J[3] * dNe[i] - J[1] * dNn[i]) / Jdet;
		dNy[i] = (-J[2] * dNe[i] + J[0] * dNn[i]) / Jdet;
	}

	for(s = 0;s < SIZE;s++) {
		r = (s % 2);
		i = (s / 2);
		if(r) {
			B[1][s] = B[2][s - 1] = dNy[i];
		} else {
			B[0][s] = B[2][s + 1] = dNx[i];
		}
	}
}
void SLAB::GetQDStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR b) {
	RPoint p[NPOLY];
	DOUBLE dNe[9],dNn[9],J[4],dNx[9],dNy[9],N[9];
	UINT i,j,k,l,m,r,s;
	DOUBLE lij,lik,cij,cik,sij,sik;
	RPoint pij,pik;

	/*displacement*/
    section->type = ASECTION::Q4;
	GetShapeFunction(N,e,n,ASECTION::Q4);
	CalculateJacobian(J,Jdet,e,n,dNe,dNn);
	for(i = 0;i < 4;i++) {
		dNx[i] = (J[3] * dNe[i] - J[1] * dNn[i]) / Jdet;
		dNy[i] = (-J[2] * dNe[i] + J[0] * dNn[i]) / Jdet;
	}
	/*Q8*/
	section->type = ASECTION::Q4D;
	CalculateJacobian(J,Jdet,e,n,dNe,dNn,p);
	for(i = 4;i < 8;i++) {
		dNx[i] = (J[3] * dNe[i] - J[1] * dNn[i]) / Jdet;
		dNy[i] = (-J[2] * dNe[i] + J[0] * dNn[i]) / Jdet;
	}

	/*fill matrix*/
	for(s = 0;s < 12;s++) {
		r = (s % 3);
		i = (s / 3);
		if(r == 0) {
			B[0][s] = B[2][s + 1] = dNx[i];
			if(b) {
				b[s] = -0.5 * dNy[i];
			}
		} else if(r == 1) {
			B[1][s] = B[2][s - 1] = dNy[i];
			if(b) {
				b[s] = 0.5 * dNx[i];
			}
		} else {
			j = ((i == 0) ? 3 : i - 1);
			k = ((i == 3) ? 0 : i + 1);
			m = (i + 4);
			l = ((m == 4) ? 7 : m - 1);
			
			pij = (p[j] - p[i]) * 1;
			pik = (p[k] - p[i]) * 1;
			lij = pij.magnitude();
			lik = pik.magnitude();
			cij = -pij.x / lij;
			cik = -pik.x / lik;
			sij = -pij.y / lij;
			sik = -pik.y / lik;
			
			B[0][s] = 0.125 * (lij * cij * dNx[l] - lik * cik * dNx[m]);
			B[1][s] = 0.125 * (lij * sij * dNy[l] - lik * sik * dNy[m]);
			B[2][s] = 0.125 * (lij * cij * dNy[l] - lik * cik * dNy[m] +
				lij * sij * dNx[l] - lik * sik * dNx[m] );
			if(b) {
				b[s] = -0.0625 * (lij * cij * dNy[l] - lik * cik * dNy[m]) +
					0.0625 * (lij * sij * dNx[l] - lik * sik * dNx[m]) -
					N[i];
			}
		}
	}
}
void SLAB::GetStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR b) {
	if(section->type == ASECTION::PLATE4) {
		GetDKTQStrainDisplacementMatrix(B,Jdet,e,n);
	} else if(section->type == ASECTION::Q4D) {
		GetQDStrainDisplacementMatrix(B,Jdet,e,n,b);
	} else {
		GetQXStrainDisplacementMatrix(B,Jdet,e,n);
	}
}
/*
Gauss integration rules
*/
#define MEMBRANE_SHELL4  ASECTION::Q4D

const UINT SLAB::Rule[] = {
	0,0,0,0, 2,3,3,3,3,3, 3, 2,3
};
static const int membrane_dir[] =  {
	IUX,IUY,IRZ
};
static const int plate_dir[] =  {
	IUZ,IRX,IRY
};
static const DOUBLE C1 = 1 / sqrt(3);
static const DOUBLE C2 = sqrt(3 / 5.0);
static const DOUBLE C3 = sqrt((3 - 2 * sqrt(6.0 / 5)) / 7.0);
static const DOUBLE C4 = sqrt((3 + 2 * sqrt(6.0 / 5)) / 7.0);
static const DOUBLE C5 = sqrt(5 + 2 * sqrt(10.0 / 7)) / 3;
static const DOUBLE C6 = sqrt(5 - 2 * sqrt(10.0 / 7)) / 3;

static const DOUBLE W1 = 5.0 / 9;
static const DOUBLE W2 = 8.0 / 9;
static const DOUBLE W3 = 1 / 2.0 - sqrt(5.0 / 6) / 6;
static const DOUBLE W4 = 1 / 2.0 + sqrt(5.0 / 6) / 6;
static const DOUBLE W5 = (322 - 13 * sqrt(70)) / 900;
static const DOUBLE W6 = (322 + 13 * sqrt(70)) / 900;
static const DOUBLE W7 = (512.0) / 900;

static DOUBLE gauss_points[5][5] = {
	  0,      0,      0,      0,      0,
    -C1,     C1,      0,      0,      0,
	-C2,      0,     C2,      0,      0,
	-C4,    -C3,     C3,     C4,      0,
	-C5,    -C6,      0,     C6,     C5
};
static DOUBLE gauss_weights[5][5] = {
	  2,      0,      0,      0,      0,
      1,      1,      0,      0,      0,
	 W1,     W2,     W1,      0,      0,
	 W3,     W4,     W4,     W3,      0,
	 W5,     W6,     W7,     W6,     W5
};

void SLAB::__GetLocalStiffnessMatrix(MATRIX Kt) {
	DOUBLE h = section->h;
	DOUBLE b[12],Jdet,e,n,ew,nw,c;
	MATRIX B,Bt,E,temp;
	MATRIX K;
	UINT SIZE;
	
	SIZE = 2 * NJ;
	if(section->type == ASECTION::PLATE4 || section->type == ASECTION::Q4D)
		SIZE = 3 * NJ;

	UINT RULE = Rule[section->type];

	clear(Kt,SIZE,SIZE);

	mat_alloc(K,SIZE,SIZE);
	mat_alloc(B,3,SIZE);
	mat_alloc(Bt,SIZE,3);
	mat_alloc(temp,SIZE,3);
	mat_alloc(E,3,3);
	
	GetE(E);

	for(UINT l = 0;l < RULE;l++) {
		e = gauss_points[RULE - 1][l];
		ew = gauss_weights[RULE - 1][l];

        for(UINT m = 0;m < RULE;m++) { 
			n = gauss_points[RULE - 1][m];
			nw = gauss_weights[RULE - 1][m];

			GetStrainDisplacementMatrix(B,Jdet,e,n,b);
			c = ew * nw * Jdet;
			if(section->type <= ASECTION::Q4D) 
				c *= h;
			
			transpose(B,Bt,3,SIZE);
			multiply(Bt,E,temp,SIZE,3,3);
			multiply(temp,B,K,SIZE,3,SIZE);
            multiply(K,c,SIZE,SIZE);
			add(Kt,K,SIZE,SIZE);

			/*Q4D penalty matrix*/
            if(section->type == ASECTION::Q4D) {
                multiply(b,b,K,SIZE,SIZE);
				multiply(K,c,SIZE,SIZE);
				multiply(K,section->material->G,SIZE,SIZE);
                add(Kt,K,SIZE,SIZE);
			}
		}
	}

	vec_free(B);
	vec_free(Bt);
	vec_free(E);
	vec_free(temp);
	vec_free(K);
}
void SLAB::GetLocalStiffnessMatrix(MATRIX Kt) {
	if(section->type == ASECTION::SHELL4) {
		
		UINT state_type = section->state_type;
		
		if(state_type == ASECTION::PLANE_STRESS 
			|| state_type == ASECTION::PLANE_STRAIN 
			|| state_type == ASECTION::SHELL_ACTION) {

			section->type = MEMBRANE_SHELL4;
			section->state_type = (state_type == ASECTION::PLANE_STRAIN) 
				? ASECTION::PLANE_STRAIN : ASECTION::PLANE_STRESS;
			GetLocalStiffnessMatrix(Kt);
		}

		if(state_type == ASECTION::PLATE_ACTION 
			|| state_type == ASECTION::SHELL_ACTION) {
			section->type = ASECTION::PLATE4;
			section->state_type = ASECTION::PLATE_ACTION;
			GetLocalStiffnessMatrix(Kt);
		}

		section->type = ASECTION::SHELL4;
		section->state_type = state_type;
		return;
	} 
	UINT i,j,ig,jg,SIZE;

	SIZE = 2 * NJ;
	if(section->type == ASECTION::PLATE4 || section->type == ASECTION::Q4D)
		SIZE = 3 * NJ;
	
	MATRIX Kl;
	mat_alloc(Kl,SIZE,SIZE);
	
	__GetLocalStiffnessMatrix(Kl);
	if(section->type == ASECTION::PLATE4) {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 3) * 6 + plate_dir[(i % 3)];
			for(j = 0;j < SIZE;j++) {
				jg = (j / 3) * 6  + plate_dir[(j % 3)];
				Kt[ig][jg] = Kl[i][j];
			}
		}
	} else if(section->type == ASECTION::Q4D) {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 3) * 6 + membrane_dir[(i % 3)];
			for(j = 0;j < SIZE;j++) {
				jg = (j / 3) * 6  + membrane_dir[(j % 3)];
				Kt[ig][jg] = Kl[i][j];
			}
		}
	} else {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 2) * 6 + membrane_dir[(i % 2)];
			for(j = 0;j < SIZE;j++) {
				jg = (j / 2) * 6  + membrane_dir[(j % 2)];
				Kt[ig][jg] = Kl[i][j];
			}
		}
	}
	
	vec_free(Kl);
}
/*
Assemble joint stiffness matrix
*/
void SLAB::Assemble(MATRAN K,MATRIX Kl,int NB) {
	ELEMENT::AssembleJoints(jt,NJ,K,Kl,NB);
}
void SLAB::AddToK(MATRAN K,int NB) {
	MATRIX Kl;
	MATRIX Kg;
	MATRIX T;
	MATRIX Tt;
	MATRIX temp;
	UINT TOTAL = 6 * NJ;
	
	mat_alloc(Kl,TOTAL,TOTAL);
    mat_alloc(Kg,TOTAL,TOTAL);
	mat_alloc(T,TOTAL,TOTAL);
	mat_alloc(Tt,TOTAL,TOTAL);
	mat_alloc(temp,TOTAL,TOTAL);

	GetLocalStiffnessMatrix(Kl);
	GetTransformationMatrix(T);
	transpose(T,Tt,TOTAL,TOTAL);
	multiply(Tt,Kl,temp,TOTAL,TOTAL,TOTAL);
	multiply(temp,T,Kg,TOTAL,TOTAL,TOTAL);

    Assemble(K,Kg,NB);
	 
	vec_free(temp);
	vec_free(Tt);
	vec_free(T);
	vec_free(Kg);
	vec_free(Kl);
}
static DOUBLE value_en[9][2] = {
	-1,-1,
     1,-1,
	 1, 1,
	-1, 1,
	 0,-1,
	 1, 0,
	 0, 1,
	-1, 0,
	 0, 0
};
void SLAB::__CalculateStresses(VECTOR d) {
	DOUBLE sigma[3],sigma0[3],e,n,Jdet,factor,s,h;
	MATRIX E,B,temp;

	UINT SIZE = 2 * NJ;
	if(section->type == ASECTION::PLATE4 || section->type == ASECTION::Q4D)
		SIZE = 3 * NJ;

	h = section->h;
	factor = 6.0 / (h * h);
	mat_alloc(E,3,3);
	mat_alloc(B,3,SIZE);
	mat_alloc(temp,3,SIZE);

	/*constitutive matrix*/
	GetE(E);

	/*internal stresses*/
	for(UINT j = 0;j < NJ;j++) {
		e = value_en[j][0];
		n = value_en[j][1];

		GetStrainDisplacementMatrix(B,Jdet,e,n);
		multiply(E,B,temp,3,3,SIZE);
		multiply(temp,d,sigma,3,SIZE);

		/*Initial stress/moment*/
	    CalculateInitialStress(sigma0,E);
		add(sigma,sigma0,3);

		/*Intrernal forces/moments*/
		if(section->type == ASECTION::PLATE4) {
			multiply(sigma,-1,3);
			forces[NSTRESS * j + IMX - ISX] += sigma[0];
			forces[NSTRESS * j + IMY - ISX] += sigma[1];
			forces[NSTRESS * j + IMXY - ISX] += sigma[2];
			s = (sigma[0] + sigma[1]) / 2 + sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
			forces[NSTRESS * j + IMMAX - ISX] += s;
			s = (sigma[0] + sigma[1]) / 2 - sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
			forces[NSTRESS * j + IMMIN - ISX] += s;
            multiply(sigma,factor,3);
		} else {
			multiply(sigma,h,3);
            forces[NSTRESS * j + IFX - ISX] += sigma[0];
			forces[NSTRESS * j + IFY - ISX] += sigma[1];
			forces[NSTRESS * j + IFXY - ISX] += sigma[2];
			s = (sigma[0] + sigma[1]) / 2 + sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
			forces[NSTRESS * j + IFMAX - ISX] += s;
			s = (sigma[0] + sigma[1]) / 2 - sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
			forces[NSTRESS * j + IFMIN - ISX] += s;
			s = sqrt(pow(sigma[0] - sigma[1],2) + pow(sigma[0],2) + pow(sigma[1],2) + 6 * pow(sigma[2],2)) / sqrt(2);
			forces[NSTRESS * j + IFVM - ISX] += s;
			multiply(sigma,1 / h,3);
		}

		forces[NSTRESS * j + ISX - ISX] += sigma[0];
		forces[NSTRESS * j + ISY - ISX] += sigma[1];
		forces[NSTRESS * j + ISXY - ISX] += sigma[2];
        s = (sigma[0] + sigma[1]) / 2 + sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
		forces[NSTRESS * j + ISMAX - ISX] += s;
		s = (sigma[0] + sigma[1]) / 2 - sqrt(pow((sigma[0] - sigma[1]) / 2,2) + pow(sigma[2],2));
		forces[NSTRESS * j + ISMIN - ISX] += s;
		s = sqrt(pow(sigma[0] - sigma[1],2) + pow(sigma[0],2) + pow(sigma[1],2) + 6 * pow(sigma[2],2)) / sqrt(2);
		forces[NSTRESS * j + ISVM - ISX] += s;
	}

	vec_free(E);
	vec_free(B);
	vec_free(temp);
}
void SLAB::CalculateStresses(VECTOR Qg,VECTOR Dg) {

	/*shell*/
	if(section->type == ASECTION::SHELL4) {
		UINT state_type = section->state_type;

		if(state_type == ASECTION::PLANE_STRESS 
			|| state_type == ASECTION::PLANE_STRAIN 
			|| state_type == ASECTION::SHELL_ACTION) {
			
			section->type = MEMBRANE_SHELL4;
			section->state_type = (state_type == ASECTION::PLANE_STRAIN) 
				? ASECTION::PLANE_STRAIN : ASECTION::PLANE_STRESS;
			CalculateStresses(Qg,Dg);
		}
		
		if(state_type == ASECTION::PLATE_ACTION 
			|| state_type == ASECTION::SHELL_ACTION) {
			section->type = ASECTION::PLATE4;
			section->state_type = ASECTION::PLATE_ACTION;
			CalculateStresses(Qg,Dg);
		}

		section->type = ASECTION::SHELL4;
		section->state_type = state_type;
		return;
	} 

	/*calculate stresses*/
	UINT i,j,TOTAL = 6 * NJ;
	VECTOR D,d;
	MATRIX T;
	
	mat_alloc(T,TOTAL,TOTAL);
	vec_alloc(D,TOTAL);
	vec_alloc(d,TOTAL);

	/*joint displacements*/
	for(i = 0;i < NJ;i++) {
		for(j = IUX;j <= IRZ;j++) {
			if(jt[i]->number[j] != INVALID) 
				disps[6 * i + j] = Dg[jt[i]->number[j]];
		}
	}

	/*to local area coordinate*/
	GetTransformationMatrix(T);
	multiply(T,disps,D,TOTAL,TOTAL);
	for(i = 0;i < NJ;i++) {
		if(section->type == ASECTION::PLATE4) {
			d[3 * i + 0] = D[6 * i + IUZ];
			d[3 * i + 1] = D[6 * i + IRX];
			d[3 * i + 2] = D[6 * i + IRY];
		} else if(section->type == ASECTION::Q4D) {
			d[3 * i + 0] = D[6 * i + IUX];
			d[3 * i + 1] = D[6 * i + IUY];
			d[3 * i + 2] = D[6 * i + IRZ];
		} else {
			d[2 * i + 0] = D[6 * i + IUX];
			d[2 * i + 1] = D[6 * i + IUY];
		}
	}

	/*calculate stress*/
	__CalculateStresses(d);

	/*free*/
	vec_free(D);
	vec_free(d);
	vec_free(T);
}
void SLAB::CalculateInitialStress(VECTOR sigma0,MATRIX E) {
    LOAD* pload;
	DOUBLE sigma[3],strain[3],value;
	
	clear(sigma0,3);

	POSITION pos = load.GetHeadPosition();
	while(pos) {
		pload = &load.GetNext(pos);
		
		if(!pload->is_considered())
			continue;

		if(!pload->is_strain_load())
			continue;

		pload->ApplyFactor();

		/*calculate strain loading*/
		clear(strain,3);
        if(section->type == ASECTION::PLATE4) {
			if(pload->type == TEMPERATURE) {
				if(pload->dir == IRZ) {
					value = section->material->alphac * pload->P;
					strain[0] = value; 
					strain[1] = value;	
				}
			} else {
				value = pload->P;
				if(pload->dir == IRX) strain[0] = value; 
				if(pload->dir == IRY) strain[1] = value; 
				if(pload->dir == IRZ) strain[2] = value; 
			}
		} else {
			if(pload->type == TEMPERATURE) {
                if(pload->dir == IUX) {
					value = section->material->alphac * pload->P;
					strain[0] = value; 
					strain[1] = value;			
				}
			} else {
				value = pload->P;
				if(pload->dir == IUX) strain[0] = value; 
				if(pload->dir == IUY) strain[1] = value; 
				if(pload->dir == IUZ) strain[2] = value; 
			}
		}
		multiply(E,strain,sigma,3,3);	
		multiply(sigma,-1,3);
		add(sigma0,sigma,3);
		/*end*/
		pload->RemoveFactor();
	} 
}
void SLAB::CalculateInitialFEM(VECTOR FEM) {
	
	/*shell*/
	if(section->type == ASECTION::SHELL4) {
		UINT state_type = section->state_type;

		if(state_type == ASECTION::PLANE_STRESS 
			|| state_type == ASECTION::PLANE_STRAIN 
			|| state_type == ASECTION::SHELL_ACTION) {
			
			section->type = MEMBRANE_SHELL4;
			section->state_type = (state_type == ASECTION::PLANE_STRAIN) 
				? ASECTION::PLANE_STRAIN : ASECTION::PLANE_STRESS;
			CalculateInitialFEM(FEM);
		}
		
		if(state_type == ASECTION::PLATE_ACTION 
			|| state_type == ASECTION::SHELL_ACTION) {
			section->type = ASECTION::PLATE4;
			section->state_type = ASECTION::PLATE_ACTION;
			CalculateInitialFEM(FEM);
		}

		section->type = ASECTION::SHELL4;
		section->state_type = state_type;
		return;
	} 

	/*calculate inital FEM*/
	UINT i,j,ig;
	DOUBLE sigma0[3];
	MATRIX Bavg,BavgT,E;
	UINT TOTAL = 6 * NJ;
	UINT SIZE = 2 * NJ;
	if(section->type == ASECTION::PLATE4 || section->type == ASECTION::Q4D)
		SIZE = 3 * NJ;

	mat_alloc(Bavg,3,SIZE);
	mat_alloc(BavgT,TOTAL,3);
	mat_alloc(E,3,3);

	/*calcule*/
	CalculateUniformFactors(Bavg[0],STRESS_F);
	GetE(E);
	if(section->type == ASECTION::PLATE4) {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 3) * 6 + plate_dir[(i % 3)];
			for(j = 0;j < 3;j++)		
				BavgT[ig][j] = Bavg[j][i];
		}
	} else if(section->type == ASECTION::Q4D) {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 3) * 6 + membrane_dir[(i % 3)];
			for(j = 0;j < 3;j++)
				BavgT[ig][j] = Bavg[j][i];
		}
	} else {
		for(i = 0;i < SIZE;i++) {
			ig = (i / 2) * 6 + membrane_dir[(i % 2)];
			for(j = 0;j < 3;j++)
				BavgT[ig][j] = Bavg[j][i];
		}
	}
	 /*Bt * sigma*/
	CalculateInitialStress(sigma0,E);
	VECTOR FEMt;
	vec_alloc(FEMt,TOTAL);
	multiply(BavgT,sigma0,FEMt,TOTAL,3);
	add(FEM,FEMt,TOTAL);
	vec_free(FEMt);

	/*free*/
	vec_free(Bavg);
	vec_free(BavgT);
	vec_free(E);

}
void SLAB::CalculateUniformFactors(VECTOR Factors,int type) {
	UINT i,l,m;
	DOUBLE N[9],J[4];
	DOUBLE h = section->h;
	DOUBLE Jdet,e,n,ew,nw,c;
	UINT RULE = Rule[NJ];
	UINT FSIZE;

	UINT SIZE = 2 * NJ;
	if(section->type == ASECTION::PLATE4 || section->type == ASECTION::Q4D)
		SIZE = 3 * NJ;

	MATRIX B,Bavg;
	mat_alloc(B,3,SIZE);
	mat_alloc(Bavg,3,SIZE);

	if(type == STRESS_F) FSIZE = 3 * SIZE;
	else FSIZE = 9;

	clear(Factors,FSIZE);

	for(l = 0;l < RULE;l++) {
		e = gauss_points[RULE - 1][l];
		ew = gauss_weights[RULE - 1][l];

        for(m = 0;m < RULE;m++) { 
			n = gauss_points[RULE - 1][m];
			nw = gauss_weights[RULE - 1][m];
			
			if(type == STRESS_F) {
				GetStrainDisplacementMatrix(B,Jdet,e,n);
				c = ew * nw * Jdet;
				if(section->type <= ASECTION::Q4D) 
					c *= h;
				multiply(B,c,3,SIZE);
				add(Bavg,B,3,SIZE);

			} else {
				GetShapeFunction(N,e,n,NJ);
				if(type == MASS_F) {
					for(i = 0;i < NJ;i++) 
						N[i] *= N[i];
				}
				
				CalculateJacobian(J,Jdet,e,n);
				c = ew * nw * Jdet * h;
				multiply(N,c,NJ);
				add(Factors,N,NJ);
			}
		}
	}

	if(type == STRESS_F) {
		equ(Factors,Bavg[0],FSIZE);
	} else if(type == MASS_F) {
		DOUBLE total = 0;
		for(i = 0;i < NJ;i++) 
			total += Factors[i];
		for(i = 0;i < NJ;i++)
			Factors[i] = Factors[i] / total;
	} else if(type == FORCE_F) {
		DOUBLE factor = 1.0 / (GetArea() * h);
		multiply(Factors,factor,NJ);
	}

	vec_free(B);
	vec_free(Bavg);
}
void SLAB::CalculateLocalFEM(LOAD* pload,VECTOR Aml,VECTOR Factor) {
	UINT TOTAL = 6 * NJ;
	int dir = pload->dir;
	DOUBLE P = pload->P * GetArea();

	clear(Aml,TOTAL);
	for(UINT i = 0;i < NJ;i++) {
		Aml[6 * i + dir] = -P  * Factor[i];
	}
}
void SLAB::CalculateFEM(VECTOR FEM) {
	DOUBLE Factor[9];
	UINT i,TOTAL = 6 * NJ;

	/*calculate uniform factors*/
	CalculateUniformFactors(Factor,FORCE_F);

	/*load rotation matrices*/
	LOAD* pload;
	VECTOR Aml;
	vec_alloc(Aml,TOTAL);

	POSITION pos = load.GetHeadPosition();
	while(pos) {
		pload = &load.GetNext(pos);

		/*slab load*/
		if(!pload->is_slab_load())
			continue;

		/*apply factor*/
		if(!pload->is_considered())
			continue;

		pload->ApplyFactor();

		/*transform global load in to local load*/
		if(pload->system) {
			MATRIX T;
			VECTOR Am,Am11,Am12,Am2;
			vec_alloc(Am,TOTAL);
			vec_alloc(Am11,TOTAL);
			vec_alloc(Am12,TOTAL);
			vec_alloc(Am2,TOTAL);
			mat_alloc(T,TOTAL,TOTAL);


			GetLoadTransMatrix(T,pload->system);
			Am[pload->dir] = pload->P;
			multiply(T,Am,Am11,TOTAL,TOTAL);
			Am[pload->dir] = pload->P1;
			multiply(T,Am,Am12,TOTAL,TOTAL);
			
			for(i = 0;i < TOTAL;i++) {
				if(!EQUAL(Am11[i],0)) {
					LOAD mload = *pload;
					mload.P = Am11[i];
					mload.P1 = Am12[i];
					mload.dir = IUX + (i % 6);
					mload.system = 0;
				    CalculateLocalFEM(&mload,Aml,Factor);
					add(FEM,Aml,TOTAL);
				}
			}

			vec_free(Am2);
			vec_free(Am11);
			vec_free(Am12);
			vec_free(Am);
			vec_free(T);
		} else {
			VECTOR Am;
			vec_alloc(Am,TOTAL);
 		    CalculateLocalFEM(pload,Aml,Factor);
			add(FEM,Aml,TOTAL);
			vec_free(Am);
		}
		/*remove factor*/
		pload->RemoveFactor();
		/*end*/
	}

	/*temperature and fabrication errors*/
	CalculateInitialFEM(FEM);
		
	/*free*/
	vec_free(Aml);
}
void SLAB::ApplyLoad(VECTOR Q) {
	UINT i,j;
	UINT TOTAL = 6 * NJ;
	VECTOR FEM,temp;
	MATRIX T;

    vec_alloc(FEM,TOTAL);
	vec_alloc(temp,TOTAL);
	mat_alloc(T,TOTAL,TOTAL);

	CalculateFEM(temp);
	GetTransformationMatrix(T);
	transpose(T,TOTAL,TOTAL);
	multiply(T,temp,FEM,TOTAL,TOTAL);

	for(j = 0;j < NJ;j++) {
		for(i = IUX;i <= IRZ;i++) {
			if(jt[j]->number[i] != INVALID)
				  Q[jt[j]->number[i]] -= FEM[6 * j + i];
		}
	}

	vec_free(FEM);
	vec_free(temp);
	vec_free(T);
}