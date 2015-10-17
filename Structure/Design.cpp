#include "common.h"

void DESIGN::InitBeam(SECTION* psec) {

	/*section dimension*/
	REBAR* rebar = &psec->rebar;
	design = rebar->design;
	b = psec->w;
	h = psec->h;
	if(psec->type == CIRCULAR) {
		r = psec->r;
		b = h = 2 * r;
	}
	hp = rebar->cover + rebar->stirrup_barsize + rebar->barsize / 2;
	bp = rebar->cover + rebar->stirrup_barsize + rebar->barsize / 2;
	sectype = psec->type;
	rebtype = rebar->type;

	/*steel and concrete*/
	Es = psec->material->Es;
	fcd = 0.85 * (psec->material->fck) / detail->fs_concrete;
	fctd = (psec->material->fctk) / detail->fs_concrete;
	fyd = psec->material->fyk / detail->fs_steel;

}

void DESIGN::InitColumn(SECTION* psec) {

	/*common initialization*/
	InitBeam(psec);

	/*reinforcement arrangement*/
	REBAR* rebar = &psec->rebar;
	DOUBLE nz = rebar->nz,ny = rebar->ny,nt = rebar->nt,delta;
	int i,cc;

	cc = 0;
	if(rebtype == RECTANGULAR) {
		for(i = 0;i < nz;i++) {
			delta = (b - 2 * bp) / (nz - 1);
			bar[cc].x = bp + i * delta;
			bar[cc].y = hp;
			cc++;
			bar[cc].x = bp + i * delta;
			bar[cc].y = h - hp;
			cc++;
		}
		for(i = 1;i < ny - 1;i++) {
			delta = (h - 2 * hp) / (ny - 1);
			bar[cc].x = bp;
			bar[cc].y = hp + i * delta;
			cc++;
			bar[cc].x = b - bp;
			bar[cc].y = hp + i * delta;
			cc++;
		}
	} else if(rebtype == CIRCULAR) {
		DOUBLE m;
		for(i = 0;i < nt;i++) {
			m = ((2.0f * i) / nt) * PI;
			bar[cc].x = (r - hp) * sin(m);
			bar[cc].y = (r - hp) * cos(m);
			cc++;
		}
	}
	bcount = cc;
}

/*Do shear design*/
void DESIGN::DoShear(DOUBLE As,BOOL invert) {
	DOUBLE Vc,k1,k2,Vrd,Av,Ac,d,bw;
	UINT nv,smax,sv;
	BOOL success;

	success = FALSE;
	Av = 0;
	
	/*major direction*/
	bw = b;
	d = h - hp;
	if(sectype == RECTANGULAR) {
		Ac = b * h;
	} else if(sectype == CIRCULAR) {
		Ac = PI * r * r;
	}
		
	if(design == COULMN || !invert) {
		k1 = (1 + 50 * (As) / (bw * d));
		k1 = min(k1 , 2);
		k2 = 1.6 - d;
	    k2 = max(k2, 1.0);
		Vrd = 0.25 * fcd * bw * d;
		if(vh >= Vrd) goto END;
		Vc = 0.25 * fctd * k1 * k2 * bw * d;
		Vc += 0.10 * bw * d * (-vt) / Ac;
		if(vh >= Vc) Av += (vh - Vc) / (d * fyd);
	}

	/*minor direction*/
	bw = h;
	d = b - bp;
	if(sectype == RECTANGULAR) {
		Ac = b * h;
	} else if(sectype == CIRCULAR) {
		Ac = PI * r * r;
	}

	if(design == COULMN || invert) {
		k1 = (1 + 50 * (As) / (bw * d));
		k1 = min(k1 , 2);
		k2 = 1.6 - d;
	    k2 = max(k2, 1.0);
		Vrd = 0.25 * fcd * bw * d;
		if(vb >= Vrd) goto END;
		Vc = 0.25 * fctd * k1 * k2 * bw * d;
		Vc += 0.10 * bw * d * (-vt) / Ac;
		if(vb >= Vc) Av += (vb - Vc) / (d * fyd);
	} 
	
	/*determine stirrup spacing*/
	DOUBLE barsize;
	POSITION pos;
    pos = detail->shearbarlist.GetHeadPosition();
    while(pos) {
		barsize = detail->shearbarlist.GetNext(pos);

		nv = UINT(ceil((4 * Av) / (PI * barsize * barsize)));
		if(invert) {
			if(vb <= (2 / 3.0) * Vrd) smax = min(300, int(500 * d));
			else smax = min(200, int(300 * d));
		} else {
			if(vh <= (2 / 3.0) * Vrd) smax = min(300, int(500 * d));
			else smax = min(200, int(300 * d));
		}
		if(nv) {
			sv = 1000 / nv;
		} else {
			sv = smax;
            Av = (1000 / sv) * (PI * barsize * barsize) / 4;
		}
		sv = min(sv, smax);
		sv = (sv / 5) * 5;

		if(sv < detail->min_spacing) {
			continue;
		}

		success = TRUE;
		pdesigndata->As = Av;
		pdesigndata->sd = barsize;
		pdesigndata->s = sv;
		break;
	}
END:
	if(!success)
		pdesigndata->flags |= DESIGNDATA::INSUFF_SHEAR;
}
/*
Rectangular beam design
*/
void DESIGN::DoBeam(DOUBLE red_factor) {
	DOUBLE tvt,tuht,tubt;
	DOUBLE As,Amin,Amax,dmax;

	vt = vt / (fcd * b * h);
	uht = uht / (fcd * b * h * h);
	ubt = ubt / (fcd * b * h * b);
	tvt = vt;
	tuht = uht;
	tubt = ubt;
	vt = fabs(vt);
	uht = fabs(uht);
	ubt = fabs(ubt);
	vh = fabs(vh);
	vb = fabs(vb);
	
	DOUBLE kx,kz,z,mu,mup,ac,bc,
		kmax,Nsd,Msds,Mups,
		Asp,esp,fsp,red = 1 - red_factor;

	
	BOOL invert = (ubt > uht);
	BOOL negative = (tuht < 0);
	if(invert) negative = (tubt < 0);
	if(invert) dmax = b - bp;
	else dmax = h - hp;

	if(invert) beamw = h;
	else beamw = b;

	/*longitudinal reinforcement*/
	Amin = (0.6e+3 / fyd) * b * h;
	Amax = 0.04 * b * h;
	kmax = 0.8 * (red - 0.44);
	
	if(invert) kx = FindKx(ubt,dmax);
	else kx = FindKx(uht,dmax);
	
	if(kx == 0) {
		As = Amin;
		Asp = Amin;
	} else {
		CalcAB(kx,dmax,ac,bc,esp);
		kz = (1 - bc);
		z = kz * dmax;
		mu = ac * (1 - bc);
		
		if(invert) Msds = (negative ? -1 : 1) * mu * fcd * b * h * b;
		else Msds = (negative ? -1 : 1) * mu * fcd * b * h * h;
		Nsd = tvt * fcd * b * h;
		
		if(kx <= kmax) {
			As = fabs(Msds / (z * fyd) + Nsd / fyd);
			As = max(Amin,As);
			Asp = Amin;
		} else {
			fsp = Es * esp / 1000;
			if(esp < 0) fsp = max(-fyd,fsp);
			else fsp = min(fyd,fsp);
			
			mup = CalcMu(kmax,dmax);
			if(invert) {
				Mups = mup * fcd * b * h * b;
				As = fabs(Mups / (z * fyd) + (Msds - Mups) / (fyd * (dmax - bp)) + Nsd / fyd);
				Asp = fabs((Msds - Mups) / (fsp * (dmax - bp)));
			} else {
				Mups = mup * fcd * b * h * h;
				As = fabs(Mups / (z * fyd) + (Msds - Mups) / (fyd * (dmax - hp)) + Nsd / fyd);
				Asp = fabs((Msds - Mups) / (fsp * (dmax - hp)));
			}
			As = max(Amin,As);
			Asp = max(Amin,Asp);
		}
	}
		
	/*shear reinforcement*/
	DoShear(As + Asp,invert);

	/*save areas*/
	for(int j = 0;j < 2;j++) {
		if(j == 0) pdesigndata->Area[j] = negative ? Asp : As;
		else pdesigndata->Area[j] = negative ? As : Asp;
	}

	if(pdesigndata->Area[1] > Amax) pdesigndata->flags |= DESIGNDATA::INSUFF_BOTTOM;
    if(pdesigndata->Area[0] > Amax) pdesigndata->flags |= DESIGNDATA::INSUFF_TOP;
}
/*calculate alphac,betac,esp from neurtal axis depth ratio x / d*/
void DESIGN::CalcAB(DOUBLE kx,DOUBLE dmax,DOUBLE& ac,DOUBLE& bc,DOUBLE& esp) {
    DOUBLE ecm,es,x;

	x = kx * dmax;
	if(kx < (3.5 / 13.5)) {
		es = -10;
		ecm = -es * x / (dmax - x);
		if(ecm <= 2.0) {
			ac = ecm * (6 - ecm) * kx / 12;
			bc = (8 - ecm) * kx / (4 * (6 - ecm));
		} else {
			ac = (3 * ecm - 2) * kx / (3 * ecm);
			bc = (ecm * (3 * ecm - 4) + 2) * kx / (2 * ecm * (3 * ecm - 2));
		}
	} else {
		ecm = 3.5;
		ac = (3 * ecm - 2) * kx / (3 * ecm);
		bc = (ecm * (3 * ecm - 4) + 2) * kx / (2 * ecm * (3 * ecm - 2));
	}
	esp = ecm * (x - hp) / x;
}
/*calculate mu*/
DOUBLE DESIGN::CalcMu(DOUBLE kx,DOUBLE dmax) {
	DOUBLE ac,bc,esp;
	CalcAB(kx,dmax,ac,bc,esp);
	return (ac * (1 - bc));
}
/*find kx*/
DOUBLE DESIGN::FindKx(DOUBLE u,DOUBLE dmax) {
	DOUBLE x,y,c,fx,fy,fc,error;

	x = 0.0;
	y = 1.2;

	fx = CalcMu(x,dmax) - u;
	if(fabs(fx) <= 0.00001)
		return x;

	fy = CalcMu(y,dmax) - u;
	if(fabs(fy) <= 0.00001)
		return y;

	do {

		c = (x + y) / 2;
		fc = CalcMu(c,dmax) - u;
		if(fabs(fc) <= 0.00001)
			return c;
		
		error = fabs(y - x) / 2;

		if(fx * fc < 0) {
			y = c;
			fy = fc;
		} else if(fx * fc > 0) {
			x = c;
			fx = fc;
		} else {
			break;
		}

	} while(error >= 0.000001);

	return c;

}
/*calculate distance from neutral axis*/
DOUBLE DESIGN::CalcD(DOUBLE x,DOUBLE y,DOUBLE q) {
	DOUBLE d = 0;
	if(sectype == CIRCULAR) {
		DOUBLE cx,cy;
		cx = r * sin(q);
		cy = r * cos(q);
		x = cx - x;
		y = cy - y;
        d = x * sin(q) + y * cos(q); 
		return d;
	} else if(sectype == RECTANGULAR) {
		DOUBLE q1;
		q1 = atan(x / y);
		d = sqrt(x * x + y * y);
		d = d * cos(q1 - q);
		return d;
	}
	return d;
}
/*
Rectangular/Circular column design
*/
void DESIGN::DoColumn(DOUBLE* k2) {
	DOUBLE tvt,tuht,tubt;
	DOUBLE w,Amin,Amax,As,q;

	/*get ready*/
	DOUBLE Ag;
	if(sectype == RECTANGULAR) {
		Ag = b * h;
	} else if(sectype == CIRCULAR) {
		Ag = PI * r * r;
	}
	vt = vt / (fcd * Ag);
	uht = uht / (fcd * Ag * h);
	ubt = ubt / (fcd * Ag * b);
	

	tvt = vt;
	tuht = uht;
	tubt = ubt;
	vt = fabs(vt);
	uht = fabs(uht);
	ubt = fabs(ubt);

	/*find area of rebar*/
	if(EQUAL(uht,0)) q = PI / 2;
	else if(EQUAL(ubt,0)) q = 0;
	else q = atan(ubt / uht);

	w = FindW(q);
	Amin = 0.008 * Ag;
	Amax = 0.08 * Ag;
	As = w * Ag * fcd / fyd;
	if(As < Amin) {
		As = Amin;
		w = As * fyd / (Ag * fcd);
	}

	/*balanced moment*/
	DOUBLE d,dmax,x;
	dmax = 0;
	for(int xx = 0;xx < bcount;xx++) {
		d = CalcD(bar[xx].x,bar[xx].y,q);
		if(d > dmax) dmax = d;
	}

	x = (3.5 / 5.5) * dmax;

	DOUBLE Pu,Muz,Muy;
	FindMatch(x,q,w,Pu,Muz,Muy);
	
	if(uht == 0 || Muz == 0) k2[0] = 1.0;
	else k2[0] = uht / Muz;
	if(ubt == 0 || Muy == 0) k2[1] = 1.0;
	else k2[1] = ubt / Muy;

	/*save areas*/
	pdesigndata->Area[0] = As;
    if(pdesigndata->Area[0] > Amax) 
		pdesigndata->flags |= DESIGNDATA::INSUFF_TOP;

	/*ratios*/
	pdesigndata->ratios[0] = vt;
	pdesigndata->ratios[1] = uht;
	pdesigndata->ratios[2] = ubt;
	pdesigndata->ratios[3] = w;
}
/*
caluclate biaxial capacity
*/
void DESIGN::CalculateBiaxial(DOUBLE x,DOUBLE q,DOUBLE w,DOUBLE* Pu,DOUBLE* Muz,DOUBLE* Muy) {
	DOUBLE es,fs,ecm,
	       Ac,xc,yc,Ag,dr,
		   kx,As,xmax;

	/*xmax*/
	if(sectype == RECTANGULAR) {
		xmax = CalcD(b,h,q);
	} else if(sectype == CIRCULAR) {
		xmax = b;
	}
	
	/*find furthest rebar distance*/
	DOUBLE d,dmax;
	dmax = 0;
	for(int i = 0;i < bcount;i++) {
		d = CalcD(bar[i].x,bar[i].y,q);
		if(d > dmax) dmax = d;
	}

	/*find ecm*/
	kx = x / dmax;
	
	if(kx < (3.5 / 13.5)) {
		ecm = 10 * x / (dmax - x);
	} else {
		if(x < xmax) {
			ecm = 3.5;
		} else {
			ecm = 2.0 / (1 - (3 * xmax) / (7 * x));
		}
	}
	/*Compressed concrete area [Rectangular approximation]*/
	DOUBLE a = 0.8 * x;

	if(sectype == RECTANGULAR) {

		Ag = b * h;
		As = w * Ag * fcd / (bcount * fyd);

		if(a > xmax) {
			Ac = Ag;
			xc = b / 2;
			yc = h / 2;
		} else {
			if(q == 0) {
				Ac = b * a;
				xc = b / 2;
				yc = a / 2;
			} else if(q == (PI / 2)) {
				Ac = h * a;
				xc = a / 2;
				yc = h / 2;
			} else {
				DOUBLE xh,xb,A,Acx,Acy,xc1,yc1;
				xh = a / cos(q);
				xb = a / sin(q);
				Ac = xh * xb / 2;
				xc = xb / 3;
				yc = xh / 3;
				Acx = Ac * xc;
				Acy = Ac * yc;
				if(xh > h) {
					A = -xb * pow(xh - h, 2) / (2 * xh);
					xc1 = xb * (xh - h) / (3 * xh);
					yc1 = h + (xh - h) / 3;
					Ac += A;
					Acx += A * xc1;
					Acy += A * yc1;
				}
				if(xb > b) {
					A = -xh * pow(xb - b, 2) / (2 * xb);
					xc1 = b + (xb - b) / 3;
					yc1 = xh * (xb - b) / (3 * xb);
					Ac += A;
					Acx += A * xc1;
					Acy += A * yc1;
				}
				xc = Acx / Ac;
				yc = Acy / Ac;
			}
		}
		
		*Pu = fcd * Ac;
		*Muz = fcd * Ac * (h / 2 - yc);
		*Muy = fcd * Ac * (b / 2 - xc);
		
		/*contribution of rebars*/
		for(int i = 0;i < bcount;i++) {
			dr = CalcD(bar[i].x,bar[i].y,q);
			es = ecm * (x - dr) / x;
			fs = Es * es / 1000;
			if(es < 0) fs = max(-fyd,fs);
			else fs = min(fyd,fs);
			
			*Pu += As * fs;
			*Muz += As * fs * (h / 2 - bar[i].y); 
			*Muy += As * fs * (b / 2 - bar[i].x); 
		}
		
		*Pu = *Pu / (Ag * fcd);
		*Muz = *Muz / (Ag * fcd * h);
		*Muy = *Muy / (Ag * fcd * b);
	} else {
		Ag = PI * r * r;
		As = w * Ag * fcd / (bcount * fyd);

		if(a > xmax) {
			Ac = Ag;
			xc = 0;
			yc = 0;
		} else {
			DOUBLE q1,rc;
			q1 = acos((r - a) / r);
			rc = r * (1 - 2 * pow(sin(q1), 3) / (3 * q1 - 3 * sin(q1) * cos(q1)));
			Ac = r * r * (q1 - sin(q1) * cos(q1));
			xc = (r - rc) * sin(q);
			yc = (r - rc) * cos(q);
		}

		*Pu = fcd * Ac;
		*Muz = fcd * Ac * yc;
		*Muy = fcd * Ac * xc;
		
		/*contribution of rebars*/
		for(int i = 0;i < bcount;i++) {
			dr = CalcD(bar[i].x,bar[i].y,q);
			es = ecm * (x - dr) / x;
			fs = Es * es / 1000;
			if(es < 0) fs = max(-fyd,fs);
			else fs = min(fyd,fs);
			
			*Pu += As * fs;
			*Muz += As * fs * bar[i].y; 
			*Muy += As * fs * bar[i].x; 
		}
		
		*Pu = *Pu / (Ag * fcd);
		*Muz = *Muz / (Ag * fcd * h);
		*Muy = *Muy / (Ag * fcd * b);
	}

    /*truncate after 8 digits (helps convergence)*/
	*Pu = int(*Pu * 1e8) / 1.0e8;
	*Muz = int(*Muz * 1e8) / 1.0e8;
	*Muy = int(*Muy * 1e8) / 1.0e8;
}

void DESIGN::FindMatch(DOUBLE x,DOUBLE q,DOUBLE w,DOUBLE& v,DOUBLE& uh,DOUBLE& ub) {
	
	DOUBLE l, m, c, fl, fm, fc, error;
	DOUBLE vl,uhl,ubl;
	DOUBLE vm,uhm,ubm;
	DOUBLE vc,uhc,ubc;

	l = q;
	CalculateBiaxial(x,l,w,&vl, &uhl, &ubl);
	if(uhl == 0) fl = PI / 2 - q;
	else if(ubl == 0) fl = 0 - q;
	else fl = atan(ubl / uhl) - q;
	if(fabs(fl) <= 0.0001) {
		v = vl;
		uh = uhl;
		ub = ubl;
		return;
	}

	for(int i = 0;i < 16;i++) {
	    m = (i / 32.0f) * PI;
        CalculateBiaxial(x,m,w,&vm, &uhm, &ubm);
		if(uhm == 0) fm = PI / 2 - q;
		else if(ubm == 0) fm = 0 - q;
		else fm = atan(ubm / uhm) - q;
		if(fabs(fm) <= 0.0001) {
			v = vm;
			uh = uhm;
			ub = ubm;
			return;
		}
		if(fl * fm < 0) {
			break;
		}
	}

	do {
		c = (l + m) / 2;
		CalculateBiaxial(x,c,w,&vc, &uhc, &ubc);
	    if(uhc == 0) fc = PI / 2 - q;
	    else if(ubc == 0) fc = 0 - q;
	    else fc = atan(ubc / uhc) - q;
		if(fabs(fc) <= 0.0001) {
			break;
		}

		error = fabs(m - l) / 2;

		if(fl * fc < 0) {
			m = c;
			fm = fc;
		} else if(fl * fc > 0) {
			l = c;
			fl = fc;
		} else {
			break;
		}

	} while(error >= 0.000001);

    v = vc;
	uh = uhc;
	ub = ubc;
}
void DESIGN::FindNearest(DOUBLE q,DOUBLE w,DOUBLE& dmin,DOUBLE& error) {
	DOUBLE v,uh,ub,x,y,c,qx,qy,qc,qt,d1,dt,err;

	dt = sqrt(vt * vt + uht * uht + ubt * ubt);
    qt = atan(vt / sqrt(uht * uht + ubt * ubt));

	y = 1e6;
	FindMatch(y,q,w,v,uh,ub);
	if(uh + ub == 0) qy = qc =  PI / 2 - qt;
	else qy = qc = atan(v / sqrt(uh * uh + ub * ub)) - qt;
	if(fabs(qc) <= 0.001) 
		goto END;

	x = 0.0001;
	FindMatch(x,q,w,v,uh,ub);
	if(uh + ub == 0) qx = qc =  PI / 2 - qt;
	else qx = qc = atan(v / sqrt(uh * uh + ub * ub)) - qt;
	if(fabs(qc) <= 0.001) 
		goto END;

	if(qx * qy > 0) 
		goto END;

	do {
		c = (x + y) / 2;
		FindMatch(c,q,w,v,uh,ub);
		if(uh + ub == 0) qc =  PI / 2 - qt;
	    else qc = atan(v / sqrt(uh * uh + ub * ub)) - qt;
		if(fabs(qc) <= 0.001) 
			goto END;

		err = fabs(y - x) / 2;
        
		if(qx * qc < 0) {
			y = c;
			qy = qc;
		} else if(qx * qc > 0) {
			x = c;
			qx = qc;
		} else {
			goto END;
		}

	} while(err >= 0.000001);

END:
	d1 = sqrt(v * v + uh * uh + ub * ub);
	dmin = d1 * cos(qc) - dt;
	error = sqrt(pow(v - vt, 2) + pow(uh - uht, 2)  + pow(ub - ubt, 2));
}

DOUBLE DESIGN::FindW(DOUBLE q) {
	DOUBLE x,y,c,fx,fy,fc,dx,dy,dc,error;
	
	x = 0.0;
	y = 3.0;

	FindNearest(q,x,fx,dx);
	if(dx <= 0.00001)
		return x;
	FindNearest(q,y,fy,dy);
	if(dy <= 0.00001)
		return y;

	if(fx * fy > 0) {
		if(fx < fy) return x;
		else return y; 
	}

	do {
		c = (x * dy + y * dx) / (dx + dy);
		FindNearest(q,c,fc,dc);
		if(dc <= 0.00001)
			return c;
		
		error = fabs(y - x) / 2;

		if(fx * fc < 0) {
			y = c;
			fy = fc;
			dy = dc;
		} else if(fx * fc > 0) {
			x = c;
			fx = fc;
			dx = dc;
		} else {
			break;
		}

	} while(error >= 0.000001);

    return c;
}
/*
select bar diameters for the give area
*/
void DetermineBars(BARSIZELIST* plist,DOUBLE As,DESIGNDATA* pdata,int which,UINT nt)  {
  
	POSITION pos;
	DOUBLE barsize,Area,found;
	UINT i,curr;

	found = FALSE;
	for(i = 0;i < 50; i++) {
		pos = plist->GetHeadPosition();
		while(pos) {
			barsize = plist->GetNext(pos);

			Area = (nt * pdata->count[which][0] * PI * pdata->diam[which][0] * pdata->diam[which][0]) / 4;
            if(pdata->total[which] == 2)   
				Area += (nt * pdata->count[which][1] * PI * pdata->diam[which][1] * pdata->diam[which][1]) / 4;
			if(Area >= As)
				return;
			Area += (nt * PI * barsize * barsize) / 4;
			if(Area >= As) {
				found = TRUE;
				break;
			}
		}

		if(!pdata->total[which]) {
			pdata->count[which][0] = 1;
			pdata->diam[which][0] = barsize;
			pdata->total[which] = 1;
			if(found) break;
			continue;
		}

		curr = pdata->total[which] - 1;
		
		if(EQUAL(barsize,pdata->diam[which][curr])) {
			pdata->count[which][curr]++;
		} else {
			if(curr == 0) {
				pdata->count[which][1]++;
				pdata->diam[which][1] = barsize;
				pdata->total[which] = 2;
			} else {
				pdata->count[which][1]++;
			}
		}
		if(found) break;
    }
}
/*
member design
*/
void MEMBER::Design(DETAILING* detail) {
	if(e_section) return;
	if(section->material->type == MATERIAL::STEEL) {
		DesignSteelColumn(detail);
	} else if(section->material->type == MATERIAL::CONCRETE) {
		if(section->rebar.design == COULMN) {
			if(section->type == RECTANGULAR || section->type == CIRCULAR)
				DesignColumn(detail);
		} else {
			if(section->type == RECTANGULAR)
				DesignBeam(detail);
		}
	}
}
void MEMBER::DesignBeam(DETAILING* detail) {
	UINT i,j;
	DESIGN design;
	design.detail = detail;
	design.InitBeam(section);

	/*zero*/
	for(i = 0;i < nDiv;i++)
		design_data[i].flags = 0;

	/*minimum reinforcement*/
	DESIGNDATA* pdata = &design_data[0];
	DOUBLE Amin = (0.6e+3 / design.fyd) * design.b * design.h;
	POSITION pos;
	DOUBLE barsize = 0;
	BOOL found = FALSE;
	for(i = 1;i < 50;i++) {
		pos = detail->beambarlist.GetHeadPosition();
		while(pos) {
			barsize = detail->beambarlist.GetNext(pos);
			if(Amin <= 2 * i * (PI * barsize * barsize) / 4) {
				found = TRUE;
				break;
			}
		}
		if(found) break;
	}
	pdata->total[0] = 1;
	pdata->total[1] = 1;
	pdata->count[0][0] = pdata->count[1][0] = 2 * i;
	pdata->diam[0][0] = pdata->diam[1][0] = barsize;

	for(i = 1;i < nDiv;i++) {
		design_data[i] = design_data[0];
	}
	/*design for each section*/
	for(i = 0;i < nDiv;i++) {
		design.pdesigndata = &design_data[i];
		design.vt = forces[i][IUX];
		design.uht = forces[i][IRY];
		design.ubt = forces[i][IRZ];
		design.vh = forces[i][IUZ];
		design.vb = forces[i][IUY];
		design.DoBeam(redistribution_factor);
	}
	/*calculate bars*/
	int n;
	DOUBLE width;
	DOUBLE Area;
	
	found = FALSE;
	pos = detail->beambarlist.GetHeadPosition();
	while(pos) {
		barsize = detail->beambarlist.GetNext(pos);

		for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			for(j = 0;j < 2;j++) {
				Area = pdata->Area[j];
				Area -= (pdata->count[j][0] * PI * pdata->diam[j][0] * pdata->diam[j][0]) / 4;
				if(Area <= 0)
					continue;
				n = int(ceil((4 * Area) / (PI * barsize * barsize)));
				
				/*check if it fits in the beam's width*/
				width = section->rebar.cover * (n + pdata->count[j][0] - 1) 
					+ pdata->count[j][0] * pdata->diam[j][0] + n * barsize;
				if(width > (detail->max_layers + 1) * design.beamw) {
					goto BOTTOM;
				}
			}
		}

		found = TRUE;
		for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			for(j = 0;j < 2;j++) {
				Area = pdata->Area[j];
				Area -= (pdata->count[j][0] * PI * pdata->diam[j][0] * pdata->diam[j][0]) / 4;
				if(Area <= 0)
					continue;
				n = int(ceil((4 * Area) / (PI * barsize * barsize)));
				if(EQUAL(barsize,pdata->diam[j][0])) {
					pdata->count[j][0] += n;
				} else {
					pdata->total[j] = 2;
					pdata->count[j][1] = n;
					pdata->diam[j][1] = barsize;
				}
			}
		}
		break;
BOTTOM:;
	}

	/*width insufficient*/
	if(!found) {
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			pdata->flags |= DESIGNDATA::TOP_DOUBLE_ROW;
		}
	}

	/*cut at L/3*/
	DOUBLE x;
	DOUBLE L = GetLength();
	if(detail->cut_at_L3) {
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
	        x = station[i];
			if(x <= L / 3) *pdata = design_data[0];
			else if(x <= 2 * L / 3) *pdata = design_data[nDiv / 2];
			else *pdata = design_data[nDiv - 1];
		}
	}
	/*single shear reinforcement*/
	if(detail->single_shear) {
		DOUBLE Amax = 0;
		int imax;
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			if(pdata->As > Amax)  {
				Amax = pdata->As;
				imax = i;
			}
		}
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			pdata->As = design_data[imax].As;
            pdata->s = design_data[imax].s;
			pdata->sd = design_data[imax].sd;
		}
	}
	/*print result*/
	print(" X                   Top                             Bottom                      Shear         \n"
		"===============================================================================================\n");
	for(j = 0;j < nDiv;j++) {
		pdata = &design_data[j];
		print("%4.2f ",station[j]);

		if(pdata->flags) {
			print("\tInsufficient!\n");
		} else {
			CString str;
			for(i = 0;i < 2;i++) {
				if(pdata->total[i] > 1) {
					str.Format("||  %3dD%-2d and %3dD%-2d [%-4dmm2]  ",
						pdata->count[i][0],int(pdata->diam[i][0] * 1000),
						pdata->count[i][1],int(pdata->diam[i][1] * 1000),
						int(pdata->Area[i] * 1e6));
				} else {
					str.Format("||      %3dD%-2d        [%-4dmm2]  ",
						pdata->count[i][0],int(pdata->diam[i][0] * 1000),
						int(pdata->Area[i] * 1e6));
				}
				print(str.GetBuffer(0));
			}
			print("||  D%-2dc/c%4d [%-4dmm2]\n",
				int(pdata->sd * 1000),pdata->s,int(pdata->As * 1e6));
		}
	}
	print("\n");
}
void MEMBER::DesignColumn(DETAILING* detail) {
	UINT i,j,nt;
    DESIGN design;
	BOOL rectangular = section->rebar.type == RECTANGULAR;
	design.detail = detail;
	design.InitColumn(section);

	/*zero*/
	for(i = 0;i < nDiv;i++)
		design_data[i].flags = 0;

    /*reinforcement*/
    DOUBLE r,e[2],k2[2],tk2[2],temp[3];
	for(j = 0;j < 2;j++) {
		e[j] = Le[j] / 300;
		e[j] = max(e[j] , 0.02);
	}

    if(rectangular) nt = 2 * (section->rebar.nz + section->rebar.ny - 2);
    else nt = section->rebar.nt;

	for(j = 0;j < nDiv;j++) {
		design.pdesigndata = &design_data[j];
		design.vt = -forces[j][IUX];
		design.uht = forces[j][IRY];
		design.ubt = forces[j][IRZ];
		design.vh = forces[j][IUZ];
		design.vb = forces[j][IUY];

		/*accidental eccentricity*/
		if(detail->accidental_eccentricity) {
			if(e[0] * design.vt * design.uht > 0) design.uht += e[0] * design.vt;
			else design.uht -= e[0] * design.vt;
			if(e[1] * design.vt * design.ubt > 0) design.ubt += e[1] * design.vt;
			else design.ubt -= e[1] * design.vt;
		}

		/*individual buckling effect*/
		if(detail->individual_buckling) {
			temp[0] = design.uht;
			temp[1] = design.ubt;
			temp[2] = design.vt;
			
			k2[0] = k2[1] = 1.0;
			tk2[0] = tk2[1] = 1.0;

			for(int p = 0;p < 5;p++) {
				
				DOUBLE e2[2] = {0,0};
				for(i = 0;i < 2;i++) {
					DOUBLE k1,kappa;
					if(i == 0) {
						r = Le[i] / section->ry;
						if(rectangular) kappa = k2[i] * (5 / (section->h - design.hp)) * 1e-3;
						else kappa = k2[i] * (5 / (2 * section->r - design.hp)) * 1e-3;
					} else {
						r = Le[i] / section->rz;
						if(rectangular) kappa = k2[i] * (5 / (section->w - design.hp)) * 1e-3;
						else kappa = k2[i] * (5 / (2 * section->r - design.hp)) * 1e-3;
					}
					
					if(r >= 15 && r <= 35) {
						k1 = r / 20 - 0.75;
					} else {
						k1 = 1.0;
					}
					
					e2[i] = k1 * kappa * pow(Le[i], 2) / 10;
		
					if(i == 0) {
						if(e2[i] * temp[2] * temp[i] > 0) design.uht = temp[i] + e2[i] * temp[2];
						else design.uht = temp[i] - e2[i] * temp[2];
					} else {
						if(e2[i] * temp[2] * temp[i] > 0) design.ubt = temp[i] + e2[i] * temp[2];
						else design.ubt = temp[i] - e2[i] * temp[2];
					}
				}
				design.vt = temp[2];
				
				tk2[0] = k2[0];
				tk2[1] = k2[1];
				
				design.DoColumn(k2);

				if(design.pdesigndata->flags)
					break;
				
				if((fabs(k2[0] - tk2[0]) <= 0.01)
					&& (fabs(k2[1] - tk2[1]) <= 0.01)) {
					break;
				}
			}
		} else {
			design.DoColumn(k2);
		}
		/*shear and rebar*/
		design.DoShear(design.pdesigndata->Area[0]);
		DetermineBars(&detail->columnbarlist,design.pdesigndata->Area[0],design.pdesigndata,0,nt);
	}

	/*single shear reinforcement*/
	DESIGNDATA* pdata;

	if(detail->single_shear) {
		DOUBLE Amax = 0;
		int imax;
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			if(pdata->As > Amax)  {
				Amax = pdata->As;
				imax = i;
			}
		}
        for(i = 0;i < nDiv;i++) {
			pdata = &design_data[i];
			pdata->As = design_data[imax].As;
            pdata->s = design_data[imax].s;
			pdata->sd = design_data[imax].sd;
		}
	}

	/*print result*/
	print(" X                   Longtiudinal                       Shear                    v              uh              ub              w\n"
		"=====================================================================================================================================\n");
	for(j = 0;j < nDiv;j++) {
		pdata = &design_data[j];
		print("%4.2f ",station[j]);

		if(pdata->flags) {
			print("\tInsufficient!\n");
		} else {
			CString str;
			if(pdata->total[0] > 1) {
				str.Format("||  %d X (%3dD%-2d and %3dD%-2d)  [%-4dmm2]  ",nt,
					pdata->count[0][0],int(pdata->diam[0][0] * 1000),
					pdata->count[0][1],int(pdata->diam[0][1] * 1000),
					int(pdata->Area[0] * 1e6));
			} else {
				str.Format("||          %3dD%-2d           [%-4dmm2]  ",
					nt * pdata->count[0][0],
					int(pdata->diam[0][0] * 1000),
					int(pdata->Area[0] * 1e6));
			}
			print(str.GetBuffer(0));
			print("||  D%-2dc/c%4d [%-4dmm2]  ||%12.6f  ||%12.6f  ||%12.6f  ||%12.6f  \n",
				int(pdata->sd * 1000),pdata->s,int(pdata->As * 1e6),
				pdata->ratios[0],pdata->ratios[1],pdata->ratios[2],pdata->ratios[3]);
		}
	}
	print("\n");
}
/*
STEEL design
*/
static DOUBLE imperfection[] = {0.21, 0.34, 0.49, 0.76};

void MEMBER::DesignSteelColumn(DETAILING* detail) {
    DOUBLE s1,s2,Nsd,Msdc[2];
	DOUBLE fy = section->material->fyk;
	DOUBLE fyd = fy / detail->fs_steel;
	UINT i,j;

	/*zero*/
	for(i = 0;i < nDiv;i++)
		design_data[i].flags = 0;

	/*varialbes*/
	DOUBLE X[2],Xmin,al,q,l,lp[2],ba,si[2];
	DOUBLE e = sqrt(235e+3 / fy);
	DOUBLE bm[2],x,tm1,tm2,tm3,mmaxl,mmax,mmin;
	DOUBLE L = GetLength();
	
	/*si*/
	for(i = 0;i < 2;i++) {
		if(fabs(forces[0][IRY + i]) > fabs(forces[nDiv - 1][IRY + i]))
			si[i] = forces[nDiv - 1][IRY + i] / forces[0][IRY + i];
		else if(EQUAL(forces[nDiv - 1][IRY + i],0)) {
			if(EQUAL(forces[0][IRY + i],0)) si[i] = 1;
            else if(forces[0][IRY + i] < 0) si[i] = -1;
			else si[i] = 1;
		} else
			si[i] = forces[0][IRY + i] / forces[nDiv - 1][IRY + i];
	}

	/*equivalent uniform moment factor*/
	for(i = 0;i < 2;i++) {
		mmaxl = 0;
		mmax = 0;
		mmin = 1e16;
		for(j = 0;j < nDiv;j++) {
			x = station[j];
			tm1 = forces[0][IRY + i] + (forces[nDiv - 1][IRY + i] - forces[0][IRY + i]) * (x / L);
			tm2 = fabs(forces[j][IRY + i] - tm1);
			tm3 = forces[j][IRY + i];
			if(tm2 > mmaxl) mmaxl = tm2;
			if(tm3 > mmax) mmax = tm3;
			if(tm3 < mmin) mmin = tm3;
		}
		
		if(mmax * mmin < 0)
			mmax = fabs(mmax) + fabs(mmin);
		else {
			mmax = max(fabs(mmax),fabs(mmin));
		}
		
		bm[i] = 1.8 - 0.7 * si[i];
		if(mmax)
			bm[i] = bm[i] + (mmaxl / mmax) * (1.35 - bm[i]);
	}

	print("Effective Length:\n\tY = %9.3f Z = %9.3f\n",Le[0],Le[1]);
    print("Equivalent Uniform Moment Factor:\n\tY = %9.3f Z = %9.3f\n",bm[0],bm[1]);
	print("\n");
	/*
	CHECK each section
	*/
	print("\n");
	print("  Loc       CLASS        Xy        Xz       Xmin       Axial        Major        Minor        Total      LT Axial     LT Major     LT Minor     LT Total\n");
	print("========================================================================================================================================================\n");

    DESIGNDATA* pdata;
	for(j = 0;j < nDiv;j++) {
		pdata = &design_data[j];
		Nsd = fabs(forces[j][IUX]);
		Msdc[0] = fabs(forces[j][IRY]);
		Msdc[1] = fabs(forces[j][IRZ]);
		
		/*classify section*/
		s1 = Msdc[0] / section->Zy + Nsd / section->A;
		s2 = Msdc[0] / section->Zy - Nsd / section->A;
		section->si[0] = s2 / s1;
		section->fm[0] = (s1 + s2) / 2;
		s1 = Msdc[1] / section->Zz + Nsd / section->A;
		s2 = Msdc[1] / section->Zz - Nsd / section->A;
		section->si[1] = s2 / s1;
		section->fm[1] = (s1 + s2) / 2;
		section->Classify(this);
		
		/*reduction factor*/
		ba = (1 - section->dA / section->A);
		for(i = 0;i < 2;i++) {
			if(i == 0) l = Le[i] / section->ry;
			else l = Le[i] / section->rz;
			lp[i] = (l / (93.9 * e)) * sqrt(ba);
			al = imperfection[section->SelectBucklingCurve(IUY + i)];
			q = 0.5 * (1 + al * (lp[i] - 0.2) + lp[i] * lp[i]);
			X[i] = 1 / (q + sqrt(q * q - lp[i] * lp[i]));
		}
		Xmin = min(X[0],X[1]);
		
		/*resistance to bending and axial force*/
		DOUBLE k[2],mu[2];
		for(i = 0;i < 2;i++) {
			if(section->sclass == SECTION::CLASS3 || section->sclass == SECTION::CLASS4) {
				mu[i] = lp[i] * (2 * bm[i] - 4);
				if(mu[i] > 0.9) mu[i] = 0.9;
				k[i] = 1 - mu[i] * Nsd / (X[i] * ba * section->A * fy); 
				if(k[i] > 1.5) k[i] = 1.5;
			} else {
				if(i == 0) mu[i] = lp[i] * (2 * bm[i] - 4) + (section->Sy - section->Zy) / section->Zy;
				else mu[i] = lp[i] * (2 * bm[i] - 4) + (section->Sz - section->Zz) / section->Zz;
				if(mu[i] > 0.9) mu[i] = 0.9;
				k[i] = 1 - mu[i] * Nsd / (X[i] * section->A * fy);
				if(k[i] > 1.5) k[i] = 1.5;
			}
		}
		
		if(section->sclass == SECTION::CLASS1 || section->sclass == SECTION::CLASS2) {
			pdata->ratios[0] = Nsd / (Xmin * section->A * fyd);
			pdata->ratios[1] = k[0] * Msdc[0] / (section->Sy * fyd);
			pdata->ratios[2] = k[1] * Msdc[1] / (section->Sz * fyd);
		} else if(section->sclass == SECTION::CLASS3) {
			pdata->ratios[0] = Nsd / (Xmin * section->A * fyd);
			pdata->ratios[1] = k[0] * Msdc[0] / (section->Zy * fyd);
			pdata->ratios[2] = k[1] * Msdc[1] / (section->Zz * fyd);
		} else if(section->sclass == SECTION::CLASS4) {
			pdata->ratios[0] = Nsd / (Xmin * ba * section->A * fyd);
			pdata->ratios[1] = (k[0] * Msdc[0] + Nsd * section->eA[0]) / ((section->Zy - section->dZ[0]) * fyd);
			pdata->ratios[2] = (k[1] * Msdc[1] + Nsd * section->eA[1]) / ((section->Zz - section->dZ[1]) * fyd);
		}
		
		DOUBLE sum = 0;
		for(i = 0;i < 3;i++) sum += pdata->ratios[i];
		pdata->ratios[3] = sum;
		if(sum > 1)
			pdata->flags |= DESIGNDATA::INSUFF_TOP;
		
		/*resistance to lateral torsional buckling and axial force*/
		DOUBLE klt,mult,Xlt = X[0];
		mult = 0.15 * lp[1] * bm[0] - 0.15;
		if(mult > 0.9) mult = 0.9;
		klt = 1 - mult * Nsd / (X[1] * ba * section->A * fy); 
		if(klt > 1) klt = 1;
		
		if(section->sclass == SECTION::CLASS1 || section->sclass == SECTION::CLASS2) {
			pdata->ratios[4] = Nsd / (X[1] * section->A * fyd);
			pdata->ratios[5] = klt * Msdc[0] / (Xlt * section->Sy * fyd);
			pdata->ratios[6] = k[1] * Msdc[1] / (section->Sz * fyd);
		} else if(section->sclass == SECTION::CLASS3) {
			pdata->ratios[4] = Nsd / (X[1] * section->A * fyd);
			pdata->ratios[5] = klt * Msdc[0] / (Xlt * section->Zy * fyd);
			pdata->ratios[6] = k[1] * Msdc[1] / (section->Zz * fyd);
		} else if(section->sclass == SECTION::CLASS4) {
			pdata->ratios[4] = Nsd / (X[1] * ba * section->A * fyd);
			pdata->ratios[5] = (klt * Msdc[0] + Nsd * section->eA[0]) / (Xlt * (section->Sy - section->dZ[0]) * fyd);
			pdata->ratios[6] = (k[1] * Msdc[1] + Nsd * section->eA[1]) / ((section->Sz - section->dZ[1]) * fyd);
		}
		
		sum = 0;
		for(i = 0;i < 3;i++) sum += pdata->ratios[i + 4];
		pdata->ratios[7] = sum;
		if(sum > 1)
			pdata->flags |= DESIGNDATA::INSUFF_BOTTOM;

		/*print*/
		print("%6.3f      CLASS%d %9.3f %9.3f %9.3f ",
			station[j],section->sclass + 1,X[0],X[1],Xmin);
		for(i = 0;i < 8;i++) {
			print("%12.6f ",pdata->ratios[i]);
		}
		print("\n");

	}
}
/*
calculate buckling factor
*/
static DOUBLE GetKa(DOUBLE si) {
    if(EQUAL(si,1)) {
		return 0.43;
	} else if(si > 0) {
		return (0.578 / (si + 0.34));
	} else if(EQUAL(si,0)) {
		return 1.70;
	} else if(si > -1) {
		return (1.7 - 5 * si + 17.1 * si * si);
	} else if(EQUAL(si,-1)) {
		return 23.8;
	}
	return 0;
}

/*
select buckling curve
*/
int SECTION::SelectBucklingCurve(int dir) {
	if(type == WIDEFLANGE) {
		if(built_up) {
			if(dir == IUY) {
				if(tf <= 0.04) return BB;
				else return BC;
			} else {
				if(tf <= 0.04) return BC;
				else return BD;
			}
		} else {
			if(dir == IUY) {
				if(h / w > 1.2) {
					if(tf <= 0.04) return BA;
					else return BB;
				} else {
					if(tf <= 0.1) return BB;
					else return BD;
				}
			} else {
				if(h / w > 1.2) {
					if(tf <= 0.04) return BB;
					else return BC;
				} else {
					if(tf <= 0.1) return BC;
					else return BD;
				}
			}
		}
	} else if(type == BOX || type == PIPE) {
		if(type == BOX && built_up) {
			if(w / tf < 30 && h / tw < 30) return BC;
			else return BB;
		} else {
			return BA; //NB:cold formed BB; BC;
		}
    } else {
		return BC;
	}
}
/*
classify section
*/
void SECTION::Classify(MEMBER* member) {
	DOUBLE e = sqrt(235e+3 / material->fyk);
	DOUBLE b,c,d,t,ratio,ratio1;
	int cl;

	dA = 0;
	dZ[0] = dZ[1] = 0;
	eA[0] = eA[1] = 0;

	if(type == WIDEFLANGE ||
		type == CHANNEL ||
		type == DOUBLE_CHANNEL ||
		type == TEE
		) {

		if(type == CHANNEL || type == DOUBLE_CHANNEL ) c = w;
		else c = w / 2;
		ratio = c / tf;

		/*flange*/
		if(built_up) {
			if(ratio < 9 * e) cl = CLASS1;
			else if(ratio < 10 * e) cl = CLASS2;
			else if(ratio < 14 * e) cl = CLASS3;
			else cl = CLASS4;
		} else {
			if(ratio < 10 * e) cl = CLASS1;
			else if(ratio < 11 * e) cl = CLASS2;
			else if(ratio < 15 * e) cl = CLASS3;
			else cl = CLASS4;
		}
		sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			DOUBLE ka,red,l;
			ka = GetKa(si[0]);
			l = (c / tf) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
            if(type == WIDEFLANGE || type == DOUBLE_CHANNEL) {
				dA += 4 * red * c * tf;
				dZ[0] += 4 * red * c * tf * pow(h / 2,2);
                dZ[1] += 4 * (pow(red * c,3) * tf / 3);
				if(type == DOUBLE_CHANNEL) dZ[1] += 4 * red * c * tf * pow(bd / 2, 2);
			} else  {
				dA += 2 * red * c * tf;
                if(type == CHANNEL) {
					dZ[0] += 2 * red * c * tf * pow(h / 2, 2);
                    dZ[1] += 2 * (pow(red * c, 3) * tf / 12);
				} else {
					DOUBLE y,hm,bm;
					hm = h + tf / 2;
					bm = w;
					y = hm - (pow(hm, 2) * tw + pow(tf,2) * (bm - tw)) / (2 * A);
					dZ[0] += 2 * red * c * tf * pow(y,2);
					dZ[1] += 2 * (pow(red * c, 3) * tf / 12);
				}
			}
		}

		/*web*/
		if(type == TEE) {
			ratio = (h + tf / 2) / tw;
			if(ratio < 9 * e) cl = CLASS1;
			else if(ratio < 10 * e) cl = CLASS2;
			else if(ratio < 16 * e) cl = CLASS3;
			else cl = CLASS4;
		} else {
			d = h - tf;
			ratio = d / tw;
			DOUBLE ac = 2 / (1 + si[0]);
			
			if(ratio < 66 * e / (0.4 + 0.6 * ac)) cl = CLASS1;
			else if(ratio < 75 * e / ac) cl = CLASS2;
			else  {
				DOUBLE p = fm[0] / material->fyk;
				if(p > 0) {
					if(built_up) {
						if(ratio <= 120 * e / (1 + 1.5 * p) && ratio <= (41 / p - 13) * e) cl = CLASS3;
						else cl = CLASS4;
					} else {
						if(ratio <= 120 * e / (1 + 1.5 * p) && ratio <= (41 / p - 2) * e) cl = CLASS3;
						else cl = CLASS4;
					}
				} else {
					if(ratio <= 120 * e / pow(1 + p,2) && ratio <= 250 * e) cl = CLASS3;
					else cl = CLASS4;
				}
			}

		}
		if(sclass < cl) sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			c = ratio * tw;
			DOUBLE ka,red,l;
			ka = GetKa(si[0]);
			l = (c / tw) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
            if(type == DOUBLE_CHANNEL) {
				dA += 2 * red * c * tw;
				dZ[0] += 2 * (pow(red * c, 3) * tw / 12);
			} else {
				dA += red * c * tw;
				dZ[0] += (pow(red * c, 3) * tw / 12);
				if(type == CHANNEL) {
					DOUBLE x,hm,bm;
					hm = h + tf;
					bm = w + tw / 2;
					x = (2 * bm * bm * tf + (hm - 2 * tf) * tw * tw) / (2 * A);
					eA[1] = (A * x) / (A - 2 * red * c * tf) - x;
				}
			}
		}
	} else if(type == ANGLE || type == DOUBLE_ANGLE) {
		t = (tw + tf) / 2;
		ratio = (w) / t;
		ratio1 = (h + tf/2) / t;
		
		if(ratio < 9 * e && ratio1 < 9 * e) cl = CLASS1;
		else if(ratio < 10 * e && ratio1 < 10 * e) cl = CLASS2;
		else if(ratio < 15 * e && ratio1 < 15 * e && ratio + ratio1 < 23 * e) cl = CLASS3;
		else cl = CLASS4;

		sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			c = (w + h + tf/2) / 2;
			if(c > h) c = h;

			DOUBLE ka,red,l;
			ka = (GetKa(si[0]) + GetKa(si[1])) / 2;
			l = (c / t) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
			dA += red * A;

			if(type == ANGLE) {
				DOUBLE x,y,hm,bm,l;
				hm = h + tf / 2;
				bm = w + tw / 2;
				x = bm - (tw * (2 * (bm - tw) + hm) + pow(bm - tw,2)) / (2 * (bm - tw + hm));
				y = hm - (tf * (2 * (hm - tf) + bm) + pow(hm - tf,2)) / (2 * (hm - tf + bm)); 
				l = (red * A * h / (h + w)) / tw;
				dZ[0] += tw * pow(l,3) / 12;
				l = (red * A * w / (h + w)) / tf;
                dZ[1] += tf * pow(l,3) / 12;
				eA[0] = (A * y) / (A - red * A) - y;
				eA[1] = (A * x) / (A - red * A) - x;
			} else {
				DOUBLE y,hm,bm;
				hm = h + tf / 2;
				bm = w;
				y = hm - (tf * (2 * (hm - tf) + bm) + pow(hm - tf,2)) / (2 * (hm - tf + bm)); 
				l = (red * A * h / (h + w)) / (2 * tw);
				dZ[0] += 2 * tw * pow(l,3) / 12;
				l = (red * A * w / (h + w)) / (tf);
				dZ[1] += tf * pow(l,3) / 12;
				dZ[1] += (red * A * h / (h + w)) * pow(bd / 2, 2);
			}
			
		}
	} else if(type == BOX) {
		b = w - 2 * tw; 
		ratio = b / tf;
		if(ratio < 28 * e) cl = CLASS1;
		else if(ratio < 35 * e) cl = CLASS2;
		else if(ratio < 42 * e) cl = CLASS3;
		else cl = CLASS4;
		sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			DOUBLE ka,red,l;
			ka = GetKa(si[0]);
			l = (b / tf) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
			dA += 2 * red * b * tf;
			dZ[0] += 2 * red * b * tf * pow(h / 2, 2);
			dZ[1] += 2 * (pow(red * b,3) * tf / 12);
		}

		b = h - 2 * tf; 
		ratio = b / tf;
		if(ratio < 28 * e) cl = CLASS1;
		else if(ratio < 35 * e) cl = CLASS2;
		else if(ratio < 42 * e) cl = CLASS3;
		else cl = CLASS4;
		if(sclass < cl) sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			DOUBLE ka,red,l;
			ka = GetKa(si[1]);
			l = (b / tw) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
			dA += 2 * red * b * tw;
			dZ[0] += 2 * (pow(red * b, 3) * tw / 12);
			dZ[1] += 2 * red * b * tw * pow(w / 2,2);
		}

	} else if(type == PIPE) {
		ratio = (2 * r + tw) / tw;
		if(ratio < 47 * e * e) cl = CLASS1;
		else if(ratio < 67 * e * e) cl = CLASS2;
		else if(ratio < 90 * e * e) cl = CLASS3;
		else cl = CLASS4;
		sclass = cl;

		/*reduce area for class 4*/
		if(cl == CLASS4) {
			DOUBLE ka,red,l;
			ka = (GetKa(si[0]) + GetKa(si[1])) / 2;
			l = ((2 * r) / tw) * 28.4 * e * sqrt(ka);
			if(l <= 0.673) red = 1;
			else red = (l - 0.22) / (l * l);
			dA += red * A;
			DOUBLE y;
			y = (red * A / (4 * tw)) / r;
			y = r * sin(l);
			l = (red * A / 4) / tw;
			dZ[0] += (red * A / 2) * (pow(y,2) + l * l / 12);
			dZ[1] += (red * A / 2) * (pow(y,2) + l * l / 12);
		}
	} 
}