#include "common.h"
#include "dxf.h"
#include "GUI.h"

void CmyDocument::ExportToAutocad(CString sPath) {
	DXF dxf(sPath);
	
	/*definitions*/
	DXFLTYPE solidline;
	solidline.name = "solid";
	solidline.description = "solid line ------";
	
	
	DXFLAYER layer1;
	layer1.name = "layer1";
	layer1.color = 255;
	layer1.flags = 0;
	layer1.line = &solidline;
	
	dxf.set_layer(&layer1);
	
	/*drawing file*/
	dxf.init();	
	if(!dxf.f) {
		AfxMessageBox("File is in use! Close that and try again.");
		return;
	}
	
	/*tables*/
	dxf.begin_section(DXF::TABLES);
	dxf.add_line(&solidline);
	dxf.add_layer(&layer1);
	dxf.end_section();
	
	/*entities*/
	dxf.begin_section(DXF::ENTITIES);
	
	POSITION pos;
	MEMBER* member;
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		dxf.draw_line(member->j1->p,member->j2->p);
	}
	
	dxf.end_section();
	
	dxf.finish();
	/*end*/
}
void CmyDocument::ImportFromAutocad(CString sPath) {
	RPOINTLIST list;
	DXF dxf(sPath);
	dxf.init(true);
	dxf.read_lines(list);
	
	SYSTEM system;
	system.name = "GLOBAL";
	system.rank = LEVEL3;
	JOINT j1,j2;
	POSITION pos = list.GetHeadPosition();
	while(pos) {
		j1.p = list.GetNext(pos);
        j2.p = list.GetNext(pos);
		system.grid[IUX].AddSorted(j1.p.x);
        system.grid[IUY].AddSorted(j1.p.y);
		system.grid[IUZ].AddSorted(j1.p.z);
		system.grid[IUX].AddSorted(j2.p.x);
        system.grid[IUY].AddSorted(j2.p.y);
		system.grid[IUZ].AddSorted(j2.p.z);
		AddMember(&j1,&j2);
	}
	*global = system;
	
	/*views*/
	CmyView* pView;
	pos = GetFirstViewPosition();
	while(pos) {
		pView = (CmyView*) GetNextView(pos);
		pView->force_diagram = INONE;
		if(pView->view != VIEW_3D) {
			GRIDLIST* pgrid = &global->grid[NORMGRID[pView->view]];
			if(!pgrid->IsEmpty()) {
				pView->position = pgrid->GetCount() - 1;
				pView->pvalue = pgrid->GetAt(pgrid->FindIndex(pView->position));
			} 
		}
		pView->InitView();
	}
	UpdateAllViews(NULL);
	
	/*close file*/
	dxf.finish(true);
}
struct CONNECT {
	UBMP8 flags;
	DOUBLE width[6];
	CONNECT() {
		flags = 0;
	}
};
void CmyDocument::OnPrepareDrawing(UINT nID) {
	/*open file*/
	CString sPath = GetPathName();
	if(sPath.IsEmpty()) {
		AfxMessageBox("Save file before continuing");
		return;
	}
	int len = sPath.GetLength();
	sPath.SetAt(len - 4,'_');
	sPath += ".dxf";
	DXF dxf(sPath);
	dxf.init();	
	if(!dxf.f) {
		AfxMessageBox("File is in use! Close that and try again.");
		return;
	}

	/*fill joint data*/
	POSITION pos,pos1;
	MEMBER* member;
	JOINT* joint,*joint1;
	RPoint p1,p2;
	GRIDINFO ginfo;
	GetGridInfo(&ginfo);
	CONNECT* connect;
	UINT design_selected = (nID == IDM_PREPARE_DRAWING_SELECTED);

	/*allocate joint connection data*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->ptr = new CONNECT;
	}
	
#define AddConnect(x) {\
	connect = (CONNECT*) joint->ptr;\
	connect->flags |= (1 << x);\
	if(member->section->type == RECTANGULAR) connect->width[x] = member->section->h;\
	else connect->width[x] = member->section->r;\
};

	/*form connection*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		pos1 = members.GetHeadPosition();
		while(pos1) {
			member = &members.GetNext(pos1);
			
			p1 = member->j1->p;
            p2 = member->j2->p;
			if(member->j1 == joint || member->j2 == joint) {
                if(member->j1 == joint) joint1 = member->j2;
				else joint1 = member->j1;
				
				if(EQUAL(p1.z,p2.z)) {
					if(EQUAL(p1.y,p2.y)) {
						if(joint->p.x < member->j2->p.x) {
							AddConnect(IUZ);
						} else {
							AddConnect(IRZ);
						}
					} else if(EQUAL(p1.x,p2.x)) {
						if(joint->p.y < member->j2->p.y) {
							AddConnect(IUZ);
						} else { 
							AddConnect(IRZ);
						}
					}
				} else if(EQUAL(p1.y,p2.y) && EQUAL(p1.x,p2.x)) {
					if(joint->p.z < member->j2->p.z) {
						AddConnect(IUX);
                        AddConnect(IUY);
					} else {
						AddConnect(IRX);
						AddConnect(IRY);
					}
				}
			}
		}
	}
	
	
	/*definitions*/
	DXFLTYPE solidline;
	solidline.name = "solid";
	solidline.description = "solid line ------";
	
	DXFLTYPE dotline;
	dotline.name = "DOTTED";
	dotline.description = "dotted line -.-.-.";
	dotline.ecount = 4;
	dotline.length = 0.5;
	dotline.elem[0] = 0.25;
	dotline.elem[1] = -0.125;
	dotline.elem[2] = 0.0;
	dotline.elem[3] = -0.125;
	
	DXFLAYER layer1;
	layer1.name = "layer1";
	layer1.color = 255;
	layer1.flags = 0;
	layer1.line = &solidline;
	
	DXFLAYER layer2;
	layer2.name = "layer2";
	layer2.color = 100;
	layer2.flags = 0;
	layer2.line = &solidline;
	
	
	DXFLAYER layer3;
	layer3.name = "layer3";
	layer3.color = 253;
	layer3.flags = 0;
	layer3.line = &dotline;
	
	DXFLAYER layer4;
	layer4.name = "layer4";
	layer4.color = 9;
	layer4.flags = 0;
	layer4.line = &solidline;
	
	DXFLAYER layer5;
	layer5.name = "layer5";
	layer5.color = 6;
	layer5.flags = 0;
	layer5.line = &solidline;
	
	DXFLAYER layer6;
	layer6.name = "layer6";
	layer6.color = 1;
	layer6.flags = 0;
	layer6.line = &solidline;
	
	
	DXFSTYLE style1;
	style1.name = "style1";
	style1.font = "TIMES.TIF";
	style1.height = 0.1;
	style1.widthfactor = 1;
	style1.flags = 0;
	style1.textflags = 0;
	
	DXFSTYLE style2;
	style2.name = "style2";
	style2.font = "TIMES.TIF";
	style2.height = 0.3;
	style2.widthfactor = 1;
	style2.flags = 0;
	style2.textflags = 0;
	
	DXFSTYLE style3;
	style3.name = "style3";
	style3.font = "TIMES.TIF";
	style3.height = 0.2;
	style3.widthfactor = 1;
	style3.flags = 0;
	style3.textflags = 0;
	
	dxf.set_style(&style1);
	dxf.set_layer(&layer1);
	
	/*header*/
	dxf.begin_section(DXF::HEADER);
	dxf.end_section();
	
	/*tables*/
	dxf.begin_section(DXF::TABLES);
	dxf.add_line(&solidline);
    dxf.add_line(&dotline);
	dxf.add_layer(&layer1);
    dxf.add_layer(&layer2);
	dxf.add_layer(&layer3);
	dxf.add_layer(&layer4);
	dxf.add_layer(&layer5);
	dxf.add_layer(&layer6);
	dxf.add_style(&style1);
    dxf.add_style(&style2);
    dxf.add_style(&style3);
	dxf.end_section();
	
	/*entities*/
	dxf.begin_section(DXF::ENTITIES);
	
#define FindGrid(index,item,data) {\
	index = 0;\
    POSITION mpos = (data).GetHeadPosition();\
	while(mpos) {\
	     if(EQUAL(item,(data).GetNext(mpos)))\
	         break;\
	     index++;\
	}\
};
	
#define DrawLine(x1,y1,CD,D) {\
	DOUBLE st,en;\
	connect = (CONNECT*) member->j1->ptr;\
	if(connect->flags & CD) {\
	   st = x1 + connect->width[D];\
	   if(D <= IUZ) rp1.set(st,y1 + (st - x1),0);\
	   else rp1.set(st,y1 - (st - x1),0);\
	   rp2.set(st,y1,0);\
	   dxf.draw_line(rp1,rp2);\
	} else st = x1;\
	connect = (CONNECT*) member->j2->ptr;\
    if(connect->flags & CD) {\
	   en = x1 + L - connect->width[D];\
	   if(D <= IUZ) rp1.set(en,y1 + (L + x1 - en),0);\
	   else rp1.set(en,y1 - (L + x1 - en),0);\
	   rp2.set(en,y1,0);\
	   dxf.draw_line(rp1,rp2);\
	} else en = x1 + L;\
	rp1.set(st,y1,0);\
    rp2.set(en,y1,0);\
	dxf.draw_line(rp1,rp2);\
};
	
#define DrawAxisLabel(I,v1,v2) {\
	FindGrid(index,v1,global->grid[I]);\
	if(index != maxn[I]) {\
	   if(index < 10) str1.Format(" %c%d",'X' + I,index);\
	   else str1.Format("c%d",'X' + I,index);\
	}\
	FindGrid(index,v2,global->grid[I]);\
	if(index != maxn[I]) {\
	   if(index < 10) str2.Format(" %c%d",'X' + I,index);\
	   else str2.Format("c%d",'X' + I,index);\
	}\
};
	
    CString str,str1,str2;
	int i,j,index,cdir,inc;
	DOUBLE w,h,r,c,L,x,y,spacing,xd,scale,barsize,cp;
	int nz,ny,nt;
	RPoint rp1,rp2,center, end;
	int maxn[3];
	
	for(i = IUX;i <= IUZ;i++)
		maxn[i] = global->grid[i].GetCount();
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		if(design_selected && !member->sel) {
			continue;
		}
		if(member->design_data[0].flags)
			continue;
		
		w = member->section->w;
		h = member->section->h;
		r = member->section->r;
		c = member->section->rebar.cover;
		barsize = member->section->rebar.barsize;
		nz = member->section->rebar.nz;
		ny = member->section->rebar.ny;
		if(member->section->type == CIRCULAR) nt = member->section->rebar.nt;
		else nt = 2 * (nz + ny - 2);
		p1 = member->j1->p;
        p2 = member->j2->p;
		L = member->GetLength();
		/*draw axis*/
		str1 = "";
		str2 = "";
		cdir = -1;
		if(EQUAL(p1.z,p2.z)) {
			if(EQUAL(p1.y,p2.y)) {
				x = p1.x;
				y = 5 * p1.z + p1.y * (ginfo.ymax - ginfo.ymin);
				DrawAxisLabel(IUX,p1.x,p2.x);
				cdir = IUX;
			} else if(EQUAL(p1.x,p2.x)) {
				x = (ginfo.xmax - ginfo.xmin) + 5 + p1.y;
				y = 5 * p1.z + p1.x * (ginfo.xmax - ginfo.xmin);
				DrawAxisLabel(IUY,p1.y,p2.y);
				cdir = IUY;
			}
		} else if(EQUAL(p1.y,p2.y)) {
			if(EQUAL(p1.x,p2.x)) {
				x = (ginfo.xmax - ginfo.xmin + ginfo.ymax - ginfo.ymin) + 5 + p1.z;
				y = 5 * p1.y + p1.x * (ginfo.xmax - ginfo.xmin);
				DrawAxisLabel(IUZ,p1.z,p2.z);
				cdir = IUZ;
			} 
		}
		if(cdir == -1) {
			print("Can't Draw member %s\n",member->name);
			continue;
		}
		
		dxf.set_layer(&layer5);
		rp1.set(x,y - 1,0);
		rp2.set(x,y + 1.5,0);
		dxf.draw_line(rp1,rp2);
		rp1.set(x,y + 1.5 + 0.4,0);
		dxf.draw_circle(rp1,0.4);
		
		dxf.set_style(&style2);
		rp1.x -= 0.35;
		rp1.y -= 0.15;
		rp2 = rp1 + RPoint(0.4,0.3,0);
		dxf.draw_text(str1,rp1,rp2);
		
		rp1.set(x + L,y - 1,0);
		rp2.set(x + L,y + 1.5,0);
		dxf.draw_line(rp1,rp2);
		rp1.set(x + L,y + 1.5 + 0.4,0);
		dxf.draw_circle(rp1,0.4);
		
		dxf.set_style(&style2);
		rp1.x -= 0.35;
		rp1.y -= 0.15;
		rp2 = rp1 + RPoint(0.4,0.3,0);
		dxf.draw_text(str2,rp1,rp2);
		
		/*draw border*/
		dxf.set_layer(&layer1);
		rp1.set(x,y - h/2 + c,0);
        rp2.set(x + L,y - h/2 + c,0);
		dxf.draw_line(rp1,rp2);
		rp1.set(x,y + h/2 - c,0);
        rp2.set(x + L,y + h/2 - c,0);
		dxf.draw_line(rp1,rp2);
		
		dxf.set_layer(&layer2);
		if(cdir == IUX) {
			DrawLine(x, y - h/2,RX,IRX);
			DrawLine(x, y + h/2,UX,IUX); 
		} else if(cdir == IUY) {
			DrawLine(x, y - h/2,RY,IRY);
			DrawLine(x, y + h/2,UZ,IUZ); 
		} else {
			DrawLine(x, y - h/2,RZ,IRZ);
			DrawLine(x, y + h/2,UZ,IUZ); 
		}
		
		dxf.set_layer(&layer4);
		rp1.set(x,y + 0.9,0);
        rp2.set(x + L,y + 0.9,0);
		dxf.draw_line(rp1,rp2);
		
        dxf.set_style(&style1);
        str.Format("%.2fm",L);
		rp1.x = (rp1.x + rp2.x) * 0.5;
		rp1.y = rp2.y = rp1.y + 0.05;
		dxf.draw_text(str,rp1,rp2);
		
		/*shear reinforcement*/
		dxf.set_layer(&layer6);
		spacing = member->design_data[0].s / 1000.0f;
		for(i = 0;i < 7;i ++) {
			if(spacing * i >= L / 5) break;
			rp1.set(x + spacing * i,y - h/2 + c,0);
			rp2.set(x + spacing * i,y + h/2 - c,0);
			dxf.draw_line(rp1,rp2);
		}
		spacing = member->design_data[member->nDiv / 2].s / 1000.0f;
		for(i = 0;i < 3;i ++) {
			if(spacing * i >= L / 8) break;
			rp1.set(x + L / 2 + spacing * i,y - h/2 + c,0);
			rp2.set(x + L / 2 + spacing * i,y + h/2 - c,0);
			dxf.draw_line(rp1,rp2);
			rp1.set(x + L / 2 - spacing * i,y - h/2 + c,0);
			rp2.set(x + L / 2 - spacing * i,y + h/2 - c,0);
			dxf.draw_line(rp1,rp2);
		}
		spacing = member->design_data[member->nDiv - 1].s / 1000.0f;
		for(i = 0;i < 7;i ++) {
			if(spacing * i >= L / 5) break;
			rp1.set(x + L - spacing * i,y - h/2 + c,0);
			rp2.set(x + L - spacing * i,y + h/2 - c,0);
			dxf.draw_line(rp1,rp2);
		}
        /*shear rebar spacing text*/
		dxf.set_layer(&layer1);
		if(detailing.single_shear) {
			dxf.set_style(&style3);
			str.Format("%s%dc/c%d","%%C",int(member->design_data[0].sd * 1000),member->design_data[0].s);
			rp1.x = x + L / 2 - 0.45;
			rp2.x = rp1.x + 0.9;
			rp1.y = rp2.y = y - 0.1;
			dxf.draw_text(str,rp1,rp2);
		} else {
			dxf.set_style(&style1);
			str.Format("%s%dc/c%d","%%C",int(member->design_data[0].sd * 1000),member->design_data[0].s);
			rp1.x = x + L / 20;
			rp2.x = rp1.x + 0.5;
			rp1.y = rp2.y = y - 0.07;
			dxf.draw_text(str,rp1,rp2);
			str.Format("%s%dc/c%d","%%C",int(member->design_data[member->nDiv / 2].sd * 1000),member->design_data[member->nDiv / 2].s);
			rp1.x = x + L / 2 - 0.25;
			rp2.x = rp1.x + 0.5;
			rp1.y = rp2.y = y - 0.07;
			dxf.draw_text(str,rp1,rp2);
			str.Format("%s%dc/c%d","%%C",int(member->design_data[member->nDiv - 1].sd * 1000),member->design_data[member->nDiv - 1].s);
			rp1.x = x + L - 0.5 - L / 20;
			rp2.x = rp1.x + 0.5;
			rp1.y = rp2.y = y - 0.07;
			dxf.draw_text(str,rp1,rp2);
		}
		/*section*/
		if(member->section->rebar.design == BEAM) inc = 1;
		else inc = 3;
		const DOUBLE DSECTION = 3.0;
		
		dxf.set_style(&style3);
		for(int xx = 0;xx < 3;xx += inc) {
			if(xx == 0) xd = x + L / 2;
			else if(xx == 1) xd = x + 0.3 + L / 20;
			else xd = x + L - L / 20 - 0.3;
			
			dxf.set_layer(&layer6);
			rp1.set(xd,y - h/2,0);
			rp2.set(xd,y - h,0);
			dxf.draw_line(rp1,rp2);
			rp1.set(xd,y + h/2,0);
			rp2.set(xd,y + h,0);
			dxf.draw_line(rp1,rp2);
			
			
			dxf.set_layer(&layer1);
			str.Format("%s%d",str1,xx + 1);
			rp1.set(xd - 0.3,y - h - 0.20,0);
			rp2.set(xd + 0.3,y - h - 0.20,0);
			dxf.draw_text(str,rp1,rp2);
			rp1.set(xd - 0.3,y + h + 0.05,0);
			rp2.set(xd + 0.3,y + h + 0.05,0);
			dxf.draw_text(str,rp1,rp2);
			
			scale = 3;
			if(member->section->type == RECTANGULAR) {
				dxf.set_layer(&layer2);
				rp1.set(xd - (scale * w) / 2,y - DSECTION,0);
				dxf.draw_rectangle(rp1,scale * w,scale * h);
                dxf.set_layer(&layer1);
				rp1.x += scale * c;
				rp1.y -= scale * c;
				dxf.draw_rectangle(rp1,scale * (w - 2 * c) ,scale * (h - 2 * c));
				
				/*rebar*/
				if(member->section->rebar.design == COULMN) {
					center.set(xd,y - DSECTION - scale * h / 2,0);
					cp = c + barsize / 2;
					
					end.set(xd ,y - DSECTION - (scale * h * 1.5),0);
					
					for(i = 0;i < nz;i++) {
						dxf.set_layer(&layer1);
						rp1.x = center.x - (scale * (w / 2 - i * (w - 2 * cp) / (nz - 1) - cp));
						rp1.y = center.y - (scale * (h / 2 - cp));
						dxf.draw_circle(rp1,scale * barsize / 2);
						dxf.set_layer(&layer4);
						dxf.draw_line(rp1,end);
						
						dxf.set_layer(&layer1);
						rp1.x = center.x - (scale * (w / 2 - i * (w - 2 * cp) / (nz - 1) - cp));
						rp1.y = center.y - (scale * ( -h / 2 + cp));
						dxf.draw_circle(rp1,scale * barsize / 2);
						dxf.set_layer(&layer4);
						dxf.draw_line(rp1,end);
					}
					for(i = 1;i < ny - 1;i++) {
						dxf.set_layer(&layer1);
						rp1.x = center.x - (scale * (w / 2 - cp));
						rp1.y = center.y - (scale * (h / 2 - i * (h - 2 * cp) / (ny - 1) - cp));
						dxf.draw_circle(rp1,scale * barsize / 2);
						dxf.set_layer(&layer4);
						dxf.draw_line(rp1,end);
						
						dxf.set_layer(&layer1);
						rp1.x = center.x - (scale * (-w / 2 + cp));
						rp1.y = center.y - (scale * (h / 2 - i * (h - 2 * cp) / (ny - 1) - cp));
						dxf.draw_circle(rp1,scale * barsize / 2);
						dxf.set_layer(&layer4);
						dxf.draw_line(rp1,end);
					}
				} else {
					DESIGNDATA* pdata;
					if(xx == 0) pdata = &member->design_data[member->nDiv / 2];
					else if(xx == 1)  pdata = &member->design_data[0];
					else pdata = &member->design_data[member->nDiv - 1];
					
					cp = c + barsize / 2;
					for(j = 0;j < 2;j++) {
						
						/*top/bottom reinforcement*/
						if(j == 0) {
							center.set(xd,y - DSECTION - scale * cp,0);
							end.set(xd ,y - DSECTION + 0.8,0);
						} else {
							center.set(xd,y - DSECTION - scale * (h - cp),0);
							end.set(xd ,y - DSECTION - (scale * h * 1.5),0);
						}
						if(pdata->total[j]) {
							if(pdata->total[j] > 1)
								nz = pdata->count[j][0] + pdata->count[j][1];
							else 
								nz = pdata->count[j][0];
							barsize = pdata->diam[j][0];
							for(i = 0;i < nz;i++) {
								dxf.set_layer(&layer1);
								rp1.x = center.x - (scale * (w / 2 - i * (w - 2 * cp) / (nz - 1) - cp));
								rp1.y = center.y;
								dxf.draw_circle(rp1,scale * barsize / 2);
								dxf.set_layer(&layer4);
								dxf.draw_line(rp1,end);
							}
							
							/*number of bars*/
							rp1 = end;
							if(j) rp1.y -= 0.1;
							rp1.x -= 0.4;
							rp2 = rp1;
							rp2.x += 0.5;
							dxf.set_layer(&layer1);
							dxf.set_style(&style1);
							if(pdata->total[j] > 1) {
								str.Format("%d%s%d and %d%s%d",
									pdata->count[j][0],"%%C",int(pdata->diam[j][0] * 1000),
									pdata->count[j][1],"%%C",int(pdata->diam[j][1] * 1000)
									);
							} else {
								str.Format("%d%s%d",pdata->count[j][0],"%%C",int(pdata->diam[j][0] * 1000));
							}
							dxf.draw_text(str,rp1,rp2);
						}
					}
				}
			} else {
                dxf.set_layer(&layer2);
				rp1.set(xd,y - DSECTION - scale * r,0);
				dxf.draw_circle(rp1,scale * r);
				dxf.set_layer(&layer1);
				dxf.draw_circle(rp1,scale * (r - c));
				end.set(xd - (scale * r) + 0.1,y - DSECTION - (scale * r * 3),0);
				/*rebar*/
				if(member->section->rebar.design == COULMN) {
					center.set(xd,y - DSECTION - scale * r,0);
					cp = c + barsize / 2;
					DOUBLE q = 2 * PI / nt;
					for(i = 0;i < nt;i++) {
						dxf.set_layer(&layer1);
						rp1.x = center.x - scale * (r - cp) * sin(q * i);
						rp1.y = center.y - scale * (r - cp) * cos(q * i);
						dxf.draw_circle(rp1,scale * barsize / 2);
						dxf.set_layer(&layer4);
						dxf.draw_line(rp1,end);
					}
				}
			}
			/*total reinforement*/
			if(member->section->rebar.design == COULMN) {
                DESIGNDATA* pdata = &member->design_data[0];
				rp1 = end;
				rp1.y -= 0.1;
				rp2 = rp1;
				rp2.x += 0.5;
				dxf.set_layer(&layer1);
				dxf.set_style(&style1);
				if(pdata->total[0] > 1) {
					str.Format("%d X (%d%s%d and %d%s%d)",nt,
						pdata->count[0][0],"%%C",int(pdata->diam[0][0] * 1000),
						pdata->count[0][1],"%%C",int(pdata->diam[0][1] * 1000)
						);
				} else {
					str.Format("%d%s%d",
						nt * pdata->count[0][0],"%%C",
						int(pdata->diam[0][0] * 1000));
				}
				dxf.draw_text(str,rp1,rp2);
			}
			/*section text*/
			dxf.set_layer(&layer1);
			dxf.set_style(&style1);
			rp1 = end;
			if(member->section->rebar.design != COULMN) rp1.x -= 0.4;
			rp1.y -= 0.4;
			rp2 = rp1;
			rp2.x += 1;
			str.Format("Section %s%d",str1,xx + 1);
			dxf.draw_text(str,rp1,rp2);
			
			/*shear bars*/
			if(!detailing.single_shear || xx == 0) {
				const DOUBLE DSHEAR = 10.0;
				DOUBLE length;
				if(member->section->type == RECTANGULAR) {
					length = 2 * w + 2 * h;
					rp1.set(xd - (scale * w) / 2,y - DSHEAR,0);
					rp2.set(xd + (scale * w) / 2,y - DSHEAR,0);
					dxf.draw_line(rp1,rp2);
					rp1.set(xd - (scale * w) / 2,y - DSHEAR,0);
					rp2.set(xd - (scale * w) / 2,y - DSHEAR - scale * h,0);
					dxf.draw_line(rp1,rp2);
					rp1.set(xd - (scale * w) / 2,y - DSHEAR - scale * h,0);
					rp2.set(xd + (scale * w) / 2,y - DSHEAR - scale * h,0);
					dxf.draw_line(rp1,rp2);
					rp1.set(xd + (scale * w) / 2,y - DSHEAR - scale * h,0);
					rp2.set(xd + (scale * w * 3) / 4,y - DSHEAR - scale * h / 4,0);
					dxf.draw_line(rp1,rp2);
					rp1.set(xd + (scale * w) / 2,y - DSHEAR,0);
					rp2.set(xd + (scale * w) / 4,y - DSHEAR - scale * h / 8,0);
					dxf.draw_line(rp1,rp2);
					rp1.set(xd + (scale * w * 3) / 4,y - DSHEAR - scale * h / 4,0);
					rp2.set(xd + (scale * w) / 2,y - DSHEAR - scale * h / 4,0);
					dxf.draw_line(rp1,rp2); 
				} else {
                    length = 2 * PI * r;
					rp1.set(xd,y - DSHEAR - scale * r,0);
					dxf.draw_circle(rp1,scale * r);
					rp1.set(xd + scale * r,y - DSHEAR - scale * r,0);
                    rp2.set(xd + scale * r * 4 / 5,y - DSHEAR - scale * r + 0.1,0);
					dxf.draw_line(rp1,rp2);
					rp2.set(xd + scale * r * 4 / 5,y - DSHEAR - scale * r - 0.1,0);
					dxf.draw_line(rp1,rp2);
				}
				/*text*/
				length += 0.1;
				str.Format("%s%d L = %d","%%C",int(member->design_data[0].sd * 1000),int(length * 1000));
				rp1.set(xd - 0.4,y - DSHEAR + 0.1,0);
				rp2 = rp1;
				rp2.x += 0.5;
				dxf.draw_text(str,rp1,rp2);
			}
		}
		
		/*longituidinal bars*/
		const DOUBLE DBAR = 7.0;
		if(member->section->rebar.design == COULMN) {
			DESIGNDATA* pdata = &member->design_data[0];
			rp1.set(x,y - DBAR,0);
			rp2.set(x + L,y - DBAR,0);
			dxf.draw_line(rp1,rp2);
			
			if(pdata->total[0] > 1) {
				str.Format("%d X (%d%s%d and %d%s%d)",nt,
					pdata->count[0][0],"%%C",int(pdata->diam[0][0] * 1000),
					pdata->count[0][1],"%%C",int(pdata->diam[0][1] * 1000)
					);
			} else {
				str.Format("%d%s%d",
					nt * pdata->count[0][0],"%%C",
					int(pdata->diam[0][0] * 1000));
			}
			int length = int(L * 1000) + 800;
			str.Format("%s L = %d",str,length);
			rp1.set(x + L / 2 - 0.6,y - DBAR + 0.1,0);
			rp2 = rp1;
			rp2.x += 0.5;
			dxf.draw_text(str,rp1,rp2);
		}
		if(detailing.cut_at_L3 && member->section->rebar.design == BEAM) {
			DESIGNDATA* pdata;
			UINT minb,maxb;
			
			minb =  100;
			maxb = 0;
			for(UINT i = 0;i < member->nDiv;i++) {
				if(member->design_data[i].count[0][0] < minb)
					minb = member->design_data[i].count[0][0];
                if(member->design_data[i].count[1][0] > maxb)
					maxb = member->design_data[i].count[1][0];
			}
			
			pdata = &member->design_data[0];
			dxf.set_style(&style1);
			dxf.set_layer(&layer1);
			
			/*common top rebar*/
			rp1.set(x,y - DBAR,0);
			rp2.set(x + L,y - DBAR,0);
			dxf.draw_line(rp1,rp2);
			
			str.Format("%d%s%d L = %d(T)",minb,"%%C",int(pdata->diam[0][0] * 1000),int(L * 1000));
			rp1.set(x + L / 2 - 0.6,y - DBAR + 0.1,0);
			rp2 = rp1;
			rp2.x += 0.5;
			dxf.draw_text(str,rp1,rp2);
			
			/*top 1 rebar*/
			pdata = &member->design_data[0];
			if(pdata->count[0][1] || pdata->count[0][0] - minb > 0) { 
				rp1.set(x,y - DBAR - 1,0);
				rp2.set(x + L / 3,y - DBAR - 1,0);
				dxf.draw_line(rp1,rp2);
				
				if(pdata->count[0][1]) 
					str.Format("%d%s%d L = %d(T)",pdata->count[0][1],"%%C",int(pdata->diam[0][1] * 1000),int((L / 3) * 1000));
				else
					str.Format("%d%s%d L = %d(T)",pdata->count[0][0] - minb,"%%C",int(pdata->diam[0][0] * 1000),int((L / 3) * 1000));
				rp1.set(x + 0.8 ,y - DBAR - 1 + 0.1,0);
				rp2 = rp1;
				rp2.x += 0.5;
				dxf.draw_text(str,rp1,rp2);
			}
			
			/*top 2 rebar*/
			pdata = &member->design_data[member->nDiv - 1];
			if(pdata->count[0][1] || pdata->count[0][0] - minb > 0) { 
				rp1.set(x + 2 * L / 3,y - DBAR - 1,0);
				rp2.set(x + L,y - DBAR - 1,0);
				dxf.draw_line(rp1,rp2);
				
				if(pdata->count[0][1]) 
					str.Format("%d%s%d L = %d(T)",pdata->count[0][1],"%%C",int(pdata->diam[0][1] * 1000),int((L / 3) * 1000));
				else
					str.Format("%d%s%d L = %d(T)",pdata->count[0][0] - minb,"%%C",int(pdata->diam[0][0] * 1000),int((L / 3) * 1000));
				rp1.set(x + 2 * L / 3 ,y - DBAR - 1 + 0.1,0);
				rp2 = rp1;
				rp2.x += 0.5;
				dxf.draw_text(str,rp1,rp2);
			}
			
			/*bottom rebar*/
			rp1.set(x,y - DBAR - 2,0);
			rp2.set(x + L,y - DBAR - 2,0);
			dxf.draw_line(rp1,rp2);
            
			str.Format("%d%s%d L = %d(B)",maxb,"%%C",int(pdata->diam[1][0] * 1000),int(L * 1000));
			rp1.set(x + L / 2 - 0.6,y - DBAR - 2 + 0.1,0);
			rp2 = rp1;
			rp2.x += 0.5;
			dxf.draw_text(str,rp1,rp2);
		}
	}
	
#undef DrawLine
#undef FindGrid
#undef DrawAxisLabel
    
	dxf.end_section();
	
    /*end*/
	dxf.finish();

	/*free joint connection data*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		delete joint->ptr;
	}
}