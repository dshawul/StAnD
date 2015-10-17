#include "common.h"
#include "dxf.h"

DXF::DXF(CString name) {
	fname = name;
}
void DXF::init(bool read) {
	if(read) f = fopen(fname,"r");
	else f = fopen(fname,"w");
}
void DXF::finish(bool read) {
	if(!read) fprintf(f,"0\nEOF\n");
	fclose(f);
}
void DXF::begin_section(int section) {
	CString str;
	switch(section) {
	case HEADER:str = "HEADER";break;
	case CLASSES:str = "CLASSES";break;
	case TABLES:str = "TABLES";break;
	case BLOCKS:str = "BLOCKS";break;
	case ENTITIES:str = "ENTITIES";break;
	case OBJECTS:str = "OBJECTS";break;
	case THUMBNAILIMAGE:str = "THUMBNAILIMAGE";break;
	}
	fprintf(f,"0\nSECTION\n2\n%s\n",str);
}
void DXF::end_section() {
	fprintf(f,"0\nENDSEC\n");
}
void DXF::begin_table(int table) {
	CString str;
	switch(table) {
	case LTYPE:str = "LTYPE";break;
	case LAYER:str = "LAYER";break;
	case STYLE:str = "STYLE";break;
	case APPID:str = "APPID";break;
	}
	fprintf(f,"0\nTABLE\n2\n%s\n0\n%s\n",str,str);
}
void DXF::end_table() {
	fprintf(f,"0\nENDTAB\n");
}
void DXF::draw_circle(RPoint c,DOUBLE r) {
	fprintf(f,"0\nCIRCLE\n8\n%s\n62\n%d\n6\n%s\n10\n%f\n20\n%f\n30\n%f\n40\n%f\n",
		current_layer.name,current_layer.color,current_layer.line->name,c.x,c.y,c.z,r);
}
void DXF::draw_line(RPoint p1,RPoint p2) {
	fprintf(f,"0\nLINE\n8\n%s\n62\n%d\n6\n%s\n39\n%d\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n",
		current_layer.name,current_layer.color,current_layer.line->name,0,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
}
void DXF::draw_rectangle(RPoint p1,DOUBLE w,DOUBLE h) {
	RPoint p2,p3;
	p2 = p1;
	p2.x += w;
	draw_line(p1,p2);
	p3 = p2;
	p3.y -= h;
	draw_line(p2,p3);
	p2 = p3;
	p2.x -= w;
    draw_line(p3,p2);
    draw_line(p2,p1);
}
void DXF::draw_text(CString str,RPoint p1,RPoint p2) {
	fprintf(f,"0\nTEXT\n8\n%s\n62\n%d\n6\n%s\n7\n%s\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n40\n%f\n41\n%f\n50\n%f\n72\n%d\n73\n%d\n1\n%s\n",
        current_layer.name,current_layer.color,current_layer.line->name,current_style.name,
		p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,
		current_style.height,current_style.widthfactor,current_style.obliqueangle,
		(current_style.textflags & 0xff),((current_style.textflags >> 8) & 0xff),str);
}
void DXF::add_line(DXFLTYPE* line) {
	begin_table(LTYPE);
    fprintf(f,"2\n%s\n70\n%d\n3\n%s\n72\n65\n73\n%d\n40\n%f\n",
		line->name,line->flags,line->description,line->ecount,line->length);
	for(int i = 0;i < line->ecount;i++) {
		fprintf(f,"49\n%f\n",line->elem[i]);
	}
	end_table();
}
void DXF::add_layer(DXFLAYER* layer) {
	begin_table(LAYER);
    fprintf(f,"2\n%s\n70\n%d\n62\n%d\n6\n%s\n",
		   layer->name,layer->flags,layer->color,layer->line->name);
	end_table();
}
void DXF::add_style(DXFSTYLE* style) {
	begin_table(STYLE);
    fprintf(f,"2\n%s\n3\n%s\n70\n%d\n71\n%d\n40\n%f\n41\n%f\n42\n%f\n50\n%f\n",
		   style->name,style->font,style->flags,style->textflags,style->fixedheight,
		   style->widthfactor,style->height,style->obliqueangle);
	end_table();
}
void DXF::add_variable(CString name) {
	fprintf(f,"9\n$%s\n",name);
}

#define MAX_STR            256

bool DXF::read_codes(char* code,char* value) {
	if(!fgets(code,MAX_STR,f)) return false;
	if(!fgets(value,MAX_STR,f)) return false;
	value[strlen(value) - 1] = 0;
	code[strlen(code) - 1] = 0;
	return true;
}

void DXF::read_lines(RPOINTLIST& list) {
	RPoint rp1,rp2;
	char   code[MAX_STR];
    char   value[MAX_STR];
	int val;
	while(true) {
		if(!read_codes(code,value)) return;
		val = atoi(code);
		if(val == 0 && !strcmp(value,"SECTION")) {
			if(!read_codes(code,value)) return;
			
			if(!strcmp(value,"ENTITIES")) {
				while(true) {
					if(!read_codes(code,value)) return;
					if(!strcmp(value,"ENDSEC")) return;
					
					if(!strcmp(value,"LINE")) {
						while(true) {
							if(!read_codes(code,value)) return;
							int val = atoi(code);
							switch(val) {
							case 10: rp1.x = atof(value);break;
							case 20: rp1.y = atof(value);break;
							case 30: rp1.z = atof(value);break;
							case 11: rp2.x = atof(value);break;
							case 21: rp2.y = atof(value);break;
							case 31: rp2.z = atof(value);break;
							}
							if(val == 0) {
								list.AddTail(rp1);
								list.AddTail(rp2);
								if(strcmp(value,"LINE")) 
									break;
							}
						}
					}
				}
			}
		}
	}
}