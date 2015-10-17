#ifndef __DXF_H
#define __DXF_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

/*
DXF file export
*/
struct DXFSTYLE {
	CString name;
	CString font;
	int flags;
	int textflags;
	DOUBLE fixedheight;
	DOUBLE widthfactor;
	DOUBLE height;
	DOUBLE obliqueangle;
	DXFSTYLE() {
		name = "";
		font = "";
		flags = 0;
		textflags = 0;
		fixedheight = 0;
		widthfactor = 0;
		height = 0;
		obliqueangle = 0;
	}
};
struct DXFLTYPE {
	CString name;
	CString description;
	int		flags;
	int		ecount;		
	double	length;	
	double	elem[30];	
	DXFLTYPE() {
		name = "";
		description = "";
		flags = 0;
		ecount = 0;
		length = 0.0;
	}
};
struct DXFLAYER {
	CString name;
	int	flags;
	int color;
	DXFLTYPE* line;
	DXFLAYER() {
		name = "";
		flags = 0;
		color = 256;
		line = 0;
	}
};

class DXF {
public:
	CString fname;
	FILE* f;
	enum {HEADER,CLASSES,TABLES,BLOCKS,ENTITIES,OBJECTS,THUMBNAILIMAGE};
	enum {LTYPE,LAYER,STYLE,APPID};
	DXFSTYLE current_style;
	DXFLAYER current_layer;
public:
	DXF(CString);
	void init(bool read = false);
	void finish(bool read = false);
	void begin_section(int);
	void end_section();
	void begin_table(int);
	void end_table();
	void draw_line(RPoint,RPoint);
	void draw_circle(RPoint,DOUBLE);
	void draw_rectangle(RPoint,DOUBLE,DOUBLE);
	void draw_text(CString,RPoint,RPoint);
	void add_variable(CString); 
	void add_line(DXFLTYPE* line);
    void add_layer(DXFLAYER* layer);
	void add_style(DXFSTYLE* style);
	void set_layer(DXFLAYER* layer) {current_layer = *layer;}
	void set_style(DXFSTYLE* style) {current_style = *style;}
	bool read_codes(char*,char*);
	void read_lines(RPOINTLIST&);
};

#endif