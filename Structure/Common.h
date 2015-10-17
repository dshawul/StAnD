#ifndef __COMMON_H
#define __COMMON_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

#include <math.h>
#include "matrix.h"
#include "rpoint.h"
#include "utils.h"
#include "stdafx.h"


#undef pow
#define pow(x,y) pow((double)x,(double)y)
#undef log
#define log(x) log((double)x)
#undef sqrt
#define sqrt(x) sqrt((double)x)

typedef unsigned char UBMP8;
typedef unsigned int  UBMP32;

/*
Data types
*/
#define DATALIST            CList<DATA,DATA&>
#define PDATALIST           CList<DATA,DATA&>*

#define LOADCOMBOLIST       CList<LOADCOMBO,LOADCOMBO&>
#define PLOADCOMBOLIST      CList<LOADCOMBO,LOADCOMBO&>*
#define COMBOLIST           CList<COMBINATION,COMBINATION&>
#define PCOMBOLIST          CList<COMBINATION,COMBINATION&>*
#define NLCASELIST          CList<NLCASE,NLCASE&>
#define PNLCASELIST         CList<NLCASE,NLCASE&>*
#define BUCKLINGCASELIST    CList<BUCKLINGCASE,BUCKLINGCASE&>
#define PBUCKLINGCASELIST   CList<BUCKLINGCASE,BUCKLINGCASE&>*
#define LOADCASELIST        CList<LOADCASE,LOADCASE&>
#define PLOADCASELIST       CList<LOADCASE,LOADCASE&>*
#define MODALCASELIST       CList<MODALCASE,MODALCASE&>
#define PMODALCASELIST      CList<MODALCASE,MODALCASE&>*
#define RESPONSECASELIST    CList<RESPONSEHIST,RESPONSEHIST&>
#define PRESPONSECASELIST   CList<RESPONSEHIST,RESPONSEHIST&>*
#define RESPONSESPECLIST    CList<RESPONSESPEC,RESPONSESPEC&>
#define PRESPONSESPECLIST   CList<RESPONSESPEC,RESPONSESPEC&>*
#define ANALYSISCASELIST    CList<ANALYSISCASE,ANALYSISCASE&>
#define PANALYSISCASELIST   CList<ANALYSISCASE,ANALYSISCASE&>*

#define BLOADLIST           CList<BASE_LOAD,BASE_LOAD&>
#define PBLOADLIST          CList<BASE_LOAD,BASE_LOAD&>*
#define LOADLIST            CList<LOAD,LOAD&>
#define PLOADLIST           CList<LOAD,LOAD&>*
#define JLOADLIST           CList<JLOAD,JLOAD&>
#define PJLOADLIST          CList<JLOAD,JLOAD&>*

#define ENTITYLIST          CList<ENTITY,ENTITY&>
#define PENTITYLIST         CList<ENTITY,ENTITY&>*
#define ENTITYPLIST         CList<ENTITY*,ENTITY*&>
#define ENTITYPLIST_PLIST   CList<ENTITYPLIST*,ENTITYPLIST*&>
#define ELEMENTLIST         CList<ELEMENT,ELEMENT&>
#define PELEMENTLIST        CList<ELEMENT,ELEMENT&>*
#define ELEMENTPLIST        CList<ELEMENT*,ELEMENT*&>
#define ELEMENTPLIST_PLIST  CList<ELEMENTPLIST*,ELEMENTPLIST*&>
#define JOINTLIST           CList<JOINT,JOINT&>
#define PJOINTLIST          CList<JOINT,JOINT&>*
#define JOINTPLIST          CList<JOINT*,JOINT*&>
#define JOINTPLIST_PLIST    CList<JOINTPLIST*,JOINTPLIST*&>
#define MEMBERLIST          CList<MEMBER,MEMBER&>
#define PMEMBERLIST         CList<MEMBER,MEMBER&>*
#define MEMBERPLIST         CList<MEMBER*,MEMBER*&>
#define MEMBERPLIST_PLIST   CList<MEMBERPLIST*,MEMBERPLIST*&>
#define SLABLIST            CList<SLAB,SLAB&>
#define PSLABLIST           CList<SLAB,SLAB&>*
#define SLABPLIST           CList<SLAB*,SLAB*&>
#define SLABPLIST_PLIST     CList<SLABPLIST*,SLABPLIST*&>

#define MATERIALLIST        CList<MATERIAL,MATERIAL&>
#define PMATERIALLIST       CList<MATERIAL,MATERIAL&>*
#define SECTIONLIST         CList<SECTION,SECTION&>
#define PSECTIONLIST        CList<SECTION,SECTION&>*
#define ASECTIONLIST        CList<ASECTION,ASECTION&>
#define PASECTIONLIST       CList<ASECTION,ASECTION&>*
#define GRIDLIST            CSortedList<DOUBLE,DOUBLE&>
#define PGRIDLIST           CSortedList<DOUBLE,DOUBLE&>*
#define SYSTEMLIST          CList<SYSTEM,SYSTEM&>
#define PSYSTEMLIST         CList<SYSTEM,SYSTEM&>*
#define CONSTRAINTLIST      CList<CONSTRAINT,CONSTRAINT&>
#define PCONSTRAINTLIST     CList<CONSTRAINT,CONSTRAINT&>*
#define RPOINTLIST          CList<RPoint,RPoint&>
#define PRPOINTLIST         CList<RPoint,RPoint&>*
#define SRPOINTLIST         CSortedList<RPoint,RPoint&>
#define PSRPOINTLIST        CSortedList<RPoint,RPoint&>*
#define FUNCTIONLIST        CList<FUNCTION,FUNCTION&>
#define PFUNCTIONLIST       CList<FUNCTION,FUNCTION&>*
#define SPECFUNCLIST        CList<SPECFUNC,SPECFUNC&>
#define PSPECFUNCLIST       CList<SPECFUNC,SPECFUNC&>*
#define GROUPLIST           CList<GROUP,GROUP&>
#define PGROUPLIST          CList<GROUP,GROUP&>*

/*
Definitions
*/
enum {
	NOUR  =  0,
	UX    =  1,  UY =  2,  UZ =  4,
	RX    =  8,  RY = 16,  RZ = 32,  
    ALLU  =  7,
	ALLR  = 56,
	ALLUR = 63
};
enum {
    IUX,IUY,IUZ,IRX,IRY,IRZ,INONE,IDEFLECTED,IREACTIONS,IDISPLACEMENT,
	ISX,ISY,ISXY,ISMAX,ISMIN,ISVM,ISYZ,ISXZ,ISVMAX,
	IFX,IFY,IFXY,IFMAX,IFMIN,IFVM,
    IMX,IMY,IMXY,IMMAX,IMMIN,
	IVYZ,IVXZ,IVMAX,
};
enum {
	CONCENTRATED, UNIFORM, TRAPEZOIDAL, 
	FORCE, DISPLACEMENT, MASS, 
	TEMPERATURE, STRAIN
};
enum {
	F_ADD, F_REPLACE, F_DELETE
};
enum {
	RECTANGULAR, CIRCULAR,WIDEFLANGE,CHANNEL,DOUBLE_CHANNEL,TEE,ANGLE,DOUBLE_ANGLE,BOX,PIPE,GENERAL
};
enum {
	GLOBAL, LOCAL
};
enum {
	COULMN, BEAM, ASLAB
};
enum {
	PERIODIC, TRANSIENT
};
enum {
	MODAL, DIRECTINT
};
enum {
	PERMANENT, TEMPORARY
};
enum {
	LEVEL1, LEVEL2, LEVEL3,LEVEL4
};
enum {
	LINEARCOMBO, ENVELOPECOMBO, ABSSUMCOMBO, SRSSCOMBO
};
enum {
	ABSSUM, SRSS, CQC
};
enum {
	CARTESIAN, RADIAL
};

/*
Serialize controller
*/
class CmyDatabase;

class CSerializer {
public:
	int type;
	BOOL is_storing;
	CArchive* par;
	CDatabase* pdb;
	CmyDatabase* pmydb;

	CString str_values;
	CString str_name;
	CString str_name_type;
	int table_index;
	int skip;
	int moving_index;

	enum {
		ARCHIVE, DB, MYDB
	};
public:
	CSerializer() {
		skip = false;
	}
	void Set(CArchive& ar,BOOL store) {
		par = &ar;
		type = ARCHIVE;
		is_storing = store;
	}
	void Set(CDatabase& db,BOOL store) {
		pdb = &db;
		type = DB;
		is_storing = store;
	}
	void Set(CmyDatabase& mydb,BOOL store) {
		pmydb = &mydb;
		type = MYDB;
		is_storing = store;
	}
	BOOL is_db() {
		return (type != ARCHIVE);
	}
	void CreateTable(CString);
	void SetStart();
	void SetEnd();
};

/*
Data
*/
class DATA {
public:
	CString name;
	int rank;
public:
    friend int operator == (const DATA& left,const DATA& right) {
		return(left.name == right.name);
	}
	int Modify(void*,void*) {
		return IDCANCEL;
	};
	void Serialize(CSerializer&) {
	}
};
/*
Database
*/
class CmyData {
public:
	CString name;
	int nFields;
	CString field_names;
	CList<CString,CString&> data;
public:
	void GetFieldsInString(CString,CString&,int&,int,BOOL);
    void operator = (const CmyData& right) {
		name = right.name;
		nFields = right.nFields;
		field_names = right.field_names;
		data.RemoveAll();
		CString v;
		POSITION pos = right.data.GetHeadPosition();
		while(pos) {
			v = right.data.GetNext(pos);
			data.AddTail(v);
		}
	}
	friend int operator == (const CmyData& left,const CmyData& right) {
		return(left.name == right.name);
	}
	static int tokenize(CString str, CString* tokens, const CString str2);
};

class CmyDatabase {
public:
	CList<CmyData,CmyData&> database;
public:
	void Create(CString str);
	void Insert(CString str);
	void Clear();
};
/*
Template sorted list
*/
template <class T,class Tr>
class CSortedList : public CList<T,Tr> {
public:
	CSortedList() {
	}
	CSortedList(const CSortedList& right) {
		*this = right;
	}
	void operator = (const CSortedList& right) {
		T v,v1;
		POSITION pos,pos1,pos2;
		BOOL skip;
		RemoveAll();
		pos = right.GetHeadPosition();
		while(pos) {
			v = right.GetNext(pos);
			
			/*sorted copy*/
			skip = FALSE;
			pos1 = GetHeadPosition();
			while(pos1) {
				pos2 = pos1;
				v1 = GetNext(pos1);
				if(v1 == v) {
					skip = TRUE;
					continue;
				}
				if(v1 > v)  {
					pos1 = pos2;
					break;
				}
			}
			if(skip) continue;
			if(pos1) InsertBefore(pos1,v);
			else AddTail(v);
		}
	}
	void AddSorted(T& right) {
		T v1;
		POSITION pos = GetHeadPosition(),ipos = NULL,ppos;
		while(pos) {
			ppos = pos;
			v1 = GetNext(pos);
			if(right == v1)
				return;
			if(v1 > right) {
				ipos = ppos;
				break;
			}
		}
		if(ipos) InsertBefore(ipos,right);
		else AddTail(right);
	}
};
/*
Grids and coordinate system
*/
class SYSTEM : public DATA {
public:
    GRIDLIST grid[3];
	RPoint   origin;
	RPoint   rotation;
	int coordinate;
public:
	SYSTEM() {
		coordinate = CARTESIAN;
		name = "CS";
		rank = LEVEL1;
	}
	SYSTEM(const SYSTEM& right) {
		*this = right;
	}
	void operator = (const SYSTEM& right) {
        name = right.name;
		rank = right.rank;
		coordinate = right.coordinate;
		origin = right.origin;
		rotation = right.rotation;
		for(int i = IUX;i <= IUZ;i++)
			grid[i] = right.grid[i];
	}
	
    void Serialize(CSerializer&);
	int Modify(void*,void*);
};

class GRIDINFO {
public:
	DOUBLE xmin;
	DOUBLE xmax;
	DOUBLE ymin;
	DOUBLE ymax;
	DOUBLE zmin;
	DOUBLE zmax;
public:
	GRIDINFO() {
		Reset();
	}
	void Reset() {
		xmin = 1e6;
		xmax = -1e6;
		ymin = 1e6;
		ymax = -1e6;
		zmin = 1e6;
		zmax = -1e6;
	}
};
/*
load cases
*/
class NLCASE;

#define CASETYPES 7

class ANALYSISCASE  : public DATA {
public:
	int index;
	BOOL finished;
	BOOL run;
	int acase_type;
	CString start_case;
	BOOL selected; 
	enum {
		LOAD_CASE,MODAL_CASE,RESPONSEH_CASE,RESPONSES_CASE,COMBO_CASE,NL_CASE,BUCKLING_CASE
	};
public:
	ANALYSISCASE() {
		name = "";
		rank = LEVEL1;
		index = 0;
		finished = FALSE;
		run = TRUE;
		acase_type = LOAD_CASE;
		start_case = "";
	}
	ANALYSISCASE(ANALYSISCASE& right) {
		*this = right;
	}
	void operator = (const ANALYSISCASE& right) {
		name = right.name;
		rank = right.rank;
	    index = right.index;
	    finished = right.finished;
	    run = right.run;
		acase_type = right.acase_type;
		start_case = right.start_case;
	}
	void Serialize(CSerializer&);
};
/*
Static Load case
*/
class LOADCASE : public ANALYSISCASE {
public:
	DOUBLE swm;
	int type;
	enum {DEAD, LIVE, WIND, QUAKE};
public:
	LOADCASE() {
		name = "LDC";
		swm = 0;
		type = DEAD;
		acase_type = LOAD_CASE;
	}
	void Serialize(CSerializer& ar);
	int Modify(void*,void*);
};
/*
Modal anaysis case
*/
class MODALCASE : public ANALYSISCASE {
public:
	UINT minm;
	UINT maxm;
	UINT selmode;
	UINT runmodes;
	DOUBLE tolerance;
	DOUBLE shift;
	VECTOR eigvalue;
	VECTOR eigvector;
public:
	MODALCASE() {
		name = "MODAL";
		minm = 1;
		maxm = 12;
		selmode = 0;
		runmodes = 0;
		tolerance = 1e-7;
		shift = 0;
		eigvalue = 0;
		eigvector = 0;
		acase_type = MODAL_CASE;
	}
	void Serialize(CSerializer& ar);
	int Modify(void*,void*);
};
/*
Time History function
*/
enum {
	F_EBCS_DESIGN,F_EBCS_ELASTIC,F_USER,F_SINE,F_COSINE,F_TRIANGULAR,F_PERIODIC
};
enum {
	SOILA, SOILB, SOILC
};
class FUNCTION  : public DATA {
public:
	DOUBLE period;
	DOUBLE amplitude;
	int divisions;
	int cycles;
	int type;
	DOUBLE Ag;
	int soil_type;
    SRPOINTLIST points;
	static BOOL defspectrumtype;
public:
	FUNCTION() {
		name = "FUNC";
		rank = LEVEL1;
		type = F_USER;
		period = 1;
		amplitude = 1;
		cycles = 1;
		divisions = 20;
		Ag = 0.1;
		soil_type = SOILA;
	}
	FUNCTION(FUNCTION& right) {
		*this = right;
	}
	void operator = (const FUNCTION& right) {
        name = right.name;
		rank = right.rank;
		cycles = right.cycles;
		type = right.type;
		period = right.period;
		amplitude = right.amplitude;
		cycles = right.cycles;
		divisions = right.divisions;
		Ag = right.Ag;
		soil_type = right.soil_type;
		points = right.points;
	}
	DOUBLE Get(DOUBLE x,int type,DOUBLE scale,BOOL = FALSE);

	void Serialize(CSerializer& ar);
	int Modify(void*,void*);
};
/*
spectrum function
*/
class SPECFUNC {
public:
	FUNCTION* function;
	int dir;
	DOUBLE scale;
public:
	SPECFUNC() {
		function = NULL;
		dir = IUX;
		scale = 1.0;
	}
	friend int operator == (const SPECFUNC& left,const SPECFUNC& right) {
		return(left.function == right.function && left.dir == right.dir && left.scale == right.scale);
	}
	void Serialize(CSerializer& ar,PFUNCTIONLIST flist);
};
class DAMPING {
public:
	DOUBLE w1,w2,e1,e2,fm,fk;
	DOUBLE cdamping;
	int type;
	FUNCTION dvalues;
	enum {RAYLEIGH, CONSTANT,NONE};
public:
	DAMPING() {
		w1 = 1;
		w2 = 10;
		e1 = 0.05;
		e2 = 0.05;
		fm = 0.0909090909;
		fk = 0.00909090909;
		cdamping = 0.05;
		type = NONE;
	}
	DAMPING(DAMPING& right) {
		*this = right;
	}
	void operator = (const DAMPING& right) {
		w1 = right.w1;
		w2 = right.w2;
		e1 = right.e1;
		e2 = right.e2;
		fm = right.fm;
		fk = right.fk;
		cdamping = right.cdamping;
		type = right.type;
		dvalues = right.dvalues;
	}
	void Calculate() {
		RPoint v;
		DOUBLE value;
		UINT i,maxm;

		dvalues.points.RemoveAll();
		
		if(type == RAYLEIGH) {
			maxm = UINT(w2);
			value = maxm / w2;

			fm = 2.0 * ((w2 * w2 * w1 * e1) - (w1 * w1 * w2 * e2)) / (w2 * w2 - w1 * w1);
			fk = 2.0 * ((w2 * e2) - (w1 * e1)) / (w2 * w2 - w1 * w1);
			
			for(i = 1;i < 5;i++) {
				v.x = i * 0.2;
				v.y = fm / (2 * v.x) + fk * v.x / 2;
				dvalues.points.AddTail(v);
			}
			for(i = 1;i < maxm;i++) {
				v.x = i * value;
				v.y = fm / (2 * v.x) + fk * v.x / 2;
				dvalues.points.AddTail(v);
			}
		} else if(type == CONSTANT) {
			for(i = 0;i <= 10;i+=2) {
				v.x = i;
				v.y = cdamping;
				dvalues.points.AddTail(v);
			}
		}
	}
	void Serialize(CSerializer& ar);
};
class RESPONSE : public ANALYSISCASE {
public:
	DAMPING damping;
	MODALCASE* modalcase;
	SPECFUNCLIST funclist;
	static MODALCASE* defmodalcase;
	static FUNCTION* deffunction;
public:
	void ConstructDampingMatrix(MATRAN C,UINT SJ);
	RESPONSE() {
		modalcase = defmodalcase;
	}
	~RESPONSE() {
	}
	void operator = (const RESPONSE& right) {
		ANALYSISCASE::operator = (right);
		damping = right.damping;
		modalcase = right.modalcase;
		POSITION pos;
		SPECFUNC v;
		funclist.RemoveAll();
		pos = right.funclist.GetHeadPosition();
		while(pos) {
             v = right.funclist.GetNext(pos);
             funclist.AddTail(v);
		}
	}
	void Serialize(CSerializer& ar);
	void Serialize(CSerializer& ar,PFUNCTIONLIST);
};
/*
Response analysis case
*/
class RESPONSEHIST : public RESPONSE {
public:
	DOUBLE dt;
	UINT N;
	UINT seltime;
	int type;
	int ana_type;
public:
	RESPONSEHIST(RESPONSEHIST& right) {
		*this = right;
	}
	RESPONSEHIST() {
		name = "HIST";
		dt = 0.1;
		N = 10;
		type = TRANSIENT;
		ana_type = MODAL;
		seltime = 0;
		acase_type = RESPONSEH_CASE;
	}
	void operator = (const RESPONSEHIST& right) {
		RESPONSE::operator = (right);
		dt = right.dt;
		N = right.N;
		seltime = right.seltime;
		type = right.type;
		ana_type = right.ana_type;
	}
	void Serialize(CSerializer& ar);
	int Modify(void*,void*);
};
/*
Response Spectrum Case
*/
class RESPONSESPEC: public RESPONSE {
public:
	int modal_comb;
	int dir_comb;
public:
	RESPONSESPEC(RESPONSESPEC& right) {
		*this = right;
	}
	RESPONSESPEC() {
		name = "SPEC";
		modal_comb = CQC;
		dir_comb = ABSSUM;
		acase_type = RESPONSES_CASE;
	}
	void operator = (const RESPONSESPEC& right) {
		RESPONSE::operator = (right);
		modal_comb = right.modal_comb;
		dir_comb = right.dir_comb;
	}
	void Serialize(CSerializer& ar);
	int Modify(void*,void*);
};
/*
Load Case Combo
*/
class LOADCOMBO {
public:
    ANALYSISCASE* loadcase;
	DOUBLE FS;
public:
	friend int operator == (const LOADCOMBO& left,const LOADCOMBO& right) {
		return(left.loadcase->name == right.loadcase->name);
	}
	void Serialize(CSerializer&);
};
/*
Comb_Type
*/
class COMB_TYPE : public ANALYSISCASE {
public:
	LOADCOMBOLIST loadlist;
	int type;
public:
	friend int operator == (const COMB_TYPE& left,const COMB_TYPE& right) {
		return(left.name == right.name);
	}
	COMB_TYPE() {
	}
	COMB_TYPE(COMB_TYPE& right) {
		*this = right;
	}
	void operator = (const COMB_TYPE& right) {
		ANALYSISCASE::operator = (right);
		type = right.type;
		LOADCOMBO v;
		loadlist.RemoveAll();
		POSITION pos = right.loadlist.GetHeadPosition();
		while(pos) {
             v = right.loadlist.GetNext(pos);
             loadlist.AddTail(v);
		}
	}
	void Serialize(CSerializer&);
};
/*
Load combinations
*/
class COMBINATION : public COMB_TYPE {
public:
	BOOL designcombo;
	int selcombo;
	enum {CMAX, CMIN};
public:
	COMBINATION() {
		name = "COMB";
		type = LINEARCOMBO;
		acase_type = COMBO_CASE;
		selcombo = CMAX;
		designcombo = FALSE;
	}
	void operator = (const COMBINATION& right) {
		COMB_TYPE::operator = (right);
		selcombo = right.selcombo;
		designcombo = right.designcombo;
	}
	int Modify(void*,void*);
	void Serialize(CSerializer&);
};
/*
Non-linear analysis case
*/
class NLCASE : public COMB_TYPE {
public:
	MATRAN end_stiffness;
	UINT nSteps;
	UINT nIteration;
	DOUBLE tolerance;
	BOOL save_stiffness;
	enum {PDELTA_LINEARIZED,PDELTA_CONSISTENT};
public:
	NLCASE() {
		nSteps = 200;
		nIteration = 10;
		tolerance = 1e-4;
		name = "NLCASE";
		type = PDELTA_LINEARIZED;
		acase_type = NL_CASE;
		end_stiffness = NULL;
		save_stiffness = FALSE;
	}
	int Modify(void*,void*);
    void Serialize(CSerializer&);
};
class BUCKLINGCASE : public COMB_TYPE {
public:
	UINT nmodes;
	UINT runmodes;
	int selmode;
	VECTOR eigvalue;
	VECTOR eigvector;
	DOUBLE tolerance;
	enum {LINEAR};
public:
	BUCKLINGCASE() {
		nmodes = 6;
		runmodes = 0;
		selmode = 0;
		tolerance = 1e-7;
		type = LINEAR;
		name = "BUCKL";
		eigvector = 0;
		eigvalue = 0;
		acase_type = BUCKLING_CASE;
	}
	int Modify(void*,void*);
    void Serialize(CSerializer&);
};
/*
Load
*/
class BASE_LOAD {
public:
    UBMP8 type;
	SYSTEM* system;
	int casetype;
	LOADCASE* loadcase;
	static DOUBLE factor;
	static DOUBLE load_step_factor;
	static PLOADCOMBOLIST cloadcombolist;
public:
	BASE_LOAD() {
		system = NULL;
		loadcase = NULL;
		casetype = PERMANENT;
	}
	virtual void ApplyFactor() = 0;
	virtual void RemoveFactor() = 0;
	BOOL is_joint_load() {
		return (type == FORCE || type == DISPLACEMENT);
	}
	BOOL is_member_load() {
		return (type == CONCENTRATED || type == UNIFORM || type == TRAPEZOIDAL);
	}
	BOOL is_slab_load() {
		return (type == UNIFORM);
	}
	BOOL is_strain_load() {
		return (type == TEMPERATURE || type == STRAIN);
	}
	BOOL is_considered();
	void Serialize(CSerializer& ar);
};
/*
Joint Load
*/
class JLOAD : public BASE_LOAD {

public:
	DOUBLE Q[6];
public:
	virtual void ApplyFactor();
	virtual void RemoveFactor();
	JLOAD() {
		for(int i = IUX;i <= IRZ;i++) 
			Q[i] = 0;
		type = FORCE;
	}
	void Serialize(CSerializer& ar);
};
/*
Member Load
*/
class LOAD : public BASE_LOAD {
public:
	UBMP8 dir;
	DOUBLE P;
	DOUBLE x;
	DOUBLE P1;
	DOUBLE x1;
public:
	virtual void ApplyFactor();
	virtual void RemoveFactor();
	LOAD() {
		type = CONCENTRATED;
		dir = IUY;
		P = 0;
		x = 0;
		P1 = 0;
		x1 = 0;
	}
	void Serialize(CSerializer& ar);
};
/*
ENTITY
*/
#define  NENTITIES   3

class ENTITY {
public:
	int sel;
	CString name;
	void* ptr;

	static  UBMP32 TotalLoadCases;
public:
	ENTITY() {
		sel = false;
		name = "";
		ptr = 0;
	}
	void operator = (const ENTITY& right) {
		sel = right.sel;
		name = right.name;
		ptr = right.ptr;
	}
	virtual void AllocVectors() {};
	virtual void FreeVectors() {};
    virtual void AllocAnalysis(int) {};
	virtual void FreeAnalysis(int) {};
	virtual void SerializeAnalysis(CSerializer&,int,ANALYSISCASE*) {};
};

/*
JOINT
*/
class JOINT : public ENTITY {
public:
	RPoint  p;
	UBMP8   restraint;
	UBMP8   constraint;
	int     number[6];
	JLOADLIST load;
	VECTOR   forces;
	VECTOR   disps;
	MATRIX   forces_all;
	MATRIX   disps_all;
	RPoint   rotation;
	UINT     mconnect;
	UINT     cconnect;
	UINT     sconnect;
	DOUBLE   mself;
	DOUBLE   massembled[6];
	int      NodeNumber;
	static  UBMP32 TotalNumber;
public:
    JOINT() {
		restraint = NOUR;
		constraint = NOUR;
		for(int i = IUX;i <= IRZ;i++) {
			number[i] = INVALID;
			massembled[i] = 0;
		}
		mconnect = 0;
		cconnect = 0;
		sconnect = 0;
		forces = NULL;
		disps = NULL;
		forces_all = NULL;
		disps_all = NULL;
	}
	~JOINT() {
	}
	void AllocVectors();
	void FreeVectors();
	void AllocAnalysis(int);
	void FreeAnalysis(int);
	void ApplyLoad(VECTOR,VECTOR);
    void GetRotMatrix(MATRIX);

	friend int operator == (const JOINT& left,const JOINT& right) {
		return(left.p == right.p);
	}
	JOINT(JOINT& right) {
		*this = right;
	}
	void operator = (const JOINT& right) {
		ENTITY::operator = (right);
		for(int i = IUX;i <= IRZ;i++) {
			number[i] = right.number[i];
			massembled[i] = right.massembled[i];
		}
		rotation = right.rotation;
		mself = right.mself;
        p = right.p;
		restraint = right.restraint;
		constraint = right.constraint;
		mconnect = right.mconnect;
		cconnect = right.cconnect;
		sconnect = right.sconnect;
	
		JLOAD v;
		load.RemoveAll();
		POSITION pos = right.load.GetHeadPosition();
		while(pos) {
             v = right.load.GetNext(pos);
             load.AddTail(v);
		}
	}

	void Serialize(CSerializer& ar);
	void SerializeAnalysis(CSerializer&,int,ANALYSISCASE*);
};
/*
Rigid Body
*/
class CONSTRAINT : public DATA {
public:
	int type;
	BOOL isauto;
	UBMP8 constraint;
	JOINTPLIST jointlist;
	SYSTEM* axis;
	DOUBLE tolerance;
	static DOUBLE PenaltyWeight;
	enum {BODY,DIAPHRAGM,EQUAL,PLATE,BEAM,ROD,LINE,LOCAL,WELD};
public:
	CONSTRAINT() {
		name = "BODY";
		rank = LEVEL1;
		type = EQUAL;
		constraint = ALLUR;
		axis = 0;
		isauto = FALSE;
		tolerance = 0.01;
	}
	~CONSTRAINT() {
	}
	CONSTRAINT(CONSTRAINT& right) {
		*this = right;
	}
	void operator = (const CONSTRAINT& right) {
        name = right.name;
		rank = right.rank;
		type = right.type;
		isauto = right.isauto;
		constraint = right.constraint;
		axis = right.axis;
		tolerance = right.tolerance;
		JOINT* v;
		jointlist.RemoveAll();
		POSITION pos = right.jointlist.GetHeadPosition();
		while(pos) {
             v = right.jointlist.GetNext(pos);
             jointlist.AddTail(v);
		}
	}
	bool DetermineAxis(RPoint&,RPoint);
	int Modify(void*,void*);
	void Serialize(CSerializer& ar);
	void AddToK(MATRAN,UINT);
};
/*
Reinforcement property
*/
class REBAR {
public:
	int design;
	int type;
	UBMP8 nz;
	UBMP8 ny;
	UBMP8 nt;
	DOUBLE barsize;
	DOUBLE stirrup_barsize;
	DOUBLE cover;
public:
	REBAR() {
		design = COULMN;
		type = RECTANGULAR;
		nz = 2;
		ny = 2;
		nt = 8;
		cover = 0.025;
		barsize = 0.012;
		stirrup_barsize = 0.006;
	}
	~REBAR() {
	}
	friend int operator == (const REBAR& left,const REBAR& right) {
		return( (left.design == right.design)
			 && (left.type == right.type)
			 && (left.nz == right.nz)
			 && (left.ny == right.ny)
			 && (left.nt == right.nt)
			 && (EQUAL(left.barsize,right.barsize))
			 && (EQUAL(left.stirrup_barsize,right.stirrup_barsize))
			 && (EQUAL(left.cover,right.cover))
			);
	}
	void Serialize(CSerializer& ar);
};
/*
Material property
*/
class MATERIAL  : public DATA {
public:
	DOUBLE E;
	DOUBLE nu;
	DOUBLE unitweight;
	DOUBLE density;
	DOUBLE alphac;
	DOUBLE G;
    UBMP8 type;
	DOUBLE fck;
	DOUBLE fctk;
	DOUBLE fyk;
	DOUBLE fyks;
	DOUBLE ftk;
	DOUBLE ftks;
	DOUBLE Es;
	enum {CONCRETE, STEEL, FILLER};
public:
	MATERIAL() {
		name = "MAT";
		rank = LEVEL1;
		E = 32.0e+6;
		nu = 0.25;
		unitweight = 25;
		density = 2.4;
		alphac = 10.8e-6;
		G = E / (2 * (1 + nu)) ;
		type = CONCRETE;

		fck = 27579.032;
		fctk = 17000;
		fyk = 413685.8;
		fyks = 275790.32;
		ftk = 413685.8;
		ftks = 275790.32;
		Es = 200.0e+6;
	}
	~MATERIAL() {
	}
	int Modify(void*,void*);
	void Serialize(CSerializer& ar);
};
/*
Frame Section property
*/
class MEMBER;
class DETAILING;

class SECTION : public DATA {
public:
	UBMP8 type;
	REBAR rebar;
	MATERIAL* material;
	/*values*/
	DOUBLE w;
	DOUBLE h;
	DOUBLE r;
	DOUBLE tw;
	DOUBLE tf;
	DOUBLE bd;
	DOUBLE A;
	DOUBLE Ay;
	DOUBLE Az;
	DOUBLE Ix;
	DOUBLE Iy;
	DOUBLE Iz;
	DOUBLE ry;
	DOUBLE rz;
	DOUBLE Zy;
	DOUBLE Zz;
	DOUBLE Sy;
	DOUBLE Sz;
	/*stiffness modification factors*/
	DOUBLE fA;
	DOUBLE fAy;
	DOUBLE fAz;
	DOUBLE fIx;
	DOUBLE fIy;
	DOUBLE fIz;
	static MATERIAL* defmaterial;
	static DETAILING* pDocDetailing; 
	/*steel design*/
	int sclass;
	BOOL built_up;
	DOUBLE dA;
	DOUBLE dZ[2];
	DOUBLE eA[2];
	DOUBLE fm[2];
	DOUBLE si[2];
	enum {CLASS1,CLASS2,CLASS3,CLASS4};
	enum {BA,BB,BC,BD};
public:
	SECTION() {
		name = "SEC";
		rank = LEVEL1;
		type = RECTANGULAR;
		w = 0.4;
		h = 0.4;
		r = 0.2;
		tw = 0.002;
		tf = 0.002;
		bd = 0;
		A  = 0.16;
		Ay = 0.13333333333333;
		Az = 0.13333333333333;
		Ix = 0.00360533348719;
		Iy = 0.00213333333333;
		Iz = 0.00213333333333;
		ry = 0.1154701;
		rz = 0.1154701;
		Zy = 0.01066666666666;
        Zz = 0.01066666666666;
		Sy = 0.016;
		Sz = 0.016;
		fA = 1.0;
		fAy = 1.0;
		fAz = 1.0;
		fIx = 1.0;
		fIy = 1.0;
		fIz = 1.0;
		material = defmaterial;
		built_up = FALSE;
	}
	~SECTION() {
	}
	void Classify(MEMBER*);
	int SelectBucklingCurve(int);
	void GetData();
	
	int Modify(void*,void*);
	void Serialize(CSerializer& ar);
};
/*
General Abstract Member Class from which all
other types of members are derived
*/
struct DESIGNDATA {
	UINT flags;
	DOUBLE As;
	DOUBLE sd;
	UINT s;
	DOUBLE Area[2];
	DOUBLE diam[2][2];
	UINT count[2][2];
	UINT total[2];
	DOUBLE ratios[8];
	enum {TOP,BOTTOM};
	enum {INSUFF_TOP = 1,INSUFF_BOTTOM = 2,INSUFF_SHEAR = 4,TOP_DOUBLE_ROW = 8,BOT_DOUBLE_ROW = 16,SKIPPED = 32};
	DESIGNDATA() {
		flags = SKIPPED;
		s = 0;
		sd = 0;
		As = 0;
		total[TOP] = total[BOTTOM] = 0;
		Area[TOP] = Area[BOTTOM] = 0;
		total[TOP] = total[BOTTOM] = 0;
		count[TOP][0] = count[TOP][1] = 0;
		count[BOTTOM][0] = count[BOTTOM][1] = 0;
		diam[TOP][0] = diam[TOP][1] = 0;
		diam[BOTTOM][0] = diam[BOTTOM][1] = 0;
		for(int i = 0;i < 8;i++) ratios[i] = 0;
	}
	void Serialize(CSerializer&);
};

#define BARSIZELIST CSortedList<DOUBLE,DOUBLE&> 

class DETAILING {
public:
	BOOL single_shear;
	BOOL cut_at_L3;
	int  le_option;
	int  frame_type;
	BOOL accidental_eccentricity;
	BOOL individual_buckling;
	UINT nWDIV;
	UINT nQDIV;
	UINT min_spacing;
	int max_layers;
	DOUBLE fs_concrete;
	DOUBLE fs_steel;
	BARSIZELIST shearbarlist;
	BARSIZELIST beambarlist;
	BARSIZELIST columnbarlist;
	BARSIZELIST slabbarlist;
	enum {LORIGINAL,LCALCULATED,LASSIGNED,LPREF};
	enum {SWAY,NON_SWAY,SPREF};
public:
	DETAILING() {
		single_shear = TRUE;
		cut_at_L3 = TRUE;
		min_spacing = 40;
		max_layers = 1;
		fs_concrete = 1.5;
		fs_steel = 1.15;
		nWDIV = 12;
		nQDIV = 16;
		le_option = LCALCULATED;
		frame_type = SWAY;
		individual_buckling = TRUE;
		accidental_eccentricity = TRUE;
	}
    void Serialize(CSerializer& ar);
};

class ELEMENT : public ENTITY {
public:
	LOADLIST load;
public:
	ELEMENT() {
	}
	void operator = (const ELEMENT& right) {
		ENTITY::operator = (right);
		LOAD v;
		load.RemoveAll();
		POSITION pos = right.load.GetHeadPosition();
		while(pos) {
			v = right.load.GetNext(pos);
			load.AddTail(v);
		}
	}
	static void AssembleJoints(JOINT** jt,int NJ,MATRAN K,MATRIX Kl,int NB);
};

class MEMBER : public ELEMENT {
public:
	JOINT* j1;
    JOINT* j2;
	DOUBLE alpha;
	SECTION* section;
	SECTION* e_section;
	int EIyy;
	int EIzz;
	UBMP8 nrelease;
	UBMP8 frelease;
	MATRIX forces;
	MATRIX disps;
	MATRIX* forces_all;
	MATRIX* disps_all;
	VECTOR station;
	DESIGNDATA* design_data;
	UBMP32 nDiv;
	UBMP32 nMinDiv;
	UBMP32 nMinFrameDiv;
	BOOL   DivAtInterim;
	DOUBLE start_offset;
	DOUBLE end_offset;
	MATRIX K_fem;
	BOOL is_curved;
	RPOINTLIST curved_list;

	/*effective length*/
	int le_option;
	DOUBLE Le[2];
	DOUBLE LeCalculated[2];
	/*moment redistribution factor*/
	DOUBLE redistribution_factor;
	/*sway/non-sway*/
	int frame_type;

	static  BOOL add_geometric_stiffness;
	static  BOOL only_geometric_stiffness;
	static  int p_delta;
	static  UBMP32 TotalNumber;
	static  SECTION* defsection;
	enum {LINEAR,PARABOLIC,CUBIC};
	enum {LINEARIZED, CONSISTENT};
public:
	MEMBER() {
		j1 = NULL;
		j2 = NULL;
		nrelease = NOUR;
		frelease = NOUR;
		forces = NULL;
		forces_all = NULL;
		disps = NULL;
		disps_all = NULL;
		station = NULL;
		alpha = 0;
		section = defsection;
		e_section = NULL;
		nMinDiv = 9;
		nMinFrameDiv = 1;
		DivAtInterim = TRUE;
		nDiv = 1;
		EIyy = PARABOLIC;
		EIzz = LINEAR;
		start_offset = 0;
		end_offset = 0;
		is_curved = FALSE;

		frame_type = DETAILING::SPREF;
		le_option = DETAILING::LPREF;
		for(UINT i = 0;i < 2;i++) {
			Le[i] = 0;
			LeCalculated[i] = 0;
		}
		redistribution_factor = 0;
	}

	~MEMBER() {
	}

	DOUBLE GetLength();
	void AllocVectors();
	void FreeVectors();
	void AllocAnalysis(int);
	void FreeAnalysis(int);
	void AssignSections();
	void GetTransformationMatrix(MATRIX);
	void GetLoadTransMatrix(MATRIX,SYSTEM*);
	void GetLocalStiffnessMatrix(MATRIX,int = 0);
	void GetRotMatrix(MATRIX,int,int = -1,int = -1);
	void ApplyLoad(VECTOR);
	void CalcMF(VECTOR Qg,VECTOR Dg,VECTOR q,VECTOR d);
	void CalculateLocalConcentratedFEM(VECTOR,DOUBLE,DOUBLE,DOUBLE,int);
	void CalculateLocalFEM(LOAD*,VECTOR);
	void CalculateFEM(VECTOR);
    void CalculateLocalConcentratedIF(VECTOR,DOUBLE,DOUBLE,DOUBLE,int);
	void CalculateLocalIF(LOAD* ,DOUBLE, VECTOR);
	void CalculateInternalForceAndDisplacement(VECTOR,VECTOR,BOOL = FALSE);
    DOUBLE _CalculateInternal(BOOL bforce,DOUBLE x,int dir,int index);
	DOUBLE CalculateForce(DOUBLE x, int dir,int index = -1);
	DOUBLE CalculateDeflection(DOUBLE x,int dir,int index = -1);
	void CalculateDeflection(DOUBLE x,RPoint& v,int index = -1);
	BOOL Mesh(PMEMBERLIST,PJOINTLIST,MEMBERPLIST_PLIST*);
    BOOL MeshCurved(PMEMBERLIST,PJOINTLIST,MEMBERPLIST_PLIST*);
	void Assemble(MATRAN,MATRIX,int);
	void AddToK(MATRAN,int);
	void DetermineFlexibility(MATRIX,DOUBLE,BOOL);
	void Design(DETAILING*);
	void DesignBeam(DETAILING*);
	void DesignColumn(DETAILING*);
    void DesignSteelColumn(DETAILING*);
	void ApplyGeometricStiffness(MATRIX K);
	void ApplyReleases(MATRIX K);
	void ApplyReleases(VECTOR F);
	void Break(JOINTPLIST*,PMEMBERLIST,MEMBERPLIST_PLIST* = NULL);
	void SwapEnds();
	void GetRestraint(VECTOR);
	int FindCurvedPoint(DOUBLE,DOUBLE* = 0);
	RPoint FindPoint(DOUBLE);
	static void Join(MEMBERPLIST*,PMEMBERLIST);

	friend int operator == (const MEMBER& left,const MEMBER& right) {
		return(((left.j1->p == right.j1->p) && (left.j2->p == right.j2->p)) ||
			   ((left.j1->p == right.j2->p) && (left.j2->p == right.j1->p)));
	}
	MEMBER(MEMBER& right) {
		*this = right;
	}
	void operator = (const MEMBER& right) {
		ELEMENT::operator = (right);
		CopyNormalData(&right);
		CopyCurveData(&right);
	}
	void CopyNormalData(const MEMBER* right) {
		j1 = right->j1;
		j2 = right->j2;
		alpha = right->alpha;
		section = right->section;
		e_section = right->e_section;
		EIyy = right->EIyy;
		EIzz = right->EIzz;
		nrelease = right->nrelease;
		frelease = right->frelease;
		nMinDiv = right->nMinDiv;
        nMinFrameDiv = right->nMinFrameDiv;
		DivAtInterim = right->DivAtInterim;
		nDiv = right->nDiv;
		start_offset = right->start_offset;
		end_offset = right->end_offset;

		le_option = right->le_option;
		for(UINT i = 0;i < 2;i++) {
			Le[i] = right->Le[i];
			LeCalculated[i] = right->LeCalculated[i];
		}
		redistribution_factor = right->redistribution_factor;
	}
	void CopyCurveData(const MEMBER* right) {
		is_curved = right->is_curved;
        if(is_curved) {
			RPoint rp;
			curved_list.RemoveAll();
			POSITION pos = right->curved_list.GetHeadPosition();
			while(pos) {
				rp = right->curved_list.GetNext(pos);
				curved_list.AddTail(rp);
			}
		}
	}
	void Serialize(CSerializer& ar);
	void SerializeAnalysis(CSerializer& ar,int,ANALYSISCASE*);
};
/*
Area Section property
*/
class ASECTION : public DATA{
public:
	DOUBLE h;
	int type;
	int state_type;
	REBAR rebar;
	MATERIAL* material;
	static MATERIAL* defmaterial;
	enum {Q4 = 4,Q5,Q6,Q7,Q8,Q9,Q4D,PLATE4,SHELL4};
	enum {PLANE_STRESS,PLANE_STRAIN,PLATE_ACTION,SHELL_ACTION};
public:
	ASECTION() {
		name = "ASEC";
		rank = LEVEL1;
		type = SHELL4;
		state_type = SHELL_ACTION;
		material = defmaterial;
		h = 0.25;
	}
	~ASECTION() {
	}
	
	int Modify(void*,void*);
	void Serialize(CSerializer& ar);
};
/*
Area
*/
const int NPOLY = 20;

class SLAB : public ELEMENT {
public:
    JOINT* jt[NPOLY];
	UINT NJ;
	DOUBLE alpha;
	ASECTION* section;
	VECTOR forces;
	VECTOR disps;
	MATRIX forces_all;
	MATRIX disps_all;
	UINT nDivx;
	UINT nDivy;
	BOOL edge_constraint;

	static  UBMP32 TotalNumber;
	static  ASECTION* defsection;
	static  const UINT NSTRESS;
	static  const UINT Rule[];
	enum {FORCE_F,MASS_F,STRESS_F};
public:
	SLAB() {
		NJ = 4;
		for(UINT i = 0;i < NJ;i++)
			jt[i] = NULL;
		alpha = 0;
		section = defsection;
		forces = NULL;
		forces_all = NULL;
		disps = NULL;
		disps_all = NULL;
		nDivx = 1;
		nDivy = 1;
		edge_constraint = FALSE;
	}
	~SLAB() {
	}
	friend int operator == (const SLAB& left,const SLAB& right) {
		BOOL found;
		UINT NJ = left.NJ;
		if(right.NJ != NJ) 
			return 0;
		for(UINT i = 0;i < NJ;i++) {
			found = FALSE;
			for(UINT j = 0;j < NJ;j++) {
				if(left.jt[i]->p == right.jt[j]->p) {
					found = TRUE;
					break;
				}
			}
			if(!found) return 0;
		}
		return 1;
	}
	SLAB(SLAB& right) {
		*this = right;
	}
	void operator = (const SLAB& right) {
		ELEMENT::operator = (right);
		NJ = right.NJ;
		for(UINT i = 0;i < NJ;i++)
			jt[i] = right.jt[i];
		alpha = right.alpha;
		section = right.section;
		nDivx = right.nDivx;
		nDivy = right.nDivy;
		edge_constraint = right.edge_constraint;
	}
	DOUBLE GetArea();
	void GetRotMatrix(MATRIX T,int);
	void GetLoadTransMatrix(MATRIX,SYSTEM*);
	void GetTransformationMatrix(MATRIX);
	void __GetLocalStiffnessMatrix(MATRIX);
	void GetLocalStiffnessMatrix(MATRIX);
	void Assemble(MATRAN,MATRIX,int);
	void AddToK(MATRAN,int);
	void __CalculateStresses(VECTOR d);
	void CalculateStresses(VECTOR,VECTOR);
	void CalculateLocalFEM(LOAD*,VECTOR,VECTOR);
	void CalculateFEM(VECTOR);
	void CalculateInitialStress(VECTOR,MATRIX);
	void CalculateInitialFEM(VECTOR);
	void CalculateUniformFactors(VECTOR,int);
	void __CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn,RPoint* p);
	void CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn,RPoint* p);
	void CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR dNe,VECTOR dNn);
    void CalculateJacobian(VECTOR J,DOUBLE& Jdet,DOUBLE e,DOUBLE n);
	void ApplyLoad(VECTOR);
	void AllocVectors();
	void FreeVectors();
	void AllocAnalysis(int);
	void FreeAnalysis(int);
	void GetDKTQStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n);
	void GetQXStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n);
	void GetQDStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR b);
	void GetStrainDisplacementMatrix(MATRIX B,DOUBLE& Jdet,DOUBLE e,DOUBLE n,VECTOR = 0);
	void GetE(MATRIX E);
	void GetLocalJointCoordinates(RPoint*,UINT,RPoint* = 0);

	void MeshRect(PSLABLIST slabs,PJOINTLIST joints,SLABPLIST* mlist);
	void MeshPoly(RPoint* rp,UINT np,UINT n,PJOINTLIST joints,PSLABLIST slabs,SLABPLIST* mlist);

	static void GetShapeFunction(VECTOR Nf,DOUBLE e,DOUBLE n,UINT);
	static void GetShapeFunctionDerivative(VECTOR dN,bool,DOUBLE e,DOUBLE n,UINT);
	
	BOOL Mesh(PSLABLIST,PJOINTLIST,SLABPLIST_PLIST*);
	void AddEdgeConstraint(PJOINTLIST joints,SYSTEM* global,MATRAN K,UINT NB);
	void Serialize(CSerializer& ar);
	void SerializeAnalysis(CSerializer& ar,int,ANALYSISCASE*);
};
/*
Group property
*/
class GROUP : public DATA {
public:
	COLORREF color;
	JOINTPLIST jointlist;
	MEMBERPLIST memberlist;
	SLABPLIST slablist;
public:
	GROUP() {
		name = "GRP";
		rank = LEVEL1;
		color = RGB(255,255,0);
	}
	GROUP(GROUP& right) {
		*this = right;
	}
	void operator = (const GROUP& right) {
        name = right.name;
		rank = right.rank;
		color = right.color;
		JOINT* joint;
		MEMBER* mem;
		SLAB* slab;
		POSITION pos;

		jointlist.RemoveAll();
		pos = right.jointlist.GetHeadPosition();
		while(pos) {
             joint = right.jointlist.GetNext(pos);
             jointlist.AddTail(joint);
		}

		memberlist.RemoveAll();
		pos = right.memberlist.GetHeadPosition();
		while(pos) {
             mem = right.memberlist.GetNext(pos);
             memberlist.AddTail(mem);
		}

		slablist.RemoveAll();
		pos = right.slablist.GetHeadPosition();
		while(pos) {
             slab = right.slablist.GetNext(pos);
             slablist.AddTail(slab);
		}
	}
	
	int Modify(void*,void*);
	void Serialize(CSerializer& ar);
};
/*
DESIGN
*/
struct REINF {
	DOUBLE x;
	DOUBLE y;
};

struct DESIGN {
	DETAILING* detail;
	DESIGNDATA* pdesigndata;

	REINF bar[64];
	int bcount;
	UBMP8 design;
	UBMP8 sectype;
	UBMP8 rebtype;
	DOUBLE h;
	DOUBLE b;
	DOUBLE beamw;
	DOUBLE r;
	DOUBLE hp;
	DOUBLE bp;
	DOUBLE Es;
	DOUBLE fcd;
	DOUBLE fctd;
	DOUBLE fyd;

	DOUBLE vt;
	DOUBLE uht;
	DOUBLE ubt;
	DOUBLE vh;
	DOUBLE vb;
	
	void InitColumn(SECTION*);
	void InitBeam(SECTION*);
	void DoBeam(DOUBLE);
	void DoColumn(DOUBLE*);
    void DoShear(DOUBLE,BOOL = FALSE);

	/*coulmn*/
	DOUBLE FindW(DOUBLE q);
	void FindNearest(DOUBLE q,DOUBLE w,DOUBLE& dmin,DOUBLE& error);
	void CalculateBiaxial(DOUBLE x,DOUBLE q,DOUBLE w,DOUBLE* Pu,DOUBLE* Muz,DOUBLE* Muy);
	void FindMatch(DOUBLE x,DOUBLE q,DOUBLE w,DOUBLE& v,DOUBLE& uh,DOUBLE& ub);
	DOUBLE CalcD(DOUBLE x,DOUBLE y,DOUBLE q);

	/*beam*/
	void CalcAB(DOUBLE kx,DOUBLE dmax,DOUBLE& ac,DOUBLE& bc,DOUBLE& esp);
	DOUBLE CalcMu(DOUBLE kx,DOUBLE dmax);
	DOUBLE FindKx(DOUBLE u,DOUBLE dmax);
};

/*
Tables
*/
struct TABLE {
	char name[128];
	BOOL first_record;
	CRecordset* prs;
	BOOL save;
	int view;
	int parent;

	static int current_view;

	enum {
	     GSettings
		,GColors
		,PMaterial
		,PMSection
		,PASection
		,GSystems
		,GGridLines
		,ALoad_cases
		,AModal_cases
		,RH_function
		,RS_function
		,RH_cases
        ,RH_damping
        ,RH_loading
		,RS_cases
        ,RS_damping
        ,RS_loading
		,ACombinations
		,ANonLinear_cases
		,ABuckling_cases
		,JGeneral
		,JCoordinates
		,JRestraints
		,JLocal_axis
		,JLoads
		,MGeneral
		,MConnectivity
		,MSection
		,MLocal_axis
		,MLoads
		,MReleases
		,MDivisions
		,MCurvedMember
		,SGeneral
		,SConnectivity
		,SSection
		,SLocal_axis
		,SLoads
		,GConstraints
		,GC_Joints
		,GGroups
		,GR_Joints
		,GR_Members
		,GR_Slabs
		,GDetailing
		,GDShear
		,GDBeam
		,GDColumn
		,GDSlab
		,JForces
		,JDisplacements
		,MForces
		,MDisplacements
		,SForces
		,SDisplacements
		,MModal
		,MBuckling

		,Model
		,MSettings
		,MProperties
		,MFunctions
		,MAnalysis_Cases
		,ARH
		,ARS

		,Analysis
		,AJoint
		,AMember
		,ASlab
		,AModes

		,Design
		,MDesign
	};
	enum {
        MODEL,ANALYSIS,DESIGN
	};
};
const int NTABLES = 71;
extern TABLE all_tables[NTABLES];
void clear_tables();

#endif
