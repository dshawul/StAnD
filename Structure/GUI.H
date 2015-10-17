#ifndef __GUI_H
#define __GUI_H

#ifdef _MSC_VER
#    define _CRT_SECURE_NO_DEPRECATE
#    define _SCL_SECURE_NO_DEPRECATE
#endif

#include "stdafx.h"
#include "resource.h"

/*
Enums
*/
enum STATE {
	SELECT, ADD_JOINT, ADD_MEMBER, ADD_CURVED_MEMBER, ADD_RECT_AREA, ADD_QUAD_AREA, QUICK_DRAW_MEMBER, LINE_SELECT, RUBBER_ZOOM, PAN
};
enum VIEW {
	VIEW_XY, VIEW_XZ, VIEW_YZ, VIEW_3D
};
enum {
	NOARROW, SARROW, DARROW
};
const int NORMGRID[3] = {
	IUZ, IUY, IUX
};
enum {
	SURFACE, UNIAXIAL, BIAXIAL
};
enum {
	EXCEL, ACCESS
};
enum {
	MODELBLANK, MODELGRID, MODELBEAM, MODEL2DTRUSS, MODEL3DTRUSS, MODEL2DFRAME, MODEL3DFRAME, MODELSLAB, MODELWALL, MODELSHELL
};
enum {
	MAJOR, MINOR, AXIAL
};
/*
Numeric Edit
*/
class CNumericEdit : public CEdit {
public:
	CNumericEdit();
	CNumericEdit(LPCTSTR fmt);
	void SetWindowText(DOUBLE val, LPCTSTR fmt = NULL);
	void SetWindowText(LPCTSTR val);
protected:
	LPCTSTR deffmt;
};
/*
Gird control
*/
class GRIDCTRLDATA {
public:
	int type;
	int data_type;
	void* pdata;
	CList<CString,CString&> list;
	enum {NONE,EDIT,COMBO,COLOR};
	enum {TINT,TDBL,TDBL_F,TTEXT,TBOOL,TINDEX};
public:
	GRIDCTRLDATA() {
		type = NONE;
		pdata = NULL;
		data_type = TTEXT;
	}
	void operator = (const GRIDCTRLDATA& right) {
		type = right.type;
		data_type = right.data_type;
		pdata = right.pdata;
		if(type == COMBO) {
			CString v;
			list.RemoveAll();
			POSITION pos = right.list.GetHeadPosition();
			while(pos) {
				v = right.list.GetNext(pos);
				list.AddTail(v);
			}
		}
	}
};

class CGridCtrl : public CListCtrl {
public:
	int type;
	CWnd* pControl;
	int selItem;
	int selSubItem;
	int MyItemCount;
	enum {TYPE1,TYPE2,TYPE3};
public:
	CGridCtrl() {
		pControl = NULL;
		MyItemCount = 0;
	}
	CGridCtrl(int tp) {
		type = tp;
		pControl = NULL;
		MyItemCount = 0;
	}
	~CGridCtrl();
	void GetData();
	void PreSubclassWindow();
	void SetTwoHeaders();
	void InsertMyItem(int type,int data_type,void* pdata,CString,CString* = 0);
	int InsertItem(const LVITEM* pItem);
	int InsertItem(const LVITEM* pItem,GRIDCTRLDATA*);
	afx_msg void OnCustomDraw(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg BOOL OnClick(LPNMHDR pnmhdr, LRESULT* /*pResult*/);
    DECLARE_MESSAGE_MAP(); 
};
/*
Interaction chart
*/
class CChartDia;

class CChart : public CStatic {
public:
	CChartDia* parent;
public:
	afx_msg void OnMouseMove(UINT,CPoint);
	DECLARE_MESSAGE_MAP();
};

class CChartDia:public CDialog {
public:
	CPoint  Center;
	CPoint  Corner;
	DOUBLE  Scale;
	RPoint  eye; 
	int type;
	UINT isection;
	UINT WDIV;
	UINT QDIV;
	DESIGN* design;
	RPOINTLIST** allrplist;
	RPOINTLIST* rplist;
public:
	void Fill(DESIGN*);
	void DrawGrid(CDC*,int);
	void DrawChart(BOOL = FALSE);
	void DrawOrigin(CDC*);
	void DrawArrow(CDC* pDC,RPoint& rp1,RPoint& rp2);
	CPoint RPtoCP(RPoint);
	RPoint CPtoRP(CPoint);
public:
	CChart c_Drawing;
	CEdit c_Index;
	CEdit c_Current;
public:
	CChartDia(CWnd* parent,DESIGN*);
	~CChartDia();
	afx_msg void OnPlot(UINT nID);
	afx_msg void OnPaint();
	void OnCoordinate(CPoint&);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Template Define Database class
*/
template <class T>
class CDefineDia : public CDialog {
public:
	T* pT;
    CList<T,T&>* list;
	CListBox c_List;
	CButton c_Delete;
	CButton c_Modify;
	CButton c_Import;
	CButton c_Export;
	CString title;
	BOOL import;
	void* other;
	void* other1;
	BOOL modify;
public:
	CDefineDia(CWnd* parent,CList<T,T&>*,CString,void* = NULL,void* = NULL,BOOL = FALSE);
	void Clear();
protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnOK();
	afx_msg void OnCancel();
	afx_msg void OnAdd();
	afx_msg void OnDelete();
	afx_msg void OnModify();
	afx_msg void OnImport();
	afx_msg void OnExport();
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP();
};
/*
Template Select Database class
*/
template <class T>
class CSelectDia : public CDialog {
public:
	CList<CString,CString&> selected_list;
    CList<T,T&>* list;
	CListBox c_List;
	CString title;
public:
	CSelectDia(CWnd* parent,CList<T,T&>*,CString);
protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnOK();
};
/*
Group property Dialog
*/
class CGroupPropertyDia:public CDialog {
public:
    GROUP* group;
	CEdit c_Name;
	CListBox c_JList;
	CListBox c_MList;
	CListBox c_SList;
public:
	CGroupPropertyDia(CWnd* parent,GROUP*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Curved member
*/
class CCurvedMemberDia : public CDialog {
public:
	int nDiv;
	int type;
	BOOL single;
	RPoint rp;
	CComboBox c_Type;
	enum {STRAIGHT,CIRCULAR,PARABOLIC};
public:
	CCurvedMemberDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
};
/*
New model Dialog
*/
class CNewDia : public CDialog {
public:
	int type;
public:
	CNewDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Frame Dialog
*/
class CmyDocument;
class CmyView;

class CNewFrameDia : public CDialog {
public:
	int type;
	int frame_type;
	int N[3];
	DOUBLE D[3];
	DOUBLE gap;
	CEdit c_Cover;
	CStatic c_GapText;
	CEdit c_Gap;
	BOOL restraint;
	CComboBox c_FrameType;
	enum {PORTAL,BRACED_CONCENTRIC,BRACED_ECCENTRIC};
	enum {OPEN_FRAME,PERIMETER_FRAME,BEAM_SLAB,FLAT_PLATE};
public:
	CNewFrameDia(CWnd* parent,int);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP();
};
/*
Slab Dialog
*/
class CNewSlabDia : public CDialog {
public:
	int N[2];
	DOUBLE D[2];
	DOUBLE M[2];
	BOOL restraint;
public:
	CNewSlabDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP();
};
/*
Wall Dialog
*/
class CNewWallDia : public CDialog {
public:
	int N[2];
	DOUBLE D[2];
	BOOL restraint;
public:
	CNewWallDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP();
};
/*
Shell Dialog
*/
class CNewShellDia : public CDialog {
public:
	CComboBox c_Type;
	int type;
	int N[2];
	DOUBLE D[2];
	BOOL restraint;
	enum {CYLINDER};
public:
	CNewShellDia(CWnd* parent);
	afx_msg void OnUpdateControls();
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP();
};
/*
Beam Dialog
*/
class CNewBeamDia : public CDialog {
public:
	int nspan;
	DOUBLE spanlength;
	BOOL restraint;
public:
	CNewBeamDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual void OnOK();
};
/*
2D truss Dialog
*/
class CNewTruss2DDia : public CDialog {
public:
	/*for pratt truss*/
	DOUBLE L1;
	DOUBLE L2;
	int N1;
	int N2;
	DOUBLE H1;
	DOUBLE H2;
	DOUBLE H3;
	/*other*/
	int ndiv;
	DOUBLE divlength;
	DOUBLE height;
	/*other*/
	enum {SLOPED,VERTICAL,PRATT};
	BOOL restraint;
	int type;
	CComboBox c_Type;
	CEdit c_Cover1;
	CEdit c_Cover2;
public:
	CNewTruss2DDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	afx_msg void OnType();
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP();
};
/*
3D truss
*/
class TOWER {
public:
    DOUBLE d[4];
public:
	TOWER() {
		for(int i = 0;i < 4;i++)
			d[i] = 0;
	}
	friend int operator == (const TOWER& left,const TOWER& right) {
		return (left.d[0] == right.d[0]);
	}
};

#define TOWERLIST     CList<TOWER,TOWER&>
#define PTOWERLIST    CList<TOWER,TOWER&>*

class CNewTruss3DDia : public CDialog {
public:
	/*for roof truss*/
	int Nx;
	int Ny;
	int Ndx;
	int Ndy;
	DOUBLE Sx;
	DOUBLE Sy;
	DOUBLE H;
	/*tower*/
	TOWERLIST tower;
	/*other*/
	BOOL restraint;
	int type;
	int tower_shape;
	CComboBox c_Type;
	CComboBox c_Tower;
	CEdit c_Cover1;
	CEdit c_Cover2;
	enum {ROOF,TELE_TOWER};
	enum {RECTANGULAR, TRIANGULAR};
	CListBox c_List[4];
	CEdit c_Value[4];
public:
	CNewTruss3DDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	afx_msg void OnType();
	afx_msg void OnUpdateControls();
	afx_msg void OnAdd();
	afx_msg void OnMove();
	afx_msg void OnDelete();
	afx_msg void OnDeleteAll();
	afx_msg void OnUpdateLocTh();
	afx_msg void OnUpdateLocTw();
	afx_msg void OnUpdateLocTah();
	afx_msg void OnUpdateLocTaw();
	afx_msg void OnUpdateLoc(UINT);
	afx_msg void OnTower();
	DECLARE_MESSAGE_MAP();
};
/*
Grid Dialog
*/
class CSystemPage:public CPropertyPage {
public:
	SYSTEM* system;
	BOOL cartesian;
	DOUBLE grids[3];
	int ngrid[3];
public:
	CSystemPage(UINT,SYSTEM*);
protected:
	virtual void OnOK();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
property sheet
*/
class CNewSystemSht : public CPropertySheet {
public:
    CSystemPage cartesian;
    CSystemPage radial;
public:
    CNewSystemSht(CWnd*, SYSTEM*);
};
/*
System Origin Dialog
*/
class CSystemOriginDia:public CDialog {
public:
    SYSTEM* system;
	RPoint origin;
	RPoint rotation;
public:
	CSystemOriginDia(CWnd* parent,SYSTEM*);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual void OnOK();
};
/*
Edit Grid Dialog
*/
class CSystemPropertyDia:public CDialog {
public:
	SYSTEM system;
	SYSTEM* psystem;
	CListBox c_List;
	CEdit c_Location;
	CEdit c_Name;
	CButton c_Z;
	CButton c_Y;
	CButton c_X;
	CButton c_Origin;
	int selected;
public:
	CSystemPropertyDia(CWnd* parent,SYSTEM*);
	afx_msg void OnUpdateLoc();
	afx_msg void OnUpdateControls();
	afx_msg void OnAddGrid();
	afx_msg void OnMoveGrid();
	afx_msg void OnDeleteGrid();
	afx_msg void OnDeleteAll();
	afx_msg void OnNew();
	afx_msg void OnOrigin();
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Function Property Dialog
*/
class CFunctionDrawing : public CStatic {
public:
	CListBox c_Lx;
    CListBox c_Ly;
	CEdit c_Ex;
	CEdit c_Ey;
	CButton c_Add;
	CButton c_Move;
	CButton c_Delete;
	CButton c_DeleteAll;
	CEdit c_FunctionValue;
	FUNCTION* function;
	DOUBLE ScaleX;
	DOUBLE section;
	int type;
	enum {VIEWONLY,EDITABLE};
public:
	CFunctionDrawing(int type);
	void Plot(BOOL);
	void Init();
	void AcceptMessage(UINT nID);
	afx_msg void OnUpdateLoc(UINT nID);
	afx_msg void OnUpdateControls();
	afx_msg void OnPaint();
	afx_msg void OnMouseMove(UINT,CPoint);
	afx_msg void OnAdd();
	afx_msg void OnMove();
	afx_msg void OnDelete();
	afx_msg void OnDeleteAll();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};

class CFunctionPropertyDia:public CDialog {
public:
	FUNCTION function;
	FUNCTION* pfunc;
	CEdit c_Name;
	CFunctionDrawing c_Drawing;
	CComboBox c_Type;
	CEdit c_Amplitude;
    CEdit c_Period;
	CEdit c_Cycles;
	CEdit c_Divisions;
	CEdit c_Ag;
	CEdit c_Cover1;
	CEdit c_Cover2;
	CButton c_Browse;
	CEdit c_File;
	CComboBox c_Soil;
public:
	CFunctionPropertyDia(CWnd* parent,FUNCTION*);
	void RelayMessage(UINT nID);
	afx_msg void OnUpdateControls();
	afx_msg void OnType();
	afx_msg void OnPaint();
	afx_msg void OnBrowse();
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Function Plot Dialog
*/
class CFunctionPlotDia:public CDialog {
public:
	FUNCTION* function;
	CFunctionDrawing c_Drawing;
	CEdit c_FuncValue;
public:
	CFunctionPlotDia(CWnd* parent,FUNCTION*);
	afx_msg void OnPaint();
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Result Plot dialog
*/
class CResultPlotDia:public CDialog {
public:
	FUNCTION function;
	CComboBox c_Dir;
	CComboBox c_Type;
	CListBox c_JointL;
	CEdit c_Section;
	CEdit c_Cover1;
	enum {JOINT_DISPLACEMENT,JOINT_FORCE,MEMBER_FORCE};
public:
	CResultPlotDia(CWnd* parent);
protected:
	void OnDisplay();
	afx_msg void OnUpdateControls();
	afx_msg void OnType();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Loadcase property Dialog
*/
class CLoadCasePropertyDia:public CDialog {
public:
    LOADCASE* loadcase;
	CComboBox c_List;
	CEdit c_Name;
public:
	CLoadCasePropertyDia(CWnd* parent,LOADCASE*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Modalcase property Dialog
*/
class CModalCasePropertyDia:public CDialog {
public:
    MODALCASE* modalcase;
	CEdit c_Name;
public:
	CModalCasePropertyDia(CWnd* parent,MODALCASE*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Damping Dialog
*/
class CDampingDia:public CDialog {
public:
	CComboBox c_Type;
	CEdit c_Cover1;
	CEdit c_Cover2;
	DAMPING damping;
	CEdit c_W2;
	CEdit c_Constant;
	CFunctionDrawing c_Drawing;
	CEdit c_FuncValue;
public:
	CDampingDia(CWnd* parent,DAMPING*);
	void RelayMessage(UINT nID);
protected:
	afx_msg void OnUpdateControls();
	afx_msg void OnType();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Response Property Dialog
*/
class CResponsePropertyDia:public CDialog {
public:
    RESPONSE* response;
	DAMPING damping;
	PMODALCASELIST modalcases;
	PFUNCTIONLIST functions;
	CEdit c_Name;
	CComboBox c_ModalList;
	CComboBox c_Func;
	CComboBox c_Dir;
	CEdit c_Scale;
	CListBox c_DirList;
	CListBox c_FuncList;
	CListBox c_ScaleList;
	UINT rID;
public:
	CResponsePropertyDia(UINT mID,CWnd* parent,RESPONSE*,PMODALCASELIST,PFUNCTIONLIST);
protected:
	afx_msg void OnDamping();
	afx_msg void OnAdd();
	afx_msg void OnMove();
	afx_msg void OnDelete();
	afx_msg void OnDeleteAll();
	void OnUpdateLoc(UINT nID);
	afx_msg void OnUpdateLocDir();
    afx_msg void OnUpdateLocFunc();
	afx_msg void OnUpdateLocScale();
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Loadcombo property Dialog
*/
class CLoadComboPropertyDia:public CDialog {
public:
    COMB_TYPE* combo;
	PANALYSISCASELIST loadcase[CASETYPES];
	CComboBox c_CaseList;
	CComboBox c_Type;
	CListBox c_List;
	CEdit c_Name;
	CEdit c_FS;
	int type;
	enum {TCOMBO,TNL,TBUCKLING};
public:
	CLoadComboPropertyDia(CWnd* parent,UINT,COMB_TYPE*,PANALYSISCASELIST*,int);
	afx_msg void OnAdd();
	afx_msg void OnMove();
	afx_msg void OnDelete();
	afx_msg void OnDeleteAll();
	afx_msg void OnUpdateLoc();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};

class CCombinationPropertyDia : public CLoadComboPropertyDia {
public:
    CCombinationPropertyDia(CWnd* parent,COMB_TYPE*,PANALYSISCASELIST*);
};
class CNLcasePropertyDia : public CLoadComboPropertyDia {
public:
	NLCASE* pnl_case;
public:
    CNLcasePropertyDia(CWnd* parent,NLCASE*,PANALYSISCASELIST*);
	virtual void DoDataExchange(CDataExchange* pDX);
};
class CBucklingPropertyDia : public CLoadComboPropertyDia {
public:
	BUCKLINGCASE* pbkl_case;
public:
    CBucklingPropertyDia(CWnd* parent,BUCKLINGCASE*,PANALYSISCASELIST*);
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Rigid body property Dialog
*/
class CConstraintPropertyDia:public CDialog {
public:
	int dof[6];
    CONSTRAINT constraint;
	CONSTRAINT* pconst;
	CListBox c_List;
    CComboBox c_Type;
	CComboBox c_System;
	CEdit c_Cover;
	CButton c_Check[6];
	CButton c_Radio[4];
    PSYSTEMLIST systems;
public:
	CConstraintPropertyDia(CWnd* parent,CONSTRAINT*,PSYSTEMLIST);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnUpdateControls();
	afx_msg void OnType();
	DECLARE_MESSAGE_MAP();
};
/*
Joint Restraint Dialog
*/
class CJointRestraintDia:public CDialog {
public:
	BOOL res[6];
	BOOL con[6];
public:
	CJointRestraintDia(CWnd* parent,JOINT* pjoint);
	afx_msg void OnFastRestraint(UINT nID);
protected:
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Joint Load Dialog
*/
class CJointLoadDia:public CDialog {
public:
	LOADCASE* ploadcase;
	SYSTEM* system;
	PLOADCASELIST loadcases;
	PSYSTEMLIST systems;
	JOINT* joint;
	int type;
	DOUBLE Q[6];
	int mode;
	CComboBox c_List;
	CComboBox c_System;
	CEdit c_Cover;
public:
	CJointLoadDia(CWnd* parent,PLOADCASELIST,PSYSTEMLIST,JOINT* pjoint,int,LOADCASE*);
protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual void OnOK();
};
/*
Joint Local Axis Dialog
*/
class CJointAxisDia:public CDialog {
public:
	JOINT* joint;
	RPoint rotation;
public:
	CJointAxisDia(CWnd* parent,JOINT* pjoint);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual void OnOK();
};
/*
Joint Load Dialog
*/
class CJointReactionsDia:public CDialog {
public:
	int type;
	JOINT* joint;
	CGridCtrl c_List;
public:
	CJointReactionsDia(CWnd* parent,JOINT* pjoint,int);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
};
/*
Member Local Axis Dialog
*/
class CAxisRotationDia : public CDialog {
public:
	DOUBLE alpha;
	CString title;
public:
	CAxisRotationDia(CWnd* parent,DOUBLE,CString);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
};
/*
Member Load Dialog
*/
enum {
	GRAVITYD = 3
};
class CMemberLoadDia:public CDialog {
public:
	LOADCASE* ploadcase;
	PLOADCASELIST loadcases;
	SYSTEM* system;
	PSYSTEMLIST systems;
	DOUBLE d[4];
	DOUBLE Q[4];
	DOUBLE Qu;
	int mode;
	int dir;
	int loadtype;
	CComboBox c_List;
	CComboBox c_LoadDir;
	CComboBox c_System;
	CEdit c_Cover;
	BOOL use_relative;
	enum {TFORCE,TMOMENT};
	BOOL distributed;
public:
	CMemberLoadDia(CWnd* parent,PLOADCASELIST,PSYSTEMLIST,BOOL,LOADCASE*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP();
};
/*
Slab Distributed Load Dialog
*/
class CSlabDistLoadDia:public CDialog {
public:
	LOADCASE* ploadcase;
	PLOADCASELIST loadcases;
	SYSTEM* system;
	PSYSTEMLIST systems;
	DOUBLE Qu;
	int mode;
	int dir;
	int loadtype;
	CComboBox c_List;
	CComboBox c_System;
	CComboBox c_LoadDir;
	enum {TFORCE,TMOMENT};
public:
	CSlabDistLoadDia(CWnd* parent,PLOADCASELIST,PSYSTEMLIST,LOADCASE*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP();
};
/*
Strain loading
*/
class CStrainLoadDia:public CDialog {
public:
	LOADCASE* ploadcase;
	PLOADCASELIST loadcase;
	DOUBLE Value;
	int type;
	int dir;
	int mode;
	int stype;
	CComboBox c_List;
	CComboBox c_Type;
	enum {TEMP_STRAIN,FAB_STRAIN};
	enum {TMEMBER,TSLAB};
public:
	CStrainLoadDia(CWnd* parent,PLOADCASELIST,LOADCASE*,UINT);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Reinforcement Dialog
*/
class CRebarDia:public CDialog {
public:
    REBAR* rebar;
public:
	CButton c_Rectangular;
	CButton c_Circular;
	CButton c_Coulmn;
	CEdit c_Ny;
	CEdit c_Nz;
	CEdit c_Total;
	CEdit c_BarSize;
	CEdit c_SBarSize;
	CRebarDia(CWnd* parent,REBAR*);
protected:
	afx_msg void OnUpdateControls();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Material Property Dialog
*/
class CMaterialPropertyDia:public CDialog {
public:
    MATERIAL* material;
	CComboBox c_List;
	CEdit c_Elasticity;
	CEdit c_Poisson;
	CEdit c_Thermal;
	CEdit c_Fctk;
	CEdit c_Fck;
	CEdit c_Fyks;
	CEdit c_Fyk;
	CEdit c_Ftks;
	CEdit c_Ftk;
	CEdit c_Es;
	CEdit c_Name;
public:
	CMaterialPropertyDia(CWnd* parent,MATERIAL*);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP()
};
/*
Section Stiffness Modifiers dialog
*/
class CSectionModifiersDia : public CDialog {
public:
	SECTION* section;
	CGridCtrl c_Table;
public:
	CSectionModifiersDia(CWnd* parent,SECTION*);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Section Property Dialog
*/
class CSectionPropertyDia:public CDialog {
public:
	PMATERIALLIST materials;
    SECTION* section;
	DETAILING* detailing;
	CComboBox c_List;
	CComboBox c_TypeList;
	CStatic c_Drawing;
	CButton c_Rebar;
	CButton c_Chart;
	CEdit c_Name;
	CEdit c_Radius;
	CEdit c_Width;
	CEdit c_Height;
	CEdit c_Tw;
	CEdit c_Th;
	CEdit c_Bd;
	CNumericEdit c_A;
	CNumericEdit c_Ay;
	CNumericEdit c_Az;
    CNumericEdit c_Zy;
	CNumericEdit c_Zz;
	CNumericEdit c_Sy;
	CNumericEdit c_Sz;
	CNumericEdit c_ry;
	CNumericEdit c_rz;
	CNumericEdit c_Ix;
	CNumericEdit c_Iy;
	CNumericEdit c_Iz;
public:
	CSectionPropertyDia(CWnd* parent,SECTION*,PMATERIALLIST);
protected:
	afx_msg void OnEditMaterial();
    afx_msg void OnEditRebar();
	afx_msg void OnSetModifiers();
	afx_msg void OnChart();
	afx_msg void OnPaint();
	afx_msg void OnUpdateControls();
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP()
};
/*
Section Property Dialog
*/
class CASectionPropertyDia:public CDialog {
public:
	PMATERIALLIST materials;
    ASECTION* section;
	CComboBox c_List;
	CComboBox c_TypeList;
	CComboBox c_StateTypeList;
	CButton c_Rebar;
	CEdit c_Name;
	CEdit c_Height;
public:
	CASectionPropertyDia(CWnd* parent,ASECTION*,PMATERIALLIST);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnEditMaterial();
    afx_msg void OnEditRebar();
	afx_msg void OnUpdateControls();
	DECLARE_MESSAGE_MAP()
};
/*
Joint Information Dialog
*/
class CJointInfoDia:public CDialog {
public:
	JOINT* joint;
	PLOADCASELIST loadcase;
	CListBox c_FList;
	CGridCtrl c_Jpoint;
public:
	CJointInfoDia(CWnd* parent,JOINT* pjoint,PLOADCASELIST);
protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnEditRestraint();
	afx_msg void OnEditAxis();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP();
};
/*
Member Information Dialog
*/
class CMemberInfoDia:public CDialog {
public:
	MEMBER* member;
	PSECTIONLIST sections;
	PMATERIALLIST materials;
	PLOADCASELIST loadcases;
	CListBox c_FList;
public:
	CMemberInfoDia(CWnd* parent,MEMBER* pmember,
		PSECTIONLIST psection,
		PMATERIALLIST pmaterial,
		PLOADCASELIST ploadcase);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnEditSection();
    afx_msg void OnEditRelease();
	afx_msg void OnEditJoint(UINT nID);
	afx_msg void OnEditAxis();
	afx_msg void OnEditDivisions();
	DECLARE_MESSAGE_MAP();
};
/*
Member break dialog
*/
class CMemberBreakDia:public CDialog {
public:
	UINT breakn;
	DOUBLE breakr;
public:
	CMemberBreakDia(CWnd* parent);
protected:
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Member Tapered Section
*/
class CTaperedSectionDia:public CDialog {
public:
	MEMBER* member;
	SECTION* s_sec;
	SECTION* e_sec;
	int EIyy;
	int EIzz;
	DOUBLE soff;
	DOUBLE eoff;
	PSECTIONLIST sections;
	PMATERIALLIST materials;
	CComboBox c_Start;
    CComboBox c_End;
	CComboBox c_EIyy;
	CComboBox c_EIzz;
public:
	CTaperedSectionDia(CWnd* parent,PSECTIONLIST,PMATERIALLIST,MEMBER* = 0);
protected:
	virtual void OnOK();
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	afx_msg void OnEditSection(UINT nID);
	DECLARE_MESSAGE_MAP();
};

/*
Slab Information Dialog
*/
class CSlabInfoDia:public CDialog {
public:
	SLAB* slab;
	PASECTIONLIST section;
	PMATERIALLIST material;
	PLOADCASELIST loadcase;
	CComboBox c_List;
	CListBox c_FList;
public:
	CSlabInfoDia(CWnd* parent,SLAB* pslab,
		PASECTIONLIST psection,
		PMATERIALLIST pmaterial,
		PLOADCASELIST ploadcase);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnEditSection();
	afx_msg void OnEditAxis();
	afx_msg void OnEditDivisions();
	DECLARE_MESSAGE_MAP();
};
/*
Internal Member force diagrams
*/
class CForceDiagramDia:public CDialog {
public:
	int dir;
	int fill;
	int value;
public:
	CForceDiagramDia(CWnd* parent);
protected:
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Slab Force diagrams
*/
class CSlabForceDiagramDia:public CDialog {
public:
	int dir;
	int average;
	CComboBox c_Type;
	CComboBox c_Stress;
public:
	CSlabForceDiagramDia(CWnd* parent);
protected:
	virtual void OnOK();
	virtual BOOL OnInitDialog();
    afx_msg void OnUpdateControls();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Member Releases Dialog
*/
class CMemberReleaseDia:public CDialog {
public:
	int r[12];
	CButton c_R[12];
public:
	CMemberReleaseDia(CWnd*,MEMBER*);
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnUpdateControls(UINT nID);
	DECLARE_MESSAGE_MAP();
};
/*
Member Divisons
*/
class CMemberDivisionsDia:public CDialog {
public:
	int nMinDiv;
	int nMinFrameDiv;
	BOOL DivAtInterim;
public:
	CMemberDivisionsDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Slab Divisons
*/
class CSlabDivisionsDia:public CDialog {
public:
	int nDivx;
	int nDivy;
	BOOL constraint;
public:
	CSlabDivisionsDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
};

/*
Select Tables Dialog
*/
class CmyTreeCtrl : public CTreeCtrl {
public:
	afx_msg void OnLButtonDown(UINT,CPoint);
	void SelectSiblings(HTREEITEM,BOOL);
    void SelectParents(HTREEITEM,BOOL);
	DECLARE_MESSAGE_MAP();
};
class CSelectTableDia:public CDialog {
public:
	PANALYSISCASELIST* analysis_cases;
	CListBox c_Cases;
	CmyTreeCtrl c_Table;
	HTREEITEM  m_rghItem[NTABLES];
	BOOL state[NTABLES];
//*/
	enum {
		     Model
			,MGeneral
			,GSettings
			,GCoordinate_systems
			,Property_definitions
			,DMaterial
			,DSection
			,Function_definitions
			,LR_history
			,LR_spectrum
			,Joint
			,JCoordinates
			,JRestraints
			,JLoads
			,JAddedMass
			,JLocal_axis
			,JConstraints
			,Member
			,MConnectivity
			,MSection
			,MLoads
			,MStrain_loads
			,MLocal_axis
			,MRelease
			,MDivisions
			,Slab
			,SConnectivity
			,SSection
			,SLoads
			,SLocal_axis
			,Analysis_cases
			,ADefinitions
			,ALoad_cases
			,ANonLinear_cases
            ,AModal_cases
			,ABuckling_cases
			,ARH_cases
			,ARS_cases
			,ACombinations
			,Analysis
			,AGeneral
			,GModes
			,GBucklingModes
			,AJoint
			,JReactions
			,JDisplacements
			,JMasses
			,AMember
			,MInternal_forces
			,ASlab
			,SMoments
			,Design
			,DResult
	};
	//*/
public:
	CSelectTableDia(CWnd* parent,PANALYSISCASELIST*);
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	ANALYSISCASE* find_case();
};
/*
Table view
*/
class CTableViewDia : public CDialog {
public:
	CGridCtrl c_Table;
    CComboBox c_List;
	CEdit c_Index;
	CEdit c_Max;
	CmyDatabase* pdatabase;
public:
    CTableViewDia(CWnd* parent,CmyDatabase*);
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnTableChange();
	afx_msg void OnSelectChange(LPNMHDR pnmhdr, LRESULT* /*pResult*/);
	DECLARE_MESSAGE_MAP();
};
/*
Analysis Dialog
*/
class CAnalysisDia:public CDialog {
public:
	PANALYSISCASELIST* analysis_cases;
	CListBox c_List;
	CStatic c_Action;
	CStatic c_Status;
	CComboBox c_Stiffness;
public:
	CAnalysisDia(CWnd*,PANALYSISCASELIST*);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	afx_msg void OnRun();
	afx_msg void OnUpdateLoc();
	afx_msg void OnChangeAction();
	afx_msg void OnUpdateStiffness();
	ANALYSISCASE* find_case();
	DECLARE_MESSAGE_MAP();
};
/*
Timer update dialog
*/
class CmyView;

class CTimerUpdateDia : public CDialog {
public:
	UINT_PTR myTimer;
	MEMBER* member;
	CmyView* pView;
	BOOL blink;
public:
	CTimerUpdateDia(CWnd*,UINT,MEMBER*,CmyView*);
	virtual BOOL OnInitDialog();
	afx_msg void OnTimer(UINT);
	DECLARE_MESSAGE_MAP();
};
/*
Analysis Update Dialog
*/
class CmyListBox : public CListBox {
public:
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDIS);
};

class CAnalysisProgressDia:public CDialog {
public:
	CmyListBox c_List;
	CEvent* pEvent;
	CTime beginTime;
	CStatic c_Begin;
	CStatic c_End;
	CStatic c_Elapsed;
	CButton c_Ok;
    CButton c_Cancel;
	UINT_PTR myTimer;
public:
	CAnalysisProgressDia(CWnd*,CEvent*);
	void AnalysisEnd();
	afx_msg void OnTimer(UINT);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	DECLARE_MESSAGE_MAP();
};
/*
Display mode
*/
class DISPLAY : public CObject{
	DECLARE_SERIAL(DISPLAY)
public:
	BOOL JointLabel;
	BOOL JointRestraint;
	BOOL JointShow;
	BOOL JointLoad;
	BOOL JointAxis;
	BOOL MemberLabel;
	BOOL MemberSection;
    BOOL MemberMaterial;
	BOOL MemberLoad;
	BOOL MemberAxis;
	BOOL MemberRelease;
	BOOL SlabLabel;
	BOOL SlabSection;
    BOOL SlabMaterial;
	BOOL SlabLoad;
	BOOL SlabAxis;
	BOOL Origin;
	BOOL Grid;
	BOOL FillDiagram;
	BOOL ShowValues;
	BOOL Constraints;
	BOOL Numbering;
	BOOL Extrude;
	DISPLAY() {
        JointLabel = FALSE; 
        JointRestraint = TRUE;  
		JointShow = FALSE;
		JointLoad = TRUE;
		JointAxis = FALSE;
		MemberLabel = FALSE;
		MemberSection = FALSE;
		MemberMaterial = FALSE;
		MemberLoad = TRUE;
		MemberAxis = FALSE;
		MemberRelease = TRUE;
		SlabLabel = FALSE;
		SlabSection = FALSE;
		SlabMaterial = FALSE;
		SlabLoad = TRUE;
		SlabAxis = FALSE;
		Origin = TRUE;
		Grid = TRUE;
		FillDiagram = TRUE;
		Constraints = TRUE;
		Numbering = FALSE;
		ShowValues = FALSE;
		Extrude = FALSE;
	}
	void operator = (const DISPLAY& right) {
		memcpy(this,&right,sizeof(right));
	}
	void Serialize(CSerializer& ar);
};
/*
Display Options Dialog
*/
class CDisplayOptionsDia:public CDialog {
public:
	DISPLAY mydisplay;
	BOOL ApplyToAll;
public:
	CDisplayOptionsDia(CWnd* parent,DISPLAY* pdisplay);
protected:
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Color Options Dialog
*/
class CColorsDia:public CDialog {
public:
	CGridCtrl c_List1;
    CGridCtrl c_List2;
    CGridCtrl c_List3;
public:
	CColorsDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
};
/*
Paste Dialog
*/
class CPasteDia:public CDialog {
public:
	RPoint v;
public:
	CPasteDia(CWnd* parent);
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Result Dialog
*/
class CResultDia:public CDialog {
public:
	UINT Total;
	int selected;
	CComboBox c_List;
	int type;
	enum {NORMAL,MINMAX};
public:
	CResultDia(CWnd* parent,UINT,UINT,int = NORMAL);
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
};
/*
Select Labels Dialog
*/
class CmyDocument;

class CSelectLabelsDia:public CDialog {
public:
	CComboBox c_Type;
	CListBox c_List;
	BOOL select;
	int minid;
	int maxid;
	int inc;
	CEdit c_Min;
	CEdit c_Max;
	CEdit c_Inc;
	enum {JOINTT,MEMBERT,SLABT};
public:
	CSelectLabelsDia(CWnd* parent,BOOL);
	afx_msg void OnType();
	afx_msg void OnSelect(UINT nID);
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Display Analysis Options Dialog
*/
class CAnalysisOptionsDia:public CDialog {
public:
	int deflection;
	int divisions;
	BOOL sdof[6];
	BOOL ddof[3];
	UBMP8* sdofs;
	UBMP8* ddofs;
	CEdit c_Div;
public:
	CAnalysisOptionsDia(CWnd* parent,UBMP8*,UBMP8*);
	afx_msg void OnFastDof(UINT nID);
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Display Design Options Dialog
*/
class CDesignOptionsDia:public CDialog {
public:
	DETAILING detail;
	DETAILING* pdetail;
	CListBox c_BarList;
	CListBox c_List[4];
	CComboBox c_Type;
    CComboBox c_Layers;
	CGridCtrl c_Options;
	enum {SHEART,BEAMT,COLUMNT,SLABT};
public:
	CDesignOptionsDia(CWnd* parent,DETAILING*);
	afx_msg void OnType();
	afx_msg void OnAdd();
    afx_msg void OnDelete();
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Misc Options Dialog
*/
class CmyDocument;

class CMiscOptionsDia:public CDialog {
public:
	BOOL AnalysisResult;
	int deflection;
	DOUBLE scale_deflection;
	DOUBLE scale_stress;
	DOUBLE scale_load;
public:
	CMiscOptionsDia(CWnd* parent);
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Stress Diagram Dialog
*/
class CStressDiagramDia;

class CStressControl : public CStatic {
public:
	CStressDiagramDia* parent;
public:
	void UpdateParent(CPoint&);
	afx_msg void OnLButtonDown(UINT,CPoint);
    afx_msg void OnMouseMove(UINT,CPoint);
	DECLARE_MESSAGE_MAP();
};
class CStressDiagramDia : public CTimerUpdateDia {
public:
	int dir;
	int scroll;
	MEMBER* member;
	CComboBox c_Dir;
    CComboBox c_Scroll;
	ANALYSISCASE* analysis_case;
	CStressControl c_Drawing[3];
	CStressControl c_Value[3];
	CEdit c_Section;
	DOUBLE xscale;
	DOUBLE yscale;
	int Width;
	int Height;
	DOUBLE section;
	int deflection_type;
	static const int OFFSET;
	enum {SCROLL,SHOWMAX,SHOWMIN};
	enum {DABSOLUTE,DRELATIVE};
public:
	void DrawDiagram(int,int,BOOL);
	void UpdateDrawing(BOOL = FALSE);
	void DrawRectangle(CDC*,DOUBLE,DOUBLE,DOUBLE,DOUBLE);
	CStressDiagramDia(CWnd* parent,MEMBER*,ANALYSISCASE*,CmyView*);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg void OnChangeDir();
	afx_msg void OnScroll();
	afx_msg void OnSection();
	afx_msg void OnDeflectionType(UINT nID);
	DECLARE_MESSAGE_MAP();
};
/*
Design Overwrites dialog
*/
class CDesignOverwritesDia : public CDialog {
public:
	MEMBER* member;
	CGridCtrl c_Table;
public:
	CDesignOverwritesDia(CWnd* parent,MEMBER*);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP();
};
/*
Design result dialog
*/
class CDesignResultDia : public CTimerUpdateDia {
public:
	MEMBER* member;
	ANALYSISCASE* analysis_case;
	CGridCtrl c_Table;
public:
	CDesignResultDia(CWnd* parent,MEMBER*,ANALYSISCASE*,CmyView*);
	virtual void DoDataExchange(CDataExchange* pDX);
	virtual BOOL OnInitDialog();
	afx_msg void OnDesignOverwrites();
	DECLARE_MESSAGE_MAP();
};
/*
Animation Dialog
*/
class CAnimationDia:public CDialog {
public:
	CEdit c_Time;
	BOOL sound;
	UINT interval;
	BOOL stop;
	BOOL finished;
	CmyView* pView;
	CButton c_Redo;
public:
	CAnimationDia(CWnd* parent,CmyView*);
	void Finish();
public:
	virtual void OnOK();
	afx_msg void OnRedo();
	afx_msg void OnSpeed(UINT nID);
	afx_msg void OnSound();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();
};
/*
Document
*/
class UndoData {
public:
	BOOL adding;
	int type;
	CString name;
	JOINTLIST joints;
	enum {TJOINT,TMEMBER,TSLAB};
	enum {ADD, REMOVE};
public:
	UndoData() {
	}
	UndoData(UndoData& right) {
		*this = right;
	}
	void operator = (const UndoData& right) {
		adding = right.adding;
		type = right.type;
		name = right.name;
		JOINT v;
		joints.RemoveAll();
		POSITION pos = right.joints.GetHeadPosition();
		while(pos) {
			v = right.joints.GetNext(pos);
			joints.AddTail(v);
		}
	}
};

class CUndoRedo {
public:
	int current;
	int count;
	CList<UndoData,UndoData&> data;
public:
	CUndoRedo();
    void Do(CmyDocument*,int);
	void Save(JOINT*,BOOL = TRUE);
    void Save(MEMBER*,BOOL = TRUE);
	void Save(SLAB*,BOOL = TRUE);
	void Clear();
	void Add(UndoData&);
};

class CmyDocument : public CDocument {
	DECLARE_DYNCREATE(CmyDocument)
public:
	static CmyDocument* Current_Document;

	SYSTEM* global;
	SYSTEMLIST systems;
	/*objects*/
	JOINTLIST joints;
	MEMBERLIST members;
	SLABLIST slabs;
	PENTITYLIST entities[NENTITIES];
	/*removed objects*/
	MEMBERLIST removed_members;
	SLABLIST removed_slabs;
	/*meshed objects*/
	SLABPLIST_PLIST slab_elements;
	MEMBERPLIST_PLIST member_elements;
	/*others*/
	MATERIALLIST materials;
	SECTIONLIST sections;
	ASECTIONLIST asections;
	GROUPLIST groups;
	LOADCASELIST loadcases;
	MODALCASELIST modalcases;
	RESPONSECASELIST responsecases;
	RESPONSESPECLIST responsespecs;
	COMBOLIST combinations;
	NLCASELIST nlcases;
	BUCKLINGCASELIST bucklingcases;
	PANALYSISCASELIST all_analysis_cases[CASETYPES];
	CONSTRAINTLIST constraints;
	FUNCTIONLIST rhfunctions;
	FUNCTIONLIST rsfunctions;
	MEMBERLIST copyMembers;
	SLABLIST copySlabs;
	JOINTLIST copyJoints;
	BOOL cutStatus;
	UBMP8 static_dofs;
	UBMP8 dynamic_dofs;
	DETAILING detailing;
	CmyDatabase myTables;
	
	DOUBLE load_scale;
	DOUBLE force_scale[6];
	DOUBLE disp_scale;
	DOUBLE user_load_scale;
	DOUBLE user_force_scale;
	DOUBLE user_disp_scale;
	BOOL select_all;
	STATE state;
	UINT nClicks;
	RPoint selrp[20];
	BOOL Lock;
	BOOL AnalysisResult;
	BOOL DesignResult;
	BOOL DisplacementResult;
	CComboBox* LoadCaseListBox;
	ANALYSISCASE* c_AnalysisCase;
	static UINT nViews;
	CEvent threadStart;
	CAnalysisProgressDia* volatile pAnalysisProgressDia;
	BOOL Snap_To_Grid;
	BOOL Snap_To_Joint;
	BOOL Snap_To_MemberMid;
	BOOL Snap_To_Member;
    CUndoRedo UndoRedo;
public:
	CmyDocument();
	~CmyDocument();
	BOOL OnNewDocument();
	BOOL OnSaveDocument(LPCTSTR);
	void InitViews(LPCSTR);
	void DeleteContents();
	void Serialize(CArchive&);
	void Serialize(CSerializer&);
	void SerializeAnalysisResult(CSerializer&);
	void SerializeGeneral(CSerializer&);
	void DeleteSlab(SLAB*);
	SLAB* AddSlab(JOINT*,JOINT*,JOINT*,JOINT*);
	SLAB* AddSlab(JOINT**,UINT);
	SLAB* AddSlab(JOINT*,UINT);
    SLAB* AddSlab(SLAB*);
	void DeleteMember(MEMBER*);
	MEMBER* AddMember(JOINT*,JOINT*);
	void DeleteJoint(JOINT*);
	JOINT* AddJoint(JOINT*);
	ENTITY* FindEntity(CString&,int);
	void RefreshConstraints();
	void CombineLoad(COMB_TYPE*);
	void CombineModes(RESPONSESPEC*);
	void ApplyDeadLoad();
	void RemoveTemporaryLoads();
	void LumpMasses();
	void AllocMemory(int);
	void FreeMemory(int);
	void PrintOutput();
	void CalculateForces(VECTOR,VECTOR,UINT,int,BOOL = FALSE);
	void CalculateEffectiveLength(MEMBER*);
	void ShowFile(int);
	void ExportToAutocad(CString);
	void ImportFromAutocad(CString);
	void CheckCaseInCombo(int);
	void CheckCaseInResp();
	void FillCases();
	void FillAnalysisCases();
	void PrintProgress(CString,BOOL = FALSE);
	void PrintStylish(CString str,BOOL year = FALSE);
	void RunAnalysis();
	void AssignIndex(int);
	void DetermineScale();
	void DetermineHistoryScale();
	void FreeAnalysis();
	void FreeAnalysis(ENTITY*);
	BOOL Export(UINT nID,CmyDatabase* = 0);
	BOOL Import(UINT nID,CmyDatabase* = 0);
	void NumberJoint(JOINT*,UINT& N,UINT& SUN,UINT& DUN,UINT& RUN,INTEGER*);
	void NumberDofs(UINT& N,UINT& SUN,UINT& DUN,UINT& RUN,INTEGER*);
	void RestrainK(MATRAN K,UINT N,UINT NB,INTEGER* dof_flags,DOUBLE = 1,DOUBLE = SMALLNUMBER);
	void FormK(MATRAN K,UINT N,UINT NB,BOOL = FALSE);
	BOOL FactorK(MATRAN K,MATRAN FAC,UINT N,UINT NB,INTEGER* dof_flags,INTEGER* IPVT = NULL);
	void DetermineBandwidth(UINT N,UINT& NB);
	void GetExcitationLoad(RESPONSEHIST* responsecase,VECTOR P,VECTOR M,INTEGER* Dir,UINT N,UINT SJ,VECTOR E,DOUBLE step);
	void SolveStaticAnalysisCases(MATRAN K,MATRAN FAC,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index);
	void SolveStaticNonLinearAnalysisCases(UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index);
	void SolveModalAnalysisCases(MATRAN K,MATRAN M,UINT N,UINT DUN,UINT NB,INTEGER* dof_flags,INTEGER& index);
	void SolveBucklingAnalysisCases(MATRAN K,UINT N,UINT DUN,UINT NB,INTEGER* dof_flags,INTEGER& index);
	void SolveTimeHistoryAnalysisCases(MATRAN K,MATRAN M,INTEGER* Dir,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index);
    void SolveResponseSpectrumAnalysisCases(MATRAN K,MATRAN M,MATRAN FAC,INTEGER* Dir,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index);
	void SolveCombination(INTEGER& index);
	void ConvertToEquivalentFrame();
	void UpdateDrawing();
	void RemoveDesignCombos();
	void AddDesignCombos();
	void RestoreObjects();
	void GetGridInfo(GRIDINFO*);

	afx_msg void OnRunAnalysis();
	afx_msg void OnAnalysisOptions();
	afx_msg void OnDesignOptions();
	afx_msg void OnDesignVerify();
	afx_msg void OnMiscOptions();
	afx_msg void OnDesign(UINT nID);
	afx_msg void OnJointLoad(UINT nID);
	afx_msg void OnJointAxis();
	afx_msg void OnJointRestraint();
    afx_msg void OnJointConstraint(UINT nID);
	afx_msg void OnMemberPointLoad();
    afx_msg void OnMemberDistLoad();
	afx_msg void OnMemberRelease();
	afx_msg void OnMemberAxis();
	afx_msg void OnMemberDivisions();
	afx_msg void OnSlabDistLoad();
	afx_msg void OnSlabAxis();
	afx_msg void OnSlabDivisions();
	afx_msg void OnDefineMaterial();
    afx_msg void OnDefineFrameSection();
    afx_msg void OnDefineSlabSection();
	afx_msg void OnDefineSystem();
    afx_msg void OnDefineLoadCase();
	afx_msg void OnDefineModalCase();
	afx_msg void OnDefineResponseHistory();
    afx_msg void OnDefineResponseSpectrum();
    afx_msg void OnDefineLoadCombo();
	afx_msg void OnDefineNonLinearCase();
    afx_msg void OnDefineBucklingCase();
	afx_msg void OnDefineConstraint();
	afx_msg void OnDefineTimeHistoryFunction();
	afx_msg void OnDefineSpectrumFunction();
	afx_msg void OnDefineGroup();
	afx_msg void OnSelectAllJoints();
    afx_msg void OnSelectAllMembers();
	afx_msg void OnSelectAllSlabs();
	afx_msg void OnSelectAll();
    afx_msg void OnDeselectAllJoints();
    afx_msg void OnDeselectAllMembers();
	afx_msg void OnDeselectAllSlabs();
	afx_msg void OnDeselectAll();
	afx_msg void OnSelectLabels(UINT nID);
	afx_msg void OnSelectSection(UINT nID);
	afx_msg void OnSelectSlabSection(UINT nID);
	afx_msg void OnSelectConstraint(UINT nID);
	afx_msg void OnSelectGroup(UINT nID);
	afx_msg void OnSelectFromAll();
	afx_msg void OnUpdateSelectFromAllUI(CCmdUI* pCmdUI);
    afx_msg void OnInvertSelection();
	afx_msg void OnAssignSection();
	afx_msg void OnAssignTaperedSection();
	afx_msg void OnAssignSlabSection();
	afx_msg void OnAssignGroup();
	afx_msg void OnJoinMembers();
    afx_msg void OnBreakMember();
	afx_msg void OnJoinSlabs();
    afx_msg void OnBreakSlab();
	afx_msg void OnUpdateDeselectUI(CCmdUI* pCmdUI);
	afx_msg void OnUpdateAssignJointUI(CCmdUI* pCmdUI);
	afx_msg void OnUpdateAssignMemberUI(CCmdUI* pCmdUI);
	afx_msg void OnUpdateAssignSlabUI(CCmdUI* pCmdUI);
	afx_msg void OnLock();
	afx_msg void OnUpdateLockUI(CCmdUI* pCmdUI);
	afx_msg void OnShowResult();
    afx_msg void OnUpdateAnalysisResultUI(CCmdUI* pCmdUI);
	afx_msg void OnShowDesignResult();
	afx_msg void OnUpdateDesignResultUI(CCmdUI* pCmdUI);
	afx_msg void OnUpdateDesignUI(CCmdUI* pCmdUI);
	afx_msg void OnMemberStrain();
	afx_msg void OnSlabStrain();
	afx_msg void OnStrain(UINT nID);
	afx_msg void OnSelectLoadCase();
	afx_msg void OnPaste();
	afx_msg void OnCopy(UINT nID);
	afx_msg void OnDelete();
	afx_msg void OnEnter();
	afx_msg void OnUpdateLockMenuUI(CCmdUI* pCmdUI);
	afx_msg void OnExport(UINT nID);
	afx_msg void OnImport(UINT nID);
	afx_msg void OnPrepareDrawing(UINT nID);
	afx_msg void OnPlotFunctions();
	afx_msg void OnUpdatePlotFunctionsUI(CCmdUI* pCmdUI);
	afx_msg void OnSnapTo(UINT nID);
	afx_msg void OnUpdateSnapToUI(CCmdUI* pCmdUI);
	afx_msg void OnShowTables();
	afx_msg void OnUndoRedo(UINT nID);
	afx_msg void OnUpdateUndoRedoUI(CCmdUI* pCmdUI);
	afx_msg void OnAddDesignCombos();
	afx_msg void OnReverseSlabAxis();
	afx_msg void OnColors();
	DECLARE_MESSAGE_MAP()
};
/*
View
*/

class CmyView : public CView {
	DECLARE_DYNCREATE(CmyView)
public:
	VIEW view;
	VIEW perspective_toggle;
	DOUBLE  pvalue;
	int position;
	CPoint  Center;
	CPoint  Corner;
	DOUBLE  Scale;
	RPoint  eye; 
	RPoint  reye;
	RPoint  gcenter;
	static UINT nViews;
    	
	VIEW oldview;
	DOUBLE oldpvalue;
	BOOL hasGridPoint;
	BOOL hadGridPoint;
	RPoint oldGridPoint;
	
	BOOL   has_tracker;
	BOOL   had_tracker;
	CPoint prev_mouse_cp;
	BOOL   has_plane;
	BOOL   had_plane;
	DOUBLE SELSIZE;
	DISPLAY mydisplay;
	int force_diagram;
	BOOL average_stress;
	BOOL DrawingHasChanged;

	static VIEW CREATION[4];
	
	CDC memDC;
	CBitmap memBmp;
	int memX;
	int memY;

	DOUBLE min_stress;
	DOUBLE max_stress;


	TOOLINFO m_ti;
	CToolTipCtrl m_toolTip;
	CAnimationDia* pAnimationDia;

	enum {SNAP_GRID = 1,SNAP_JOINT,SNAP_MEMBERMID,SNAP_MEMBER};
	enum {LINE_TRACK,RECT_TRACK,MEMBER_TRACK};

	static COLORREF BackColor;
	static COLORREF JointColor;
	static COLORREF MemberColor;
	static COLORREF SlabColor;
	static COLORREF HighlightColor;
	static COLORREF SelectionColor;
	static COLORREF RestraintColor;
	static COLORREF ConstraintColor;
	static COLORREF AxisColor;
	static COLORREF ReleaseColor;
	static COLORREF TextColor;
	static COLORREF GridColor;
	static COLORREF PositiveColor;
	static COLORREF NegativeColor;
	static COLORREF ContorColors[];
	static UINT NCONTORS;
public:
	CmyView();
	~CmyView();
	CmyDocument* GetDocument() {
		return (CmyDocument*) m_pDocument;
	}
	static void ResetColors();
	void Serialize(CSerializer&);
	static void SerializeColors(CSerializer&);
	BOOL InPlane(RPoint);
	RPoint RPtoRP(RPoint,SYSTEM*,BOOL = FALSE);
	CPoint RPtoCP(RPoint,SYSTEM* = 0,BOOL = FALSE); 
	RPoint CPtoRP(CPoint);
	RPoint SnapToGrid(CPoint&,int&,JOINT** = NULL);
	RPoint SnapToGridLocal(SYSTEM*,CPoint&,BOOL&);
	void InitView();
	void OnDraw(CDC*);
	BOOL OnPreparePrinting(CPrintInfo* pInfo);
	void OnPrint(CDC* pDC,CPrintInfo* pInfo);
	void AdjustCenter();
	void AdjustScale();
	void AdjustPosition();
	void InitMinMaxStress();
	DOUBLE FindStress(SLAB* slab,RPoint);
	void GetEdgeStress(SLAB* slab,DOUBLE* V);
	void DrawCircle(CDC*,RPoint&,int,COLORREF);
	void DrawGridPoint(RPoint&,BOOL);
	void DrawJoint(JOINT&,BOOL = FALSE,CDC* pDC = NULL);
	void DrawMember(MEMBER&,CDC* pDC = NULL);
	void DrawColorTable(CDC* pDC);
	void DrawSlab(SLAB&,CDC* pDC = NULL);
	void DrawSlabMiniBlockStress(CDC*,DOUBLE*,RPoint*);
	void DrawSlabBlockStress(CDC*,SLAB&,RPoint*);
	void DrawSlabStress(CDC*,SLAB&,MATRIX);
	void DrawDeflectedShape(CDC*,MEMBER&);
	void DrawTracker(CPoint&,int);
	void DrawPlaneGrid(CDC*,SYSTEM*,DOUBLE*,int,int);
	void DrawGridOnView(CDC*,SYSTEM*,int,DOUBLE);
	void DrawGrids(CDC*);
	void DrawView(CDC*);
	void DrawPlane(CmyView*);
	void DrawOrigin(CDC*);
	void DrawAxis(CDC*,RPoint&,MATRIX,DOUBLE);
    void DrawRestraint(CDC*,JOINT*,MATRIX);
	void DrawArrow(CDC*,RPoint&,RPoint&,int = NOARROW);
	void DrawText(CDC*,CString,CPoint*,CPoint* = NULL);
	void DrawRectangle(CDC*,RPoint&,RPoint&,RPoint&,RPoint&);
	void DrawTriangle(CDC*,RPoint&,RPoint&,RPoint&);
	void UpdatePlanes(BOOL);
	BOOL IsNear(RPoint,RPoint,RPoint);
	JOINT* FindJoint(RPoint,BOOL OriginalPos = FALSE);
	MEMBER* FindMember(RPoint,RPoint* = 0,DOUBLE* = 0);
	SLAB* FindSlab(RPoint);
	BOOL SelectJoint(RPoint);
	BOOL SelectMember(RPoint);
    BOOL SelectSlab(RPoint);
	void Animate();
	RPoint GetJointPoint(JOINT* joint,CmyDocument* pDoc);
	void SetTitle();
	void CalculateReye();
public:
	BOOL PreCreateWindow(CREATESTRUCT&);
	afx_msg  int OnCreate(LPCREATESTRUCT);
	afx_msg void OnSize(UINT nType,int cx,int cy);
	afx_msg void OnLButtonDown(UINT,CPoint);
	afx_msg void OnRButtonDown(UINT,CPoint);
	afx_msg void OnLButtonUp(UINT,CPoint);
	afx_msg void OnMouseMove(UINT,CPoint);
	afx_msg BOOL OnMouseWheel(UINT,short,CPoint);
	afx_msg LRESULT OnMouseLeave(WPARAM wPawam, LPARAM lParam);
    afx_msg void OnChangeView(UINT nID);
	afx_msg void OnCoordinate(CCmdUI *pCmdUI);
	afx_msg BOOL OnSetCursor(CWnd* pWnd,UINT nHitTest,UINT message);
	afx_msg void OnSetFocus(CWnd* pOldWnd);
	afx_msg void OnKillFocus(CWnd* pOldWnd);
	afx_msg void OnEscape();
	afx_msg void OnAdd(UINT nID);
	afx_msg void OnUpdateAddUI(CCmdUI* pCmdUI);
	afx_msg void OnMoveView(UINT nID);
	afx_msg void OnUpdateMoveViewUI(CCmdUI* pCmdUI);
	afx_msg void OnZoomIn();
    afx_msg void OnZoomOut();
	afx_msg void OnRestoreView();
    afx_msg void OnRubberZoom();
    afx_msg void OnPan();
	afx_msg void OnDisplayOptions();
	afx_msg void OnDiagrams();
	afx_msg void OnSlabDiagrams();
	afx_msg void OnUpdateDiagramsUI(CCmdUI* pCmdUI);
	afx_msg void OnDeflected();
	afx_msg void OnReactions();
	afx_msg void OnDisplacement();
	afx_msg void OnRotate(UINT nID);
	afx_msg void OnUpdateRotateUI(CCmdUI* pCmdUI);
	afx_msg void OnShowLoads();
    afx_msg void OnShowGrids();
	afx_msg void OnShowOrigin();
	afx_msg void OnRefresh();
	afx_msg void OnAnimate();
    afx_msg void OnUpdateAnimateUI(CCmdUI* pCmdUI);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnPerspective();
	DECLARE_MESSAGE_MAP()
};
/*
Child Frame
*/
class CChildFrame:public CMDIChildWnd  {
	DECLARE_DYNCREATE(CChildFrame)
public:
	CChildFrame() {};
	~CChildFrame() {};
	afx_msg void OnClose();
	BOOL PreCreateWindow(CREATESTRUCT&);
	DECLARE_MESSAGE_MAP()
};
/*
Custom tool bar
*/
class CCommonToolBar : public CToolBar {
};

class CMyToolBar : public CCommonToolBar {
public:
	CComboBox c_List;
public:
    afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
    DECLARE_MESSAGE_MAP()
};
/*
Main Frame
*/
class CMainFrame : public CMDIFrameWnd {
	DECLARE_DYNCREATE(CMainFrame)
public:
	CStatusBar myStatusBar;
	CMyToolBar myToolBar;
	CCommonToolBar myAnalysisToolBar;
	CCommonToolBar myDesignToolBar;
	CCommonToolBar myDrawToolBar;
	CCommonToolBar mySelectToolBar;
	CCommonToolBar myViewToolBar;
	CCommonToolBar mySnapToolBar;
	CCommonToolBar myEditToolBar;
	CCommonToolBar myAssignToolBar;
public:
	CMainFrame() {};
	~CMainFrame() {};
	CToolBar* GetTool(UINT nID);
	afx_msg void OnToolbar(UINT nID);
    afx_msg void OnUpdateToolbarUI(CCmdUI* pCmdUI);
	afx_msg void ResetToolBars();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnUpdateLoadcaseUI(CCmdUI* pCmdUI);
	afx_msg void OnSize(UINT nType,int cx,int cy);
	void DockControlBarBottomOf(CToolBar* Bar,CToolBar* LeftOf);
	DECLARE_MESSAGE_MAP()
};
/*
Application
*/
class CmyApplication : public CWinApp {
public:
	BOOL skip_model;
	char app_path[256];
public:
	CmyApplication() {
		skip_model = FALSE;
	};
	virtual BOOL InitInstance();
	
	void CreateMyView();
	void DestroyMyView();
	void OpenViews(UINT count);
	afx_msg void OnAppAbout();
	afx_msg void OnFileNew();
	afx_msg void OnFileOpen();
	afx_msg void OnWindows(UINT nID);
	afx_msg void OnImport(UINT nID);
protected:
	DECLARE_MESSAGE_MAP()
};
/*
Get Document
*/
__inline CmyDocument* GetMyDocument() {
	return CmyDocument::Current_Document;
}
#endif