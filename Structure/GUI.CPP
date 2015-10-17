#include "common.h"
#include "GUI.h"

MATERIAL* SECTION::defmaterial;
MATERIAL* ASECTION::defmaterial;
SECTION* MEMBER::defsection;
ASECTION* SLAB::defsection;
MODALCASE* RESPONSE::defmodalcase;
FUNCTION* RESPONSE::deffunction;
BOOL FUNCTION::defspectrumtype;
VIEW CmyView::CREATION[4] = {
	VIEW_3D,VIEW_XY,VIEW_XZ,VIEW_YZ
};

static CString DirString[]   = {"UX","UY","UZ","RX","RY","RZ"};
static CString ForceString[] = {"FX","FY","FZ","MX","MY","MZ"};
static CString CoordinateString[] = {"GLOBAL","LOCAL"};
static CString StressString[] = {
	"SX","SY","SXY","SMAX","SMIN","SVM","SYZ","SXZ","SVMAX",
	"FX","FY","FXY","FMAX","FMIN","FVM",
    "MX","MY","MXY","MMAX","MMIN",
    "VYZ","VXZ","VMAX"
};
static CString LoadTypeStr[] = {
	"CONCENTRATED","UNIFORM","TRAPEZOIDAL","FORCE","DISPLACEMENT",
	"MASS","TEMPERATURE","STRAIN"
};
/*
Function
*/
DOUBLE FUNCTION::Get(DOUBLE x,int type,DOUBLE scale,BOOL extrapolate) {
	RPoint v,p1,p2;
	DOUBLE xmin = points.GetHead().x,
		   ymin = points.GetHead().y,
		   xmax = points.GetTail().x,
		   ymax = points.GetTail().y;

	if(type != PERIODIC) {
		if(extrapolate) {
			POSITION pos;
			if(points.GetCount() < 2) 
				return (scale * ymax);
            if(x < xmin) {
				p1 = points.GetHead();
                pos = points.GetHeadPosition();
				points.GetNext(pos);
				p2 = points.GetNext(pos);
				return scale * (p1.y + (p2.y - p1.y) * (x - p1.x) / (p2.x - p1.x));
			} else if(x > xmin + cycles * (xmax - xmin)) {
				p1 = points.GetTail();
                pos = points.GetTailPosition();
				points.GetPrev(pos);
				p2 = points.GetPrev(pos);
				return scale * (p1.y + (p2.y - p1.y) * (x - p1.x) / (p2.x - p1.x));
			}
		} else {
		    if(x < xmin) return 0;
			else if(x > xmin + cycles * (xmax - xmin)) return 0;
		}
	}
	
	POSITION pos = points.GetHeadPosition(),ppos = pos,ppos1 = ppos;
	while(pos) {
		ppos1 = ppos;
		ppos = pos;
		v = points.GetNext(pos);
		if(fabs(v.x - x) <= 0.001)
			return scale * v.y;
		if(v.x > x)
			break;
	}
	p1 = points.GetAt(ppos1);
	p2 = points.GetAt(ppos);
	return scale * (p1.y + (p2.y - p1.y) * (x - p1.x) / (p2.x - p1.x));
}
/*
Numeric Edit
*/
CNumericEdit::CNumericEdit() {
	deffmt = NULL;
}
CNumericEdit::CNumericEdit(LPCTSTR fmt) {
	deffmt = fmt;
}
void CNumericEdit::SetWindowText(DOUBLE val, LPCTSTR fmt) {
	if(fmt == NULL) {
		if(deffmt == NULL) fmt = _T("%.6e");
		else fmt = deffmt;
	}
    CString str;
	str.Format(fmt,val);

	CString oldstr;
	CEdit::GetWindowText(oldstr);
	if(lstrcmp(oldstr,str) != 0)
		CEdit::SetWindowText(str);

}
void CNumericEdit::SetWindowText(LPCTSTR val) {
	CEdit::SetWindowText(val);
}
/*
Grid control
*/
BEGIN_MESSAGE_MAP(CGridCtrl, CListCtrl)
    ON_NOTIFY_REFLECT(NM_CUSTOMDRAW, OnCustomDraw)
	ON_NOTIFY_REFLECT_EX(NM_CLICK,OnClick)
END_MESSAGE_MAP()

CGridCtrl::~CGridCtrl() {
	if(type != TYPE1) {
		GRIDCTRLDATA* pXLCD;
		UINT count = GetItemCount();
		for(UINT i = 0;i < count;i++) {
			pXLCD = (GRIDCTRLDATA*) CListCtrl::GetItemData(i);
			if(!pXLCD)
				continue;
			delete pXLCD;
		}
	}
}
void CGridCtrl::PreSubclassWindow() {
	CListCtrl::PreSubclassWindow();

	DWORD style;
	style = GetStyle();
	style |= LVS_REPORT | LVS_SINGLESEL;
	if(type == TYPE3)
		style |= LVS_NOCOLUMNHEADER | LVS_NOSCROLL;
	SetWindowLong(m_hWnd,GWL_STYLE,style);

	style = GetExtendedStyle();
	style |= LVS_EX_GRIDLINES | LVS_EX_FULLROWSELECT;
	SetExtendedStyle(style);
}
void CGridCtrl::SetTwoHeaders() {
	CRect rect;
	GetClientRect(&rect);
	UINT nInterval = (rect.Width()) / 2;
	InsertColumn(0,"Item",LVCFMT_CENTER,nInterval);
	InsertColumn(1,"Value",LVCFMT_CENTER,nInterval);
}
int CGridCtrl::InsertItem(const LVITEM* pItem) {
	if(type == TYPE1) return CListCtrl::InsertItem(pItem);
	else {
		GRIDCTRLDATA data;
		return InsertItem(pItem,&data);
	}
}
int CGridCtrl::InsertItem(const LVITEM* pItem,GRIDCTRLDATA* pData) {
	if (pItem->iItem < 0)
		return -1;

	int index = CListCtrl::InsertItem(pItem);

	if (index < 0)
		return index;

	GRIDCTRLDATA *pXLCD = new GRIDCTRLDATA;
	if (!pXLCD)
		return -1;

    *pXLCD = *pData;
	CListCtrl::SetItemData(index, (DWORD) pXLCD);

	CString str;
	switch(pXLCD->data_type) {
	case GRIDCTRLDATA::TINT: str.Format("%d",*(int*)(pXLCD->pdata)); break;
	case GRIDCTRLDATA::TDBL: str.Format("%.2f",*(double*)(pXLCD->pdata)); break;
	case GRIDCTRLDATA::TDBL_F: str.Format("%.6g",*(double*)(pXLCD->pdata)); break;
	case GRIDCTRLDATA::TTEXT: str = *(CString*)(pXLCD->pdata); break;
	case GRIDCTRLDATA::TBOOL: str = *(BOOL*)(pXLCD->pdata) ? "Yes" : "No"; break;
	case GRIDCTRLDATA::TINDEX: str = pXLCD->list.GetAt(pXLCD->list.FindIndex(*(int*)pXLCD->pdata)); break;
	}
	SetItemText(index,1,str);

	return index;
}
void CGridCtrl::InsertMyItem(int type,int data_type,void* pdata,CString ItemText,CString* list) {
	LV_ITEM lvItem;
	lvItem.mask = LVIF_TEXT | LVIF_STATE;
	lvItem.state = 0;
	lvItem.stateMask = 0;
	lvItem.iSubItem = 0;

	GRIDCTRLDATA data;
	data.type = type;
	data.data_type = data_type;
	data.pdata = pdata;

	/*add string*/
	CString str;
	if(type == GRIDCTRLDATA::COMBO) {
		if(list) {
			UINT i = 0;
			while(1) {
				str = list[i++];
				if(str == "0") break;
				data.list.AddTail(str);
			};
		} else if(data_type == GRIDCTRLDATA::TBOOL) {
			str = "Yes";
			data.list.AddTail(str);
			str = "No";
			data.list.AddTail(str);
		}
	}

	lvItem.iItem = MyItemCount;
	lvItem.pszText = ItemText.GetBuffer(0);
	InsertItem(&lvItem,&data);
	MyItemCount++;
}
void CGridCtrl::OnCustomDraw(NMHDR* pNMHDR, LRESULT* pResult) {
	NMLVCUSTOMDRAW* pLVCD = reinterpret_cast<NMLVCUSTOMDRAW*>(pNMHDR);
	NMCUSTOMDRAW* pCD = reinterpret_cast<NMCUSTOMDRAW*>(pNMHDR);
	*pResult = CDRF_DODEFAULT;
	if (pLVCD->nmcd.dwDrawStage == CDDS_PREPAINT) {
		*pResult = CDRF_NOTIFYITEMDRAW;
	} else if (pLVCD->nmcd.dwDrawStage == CDDS_ITEMPREPAINT) {
		*pResult = CDRF_NOTIFYSUBITEMDRAW;
	} else if (pLVCD->nmcd.dwDrawStage == (CDDS_ITEMPREPAINT | CDDS_SUBITEM)) {
		int nItem = static_cast<int> (pLVCD->nmcd.dwItemSpec);
		int nSubItem = pLVCD->iSubItem;
		
		if(nSubItem == 0) 
			return;

		GRIDCTRLDATA *pXLCD = (GRIDCTRLDATA *) CListCtrl::GetItemData(nItem);
		CDC* pDC = CDC::FromHandle(pLVCD->nmcd.hdc);
		CRect rect;
		GetSubItemRect(nItem, nSubItem, LVIR_BOUNDS, rect);

		if(pXLCD && pXLCD->type == GRIDCTRLDATA::COLOR) {
			CString str = GetItemText(nItem, nSubItem);
			COLORREF col = atoi(str);
			pDC->FillSolidRect(&rect, col);
		} else {
			COLORREF crText  = RGB(255,0,0);
			COLORREF crBkgnd = RGB(100,255,255);
			
			if(pCD->uItemState == CDIS_SELECTED)
				pDC->FillSolidRect(&rect, RGB(200,200,200));
			else
				pDC->FillSolidRect(&rect, crBkgnd);
			
			CString str;
			str = GetItemText(nItem, nSubItem);
			if(!str.IsEmpty()) {
				pDC->SetBkMode(TRANSPARENT);
				pDC->SetTextColor(crText);
				pDC->SetBkColor(crBkgnd);
				pDC->DrawText(str, &rect, DT_VCENTER | DT_SINGLELINE | DT_CENTER);
			}
		}

		if(pCD->uItemState != CDIS_SELECTED) 
			*pResult = CDRF_SKIPDEFAULT;
		else
			*pResult = CDRF_SKIPDEFAULT;
	}
}
BOOL CGridCtrl::OnClick(LPNMHDR pNMHDR, LRESULT* /*pResult*/) {
	if(type == TYPE1)
		return FALSE;

	NMLISTVIEW* pnmlv = (NMLISTVIEW*)pNMHDR;
	int nItem = pnmlv->iItem;
	int nSubItem = pnmlv->iSubItem;
	GRIDCTRLDATA *pXLCD = (GRIDCTRLDATA *) CListCtrl::GetItemData(nItem);

	if(nSubItem == 0)
		nSubItem++;
	if(pControl) {
		CString str;
        pControl->GetWindowText(str);
        SetItemText(selItem,selSubItem,str);
		delete pControl;
		pControl = NULL;
	}
	if(nItem < 0)
		return FALSE;
	if(pXLCD->type == GRIDCTRLDATA::NONE)
		return FALSE;

	COLORREF crBkgnd = RGB(100,255,255);
	CString str,str1;
	CRect rect;
	GetSubItemRect(nItem, nSubItem, LVIR_BOUNDS, rect);
	str1 = GetItemText(nItem,nSubItem);

	if(pXLCD->type == GRIDCTRLDATA::COMBO) {
		CComboBox* cBox = new CComboBox();
		rect.bottom = 300;
		cBox->Create(WS_CHILD|WS_VISIBLE|WS_VSCROLL|CBS_DROPDOWNLIST,rect, this, 1);
	
		POSITION pos = pXLCD->list.GetHeadPosition();
		while(pos) {
			str = pXLCD->list.GetNext(pos);
			cBox->AddString(str);
		}

		cBox->SetCurSel(cBox->FindString(-1,str1));
		cBox->SetFocus();
		pControl = cBox;
		selItem = nItem;
		selSubItem = nSubItem;
	} else if(pXLCD->type == GRIDCTRLDATA::EDIT) {
		CEdit* cBox = new CEdit();
		cBox->Create(WS_CHILD|WS_VISIBLE|ES_CENTER,rect, this, 1);
		cBox->SetWindowText(str1);
		cBox->SetFocus();
		pControl = cBox;
		selItem = nItem;
		selSubItem = nSubItem;
	} else if(pXLCD->type == GRIDCTRLDATA::COLOR) {
		COLORREF color = *((COLORREF*)(pXLCD->pdata));
		CColorDialog Dia(color,CC_FULLOPEN);
		if(Dia.DoModal() == IDOK) {
			color = Dia.GetColor();
			CString str;
			str.Format("%d",color);
			SetItemText(nItem,1,str);
		}
	}
	return TRUE;
}
void CGridCtrl::GetData() {
	GRIDCTRLDATA *pXLCD;
	CString str,str1;
	POSITION pos;
	UINT count = GetItemCount(),index;

	if(pControl) {
		CString str;
        pControl->GetWindowText(str);
        SetItemText(selItem,selSubItem,str);
		delete pControl;
		pControl = NULL;
	}

	for(UINT i = 0;i < count;i++) {
		pXLCD = (GRIDCTRLDATA *) CListCtrl::GetItemData(i);
		if(!pXLCD->pdata)
			continue;
		str = GetItemText(i,1);
		switch(pXLCD->data_type) {
		case GRIDCTRLDATA::TINT: *(int*)(pXLCD->pdata) = atoi(str); break;
		case GRIDCTRLDATA::TDBL: *(double*)(pXLCD->pdata) = atof(str); break;
		case GRIDCTRLDATA::TTEXT: *(CString*)(pXLCD->pdata) = str; break;
		case GRIDCTRLDATA::TBOOL: *(BOOL*)(pXLCD->pdata) = (str == "Yes"); break;
		case GRIDCTRLDATA::TINDEX: 
			index = 0;
            pos = pXLCD->list.GetHeadPosition();
			while(pos) {
				str1 = pXLCD->list.GetNext(pos);
				if(str == str1) {
					*(int*)(pXLCD->pdata) = index;
					break;
				}
				index++;
			}
			break;
		}
	}
}
/*
Database modifiers
*/
int MATERIAL::Modify(void* other,void* other1) {
	CMaterialPropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}
int SECTION::Modify(void* other,void* other1) {
	CSectionPropertyDia PropertyDia(AfxGetMainWnd(),this,(MATERIALLIST*)other);
	return PropertyDia.DoModal();
}
int ASECTION::Modify(void* other,void* other1) {
	CASectionPropertyDia PropertyDia(AfxGetMainWnd(),this,(MATERIALLIST*)other);
	return PropertyDia.DoModal();
}
int GROUP::Modify(void* other,void* other1) {
	CGroupPropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}
int SYSTEM::Modify(void* other,void* other1) {
	CSystemPropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}
int LOADCASE::Modify(void* other,void* other1) {
	CLoadCasePropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}
int MODALCASE::Modify(void* other,void* other1) {
	CModalCasePropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}
int RESPONSEHIST::Modify(void* other,void* other1) {
	CResponsePropertyDia PropertyDia(IDD_RESPONSECASE,AfxGetMainWnd(),this,(MODALCASELIST*)other,(FUNCTIONLIST*)other1);
	return PropertyDia.DoModal();
}
int RESPONSESPEC::Modify(void* other,void* other1) {
	CResponsePropertyDia PropertyDia(IDD_RESPONSESPEC,AfxGetMainWnd(),this,(MODALCASELIST*)other,(FUNCTIONLIST*)other1);
	return PropertyDia.DoModal();
}
int COMBINATION::Modify(void* other,void* other1) {
	CCombinationPropertyDia PropertyDia(AfxGetMainWnd(),this,(PANALYSISCASELIST*)other);
	return PropertyDia.DoModal();
}
int NLCASE::Modify(void* other,void* other1) {
	CNLcasePropertyDia PropertyDia(AfxGetMainWnd(),this,(PANALYSISCASELIST*)other);
	return PropertyDia.DoModal();
}
int BUCKLINGCASE::Modify(void* other,void* other1) {
	CBucklingPropertyDia PropertyDia(AfxGetMainWnd(),this,(PANALYSISCASELIST*)other);
	return PropertyDia.DoModal();
}
int CONSTRAINT::Modify(void* other,void* other1) {
	CConstraintPropertyDia PropertyDia(AfxGetMainWnd(),this,(PSYSTEMLIST)other);
	return PropertyDia.DoModal();
}
int FUNCTION::Modify(void* other,void* other1) {
	CFunctionPropertyDia PropertyDia(AfxGetMainWnd(),this);
	return PropertyDia.DoModal();
}

/*
Template Define Database class
*/
BEGIN_TEMPLATE_MESSAGE_MAP(CDefineDia, T , CDialog)
    ON_LBN_SELCHANGE(IDC_DATA_LIST,OnUpdateControls)
    ON_COMMAND(IDC_DATA_ADD,OnAdd)
    ON_COMMAND(IDC_DATA_DELETE,OnDelete)
	ON_COMMAND(IDC_DATA_MODIFY,OnModify)
	ON_COMMAND(IDC_DATA_IMPORT,OnImport)
	ON_COMMAND(IDC_DATA_EXPORT,OnExport)
END_MESSAGE_MAP()

template <class T>
CDefineDia<T>::CDefineDia(CWnd* parent,CList<T,T&>* plist,CString text,void* pother,void* pother1,BOOL timport) : 
                    CDialog(IDD_DATABASE,parent) {
	list = plist;
	title = text;
	other = pother;
	other1 = pother1;
	import = timport;
}
template <class T>
BOOL CDefineDia<T>::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	SetWindowText(title);
	POSITION pos = list->GetHeadPosition();
	while(pos) {
		c_List.AddString(list->GetNext(pos).name);
	}
	c_List.SetCurSel(0);
	OnUpdateControls();
	
	return TRUE;
}
template <class T>
void CDefineDia<T>::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
    DDX_Control(pDX,IDC_DATA_LIST,c_List);
	DDX_Control(pDX,IDC_DATA_DELETE,c_Delete);
	DDX_Control(pDX,IDC_DATA_MODIFY,c_Modify);
	DDX_Control(pDX,IDC_DATA_IMPORT,c_Import);
	DDX_Control(pDX,IDC_DATA_EXPORT,c_Export);
}
template <class T>
void CDefineDia<T>::OnAdd() {
	T t;
	int index;

	CString num;
	num.Format("%d",list->GetCount() + 1);
	t.name += num;
	list->AddTail(t);
	index = c_List.SetCurSel(c_List.AddString(t.name));

	OnModify();
	if(!modify) {
		list->RemoveTail();
		c_List.DeleteString(index);
		c_List.SetCurSel(index - 1); 
		OnUpdateControls();
	}
}
template <class T>
void CDefineDia<T>::OnDelete() {
	T t;
	POSITION pos;
	int index = c_List.GetCurSel();
	c_List.GetText(index,t.name);
	if(pos = list->Find(t)) {
		list->RemoveAt(pos);
		c_List.DeleteString(index);
		if(index) c_List.SetCurSel(index - 1);
		else c_List.SetCurSel(index);
		OnUpdateControls();
	}
}
template <class T>
void CDefineDia<T>::OnModify() {
	T t;
	POSITION pos;
    int index = c_List.GetCurSel();
	c_List.GetText(index,t.name);

	if(pos = list->Find(t)) {
		T* pt = &list->GetAt(pos);
		t = *pt;
		
		if(t.Modify(other,other1) == IDOK) {
			*pt = t;
			c_List.DeleteString(index);
			c_List.SetCurSel(c_List.AddString(t.name));
			OnUpdateControls();
			modify = TRUE;
			return;
		}
		
	}
	modify = FALSE;
}
template <class T>
void CDefineDia<T>::OnOK() {
	T t;
	POSITION pos;
	c_List.GetText(c_List.GetCurSel(),t.name);
	if(pos = list->Find(t)) {
        pT = &list->GetAt(pos);
	}
	Clear();
	CDialog::OnOK();
}
template <class T>
void CDefineDia<T>::OnExport() {
	char strFilter[] = { "Database files (*.dat)|*.dat|All Files (*.*)|*.*||" };
	CFileDialog FileDlg(FALSE,".dat",NULL,0,strFilter);
	
	if(FileDlg.DoModal() == IDOK) {
		CString name = FileDlg.GetPathName();
		CFile viewfile;
		if(viewfile.Open(name, CFile::modeCreate | CFile::modeWrite)) {
			CArchive ar(&viewfile, CArchive::store);
			CSerializer sr;
			sr.Set(ar,TRUE);
			
			POSITION pos;
			UINT i,count;
			
			count = list->GetCount();                     
			ar << count;                                
			pos = list->GetHeadPosition();                
			for(i = 0; i < count;i++) {                    
				T* mydat;                         
				mydat = &(list->GetNext(pos));            
				mydat->Serialize(sr);                      
			}     
			

			ar.Close();
			viewfile.Close();
		}
	}
}
template <class T>
void CDefineDia<T>::OnImport() {
	char strFilter[] = { "Database files (*.dat)|*.dat|All Files (*.*)|*.*||" };
	CFileDialog FileDlg(TRUE,".dat",NULL,0,strFilter);
	
	if(FileDlg.DoModal() == IDOK) {
		CString name = FileDlg.GetPathName();
		CFile viewfile;
		if(viewfile.Open(name, CFile::modeRead)) {
			CArchive ar(&viewfile, CArchive::load);
			CSerializer sr;
			sr.Set(ar,FALSE);

			POSITION pos;
			UINT i,count;
			
			ar >> count;
			for(i = 0;i < count;i++) {                     
				T mydat;                          
				mydat.Serialize(sr);  
				pos = list->Find(mydat);
				if(!pos) {
					list->AddTail(mydat);  
					c_List.AddString(mydat.name);
				} else {
					list->SetAt(pos,mydat);  
				}
			}   
			
			ar.Close();
			viewfile.Close();
		}
	}
}
template <class T>
void CDefineDia<T>::OnCancel() {
	Clear();
	CDialog::OnCancel();
}
template <class T>
void CDefineDia<T>::Clear() {
	POSITION pos;
	T* t;
	pos = list->GetHeadPosition();
	while(pos) {
		t = &list->GetNext(pos);
		if(t->rank == LEVEL2)
            t->rank = LEVEL1;
	}
}
template <class T>
void CDefineDia<T>::OnUpdateControls() {
	T t;
	POSITION pos;
	T* pt;
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		c_List.GetText(index,t.name);
		if(pos = list->Find(t)) {
			pt = &list->GetAt(pos);
			c_Delete.EnableWindow(pt->rank == LEVEL1 || pt->rank == LEVEL4);
		}
		c_Modify.EnableWindow(TRUE);
	} else {
		c_Delete.EnableWindow(FALSE);
		c_Modify.EnableWindow(FALSE);
	}
	c_Import.ShowWindow(import);
	c_Export.ShowWindow(import);
}
/*
Template Define Database class
*/
template <class T>
CSelectDia<T>::CSelectDia(CWnd* parent,CList<T,T&>* plist,CString text) : 
                    CDialog(IDD_SELECT_DATABASE,parent) {
	list = plist;
	title = text;
}
template <class T>
BOOL CSelectDia<T>::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	SetWindowText(title);
	POSITION pos = list->GetHeadPosition();
	while(pos) {
		c_List.AddString(list->GetNext(pos).name);
	}
	c_List.SetCurSel(0);
	
	return TRUE;
}
template <class T>
void CSelectDia<T>::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
    DDX_Control(pDX,IDC_DATA_LIST,c_List);
}
template <class T>
void CSelectDia<T>::OnOK() {
	int nCount = c_List.GetSelCount();
	CArray<int,int> sel_index;
	sel_index.SetSize(nCount);
	c_List.GetSelItems(nCount,sel_index.GetData());

	selected_list.RemoveAll();
	CString name;
	for(int i = 0;i < nCount;i++) {
		c_List.GetText(sel_index[i],name);
		selected_list.AddTail(name);
	}

	CDialog::OnOK();
}
/*
Curved Member Dialog
*/
CCurvedMemberDia::CCurvedMemberDia(CWnd* parent) : 
               CDialog(IDD_CURVED_MEMBER,parent) {
	nDiv = 10;
	type = CIRCULAR;
	single = TRUE;
}
void CCurvedMemberDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Text(pDX,IDC_DIVISIONS,nDiv);
	DDV_MinMaxInt(pDX,nDiv,1,10000);
	DDX_Text(pDX,IDC_IX,rp.x);
    DDX_Text(pDX,IDC_IY,rp.y);
	DDX_Text(pDX,IDC_IZ,rp.z);
	DDX_CBIndex(pDX,IDC_TYPE,type);
	DDX_Check(pDX,IDC_SINGLE,single);
}
BOOL CCurvedMemberDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
    c_Type.AddString("Straight");
	c_Type.AddString("Circular");
	c_Type.AddString("Parabolic");
	c_Type.SetCurSel(type);
	return TRUE;
}
/*
New model Dialog
*/
CNewDia::CNewDia(CWnd* parent) : 
               CDialog(IDD_NEW,parent) {
	type = MODELGRID;
}
void CNewDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_BLANK, type);
}
/*
New truss Dialog
*/
CNewBeamDia::CNewBeamDia(CWnd* parent) : 
               CDialog(IDD_NEW_BEAM,parent) {
	
	nspan = 2;
	spanlength = 6;
	restraint = TRUE;
}
void CNewBeamDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NSPAN, nspan);
	DDV_MinMaxInt(pDX,nspan,0,1000);
	DDX_Text(pDX, IDC_SPANLENGTH, spanlength);
	DDV_MinMaxDouble(pDX,spanlength,0,1000000);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
}

void CNewBeamDia::OnOK() {
	CmyDocument* pDoc = GetMyDocument();
	CDialog::OnOK();
	register int i,j;
	DOUBLE v;
	for(j = 0;j <= nspan;j++) {
		v = j * spanlength - nspan * spanlength / 2;
		pDoc->global->grid[IUX].AddTail(v);
	}
	v = 0;
	pDoc->global->grid[IUY].AddTail(v);
	pDoc->global->grid[IUZ].AddTail(v);
	pDoc->global->coordinate = CARTESIAN;

	/*members and restraints*/
	for(i = 0;i < nspan;i++) {
		JOINT j1,j2;
		j1.p.x = i * spanlength - nspan * spanlength / 2;
		j2.p.x = j1.p.x + spanlength;
		pDoc->AddMember(&j1,&j2);
	}
	if(restraint) {
		JOINT* joint;
		POSITION pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			joint->restraint = UY | UZ | RX | RZ;
		}
		pDoc->joints.GetHead().restraint |= UX;
	}
	CmyView::CREATION[1] = VIEW_XZ;
	CmyView::CREATION[2] = VIEW_XY;
}
/*
New truss Dialog
*/
BEGIN_MESSAGE_MAP(CNewTruss2DDia , CDialog)
    ON_CBN_SELCHANGE(IDC_TYPE,OnType)
END_MESSAGE_MAP()

CNewTruss2DDia::CNewTruss2DDia(CWnd* parent) : 
               CDialog(IDD_NEW_2DTRUSS,parent) {
	
	/*pratt*/
	L1 = 5;
	L2 = 9;
	N1 = 4;
	N2 = 6;
	H1 = 1;
	H2 = 4;
	H3 = 2;
	/*other*/
	ndiv = 3;
	divlength = 3;
	height = 3;
    type = SLOPED;
	restraint = TRUE;
}
void CNewTruss2DDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NDIV, ndiv);
	DDV_MinMaxInt(pDX,ndiv,0,1000);
	DDX_Text(pDX, IDC_DIVLENGTH, divlength);
	DDV_MinMaxDouble(pDX,divlength,0,1000000);
	DDX_Text(pDX, IDC_HEIGHT, height);
	DDV_MinMaxDouble(pDX,height,0,1000000);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
	DDX_Text(pDX,IDC_L1,L1);
	DDV_MinMaxDouble(pDX,L1,0,1000000);
	DDX_Text(pDX,IDC_L2,L2);
	DDV_MinMaxDouble(pDX,L2,0,1000000);
	DDX_Text(pDX,IDC_N1,N1);
	DDV_MinMaxInt(pDX,N1,0,1000);
	DDX_Text(pDX,IDC_N2,N2);
	DDV_MinMaxInt(pDX,N2,0,1000);
	DDX_Text(pDX,IDC_H1,H1);
	DDV_MinMaxDouble(pDX,H1,0,1000000);
	DDX_Text(pDX,IDC_H2,H2);
	DDV_MinMaxDouble(pDX,H2,0,1000000);
	DDX_Text(pDX,IDC_H3,H3);
	DDV_MinMaxDouble(pDX,H3,0,1000000);
	DDX_Control(pDX, IDC_TYPE, c_Type);
    DDX_Control(pDX, IDC_COVER1, c_Cover1);
	DDX_Control(pDX, IDC_COVER2, c_Cover2);
}
BOOL CNewTruss2DDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Type.AddString("Sloped truss");
	c_Type.AddString("Vertica truss");
	c_Type.AddString("Pratt truss");
	c_Type.SetCurSel(type);
	OnUpdateControls();
	return TRUE;
}
void CNewTruss2DDia::OnType() {
	type = c_Type.GetCurSel();
    OnUpdateControls();
}
void CNewTruss2DDia::OnUpdateControls() {
    c_Cover1.ShowWindow(type == PRATT);
    c_Cover2.ShowWindow(type != PRATT);
}
void CNewTruss2DDia::OnOK() {
	CDialog::OnOK();
	CmyDocument* pDoc = GetMyDocument();
	register int i,j;
	DOUBLE v; 

	if(type == SLOPED) {
		for(j = 0;j <= 2 * ndiv;j++) {
			v = j * divlength / 2 - ndiv * divlength / 2;
			pDoc->global->grid[IUX].AddTail(v);
		}
		v = 0;
		pDoc->global->grid[IUY].AddTail(v);
		pDoc->global->grid[IUZ].AddTail(v);
		v = height;
        pDoc->global->grid[IUZ].AddTail(v); 
		
		/*members and restraints*/
		JOINT j1,j2,j3,j4;
		for(i = 0;i < ndiv;i++) {
			j1.p.set(0,0,0);
			j1.p.x = i * divlength - ndiv * divlength / 2;
			j2 = j1;
			j2.p.x += divlength;
			j3 = j1;
			j3.p.x += divlength / 2;
			j3.p.z += height;
			j4 = j3;
			j4.p.x += divlength; 
			pDoc->AddMember(&j1,&j2);
			pDoc->AddMember(&j1,&j3);
			pDoc->AddMember(&j3,&j2);
			if(i != ndiv - 1) 
				pDoc->AddMember(&j3,&j4);
		}
	} else if(type == VERTICAL) {
		for(j = 0;j <= ndiv;j++) {
			v = j * divlength - ndiv * divlength / 2;
			pDoc->global->grid[IUX].AddTail(v);
		}
		v = 0;
		pDoc->global->grid[IUY].AddTail(v);
		pDoc->global->grid[IUZ].AddTail(v);
        pDoc->global->grid[IUZ].AddTail(height); 
		
		/*members and restraints*/
		JOINT j1,j2,j3,j4;
		for(i = 0;i < ndiv;i++) {
			j1.p.set(0,0,0);
			j1.p.x = i * divlength - ndiv * divlength / 2;
			j2 = j1;
			j2.p.x += divlength;
			j3 = j1;
			j3.p.z += height;
			j4 = j3;
			j4.p.x += divlength; 
			pDoc->AddMember(&j1,&j2);
			pDoc->AddMember(&j1,&j3);
			pDoc->AddMember(&j1,&j4);
			pDoc->AddMember(&j3,&j4);
			pDoc->AddMember(&j3,&j2);
		}
		pDoc->AddMember(&j2,&j4);
	} else if(type == PRATT) {
		for(j = 0;j <= N1;j++) {
			v = (j - N1) * L1 / N1;
			pDoc->global->grid[IUX].AddTail(v);
		}
		for(j = 0;j <= N2;j++) {
			v = j * L2 / N2;
			pDoc->global->grid[IUX].AddTail(v);
		}
		
		v = 0;
		pDoc->global->grid[IUY].AddTail(v);
		pDoc->global->grid[IUZ].AddTail(v);
        pDoc->global->grid[IUZ].AddTail(H1); 
        pDoc->global->grid[IUZ].AddTail(H2); 
        pDoc->global->grid[IUZ].AddTail(H3); 
		/*members and restraints*/
		JOINT j1,j2,j3,j4;
		for(i = 0;i < N1;i++) {
			j1.p.set(0,0,0);
			j1.p.x = -i * L1 / N1;
			j2 = j1;
			j2.p.x -= L1 / N1;
			j3 = j1;
			j3.p.z += H2 - i * (H2 - H1) / N1;
			j4 = j2;
			j4.p.z += H2 - (i + 1) * (H2 - H1) / N1;
			pDoc->AddMember(&j1,&j2);
			pDoc->AddMember(&j1,&j3);
			pDoc->AddMember(&j3,&j2);
			pDoc->AddMember(&j3,&j4);
		}
		pDoc->AddMember(&j2,&j4);
		/*right wing*/
		for(i = 0;i < N2;i++) {
			j1.p.set(0,0,0);
			j1.p.x = i * L2 / N2;
			j2 = j1;
			j2.p.x += L2 / N2;
			j3 = j1;
			j3.p.z += H2 - i * (H2 - H3) / N2;
			j4 = j2;
			j4.p.z += H2 - (i + 1) * (H2 - H3) / N2;
			pDoc->AddMember(&j1,&j2);
			pDoc->AddMember(&j1,&j3);
			pDoc->AddMember(&j3,&j2);
			pDoc->AddMember(&j3,&j4);
		}
		pDoc->AddMember(&j2,&j4);
	}
	pDoc->global->coordinate = CARTESIAN;

	if(restraint) {
		JOINT* joint;
		DOUBLE minx = pDoc->global->grid[IUX].GetHead(),
			   maxx = pDoc->global->grid[IUX].GetTail();
		POSITION pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(joint->p.z == 0) {
				if(joint->p.x == minx) joint->restraint = UX | UY | UZ | RX | RZ;
				else if(joint->p.x == maxx) joint->restraint = UY | UZ | RX | RZ;
			}
		}
	}
	CmyView::CREATION[1] = VIEW_XZ;
	CmyView::CREATION[2] = VIEW_XY;
}
/*
New truss 3D Dialog
*/
BEGIN_MESSAGE_MAP(CNewTruss3DDia , CDialog)
    ON_CBN_SELCHANGE(IDC_TYPE,OnType)
	ON_CBN_SELCHANGE(IDC_TOWERTYPE,OnTower)
    ON_LBN_SELCHANGE(IDC_LISTTH,OnUpdateLocTh)
	ON_LBN_SELCHANGE(IDC_LISTTW,OnUpdateLocTw)
	ON_LBN_SELCHANGE(IDC_LISTTAH,OnUpdateLocTah)
	ON_LBN_SELCHANGE(IDC_LISTTAW,OnUpdateLocTaw)
	ON_COMMAND(IDC_ADD,OnAdd)
    ON_COMMAND(IDC_MOVE,OnMove)
	ON_COMMAND(IDC_DELETE,OnDelete)
	ON_COMMAND(IDC_DELETEALL,OnDeleteAll)
END_MESSAGE_MAP()

CNewTruss3DDia::CNewTruss3DDia(CWnd* parent) : 
               CDialog(IDD_NEW_3DTRUSS,parent) {
	
    type = TELE_TOWER;
	tower_shape = RECTANGULAR;
	restraint = TRUE;
	Nx = 3;
	Ny = 3;
	Ndx = 6;
	Ndy = 6;
	Sx = 9;
	Sy = 9;
	H = 1.5;
}
void CNewTruss3DDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
	DDX_Text(pDX,IDC_NBAYX,Nx);
	DDV_MinMaxInt(pDX,Nx,0,1000);
	DDX_Text(pDX,IDC_NBAYY,Ny);
	DDV_MinMaxInt(pDX,Ny,0,1000);
	DDX_Text(pDX,IDC_NDIVX,Ndx);
	DDV_MinMaxInt(pDX,Ndx,0,1000);
	DDX_Text(pDX,IDC_NDIVY,Ndy);
	DDV_MinMaxInt(pDX,Ndy,0,1000);
	DDX_Text(pDX,IDC_SPANX,Sx);
	DDV_MinMaxDouble(pDX,Sx,0,1000000);
	DDX_Text(pDX,IDC_SPANY,Sy);
	DDV_MinMaxDouble(pDX,Sy,0,1000000);
    DDX_Text(pDX,IDC_HEIGHT,H);
	DDV_MinMaxDouble(pDX,H,0,1000000); 
	DDX_Control(pDX, IDC_TYPE, c_Type);
	DDX_Control(pDX, IDC_TOWERTYPE, c_Tower);
	DDX_Control(pDX, IDC_COVER1, c_Cover1);
	DDX_Control(pDX, IDC_COVER2, c_Cover2);
	for(int i = 0;i < 4;i++) {
		DDX_Control(pDX, IDC_TH + i, c_Value[i]);
        DDX_Control(pDX, IDC_LISTTH + i, c_List[i]);
	}
}
BOOL CNewTruss3DDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Type.AddString("Roof");
	c_Type.AddString("Tower");
	c_Type.SetCurSel(type);
	c_Value[0].SetWindowText("0.00");
	c_Value[1].SetWindowText("10.00");
	c_Value[2].SetWindowText("0.00");
	c_Value[3].SetWindowText("0.00");
	c_Tower.AddString("Tower type 1");
    c_Tower.AddString("Tower type 2");
	c_Tower.AddString("Tower type 3");
	c_Tower.AddString("Tower type 4");
	c_Tower.SetCurSel(0);
	OnUpdateControls();
	return TRUE;
}
void CNewTruss3DDia::OnType() {
	type = c_Type.GetCurSel();
    OnUpdateControls();
}
void CNewTruss3DDia::OnTower() {
	int i,k;
	int towertype = c_Tower.GetCurSel();
	tower.RemoveAll();
    for(i = 0;i < 4;i++) {
		c_List[i].ResetContent();
	}

	CString str;
	TOWER t;
	if(towertype == 0) {
		for(k = 0; k < 8; k++) {
			if(k <= 4)  {
				t.d[0] = k * 4;
				t.d[1] = 10 - k * 2;
				if(k == 0) {
					t.d[2] = 0;
					t.d[3] = 0;
				} else if(k == 4) {
					t.d[2] = 0.3;
					t.d[3] = 1.0;
				} else {
					t.d[2] = 0.3;
					t.d[3] = 0.3;
				}
			} else {
				t.d[0] = 16 + (k - 4) * 2;
				t.d[1] = 2;
				t.d[2] = 0;
				t.d[3] = 1;
			}
			tower.AddTail(t);
			for(i = 0;i < 4;i++) {
				str.Format("%.2f",t.d[i]);
				c_List[i].AddString(str);
			}
		}
	} else if(towertype == 1) {
		for(k = 0; k < 5; k++) {
			t.d[0] = k * 5;
			t.d[1] = 10 - k * 2;
			if(k == 0) t.d[3] = 0;
			else t.d[3] = 1.0;
			tower.AddTail(t);
			for(i = 0;i < 4;i++) {
				str.Format("%.2f",t.d[i]);
				c_List[i].AddString(str);
			}
		}
	} else if(towertype == 2) {
		for(k = 0; k < 5; k++) {
			t.d[0] = k * 5;
			t.d[1] = 4;
			if(k == 0) t.d[3] = 0;
			else t.d[3] = 0.5;
			tower.AddTail(t);
			for(i = 0;i < 4;i++) {
				str.Format("%.2f",t.d[i]);
				c_List[i].AddString(str);
			}
		}
	} else if(towertype == 3) {
		tower_shape = TRIANGULAR;
		for(k = 0; k < 5; k++) {
			t.d[0] = k * 5;
			t.d[1] = 4;
			t.d[2] = 0.3;
			if(k == 0) t.d[3] = 0;
			else t.d[3] = 0.5;
			tower.AddTail(t);
			for(i = 0;i < 4;i++) {
				str.Format("%.2f",t.d[i]);
				c_List[i].AddString(str);
			}
		}
	}
}
void CNewTruss3DDia::OnUpdateControls() {
	c_Cover1.ShowWindow(type == ROOF);
    c_Cover2.ShowWindow(type == TELE_TOWER);
}
void CNewTruss3DDia::OnOK() {
	CDialog::OnOK();

	CmyDocument* pDoc = GetMyDocument();
	int i,j,k;

	if(type == TELE_TOWER) {
		/*sort tower data*/
		POSITION pos,pos1;
		TOWER *v, *v1, temp;
		pos = tower.GetHeadPosition();
		if(!pos)
			return;
		while(pos) {
			v = &tower.GetNext(pos);

			pos1 = pos;
			while(pos1) {
				v1 = &tower.GetNext(pos1);
				if(v1->d[0] < v->d[0]) {
					temp = *v;
					*v = *v1;
					*v1 = temp;
				}
			}
		}
		/*grids*/
		SYSTEM system = *pDoc->global;
		DOUBLE g;

		if(tower_shape == TRIANGULAR) {
			pos = tower.GetHeadPosition();
			while(pos) {
				v = &tower.GetNext(pos);
				g = v->d[0];
				system.grid[IUZ].AddTail(g);
				g = v->d[1] / sqrt(3.0);
				system.grid[IUX].AddTail(g);
				for(i = 0;i <= 12;i++) {
					g = i * 30;
					system.grid[IUY].AddTail(g);
				}
			}
			system.coordinate = RADIAL;
		} else {
			pos = tower.GetHeadPosition();
			while(pos) {
				v = &tower.GetNext(pos);
				g = v->d[0];
				system.grid[IUZ].AddTail(g);
				g = -v->d[1] / 2;
				system.grid[IUX].AddTail(g);
				system.grid[IUY].AddTail(g);
				g = v->d[1] / 2;
				system.grid[IUX].AddTail(g);
				system.grid[IUY].AddTail(g);
			}
			system.coordinate = CARTESIAN;
		}
		*pDoc->global = system;
        /*Tower*/
		JOINT j[4],tj[4],bj[8],bhj[4];
		int N = (tower_shape == TRIANGULAR) ? 3 : 4;
		BOOL first = TRUE;

		pos = tower.GetHeadPosition();
		while(pos) {
			v = &tower.GetNext(pos);

			if(tower_shape == TRIANGULAR) {
				j[0].p.z = v->d[0];
				j[0].p.x = -v->d[1] / 2;
				j[0].p.y = -sqrt(3.0) * v->d[1] / 6;
				j[1].p = j[0].p;
				j[1].p.x += v->d[1];
				j[2].p = j[1].p;
				j[2].p.x = 0;
				j[2].p.y = v->d[1] / sqrt(3.0);
			} else {
				j[0].p.z = v->d[0];
				j[0].p.x = -v->d[1] / 2;
				j[0].p.y = -v->d[1] / 2;
				j[1].p = j[0].p;
				j[1].p.x += v->d[1];
				j[2].p = j[1].p;
				j[2].p.y += v->d[1];
				j[3].p = j[0].p;
				j[3].p.y += v->d[1];
			}
			if(!first) {
				for(i = 0;i < N;i++) {
					pDoc->AddMember(&tj[i],&j[i]);
					bhj[i].p = tj[i].p + (j[i].p - tj[i].p) * v->d[2];
				}
				for(i = 0;i < N;i++) { 
					k = (i == N - 1) ? 0 : (i + 1);
					bj[2 * i].p = j[i].p + (j[k].p - j[i].p) * v->d[3];
					pDoc->AddMember(&bhj[i],&bj[2 * i]);
					
					bj[2 * i + 1].p = j[k].p + (j[i].p - j[k].p) * v->d[3];
					pDoc->AddMember(&bhj[k],&bj[2 * i + 1]);
				}
				for(i = 0;i < N;i++) {
					k = (i == N - 1) ? 0 : (i + 1);
					pDoc->AddMember(&bj[2 * i + 1],&bj[2 * k]);
				}

				pDoc->AddMember(&j[0],&j[1]);
				pDoc->AddMember(&j[1],&j[2]);
				if(tower_shape == TRIANGULAR) {
					pDoc->AddMember(&j[2],&j[0]);
				} else {
					pDoc->AddMember(&j[2],&j[3]);
					pDoc->AddMember(&j[3],&j[0]);
				}
			}

            for(i = 0;i < N;i++)
				tj[i] = j[i];
			first = FALSE;
		}

		/*restraint*/
		DOUBLE base = tower.GetHead().d[0];
		if(restraint) {
			JOINT* joint;
			POSITION pos = pDoc->joints.GetHeadPosition();
			while(pos) {
				joint = &pDoc->joints.GetNext(pos);
				if(EQUAL(joint->p.z,base)) joint->restraint = UX | UY | UZ;
			}
		}
	} else if(type == ROOF) {
		/*grids*/
		DOUBLE v;
		for(i = 0;i <= (Nx * Ndx); i++) {
			v = i * Sx / Ndx - Nx * Sx / 2;
			pDoc->global->grid[IUX].AddTail(v);
		}
		for(j = 0;j <= (Ny * Ndy); j++) {
			v = j * Sy / Ndy - Ny * Sy / 2;
			pDoc->global->grid[IUY].AddTail(v);
		}
		v = 0;
		pDoc->global->grid[IUZ].AddTail(v);
		pDoc->global->coordinate = CARTESIAN;

		/*truss*/
		JOINT j1,j2,j3,j4,j5;
		for(i = 0;i <= (Nx * Ndx); i++) {
			for(j = 0;j <= (Ny * Ndy); j++) {
				j1.p.x = i * Sx / Ndx - Nx * Sx / 2;
				j1.p.y = j * Sy / Ndy - Ny * Sy / 2;
				j1.p.z = 0;
				if(i < (Nx * Ndx)) {
					j2.p = j1.p;
					j2.p.x += Sx / Ndx;
					pDoc->AddMember(&j1,&j2);
				}
				if(j < (Ny * Ndy)) {
					j2.p = j1.p;
					j2.p.y += Sy / Ndy;
					pDoc->AddMember(&j1,&j2);
				}
				j2.p = j1.p;
				j2.p.x -= Sx / (2 * Ndx);
				j2.p.y -= Sy / (2 * Ndy);
				j2.p.z = H;

				j3.p = j1.p;
				j3.p.x -= Sx / (2 * Ndx);
				j3.p.y += Sy / (2 * Ndy);
				j3.p.z = H;

				j4.p = j1.p;
				j4.p.x += Sx / (2 * Ndx);
				j4.p.y += Sy / (2 * Ndy);
				j4.p.z = H;

				j5.p = j1.p;
				j5.p.x += Sx / (2 * Ndx);
				j5.p.y -= Sy / (2 * Ndy);
				j5.p.z = H;

				pDoc->AddMember(&j1,&j2);
				pDoc->AddMember(&j1,&j3);
				pDoc->AddMember(&j1,&j4);
				pDoc->AddMember(&j1,&j5);

				pDoc->AddMember(&j2,&j3);
				pDoc->AddMember(&j3,&j4);
				pDoc->AddMember(&j4,&j5);
				pDoc->AddMember(&j5,&j2);
			}
		}

		/*restraint*/
		if(restraint) {
			JOINT j1,*joint;
			POSITION pos;
			for(i = 0;i <= Nx; i++) {
				for(j = 0;j <= Ny; j++) {
					j1.p.x = i * Sx - Nx * Sx / 2;
					j1.p.y = j * Sy - Ny * Sy / 2;
					j1.p.z = 0;
					pos = pDoc->joints.Find(j1);
					if(pos) {
						joint = &pDoc->joints.GetAt(pos);
						joint->restraint = UX | UY | UZ;
					}
				}
			}
		}
	}
}
void CNewTruss3DDia::OnAdd() {
    int i;
	TOWER t;
	CString str;
	for(i = 0;i < 4;i++) {
		c_Value[i].GetWindowText(str);
		t.d[i] = atof(str);
	}
	if(tower.Find(t))
		return;
	tower.AddTail(t);
	for(i = 0;i < 4;i++) {
		str.Format("%.2f",t.d[i]);
		c_List[i].SetCurSel(c_List[i].AddString(str));
	}
}
void CNewTruss3DDia::OnMove() {
	TOWER *pt;
	CString str;
	POSITION pos;
	int index = c_List[0].GetCurSel();
	if(index == LB_ERR)
		return;
	pos = tower.FindIndex(index);
	if(pos) {
		pt = &tower.GetAt(pos);
		for(int i = 0;i < 4;i++) {
			c_Value[i].GetWindowText(str);
			pt->d[i] = atof(str);
			str.Format("%.2f",pt->d[i]);
			c_List[i].DeleteString(index);
            c_List[i].InsertString(index,str);
            c_List[i].SetCurSel(index);
		}
	}
}
void CNewTruss3DDia::OnDelete() {
	CString str;
	POSITION pos;
	int index = c_List[0].GetCurSel();
	if(index == LB_ERR)
		return;
	pos = tower.FindIndex(index);
	if(pos) {
		tower.RemoveAt(pos);
		for(int i = 0;i < 4;i++) {
			c_List[i].DeleteString(index);
            c_List[i].SetCurSel(index);
		}
	}
}
void CNewTruss3DDia::OnDeleteAll() {
	tower.RemoveAll();
    for(int i = 0;i < 4;i++) {
		c_List[i].ResetContent();
	}
}
void CNewTruss3DDia::OnUpdateLoc(UINT nID) {
	CString str;
	int index = c_List[nID - IDC_TH].GetCurSel();
	if(index == LB_ERR)
		return;
	for(int i = 0;i < 4;i++) {
		c_List[i].SetCurSel(index);
		c_List[i].GetText(index,str);
		c_Value[i].SetWindowText(str);
	}
}
void CNewTruss3DDia::OnUpdateLocTh() {
	OnUpdateLoc(IDC_TH);
}
void CNewTruss3DDia::OnUpdateLocTw() {
	OnUpdateLoc(IDC_TW);
}
void CNewTruss3DDia::OnUpdateLocTah() {
	OnUpdateLoc(IDC_TAH);
}
void CNewTruss3DDia::OnUpdateLocTaw() {
	OnUpdateLoc(IDC_TAW);
}
/*
my add slab
*/
static SLAB* AddMySlab(CmyDocument* pDoc,RPoint p,DOUBLE w,DOUBLE h,BOOL usexz = FALSE) {
	JOINT j1,j2,j3,j4;
	j1.p = p;
	j2.p = j1.p;
	j2.p.x += w;
	j3.p = j2.p;
    if(usexz) j3.p.z += h;
	else j3.p.y += h;
	j4.p = j1.p;
	if(usexz) j4.p.z += h;
	else j4.p.y += h;
	SLAB* slab = pDoc->AddSlab(&j1,&j2,&j3,&j4);
	return slab;
}
/*
New frame Dialog
*/
BEGIN_MESSAGE_MAP(CNewFrameDia , CDialog)
    ON_LBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
END_MESSAGE_MAP()

CNewFrameDia::CNewFrameDia(CWnd* parent,int ttype) : 
               CDialog(IDD_NEW_3DFRAME,parent) {
	
	type = ttype;
	if(type == MODEL2DFRAME) {
		N[IUX] = 3;
	    D[IUX] = 6;
		N[IUY] = 0;
		D[IUY] = 0;
		N[IUZ] = 2;
	    D[IUZ] = 3;
		frame_type = PORTAL;
	} else {
		N[IUX] = 3;
	    D[IUX] = 6;
		N[IUY] = 2;
		D[IUY] = 6;
		N[IUZ] = 2;
	    D[IUZ] = 3;
		frame_type = OPEN_FRAME;
	}
	restraint = TRUE;
	gap = 1;
}
void CNewFrameDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NSTOREY, N[IUZ]);
	DDV_MinMaxInt(pDX,N[IUZ],0,1000);
	DDX_Text(pDX, IDC_SHEIGHT, D[IUZ]);
	DDV_MinMaxDouble(pDX,D[IUZ],0,1000000);
	DDX_Text(pDX, IDC_NBAYX, N[IUX]);
	DDV_MinMaxInt(pDX,N[IUX],0,1000);
	DDX_Text(pDX, IDC_BWIDTHX, D[IUX]);
	DDV_MinMaxDouble(pDX,D[IUX],0,1000000);
	if(type == MODEL3DFRAME) {
		DDX_Text(pDX, IDC_NBAYY, N[IUY]);
		DDV_MinMaxInt(pDX,N[IUY],0,1000);
		DDX_Text(pDX, IDC_BWIDTHY, D[IUY]);
		DDV_MinMaxDouble(pDX,D[IUY],0,1000000);
	}
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
	DDX_Control(pDX,IDC_COVER,c_Cover);
	DDX_Control(pDX,IDC_TYPE,c_FrameType);
	DDX_Control(pDX,IDC_GAPWIDTH,c_Gap);
	DDX_Text(pDX,IDC_GAPWIDTH,gap);
	DDV_MinMaxDouble(pDX,gap,0,1000000);
	DDX_Control(pDX,IDC_GAPTEXT,c_GapText);
}
BOOL CNewFrameDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Cover.ShowWindow(type == MODEL2DFRAME);
	if(type == MODEL2DFRAME) {
		c_FrameType.AddString("Portal");
		c_FrameType.AddString("Braced (Concentric)");
		c_FrameType.AddString("Braced (Eccentric)");
	} else {
		c_FrameType.AddString("Open Frame Building");
		c_FrameType.AddString("Perimeter Frame Building");
		c_FrameType.AddString("Beam-slab building");
		c_FrameType.AddString("Flat-plate building");
	}
	c_FrameType.SetCurSel(frame_type);
	OnUpdateControls();
	return TRUE;
}
void CNewFrameDia::OnUpdateControls() {
	frame_type = c_FrameType.GetCurSel();
	c_Gap.ShowWindow(type == MODEL2DFRAME && frame_type == BRACED_ECCENTRIC);
    c_GapText.ShowWindow(type == MODEL2DFRAME && frame_type == BRACED_ECCENTRIC);
}
void CNewFrameDia::OnOK() {
	frame_type = c_FrameType.GetCurSel();
	CDialog::OnOK();
	CmyDocument* pDoc = GetMyDocument();

	BOOL concentric = (type == MODEL2DFRAME && frame_type == BRACED_CONCENTRIC); 
    BOOL eccentric = (type == MODEL2DFRAME && frame_type == BRACED_ECCENTRIC); 
	BOOL beamslabf = (type == MODEL3DFRAME && frame_type == BEAM_SLAB);
	BOOL perimeterf = (type == MODEL3DFRAME && frame_type == PERIMETER_FRAME);
	BOOL flatplatef = (type == MODEL3DFRAME && frame_type == FLAT_PLATE);

	register int i,j,k;
	DOUBLE v;
	for(i = IUX;i <= IUZ ;i++) {
		for(j = 0;j <= N[i];j++) {
			v = j * D[i];
			if(i != IUZ) v -= N[i] * D[i] / 2;
			pDoc->global->grid[i].AddTail(v);

			if(eccentric && i == IUX && j < N[i] && (j % 2)) {
				v += (D[i] - gap) / 2; 
				pDoc->global->grid[i].AddTail(v);
				v += gap; 
				pDoc->global->grid[i].AddTail(v);
			}
		}
	}
	pDoc->global->coordinate = CARTESIAN;
	/*members and restraints*/
	for(i = 0;i <= N[IUX];i++) {
		for(j = 0;j <= N[IUY];j++) {
			if(perimeterf) {
				if(i == 0 || i == N[IUX] || j == 0 || j == N[IUY]);
				else continue;
			}
			for(k = 0;k <= N[IUZ];k++) {
				JOINT j1,j2,j3;
				j1.p.x = i * D[IUX] - N[IUX] * D[IUX] / 2;
				j1.p.y = j * D[IUY] - N[IUY] * D[IUY] / 2;
				j1.p.z = k * D[IUZ];
				
				if(k > 0 && !flatplatef) {
					if(i < N[IUX] && (!perimeterf || (j == 0 || j == N[IUY]))) {
						j2.p = j1.p;
						j2.p.x += D[IUX];
						pDoc->AddMember(&j1,&j2);
					}
					if(j < N[IUY] && (!perimeterf || (i == 0 || i == N[IUX]))) {
						j2.p = j1.p;
						j2.p.y += D[IUY];
						pDoc->AddMember(&j1,&j2);
					}
				}
				if(k < N[IUZ]) {
					j2.p = j1.p;
					j2.p.z += D[IUZ];
					pDoc->AddMember(&j1,&j2);
					if(i < N[IUX] && (i % 2)) {
						if(concentric) {
							j3.p = j1.p;
							j3.p.x += D[IUX];
							pDoc->AddMember(&j2,&j3);
							j2.p.x += D[IUX];
							pDoc->AddMember(&j1,&j2);
						} else if(eccentric) {
							j3.p = j2.p;
							j3.p.x += (D[IUX] - gap) / 2;
							pDoc->AddMember(&j1,&j3);
							j3.p.x += gap;
							j2.p = j1.p;
							j2.p.x += D[IUX];
							pDoc->AddMember(&j3,&j2);
						}
					}
					
				}
			}
		}
	}
	if(restraint) {
		JOINT* joint;
		POSITION pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(joint->p.z == 0) 
				joint->restraint = ALLUR;
		}
	}
	
	if(beamslabf || flatplatef) {
		for(k = 1;k <= N[IUZ];k++) {
			for(i = 0;i < N[IUX];i++) {
				for(j = 0;j < N[IUY];j++) {
					RPoint p;
					p.x = i * D[IUX] - N[IUX] * D[IUX] / 2;
					p.y = j * D[IUY] - N[IUY] * D[IUY] / 2;
					p.z = k * D[IUZ];
					SLAB* slab = AddMySlab(pDoc,p,D[IUX],D[IUY]);
					slab->nDivx = 4;
					slab->nDivy = 4;
				}
			}
		}
	}
	
	if(type == MODEL2DFRAME) {
		CmyView::CREATION[1] = VIEW_XZ;
		CmyView::CREATION[2] = VIEW_XY;
	}
}
/*
New Slab Dialog
*/
BEGIN_MESSAGE_MAP(CNewSlabDia , CDialog)
END_MESSAGE_MAP()

CNewSlabDia::CNewSlabDia(CWnd* parent) : 
CDialog(IDD_NEW_SLAB,parent) {
	
	N[IUX] = 3;
	D[IUX] = 8;
	N[IUY] = 3;
	D[IUY] = 8;
	M[IUX] = 4;
	M[IUY] = 4;
	restraint = TRUE;
}
void CNewSlabDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NBAYX, N[IUX]);
	DDV_MinMaxInt(pDX,N[IUX],0,1000);
	DDX_Text(pDX, IDC_BWIDTHX, D[IUX]);
	DDV_MinMaxDouble(pDX,D[IUX],0,1000000);
	DDX_Text(pDX, IDC_NBAYY, N[IUY]);
	DDV_MinMaxInt(pDX,N[IUY],0,1000);
	DDX_Text(pDX, IDC_BWIDTHY, D[IUY]);
	DDV_MinMaxDouble(pDX,D[IUY],0,1000000);
	DDX_Text(pDX, IDC_MIDX, M[IUX]);
	DDV_MinMaxDouble(pDX,M[IUX],0,1000000);
	DDX_Text(pDX, IDC_MIDY, M[IUY]);
	DDV_MinMaxDouble(pDX,M[IUY],0,1000000);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
}
BOOL CNewSlabDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	return TRUE;
}
void CNewSlabDia::OnOK() {
	CDialog::OnOK();

	CmyDocument* pDoc = GetMyDocument();
	register int i,j,k,l;
	DOUBLE v;
	for(i = IUX;i <= IUY ;i++) {
		for(j = 0;j <= N[i];j++) {
			v = j * D[i];
			v -= N[i] * D[i] / 2;
			pDoc->global->grid[i].AddTail(v);
			if(j < N[i]) {
				v += (D[i] - M[i]) / 2.0;
				pDoc->global->grid[i].AddTail(v);
				v += M[i];
				pDoc->global->grid[i].AddTail(v);
			}
		}
	}
	v = 0;
	pDoc->global->grid[IUZ].AddTail(v);
	pDoc->global->coordinate = CARTESIAN;

	/*slabs*/
	RPoint p,pt;
	DOUBLE w,h,wt,ht;
	for(i = 0;i < N[IUX];i++) {
		for(j = 0;j < N[IUY];j++) {

			p.x = i * D[IUX] - N[IUX] * D[IUX] / 2;
			p.y = j * D[IUY] - N[IUY] * D[IUY] / 2;
			p.z = 0;
			pt = p;

			wt = 0;
			for(k = 0;k < 6; k++) {
				if(k < 2) w = (D[IUX] - M[IUX]) / 4;
				else if(k < 4) w = (M[IUX]) / 2;
				else w = (D[IUX] - M[IUX]) / 4;

				pt.x = p.x + wt;
				wt += w;

				ht = 0;
				for(l = 0;l < 6;l++) {
					if(l < 2) h = (D[IUY] - M[IUY]) / 4;
					else if(l < 4) h = (M[IUY]) / 2;
					else h = (D[IUY] - M[IUY]) / 4;

				    pt.y = p.y + ht;	
					ht += h;

					AddMySlab(pDoc,pt,w,h);
				}
			}
		}
	}
	/*restraints*/
	if(restraint) {
		POSITION pos;
		JOINT* joint;
		for(i = 0;i <= N[IUX];i++) {
			for(j = 0;j <= N[IUY];j++) {
				p.x = i * D[IUX] - N[IUX] * D[IUX] / 2;
				p.y = j * D[IUY] - N[IUY] * D[IUY] / 2;
				p.z = 0;
				
				pos = pDoc->joints.GetHeadPosition();
				while(pos) {
					joint = &pDoc->joints.GetNext(pos);
					if(joint->p == p) break;
				}
				joint->restraint = ALLU;
			}
		}
	}
}
/*
New Wall Dialog
*/
BEGIN_MESSAGE_MAP(CNewWallDia , CDialog)
END_MESSAGE_MAP()

CNewWallDia::CNewWallDia(CWnd* parent) : 
CDialog(IDD_NEW_WALL,parent) {
	
	N[IUX] = 6;
	D[IUX] = 1;
	N[IUY] = 6;
	D[IUY] = 1;
	restraint = TRUE;
}
void CNewWallDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NSTOREY, N[IUY]);
	DDV_MinMaxInt(pDX,N[IUY],0,1000);
	DDX_Text(pDX, IDC_SHEIGHT, D[IUY]);
	DDV_MinMaxDouble(pDX,D[IUY],0,1000000);
	DDX_Text(pDX, IDC_NBAYX, N[IUX]);
	DDV_MinMaxInt(pDX,N[IUX],0,1000);
	DDX_Text(pDX, IDC_BWIDTHX, D[IUX]);
	DDV_MinMaxDouble(pDX,D[IUX],0,1000000);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
}
BOOL CNewWallDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	return TRUE;
}
void CNewWallDia::OnOK() {
	CDialog::OnOK();

	CmyDocument* pDoc = GetMyDocument();
	register int i,j;
	DOUBLE v;
	for(j = 0;j <= N[IUX];j++) {
		v = j * D[IUX];
		v -= N[IUX] * D[IUX] / 2;
		pDoc->global->grid[IUX].AddTail(v);
	}
    for(j = 0;j <= N[IUY];j++) {
		v = j * D[IUY];
		pDoc->global->grid[IUZ].AddTail(v);
	}
	v = 0;
	pDoc->global->grid[IUY].AddTail(v);
	pDoc->global->coordinate = CARTESIAN;


	RPoint p;
	for(i = 0;i < N[IUX];i++) {
		for(j = 0;j < N[IUY];j++) {
			p.x = i * D[IUX] - N[IUX] * D[IUX] / 2;
			p.y = 0;
			p.z = j * D[IUY];
			AddMySlab(pDoc,p,D[IUX],D[IUY],TRUE);
		}
	}

	if(restraint) {
		JOINT* joint;
		POSITION pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(joint->p.z == 0) 
				joint->restraint = ALLUR;
		}
	}

	CmyView::CREATION[1] = VIEW_XZ;
	CmyView::CREATION[2] = VIEW_XY;
}
/*
New Shell Dialog
*/
BEGIN_MESSAGE_MAP(CNewShellDia , CDialog)
    ON_LBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
END_MESSAGE_MAP()

CNewShellDia::CNewShellDia(CWnd* parent) : 
CDialog(IDD_NEW_SHELL,parent) {
	
	N[IUX] = 4;
	D[IUX] = 4;
	N[IUY] = 16;
	D[IUY] = 2;
	restraint = TRUE;
	type = CYLINDER;
}
void CNewShellDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NSTOREY, N[IUX]);
	DDV_MinMaxInt(pDX,N[IUX],0,1000);
	DDX_Text(pDX, IDC_SHEIGHT, D[IUX]);
	DDV_MinMaxDouble(pDX,D[IUX],0,1000000);
	DDX_Text(pDX, IDC_NANGULAR, N[IUY]);
	DDV_MinMaxInt(pDX,N[IUY],0,1000);
	DDX_Text(pDX, IDC_RADIUS, D[IUY]);
	DDV_MinMaxDouble(pDX,D[IUY],0,1000000);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Check(pDX, IDC_RESTRAINT, restraint);
}
BOOL CNewShellDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Type.AddString("Cylinder");
	c_Type.SetCurSel(type);

	OnUpdateControls();
	return TRUE;
}
void CNewShellDia::OnOK() {
	CDialog::OnOK();

	CmyDocument* pDoc = GetMyDocument();
	register int i,j,k;
	DOUBLE Zinc;

	/*grids*/
	if(type == CYLINDER) {
		Zinc = D[IUX] / N[IUX];
		
		DOUBLE v;
		for(j = 0;j <= N[IUX];j++) {
			v = j * Zinc;
			pDoc->global->grid[IUZ].AddTail(v);
		}
		for(j = 0;j <= N[IUY];j++) {
			v = j * (360.0 / N[IUY]);
			pDoc->global->grid[IUY].AddTail(v);
		}
		v = 0;
		pDoc->global->grid[IUX].AddTail(v);
		v = D[IUY];
		pDoc->global->grid[IUX].AddTail(v);
		pDoc->global->coordinate = RADIAL;
	} else {
	
	}

	/*add shells*/
	SLAB* slab;
	RPoint p[4];
	
	for(i = 0;i < N[IUX];i++) {
		for(j = 0;j < N[IUY];j++) {
			 JOINT jt[4];
             p[0].x = D[IUY];
			 p[0].y = j * (360.0 / N[IUY]);
			 p[0].z = i * Zinc;
			 p[1].x = D[IUY];
			 p[1].y = (j + 1) * (360.0 / N[IUY]);
			 p[1].z = i * Zinc;
			 p[2].x = D[IUY];
			 p[2].y = (j + 1) * (360.0 / N[IUY]);
			 p[2].z = (i + 1) * Zinc;
			 p[3].x = D[IUY];
			 p[3].y = j * (360.0 / N[IUY]);
			 p[3].z = (i + 1) * Zinc;
			 for(k = 0;k < 4;k++) {
				 jt[k].p = RDtoCT(p[k]);
			 }
			 
			 slab = pDoc->AddSlab(jt,4);
		}
	}

	if(restraint) {
		JOINT* joint;
		POSITION pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(joint->p.z == 0) 
				joint->restraint = ALLU;
		}
	}

	CmyView::CREATION[1] = VIEW_XZ;
	CmyView::CREATION[2] = VIEW_XY;
}
void CNewShellDia::OnUpdateControls() {
	type = c_Type.GetCurSel();
}
/*
Analysis Dialog
*/
BEGIN_MESSAGE_MAP(CAnalysisDia , CDialog)
    ON_LBN_SELCHANGE(IDC_LIST,OnUpdateLoc)
	ON_CBN_SELCHANGE(IDC_STIFFNESS,OnUpdateStiffness)
   	ON_COMMAND(IDC_RUN,OnRun)
    ON_COMMAND(IDC_CHANGEACTION,OnChangeAction)
END_MESSAGE_MAP()

CAnalysisDia::CAnalysisDia(CWnd* parent,PANALYSISCASELIST* pcase) : 
          CDialog(IDD_ANALYSIS,parent) {
	analysis_cases = pcase;
}
void CAnalysisDia::OnOK() {
    CDialog::OnOK();
}
void CAnalysisDia::OnRun() {
	EndDialog(IDC_RUN);
}
BOOL CAnalysisDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	ANALYSISCASE* pcase;
	CString str;
	POSITION pos;

	for(int i = 0;i < CASETYPES;i++) {
		if(i == ANALYSISCASE::COMBO_CASE) continue;
		pos = analysis_cases[i]->GetHeadPosition();
		while(pos) {
			pcase = &analysis_cases[i]->GetNext(pos);
			str.Format("%s",pcase->name);
			c_List.AddString(str);
		}
	}
	c_List.SetCurSel(0);

	OnUpdateLoc();
	return TRUE;
}
void CAnalysisDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
    DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_ACTION,c_Action);
	DDX_Control(pDX,IDC_STATUS,c_Status);
	DDX_Control(pDX,IDC_STIFFNESS,c_Stiffness);
}
ANALYSISCASE* CAnalysisDia::find_case() {
	POSITION pos;
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		ANALYSISCASE acase;
		ANALYSISCASE* pcase;
		c_List.GetText(index,acase.name);

		for(int i = 0;i < CASETYPES;i++) {
			if(i == ANALYSISCASE::COMBO_CASE) continue;
			pos = analysis_cases[i]->Find(acase);
			if(pos) {
                pcase = &analysis_cases[i]->GetAt(pos);
			    break;
			}
		}
		return pcase;
	}
	return NULL;
}
void CAnalysisDia::OnUpdateLoc() {
	ANALYSISCASE* pcase = find_case();
	if(pcase) {
		CString str;
		str = pcase->finished ? "Finished" : "Not finished";
		c_Status.SetWindowText(str);
		str = pcase->run ? "Run" : "Not run";
		c_Action.SetWindowText(str);

		NLCASE* pnl_case;
		c_Stiffness.ResetContent();
		c_Stiffness.AddString("Unstressed state");
		POSITION pos = analysis_cases[5]->GetHeadPosition();
		while(pos) {
			pnl_case = (NLCASE*)&analysis_cases[5]->GetNext(pos);
			if(pnl_case->name == pcase->name) continue;
			if(pnl_case->save_stiffness) {
				str.Format("%s",pnl_case->name);
				c_Stiffness.AddString(str);
			}
		}

		if(pcase->start_case == "") {
			c_Stiffness.SetCurSel(0);
		} else {
			c_Stiffness.SetCurSel(c_Stiffness.FindString(-1,pcase->start_case));
		}
        c_Stiffness.ShowWindow(pcase->acase_type != ANALYSISCASE::RESPONSES_CASE);
	}
}
void CAnalysisDia::OnChangeAction() {
	ANALYSISCASE* pcase = find_case();
	if(pcase) {
		pcase->run = !pcase->run;
		CString str;
		str = pcase->run ? "Run" : "Not run";
		c_Action.SetWindowText(str);
	}
}
void CAnalysisDia::OnUpdateStiffness() {
	ANALYSISCASE* pcase = find_case();
	if(pcase) {
		int index = c_Stiffness.GetCurSel();
		if(index == 0) pcase->start_case = "";
		else c_Stiffness.GetLBText(index,pcase->start_case);
	}
}
/*
Custom list box
*/
void CmyListBox::DrawItem(LPDRAWITEMSTRUCT lpDIS) {
	CDC* pDC = CDC::FromHandle(lpDIS->hDC);
	COLORREF cr = RGB(0,0,0);
	CString str;
	GetText(lpDIS->itemID,str);
	str.MakeUpper();

	if (lpDIS->itemAction & ODA_DRAWENTIRE) {
		CBrush br(cr);
		pDC->FillRect(&lpDIS->rcItem, &br);
		COLORREF bkmode = pDC->SetBkMode(TRANSPARENT);
		COLORREF col = pDC->SetTextColor(RGB(0,255,0));
		pDC->TextOut(lpDIS->rcItem.left,lpDIS->rcItem.top,str);
		pDC->SetTextColor(col);
		pDC->SetBkMode(bkmode);
	}

	if ((lpDIS->itemState & ODS_SELECTED) && (lpDIS->itemAction & (ODA_SELECT | ODA_DRAWENTIRE))) {
		COLORREF crHilite = RGB(0,0,255);
		CBrush br(crHilite);
		pDC->FrameRect(&lpDIS->rcItem, &br);
	}

	if (!(lpDIS->itemState & ODS_SELECTED) && (lpDIS->itemAction & ODA_SELECT)) {
		CBrush br(cr);
		pDC->FrameRect(&lpDIS->rcItem, &br);
	}
}
/*
Timer Update Dialog
*/
BEGIN_MESSAGE_MAP(CTimerUpdateDia , CDialog)
    ON_WM_TIMER()
END_MESSAGE_MAP()

CTimerUpdateDia::CTimerUpdateDia(CWnd* parent,UINT nID,MEMBER* membert,CmyView* pViewt) :
           CDialog(nID,parent)
{
	member = membert;
	pView = pViewt;
	blink = FALSE;
}
BOOL CTimerUpdateDia::OnInitDialog() {
	CDialog::OnInitDialog();
	myTimer = SetTimer(1, 500, 0);
	
	return TRUE;
}
void CTimerUpdateDia::OnTimer(UINT nIDEvent) {
	CDialog::OnTimer(nIDEvent);
    blink = !blink;
	CPoint p1,p2;
	CPen pen(PS_SOLID,2,blink ? RGB(255,0,0) : RGB(255,255,0));
	CPen* oldPen;
	CClientDC clientDC(pView);
	oldPen = clientDC.SelectObject(&pen);
	if(member->is_curved) {
		POSITION pos1 = member->curved_list.GetHeadPosition();
		p1 = pView->RPtoCP(member->curved_list.GetNext(pos1) + member->j1->p);
		while(pos1) {
			p2 = pView->RPtoCP(member->curved_list.GetNext(pos1) + member->j1->p);
			clientDC.MoveTo(p1);
			clientDC.LineTo(p2);
			p1 = p2;
		}
	} else {
		p1 = pView->RPtoCP(member->j1->p);
		p2 = pView->RPtoCP(member->j2->p);
		clientDC.MoveTo(p1);
		clientDC.LineTo(p2);
	}
	
	clientDC.SelectObject(oldPen);
}
/*
Analysis UpdateDialog
*/
BEGIN_MESSAGE_MAP(CAnalysisProgressDia , CDialog)
    ON_WM_TIMER()
END_MESSAGE_MAP()

CAnalysisProgressDia::CAnalysisProgressDia(CWnd* parent,CEvent* pEventt) : 
          CDialog(IDD_ANALYSIS_PROGRESS,parent) {
	pEvent = pEventt;
}
void CAnalysisProgressDia::AnalysisEnd() {
	KillTimer(myTimer);
	MessageBeep(0xFFFFFFFF);
	CTime finish = CTime::GetCurrentTime();
	c_End.SetWindowText(finish.Format("%Y/%m/%d    %H:%M:%S"));
    c_Ok.EnableWindow(TRUE);
	c_Cancel.EnableWindow(FALSE);
}
BOOL CAnalysisProgressDia::OnInitDialog() {
	CDialog::OnInitDialog();
	myTimer = SetTimer(1, 1000, 0);
	pEvent->SetEvent();
	beginTime = CTime::GetCurrentTime();
	c_Begin.SetWindowText(beginTime.Format("%Y/%m/%d    %H:%M:%S"));
	c_End.SetWindowText("Not applicable");
	c_Elapsed.SetWindowText("00:00:00");
	c_Ok.EnableWindow(FALSE);
	c_Cancel.EnableWindow(TRUE);
	return TRUE;
}
void CAnalysisProgressDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
    DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_BTIME,c_Begin);
	DDX_Control(pDX,IDC_ETIME,c_End);
	DDX_Control(pDX,IDC_ELAPSED,c_Elapsed);
    DDX_Control(pDX,IDOK,c_Ok);
	DDX_Control(pDX,IDCANCEL,c_Cancel);
}
void CAnalysisProgressDia::OnTimer(UINT nIDEvent) {
    CTimeSpan elapsed = CTime::GetCurrentTime() - beginTime;
	c_Elapsed.SetWindowText(elapsed.Format("%H:%M:%S"));
	CDialog::OnTimer(nIDEvent);
}
/*
Grid Dialog
*/
CSystemPage::CSystemPage(UINT nID,SYSTEM* psystem) : 
          CPropertyPage(nID) {
    system = psystem;
	cartesian = (nID == IDD_CARTESIAN_GRID);
	if(cartesian) {
		ngrid[IUX] = 4; ngrid[IUY] = 4; ngrid[IUZ] = 5;
		grids[IUX] = 6; grids[IUY] = 6; grids[IUZ] = 3;
	} else {
		ngrid[IUX] = 4; ngrid[IUY] = 7; ngrid[IUZ] = 5;
		grids[IUX] = 6; grids[IUY] = 30; grids[IUZ] = 3;
	}
}
void CSystemPage::DoDataExchange(CDataExchange* pDX) {
	CPropertyPage::DoDataExchange(pDX);

    DDX_Text(pDX, IDC_NXGRID, ngrid[IUX]);
    DDV_MinMaxInt(pDX,ngrid[IUX],0,1000);

	DDX_Text(pDX, IDC_NYGRID, ngrid[IUY]);
    DDV_MinMaxInt(pDX,ngrid[IUY],0,1000);

	DDX_Text(pDX, IDC_NZGRID, ngrid[IUZ]);
    DDV_MinMaxInt(pDX,ngrid[IUZ],0,1000);

	DDX_Text(pDX, IDC_XGRID, grids[IUX]);
    DDV_MinMaxDouble(pDX,grids[IUX],0,500000);

	DDX_Text(pDX, IDC_YGRID, grids[IUY]);
    DDV_MinMaxDouble(pDX,grids[IUY],0,500000);

	DDX_Text(pDX, IDC_ZGRID, grids[IUZ]);
    DDV_MinMaxDouble(pDX,grids[IUZ],0,500000);
}
void CSystemPage::OnOK() {
	CPropertySheet* parent = (CPropertySheet*) GetParent();
	if(parent->GetActivePage() == this) {
		register int i,j;
		SYSTEM mysystem;
		DOUBLE v;
		for(i = IUX;i <= IUZ ;i++) {
			for(j = 0;j < ngrid[i];j++) {
				v = j * grids[i];
				if(!cartesian && i == IUY && v > 360) {
					v = v - int(v / 360) * 360;
				}
				mysystem.grid[i].AddTail(v);
			}
		}
		if(cartesian) mysystem.coordinate = CARTESIAN;
		else mysystem.coordinate = RADIAL;
		*system = mysystem;
	}
	CPropertyPage::OnOK();
}
/*
property sheet
*/
CNewSystemSht::CNewSystemSht(CWnd* parent,SYSTEM* psystem)
	: CPropertySheet("New Coordinate System", parent, 0)
	, cartesian(IDD_CARTESIAN_GRID,psystem)
    , radial(IDD_RADIAL_GRID,psystem)
{
	AddPage(&cartesian);
    AddPage(&radial);
}
/*
System Origin Dialog
*/
CSystemOriginDia::CSystemOriginDia(CWnd* parent,SYSTEM* psystem) : 
                    CDialog(IDD_SYSTEM_ORIGIN,parent) {
    system = psystem;
	origin = system->origin;
	rotation = system->rotation;
}

void CSystemOriginDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX,IDC_UX,origin.x);
	DDX_Text(pDX,IDC_UY,origin.y);
	DDX_Text(pDX,IDC_UZ,origin.z);
	DDX_Text(pDX,IDC_RX,rotation.x);
	DDX_Text(pDX,IDC_RY,rotation.y);
	DDX_Text(pDX,IDC_RZ,rotation.z);
}
void CSystemOriginDia::OnOK() {
	CDialog::OnOK();
	system->origin = origin;
	system->rotation = rotation;
}
/*
Edit Grid Dialog
*/
BEGIN_MESSAGE_MAP(CSystemPropertyDia , CDialog)
    ON_LBN_SELCHANGE(IDC_LOCLIST,OnUpdateLoc)
	ON_BN_CLICKED(IDC_X,OnUpdateControls)
    ON_BN_CLICKED(IDC_Y,OnUpdateControls)
    ON_BN_CLICKED(IDC_Z,OnUpdateControls)
	ON_COMMAND(IDC_ADD,OnAddGrid)
    ON_COMMAND(IDC_MOVE,OnMoveGrid)
	ON_COMMAND(IDC_DELETE,OnDeleteGrid)
	ON_COMMAND(IDC_DELETEALL,OnDeleteAll)
    ON_COMMAND(IDC_NEW,OnNew)
	ON_COMMAND(IDC_ORIGIN,OnOrigin)
END_MESSAGE_MAP()

CSystemPropertyDia::CSystemPropertyDia(CWnd* parent,SYSTEM* tsystem) : 
              CDialog(IDD_EDITGRID,parent) {
	selected = IUX;
	system = *tsystem;
	psystem = tsystem;
}
void CSystemPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_NAME,system.name);
	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Control(pDX,IDC_Z,c_Z);
	DDX_Control(pDX,IDC_Y,c_Y);
	DDX_Control(pDX,IDC_X,c_X);
	DDX_Control(pDX,IDC_LOC,c_Location);
	DDX_Control(pDX,IDC_LOCLIST,c_List);
	DDX_Control(pDX,IDC_ORIGIN,c_Origin);
}
void CSystemPropertyDia::OnUpdateControls() {
	if(c_X.GetCheck() == BST_CHECKED) selected = IUX;
	else if(c_Y.GetCheck() == BST_CHECKED) selected = IUY;
	else selected = IUZ;
	
	CString str;
	c_List.ResetContent();
	POSITION pos = system.grid[selected].GetHeadPosition();
	while(pos) {
		str.Format("%+.2f",system.grid[selected].GetNext(pos));
		c_List.AddString(str);
	}
	c_List.SetCurSel(0);
	OnUpdateLoc();
}
void CSystemPropertyDia::OnUpdateLoc() {
	CString text;
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		c_List.GetText(index,text);
		c_Location.SetWindowText(text);
	}
}
BOOL CSystemPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	if(system.coordinate == RADIAL) {
		SetWindowText("Radial Grid");
		c_X.SetWindowText("R");
		c_Y.SetWindowText("Q");
		c_Z.SetWindowText("Z");
	} else {
		SetWindowText("Cartesian Grid");
		c_X.SetWindowText("X");
		c_Y.SetWindowText("Y");
		c_Z.SetWindowText("Z");
	}

	if(system.rank) {
		c_Name.EnableWindow(FALSE);
		c_Origin.ShowWindow(FALSE);
	}

	if(selected == IUX) c_X.SetCheck(TRUE);
	else if(selected == IUY) c_Y.SetCheck(TRUE);
	else  c_Z.SetCheck(TRUE);
	
    CString str;
	POSITION pos = system.grid[selected].GetHeadPosition();
	while(pos) {
		str.Format("%+.2f",system.grid[selected].GetNext(pos));
		c_List.AddString(str);
	}
	c_List.SetCurSel(0);
	OnUpdateLoc();
	return TRUE;
}
void CSystemPropertyDia::OnAddGrid() {
	DOUBLE d;
	CString str;
	c_Location.GetWindowText(str);
	d = atof(str);
	if(!system.grid[selected].Find(d)) {
		system.grid[selected].AddTail(d);
		str.Format("%+.2f",d);
		c_List.SetCurSel(c_List.AddString(str));
	}
}
void CSystemPropertyDia::OnMoveGrid() {
	CString str;
	POSITION pos;
	DOUBLE d;
	int index = c_List.GetCurSel();
	c_List.GetText(index,str);
	d = atof(str);
	if(pos = system.grid[selected].Find(d)) {
		c_Location.GetWindowText(str);
		d = atof(str);
		str.Format("%+.2f",d);
		system.grid[selected].GetAt(pos) = d;
		c_List.DeleteString(index);
		c_List.InsertString(index,str);
		c_List.SetCurSel(index);
	}
}
void CSystemPropertyDia::OnDeleteGrid() {
	CString str;
	POSITION pos;
	DOUBLE d;
	int index = c_List.GetCurSel();
	c_List.GetText(index,str);
	d = atof(str);
	if(pos = system.grid[selected].Find(d)) {
		system.grid[selected].RemoveAt(pos);
		c_List.DeleteString(index);
		c_List.SetCurSel(index);
	}
}
void CSystemPropertyDia::OnDeleteAll() {
	c_List.ResetContent();
	system.grid[selected].RemoveAll();
	c_List.SetCurSel(0);
}
void CSystemPropertyDia::OnNew() {
	SYSTEM sys;
	CNewSystemSht Dia(this,&sys);
	if(Dia.DoModal() == IDOK) {
		c_List.ResetContent();
		system = sys;
		OnUpdateControls();
	}
}
void CSystemPropertyDia::OnOrigin() {
	CSystemOriginDia Dia(this,&system);
	Dia.DoModal();
}
void CSystemPropertyDia::OnOK() {
	CDialog::OnOK();
	if(system.coordinate == RADIAL) {
		DOUBLE* v;
		POSITION pos = system.grid[IUY].GetHeadPosition();
        while(pos) {
			v = &system.grid[IUY].GetNext(pos);
			if(*v > 360) {
				*v = *v - int(*v / 360) * 360;
			}
		}
	}
	*psystem = system;
}
/*
Group Property Dialog
*/
CGroupPropertyDia::CGroupPropertyDia(CWnd* parent,GROUP* pgroup) : 
                    CDialog(IDD_GROUP_PROPERTY,parent) {
    group = pgroup;
}
BOOL CGroupPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	if(group->rank) {
		c_Name.EnableWindow(FALSE);
	}
	
	POSITION pos;
	JOINT* joint;
	MEMBER* mem;
	SLAB* sla;

    pos = group->jointlist.GetHeadPosition();
	while(pos) {
		joint = group->jointlist.GetNext(pos);
		c_JList.AddString(joint->name);
	}
	c_JList.SetCurSel(0);
	pos = group->memberlist.GetHeadPosition();
	while(pos) {
		mem = group->memberlist.GetNext(pos);
		c_MList.AddString(mem->name);
	}
	c_MList.SetCurSel(0);
	pos = group->slablist.GetHeadPosition();
	while(pos) {
		sla = group->slablist.GetNext(pos);
		c_SList.AddString(sla->name);
	}
	c_SList.SetCurSel(0);

	return TRUE;
}
void CGroupPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Control(pDX,IDC_JLIST,c_JList);
	DDX_Control(pDX,IDC_MLIST,c_MList);
	DDX_Control(pDX,IDC_SLIST,c_SList);
	DDX_Text(pDX,IDC_NAME,group->name);
}
void CGroupPropertyDia::OnOK() {
	CDialog::OnOK();
}
/*
Function plot
*/
BEGIN_MESSAGE_MAP(CFunctionDrawing , CStatic)
    ON_WM_MOUSEMOVE()
	ON_WM_PAINT()
END_MESSAGE_MAP()

CFunctionDrawing::CFunctionDrawing(int ttype) {
	section = 0;
	ScaleX = 1;
	type = ttype;
}
void CFunctionDrawing::OnPaint() {
	Plot(TRUE);
	CStatic::OnPaint();
}
void CFunctionDrawing::DoDataExchange(CDataExchange* pDX) {
	CStatic::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_FUNCTION_VALUE,c_FunctionValue);
	
	if(type == EDITABLE) {
		DDX_Control(pDX,IDC_LY,c_Ly);
		DDX_Control(pDX,IDC_LX,c_Lx);
		DDX_Control(pDX,IDC_EY,c_Ey);
		DDX_Control(pDX,IDC_EX,c_Ex);
		DDX_Control(pDX,IDC_ADD,c_Add);
		DDX_Control(pDX,IDC_MOVE,c_Move);
		DDX_Control(pDX,IDC_DELETE,c_Delete);
		DDX_Control(pDX,IDC_DELETEALL,c_DeleteAll);
	}
}
void CFunctionDrawing::Init() {
	if(type != EDITABLE)
		return;

    CString str;
	POSITION pos;
	RPoint v;

	c_Lx.ResetContent();
	c_Ly.ResetContent();

	pos = function->points.GetHeadPosition();
	while(pos) {
		v = function->points.GetNext(pos);
		str.Format("%+.2f",v.x);
		c_Lx.AddString(str);
        str.Format("%+.2f",v.y);
		c_Ly.AddString(str);
	}

	c_Lx.SetCurSel(0);
	c_Ly.SetCurSel(0);

	OnUpdateLoc(IDC_LX);
    OnUpdateLoc(IDC_LY);
	OnUpdateControls();
	Invalidate();
}
void CFunctionDrawing::Plot(BOOL paint) {
	if(function->points.GetCount() < 2)
		return;

    /*paint section drawing*/
	CDC* pDC;
	CPaintDC paintDC(this);
	if(paint) pDC = &paintDC;
	else pDC = GetDC();

	CRect rect;
	GetClientRect(&rect);
	CBrush wbrush(RGB(255,255,255));
	CBrush* oldBrush = paintDC.SelectObject(&wbrush);
	int H = rect.Height(),W = rect.Width();
	pDC->Rectangle(rect);
	pDC->SelectObject(oldBrush);

	CRgn rgn;
	rgn.CreateRectRgn(0,0,W,H);
	pDC->SelectClipRgn(&rgn);

	/*scale*/
	RPoint v;
	DOUBLE xmin = 10000,xmax = -10000,ymin = 10000,ymax = -10000;
	POSITION pos = function->points.GetHeadPosition();
	while(pos) {
        v = function->points.GetNext(pos);
		if(v.x < xmin) xmin = v.x;
        if(v.x > xmax) xmax = v.x;
		if(v.y < ymin) ymin = v.y;
        if(v.y > ymax) ymax = v.y;
	}
	CPoint Center;
	DOUBLE ScaleY;

	ScaleX = (W - 5) / (function->cycles * (xmax - xmin));
	Center.x = int(xmin * ScaleX);

	if(ymax != ymin) {
		ScaleY = (H - 5) / (ymax - ymin);
		Center.y = H / 2 + int((ymin + ymax) * ScaleY / 2);
	} else {
		ScaleY = 0;
		Center.y = 0;
	}
	/*grids*/
	CPen pen(PS_SOLID,1,RGB(200,200,200));
	CPen pen1(PS_SOLID,2,RGB(0,0,255));
	CPen* oldPen = pDC->SelectObject(&pen);
	CPoint p1,p2;
	const int SIZE = 5;
	int x = 0;
	while(x < W) {
		p1.x = x;
		p1.y = 0;
		p2.x = x;
		p2.y = H;
		pDC->MoveTo(p1);
		pDC->LineTo(p2);
        x += SIZE;
	}
	x = 0;
	while(x < H) {
		p1.x = 0;
		p1.y = x;
		p2.x = W;
		p2.y = x;
		pDC->MoveTo(p1);
		pDC->LineTo(p2);
        x += SIZE;
	}
	/*function*/
	pDC->SelectObject(&pen1);
	
	for(int i = 0;i < function->cycles;i++) {
		pos = function->points.GetHeadPosition();
		v = function->points.GetNext(pos);
		p1.x = Center.x + int(ScaleX * (v.x + i * (xmax - xmin)));
		p1.y = Center.y - int(ScaleY * v.y); 
		if(i == 0) pDC->MoveTo(p1);
		else pDC->LineTo(p1);
		
		while(pos) {
			v = function->points.GetNext(pos);
			p1.x = Center.x + int(ScaleX * (v.x + i * (xmax - xmin)));
			p1.y = Center.y - int(ScaleY * v.y); 
			pDC->LineTo(p1);
		}
	}
	pDC->SelectObject(oldPen);

	/*section*/
	DOUBLE y = function->Get(section,TRANSIENT,1);
	CPoint p;
	p.x = Center.x + int(section * ScaleX);
	p.y = Center.y - int(y * ScaleY);

	CBrush brush(RGB(255,0,0));
	CBrush* pOldBrush = pDC->SelectObject(&brush);
	CRect rect1(p.x - 4,p.y - 4,p.x + 4,p.y + 4);
	pDC->Ellipse(&rect1);
	pDC->SelectObject(pOldBrush);

	/*releae dc*/
	ReleaseDC(pDC);
	rgn.DeleteObject();
}
void CFunctionDrawing::OnMouseMove(UINT nFlags,CPoint point) {
	CStatic::OnMouseMove(nFlags,point);
	if(function->points.GetCount() < 2)
		return;
	section = point.x / ScaleX;
	DOUBLE y = function->Get(section,TRANSIENT,1,TRUE);
	CString str;
	str.Format("%.6g at %.3f",y,section);
	c_FunctionValue.SetWindowText(str);
	Plot(FALSE);
}
void CFunctionDrawing::OnUpdateControls() {
	BOOL acyclic = ((function->type == F_USER || function->type == F_PERIODIC) && !function->rank);

	c_Ex.EnableWindow(acyclic);
    c_Ey.EnableWindow(acyclic);
	c_Add.EnableWindow(acyclic);
    c_Move.EnableWindow(acyclic);
	c_Delete.EnableWindow(acyclic);
    c_DeleteAll.EnableWindow(acyclic);
}
void CFunctionDrawing::OnUpdateLoc(UINT nID) {
	CString text;
	int index;
	if(nID == IDC_LX) index = c_Lx.GetCurSel();
	else index = c_Ly.GetCurSel();
	if(index != LB_ERR) {
		if(nID == IDC_LX) c_Ly.SetCurSel(index);
		else c_Lx.SetCurSel(index);
		c_Lx.GetText(index,text);
		c_Ex.SetWindowText(text);
		c_Ly.GetText(index,text);
		c_Ey.SetWindowText(text);
	}
}
void CFunctionDrawing::OnAdd() {
	RPoint v,v1;
	CString str;
	c_Ex.GetWindowText(str);
	v.x = atof(str);
	c_Ey.GetWindowText(str);
	v.y = atof(str);
    function->points.AddSorted(v);
	/*add to list*/
	str.Format("%+.2f",v.x);
	c_Lx.AddString(str);
    str.Format("%+.2f",v.y);
	c_Ly.AddString(str);
	Invalidate();
}
void CFunctionDrawing::OnMove() {
	
	CString str;
	POSITION pos;
	RPoint v;
	DOUBLE d;
	int index = c_Lx.GetCurSel();
	if(index == LB_ERR)
		return;
	c_Lx.GetText(index,str);
	v.x = atof(str);
	c_Ly.GetText(index,str);
	v.y = atof(str);
	if(pos = function->points.Find(v)) {
		c_Ex.GetWindowText(str);
		d = atof(str);
		if(d != v.x)
			return;
		v.x = d;
		c_Ey.GetWindowText(str);
		v.y = atof(str);
		function->points.GetAt(pos) = v;

		str.Format("%+.2f",v.x);
		c_Lx.DeleteString(index);
		c_Lx.InsertString(index,str);
		str.Format("%+.2f",v.y);
		c_Ly.DeleteString(index);
		c_Ly.InsertString(index,str);

		c_Lx.SetCurSel(index);
        c_Ly.SetCurSel(index);
	}
	Invalidate();
}
void CFunctionDrawing::OnDelete() {
	CString str;
	POSITION pos;
	RPoint v;
	int index = c_Lx.GetCurSel();
	if(index == LB_ERR)
		return;
	c_Lx.GetText(index,str);
	v.x = atof(str);
	c_Ly.GetText(index,str);
	v.y = atof(str);
	if(pos = function->points.Find(v)) {
		function->points.RemoveAt(pos);
		c_Lx.DeleteString(index);
		c_Ly.DeleteString(index);
		c_Lx.SetCurSel(index);
		c_Ly.SetCurSel(index);
	}
	Invalidate();
}
void CFunctionDrawing::OnDeleteAll() {
	c_Lx.ResetContent();
    c_Ly.ResetContent();
	function->points.RemoveAll();
	c_Lx.SetCurSel(0);
    c_Ly.SetCurSel(0);
	Invalidate();
}
void CFunctionDrawing::AcceptMessage(UINT nID) {
	if(nID == IDC_LX) OnUpdateLoc(nID);
    else if(nID == IDC_LY) OnUpdateLoc(nID);
    else if(nID == IDC_ADD) OnAdd();
	else if(nID == IDC_MOVE) OnMove();
	else if(nID == IDC_DELETE) OnDelete();
	else if(nID == IDC_DELETEALL) OnDeleteAll();
}
/*
Function Property Dialog
*/
BEGIN_MESSAGE_MAP(CFunctionPropertyDia , CDialog)
    ON_WM_PAINT()
	ON_LBN_SELCHANGE(IDC_TYPE,OnType)
	ON_EN_CHANGE(IDC_PERIOD,OnType)
	ON_EN_CHANGE(IDC_AMPLITUDE,OnType)
	ON_EN_CHANGE(IDC_CYCLES,OnType)
	ON_EN_CHANGE(IDC_DIVISIONS,OnType)
	ON_EN_CHANGE(IDC_AG,OnType)
	ON_CBN_SELCHANGE(IDC_SOILTYPE,OnType)
	ON_COMMAND(IDC_BROWSE,OnBrowse)
	ON_CONTROL_RANGE(LBN_SELCHANGE,IDC_LX,IDC_LY,RelayMessage)
	ON_COMMAND_RANGE(IDC_ADD,IDC_DELETEALL,RelayMessage)
END_MESSAGE_MAP()

CFunctionPropertyDia::CFunctionPropertyDia(CWnd* parent,FUNCTION* pfunction) : 
              CDialog(IDD_FUNCTION,parent),
			  c_Drawing(CFunctionDrawing::EDITABLE) {
	function = *pfunction;
	pfunc = pfunction;
}
void CFunctionPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_NAME,function.name);
	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_PERIOD,c_Period);
	DDX_Control(pDX,IDC_AMPLITUDE,c_Amplitude);
	DDX_Control(pDX,IDC_CYCLES,c_Cycles);
	DDX_Control(pDX,IDC_DIVISIONS,c_Divisions);
    DDX_Control(pDX,IDC_AG,c_Ag);
	DDX_Control(pDX,IDC_SOILTYPE,c_Soil);
	DDX_Text(pDX,IDC_PERIOD,function.period);
	DDX_Text(pDX,IDC_AMPLITUDE,function.amplitude);
	DDX_Text(pDX,IDC_CYCLES,function.cycles);
	DDV_MinMaxInt(pDX,function.cycles,1,10000);
	DDX_Text(pDX,IDC_DIVISIONS,function.divisions);
	DDV_MinMaxInt(pDX,function.divisions,1,10000);
	DDX_Text(pDX,IDC_AG,function.Ag);
	DDV_MinMaxDouble(pDX,function.Ag,0,1000);
    DDX_Control(pDX,IDC_COVER1,c_Cover1);
	DDX_Control(pDX,IDC_COVER2,c_Cover2);
	DDX_Control(pDX,IDC_FILE,c_File);
    DDX_Control(pDX,IDC_BROWSE,c_Browse);
	DDX_Control(pDX,IDC_DRAWING,c_Drawing);
	c_Drawing.DoDataExchange(pDX);
}

BOOL CFunctionPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	
	if(!FUNCTION::defspectrumtype) {
		c_Type.AddString("User function");
		c_Type.AddString("Sine function");
		c_Type.AddString("Cosine function");
		c_Type.AddString("Triangular function");
		c_Type.AddString("Periodic function");
		c_Type.SetCurSel(function.type - F_USER);
		c_Cover1.ShowWindow(FALSE);
		c_Cover2.ShowWindow(TRUE);
	} else {
		c_Type.AddString("EBCS design spectrum");
		c_Type.AddString("EBCS elastric response spectrum");
		c_Type.AddString("User function");
		c_Type.SetCurSel(function.type);

		c_Soil.AddString("A");
        c_Soil.AddString("B");
		c_Soil.AddString("C");
		c_Soil.SetCurSel(function.soil_type);
		c_Cover1.ShowWindow(TRUE);
		c_Cover2.ShowWindow(FALSE);
	}

	c_Drawing.function = &function;
	c_Drawing.Init();

	OnUpdateControls();
	return TRUE;
}
void CFunctionPropertyDia::RelayMessage(UINT nID) {
	c_Drawing.AcceptMessage(nID);
}
void CFunctionPropertyDia::OnUpdateControls() {
	BOOL cyclic = ((function.type != F_USER) && !function.rank);
	if(FUNCTION::defspectrumtype) cyclic = FALSE;

	c_Name.EnableWindow(!function.rank);
	c_Type.EnableWindow(!function.rank);
	c_Amplitude.EnableWindow(cyclic && (function.type != F_PERIODIC));
	c_Period.EnableWindow(cyclic && (function.type != F_PERIODIC));
    c_Cycles.EnableWindow(cyclic);
	c_Divisions.EnableWindow(cyclic && (function.type != F_TRIANGULAR && function.type != F_PERIODIC));
	c_Ag.EnableWindow(function.type == F_EBCS_DESIGN || function.type == F_EBCS_ELASTIC);
	c_Soil.EnableWindow(function.type == F_EBCS_DESIGN || function.type == F_EBCS_ELASTIC);
	c_File.ShowWindow(function.type == F_USER && !function.rank);
    c_Browse.ShowWindow(function.type == F_USER && !function.rank);
	c_Drawing.OnUpdateControls();
}
void CFunctionPropertyDia::OnType() {

	CString str;
	function.type = c_Type.GetCurSel();
	if(!FUNCTION::defspectrumtype) 
		function.type += F_USER;

	c_Period.GetWindowText(str);
	function.period = atof(str);
	c_Amplitude.GetWindowText(str);
	function.amplitude = atof(str);
	c_Cycles.GetWindowText(str);
	function.cycles = atoi(str);
	c_Divisions.GetWindowText(str);
	function.divisions = atoi(str);
	c_Ag.GetWindowText(str);
	function.Ag = atof(str);
	function.soil_type = c_Soil.GetCurSel();

    OnUpdateControls();
	if(function.type == F_USER) {
		function.cycles = 1;
	} else if(function.type == F_PERIODIC) {
	} else if(!FUNCTION::defspectrumtype) {

		function.points.RemoveAll();

		RPoint v;
		CString str;
		int div = (function.type == F_TRIANGULAR) ? 4 : function.divisions;
		
		for(int i = 0;i < div + 1;i++) {
			v.x = i / DOUBLE(div);
			if(function.type == F_SINE) v.y = function.amplitude * sin((2 * PI * v.x) / function.period);
			else if(function.type == F_COSINE) v.y = function.amplitude * cos((2 * PI * v.x) / function.period);
			else if(function.type == F_TRIANGULAR) {
				if(i == 0 || i == 2 || i == 4) v.y = 0;
				else if(i == 1) v.y = function.amplitude;
				else v.y = -function.amplitude;
			}
			function.points.AddTail(v);
		}

		c_Drawing.Init();
	} else {
		function.cycles = 1;
		function.points.RemoveAll();
		

		DOUBLE b,s;
		RPoint v;
		CString str;

	
		switch(function.soil_type) {
		case SOILA: s = 1.0; break;
		case SOILB: s = 1.2; break;
		case SOILC: s = 1.5; break;
		}

		v.x = 0;
		while(v.x <= 10) {
			if(v.x == 0) {
				v.y = function.Ag;
			} else {
				if(function.type == F_EBCS_ELASTIC) {
					switch(function.soil_type) {
					case SOILA: b = 1 / v.x;; break;
					case SOILB: b = 1.5 / v.x;; break;
					case SOILC: b = 2.25 / v.x; break;
					}
                } else if(function.type == F_EBCS_DESIGN) {
					b = 1.2 * s / pow(v.x, 2 / 3.0);
				}
				b = min(b, 2.5);
				v.y = function.Ag * b;
			}

			function.points.AddTail(v);
		
			if(v.x < 0.2) {
				switch(function.soil_type) {
				case SOILA: v.x += 0.1; break;
				case SOILB: v.x += 0.15; break;
				case SOILC: v.x += 0.2; break;
				}
			} else if(v.x < 1) v.x += 0.05;
			else if(v.x < 3) v.x += 0.5;
			else v.x += 1;
		}

		c_Drawing.Init();
	}
	Invalidate();
}

void CFunctionPropertyDia::OnBrowse() {
	char strFilter[] = { "Ground acceleration files (*.gra)|*.gra|Response spectrum files (*.rsp)|*.rsp|Text Files(*.txt)|*.txt|All Files (*.*)|*.*||" };
	CFileDialog FileDlg(TRUE,".gra",NULL,0,strFilter);
	
	if(FileDlg.DoModal() == IDOK) {
		CString name = FileDlg.GetPathName();
		c_File.SetWindowText(name);
		FILE* file;
		if(!(file = fopen(name.GetBuffer(0),"r")))
			return;

		function.points.RemoveAll();

		int status;
		do {
			RPoint v;
			float a,b;
			status = fscanf(file,"%f %f",&a,&b);
			v.x = a;
			v.y = b;
			function.points.AddSorted(v);
		} while(status != EOF);

		c_Drawing.Init();

		fclose(file);
		Invalidate();
	}
}
void CFunctionPropertyDia::OnOK() {
	CDialog::OnOK();
	*pfunc = function;
}
void CFunctionPropertyDia::OnPaint() {
	CDialog::OnPaint();
	c_Drawing.Plot(TRUE);
}
/*
Function plot
*/
BEGIN_MESSAGE_MAP(CFunctionPlotDia , CDialog)
    ON_WM_PAINT()
END_MESSAGE_MAP()

CFunctionPlotDia::CFunctionPlotDia(CWnd* parent,FUNCTION* pfunction) : 
              CDialog(IDD_FUNCTION_PLOT,parent) ,
			  c_Drawing(CFunctionDrawing::VIEWONLY) {
	function = pfunction;
}
void CFunctionPlotDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_DRAWING,c_Drawing);
	c_Drawing.DoDataExchange(pDX);
}

BOOL CFunctionPlotDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_Drawing.function = function;
	c_Drawing.Init();

	return TRUE;
}
void CFunctionPlotDia::OnPaint() {
	CDialog::OnPaint();
	c_Drawing.Plot(TRUE);
}
void CFunctionPlotDia::OnOK() {
	CDialog::OnOK();
}
/*
Result plot dialog
*/
BEGIN_MESSAGE_MAP(CResultPlotDia , CDialog)
    ON_COMMAND(IDC_DISPLAY,OnDisplay)
	ON_CBN_SELCHANGE(IDC_TYPE,OnType)
END_MESSAGE_MAP()

CResultPlotDia::CResultPlotDia(CWnd* parent) : 
              CDialog(IDD_PLOT_RESULT,parent) {
}
void CResultPlotDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_JLIST,c_JointL);
    DDX_Control(pDX,IDC_LISTDIR,c_Dir);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_SECTION,c_Section);
	DDX_Control(pDX,IDC_COVER1,c_Cover1);
}
BOOL CResultPlotDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	for(int i = IUX;i <= IRZ;i++)
		c_Dir.AddString(DirString[i]);
	c_Dir.SetCurSel(0);

	c_Type.AddString("Joint Displacement");
    c_Type.AddString("Joint Forces");
    c_Type.AddString("Member Forces");
	c_Type.SetCurSel(0);
	OnType();

	c_Section.SetWindowText("0.00");

    OnUpdateControls();
	return TRUE;
}
void CResultPlotDia::OnUpdateControls() {
	int type = c_Type.GetCurSel();
    c_Cover1.ShowWindow(type != MEMBER_FORCE);
}
void CResultPlotDia::OnType() {
	CmyDocument* pDoc = GetMyDocument();
    int type = c_Type.GetCurSel();
	JOINT* joint;
	MEMBER* member;
	POSITION pos;

	c_JointL.ResetContent();

	if(type <= JOINT_FORCE) {
		pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			c_JointL.AddString(joint->name);
		}
		c_JointL.SetCurSel(0);
	} else {
		pos = pDoc->members.GetHeadPosition();
		while(pos) {
			member = &pDoc->members.GetNext(pos);
			c_JointL.AddString(member->name);
		}
		c_JointL.SetCurSel(0);
	}

	OnUpdateControls();
}
void CResultPlotDia::OnDisplay() {
	CmyDocument* pDoc = GetMyDocument();
	int type = c_Type.GetCurSel();
	int dir = c_Dir.GetCurSel();
	CString ename;
	c_JointL.GetText(c_JointL.GetCurSel(),ename);
	JOINT* joint;
	MEMBER* member;
	RESPONSEHIST* pcase;
	pcase = (RESPONSEHIST*) pDoc->c_AnalysisCase;
	FUNCTION function;
	RPoint rp;
	POSITION pos;
	UINT i;

	if(type <= JOINT_FORCE) {
		pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(joint->name == ename) break;
		}
		for(i = 0;i < pcase->N;i++) {
			rp.x = i * pcase->dt;
			if(type == JOINT_DISPLACEMENT) rp.y = joint->disps_all[pcase->index + i][dir];
			else if(type == JOINT_FORCE) rp.y = joint->forces_all[pcase->index + i][dir];
			function.points.AddTail(rp);
		}
	} else {
		DOUBLE x;
		CString str;
		c_Section.GetWindowText(str);
		x = atof(str);
		
		pos = pDoc->members.GetHeadPosition();
		while(pos) {
			member = &pDoc->members.GetNext(pos);
			if(member->name == ename) break;
		}
		for(i = 0;i < pcase->N;i++) {
			rp.x = i * pcase->dt;
			if(type == MEMBER_FORCE) rp.y = member->CalculateForce(x,dir,pcase->index + i);
			function.points.AddTail(rp);
		}
	}

	CFunctionPlotDia Dia(this,&function);
	Dia.DoModal();
}
/*
Loadcase Property Dialog
*/
CLoadCasePropertyDia::CLoadCasePropertyDia(CWnd* parent,LOADCASE* ploadcase) : 
                    CDialog(IDD_LOADCASE,parent) {
    loadcase = ploadcase;
}
BOOL CLoadCasePropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_List.AddString("DEAD");
    c_List.AddString("LIVE");
	c_List.AddString("WIND");
	c_List.AddString("QUAKE");
	c_List.SetCurSel(loadcase->type);
	if(loadcase->rank) {
		c_Name.EnableWindow(FALSE);
	}
	return TRUE;
}
void CLoadCasePropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Text(pDX, IDC_NAME, loadcase->name);
    DDX_Text(pDX, IDC_SWM, loadcase->swm);
	DDV_MinMaxDouble(pDX,loadcase->swm,0,10000);
	DDX_Control(pDX,IDC_TYPE,c_List);
}
void CLoadCasePropertyDia::OnOK() {
	CDialog::OnOK();
	loadcase->type = LOADCASE::DEAD + c_List.GetCurSel();
}
/*
Modalcase Property Dialog
*/
CModalCasePropertyDia::CModalCasePropertyDia(CWnd* parent,MODALCASE* pmodalcase) : 
                    CDialog(IDD_MODALCASE,parent) {
    modalcase = pmodalcase;
}
BOOL CModalCasePropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	if(modalcase->rank) {
		c_Name.EnableWindow(FALSE);
	}
	return TRUE;
}
void CModalCasePropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Text(pDX, IDC_NAME, modalcase->name);
    DDX_Text(pDX, IDC_MIN, modalcase->minm);
	DDV_MinMaxInt(pDX,modalcase->minm,1,10000);
	DDX_Text(pDX, IDC_MAX, modalcase->maxm);
	DDV_MinMaxInt(pDX,modalcase->maxm,modalcase->minm,10000);
    DDX_Text(pDX, IDC_SHIFT, modalcase->shift);
	DDX_Text(pDX, IDC_TOLERANCE, modalcase->tolerance);
}
void CModalCasePropertyDia::OnOK() {
	CDialog::OnOK();
}
/*
Damping Dialog
*/
BEGIN_MESSAGE_MAP(CDampingDia , CDialog)
    ON_COMMAND(IDC_CALCULATE,OnType)
	ON_CBN_SELCHANGE(IDC_TYPE,OnType)
	ON_CONTROL_RANGE(LBN_SELCHANGE,IDC_LX,IDC_LY,RelayMessage)
	ON_COMMAND_RANGE(IDC_ADD,IDC_DELETEALL,RelayMessage)
END_MESSAGE_MAP()

CDampingDia::CDampingDia(CWnd* parent,DAMPING* pdamping) : 
                    CDialog(IDD_DAMPING,parent),
					c_Drawing(CFunctionDrawing::EDITABLE) {
	damping = *pdamping;
}
BOOL CDampingDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

    c_Type.AddString("Rayleigh Damping");
    c_Type.AddString("Constant Damping");
	c_Type.SetCurSel(damping.type);

	c_Drawing.function = &damping.dvalues;
	c_Drawing.Init();

	OnUpdateControls();
	return TRUE;
}
void CDampingDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_FM, damping.fm);
	DDX_Text(pDX, IDC_FK, damping.fk);
	DDX_Text(pDX, IDC_W1, damping.w1);
	DDX_Text(pDX, IDC_W2, damping.w2);
	DDX_Text(pDX, IDC_E1, damping.e1);
	DDX_Text(pDX, IDC_E2, damping.e2);
	DDX_Text(pDX, IDC_CONSTANT,damping.cdamping);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_COVER1,c_Cover1);
	DDX_Control(pDX,IDC_COVER2,c_Cover2);
	DDX_Control(pDX,IDC_CONSTANT,c_Constant);
	DDX_Control(pDX,IDC_W2,c_W2);
	DDX_Control(pDX,IDC_DRAWING,c_Drawing);
	c_Drawing.DoDataExchange(pDX);
}
void CDampingDia::RelayMessage(UINT nID) {
	c_Drawing.AcceptMessage(nID);
}
void CDampingDia::OnUpdateControls() {
	int type = c_Type.GetCurSel();
	c_Cover1.ShowWindow(type != DAMPING::CONSTANT);
    c_Cover2.ShowWindow(type != DAMPING::RAYLEIGH);
}
void CDampingDia::OnType() {
	int type = c_Type.GetCurSel();
	OnUpdateControls();
	UpdateData(TRUE);
	damping.type = type;
	damping.Calculate();
	UpdateData(FALSE);
	c_Drawing.Init();
}
/*
Response Property Dialog
*/
BEGIN_MESSAGE_MAP(CResponsePropertyDia , CDialog)
    ON_LBN_SELCHANGE(IDC_LISTDIR,OnUpdateLocDir)
	ON_LBN_SELCHANGE(IDC_LISTFUNCTION,OnUpdateLocFunc)
    ON_LBN_SELCHANGE(IDC_LISTSCALE,OnUpdateLocScale)
	ON_COMMAND(IDC_ADD,OnAdd)
    ON_COMMAND(IDC_MOVE,OnMove)
	ON_COMMAND(IDC_DELETE,OnDelete)
	ON_COMMAND(IDC_DELETEALL,OnDeleteAll)
	ON_COMMAND(IDC_DAMPING,OnDamping)
END_MESSAGE_MAP()

CResponsePropertyDia::CResponsePropertyDia(UINT mID, CWnd* parent,
		            RESPONSE* presponse,PMODALCASELIST pmodalcases,PFUNCTIONLIST pfunctions) : 
                    CDialog(mID,parent) {
    response = presponse;
    modalcases = pmodalcases;
	functions = pfunctions;
	rID = mID;
	damping = response->damping;
	if(damping.type == DAMPING::NONE) {
		if(response->acase_type == ANALYSISCASE::RESPONSEH_CASE) damping.type = DAMPING::RAYLEIGH;
		else damping.type = DAMPING::CONSTANT;
        damping.Calculate();
	}
}
BOOL CResponsePropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	if(response->rank) {
		c_Name.EnableWindow(FALSE);
        c_ModalList.EnableWindow(FALSE); 
	}

	POSITION pos;

    pos = modalcases->GetHeadPosition();
	while(pos) c_ModalList.AddString(modalcases->GetNext(pos).name);
	c_ModalList.SetCurSel(c_ModalList.FindString(0,response->modalcase->name));

    pos = functions->GetHeadPosition();
	while(pos) c_Func.AddString(functions->GetNext(pos).name);
	c_Func.SetCurSel(0);

	for(int i = IUX;i <= IRZ;i++)
		c_Dir.AddString(DirString[i]);
	c_Dir.SetCurSel(0);

	c_Scale.SetWindowText("1.0");

	CString str;
	SPECFUNC* pfunc;
	pos = response->funclist.GetHeadPosition();
	while(pos) {
        pfunc = &response->funclist.GetNext(pos);
		c_DirList.AddString(DirString[pfunc->dir]);
		c_FuncList.AddString(pfunc->function->name);
		str.Format("%.2f",pfunc->scale);
		c_ScaleList.AddString(str);
	}

	OnUpdateLocDir();
	return TRUE;
}
void CResponsePropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_NAME,c_Name);
    DDX_Control(pDX,IDC_MODALCASE,c_ModalList);
	DDX_Control(pDX,IDC_FUNCTION,c_Func);
	DDX_Control(pDX,IDC_DIR,c_Dir);
    DDX_Control(pDX,IDC_SCALE,c_Scale);
	DDX_Control(pDX,IDC_LISTDIR,c_DirList);
	DDX_Control(pDX,IDC_LISTFUNCTION,c_FuncList);
	DDX_Control(pDX,IDC_LISTSCALE,c_ScaleList);
	DDX_Text(pDX, IDC_NAME, response->name);
	
	if(rID == IDD_RESPONSECASE) {
		RESPONSEHIST* resp = (RESPONSEHIST*) response;
		DDX_Text(pDX, IDC_STEP, resp->dt);
		DDX_Text(pDX, IDC_N, resp->N);
		DDX_Radio(pDX, IDC_PERIODIC,resp->type);
		DDX_Radio(pDX, IDC_MODAL,resp->ana_type);
	} else {
		RESPONSESPEC* resp = (RESPONSESPEC*) response;
		DDX_Radio(pDX, IDC_ABSSUM, resp->modal_comb);
		DDX_Radio(pDX, IDC_DIRABSSUM, resp->dir_comb);
	}
}
void CResponsePropertyDia::OnDamping() {
	CDampingDia Dia(this,&damping);
	if(Dia.DoModal() == IDOK) {
		damping = Dia.damping;
	}
}
void CResponsePropertyDia::OnOK() {
	POSITION pos;
	MODALCASE md;
	c_ModalList.GetLBText(c_ModalList.GetCurSel(),md.name);
	if(pos = modalcases->Find(md)) {
		response->modalcase = &modalcases->GetAt(pos);
		response->damping = damping;
	}
	CDialog::OnOK();
}

void CResponsePropertyDia::OnAdd() {
	SPECFUNC specfunc;
	CString str;
	c_Scale.GetWindowText(str);
	specfunc.scale = atof(str);

	specfunc.dir = c_Dir.GetCurSel();

	FUNCTION fn;
	c_Func.GetLBText(c_Func.GetCurSel(), fn.name);
	specfunc.function = &functions->GetAt(functions->Find(fn));

	if(response->funclist.Find(specfunc))
		return;

	response->funclist.AddTail(specfunc);
	c_DirList.SetCurSel(c_DirList.AddString(DirString[specfunc.dir]));
	c_FuncList.SetCurSel(c_FuncList.AddString(specfunc.function->name));
	str.Format("%.2f",specfunc.scale);
	c_ScaleList.SetCurSel(c_ScaleList.AddString(str));
}
void CResponsePropertyDia::OnMove() {
	CString str;
	POSITION pos;
	SPECFUNC fn;
	FUNCTION f;
	int index = c_DirList.GetCurSel();
	if(index == LB_ERR)
		return;

	c_ScaleList.GetText(index,str);
	fn.scale = atof(str);
	c_FuncList.GetText(index,f.name);
	fn.function = &functions->GetAt(functions->Find(f));
	c_DirList.GetText(index,str);
	fn.dir = c_Dir.FindString(0,str);

	if(pos = response->funclist.Find(fn)) {
		fn.dir = c_Dir.GetCurSel();
		c_Scale.GetWindowText(str);
		fn.scale = atof(str);
	    c_Func.GetLBText(c_Func.GetCurSel(), f.name);
	    fn.function = &functions->GetAt(functions->Find(f));
		response->funclist.GetAt(pos) = fn;

		c_DirList.DeleteString(index);
		c_DirList.InsertString(index,DirString[fn.dir]);
		c_FuncList.DeleteString(index);
		c_FuncList.InsertString(index,fn.function->name);
		str.Format("%.2f",fn.scale);
		c_ScaleList.DeleteString(index);
		c_ScaleList.InsertString(index,str);
		c_DirList.SetCurSel(index);
		c_FuncList.SetCurSel(index);
		c_ScaleList.SetCurSel(index);
	}
}
void CResponsePropertyDia::OnDelete() {
    CString str;
	POSITION pos;
	int index = c_DirList.GetCurSel();
	if(index == LB_ERR)
		return;
	SPECFUNC fn;
	c_ScaleList.GetText(index,str);
	fn.scale = atof(str);
	
	FUNCTION f;
	c_FuncList.GetText(index,f.name);
	fn.function = &functions->GetAt(functions->Find(f));

	c_DirList.GetText(index,str);
	fn.dir = c_Dir.FindString(0,str);

	if(pos = response->funclist.Find(fn)) {
		response->funclist.RemoveAt(pos);
		c_DirList.DeleteString(index);
		c_FuncList.DeleteString(index);
		c_ScaleList.DeleteString(index);
		c_DirList.SetCurSel(index);
		c_FuncList.SetCurSel(index);
		c_ScaleList.SetCurSel(index);
	}
}
void CResponsePropertyDia::OnDeleteAll() {
	response->funclist.RemoveAll();
	c_DirList.ResetContent();
    c_FuncList.ResetContent();
	c_ScaleList.ResetContent();
}
void CResponsePropertyDia::OnUpdateLocDir() {
	OnUpdateLoc(IDC_LISTDIR);
}
void CResponsePropertyDia::OnUpdateLocFunc() {
	OnUpdateLoc(IDC_LISTFUNCTION);
}
void CResponsePropertyDia::OnUpdateLocScale() {
	OnUpdateLoc(IDC_LISTSCALE);
}
void CResponsePropertyDia::OnUpdateLoc(UINT nID) {
	CString text;
	int index;
	if(nID == IDC_LISTDIR) index = c_DirList.GetCurSel();
	else if(nID == IDC_LISTFUNCTION) index = c_FuncList.GetCurSel();
	else index = c_ScaleList.GetCurSel();
	if(index != LB_ERR) {
		if(nID == IDC_LISTDIR) {
			c_FuncList.SetCurSel(index);
            c_ScaleList.SetCurSel(index);
		} else if(nID == IDC_LISTFUNCTION) {
			c_DirList.SetCurSel(index);
            c_ScaleList.SetCurSel(index);
		} else {
			c_DirList.SetCurSel(index);
            c_FuncList.SetCurSel(index);
		}
		c_DirList.GetText(index,text);
		c_Dir.SetCurSel(c_Dir.FindString(0,text));
		c_FuncList.GetText(index,text);
		c_Func.SetCurSel(c_Func.FindString(0,text));
		c_ScaleList.GetText(index,text);
		c_Scale.SetWindowText(text);
	}
}
/*
Constraint Property Dialog
*/
BEGIN_MESSAGE_MAP(CConstraintPropertyDia , CDialog)
	ON_LBN_SELCHANGE(IDC_TYPE,OnType)
END_MESSAGE_MAP()

CConstraintPropertyDia::CConstraintPropertyDia(CWnd* parent,CONSTRAINT* pconstraint,PSYSTEMLIST psystems) : 
                    CDialog(IDD_CONSTRAINT,parent) {

	pconst = pconstraint;
    constraint = *pconstraint;
	systems = psystems;

	UINT i,cons = constraint.constraint;
	BOOL plane = (
		constraint.type == CONSTRAINT::DIAPHRAGM ||
		constraint.type == CONSTRAINT::PLATE ||
		constraint.type == CONSTRAINT::BEAM ||
		constraint.type == CONSTRAINT::ROD);

	if(plane) {
        if(constraint.isauto) {
			dof[0] = IRX;
		} else {
			for(i = IUX;i <= IUZ;i++) {
				if(cons & (1 << i)) {
					dof[0] = i;
					break;
				}
			}
		}
	} else {
		for(i = IUX;i <= IRZ;i++) {
			if(cons & (1 << i)) dof[i] = TRUE;
			else dof[i] = FALSE;
		}
	}
}
BOOL CConstraintPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos;
	JOINT* joint;

    pos = constraint.jointlist.GetHeadPosition();
	while(pos) {
		joint = constraint.jointlist.GetNext(pos);
		c_List.AddString(joint->name);
	}
	c_List.SetCurSel(0);

	pos = systems->GetHeadPosition();
	while(pos)
		c_System.AddString(systems->GetNext(pos).name);
	if(constraint.axis)
		c_System.SetCurSel(c_System.FindString(-1,constraint.axis->name));
    else
		c_System.SetCurSel(0);

	c_Type.AddString("Body");
	c_Type.AddString("Diaphragm");
	c_Type.AddString("Equal");
	c_Type.AddString("Plate");
	c_Type.AddString("Beam");
	c_Type.AddString("Rod");
	c_Type.AddString("Line");
	c_Type.AddString("Local");
	c_Type.AddString("Weld");
	c_Type.SetCurSel(constraint.type);

	OnUpdateControls();
	return TRUE;
}
void CConstraintPropertyDia::OnUpdateControls() {
	BOOL plane = (
		constraint.type == CONSTRAINT::DIAPHRAGM ||
		constraint.type == CONSTRAINT::PLATE ||
		constraint.type == CONSTRAINT::BEAM ||
		constraint.type == CONSTRAINT::ROD);
	for(int i = IUX;i <= IRZ;i++) {
		c_Check[i].ShowWindow(!plane);
		if(i <= IRX)
			c_Radio[i].ShowWindow(plane);
	}
	c_System.ShowWindow(constraint.type != CONSTRAINT::LOCAL);
	c_Cover.ShowWindow(constraint.type != CONSTRAINT::WELD);
}
void CConstraintPropertyDia::OnType() {
	constraint.type = c_Type.GetCurSel();
	BOOL plane = (
		constraint.type == CONSTRAINT::DIAPHRAGM ||
		constraint.type == CONSTRAINT::PLATE ||
		constraint.type == CONSTRAINT::BEAM ||
		constraint.type == CONSTRAINT::ROD);

	
	if(plane) {
		dof[0] = IRX;
	} else {
		for(int i = IUX;i <= IRZ;i++) 
			dof[i] = FALSE;
	}
	UpdateData(FALSE);
	OnUpdateControls();
}
void CConstraintPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	register int i;
	BOOL plane = (
		constraint.type == CONSTRAINT::DIAPHRAGM ||
		constraint.type == CONSTRAINT::PLATE ||
		constraint.type == CONSTRAINT::BEAM ||
		constraint.type == CONSTRAINT::ROD);
	

	DDX_Text(pDX, IDC_NAME, constraint.name);
	DDX_Text(pDX, IDC_TOLERANCE, constraint.tolerance);
	DDX_Control(pDX,IDC_COVER,c_Cover);
	DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_SYSTEM,c_System);

	for(i = IUX;i <= IRZ;i++) {
		DDX_Control(pDX,IDC_UX + i,c_Check[i]);
		if(i <= IRX)
			DDX_Control(pDX,IDC_NUX + i,c_Radio[i]);
	}
	if(plane) {
		DDX_Radio(pDX,IDC_NUX,dof[0]);
	} else {
		for(i = IUX;i <= IRZ;i++) {
			DDX_Check(pDX,IDC_UX + i,dof[i]);
		}
	}
}
void CConstraintPropertyDia::OnOK() {
	CDialog::OnOK();

	BOOL plane = (
		constraint.type == CONSTRAINT::DIAPHRAGM ||
		constraint.type == CONSTRAINT::PLATE ||
		constraint.type == CONSTRAINT::BEAM ||
		constraint.type == CONSTRAINT::ROD);

	constraint.isauto = FALSE;
	if(plane) {
		if(dof[0] == IRX) constraint.isauto = TRUE;
		else constraint.constraint = (1 << dof[0]);
	} else {
		constraint.constraint = NOUR;
		for(int i = IUX;i <= IRZ;i++) {
			if(dof[i]) 
				constraint.constraint |= (1 << i);
		}
	} 
	
	int sel = c_System.GetCurSel();
	constraint.axis = &systems->GetAt(systems->FindIndex(sel));

	*pconst = constraint;
}
/*
Loadcombo Property Dialog
*/
BEGIN_MESSAGE_MAP(CLoadComboPropertyDia , CDialog)
    ON_LBN_SELCHANGE(IDC_LIST,OnUpdateLoc)
    ON_COMMAND(IDC_ADD,OnAdd)
    ON_COMMAND(IDC_MOVE,OnMove)
	ON_COMMAND(IDC_DELETE,OnDelete)
	ON_COMMAND(IDC_DELETEALL,OnDeleteAll)
END_MESSAGE_MAP()

CLoadComboPropertyDia::CLoadComboPropertyDia(CWnd* parent,UINT pID,COMB_TYPE* ploadcombo,PANALYSISCASELIST* ploadcase,int ttype) : 
                    CDialog(pID,parent) {
	combo = ploadcombo;
	type = ttype;
	for(int i = 0;i < CASETYPES;i++)
		loadcase[i] = ploadcase[i];
}
BOOL CLoadComboPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	if(type == TCOMBO) {
		SetWindowText("Load Combination");
	} else if(type == TNL) {
		SetWindowText("Non-Linear analysis case");
	} else {
		SetWindowText("Buckling case");
	}

	POSITION pos;
	ANALYSISCASE* pcase;
    for(int i = 0;i < CASETYPES;i++) {
		pos = loadcase[i]->GetHeadPosition();
		while(pos) {
			pcase = &loadcase[i]->GetNext(pos);
			if(type == TCOMBO) {
				if(pcase->name == combo->name) 
					continue;
				if(pcase->acase_type == ANALYSISCASE::MODAL_CASE)
					continue;
			} else {
				if(pcase->acase_type != ANALYSISCASE::LOAD_CASE)
					continue;
			}
			c_CaseList.AddString(pcase->name);
		}
	}
	c_CaseList.SetCurSel(0);

	LOADCOMBO* pload;
	CString str;
    pos = combo->loadlist.GetHeadPosition();
	while(pos) {
		pload = &combo->loadlist.GetNext(pos);
		str.Format("%s",pload->loadcase->name);
		c_List.AddString(str);
	}
	c_List.SetCurSel(0);

    if(type == TCOMBO) {
		c_Type.AddString("Linear Add");
		c_Type.AddString("Envelope");
		c_Type.AddString("Absolute Add");
		c_Type.AddString("SRSS");
	} else if(type == TNL) {
		c_Type.AddString("P-Delta - Linearized");
        c_Type.AddString("P-Delta - Consistent(cubic)");
	} else {
		c_Type.AddString("Linear Buckling");
	}

	c_Type.SetCurSel(combo->type);
	c_FS.SetWindowText("1.00");

	OnUpdateLoc();
	return TRUE;
}
void CLoadComboPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Text(pDX, IDC_NAME, combo->name);
	DDX_CBIndex(pDX,IDC_TYPE,combo->type);
	DDX_Control(pDX,IDC_LOADCASE,c_CaseList);
	DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_FS,c_FS);
	DDX_Control(pDX,IDC_TYPE,c_Type);
}
void CLoadComboPropertyDia::OnAdd() {
	LOADCOMBO ldcmb;
	LOADCASE ld;
	CString str;
	POSITION pos;

	c_FS.GetWindowText(str);
	ldcmb.FS = atof(str);
	c_CaseList.GetLBText(c_CaseList.GetCurSel(), ld.name);
	for(int i = 0;i < CASETYPES;i++) {
		if(pos = loadcase[i]->Find(ld)) {
			ldcmb.loadcase = &loadcase[i]->GetAt(pos);
			break;
		}
	}
	combo->loadlist.AddTail(ldcmb);

	str.Format("%s",ldcmb.loadcase->name);
	c_List.SetCurSel(c_List.AddString(str));
}
void CLoadComboPropertyDia::OnMove() {
	CString str;
	POSITION pos,pos1;
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		LOADCOMBO ldcmb;
		LOADCASE ld;
		ldcmb.loadcase = &ld;
		c_List.GetText(index,ldcmb.loadcase->name);

		if(pos = combo->loadlist.Find(ldcmb)) {
			c_FS.GetWindowText(str);
			ldcmb.FS = atof(str);
			c_CaseList.GetLBText(c_CaseList.GetCurSel(), ld.name);
			for(int i = 0;i < CASETYPES;i++) {
				if(pos1 = loadcase[i]->Find(ld)) {
					ldcmb.loadcase = &loadcase[i]->GetAt(pos1);
					break;
				}
			}
			combo->loadlist.GetAt(pos) = ldcmb;
			
			str.Format("%s",ldcmb.loadcase->name);
			c_List.DeleteString(index);
			c_List.InsertString(index,str);
			c_List.SetCurSel(index);
		}
	}
}
void CLoadComboPropertyDia::OnDelete() {
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		LOADCOMBO ldcmb;
		LOADCASE ldcase;
		ldcmb.loadcase = &ldcase;
		c_List.GetText(index,ldcmb.loadcase->name);
		combo->loadlist.RemoveAt(combo->loadlist.Find(ldcmb));
		c_List.DeleteString(index);
		c_List.SetCurSel(index);
	}
}
void CLoadComboPropertyDia::OnDeleteAll() {
	combo->loadlist.RemoveAll();
	c_List.ResetContent();
}
void CLoadComboPropertyDia::OnUpdateLoc() {
	int index = c_List.GetCurSel();
	if(index != LB_ERR) {
		LOADCOMBO* pldcmb;
		LOADCOMBO ldcmb;
		LOADCASE ldcase;
		ldcmb.loadcase = &ldcase;
		c_List.GetText(index,ldcmb.loadcase->name);
		pldcmb = &combo->loadlist.GetAt(combo->loadlist.Find(ldcmb));
		CString str;
		str.Format("%.2f",pldcmb->FS);
		c_FS.SetWindowText(str);
	}
}
/*
combination
*/
CCombinationPropertyDia::CCombinationPropertyDia(CWnd* parent,COMB_TYPE* ploadcombo,PANALYSISCASELIST* ploadcase) :
                         CLoadComboPropertyDia(parent,IDD_COMBINATION,ploadcombo,ploadcase,0) {
}
/*
non-linear case
*/
CNLcasePropertyDia::CNLcasePropertyDia(CWnd* parent,NLCASE* tpnl_case,PANALYSISCASELIST* ploadcase) :
                         CLoadComboPropertyDia(parent,IDD_NLCASE,tpnl_case,ploadcase,1) {
	pnl_case = tpnl_case;
}
void CNLcasePropertyDia::DoDataExchange(CDataExchange* pDX) {
	CLoadComboPropertyDia::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NSTEPS, pnl_case->nSteps);
	DDV_MinMaxInt(pDX,pnl_case->nSteps,1,100000);
	DDX_Text(pDX, IDC_NITERATION, pnl_case->nIteration);
	DDV_MinMaxInt(pDX,pnl_case->nSteps,1,100000);
	DDX_Text(pDX, IDC_TOLERANCE, pnl_case->tolerance);
    DDX_Check(pDX, IDC_STIFFNESS, pnl_case->save_stiffness);
}
/*
Buckling case
*/
CBucklingPropertyDia::CBucklingPropertyDia(CWnd* parent,BUCKLINGCASE* tpbkl_case,PANALYSISCASELIST* ploadcase) :
                         CLoadComboPropertyDia(parent,IDD_BUCKLING,tpbkl_case,ploadcase,2) {
	pbkl_case = tpbkl_case;
}
void CBucklingPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CLoadComboPropertyDia::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NMODES, pbkl_case->nmodes);
	DDV_MinMaxInt(pDX,pbkl_case->nmodes,1,100000);
	DDX_Text(pDX, IDC_TOLERANCE,pbkl_case->tolerance);
}

/*
Joint Restraint Dialog
*/
BEGIN_MESSAGE_MAP(CJointRestraintDia , CDialog)
	ON_COMMAND_RANGE(IDC_FIXED,IDC_FREE,OnFastRestraint)
END_MESSAGE_MAP()

CJointRestraintDia::CJointRestraintDia(CWnd* parent,JOINT* pjoint) : 
                  CDialog(IDD_JOINT_RESTRAINT,parent) {
	register int i;
	for(i = IUX;i <= IRZ;i++) {
		if(pjoint->restraint & (1 << i)) res[i] = TRUE;
		else res[i] = FALSE;
	}
}
void CJointRestraintDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	register int i;
	for(i = IUX;i <= IRZ;i++)
	    DDX_Check(pDX,IDC_UX + i,res[i]);
}
void CJointRestraintDia::OnFastRestraint(UINT nID) {
	UBMP8 restraint;
	if(nID == IDC_FIXED) restraint = ALLUR;
	else if(nID == IDC_PINNED) restraint = ALLU;
	else if(nID == IDC_ROLLER) restraint = UZ;
	else restraint = NOUR;

	for(int i = IUX;i <= IRZ;i++) {
		if(restraint & (1 << i)) res[i] = TRUE;
		else res[i] = FALSE;
	}
	UpdateData(FALSE);
}
/*
Joint Load Dialog
*/
CJointLoadDia::CJointLoadDia(CWnd* parent,PLOADCASELIST pldcase,PSYSTEMLIST psystems,JOINT* pjoint,int ttype,LOADCASE* ccase) : 
               CDialog(IDD_JOINT_LOAD,parent) {
	type = ttype;
	loadcases = pldcase;
	ploadcase = ccase;
	systems = psystems;
	joint = pjoint;
	for(int i = 0; i < 6;i++)
		Q[i] = 0;
    mode = F_REPLACE;
}

BOOL CJointLoadDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	if(type == MASS) {
		SetWindowText("Mass");
		c_Cover.ShowWindow(TRUE);
		c_List.ShowWindow(FALSE);
	} else {
		if(type == FORCE) SetWindowText("Force");
		else if(type == DISPLACEMENT) SetWindowText("Diplacement");
		c_Cover.ShowWindow(type == MASS);
		
		POSITION pos = loadcases->GetHeadPosition();
		while(pos) {
			c_List.AddString(loadcases->GetNext(pos).name);
		}
		c_List.SetCurSel(c_List.FindString(-1,ploadcase->name));
	}
	
	POSITION pos = systems->GetHeadPosition();
	while(pos) {
		c_System.AddString(systems->GetNext(pos).name);
	}
	c_System.AddString("Local");
	c_System.SetCurSel(0);
	
	return TRUE;
}
void CJointLoadDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_LOADCASE,c_List);
	DDX_Control(pDX,IDC_COVER,c_Cover);
	DDX_Control(pDX,IDC_SYSTEM,c_System);
    DDX_Text(pDX, IDC_UX, Q[IUX]);
 	DDX_Text(pDX, IDC_UY, Q[IUY]);
	DDX_Text(pDX, IDC_UZ, Q[IUZ]);
	DDX_Text(pDX, IDC_RX, Q[IRX]);
	DDX_Text(pDX, IDC_RY, Q[IRY]);
	DDX_Text(pDX, IDC_RZ, Q[IRZ]);
	DDX_Radio(pDX,IDC_ADD_TO_EXISTING,mode);
}
void CJointLoadDia::OnOK() {
	CDialog::OnOK();
    
	if(type != MASS) {
		LOADCASE ld;
		SYSTEM sys;
		POSITION pos;

		c_List.GetLBText(c_List.GetCurSel(),ld.name);
		if(pos = loadcases->Find(ld)) {
			ploadcase = &loadcases->GetAt(pos);
		} else {
			ploadcase = NULL;
		}

		c_System.GetLBText(c_System.GetCurSel(),sys.name);
		if(pos = systems->Find(sys)) {
			system = &systems->GetAt(pos);
		} else {
			system = NULL;
		}
	}
}
/*
Joint Local Axis Dialog
*/
CJointAxisDia::CJointAxisDia(CWnd* parent,JOINT* pjoint) : 
               CDialog(IDD_JLOCAL_AXIS,parent) {
	joint = pjoint;
	rotation = pjoint->rotation * (180 / PI);
}
void CJointAxisDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_RX, rotation.x);
	DDX_Text(pDX, IDC_RY, rotation.y);
	DDX_Text(pDX, IDC_RZ, rotation.z);
}
void CJointAxisDia::OnOK() {
	CDialog::OnOK();
	rotation = rotation * (PI / 180);
}
/*
Joint Reactions Dialogue
*/
CJointReactionsDia::CJointReactionsDia(CWnd* parent,JOINT* pjoint,int ttype) : 
        CDialog(IDD_JOINT_REACTIONS,parent) ,
		c_List(CGridCtrl::TYPE3) {
	joint = pjoint;
	type = ttype;
}
void CJointReactionsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST, c_List);
}
BOOL CJointReactionsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	VECTOR v;
	if(type == IDEFLECTED || type == IDISPLACEMENT) {
		v = joint->disps;
		SetWindowText("Joint Displacements");
	} else {
		v = joint->forces;
		SetWindowText("Joint Reactions");
	}
	CString str;
	c_List.SetTwoHeaders();
	for(int i = IUX;i <= IRZ;i++) {
		c_List.InsertMyItem(GRIDCTRLDATA::NONE,GRIDCTRLDATA::TDBL_F,&v[i],DirString[i]);
	}
	return TRUE;
}
/*
Member Local Axis Dialog
*/
CAxisRotationDia::CAxisRotationDia(CWnd* parent,DOUBLE talpha,CString ttitle) : 
               CDialog(IDD_MLOCAL_AXIS,parent) {
	alpha = (180 * talpha) / PI;
	title = ttitle;
}
void CAxisRotationDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_RX, alpha);
}
BOOL CAxisRotationDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	SetWindowText(title);
	return TRUE;
}
void CAxisRotationDia::OnOK() {
	CDialog::OnOK();
	alpha = (PI * alpha) / 180;
}
/*
Member Load Dialog
*/
BEGIN_MESSAGE_MAP(CMemberLoadDia , CDialog)
	ON_CBN_SELCHANGE(IDC_SYSTEM,OnUpdateControls)
END_MESSAGE_MAP()

CMemberLoadDia::CMemberLoadDia(CWnd* parent,PLOADCASELIST pldcase,PSYSTEMLIST psystems,BOOL dist,LOADCASE* ccase)
    : CDialog(IDD_MEMBER_LOAD,parent) {
	for(int i = 0; i < 4;i++) {
		Q[i] = 0;
		if(i < 2) d[i] = i * 0.25;
		else d[i] = (i + 1) * 0.25;
	}
	Qu = 0;
	mode = F_REPLACE;
	dir = GRAVITYD;
	loadcases = pldcase;
	ploadcase = ccase;
	systems = psystems;
	loadtype = TFORCE;
	distributed = dist;
    use_relative = TRUE;
}
BOOL CMemberLoadDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_Cover.ShowWindow(!distributed);

	POSITION pos = loadcases->GetHeadPosition();
	while(pos) {
		c_List.AddString(loadcases->GetNext(pos).name);
	}
	c_List.SetCurSel(c_List.FindString(-1,ploadcase->name));

	pos = systems->GetHeadPosition();
	while(pos) {
		c_System.AddString(systems->GetNext(pos).name);
	}
	c_System.AddString("Local");
	c_System.SetCurSel(0);
	
	OnUpdateControls();
	return TRUE;
}
void CMemberLoadDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_LOADCASE,c_List);
	DDX_Control(pDX,IDC_LOADDIR,c_LoadDir);
	DDX_Control(pDX,IDC_COVER,c_Cover);
	DDX_Control(pDX,IDC_SYSTEM,c_System);
	for(int i = 0;i < 4;i++) {
		DDX_Text(pDX, IDC_DISTANCE1 + i, d[i]);
	    DDV_MinMaxDouble(pDX,d[i],0,1000000);
		DDX_Text(pDX, IDC_LOAD1 + i, Q[i]);
	}
   	DDX_Text(pDX, IDC_UNIFORM, Qu);
	DDX_Radio(pDX,IDC_ADD_TO_EXISTING,mode);
	DDX_CBIndex(pDX,IDC_LOADDIR,dir);
	DDX_Radio(pDX,IDC_FORCETYPE,loadtype);
	DDX_Check(pDX,IDC_RELATIVE,use_relative);
}
void CMemberLoadDia::OnOK() {
	CDialog::OnOK();
       
	if(dir == GRAVITYD) {
		for(int i = 0;i < 4;i++) Q[i] = -Q[i];
		Qu = -Qu;
		dir = IUZ;
	}
	if(loadtype == TMOMENT) 
		dir += IRX;

	LOADCASE ld;
	SYSTEM sys;
	POSITION pos;
	
	c_List.GetLBText(c_List.GetCurSel(),ld.name);
	if(pos = loadcases->Find(ld)) {
        ploadcase = &loadcases->GetAt(pos);
	} else {
		ploadcase = NULL;
	}
	
	
	c_System.GetLBText(c_System.GetCurSel(),sys.name);
	if(pos = systems->Find(sys)) {
		system = &systems->GetAt(pos);
	} else {
		system = NULL;
	}
}
void CMemberLoadDia::OnUpdateControls() {
    CString name;
	c_System.GetLBText(c_System.GetCurSel(),name);

	c_LoadDir.ResetContent();
	c_LoadDir.AddString("X");
    c_LoadDir.AddString("Y");
	c_LoadDir.AddString("Z");
	if(name == "Local") {
		if(dir == GRAVITYD) dir = 0;
		c_LoadDir.SetCurSel(dir);
	} else {
		c_LoadDir.AddString("Gravity");
		c_LoadDir.SetCurSel(dir);
	}
}
/*
Slab Distributed Load Dialog
*/
BEGIN_MESSAGE_MAP(CSlabDistLoadDia , CDialog)
	ON_CBN_SELCHANGE(IDC_SYSTEM,OnUpdateControls)
END_MESSAGE_MAP()

CSlabDistLoadDia::CSlabDistLoadDia(CWnd* parent,PLOADCASELIST pldcase,PSYSTEMLIST psystems,LOADCASE* ccase)
    : CDialog(IDD_SLAB_DISTRIBUTED,parent) {
	Qu = 0;
	mode = F_ADD;
	dir = GRAVITYD;
	loadcases = pldcase;
	ploadcase = ccase;
	systems = psystems;
	loadtype = TFORCE;
}
BOOL CSlabDistLoadDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos = loadcases->GetHeadPosition();
	while(pos) {
		c_List.AddString(loadcases->GetNext(pos).name);
	}
	c_List.SetCurSel(c_List.FindString(-1,ploadcase->name));
	
	pos = systems->GetHeadPosition();
	while(pos) {
		c_System.AddString(systems->GetNext(pos).name);
	}
	c_System.AddString("Local");
	c_System.SetCurSel(0);
	
	OnUpdateControls();
	return TRUE;
}
void CSlabDistLoadDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_LOADCASE,c_List);
	DDX_Control(pDX,IDC_LOADDIR,c_LoadDir);
	DDX_Control(pDX,IDC_SYSTEM,c_System);
	DDX_Text(pDX, IDC_UNIFORM, Qu);
	DDX_Radio(pDX,IDC_ADD_TO_EXISTING,mode);
	DDX_CBIndex(pDX,IDC_LOADDIR,dir);
	DDX_Radio(pDX,IDC_FORCETYPE,loadtype);
}
void CSlabDistLoadDia::OnOK() {
	CDialog::OnOK();
       
	if(dir == GRAVITYD) {
		Qu = -Qu;
		dir = IUZ;
	}
	if(loadtype == TMOMENT) 
		dir += IRX;

	LOADCASE ld;
	SYSTEM sys;
	POSITION pos;
	
	c_List.GetLBText(c_List.GetCurSel(),ld.name);
	if(pos = loadcases->Find(ld)) {
        ploadcase = &loadcases->GetAt(pos);
	} else {
		ploadcase = NULL;
	}
	
	
	c_System.GetLBText(c_System.GetCurSel(),sys.name);
	if(pos = systems->Find(sys)) {
		system = &systems->GetAt(pos);
	} else {
		system = NULL;
	}
}
void CSlabDistLoadDia::OnUpdateControls() {
    CString name;
	c_System.GetLBText(c_System.GetCurSel(),name);
	c_LoadDir.ResetContent();
	c_LoadDir.AddString("X");
    c_LoadDir.AddString("Y");
	c_LoadDir.AddString("Z");
	if(name == "Local") {
		if(dir == GRAVITYD) dir = 0;
		c_LoadDir.SetCurSel(dir);
	} else {
		c_LoadDir.AddString("Gravity");
		c_LoadDir.SetCurSel(dir);
	}
}
/*
Strain Load Dia
*/
CStrainLoadDia::CStrainLoadDia(CWnd* parent,PLOADCASELIST pldcase,LOADCASE* ccase,UINT ttype) :
                CDialog(IDD_MEMBER_STRAIN,parent) {
	Value = 0;
	type = TEMP_STRAIN;
	mode = F_ADD;
	loadcase = pldcase;
	ploadcase = ccase;
	stype = ttype;
}
BOOL CStrainLoadDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos = loadcase->GetHeadPosition();
	while(pos) {
		c_List.AddString(loadcase->GetNext(pos).name);
	}
	c_List.SetCurSel(c_List.FindString(-1,ploadcase->name));

	if(stype == TMEMBER) {
		c_Type.AddString("Temperature");
        c_Type.AddString("Axial Strain");
	} else {
		c_Type.AddString("Temperature");
		c_Type.AddString("Temperature Gradient ZZ");
        c_Type.AddString("Strain XX");
		c_Type.AddString("Strain YY");
		c_Type.AddString("Strain XY");
		c_Type.AddString("Curvature XX");
        c_Type.AddString("Curvature YY");
		c_Type.AddString("Curvature XY");
	}
	c_Type.SetCurSel(type);

	return TRUE;
}
void CStrainLoadDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_LOADCASE,c_List);
	DDX_Control(pDX,IDC_TYPE,c_Type);
   	DDX_Text(pDX,IDC_VALUE,Value);
	DDX_Radio(pDX,IDC_ADD_TO_EXISTING,mode);
}
void CStrainLoadDia::OnOK() {
	CDialog::OnOK();
	
	LOADCASE ld;
	POSITION pos;
	c_List.GetLBText(c_List.GetCurSel(),ld.name);
	if(pos = loadcase->Find(ld)) {
        ploadcase = &loadcase->GetAt(pos);
	} else {
		ploadcase = NULL;
	}

	type = c_Type.GetCurSel();
	if(stype == TSLAB) {
		if(type >= 2) {
			dir = IUX + type - 2;		
			type = FAB_STRAIN;
		} else {
			if(type == 1) dir = IRZ;
			else dir = IUX;
			type = TEMP_STRAIN;
		}
	}
}
/*
Rebar Dialog
*/
BEGIN_MESSAGE_MAP(CRebarDia , CDialog)
    ON_BN_CLICKED(IDC_COULMN,OnUpdateControls)
	ON_BN_CLICKED(IDC_BEAM,OnUpdateControls)
	ON_BN_CLICKED(IDC_RECTANGULAR,OnUpdateControls)
	ON_BN_CLICKED(IDC_CIRCULAR,OnUpdateControls)
END_MESSAGE_MAP()

CRebarDia::CRebarDia(CWnd* parent,REBAR* prebar) : 
                    CDialog(IDD_REBAR_DATA,parent) {
    rebar = prebar;
}
BOOL CRebarDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	OnUpdateControls();
	return TRUE;
}
void CRebarDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_COULMN,c_Coulmn);
	if(rebar->design == ASLAB) {
		int design = BEAM;
	    DDX_Radio(pDX, IDC_COULMN, design);
	} else {
	    DDX_Radio(pDX, IDC_COULMN, rebar->design);
	}

	DDX_Control(pDX,IDC_RECTANGULAR,c_Rectangular);
	DDX_Radio(pDX, IDC_RECTANGULAR, rebar->type);
	DDX_Control(pDX,IDC_CIRCULAR,c_Circular);
	DDX_Text(pDX, IDC_COVER, rebar->cover);
	DDX_Control(pDX, IDC_BARSIZE, c_BarSize);
	DDX_Text(pDX, IDC_BARSIZE, rebar->barsize);
	DDX_Control(pDX, IDC_SBARSIZE, c_SBarSize);
	DDX_Text(pDX, IDC_SBARSIZE, rebar->stirrup_barsize);
	DDX_Control(pDX,IDC_NZ,c_Nz);
	DDX_Text(pDX, IDC_NZ, rebar->nz);
	DDV_MinMaxInt(pDX,rebar->nz,0,15);
	DDX_Control(pDX,IDC_NY,c_Ny);
	DDX_Text(pDX, IDC_NY, rebar->ny);
	DDV_MinMaxInt(pDX,rebar->ny,0,15);
	DDX_Control(pDX,IDC_TOTAL,c_Total);
	DDX_Text(pDX, IDC_TOTAL, rebar->nt);
	DDV_MinMaxInt(pDX,rebar->nt,0,64);
}
void CRebarDia::OnUpdateControls() {
	BOOL column = (c_Coulmn.GetCheck() == BST_CHECKED);
	BOOL rectangular = (c_Rectangular.GetCheck() == BST_CHECKED);
	c_Coulmn.EnableWindow(rebar->design != ASLAB);
    c_Total.EnableWindow(!rectangular && column);
	c_Nz.EnableWindow(rectangular && column);
	c_Ny.EnableWindow(rectangular && column);
	c_BarSize.EnableWindow(column);
	c_SBarSize.EnableWindow(column);
	c_Rectangular.EnableWindow(column);
	c_Circular.EnableWindow(column);
}
/*
Material Property Dialog
*/
BEGIN_MESSAGE_MAP(CMaterialPropertyDia , CDialog)
    ON_LBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
END_MESSAGE_MAP()

CMaterialPropertyDia::CMaterialPropertyDia(CWnd* parent,MATERIAL* pmaterial) : 
                    CDialog(IDD_MATERIAL_PROPERTY,parent) {
    material = pmaterial;
}
BOOL CMaterialPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_List.AddString("Concrete");
    c_List.AddString("Steel");
	c_List.AddString("Filler");
	c_List.SetCurSel(material->type);
	if(material->rank) {
		c_Name.EnableWindow(FALSE);
		c_List.EnableWindow(FALSE);
	}
	OnUpdateControls();
	return TRUE;
}
void CMaterialPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_NAME,c_Name);
	DDX_Text(pDX, IDC_NAME, material->name);
	DDX_Control(pDX, IDC_ELASTICITY, c_Elasticity);
	DDX_Text(pDX, IDC_ELASTICITY, material->E);
	DDX_Control(pDX, IDC_POISON, c_Poisson);
	DDX_Text(pDX, IDC_POISON, material->nu);
	DDX_Text(pDX, IDC_UNITWEIGHT, material->unitweight);
	DDX_Text(pDX, IDC_DENSITY, material->density);
	DDX_Control(pDX, IDC_THERMAL_COEFF, c_Thermal);
	DDX_Text(pDX, IDC_THERMAL_COEFF, material->alphac);
	DDX_Text(pDX, IDC_SHEAR_MODULUS, material->G);
	DDX_Control(pDX,IDC_FCK,c_Fck);
	DDX_Text(pDX, IDC_FCK, material->fck);
	DDX_Control(pDX,IDC_FCTK,c_Fctk);
	DDX_Text(pDX, IDC_FCTK, material->fctk);
    DDX_Control(pDX, IDC_FYK, c_Fyk);
	DDX_Text(pDX, IDC_FYK, material->fyk);
	DDX_Control(pDX, IDC_FYKS, c_Fyks);
	DDX_Text(pDX, IDC_FYKS, material->fyks);
	DDX_Control(pDX,IDC_FTK,c_Ftk);
	DDX_Text(pDX, IDC_FTK, material->ftk);
	DDX_Control(pDX,IDC_FTKS,c_Ftks);
	DDX_Text(pDX, IDC_FTKS, material->ftks);
	DDX_Control(pDX,IDC_ES,c_Es);
	DDX_Text(pDX, IDC_ES, material->Es);
	DDX_Control(pDX,IDC_TYPE,c_List);
}
void CMaterialPropertyDia::OnOK() {
	CDialog::OnOK();
	material->G = material->E / (2 * (material->nu + 1));
	material->type = MATERIAL::CONCRETE + c_List.GetCurSel();
}
void CMaterialPropertyDia::OnUpdateControls() {
	int type = c_List.GetCurSel();
	c_Ftk.EnableWindow(type == MATERIAL::STEEL);
	c_Ftks.EnableWindow(type == MATERIAL::STEEL);
	c_Fck.EnableWindow(type == MATERIAL::CONCRETE);
	c_Fctk.EnableWindow(type == MATERIAL::CONCRETE);
    c_Es.EnableWindow(type == MATERIAL::CONCRETE);
	c_Fyk.EnableWindow(type != MATERIAL::FILLER);
	c_Fyks.EnableWindow(type != MATERIAL::FILLER);
	c_Elasticity.EnableWindow(type != MATERIAL::FILLER);
	c_Poisson.EnableWindow(type != MATERIAL::FILLER);
	c_Thermal.EnableWindow(type != MATERIAL::FILLER);
}
/*
Section Stiffness Modifiers dialog
*/
BEGIN_MESSAGE_MAP(CSectionModifiersDia , CDialog)
END_MESSAGE_MAP()

CSectionModifiersDia::CSectionModifiersDia(CWnd* parent,SECTION* psection) : 
                    CDialog(IDD_SECTION_MODIFIERS,parent),
					c_Table(CGridCtrl::TYPE3) {
    section = psection;
}
void CSectionModifiersDia::OnOK() {
	c_Table.GetData();
	CDialog::OnOK();
}
BOOL CSectionModifiersDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_Table.SetTwoHeaders();
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fA,"Axial Area");
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fAy,"Shear Area Y");
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fAz,"Shear Area Z");
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fIx,"Torsion");
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fIy,"Moment of Inertia Y");
	c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&section->fIz,"Moment of Inertia Z");

    return TRUE;
}
void CSectionModifiersDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_LIST,c_Table);
}
/*
Section Property Dialog
*/
BEGIN_MESSAGE_MAP(CSectionPropertyDia , CDialog)
    ON_WM_PAINT()
	ON_LBN_SELCHANGE(IDC_MATERIAL_LIST,OnUpdateControls)
    ON_LBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
    ON_COMMAND(IDC_EDIT_MATERIAL,OnEditMaterial)
	ON_COMMAND(IDC_EDIT_REBAR,OnEditRebar)
    ON_COMMAND(IDC_SET_MODIFIERS,OnSetModifiers)
    ON_COMMAND(IDC_CHART,OnChart)
END_MESSAGE_MAP()

CSectionPropertyDia::CSectionPropertyDia(CWnd* parent,SECTION* psection,PMATERIALLIST pmatdata) : 
                    CDialog(IDD_SECTION_PROPERTY,parent) {
    materials = pmatdata;
    section = psection;
}
BOOL CSectionPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	/*material*/
	POSITION pos = materials->GetHeadPosition();
	while(pos) c_List.AddString(materials->GetNext(pos).name);
	c_List.SetCurSel(c_List.FindString(0,section->material->name));
    
	/*type*/
	c_TypeList.AddString("Rectangular");
	c_TypeList.AddString("Circular");
    c_TypeList.AddString("Wide Flange");
	c_TypeList.AddString("Channel");
    c_TypeList.AddString("Double Channel");
	c_TypeList.AddString("Tee");
	c_TypeList.AddString("Angle");
	c_TypeList.AddString("Double Angle");
	c_TypeList.AddString("Box");
	c_TypeList.AddString("Pipe");
	c_TypeList.AddString("General");
	c_TypeList.SetCurSel(section->type);
	
	/*controls*/
	if(section->rank) c_Name.EnableWindow(FALSE);
	OnUpdateControls();

	return TRUE;
}
void CSectionPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_SECTION_NAME,c_Name);
	DDX_Text(pDX, IDC_SECTION_NAME, section->name);
	DDX_Control(pDX, IDC_WIDTH,c_Width);
	DDX_Text(pDX, IDC_WIDTH, section->w);
	DDX_Control(pDX, IDC_HEIGHT,c_Height);
	DDX_Text(pDX, IDC_HEIGHT, section->h);
	DDX_Control(pDX, IDC_RADIUS,c_Radius);
	DDX_Text(pDX, IDC_RADIUS,section->r);
	DDX_Control(pDX, IDC_TW,c_Tw);
	DDX_Text(pDX, IDC_TH, section->tf);
    DDX_Control(pDX, IDC_TH,c_Th);
	DDX_Text(pDX, IDC_BDIST, section->bd);
    DDX_Control(pDX, IDC_BDIST,c_Bd);
	DDX_Text(pDX, IDC_TW, section->tw);
	DDX_Control(pDX, IDC_AREA, c_A);
    DDX_Control(pDX, IDC_AY, c_Ay);
	DDX_Control(pDX, IDC_AZ, c_Az);
	DDX_Control(pDX, IDC_IXX, c_Ix);
	DDX_Control(pDX, IDC_IYY, c_Iy);
	DDX_Control(pDX, IDC_IZZ, c_Iz);
	DDX_Control(pDX, IDC_ZY, c_Zy);
	DDX_Control(pDX, IDC_ZZ, c_Zz);
	DDX_Control(pDX, IDC_SY, c_Sy);
	DDX_Control(pDX, IDC_SZ, c_Sz);
	DDX_Control(pDX, IDC_RY, c_ry);
	DDX_Control(pDX, IDC_RZ, c_rz);
	DDX_Control(pDX,IDC_MATERIAL_LIST,c_List);
	DDX_Text(pDX, IDC_MATERIAL_LIST, section->material->name);
	DDX_Control(pDX,IDC_DRAWING,c_Drawing);
	DDX_Control(pDX,IDC_EDIT_REBAR,c_Rebar);
	DDX_Control(pDX,IDC_CHART,c_Chart);
	DDX_Control(pDX,IDC_TYPE,c_TypeList);
}
static void DrawRebar(CPaintDC& paintDC,CPoint& p) {
	CRect rect(p.x - 4,p.y - 4,p.x + 4,p.y + 4);
	CBrush Brush(RGB(255,0,0));
	CBrush* pOldBrush = paintDC.SelectObject(&Brush);
	paintDC.Ellipse(&rect);
	paintDC.SelectObject(pOldBrush);
}

void CSectionPropertyDia::OnPaint() {
	CDialog::OnPaint();
	DOUBLE w = section->w,h = section->h,r = section->r,tw = section->tw,tf = section->tf,
		   bd = section->bd,c = section->rebar.cover;
	int    nz = section->rebar.nz,ny = section->rebar.ny,
		   nt = section->rebar.nt;

	if(section->type == CIRCULAR || section->type == PIPE) {
		h = w = 2 * r;
	}
    /*paint section drawing*/
	CPaintDC paintDC(&c_Drawing);
	CRect rect;
	c_Drawing.GetClientRect(&rect);
	CBrush wbrush(RGB(255,255,255));
	CBrush* oldBrush = paintDC.SelectObject(&wbrush);
	CPoint Center = rect.CenterPoint();
	int H = rect.Height(),W = rect.Width();
	paintDC.Rectangle(rect);

	CRgn rgn;
	rgn.CreateRectRgn(0,0,W,H);
	paintDC.SelectClipRgn(&rgn);

	CBrush brush(RGB(100,100,100));
	CPen pen(PS_SOLID,1,RGB(200,200,200));
	CPen* oldPen = paintDC.SelectObject(&pen);
	CPoint p1,p2;
	DOUBLE Scale = (H - 40) / h,Scale1 = (W - 40) / w;
	Scale = min(Scale,Scale1);

	/*grids*/
	const int SIZE = 10;
	int x = 0;
	while(x < W) {
		p1.x = x;
		p1.y = 0;
		p2.x = x;
		p2.y = H;
		paintDC.MoveTo(p1);
		paintDC.LineTo(p2);
        x += SIZE;
	}
	x = 0;
	while(x < H) {
		p1.x = 0;
		p1.y = x;
		p2.x = W;
		p2.y = x;
		paintDC.MoveTo(p1);
		paintDC.LineTo(p2);
        x += SIZE;
	}
    
	/*section*/
	paintDC.SelectObject(&brush);

	CPen pen1(PS_SOLID,1,RGB(0,0,0));
	paintDC.SelectObject(&pen1);

#define DrawRect(xt,yt,wt,ht) {\
    p1.x = Center.x - int(Scale * (xt));\
	p1.y = Center.y - int(Scale * (yt));\
	p2.x = Center.x + int(Scale * ((wt) - (xt)));\
	p2.y = Center.y + int(Scale * ((ht) - (yt)));\
	if(p1.x == p2.x) p2.x++;\
	if(p1.y == p2.y) p2.y++;\
	paintDC.Rectangle(p1.x,p1.y,p2.x,p2.y);\
};

	if(section->type == RECTANGULAR) {
		DrawRect(w/2,h/2,w,h);
	} else if(section->type == CIRCULAR) {
		int R = int(Scale * r);
		CRect rect(Center.x - R,Center.y - R,Center.x + R,Center.y + R);
		paintDC.Ellipse(&rect);
	} else if(section->type == WIDEFLANGE) {
		DrawRect(tw/2,(h - tf)/2,tw,h - tf);
        DrawRect(w/2,(h + tf)/2,w,tf);
        DrawRect(w/2,-(h - tf)/2,w,tf);
	} else if(section->type == CHANNEL) {
		DrawRect(w/3,(h - tf)/2,tw,h - tf);
        DrawRect(w/3,(h + tf)/2,w + tw/2,tf);
        DrawRect(w/3,-(h - tf)/2,w + tw/2,tf);
	} else if(section->type == DOUBLE_CHANNEL) {
		DrawRect(-bd/2,(h - tf)/2,tw,h - tf);
        DrawRect(-bd/2,(h + tf)/2,w,tf);
        DrawRect(-bd/2,-(h - tf)/2,w,tf);
		DrawRect(tw + bd/2,(h - tf)/2,tw,h - tf);
        DrawRect(w + tw/2 + bd/2,(h + tf)/2,w + tw/2,tf);
        DrawRect(w + tw/2 + bd/2,-(h - tf)/2,w + tw/2,tf);
	} else if(section->type == TEE) {
		DrawRect(tw / 2,-h / 6 + (h - tf) / 2,tw,h - tf / 2);
        DrawRect( (w + tw) / 2,-h / 6 + (h - tf)/2,w + tw,tf);
	} else if(section->type == ANGLE) {
		DrawRect(w / 3, 2 * h / 3,tw,h - tf / 2);
        DrawRect(w / 3, -h / 3 + tf / 2,w + tw / 2,tf);
	} else if(section->type == DOUBLE_ANGLE) {
		DrawRect(-bd/2, 2 * h / 3,tw,h - tf / 2);
        DrawRect(-bd/2, -h / 3 + tf / 2,w + tw / 2,tf);
		DrawRect(tw + bd/2, 2 * h / 3,tw,h - tf / 2);
        DrawRect(w + tw / 2 + bd/2, -h / 3 + tf / 2,w + tw / 2,tf);
	} else if(section->type == BOX) {
		DrawRect((w + tw) / 2,(h + tf)/2,w + tw,tf);
        DrawRect((w + tw) / 2,-(h - tf)/2,w + tw,tf);
		DrawRect((w + tw) / 2,(h - tf)/2,tw,h - tf);
        DrawRect(-(w - tw)/2,(h - tf)/2,tw,h - tf);
	} else if(section->type == PIPE) {
		CBrush brush1;
		brush1.CreateStockObject(NULL_BRUSH);
	    paintDC.SelectObject(brush1);
		int R = int(Scale * r);
		CRect rect(Center.x - R,Center.y - R,Center.x + R,Center.y + R);
		paintDC.Ellipse(&rect);
		R = int(Scale * (r - tw));
		rect = CRect(Center.x - R,Center.y - R,Center.x + R,Center.y + R);
		paintDC.Ellipse(&rect);
	}
#undef DrawRect
	paintDC.SelectObject(oldBrush);

	/*axes*/
	CPen* pPen = new CPen(PS_SOLID,1,RGB(0,255,0));
	paintDC.SelectObject(pPen);
	paintDC.MoveTo(Center);
	p2 = Center;
	p2.x += W/2;
	paintDC.LineTo(p2);
	delete pPen;
	pPen = new CPen(PS_SOLID,1,RGB(0,0,255));
	paintDC.SelectObject(pPen);
	paintDC.MoveTo(Center);
	p2 = Center;
	p2.y -= H/2;
	paintDC.LineTo(p2);
	delete pPen;

	paintDC.SelectObject(oldPen);
	/*reinforcement*/
	if(section->material->type == MATERIAL::CONCRETE && section->rebar.design == COULMN) {
		int rebtype = section->rebar.type,i;
		if(rebtype == RECTANGULAR) {
			for(i = 0;i < nz;i++) {
				p1.x = Center.x - int(Scale * (w / 2 - i * (w - 2 * c) / (nz - 1) - c));
				p1.y = Center.y - int(Scale * (h / 2 - c));
				DrawRebar(paintDC,p1);
				p1.x = Center.x - int(Scale * (w / 2 - i * (w - 2 * c) / (nz - 1) - c));
				p1.y = Center.y - int(Scale * ( -h / 2 + c));
				DrawRebar(paintDC,p1);
			}
			for(i = 1;i < ny - 1;i++) {
				p1.x = Center.x - int(Scale * (w / 2 - c));
				p1.y = Center.y - int(Scale * (h / 2 - i * (h - 2 * c) / (ny - 1) - c));
				DrawRebar(paintDC,p1);
				p1.x = Center.x - int(Scale * (-w / 2 + c));
				p1.y = Center.y - int(Scale * (h / 2 - i * (h - 2 * c) / (ny - 1) - c));
				DrawRebar(paintDC,p1);
			}
		} else if(rebtype == CIRCULAR) {
			DOUBLE q = 2 * PI / nt;
			for(i = 0;i < nt;i++) {
				p1.x = Center.x - int(Scale * ((r - c) * sin(q * i)));
				p1.y = Center.y - int(Scale * ((r - c) * cos(q * i)));
				DrawRebar(paintDC,p1);
			}
		}
	}
	/*release*/
	rgn.DeleteObject();
}
void CSectionPropertyDia::OnOK() {
	CDialog::OnOK();
	if(section->rank != LEVEL4)
		section->GetData();
}
void CSectionPropertyDia::OnChart() {
	DESIGN design;
	design.detail = SECTION::pDocDetailing;
	design.InitColumn(section);
	CChartDia ChartDia(this,&design);
	ChartDia.DoModal();
}
void CSectionPropertyDia::OnEditMaterial() {
	MATERIAL mat;
	POSITION pos;
	int index = c_List.GetCurSel();
	c_List.GetLBText(index,mat.name);
  	if(pos = materials->Find(mat)) {
		MATERIAL* pmat = &materials->GetAt(pos);
		mat = *pmat;
		CMaterialPropertyDia MaterialDia(this,&mat);
		if(MaterialDia.DoModal() == IDOK) {
			*pmat = mat;
			c_List.DeleteString(index);
			c_List.SetCurSel(c_List.AddString(mat.name));
			UpdateData(FALSE);
		}
	}
}
void CSectionPropertyDia::OnEditRebar() {
	REBAR reb = section->rebar;
	CRebarDia RebarDia(this,&reb);
	if(RebarDia.DoModal() == IDOK) {
		section->rebar = reb;
		Invalidate();
	}
}
void CSectionPropertyDia::OnSetModifiers() {
	CSectionModifiersDia Dia(this,section);
	Dia.DoModal();
}
void CSectionPropertyDia::OnUpdateControls() {
    MATERIAL mat,*steel,*conc;
	POSITION pos;
	int prev_type = section->type;

	section->type = c_TypeList.GetCurSel();
	c_List.GetLBText(c_List.GetCurSel(),mat.name);
	if(pos = materials->Find(mat)) {
        section->material = &materials->GetAt(pos);
		c_Rebar.ShowWindow(section->material->type == MATERIAL::CONCRETE
			&& section->type <= CIRCULAR);
		Invalidate();
	}

	mat.name = "STEEL";
	if(pos = materials->Find(mat)) {
		steel = &materials->GetAt(pos);
	}
	mat.name = "CONC";
	if(pos = materials->Find(mat)) {
		conc = &materials->GetAt(pos);
	}
	
	section->rebar.type = section->type;
	if(prev_type != section->type) {
		if(section->type == WIDEFLANGE ||
			section->type == CHANNEL || 
			section->type == DOUBLE_CHANNEL || 
			section->type == TEE ||
			section->type == ANGLE || 
			section->type == DOUBLE_ANGLE ||
			section->type == BOX
			) {
			section->w = 0.1;
			section->h = 0.17;
			section->tw = 0.006;
			section->tf = 0.006;
			section->bd = 0;
			section->material = steel;
		} else if(section->type == PIPE) {
			section->r = 0.08;
			section->tw = 0.006;
			section->material = steel;
		} else {
			section->w = 0.4;
			section->h = 0.4;
            section->r = 0.2;
			section->rebar.design = COULMN;
			section->material = conc;
		}
		section->GetData();
	}

	BOOL enable = (section->type == GENERAL);
	c_A.EnableWindow(enable);
    c_Ay.EnableWindow(enable);
	c_Az.EnableWindow(enable);
	c_Ix.EnableWindow(enable);
	c_Iy.EnableWindow(enable);
	c_Iz.EnableWindow(enable);
	c_Zy.EnableWindow(enable);
    c_Zz.EnableWindow(enable);
	c_Sy.EnableWindow(enable);
	c_Sz.EnableWindow(enable);
	c_ry.EnableWindow(enable);
	c_rz.EnableWindow(enable);

	if(section->rank == LEVEL4 || section->type == GENERAL) {
		if(section->rank == LEVEL4) c_TypeList.EnableWindow(FALSE);
		c_Width.EnableWindow(FALSE);
		c_Height.EnableWindow(FALSE);
		c_Radius.EnableWindow(FALSE);
		c_Tw.EnableWindow(FALSE);
		c_Th.EnableWindow(FALSE); 
		c_Bd.EnableWindow(FALSE);
	} else {
		BOOL has_w = (section->type == RECTANGULAR || section->type == WIDEFLANGE || section->type == CHANNEL || 
			section->type == DOUBLE_CHANNEL || section->type == TEE  || section->type == ANGLE  || section->type == DOUBLE_ANGLE
			|| section->type == BOX);
		BOOL has_h = (has_w);
		BOOL is_full = (section->type == RECTANGULAR || section->type == CIRCULAR);
		c_Width.EnableWindow(has_w);
		c_Height.EnableWindow(has_h);
		c_Radius.EnableWindow(section->type == CIRCULAR || section->type == PIPE);
		c_Tw.EnableWindow(has_w && !is_full);
		c_Th.EnableWindow(has_h && !is_full); 
		c_Bd.EnableWindow(section->type == DOUBLE_CHANNEL || section->type == DOUBLE_ANGLE); 
	}
	c_Chart.ShowWindow(section->rebar.design == COULMN 
		&& section->type <= CIRCULAR 
		&& section->material->type == MATERIAL::CONCRETE);

	/*set values*/
	c_A.SetWindowText(section->A);
	c_Ay.SetWindowText(section->Ay);
	c_Az.SetWindowText(section->Az);
	c_Ix.SetWindowText(section->Ix);
	c_Iy.SetWindowText(section->Iy);
	c_Iz.SetWindowText(section->Iz);
	c_Zy.SetWindowText(section->Zy);
	c_Zz.SetWindowText(section->Zz);
	c_Sy.SetWindowText(section->Sy);
	c_Sz.SetWindowText(section->Sz);
	c_ry.SetWindowText(section->ry);
	c_rz.SetWindowText(section->rz);

	UpdateData(FALSE);
}
/*
Section Property Dialog
*/
BEGIN_MESSAGE_MAP(CASectionPropertyDia , CDialog)
	ON_LBN_SELCHANGE(IDC_MATERIAL_LIST,OnUpdateControls)
    ON_LBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
    ON_COMMAND(IDC_EDIT_MATERIAL,OnEditMaterial)
	ON_COMMAND(IDC_EDIT_REBAR,OnEditRebar)
END_MESSAGE_MAP()

CASectionPropertyDia::CASectionPropertyDia(CWnd* parent,ASECTION* psection,PMATERIALLIST pmatdata) : 
                    CDialog(IDD_ASECTION_PROPERTY,parent) {
    materials = pmatdata;
    section = psection;
}
BOOL CASectionPropertyDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	/*material*/
	POSITION pos = materials->GetHeadPosition();
	while(pos) c_List.AddString(materials->GetNext(pos).name);
	c_List.SetCurSel(c_List.FindString(0,section->material->name));
    
	/*type*/
	c_TypeList.AddString("Q4");
	c_TypeList.AddString("Q5");
	c_TypeList.AddString("Q6");
	c_TypeList.AddString("Q7");
	c_TypeList.AddString("Q8");
	c_TypeList.AddString("Q9");
	c_TypeList.AddString("Q4D");
	c_TypeList.AddString("PLATE4");
    c_TypeList.AddString("SHELL4");
	c_TypeList.SetCurSel(section->type - ASECTION::Q4);
	
	/*controls*/
	if(section->rank) c_Name.EnableWindow(FALSE);
	OnUpdateControls();

	return TRUE;
}
void CASectionPropertyDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX,IDC_SECTION_NAME,c_Name);
	DDX_Text(pDX, IDC_SECTION_NAME, section->name);
	DDX_Control(pDX, IDC_HEIGHT,c_Height);
	DDX_Text(pDX, IDC_HEIGHT, section->h);
	DDX_Control(pDX,IDC_MATERIAL_LIST,c_List);
	DDX_Text(pDX, IDC_MATERIAL_LIST, section->material->name);
	DDX_Control(pDX,IDC_EDIT_REBAR,c_Rebar);
	DDX_Control(pDX,IDC_TYPE,c_TypeList);
	DDX_Control(pDX,IDC_STATE_TYPE,c_StateTypeList);
}

void CASectionPropertyDia::OnOK() {
	section->type = c_TypeList.GetCurSel() + ASECTION::Q4;
	if(section->type == ASECTION::PLATE4)
		section->state_type = ASECTION::PLATE_ACTION;
	else
		section->state_type = c_StateTypeList.GetCurSel();
	CDialog::OnOK();
}
void CASectionPropertyDia::OnEditMaterial() {
	MATERIAL mat;
	POSITION pos;
	int index = c_List.GetCurSel();
	c_List.GetLBText(index,mat.name);
  	if(pos = materials->Find(mat)) {
		MATERIAL* pmat = &materials->GetAt(pos);
		mat = *pmat;
		CMaterialPropertyDia MaterialDia(this,&mat);
		if(MaterialDia.DoModal() == IDOK) {
			*pmat = mat;
			c_List.DeleteString(index);
			c_List.SetCurSel(c_List.AddString(mat.name));
			UpdateData(FALSE);
		}
	}
}
void CASectionPropertyDia::OnEditRebar() {
	REBAR reb = section->rebar;
	CRebarDia RebarDia(this,&reb);
	if(RebarDia.DoModal() == IDOK) {
		section->rebar = reb;
	}
}
void CASectionPropertyDia::OnUpdateControls() {
	MATERIAL mat;
	POSITION pos;

	section->type = c_TypeList.GetCurSel() + ASECTION::Q4;
	if(section->type == ASECTION::SHELL4)
		section->state_type = ASECTION::SHELL_ACTION;
	else if(section->type == ASECTION::PLATE4)
		section->state_type = ASECTION::PLATE_ACTION;
	else
		section->state_type = ASECTION::PLANE_STRESS;

	/*state type*/
	if(section->type != ASECTION::PLATE4) {
		c_StateTypeList.ResetContent();
		c_StateTypeList.AddString("Plane Stress");
		c_StateTypeList.AddString("Plane Strain");
		if(section->type == ASECTION::SHELL4) {
			c_StateTypeList.AddString("Plate action");
			c_StateTypeList.AddString("Shell action");
		}
	}
	c_StateTypeList.SetCurSel(section->state_type);
	c_StateTypeList.ShowWindow(section->type != ASECTION::PLATE4);

	/*Rebar*/
	c_List.GetLBText(c_List.GetCurSel(),mat.name);
	if(pos = materials->Find(mat)) {
        section->material = &materials->GetAt(pos);
		c_Rebar.ShowWindow(section->material->type == MATERIAL::CONCRETE);
	}
	section->rebar.type = section->type;
	section->rebar.design = ASLAB;
}
/*
Joint Information Dialog
*/
BEGIN_MESSAGE_MAP(CJointInfoDia , CDialog)
   	ON_COMMAND(IDC_EDIT_RESTRAINT,OnEditRestraint)
	ON_COMMAND(IDC_EDIT_AXIS,OnEditAxis)
END_MESSAGE_MAP()

CJointInfoDia::CJointInfoDia(CWnd* parent,JOINT* pjoint,PLOADCASELIST pldcase) : 
               CDialog(IDD_JOINT_DATA,parent),
			   c_Jpoint(CGridCtrl::TYPE3) {
	joint = pjoint;
	loadcase = pldcase;
}

void CJointInfoDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_FORCE_LIST, c_FList);
	DDX_Control(pDX, IDC_JPOINT,c_Jpoint);
}
void CJointInfoDia::OnEditRestraint() {
	CJointRestraintDia JointDia(this,joint);
	if(JointDia.DoModal() == IDOK) {
		UBMP8 restraint = 0;
		for(int i = IUX;i <= IRZ;i++) {
			if(JointDia.res[i]) restraint |= (1 << i);
		}
		joint->restraint = restraint;
	}
}
void CJointInfoDia::OnEditAxis() {
	CJointAxisDia JointDia(this,joint);
	if(JointDia.DoModal() == IDOK)
		joint->rotation = JointDia.rotation;
}
void CJointInfoDia::OnOK() {
	c_Jpoint.GetData();
	CDialog::OnOK();
}
BOOL CJointInfoDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	CString str;
	JLOAD* pload;
    POSITION pos = joint->load.GetHeadPosition();
	while(pos) {
		pload = &joint->load.GetNext(pos);
		c_FList.AddString(LoadTypeStr[pload->type]);

		if(pload->type != MASS) {
			str = "      Load case:  " + pload->loadcase->name;
			c_FList.AddString(str);
		}

		str = "      Coordinate: ";
		if(pload->system) str += pload->system->name;
		else str += "LOCAL";
		c_FList.AddString(str);

		for(int i = 0;i < 6;i++) {
			if(!EQUAL(pload->Q[i],0)) {
				if(pload->type == DISPLACEMENT) str.Format("      %s : %.3g",DirString[i],pload->Q[i]);
				else str.Format("      %s : %.2f",DirString[i],pload->Q[i]);
			    c_FList.AddString(str);
			}
		}
	}

	DOUBLE* jp;
	c_Jpoint.SetTwoHeaders();
	for(UINT i = 0;i < 3;i++) {
		if(i == 0) jp = &joint->p.x;
		else if(i == 1) jp = &joint->p.y;
		else jp = &joint->p.z;

		str.Format("%c",'X' + i);
		c_Jpoint.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,jp,str);
	}

	return TRUE;
}
/*
Member Information Dialog
*/
BEGIN_MESSAGE_MAP(CMemberInfoDia , CDialog)
   	ON_COMMAND(IDC_EDIT_SECTION,OnEditSection)
	ON_COMMAND(IDC_EDIT_RELEASE,OnEditRelease)
	ON_COMMAND_RANGE(IDC_EDIT_JI,IDC_EDIT_JII,OnEditJoint)
	ON_COMMAND(IDC_EDIT_AXIS,OnEditAxis)
	ON_COMMAND(IDC_EDIT_DIVISIONS,OnEditDivisions)
END_MESSAGE_MAP()

CMemberInfoDia::CMemberInfoDia(CWnd* parent,MEMBER* pmember,PSECTIONLIST psections,
							   PMATERIALLIST pmaterials,PLOADCASELIST pldcases) : 
               CDialog(IDD_MEMBER_DATA,parent) {
    member = pmember;
	sections = psections;
	materials = pmaterials;
	loadcases = pldcases;
}
void CMemberInfoDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_FORCE_LIST, c_FList);
}
BOOL CMemberInfoDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos;

	/*load*/
	CString str;
	LOAD* pload;
    pos = member->load.GetHeadPosition();
	while(pos) {
		pload = &member->load.GetNext(pos);
		
		c_FList.AddString(LoadTypeStr[pload->type]);
		str = "      Load case:  " + pload->loadcase->name;
		c_FList.AddString(str);

		if(pload->is_member_load()) {
			
			str = "      Coordinate: ";
			if(pload->system) str += pload->system->name;
			else str += "LOCAL";
			c_FList.AddString(str);

			str = "      Direction : " + DirString[pload->dir];
			c_FList.AddString(str);

			if(pload->type == CONCENTRATED) {
				str.Format("      P         : %.2f",pload->P);
				c_FList.AddString(str);
				str.Format("      X         : %.2f",pload->x);
				c_FList.AddString(str);
			} else if(pload->type == UNIFORM) {
				str.Format("      W         : %.2f",pload->P);
				c_FList.AddString(str);
			} else if(pload->type == TRAPEZOIDAL) {
				str.Format("      P1         : %.2f",pload->P);
				c_FList.AddString(str);
				str.Format("      X1         : %.2f",pload->x);
				c_FList.AddString(str);
				str.Format("      P2         : %.2f",pload->P1);
				c_FList.AddString(str);
				str.Format("      X2         : %.2f",pload->x1);
				c_FList.AddString(str);
			}
		} else {
			if(pload->type == TEMPERATURE) {
				str.Format("      DT      : %.2f",pload->P);
				c_FList.AddString(str);
			} else if(pload->type == STRAIN) {
				str.Format("      DL/L      : %.3g",pload->P);
				c_FList.AddString(str);
			}
		}
	}
	
	return TRUE;
}
void CMemberInfoDia::OnEditSection() {
	if(member->e_section) {
		CTaperedSectionDia Dia(this,sections,materials,member);
		Dia.DoModal();
		UpdateData(FALSE);
	} else {
		CSectionPropertyDia Dia(this,member->section,materials);
		Dia.DoModal();
		UpdateData(FALSE);
	}
}
void CMemberInfoDia::OnEditRelease() {
	CMemberReleaseDia MemberDia(AfxGetMainWnd(),member);
    if(MemberDia.DoModal() == IDOK) {

		UBMP8 nrelease = NOUR,frelease = NOUR;
		if(MemberDia.r[0]) nrelease |= UX;
		if(MemberDia.r[1]) nrelease |= UY;
		if(MemberDia.r[2]) nrelease |= UZ;
		if(MemberDia.r[3]) nrelease |= RX;
		if(MemberDia.r[4]) nrelease |= RY;
		if(MemberDia.r[5]) nrelease |= RZ;
		if(MemberDia.r[6]) frelease |= UX;
		if(MemberDia.r[7]) frelease |= UY;
		if(MemberDia.r[8]) frelease |= UZ;
		if(MemberDia.r[9]) frelease |= RX;
		if(MemberDia.r[10]) frelease |= RY;
		if(MemberDia.r[11]) frelease |= RZ;

		member->nrelease = nrelease;
		member->frelease = frelease;
	}
}
void CMemberInfoDia::OnEditJoint(UINT nID) {
	CJointInfoDia* JointDia;
	if(nID == IDC_EDIT_JI)
		JointDia = new CJointInfoDia(this,member->j1,loadcases);
	else
		JointDia = new CJointInfoDia(this,member->j2,loadcases);
	JointDia->DoModal();
	delete JointDia;
}
void CMemberInfoDia::OnEditAxis() {
	CAxisRotationDia Dia(this,member->alpha,"Member local axis");
	if(Dia.DoModal() == IDOK) {
		member->alpha = Dia.alpha;
	}
}
void CMemberInfoDia::OnEditDivisions() {
	CMemberDivisionsDia Dia(AfxGetMainWnd());
	Dia.nMinDiv = member->nMinDiv;
	Dia.nMinFrameDiv = member->nMinFrameDiv; 
	Dia.DivAtInterim = member->DivAtInterim;
    if(Dia.DoModal() == IDOK) {
		member->nMinDiv = Dia.nMinDiv;
		member->nMinFrameDiv = Dia.nMinFrameDiv;
		member->DivAtInterim = Dia.DivAtInterim;
	}
}
void CMemberInfoDia::OnOK() {
	CDialog::OnOK();
}
/*
Member break dialog
*/
CMemberBreakDia::CMemberBreakDia(CWnd* parent) : 
                    CDialog(IDD_MEMBER_BREAK,parent) {
    breakn = 2;
	breakr = 1.0;
}
void CMemberBreakDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_BREAKN, breakn);
	DDV_MinMaxInt(pDX,breakn,0,100);
	DDX_Text(pDX, IDC_BREAKR, breakr);
	DDV_MinMaxDouble(pDX,breakr,0,100);
}
/*
Tapered section dialog
*/
BEGIN_MESSAGE_MAP(CTaperedSectionDia , CDialog)
	ON_COMMAND_RANGE(IDC_EDITSS,IDC_EDITES,OnEditSection)
END_MESSAGE_MAP()

CTaperedSectionDia::CTaperedSectionDia(CWnd* parent,PSECTIONLIST psections,PMATERIALLIST pmaterials,MEMBER* pmember) : 
                    CDialog(IDD_TAPERED_SECTION,parent) {
	member = pmember;
	sections = psections;
	materials = pmaterials;
	if(member) {
		EIyy = member->EIyy;
		EIzz = member->EIzz;
		s_sec = member->section;
		e_sec = member->e_section;
		soff = member->start_offset;
		eoff = member->end_offset;
	} else {
		EIyy = MEMBER::PARABOLIC;
		EIzz = MEMBER::LINEAR;
		s_sec = &sections->GetHead();
		e_sec = &sections->GetHead();
		soff = 0;
		eoff = 0;
	}
}
void CTaperedSectionDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_START, c_Start);
	DDX_Control(pDX, IDC_END, c_End);
    DDX_Control(pDX, IDC_EIYY, c_EIyy);
	DDX_Control(pDX, IDC_EIZZ, c_EIzz);
	DDX_Text(pDX,IDC_SOFFSET,soff);
	DDV_MinMaxDouble(pDX,soff,0,1);
    DDX_Text(pDX,IDC_EOFFSET,eoff);
	DDV_MinMaxDouble(pDX,eoff,0,1);
}
BOOL CTaperedSectionDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos;

	pos = sections->GetHeadPosition();
	while(pos) c_Start.AddString(sections->GetNext(pos).name);
	c_Start.SetCurSel(c_Start.FindString(-1,s_sec->name));

	pos = sections->GetHeadPosition();
	while(pos) c_End.AddString(sections->GetNext(pos).name);
	c_End.SetCurSel(c_End.FindString(-1,e_sec->name));

	c_EIyy.AddString("Linear");
    c_EIyy.AddString("Parabolic");
	c_EIyy.AddString("Cubic");
	c_EIyy.SetCurSel(EIyy);

	c_EIzz.AddString("Linear");
    c_EIzz.AddString("Parabolic");
	c_EIzz.AddString("Cubic");
	c_EIzz.SetCurSel(EIzz);

	return TRUE;
}
void CTaperedSectionDia::OnOK() {
	SECTION sec;
	POSITION pos;

    c_Start.GetLBText(c_Start.GetCurSel(),sec.name);
  	if(pos = sections->Find(sec)) s_sec = &sections->GetAt(pos);

	c_End.GetLBText(c_End.GetCurSel(),sec.name);
  	if(pos = sections->Find(sec)) e_sec = &sections->GetAt(pos);

	EIyy = c_EIyy.GetCurSel();
    EIzz = c_EIzz.GetCurSel();

	CDialog::OnOK();

	if(member) {
		member->section = s_sec;
		member->e_section = e_sec;
		member->EIyy = EIyy;
		member->EIzz = EIzz;
		member->start_offset = soff;
		member->end_offset = eoff;
	}
}
void CTaperedSectionDia::OnEditSection(UINT nID) {
	CSectionPropertyDia* Dia;
	if(nID == IDC_EDITSS)
		Dia = new CSectionPropertyDia(this,s_sec,materials);
	else
		Dia = new CSectionPropertyDia(this,e_sec,materials);
	Dia->DoModal();
	delete Dia;
}
/*
Slab Information Dialog
*/
BEGIN_MESSAGE_MAP(CSlabInfoDia , CDialog)
   	ON_COMMAND(IDC_EDIT_SECTION,OnEditSection)
	ON_COMMAND(IDC_EDIT_AXIS,OnEditAxis)
	ON_COMMAND(IDC_EDIT_DIVISIONS,OnEditDivisions)
END_MESSAGE_MAP()

CSlabInfoDia::CSlabInfoDia(CWnd* parent,SLAB* pslab,PASECTIONLIST psection,
							   PMATERIALLIST pmaterial,PLOADCASELIST pldcase) : 
               CDialog(IDD_SLAB_DATA,parent) {
    slab = pslab;
	section = psection;
	material = pmaterial;
	loadcase = pldcase;
}
void CSlabInfoDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_SECTION_LIST, c_List);
	DDX_Control(pDX, IDC_FORCE_LIST, c_FList);
}
BOOL CSlabInfoDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	POSITION pos;

	/*section*/
	ASECTION* psec;
	pos = section->GetHeadPosition();
	while(pos) {
        psec = &section->GetNext(pos);
		c_List.AddString(psec->name);
	}
	c_List.SetCurSel(c_List.FindString(0,slab->section->name));

	/*load*/
	CString str;
	LOAD* pload;
    pos = slab->load.GetHeadPosition();
	while(pos) {
		pload = &slab->load.GetNext(pos);

		c_FList.AddString(LoadTypeStr[pload->type]);
		str = "      Load case:  " + pload->loadcase->name;
		c_FList.AddString(str);

		if(pload->type == UNIFORM) {
			str = "      Coordinate: ";
			if(pload->system) str += pload->system->name;
			else str += "LOCAL";
			c_FList.AddString(str);
			
			str = "      Direction : " + DirString[pload->dir];
			c_FList.AddString(str);
			str.Format("      W         : %.2f",pload->P);
			c_FList.AddString(str);
		} else {
			if(pload->type == TEMPERATURE) {
				if(pload->dir) {
					str.Format("      DT/DL      : %.2f",pload->P);
					c_FList.AddString(str);
					str = "      Direction : " + DirString[pload->dir];
					c_FList.AddString(str);
				} else {
					str.Format("      DT      : %.2f",pload->P);
					c_FList.AddString(str);
				}
			} else if(pload->type == STRAIN) {
				str.Format("      DL/L      : %.3g",pload->P);
				c_FList.AddString(str);
				str = "      Direction : " + DirString[pload->dir];
				c_FList.AddString(str);
			}
		}
	}

	return TRUE;
}
void CSlabInfoDia::OnEditSection() {
	ASECTION sec;
	POSITION pos;

	c_List.GetLBText(c_List.GetCurSel(),sec.name);
  	if(pos = section->Find(sec)) {
		ASECTION* psec = &section->GetAt(pos);
		sec = *psec;
		CASectionPropertyDia SectionDia(this,&sec,material);
		if(SectionDia.DoModal() == IDOK) {
			*psec = sec;
			UpdateData(FALSE);
		}
	}
}
void CSlabInfoDia::OnOK() {
	ASECTION sec;
	POSITION pos;
    c_List.GetLBText(c_List.GetCurSel(),sec.name);
  	if(pos = section->Find(sec)) {
		slab->section = &section->GetAt(pos);
	}
	CDialog::OnOK();
}
void CSlabInfoDia::OnEditAxis() {
	CAxisRotationDia Dia(this,slab->alpha,"Slab local axis");
	if(Dia.DoModal() == IDOK) {
		slab->alpha = Dia.alpha;
	}
}
void CSlabInfoDia::OnEditDivisions() {
	CSlabDivisionsDia Dia(AfxGetMainWnd());
	Dia.nDivx = slab->nDivx;
	Dia.nDivy = slab->nDivy; 
	Dia.constraint = slab->edge_constraint;
    if(Dia.DoModal() == IDOK) {
		slab->nDivx = Dia.nDivx;
		slab->nDivy = Dia.nDivy;
		slab->edge_constraint = Dia.constraint;
	}
}
/*
Internal Member force diagrams
*/
CForceDiagramDia::CForceDiagramDia(CWnd* parent) 
                    : CDialog(IDD_FORCE_DIAGRAM,parent) {
}
void CForceDiagramDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
    DDX_Radio(pDX,IDC_FX,dir);
	DDX_Check(pDX,IDC_FILLDIAGRAM,fill);
	DDX_Check(pDX,IDC_SHOWVALUE,value);
}
/*
Slab force diagrams
*/
BEGIN_MESSAGE_MAP(CSlabForceDiagramDia , CDialog)
   	ON_CBN_SELCHANGE(IDC_TYPE,OnUpdateControls)
END_MESSAGE_MAP()

CSlabForceDiagramDia::CSlabForceDiagramDia(CWnd* parent) 
                    : CDialog(IDD_SLAB_DIAGRAMS,parent) {
}
void CSlabForceDiagramDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_TYPE,c_Type);
    DDX_Control(pDX,IDC_STRESS,c_Stress);
	DDX_Check(pDX,IDC_AVERAGE,average);
}
BOOL CSlabForceDiagramDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_Type.AddString("Force");
	c_Type.AddString("Stress");

	if(dir >= ISX && dir <= ISVMAX) {
		dir -= ISX;
		c_Type.SetCurSel(1);
	} else if(dir >= IFX) {
		dir -= IFX;
		c_Type.SetCurSel(0);
	} else {
		dir = 0;
		c_Type.SetCurSel(0);
	}
		
	OnUpdateControls();
	return TRUE;
}
void CSlabForceDiagramDia::OnUpdateControls() {

	int type = c_Type.GetCurSel();

	c_Stress.ResetContent();
	if(type == 0) {
		c_Stress.AddString("FX");
		c_Stress.AddString("FY");
		c_Stress.AddString("FXY");
		c_Stress.AddString("FMAX");
		c_Stress.AddString("FMIN");
		c_Stress.AddString("FVM");
		c_Stress.AddString("MX");
		c_Stress.AddString("MY");
		c_Stress.AddString("MXY");
		c_Stress.AddString("MMAX");
		c_Stress.AddString("MMIN");
		c_Stress.AddString("VYZ");
		c_Stress.AddString("VXZ");
		c_Stress.AddString("VMAX");
	} else {
		c_Stress.AddString("SX");
		c_Stress.AddString("SY");
		c_Stress.AddString("SXY");
		c_Stress.AddString("SMAX");
		c_Stress.AddString("SMIN");
		c_Stress.AddString("SVM");
		c_Stress.AddString("SYZ");
		c_Stress.AddString("SZX");
		c_Stress.AddString("SVMAX");
	}
	c_Stress.SetCurSel(dir);
}
void CSlabForceDiagramDia::OnOK() {
    int type = c_Type.GetCurSel();
	dir = c_Stress.GetCurSel();
	if(type == 0) dir += IFX;
	else dir += ISX;
	CDialog::OnOK();
}
/*
Member Release Dialog
*/
BEGIN_MESSAGE_MAP(CMemberReleaseDia,CDialog)
    ON_COMMAND_RANGE(IDC_NUX,IDC_FRZ,OnUpdateControls)
END_MESSAGE_MAP()

CMemberReleaseDia::CMemberReleaseDia(CWnd* parent,MEMBER* member) : 
                  CDialog(IDD_MEMBER_RELEASE,parent) {
    r[0] = (member->nrelease & UX) >> 0;
	r[1] = (member->nrelease & UY) >> 1;
	r[2] = (member->nrelease & UZ) >> 2;
	r[3] = (member->nrelease & RX) >> 3;
	r[4] = (member->nrelease & RY) >> 4;
	r[5] = (member->nrelease & RZ) >> 5;
	r[6] = (member->frelease & UX) >> 0;
	r[7] = (member->frelease & UY) >> 1;
	r[8] = (member->frelease & UZ) >> 2;
	r[9] = (member->frelease & RX) >> 3;
	r[10] = (member->frelease & RY) >> 4;
	r[11] = (member->frelease & RZ) >> 5;
}
void CMemberReleaseDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	for(int i = 0;i < 12;i++) {
		DDX_Control(pDX,IDC_NUX + i,c_R[i]);
		DDX_Check(pDX,IDC_NUX + i,r[i]);
	}
}
void CMemberReleaseDia::OnUpdateControls(UINT nID) {
	int id = nID - IDC_NUX;

	if(c_R[id].GetCheck() == BST_CHECKED) {
		switch(id) {
		case 0: 
			c_R[6].SetCheck(BST_UNCHECKED); 
			break;
		case 1: 
			c_R[7].SetCheck(BST_UNCHECKED);
			if(c_R[5].GetCheck() == BST_CHECKED && c_R[11].GetCheck() == BST_CHECKED)
				c_R[1].SetCheck(BST_UNCHECKED);
			break;
		case 2: 
			c_R[8].SetCheck(BST_UNCHECKED); 
			if(c_R[4].GetCheck() == BST_CHECKED && c_R[10].GetCheck() == BST_CHECKED)
				c_R[2].SetCheck(BST_UNCHECKED);
			break;
		case 3: 
			c_R[9].SetCheck(BST_UNCHECKED); 
			break;
		case 4: 
			if(c_R[2].GetCheck() == BST_CHECKED && c_R[8].GetCheck() == BST_CHECKED)
				c_R[4].SetCheck(BST_UNCHECKED);
			else if(c_R[2].GetCheck() == BST_CHECKED || c_R[8].GetCheck() == BST_CHECKED)
			    c_R[10].SetCheck(BST_UNCHECKED);
			break;
		case 5: 
			if(c_R[1].GetCheck() == BST_CHECKED && c_R[7].GetCheck() == BST_CHECKED)
				c_R[5].SetCheck(BST_UNCHECKED);
			else if(c_R[1].GetCheck() == BST_CHECKED || c_R[7].GetCheck() == BST_CHECKED)
				c_R[11].SetCheck(BST_UNCHECKED);
			break;
		case 6: 
			c_R[0].SetCheck(BST_UNCHECKED); 
			break;
		case 7: 
			c_R[1].SetCheck(BST_UNCHECKED);
			if(c_R[5].GetCheck() == BST_CHECKED && c_R[11].GetCheck() == BST_CHECKED)
				c_R[7].SetCheck(BST_UNCHECKED);
			break;
		case 8: 
			c_R[2].SetCheck(BST_UNCHECKED);
			if(c_R[4].GetCheck() == BST_CHECKED && c_R[10].GetCheck() == BST_CHECKED)
				c_R[8].SetCheck(BST_UNCHECKED);
			break;
		case 9: 
			c_R[3].SetCheck(BST_UNCHECKED);
			break;
		case 10:
			if(c_R[2].GetCheck() == BST_CHECKED && c_R[8].GetCheck() == BST_CHECKED)
				c_R[10].SetCheck(BST_UNCHECKED);
			else if(c_R[2].GetCheck() == BST_CHECKED || c_R[8].GetCheck() == BST_CHECKED)
				c_R[4].SetCheck(BST_UNCHECKED);
			break;
		case 11:
			if(c_R[1].GetCheck() == BST_CHECKED && c_R[7].GetCheck() == BST_CHECKED)
				c_R[11].SetCheck(BST_UNCHECKED);
			else if(c_R[1].GetCheck() == BST_CHECKED || c_R[7].GetCheck() == BST_CHECKED)
				c_R[5].SetCheck(BST_UNCHECKED);
			break; 
		};
	}
}
/*
Member Divisons
*/
CMemberDivisionsDia::CMemberDivisionsDia(CWnd* parent) :
                        CDialog(IDD_MEMBER_DIVISIONS,parent) {
	nMinDiv = 9;
	nMinFrameDiv = 1;
	DivAtInterim = TRUE;
}
void CMemberDivisionsDia::DoDataExchange(CDataExchange* pDX) {
    DDX_Text(pDX,IDC_DIVISIONS,nMinDiv);
	DDV_MinMaxInt(pDX,nMinDiv,3,1000);
	DDX_Text(pDX,IDC_FRAMEDIV,nMinFrameDiv);
	DDV_MinMaxInt(pDX,nMinDiv,1,10000);
	DDX_Check(pDX,IDC_DIVATINTERIM,DivAtInterim);
}
/*
Slab Divisons
*/
CSlabDivisionsDia::CSlabDivisionsDia(CWnd* parent) :
                        CDialog(IDD_SLAB_DIVISIONS,parent) {
	nDivx = 1;
	nDivy = 1;
}
void CSlabDivisionsDia::DoDataExchange(CDataExchange* pDX) {
    DDX_Text(pDX,IDC_DIVX,nDivx);
	DDV_MinMaxInt(pDX,nDivx,1,10000);
	DDX_Text(pDX,IDC_DIVY,nDivy);
	DDV_MinMaxInt(pDX,nDivy,1,10000);
	DDX_Check(pDX,IDC_CONSTRAINT,constraint);
}
/*
Paste Dialog
*/
CPasteDia::CPasteDia(CWnd* parent) :
           CDialog(IDD_PASTE,parent) {
}
void CPasteDia::DoDataExchange(CDataExchange* pDX) {
    DDX_Text(pDX,IDC_X,v.x);
	DDX_Text(pDX,IDC_Y,v.y);
	DDX_Text(pDX,IDC_Z,v.z);
}
/*
Result Dialog
*/
CResultDia::CResultDia(CWnd* parent,UINT N, UINT sel,int ttype) :
           CDialog(IDD_RESULT,parent) {
    Total = N;
	selected = sel;
	type = ttype;
}
BOOL CResultDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
    if(type == NORMAL) {
		for(UINT i = 0;i < Total;i++) {
			CString str;
			str.Format("%d",i);
			c_List.AddString(str);
		}
	} else if(type == MINMAX) {
		c_List.AddString("Max");
        c_List.AddString("Min");
	}
	if(selected >= c_List.GetCount())
		c_List.SetCurSel(0);
	else
		c_List.SetCurSel(selected);
    return TRUE;
}
void CResultDia::DoDataExchange(CDataExchange* pDX) {
	DDX_Control(pDX,IDC_MODE,c_List);
}
void CResultDia::OnOK() {
	selected = c_List.GetCurSel();
	CDialog::OnOK();
}
/*
modified tree control
*/
BEGIN_MESSAGE_MAP(CmyTreeCtrl , CTreeCtrl)
    ON_WM_LBUTTONDOWN()
END_MESSAGE_MAP()

void CmyTreeCtrl::SelectSiblings(HTREEITEM hItem,BOOL checked) {
	HTREEITEM hChildItem = GetChildItem(hItem);
	while (hChildItem != NULL) {
		SetCheck(hChildItem,checked);
		SelectSiblings(hChildItem,checked);
		hChildItem = GetNextSiblingItem(hChildItem);
	}
}
void CmyTreeCtrl::SelectParents(HTREEITEM hItem,BOOL checked) {
	HTREEITEM hParentItem = GetParentItem(hItem);
	while (hParentItem != NULL) {
		SetCheck(hParentItem,checked);
		SelectParents(hParentItem,checked);
		hParentItem = GetParentItem(hParentItem);
	}
}
void CmyTreeCtrl::OnLButtonDown(UINT nFlags,CPoint point) {
	CTreeCtrl::OnLButtonDown(nFlags,point);
	UINT uFlags;
	HTREEITEM hItem = HitTest(point,&uFlags);
	if ((hItem != NULL) && (TVHT_ONITEM & uFlags)) {
		Select(hItem, TVGN_CARET);
		BOOL checked = !GetCheck(hItem);
		if(TVHT_ONITEMSTATEICON & uFlags)  checked = !checked;

		SetCheck(hItem,checked);
		SelectSiblings(hItem,checked);

		if(checked) {
			SelectParents(hItem,checked);
		} else {
			BOOL flag = FALSE;
			HTREEITEM hParentItem = GetParentItem(hItem);
			HTREEITEM hChildItem = GetChildItem(hParentItem);
			while (hChildItem != NULL) {
				flag = GetCheck(hChildItem);
				if(flag) break;
				hChildItem = GetNextSiblingItem(hChildItem);
			}
			if(!flag)
				SetCheck(hParentItem,FALSE);
		}
	}
}
/*
Select table dialog
*/

CSelectTableDia::CSelectTableDia(CWnd* parent,PANALYSISCASELIST* pcases) :
              CDialog(IDD_SELECT_TABLE,parent) {
	analysis_cases = pcases;
}
void CSelectTableDia::DoDataExchange(CDataExchange* pDX) {
	DDX_Control(pDX,IDC_TABLE,c_Table);
	DDX_Control(pDX,IDC_SELECT_CASES,c_Cases);
}
BOOL CSelectTableDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	
	c_Table.ModifyStyle( 0,(TVS_HASBUTTONS | TVS_HASLINES | TVS_LINESATROOT | TVS_CHECKBOXES));
	
	TABLE* item;
	TV_INSERTSTRUCT  curItem;
	HTREEITEM hItem;
	UINT i;
	int found;

	for (i = 0; i < NTABLES; i++) {
		m_rghItem[i] = 0;
	}

	do {
		found = false;
		for (i = 0; i < NTABLES; i++) {
			if(m_rghItem[i]) continue;
			item = &all_tables[i];
			if(item->view > TABLE::current_view)
				continue;
			if(item->parent == -1) {
				curItem.hParent = NULL;
			} else {
				if(!m_rghItem[item->parent]) continue;
				curItem.hParent = m_rghItem[item->parent];
			}
			curItem.hInsertAfter = TVI_LAST;
			curItem.item.pszText = item->name;
			curItem.item.mask = TVIF_TEXT;
			hItem = c_Table.InsertItem(&curItem);
			if(item->save) c_Table.SetCheck(hItem);
			m_rghItem[i] = hItem;
			found = true;
		}
	} while(found);

	c_Table.Expand(m_rghItem[TABLE::Model],TVE_EXPAND);
	if(TABLE::current_view >= TABLE::ANALYSIS) 
		c_Table.Expand(m_rghItem[TABLE::Analysis],TVE_EXPAND);
    if(TABLE::current_view >= TABLE::DESIGN) 
		c_Table.Expand(m_rghItem[TABLE::Design],TVE_EXPAND);

	c_Cases.ShowWindow(TABLE::current_view >= TABLE::ANALYSIS);
	if(TABLE::current_view >= TABLE::ANALYSIS) {
		POSITION pos;
		ANALYSISCASE* pcase;
		PANALYSISCASELIST pcaselist;
		for(i = 0;i < CASETYPES;i++) {
			pcaselist = analysis_cases[i];
			if(pcaselist) {
				pos = pcaselist->GetHeadPosition();
				while(pos) {
					pcase = &pcaselist->GetNext(pos);
					if(pcase->finished) {
						c_Cases.SetSel(c_Cases.AddString(pcase->name));
					}
				}
			}
		}
	}
    return TRUE;
}
void CSelectTableDia::OnOK() {
	CDialog::OnOK();
	UINT i;
	TABLE* item;
    for (i = 0; i < NTABLES; i ++) {
		if(!m_rghItem[i]) continue;
		item = &all_tables[i];
		item->save = c_Table.GetCheck(m_rghItem[i]);
	}

	if(TABLE::current_view >= TABLE::ANALYSIS) {
		POSITION pos;
		ANALYSISCASE* pcase;
		PANALYSISCASELIST pcaselist;
		for(i = 0;i < CASETYPES;i++) {
			pcaselist = analysis_cases[i];
			if(pcaselist) {
				pos = pcaselist->GetHeadPosition();
				while(pos) {
					pcase = &pcaselist->GetNext(pos);
					if(pcase->finished) {
						pcase->selected = c_Cases.GetSel(c_Cases.FindString(-1,pcase->name));
					}
				}
			}
		}
	}
}
/*
Display tables dialog
*/
BEGIN_MESSAGE_MAP(CTableViewDia , CDialog)
	ON_CBN_SELCHANGE(IDC_LIST,OnTableChange)
	ON_NOTIFY(LVN_ITEMCHANGED,IDC_TABLE,OnSelectChange)
END_MESSAGE_MAP()

CTableViewDia::CTableViewDia(CWnd* parent,CmyDatabase* myd) :
              CDialog(IDD_TABLE_VIEW,parent),
			  c_Table(CGridCtrl::TYPE1) {
	pdatabase = myd;
}
void CTableViewDia::DoDataExchange(CDataExchange* pDX) {
	DDX_Control(pDX,IDC_TABLE,c_Table);
	DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_MAX,c_Max);
	DDX_Control(pDX,IDC_INDEX,c_Index);
}
BOOL CTableViewDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	CmyData* ptable;
	POSITION pos = pdatabase->database.GetHeadPosition();
	while(pos) {
		ptable = &pdatabase->database.GetNext(pos);
		c_List.AddString(ptable->name);
	}
	c_List.SetCurSel(0);

	OnTableChange();
    return TRUE;
}
void CTableViewDia::OnTableChange() {
	int i;

	/*delete table*/
	c_Table.DeleteAllItems();
	int nColumnCount = c_Table.GetHeaderCtrl()->GetItemCount();
	for (i=0;i < nColumnCount;i++) {
		c_Table.DeleteColumn(0);
	}

	/*new table*/
	int index = c_List.GetCurSel();
	CmyData* ptable = &pdatabase->database.GetAt(pdatabase->database.FindIndex(index));

	CRect rect;
	c_Table.GetClientRect(&rect);
	int nInterval;
	if(ptable->nFields > 8) nInterval = (rect.Width() - 15) / 8;
	else nInterval = (rect.Width() - 15) / ptable->nFields;

	CString tokens[50];
	CmyData::tokenize(ptable->field_names,tokens," ");
	for(i = 0;i < ptable->nFields;i++) {
         c_Table.InsertColumn(i,tokens[i],LVCFMT_RIGHT,nInterval);
	}

    LV_ITEM lvItem;
	lvItem.mask = LVIF_TEXT | LVIF_STATE;
	lvItem.state = 0;
	lvItem.stateMask = 0;
	lvItem.iSubItem = 0;

	CString* pstr;
	POSITION pos = ptable->data.GetHeadPosition();
	i = 0;
	while(pos) {
		pstr = &ptable->data.GetNext(pos);
		CmyData::tokenize(*pstr,tokens,"\n");
		lvItem.iItem = i;
		lvItem.pszText = tokens[0].GetBuffer(0);
		c_Table.InsertItem(&lvItem);
		for(int j = 1;j < ptable->nFields;j++) {
			c_Table.SetItemText(i,j,tokens[j]);
		}
		i++;
	}

	c_Index.SetWindowText("");
	c_Max.SetWindowText("");
	c_Table.SetItemState(0, LVIS_SELECTED, LVIS_SELECTED | LVIS_FOCUSED);
}
void CTableViewDia::OnSelectChange(LPNMHDR pnmhdr, LRESULT* /*pResult*/) {
	CString str;
	POSITION pos = c_Table.GetFirstSelectedItemPosition();
	if(pos) {
		str.Format("%d",c_Table.GetNextSelectedItem(pos) + 1);
	    c_Index.SetWindowText(str);
	}
	str.Format("%d",c_Table.GetItemCount());
	c_Max.SetWindowText(str);
}
/*
Select Labels Dialog
*/
BEGIN_MESSAGE_MAP(CSelectLabelsDia , CDialog)
	ON_LBN_SELCHANGE(IDC_TYPE,OnType)
    ON_COMMAND_RANGE(IDC_SELECT,IDC_DESELECT,OnSelect)
END_MESSAGE_MAP()

CSelectLabelsDia::CSelectLabelsDia(CWnd* parent,BOOL tselect) :
           CDialog(IDD_SELECT_LABELS,parent) {
	select = tselect;
	minid = 1;
	maxid = 1000;
	inc = 2;
}
BOOL CSelectLabelsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Type.AddString("Joints");
    c_Type.AddString("Members");
    c_Type.AddString("Slabs");
	c_Type.SetCurSel(0);
	OnType();
    return TRUE;
}
void CSelectLabelsDia::DoDataExchange(CDataExchange* pDX) {
	DDX_Text(pDX,IDC_MIN,minid);
    DDX_Text(pDX,IDC_MAX,maxid);
	DDV_MinMaxInt(pDX,maxid,minid,10000);
	DDX_Text(pDX,IDC_INC,inc);
	DDV_MinMaxInt(pDX,inc,1,10000);
	DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_LIST,c_List);
	DDX_Control(pDX,IDC_MIN,c_Min);
	DDX_Control(pDX,IDC_MAX,c_Max);
	DDX_Control(pDX,IDC_INC,c_Inc);
}
void CSelectLabelsDia::OnSelect(UINT nID) {
	int minid,maxid,inc,index;
    CString str;
	c_Min.GetWindowText(str);
	minid = atoi(str);
    c_Max.GetWindowText(str);
	maxid = atoi(str);
	c_Inc.GetWindowText(str);
	inc = atoi(str);

	for(int i = minid;i <= maxid;i += inc) {
		str.Format("%d",i);
		index = c_List.FindString(-1,str);
		if(index != LB_ERR) {
			if(nID == IDC_SELECT) c_List.SetSel(index,TRUE);
			else c_List.SetSel(index,FALSE);
		}
	}
}
void CSelectLabelsDia::OnOK() {
	CmyDocument* pDoc = GetMyDocument();
	int type = c_Type.GetCurSel();
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	POSITION pos;

	if(type == JOINTT) {
		pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			if(c_List.GetSel(c_List.FindString(-1,joint->name)))
				joint->sel = select;
		}
	} else if(type == MEMBERT) {
		pos = pDoc->members.GetHeadPosition();
		while(pos) {
			member = &pDoc->members.GetNext(pos);
			if(c_List.GetSel(c_List.FindString(-1,member->name)))
				member->sel = select;
		}
	} else {
        pos = pDoc->slabs.GetHeadPosition();
		while(pos) {
			slab = &pDoc->slabs.GetNext(pos);
			if(c_List.GetSel(c_List.FindString(-1,slab->name)))
				slab->sel = select;
		} 
	}

	CDialog::OnOK();
}
void CSelectLabelsDia::OnType() {
	CmyDocument* pDoc = GetMyDocument();
    int type = c_Type.GetCurSel();
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	POSITION pos;

	c_List.ResetContent();

	if(type == JOINTT) {
		pos = pDoc->joints.GetHeadPosition();
		while(pos) {
			joint = &pDoc->joints.GetNext(pos);
			c_List.AddString(joint->name);
		}
		c_List.SetCurSel(0);
	} else if(type == MEMBERT) {
		pos = pDoc->members.GetHeadPosition();
		while(pos) {
			member = &pDoc->members.GetNext(pos);
			c_List.AddString(member->name);
		}
		c_List.SetCurSel(0);
	} else {
		pos = pDoc->slabs.GetHeadPosition();
		while(pos) {
			slab = &pDoc->slabs.GetNext(pos);
			c_List.AddString(slab->name);
		}
		c_List.SetCurSel(0);
	}
}
/*
Analysis Options Dialog
*/
BEGIN_MESSAGE_MAP(CAnalysisOptionsDia , CDialog)
	ON_COMMAND_RANGE(IDC_PLANETRUSS,IDC_PLANEGRID,OnFastDof)
END_MESSAGE_MAP()

CAnalysisOptionsDia::CAnalysisOptionsDia(CWnd* parent,UBMP8* psdofs,UBMP8* pddofs) : 
                  CDialog(IDD_ANALYSIS_OPTIONS,parent) {
	sdofs = psdofs;
	ddofs = pddofs;
	register int i;
	for(i = IUX;i <= IRZ;i++) {
		if(*sdofs & (1 << i)) sdof[i] = TRUE;
		else sdof[i] = FALSE;
	}
	for(i = IUX;i <= IUZ;i++) {
		if(*ddofs & (1 << i)) ddof[i] = TRUE;
		else ddof[i] = FALSE;
	}
}
void CAnalysisOptionsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	register int i;
	for(i = IUX;i <= IRZ;i++)
	    DDX_Check(pDX,IDC_NUX + i,sdof[i]);
	for(i = IUX;i <= IUZ;i++)
	    DDX_Check(pDX,IDC_FUX + i,ddof[i]);
}
void CAnalysisOptionsDia::OnOK() {
	CDialog::OnOK();
	*sdofs = NOUR;
	*ddofs = NOUR;
	register int i;
	for(i = IUX;i <= IRZ;i++)
		if(sdof[i]) *sdofs |= (1 << i);
	for(i = IUX;i <= IUZ;i++)
		if(ddof[i]) *ddofs |= (1 << i);
}
BOOL CAnalysisOptionsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
    return TRUE;
}
void CAnalysisOptionsDia::OnFastDof(UINT nID) {
	UBMP8 dofs;
	if(nID == IDC_PLANETRUSS) dofs = UX | UZ;
	else if(nID == IDC_SPACETRUSS) dofs = ALLU;
	else if(nID == IDC_PLANEFRAME) dofs = UX | UZ | RY;
    else if(nID == IDC_SPACEFRAME) dofs = ALLUR;
	else if(nID == IDC_PLANEGRID) dofs = UZ | RX | RY;
	else dofs = NOUR;

	for(int i = IUX;i <= IRZ;i++) {
		if(dofs & (1 << i)) sdof[i] = TRUE;
		else sdof[i] = FALSE;
	}
	UpdateData(FALSE);
}
/*
Design Options Dialog
*/
BEGIN_MESSAGE_MAP(CDesignOptionsDia , CDialog)
    ON_COMMAND(IDC_ADD,OnAdd)
	ON_COMMAND(IDC_DELETE,OnDelete)
	ON_CBN_SELCHANGE(IDC_TYPE,OnType)
END_MESSAGE_MAP()

CDesignOptionsDia::CDesignOptionsDia(CWnd* parent,DETAILING* tpdetail) : 
                  CDialog(IDD_DESIGN_OPTIONS,parent),
				  c_Options(CGridCtrl::TYPE3) {
	pdetail = tpdetail;
	detail = *pdetail;
}
void CDesignOptionsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX,IDC_SINGLE_SHEAR,detail.single_shear);
	DDX_Check(pDX,IDC_CUT,detail.cut_at_L3);
	DDX_Text(pDX,IDC_SPACING,detail.min_spacing);
	DDX_CBIndex(pDX,IDC_LAYERS,detail.max_layers);
	DDX_Control(pDX,IDC_LAYERS,c_Layers);
	DDX_Control(pDX,IDC_BAR_LIST,c_BarList);
	DDX_Control(pDX,IDC_SHEAR_LIST,c_List[SHEART]);
	DDX_Control(pDX,IDC_BEAM_LIST,c_List[BEAMT]);
	DDX_Control(pDX,IDC_COLUMN_LIST,c_List[COLUMNT]);
	DDX_Control(pDX,IDC_SLAB_LIST,c_List[SLABT]);
    DDX_Control(pDX,IDC_TYPE,c_Type);
	DDX_Control(pDX,IDC_OPTIONS,c_Options);
}
void CDesignOptionsDia::OnOK() {
	c_Options.GetData();
	CDialog::OnOK();
	*pdetail = detail;
}
BOOL CDesignOptionsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	c_Layers.AddString("Single");
	c_Layers.AddString("Double");
	c_Layers.AddString("Triple");
	c_Layers.SetCurSel(0);

	DOUBLE v;
	POSITION pos;
	CString str;

	pos = detail.shearbarlist.GetHeadPosition();
	while(pos) {
		v = detail.shearbarlist.GetNext(pos);
		str.Format("%d",int(v * 1000));
		c_List[SHEART].AddString(str);
	}
	pos = detail.beambarlist.GetHeadPosition();
	while(pos) {
		v = detail.beambarlist.GetNext(pos);
		str.Format("%d",int(v * 1000));
		c_List[BEAMT].AddString(str);
	}
	pos = detail.columnbarlist.GetHeadPosition();
	while(pos) {
		v = detail.columnbarlist.GetNext(pos);
		str.Format("%d",int(v * 1000));
		c_List[COLUMNT].AddString(str);
	}
	pos = detail.slabbarlist.GetHeadPosition();
	while(pos) {
		v = detail.slabbarlist.GetNext(pos);
		str.Format("%d",int(v * 1000));
		c_List[SLABT].AddString(str);
	}

	
	c_Type.AddString("Shear");
	c_Type.AddString("Beam");
	c_Type.AddString("Column");
	c_Type.AddString("Slab");
	c_Type.SetCurSel(SHEART);

	
	/*options*/
	CString length_str[] = {"Don't use factor","Use program calculated","0"};
	CString frame_str[] = {"Moment Frame","Braced frame","0"};
	c_Options.SetTwoHeaders();
	c_Options.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&detail.fs_concrete,"Concrete FS");
	c_Options.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&detail.fs_steel,"Steel FS");
	c_Options.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TINDEX,&detail.frame_type,"Framing",frame_str);
	c_Options.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TINDEX,&detail.le_option,"Effective Length",length_str);
	c_Options.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TBOOL,&detail.accidental_eccentricity,"Add accidental eccentricity");
	c_Options.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TBOOL,&detail.individual_buckling,"Add column buckling effect");
	c_Options.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TINT,&detail.nWDIV,"Number of interaction curves");
	c_Options.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TINT,&detail.nQDIV,"Interaction chart division");

	OnType();
    return TRUE;
}
void CDesignOptionsDia::OnType() {
	CString str;
	int type = c_Type.GetCurSel();
	int index,i;

	c_BarList.ResetContent();
	for(i = 6;i <= 40;i += 2) {
		str.Format("%d",i);
		index = c_List[type].FindString(-1,str);
		if(index == LB_ERR)
		   c_BarList.AddString(str);
	}
	for(i = 0;i < 4;i++) {
		c_List[i].EnableWindow(type == i);
	}
}
void CDesignOptionsDia::OnAdd() {
	CString str;
	DOUBLE v;
	int type = c_Type.GetCurSel(),i;
	int ncount = c_BarList.GetSelCount();

	BARSIZELIST* plist;
	if(type == SHEART) plist = &detail.shearbarlist;
	else if(type == BEAMT) plist = &detail.beambarlist;
	else if(type == COLUMNT) plist = &detail.columnbarlist;
	else if(type == SLABT) plist = &detail.slabbarlist;

	LPINT selected = new int[ncount];
	c_BarList.GetSelItems(ncount,selected);

	for(i = 0;i < ncount;i++) {
		c_BarList.GetText(selected[i],str);
		v = atof(str) / 1000.0;
		plist->AddTail(v);
		c_List[type].AddString(str);
	}
	for(i = ncount - 1;i >= 0;i--) {
		c_BarList.DeleteString(selected[i]);
	}
	delete[] selected;
}
void CDesignOptionsDia::OnDelete() {
	CString str;
	DOUBLE v;
	int type = c_Type.GetCurSel(),i;
	int ncount = c_List[type].GetSelCount();
	POSITION pos;

	BARSIZELIST* plist;
	if(type == SHEART) plist = &detail.shearbarlist;
	else if(type == BEAMT) plist = &detail.beambarlist;
	else if(type == COLUMNT) plist = &detail.columnbarlist;
	else if(type == SLABT) plist = &detail.slabbarlist;

	LPINT selected = new int[ncount];
	c_List[type].GetSelItems(ncount,selected);
	for(i = 0;i < ncount;i++) {
		c_List[type].GetText(selected[i],str);
		v = atof(str) / 1000.0;
		pos = plist->Find(v);
		plist->RemoveAt(pos);
		c_BarList.AddString(str);
	}
	for(i = ncount - 1;i >= 0;i--) {
		c_List[type].DeleteString(selected[i]);
	}
	delete[] selected;
}
/*
Misc Options Dialog
*/
BEGIN_MESSAGE_MAP(CMiscOptionsDia , CDialog)
END_MESSAGE_MAP()

CMiscOptionsDia::CMiscOptionsDia(CWnd* parent) : 
                  CDialog(IDD_MISC_OPTIONS,parent) {
	CmyDocument* pDoc = GetMyDocument();
	AnalysisResult = pDoc->AnalysisResult;
	scale_deflection = pDoc->user_disp_scale;
    scale_stress = pDoc->user_force_scale;
	scale_load = pDoc->user_load_scale;
}
void CMiscOptionsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX,IDC_DEFLECTION,scale_deflection);
    DDX_Text(pDX,IDC_FORCE,scale_stress);
	DDX_Text(pDX,IDC_LOAD,scale_load);
}
void CMiscOptionsDia::OnOK() {
	CDialog::OnOK();
	CmyDocument* pDoc = GetMyDocument();
	pDoc->user_disp_scale = scale_deflection;
    pDoc->user_force_scale = scale_stress;
	pDoc->user_load_scale = scale_load;
}
BOOL CMiscOptionsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
    return TRUE;
}
/*
Stress Diagram Dialog
*/
BEGIN_MESSAGE_MAP(CStressControl , CStatic)
	ON_WM_LBUTTONDOWN()
	ON_WM_MOUSEMOVE()
END_MESSAGE_MAP()

void CStressControl::UpdateParent(CPoint& point) {
	parent->section = (point.x - CStressDiagramDia::OFFSET) / parent->xscale;
	parent->scroll = parent->SCROLL;
	parent->c_Scroll.SetCurSel(parent->scroll);
	CString str;
	str.Format("%.3f",parent->section);
	parent->c_Section.SetWindowText(str);
}
void CStressControl::OnLButtonDown(UINT nFlags,CPoint point) {
	CStatic::OnLButtonDown(nFlags,point);
	UpdateParent(point);
}
void CStressControl::OnMouseMove(UINT nFlags,CPoint point) {
	CStatic::OnMouseMove(nFlags,point);
	int OFFSET = CStressDiagramDia::OFFSET;
	if(nFlags == MK_LBUTTON) {
		if(point.x < OFFSET) point.x = OFFSET;
		else if(point.x > parent->Width - OFFSET) point.x = parent->Width - OFFSET;
	    UpdateParent(point);
	}
}

const int CStressDiagramDia::OFFSET = 20;

BEGIN_MESSAGE_MAP(CStressDiagramDia,CDialog)
    ON_WM_PAINT()
	ON_LBN_SELCHANGE(IDC_ITEM_DIRECTION,OnChangeDir)
	ON_LBN_SELCHANGE(IDC_SCROLL,OnScroll)
	ON_EN_UPDATE(IDC_SECTION,OnSection)
	ON_WM_TIMER()
	ON_COMMAND_RANGE(IDC_ABSOLUTE,IDC_RELATIVE,OnDeflectionType)
END_MESSAGE_MAP()

CStressDiagramDia::CStressDiagramDia(CWnd* parent,MEMBER* tmember,ANALYSISCASE* tcase,CmyView* pView) : 
                  CTimerUpdateDia(parent,IDD_STRESS_DIAGRAM,tmember,pView) {
	member = tmember;
	analysis_case = tcase;
	dir = pView->force_diagram;
	if(dir == IUZ || dir == IRY) dir = MAJOR;
	else if(dir == IUY || dir == IRZ) dir = MINOR;
	else if(dir == IUX || dir == IRX) dir = AXIAL;
	else dir = MAJOR;
	section = 0;
	scroll = SCROLL;
	deflection_type = DABSOLUTE;
}
BOOL CStressDiagramDia::OnInitDialog() {
	if(!CTimerUpdateDia::OnInitDialog())
		return FALSE;
   	c_Dir.AddString("Major(MY & VZ)");
	c_Dir.AddString("Minor(MZ & VY)");
	c_Dir.AddString("Axial(MX & VX)");
	c_Dir.SetCurSel(dir);

	c_Scroll.AddString("Scroll for values");
	c_Scroll.AddString("Show maximum");
	c_Scroll.AddString("Show minimum");
	c_Scroll.SetCurSel(scroll);

    c_Drawing[0].parent = this;
	c_Drawing[1].parent = this;
	c_Drawing[2].parent = this;
	
	OnChangeDir();
    return TRUE;
}
void CStressDiagramDia::DoDataExchange(CDataExchange* pDX) {
	CTimerUpdateDia::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_ITEM_DIRECTION,c_Dir);
	DDX_Control(pDX,IDC_MOMENT,c_Drawing[0]);
	DDX_Control(pDX,IDC_SHEAR,c_Drawing[1]);
	DDX_Control(pDX,IDC_DEFLECTION,c_Drawing[2]);
	DDX_Control(pDX,IDC_MOMENT_VALUE,c_Value[0]);
	DDX_Control(pDX,IDC_SHEAR_VALUE,c_Value[1]);
	DDX_Control(pDX,IDC_DEFLECTION_VALUE,c_Value[2]);
    DDX_Control(pDX,IDC_SECTION,c_Section);
	DDX_Control(pDX,IDC_SCROLL,c_Scroll);
	DDX_Text(pDX,IDC_NAME,analysis_case->name);
    DDX_Text(pDX,IDC_SECTION,section);
	DDX_Radio(pDX,IDC_ABSOLUTE,deflection_type);
}
void CStressDiagramDia::OnChangeDir() {
	dir = c_Dir.GetCurSel();
	c_Drawing[2].EnableWindow(dir != AXIAL);
    c_Value[2].ShowWindow(dir != AXIAL);
	UpdateDrawing();
}
void CStressDiagramDia::OnScroll() {
	scroll = c_Scroll.GetCurSel();
	c_Section.ShowWindow(scroll == SCROLL);
	if(scroll == SCROLL) {
		CString str;
		c_Section.GetWindowText(str);
		section = atof(str);
	}
	UpdateDrawing();
}
void CStressDiagramDia::OnSection() {
	CString str;
	c_Section.GetWindowText(str);
	section = atof(str);
	UpdateDrawing();
}
void CStressDiagramDia::DrawRectangle(CDC* pDC,DOUBLE x1,DOUBLE h1,DOUBLE x2,DOUBLE h2) {
	x1 = (xscale * x1) + OFFSET;
	x2 = (xscale * x2) + OFFSET;
	h1 = (yscale * h1);
	h2 = (yscale * h2);
#define T(x,y) CPoint(int(x),int(y))
	CPoint pt[4];
	pt[0] = T(x1,Height / 2);
	pt[1] = T(x1,h1 + Height / 2);
	pt[2] = T(x2,h2 + Height / 2);
	pt[3] = T(x2,Height / 2);
#undef T
	pDC->Polygon(pt,4);
}
void CStressDiagramDia::DrawDiagram(int index,int force_diagram,BOOL paint) {
	
	/*paint drawing*/
	CRect rect;
	CStatic* pDrawing = &c_Drawing[index];
	CStatic* pText = &c_Value[index];
	CDC* pDC;
	CPaintDC paintDC(pDrawing);
	if(paint) pDC = &paintDC;
	else pDC = pDrawing->GetDC();

	CBrush wbrush(RGB(0,0,0));
	CBrush* pOldBrush = pDC->SelectObject(&wbrush);
	pDrawing->GetClientRect(&rect);
	pDC->Rectangle(rect);
	pDC->SelectObject(pOldBrush);
	Height = rect.Height();
	Width = rect.Width();
	
	/*clip*/
	CRgn rgn;
	rgn.CreateRectRgn(0,0,Width,Height);
	pDC->SelectClipRgn(&rgn);

	DOUBLE value;
	CPen* pOldPen;

	CPen pen1(PS_SOLID,1,RGB(0,255,0));
	pOldPen = pDC->SelectObject(&pen1);
	pDC->MoveTo(0,Height / 2);
	pDC->LineTo(Width,Height / 2);
	pDC->SelectObject(pOldPen);

	/*deflection*/
	DOUBLE L = member->GetLength();
	if(index == 2) {
		xscale = (Width - 2 * OFFSET) / L;
		yscale = 0;

		
		UINT i;
		DOUBLE * v,slope,y1,y2;
		UINT NDIV = member->nDiv;
		y1 = member->CalculateDeflection(0,force_diagram);
        y2 = member->CalculateDeflection(L,force_diagram);
		slope = (y2 - y1) / L;
		
		v = new DOUBLE[NDIV];
		for(i = 0;i < NDIV;i++) {
			v[i] = member->CalculateDeflection(member->station[i],force_diagram);
			if(deflection_type == DRELATIVE) v[i] -= (y1 + slope * member->station[i]);
			if(fabs(v[i]) <= 1e-8) v[i] = 0;
			if(fabs(v[i]) > yscale) yscale = fabs(v[i]);
		}
		if(yscale) yscale = (Height / 2) / yscale;
		else yscale = 1;
		
		CPoint p1;
		CPen pen(PS_SOLID,2,RGB(255,0,0));
		pOldPen = pDC->SelectObject(&pen);
		for(i = 0;i < NDIV;i++) {
			p1 = CPoint(int(xscale * member->station[i] + OFFSET),int(-yscale * v[i] + Height / 2));
			if(!i) pDC->MoveTo(p1);
			else pDC->LineTo(p1);
		}
		pDC->SelectObject(pOldPen);

		/*value*/
		if(scroll == SCROLL) {
			value = member->CalculateDeflection(section,force_diagram);
			if(deflection_type == DRELATIVE) value -= (y1 + slope * section);
			if(fabs(value) <= 1e-8) value = 0;
		} else if(scroll == SHOWMAX) {
			value = -1e16;
			for(i = 0;i < NDIV;i++) {
				if(v[i] > value) {
					section = member->station[i];
					value = v[i];
				}
			}
		} else if(scroll == SHOWMIN) {
			value = 1e16;
			for(i = 0;i < NDIV;i++) {
				if(v[i] < value) {
					section = member->station[i];
					value = v[i];
				}
			}
		}
		delete[] v;
	} else {
		/*force diagram*/
		DOUBLE x1,h1,x2,h2,hmax = 0,xmax = 0;
		CBrush BrushPositive(CmyView::PositiveColor);
		CBrush BrushNegative(CmyView::NegativeColor);
		CPen PenPositive(PS_SOLID,1,CmyView::PositiveColor);
		CPen PenNegative(PS_SOLID,1,CmyView::NegativeColor);

		CBrush* pOldBrush = pDC->SelectObject(&BrushPositive);
		pOldPen = pDC->SelectObject(&PenPositive);
		
		/*detrmine scale*/
		UBMP32 i;
		for(i = 0;i < member->nDiv;i++) {
			x1 = member->station[i];
			h1 = fabs(member->forces[i][force_diagram]);
			if(x1 > xmax) xmax = x1;
			if(h1 > hmax) hmax = h1;
		}
		xscale = (Width - 2 * OFFSET) / L;
		if(hmax) yscale = (Height / 2) / hmax;
		else yscale = 1;
		
		/*plot stress diagram*/
		for(i = 0;i < member->nDiv;i++) {
			x1 = member->station[i];
			h1 = -member->forces[i][force_diagram];
			
			if(i) {
				if(h1 * h2 >= 0) {
					if(h1 < 0 || (h1 == 0 && h2 < 0) ) {
						pDC->SelectObject(&BrushNegative);
						pDC->SelectObject(&PenNegative);
					} else if(h1 > 0 || (h1 == 0 && h2 > 0)) {
						pDC->SelectObject(&BrushPositive);
						pDC->SelectObject(&PenPositive);
					}
					DrawRectangle(pDC,x1,h1,x2,h2);
				} else {
					DOUBLE x3 = x2 + (x1 - x2) * (h2 / (h2 - h1));
					DOUBLE h3 = 0;
					
					pDC->SelectObject(&BrushPositive);
					pDC->SelectObject(&PenPositive);
					if(h1 < 0) DrawRectangle(pDC,x2,h2,x3,h3);
					else  DrawRectangle(pDC,x1,h1,x3,h3);
					
					pDC->SelectObject(&BrushNegative);
					pDC->SelectObject(&PenNegative);
					if(h1 > 0) DrawRectangle(pDC,x2,h2,x3,h3);
					else  DrawRectangle(pDC,x1,h1,x3,h3);
				}
			}
			x2 = x1;
			h2 = h1;
		}
		
		pDC->SelectObject(pOldBrush);
		pDC->SelectObject(pOldPen);
				
		/*value*/
		if(scroll == SCROLL) {
			value = member->CalculateForce(section,force_diagram);
		} else if(scroll == SHOWMAX) {
			value = -1e16;
			for(i = 0;i < member->nDiv;i++) {
				if(member->forces[i][force_diagram] > value) {
					section = member->station[i];
					value = member->forces[i][force_diagram];
				}
			}
		} else if(scroll == SHOWMIN) {
			value = 1e16;
			for(i = 0;i < member->nDiv;i++) {
				if(member->forces[i][force_diagram] < value) {
					section = member->station[i];
					value = member->forces[i][force_diagram];
				}
			}
		}
	}

	/*section*/
	CPen pen(PS_SOLID,1,RGB(0,255,0));
	pOldPen = pDC->SelectObject(&pen);
	pDC->MoveTo(int(section * xscale + OFFSET),0);
	pDC->LineTo(int(section * xscale + OFFSET),Height);
	pDC->SelectObject(pOldPen);

	/*show values*/
	CString str;
	if(index == 2) str.Format("%.6f \nat %.3f",value,section);
	else str.Format("%.4f \nat %.3f",value,section);
	pText->SetWindowText(str);

	/*releae dc*/
	pDrawing->ReleaseDC(pDC);
}
void CStressDiagramDia::UpdateDrawing(BOOL paint) {
    if(dir == MAJOR) {
		DrawDiagram(0,IRY,paint);
		DrawDiagram(1,IUZ,paint);
		DrawDiagram(2,IUZ,paint);
	} else if(dir == MINOR) {
		DrawDiagram(0,IRZ,paint);
		DrawDiagram(1,IUY,paint);
		DrawDiagram(2,IUY,paint);
	} else {
		DrawDiagram(0,IRX,paint);
		DrawDiagram(1,IUX,paint);
	}
}
void CStressDiagramDia::OnPaint() {
	CTimerUpdateDia::OnPaint();
	UpdateDrawing(TRUE);
}
void CStressDiagramDia::OnDeflectionType(UINT nID) {
	deflection_type = nID - IDC_ABSOLUTE;
	UpdateDrawing();
}
/*
Design Overwrites Dialog
*/
BEGIN_MESSAGE_MAP(CDesignOverwritesDia,CDialog)
END_MESSAGE_MAP()

CDesignOverwritesDia::CDesignOverwritesDia(CWnd* parent,MEMBER* tmember) : 
                  CDialog(IDD_DESIGN_OVERWRITES,parent),
				  c_Table(CGridCtrl::TYPE2) {
	member = tmember;
}
BOOL CDesignOverwritesDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;

	CString str;
	c_Table.SetTwoHeaders();

	if(member->section->rebar.design == BEAM &&
		member->section->material->type == MATERIAL::CONCRETE
		) {
		c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&member->redistribution_factor,"Moment redistribution factor");
	} else {
		CString length_str[] = {"Don't use factor","Use program calculated","Use Assigned","Use common preference","0"};
		c_Table.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TINDEX,&member->le_option,"Effective Length",length_str);
		for(UINT i = 0;i < 2;i++) {
			str.Format("Calculated Effective L%c",'y' + i);
			c_Table.InsertMyItem(GRIDCTRLDATA::NONE,GRIDCTRLDATA::TDBL,&member->LeCalculated[i],str);
			str.Format("Used Effective L%c",'y' + i);
			c_Table.InsertMyItem(GRIDCTRLDATA::EDIT,GRIDCTRLDATA::TDBL,&member->Le[i],str);
		}
		CString frame_str[] = {"Sway","Non sway","Use common preference","0"};
		c_Table.InsertMyItem(GRIDCTRLDATA::COMBO,GRIDCTRLDATA::TINDEX,&member->frame_type,"Framing",frame_str);
	}

    return TRUE;
}
void CDesignOverwritesDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_TABLE,c_Table);
}
void CDesignOverwritesDia::OnOK() {
	c_Table.GetData();
	CDialog::OnOK();
}
/*
Design Result Dialog
*/
BEGIN_MESSAGE_MAP(CDesignResultDia,CDialog)
    ON_WM_TIMER()
    ON_COMMAND(IDC_OVERWRITES,OnDesignOverwrites)
END_MESSAGE_MAP()

CDesignResultDia::CDesignResultDia(CWnd* parent,MEMBER* tmember,ANALYSISCASE* tcase,CmyView* pView) : 
                  CTimerUpdateDia(parent,IDD_DESIGN_RESULT,tmember,pView),
				  c_Table(CGridCtrl::TYPE1) {
	member = tmember;
	analysis_case = tcase;
}
BOOL CDesignResultDia::OnInitDialog() {
	if(!CTimerUpdateDia::OnInitDialog())
		return FALSE;

	if(member->design_data->flags & DESIGNDATA::SKIPPED)
		return TRUE;

	/*new table*/
	DESIGNDATA* pdata;
	CString str;
	UINT i,j,nt;
	int design;

	if(member->section->rebar.design != COULMN) {
		nt = 1;
	} else if(member->section->rebar.type == RECTANGULAR) {
		nt = 2 * (member->section->rebar.nz + member->section->rebar.ny - 2);
	} else {
		nt = member->section->rebar.nt;
	}
	design = member->section->rebar.design;

	CRect rect;
	c_Table.GetClientRect(&rect);
	int nInterval;

	
    LV_ITEM lvItem;
	lvItem.mask = LVIF_TEXT | LVIF_STATE;
	lvItem.state = 0;
	lvItem.stateMask = 0;
	lvItem.iSubItem = 0;

    if(member->section->material->type == MATERIAL::STEEL) {
		nInterval = (rect.Width() - 15) / 9;
		c_Table.InsertColumn(0,"Location",LVCFMT_RIGHT,nInterval);
		
		c_Table.InsertColumn(1,"Axial",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(2,"Major",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(3,"Minor",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(4,"Total",LVCFMT_RIGHT,nInterval);

		c_Table.InsertColumn(5,"LT Axial",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(6,"LT Major",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(7,"LT Minor",LVCFMT_RIGHT,nInterval);
		c_Table.InsertColumn(8,"LT Total",LVCFMT_RIGHT,nInterval);
		
		for(i = 0;i < member->nDiv;i++) {
			pdata = &member->design_data[i];
			
			str.Format("%.2f",member->station[i]);
			lvItem.iItem = i;
			lvItem.pszText = str.GetBuffer(0);
			c_Table.InsertItem(&lvItem);
			
			for(j = 0;j < 8;j++) {
				str.Format("%.6f",pdata->ratios[j]);
				c_Table.SetItemText(i,1 + j,str);
			}
		}
	} else {
		if(design == COULMN)  {
			nInterval = (rect.Width() - 15) / 3;
			c_Table.InsertColumn(0,"Location",LVCFMT_RIGHT,nInterval);
			c_Table.InsertColumn(1,"Total",LVCFMT_RIGHT,nInterval);
			c_Table.InsertColumn(2,"Shear",LVCFMT_RIGHT,nInterval);
		} else {
			nInterval = (rect.Width() - 15) / 4;
			c_Table.InsertColumn(0,"Location",LVCFMT_RIGHT,nInterval);
			c_Table.InsertColumn(1,"Top",LVCFMT_RIGHT,nInterval);
			c_Table.InsertColumn(2,"Bottom",LVCFMT_RIGHT,nInterval);
			c_Table.InsertColumn(3,"Shear",LVCFMT_RIGHT,nInterval);
		}
		
		for(i = 0;i < member->nDiv;i++) {
			pdata = &member->design_data[i];

			str.Format("%.2f",member->station[i]);
			lvItem.iItem = i;
			lvItem.pszText = str.GetBuffer(0);
			c_Table.InsertItem(&lvItem);

			if(!pdata->flags) {
				for(j = 0;j < 2;j++) {
					if(!pdata->total[j])
						continue;
					if(pdata->total[j] > 1) {
						if(nt == 1) {
							str.Format("%3dD%-2d and %3dD%-2d [%-4dmm2]  ",
								pdata->count[j][0],int(pdata->diam[j][0] * 1000),
								pdata->count[j][1],int(pdata->diam[j][1] * 1000),
								int(pdata->Area[j] * 1e6));
						} else {
							str.Format("%d X (%3dD%-2d and %3dD%-2d) [%-4dmm2]  ",
								nt,
								pdata->count[j][0],int(pdata->diam[j][0] * 1000),
								pdata->count[j][1],int(pdata->diam[j][1] * 1000),
								int(pdata->Area[j] * 1e6));
						}
					} else {
						str.Format("%3dD%-2d        [%-4dmm2]  ",
							nt * pdata->count[j][0],int(pdata->diam[j][0] * 1000),
							int(pdata->Area[j] * 1e6));
					}
					c_Table.SetItemText(i,1 + j,str);
				}
				str.Format("       D%-2dc/c%4d [%-4dmm2]",int(pdata->sd * 1000),pdata->s,int(pdata->As * 1e6));
				c_Table.SetItemText(i,(design == COULMN) ? 2 : 3,str);
			} else {
				c_Table.SetItemText(i,1,"Insufficient");
			}
		}
	}

	c_Table.SetItemState(0, LVIS_SELECTED, LVIS_SELECTED | LVIS_FOCUSED);
    return TRUE;
}
void CDesignResultDia::DoDataExchange(CDataExchange* pDX) {
	CTimerUpdateDia::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_LIST,c_Table);
	DDX_Text(pDX,IDC_MEMBER,member->name);
    DDX_Text(pDX,IDC_SECTION,member->section->name);
	DDX_Text(pDX,IDC_NAME,analysis_case->name);
}
void CDesignResultDia::OnDesignOverwrites() {
	CDesignOverwritesDia Dia(this,member);
	Dia.DoModal();
}
/*
Undo / Redo
*/
CUndoRedo::CUndoRedo() {
	count = 0;
	current = 0;
}
void CUndoRedo::Clear() {
	count = 0;
	current = 0;
	data.RemoveAll();
}
void CUndoRedo::Add(UndoData& udata) {
	for(int i = 0;i < count - current;i++)
		data.RemoveTail();
	count = current;
	if(count >= 100) {
		data.RemoveHead();
		current = --count;
	}
	data.AddTail(udata);
	current = ++count;
}
void CUndoRedo::Save(JOINT* pjoint,BOOL add) {
	UndoData udata;
	udata.name = pjoint->name;
	udata.type = UndoData::TJOINT;
	udata.joints.AddTail(*pjoint);
	udata.adding = add;
	Add(udata);
}
void CUndoRedo::Save(MEMBER* pmem,BOOL add) {
	UndoData udata;
	udata.name = pmem->name;
	udata.type = UndoData::TMEMBER;
	udata.joints.AddTail(*pmem->j1);
    udata.joints.AddTail(*pmem->j2);
	udata.adding = add;
	Add(udata);
}
void CUndoRedo::Save(SLAB* pslab,BOOL add) {
	UndoData udata;
	udata.name = pslab->name;
	udata.type = UndoData::TSLAB;
	for(UINT i = 0;i < pslab->NJ;i++)
		udata.joints.AddTail(*pslab->jt[i]);
	udata.adding = add;
	Add(udata);
}
void CUndoRedo::Do(CmyDocument* pDoc,int redo) {
	if(redo) {
		if(current == count) return;
	} else {
		if(current == 0) return;
		current--;
	}

    UndoData* pdata;
	pdata = &data.GetAt(data.FindIndex(current));
	if(pdata->type == UndoData::TJOINT) {
		JOINT* j;
		j = &pdata->joints.GetHead();
        if(pdata->adding && !redo) pDoc->DeleteJoint(j);
		else pDoc->AddJoint(j);
	} else if(pdata->type == UndoData::TMEMBER) {
		MEMBER mem;
		mem.j1 = &pdata->joints.GetHead();
		mem.j2 = &pdata->joints.GetTail();
        if(pdata->adding && !redo) pDoc->DeleteMember(&mem);
		else pDoc->AddMember(mem.j1,mem.j2);
	} else if(pdata->type == UndoData::TSLAB) {
		SLAB slab;
		POSITION pos = pdata->joints.GetHeadPosition();
		for(UINT i = 0;i < slab.NJ;i++)
			slab.jt[i] = &pdata->joints.GetNext(pos);
        if(pdata->adding && !redo) pDoc->DeleteSlab(&slab);
		else pDoc->AddSlab(slab.jt,4);
	}

	if(redo) {
		current++;
	}
}
/*
Document
*/
UINT CmyDocument::nViews = 2;
CmyDocument* CmyDocument::Current_Document = 0;

IMPLEMENT_DYNCREATE(CmyDocument,CDocument)

BEGIN_MESSAGE_MAP(CmyDocument,CDocument)
    ON_COMMAND(IDM_RUN_ANALYSIS,OnRunAnalysis)
	ON_COMMAND(IDM_ANALYSIS_OPTIONS,OnAnalysisOptions)
    ON_COMMAND(IDM_MISC_OPTIONS,OnMiscOptions)
    ON_COMMAND_RANGE(IDM_DESIGN,IDM_DESIGN_SELECTED,OnDesign)
	ON_COMMAND(IDM_DESIGN_OPTIONS,OnDesignOptions)
    ON_COMMAND(IDM_DESIGN_VERIFY,OnDesignVerify)  
	ON_COMMAND_RANGE(IDM_PREPARE_DRAWING,IDM_PREPARE_DRAWING_SELECTED,OnPrepareDrawing)
	ON_COMMAND(IDM_JOINT_RESTRAINT,OnJointRestraint)
	ON_COMMAND_RANGE(IDM_JOINT_CONSTRAINT,IDM_JOINT_CONSTRAINTZ,OnJointConstraint)
    ON_COMMAND_RANGE(IDM_JOINT_FORCE,IDM_JOINT_MASS,OnJointLoad)
	ON_COMMAND(IDM_JOINT_AXIS,OnJointAxis)
	ON_COMMAND(IDM_MEMBER_POINT,OnMemberPointLoad)
    ON_COMMAND(IDM_MEMBER_DISTRIBUTED,OnMemberDistLoad)
	ON_COMMAND(IDM_MEMBER_AXIS,OnMemberAxis)
	ON_COMMAND(IDM_MEMBER_DIVISIONS,OnMemberDivisions)
	ON_COMMAND(IDM_SLAB_DISTRIBUTED,OnSlabDistLoad)
	ON_COMMAND(IDM_SLAB_AXIS,OnSlabAxis)
	ON_COMMAND(IDM_SLAB_REVERSE_AXIS,OnReverseSlabAxis)
	ON_COMMAND(IDM_SLAB_DIVISIONS,OnSlabDivisions)
	ON_COMMAND(IDM_DEFINE_MATERIAL,OnDefineMaterial)
	ON_COMMAND(IDM_DEFINE_FRAME_SECTION,OnDefineFrameSection)
    ON_COMMAND(IDM_DEFINE_SLAB_SECTION,OnDefineSlabSection)
	ON_COMMAND(IDM_LOAD_CASE,OnDefineLoadCase)
	ON_COMMAND(IDM_MODAL_CASE,OnDefineModalCase)
	ON_COMMAND(IDM_RESPONSE_CASE,OnDefineResponseHistory)
	ON_COMMAND(IDM_RESPONSE_SPECTRUM,OnDefineResponseSpectrum)
    ON_COMMAND(IDM_LOAD_COMBO,OnDefineLoadCombo)
	ON_COMMAND(IDM_NONLINEAR_CASE,OnDefineNonLinearCase)
	ON_COMMAND(IDM_BUCKLING,OnDefineBucklingCase)
	ON_COMMAND(IDM_CONSTRAINT,OnDefineConstraint)
	ON_COMMAND(IDM_TIME_HISTORY_FUNC,OnDefineTimeHistoryFunction)
	ON_COMMAND(IDM_GROUP,OnDefineGroup)
	ON_COMMAND(IDM_SPECTRUM_FUNC,OnDefineSpectrumFunction)
	ON_COMMAND(IDM_SELECT_ALL_JOINTS,OnSelectAllJoints)
	ON_COMMAND(IDM_SELECT_ALL_MEMBERS,OnSelectAllMembers)
	ON_COMMAND(IDM_SELECT_ALL_SLABS,OnSelectAllSlabs)
	ON_COMMAND(IDM_SELECT_ALL,OnSelectAll)
	ON_COMMAND(IDM_DESELECT_ALL_JOINTS,OnDeselectAllJoints)
	ON_COMMAND(IDM_DESELECT_ALL_MEMBERS,OnDeselectAllMembers)
	ON_COMMAND(IDM_DESELECT_ALL_SLABS,OnDeselectAllSlabs)
	ON_COMMAND(IDM_DESELECT_ALL,OnDeselectAll)
	ON_COMMAND_RANGE(IDM_SELECT_LABELS,IDM_DESELECT_LABELS,OnSelectLabels)
	ON_COMMAND_RANGE(IDM_SELECT_SECTION,IDM_DESELECT_SECTION,OnSelectSection)
	ON_COMMAND_RANGE(IDM_SELECT_ASECTION,IDM_DESELECT_ASECTION,OnSelectSlabSection)
	ON_COMMAND_RANGE(IDM_SELECT_CONSTRAINT,IDM_DESELECT_CONSTRAINT,OnSelectConstraint)
	ON_COMMAND_RANGE(IDM_SELECT_GROUP,IDM_DESELECT_GROUP,OnSelectGroup)
	ON_COMMAND(IDM_INVERT_SELECTION,OnInvertSelection)
	ON_COMMAND(IDM_SELECT_FROM_ALL,OnSelectFromAll)
    ON_UPDATE_COMMAND_UI(IDM_SELECT_FROM_ALL,OnUpdateSelectFromAllUI)
	ON_COMMAND(IDM_MEMBER_SECTION,OnAssignSection)
	ON_COMMAND(IDM_TAPERED_SECTION,OnAssignTaperedSection)
	ON_COMMAND(IDM_SLAB_SECTION,OnAssignSlabSection)
    ON_COMMAND(IDM_MEMBER_JOIN,OnJoinMembers)
	ON_COMMAND(IDM_MEMBER_BREAK,OnBreakMember)
	ON_COMMAND(IDM_MEMBER_RELEASE,OnMemberRelease)
	ON_COMMAND(IDM_SLAB_JOIN,OnJoinSlabs)
	ON_COMMAND(IDM_SLAB_BREAK,OnBreakSlab)
	ON_COMMAND(IDM_ASSIGN_GROUP,OnAssignGroup)
	ON_COMMAND(IDM_DESIGN_COMBOS,OnAddDesignCombos)
    ON_UPDATE_COMMAND_UI(IDM_DESIGN_COMBOS,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_DESELECT_ALL_JOINTS,IDM_DESELECT_ALL,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_JOINT_RESTRAINT,IDM_JOINT_AXIS,OnUpdateAssignJointUI)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_MEMBER_SECTION,IDM_MEMBER_RELEASE,OnUpdateAssignMemberUI)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_SLAB_SECTION,IDM_SLAB_STRAIN,OnUpdateAssignSlabUI)
	ON_UPDATE_COMMAND_UI(IDM_TAPERED_SECTION,OnUpdateAssignMemberUI)
	ON_UPDATE_COMMAND_UI(IDM_MEMBER_DIVISIONS,OnUpdateAssignMemberUI)
	ON_CBN_SELCHANGE(IDM_SELECT_LOADCASE,OnSelectLoadCase) 
	ON_COMMAND(ID_DOCUMENT_LOCK,OnLock)
    ON_UPDATE_COMMAND_UI(ID_DOCUMENT_LOCK,OnUpdateLockUI)
	ON_COMMAND(IDM_SHOW_RESULT,OnShowResult)
	ON_COMMAND(IDM_SHOW_DESIGN_RESULT,OnShowDesignResult)
	ON_UPDATE_COMMAND_UI(IDM_SHOW_RESULT,OnUpdateAnalysisResultUI)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_DESIGN,IDM_DESIGN_SELECTED,OnUpdateDesignUI)
	ON_UPDATE_COMMAND_UI(IDM_SHOW_DESIGN_RESULT,OnUpdateDesignResultUI)
	ON_UPDATE_COMMAND_UI(IDM_DESIGN_VERIFY,OnUpdateDesignResultUI)
	ON_COMMAND(IDM_MEMBER_STRAIN,OnMemberStrain)
    ON_COMMAND(IDM_SLAB_STRAIN,OnSlabStrain)
	ON_COMMAND(IDM_SYSTEM,OnDefineSystem)
	ON_COMMAND_RANGE(IDM_COPY,IDM_CUT,OnCopy)
	ON_COMMAND(IDM_PASTE,OnPaste)
	ON_COMMAND(IDM_DELETE,OnDelete)
    ON_COMMAND(ID_DELETE,OnDelete)
	ON_COMMAND(ID_ENTER,OnEnter)
	ON_COMMAND(IDM_COLORS,OnColors)
	ON_COMMAND(IDM_PLOT_FUNCTIONS,OnPlotFunctions)
	ON_COMMAND_RANGE(ID_SNAP_TO_GRID,ID_SNAP_TO_MEMBER,OnSnapTo)
	ON_COMMAND(IDM_SHOW_TABLES,OnShowTables)
    ON_COMMAND_RANGE(IDM_UNDO,IDM_REDO,OnUndoRedo)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_COPY,IDM_CUT,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DELETE,OnUpdateDeselectUI)
    ON_UPDATE_COMMAND_UI(ID_DELETE,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DESELECT_LABELS,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DESELECT_SECTION,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DESELECT_ASECTION,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DESELECT_CONSTRAINT,OnUpdateDeselectUI)
	ON_UPDATE_COMMAND_UI(IDM_DESELECT_GROUP,OnUpdateDeselectUI)
    ON_UPDATE_COMMAND_UI_RANGE(IDM_UNDO,IDM_REDO,OnUpdateUndoRedoUI)
    ON_UPDATE_COMMAND_UI_RANGE(ID_SNAP_TO_GRID,ID_SNAP_TO_MEMBER,OnUpdateSnapToUI)
	ON_UPDATE_COMMAND_UI(IDM_PLOT_FUNCTIONS,OnUpdatePlotFunctionsUI)
	ON_UPDATE_COMMAND_UI(IDM_GROUP,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_ASSIGN_GROUP,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_PASTE,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_DEFINE_MATERIAL,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_DEFINE_FRAME_SECTION,OnUpdateLockMenuUI)
    ON_UPDATE_COMMAND_UI(IDM_DEFINE_SLAB_SECTION,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_SYSTEM,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_LOAD_CASE,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_MODAL_CASE,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_RESPONSE_CASE,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_RESPONSE_SPECTRUM,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_LOAD_COMBO,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_CONSTRAINT,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_TIME_HISTORY_FUNC,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_SPECTRUM_FUNC,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_NONLINEAR_CASE,OnUpdateLockMenuUI)
	ON_UPDATE_COMMAND_UI(IDM_BUCKLING,OnUpdateLockMenuUI)
	ON_COMMAND_RANGE(IDM_EXPORT_EXCEL,IDM_EXPORT_AUTOCAD,OnExport)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_PREPARE_DRAWING,IDM_PREPARE_DRAWING_SELECTED,OnUpdateDesignResultUI)
END_MESSAGE_MAP()

CmyDocument::CmyDocument() {
	SYSTEM sys;
	systems.AddTail(sys);
	global = &systems.GetHead();
	LoadCaseListBox = &(((CMainFrame*)AfxGetMainWnd())->myToolBar.c_List);
	Lock = FALSE;
	all_analysis_cases[0] = (PANALYSISCASELIST)&loadcases;
	all_analysis_cases[1] = (PANALYSISCASELIST)&modalcases;
	all_analysis_cases[2] = (PANALYSISCASELIST)&responsecases;
	all_analysis_cases[3] = (PANALYSISCASELIST)&responsespecs;
	all_analysis_cases[4] = (PANALYSISCASELIST)&combinations;
	all_analysis_cases[5] = (PANALYSISCASELIST)&nlcases;
	all_analysis_cases[6] = (PANALYSISCASELIST)&bucklingcases;
	entities[0] = (PENTITYLIST)&joints;
	entities[1] = (PENTITYLIST)&members;
	entities[2] = (PENTITYLIST)&slabs;
	pAnalysisProgressDia = NULL;
	SECTION::pDocDetailing = &detailing;
	Current_Document = this;
}

CmyDocument::~CmyDocument() {
}

BOOL CmyDocument::OnNewDocument() {
	if (!CDocument::OnNewDocument())
		return FALSE;

	UINT i;
	/*coordinate system*/
	SYSTEM sys;
	sys.name = "GLOBAL";
    sys.coordinate = CARTESIAN;
	sys.rank = LEVEL3;
	systems.AddTail(sys);
	global = &systems.GetHead();

	/*initialize*/
	UndoRedo.Clear();
	Snap_To_Grid = TRUE;
	Snap_To_Joint = TRUE;
	Snap_To_MemberMid = FALSE;
    Snap_To_Member = FALSE;
	select_all = FALSE;
	state = SELECT;
	nClicks = 0;
	Lock = FALSE;
	AnalysisResult = FALSE;
	DesignResult = FALSE;
	DisplacementResult = FALSE;
	JOINT::TotalNumber = 0;
	MEMBER::TotalNumber = 0;
	SLAB::TotalNumber = 0;
	static_dofs = ALLUR;
	dynamic_dofs = ALLU;
	user_force_scale = 1;
	user_load_scale = 1;
	user_disp_scale = 1;

	/*colors*/
	CmyView::ResetColors();
	
	/*materials*/
	MATERIAL conc,steel,brick,*psteel;

	conc.name = "CONC";
	conc.E = 24855578;
	conc.nu = 0.2;
	conc.unitweight = 23.5631;
	conc.density = 2.4028;
	conc.alphac = 9.9e-6;
	conc.G = conc.E / (2 * (1 + conc.nu)) ;
	conc.rank = LEVEL3;
	conc.type = MATERIAL::CONCRETE;
	conc.fck = 27579.032;
	conc.fctk = 17000;
	conc.fyk = 413685.8;
	conc.fyks = 275790.32;
	conc.ftk = 413685.8;
	conc.ftks = 275790.32;
	materials.AddTail(conc);

	steel.name = "STEEL";
	steel.E = 200e+6;
	steel.nu = 0.3;
	steel.unitweight = 76.9729;
	steel.density = 7.849;
	steel.alphac = 11.7e-6;
	steel.G = steel.E / (2 * (1 + steel.nu)) ;
	steel.rank = LEVEL3;
	steel.type = MATERIAL::STEEL;
	steel.fyk = 248211.28;
	steel.ftk = 399896;
	materials.AddTail(steel);
	psteel = &materials.GetTail();

	SECTION::defmaterial = &materials.GetHead();
	ASECTION::defmaterial = &materials.GetHead();
	
	/*sections*/
	SECTION coulmn,beam;
	coulmn.rank = LEVEL3;
	coulmn.name = "COULMN";
	coulmn.w = 0.4;
	coulmn.h = 0.4;
	coulmn.rebar.design = COULMN;
	sections.AddTail(coulmn);
	MEMBER::defsection = &sections.GetHead();

	beam.rank = LEVEL3;
	beam.name = "BEAM";
	beam.rebar.design = BEAM;
	beam.w = 0.3;
	beam.h = 0.4;
	beam.GetData();
	sections.AddTail(beam);

	/*area sections*/
	ASECTION slab;
	slab.rank = LEVEL3;
	slab.name = "SLAB";
	slab.h = 0.25;
	slab.rebar.design = ASLAB;
	asections.AddTail(slab);
	SLAB::defsection = &asections.GetTail();

	/*cases*/
	LOADCASE dead;
	dead.name = "DEAD";
	dead.swm = 1;
	dead.type = LOADCASE::DEAD;
	dead.rank = LEVEL3;
	loadcases.AddTail(dead);

	MODALCASE modal;
	modal.name = "MODAL";
	modal.rank = LEVEL3;
	modal.run = FALSE;
	modalcases.AddTail(modal);
	RESPONSE::defmodalcase = &modalcases.GetHead();

	FUNCTION unifrh;
	unifrh.name = "UNIFRH";
	unifrh.rank = LEVEL3;
	unifrh.points.AddTail(RPoint(0,1,0));
	unifrh.points.AddTail(RPoint(1,1,0));
	rhfunctions.AddTail(unifrh);
	RESPONSE::deffunction = &rhfunctions.GetHead();

	FUNCTION unifrs;
	unifrs.name = "UNIFRS";
	unifrs.rank = LEVEL3;
	unifrs.points.AddTail(RPoint(0,1,0));
	unifrs.points.AddTail(RPoint(1,1,0));
	rsfunctions.AddTail(unifrs);

	FUNCTION::defspectrumtype = FALSE;

	/*rebars*/
	DOUBLE diam;
	for(i = 6;i <= 8;i += 2) {
		diam = i * 0.001;
		detailing.shearbarlist.AddTail(diam);
	}
    for(i = 12;i <= 32;i += 2) {
		diam = i * 0.001;
		detailing.beambarlist.AddTail(diam);
	}
	for(i = 8;i <= 14;i += 2) {
		diam = i * 0.001;
		detailing.slabbarlist.AddTail(diam);
	}
	for(i = 12;i <= 32;i += 2) {
		diam = i * 0.001;
		detailing.columnbarlist.AddTail(diam);
	}
	/*New model dialog*/
	CmyApplication* pApp = (CmyApplication*) AfxGetApp();
	if(!pApp->skip_model) {
		CNewDia ModelDia(AfxGetMainWnd());
		if(ModelDia.DoModal() == IDOK) {
			CmyView::CREATION[1] = VIEW_XY;
			CmyView::CREATION[2] = VIEW_XZ;
			
			if(ModelDia.type == MODELBLANK) {
			} else if(ModelDia.type == MODELGRID) {
				CNewSystemSht Dia(AfxGetMainWnd(),global);
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODELBEAM) {
				CNewBeamDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODEL2DTRUSS) {
				CNewTruss2DDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODEL3DTRUSS) {
				CNewTruss3DDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODEL2DFRAME || ModelDia.type == MODEL3DFRAME) {
				CNewFrameDia Dia(AfxGetMainWnd(),ModelDia.type);
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODELSLAB) {
				CNewSlabDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODELWALL) {
				CNewWallDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else if(ModelDia.type == MODELSHELL) {
				CNewShellDia Dia(AfxGetMainWnd());
				if(Dia.DoModal() != IDOK)
					return FALSE;
			} else {
				return FALSE;
			}
		} else {
			return FALSE;
		}
	}
	/*fill cases*/
	FillCases();
	/*update views*/
	UpdateDrawing();
	return TRUE;
}

BOOL CmyDocument::OnSaveDocument(LPCTSTR lpszPathName) {
	if (!CDocument::OnSaveDocument(lpszPathName))
      return FALSE;

	/*view data*/
	CFile viewfile;
    CString sPath = lpszPathName;
	int len = sPath.GetLength();
	sPath.SetAt(len - 3,'V');
	sPath.SetAt(len - 2,'I');
	sPath.SetAt(len - 1,'W');
	
	if(viewfile.Open(sPath, CFile::modeCreate | CFile::modeWrite)) {
		CArchive ar(&viewfile, CArchive::store);
		CSerializer sr;
		sr.Set(ar,TRUE);

		CmyView* pView;
		POSITION pos = GetFirstViewPosition();
		while(pos) {
			pView = (CmyView*) GetNextView(pos);
			pView->Serialize(sr);
		}

		ar.Close();
		viewfile.Close();
	}
	/*end*/
	return TRUE;
}
void CmyDocument::InitViews(LPCSTR lpszPathName) {
	/*view data*/
	CFile viewfile;
	char path[256];
	strcpy(path,lpszPathName);
	int len = strlen(path);
	path[len - 3] = 'V';
	path[len - 2] = 'I';
	path[len - 1] = 'W';
	
	if(viewfile.Open(path, CFile::modeRead)) {
		CArchive ar(&viewfile, CArchive::load);
		CSerializer sr;
		sr.Set(ar,FALSE);
		
		CmyView* pView;
		POSITION pos = GetFirstViewPosition();
		while(pos) {
			pView = (CmyView*) GetNextView(pos);
			pView->Serialize(sr);
		}
		
		ar.Close();
		viewfile.Close();
	}
	/*end*/
}
void CmyDocument::UpdateDrawing() {
	CmyView* pView;
	POSITION pos = GetFirstViewPosition();
	while(pos) {
		pView = (CmyView*) GetNextView(pos);
		pView->DrawingHasChanged = TRUE;
		pView->Invalidate();
	}
}
void CmyDocument::OnColors() {
	CColorsDia Dia(AfxGetMainWnd());
    if(Dia.DoModal() == IDOK) {
		UpdateDrawing();
	}
}
/*
Undo/redo
*/
void CmyDocument::OnUndoRedo(UINT nID) {
	if(nID == IDM_UNDO) {
        UndoRedo.Do(this,0);
	} else {
		UndoRedo.Do(this,1);
	}
	UpdateDrawing();
}
void CmyDocument::OnUpdateUndoRedoUI(CCmdUI* pCmdUI) {
	if(pCmdUI->m_nID == IDM_UNDO) pCmdUI->Enable(UndoRedo.current > 0);
	else pCmdUI->Enable(UndoRedo.current < UndoRedo.count);
}
void CmyDocument::DeleteContents() {
	/*if analysis result was on*/
	if(Lock && AnalysisResult) {
		FreeAnalysis();
	}
	/*initialize*/
	UndoRedo.Clear();
	Snap_To_Grid = TRUE;
	Snap_To_Joint = TRUE;
	Snap_To_MemberMid = FALSE;
	Snap_To_Member = FALSE;
	select_all = FALSE;
	state = SELECT;
	nClicks = 0;
	Lock = FALSE;
	AnalysisResult = FALSE;
	DesignResult = FALSE;
	DisplacementResult = FALSE;
	JOINT::TotalNumber = 0;
	MEMBER::TotalNumber = 0;
	SLAB::TotalNumber = 0;
	static_dofs = ALLUR;
	dynamic_dofs = ALLU;
	user_force_scale = 1;
	user_load_scale = 1;
	user_disp_scale = 1;

	/*colors*/
	CmyView::ResetColors();

    /*clear lists*/
	members.RemoveAll();
	slabs.RemoveAll();
	joints.RemoveAll();
	constraints.RemoveAll();
	systems.RemoveAll();
	materials.RemoveAll();
    sections.RemoveAll();
	asections.RemoveAll();
	groups.RemoveAll();
	loadcases.RemoveAll();
	modalcases.RemoveAll();
	rhfunctions.RemoveAll();
    rsfunctions.RemoveAll();
	responsecases.RemoveAll();
	responsespecs.RemoveAll();
	combinations.RemoveAll();
	nlcases.RemoveAll();
	bucklingcases.RemoveAll();
	copyJoints.RemoveAll();
	copyMembers.RemoveAll();
	copySlabs.RemoveAll();
	detailing.shearbarlist.RemoveAll();
	detailing.beambarlist.RemoveAll();
	detailing.columnbarlist.RemoveAll();
	detailing.slabbarlist.RemoveAll();
}
/*
Add / Delete
*/
JOINT* CmyDocument::AddJoint(JOINT* j) {
    JOINT* joint;
	BOOL duplicate = FALSE;
	POSITION pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		if(*joint == *j) {
			duplicate = TRUE;
			break;
		}
	}
	if(duplicate)
		return joint;

	JOINT myj;
	joint = &myj;
	joint->p = j->p;
	joint->name.Format("%d",++JOINT::TotalNumber);
	joints.AddTail(*joint);
	joint = &joints.GetTail();
	return joint;
}
void CmyDocument::DeleteJoint(JOINT* j) {
	POSITION pos;
	if(!(j->mconnect + j->sconnect) && (pos = joints.Find(*j))) {
		if(AnalysisResult) FreeAnalysis(j);
		joints.RemoveAt(pos);
	}
}
MEMBER* CmyDocument::AddMember(JOINT* j1,JOINT* j2) {
	MEMBER* mem;
	MEMBER tmem;

	if(j1->p == j2->p)
		return NULL;

	/*order joints*/
	if(j1->p > j2->p) {
		JOINT* temp = j2;
		j2 = j1;
		j1 = temp;
	}
	/*end*/

	tmem.j1 = j1;
	tmem.j2 = j2;

	BOOL duplicate = FALSE;
	POSITION pos = members.GetHeadPosition();
	while(pos) {
		mem = &members.GetNext(pos);
		if(*mem == tmem) {
			duplicate = TRUE;
			break;
		}
	}
	if(duplicate)
		return mem;

	MEMBER mymem;
	mem = &mymem;

	SECTION sec;
	if(j1->p.z == j2->p.z) sec.name = "BEAM";
	else sec.name = "COULMN";
	if(pos = sections.Find(sec)) mem->section = &sections.GetAt(pos);
	if(!mem->section) mem->section = &sections.GetHead();
	mem->j1 = AddJoint(j1);
	mem->j2 = AddJoint(j2);
	mem->j1->mconnect++;
	mem->j2->mconnect++;
	mem->name.Format("%d",++MEMBER::TotalNumber);
	members.AddTail(*mem);
	mem = &members.GetTail();
	return mem;
}
void CmyDocument::DeleteMember(MEMBER* member) {
	POSITION pos;
	pos = members.Find(*member);
	if(pos) {
		member->j1->mconnect--;
        member->j2->mconnect--;
		DeleteJoint(member->j1);
		DeleteJoint(member->j2);
		if(AnalysisResult) FreeAnalysis(member);
		members.RemoveAt(pos);
	}
}
SLAB* CmyDocument::AddSlab(SLAB* slab) {
	JOINT* Jt[NPOLY];
	for(UINT i = 0;i < slab->NJ;i++) {
		Jt[i] = slab->jt[i];
	}
	return AddSlab(Jt,slab->NJ);
}
SLAB* CmyDocument::AddSlab(JOINT* J,UINT N) {
	JOINT* Jt[NPOLY];
	for(UINT i = 0;i < N;i++) {
		Jt[i] = &J[i];
	}
	return AddSlab(Jt,N);
}
SLAB* CmyDocument::AddSlab(JOINT* j1,JOINT* j2,JOINT* j3,JOINT* j4) {
	JOINT* J[4];
	J[0] = j1;
	J[1] = j2;
	J[2] = j3;
	J[3] = j4;
	return AddSlab(J,4);
}
SLAB* CmyDocument::AddSlab(JOINT** J,UINT N) {
	UINT i,j;
	SLAB* sla;
	SLAB tsla;

	for(i = 0;i < N;i++) { 
		for(j = i + 1;j < N;j++) {
		    if(J[i]->p == J[j]->p) return NULL;
		}
	}

	for(i = 0;i < N;i++) { 
		tsla.jt[i] = J[i];
	}

	BOOL duplicate = FALSE;
	POSITION pos = slabs.GetHeadPosition();
	while(pos) {
		sla = &slabs.GetNext(pos);
		if(*sla == tsla) {
			duplicate = TRUE;
			break;
		}
	}
	if(duplicate)
		return sla;

	if(clockwise(J[0]->p,J[1]->p,J[2]->p)) {
		JOINT* temp;
		for(i = 0;i < N;i++) {
			j = N - 1 - i;
			if(i >= j) continue;
			temp = J[i];
			J[i] = J[j];
		    J[j] = temp;
		}
	}
	SLAB mysla;
	sla = &mysla;

	ASECTION sec;
	sec.name = "SLAB";

	if(pos = asections.Find(sec)) sla->section = &asections.GetAt(pos);
	if(!sla->section) sla->section = &asections.GetHead();

		
	sla->NJ = N;
	for(i = 0;i < N;i++) { 
		sla->jt[i] = AddJoint(J[i]);
		sla->jt[i]->sconnect++;
	}
	
	sla->name.Format("%d",++SLAB::TotalNumber);
	slabs.AddTail(*sla);
	sla = &slabs.GetTail();
	return sla;
}
void CmyDocument::DeleteSlab(SLAB* slab) {
	POSITION pos;
	pos = slabs.Find(*slab);
	if(pos) {
		for(UINT i = 0;i < slab->NJ;i++) {
			slab->jt[i]->sconnect--;
			DeleteJoint(slab->jt[i]);
		}
		if(AnalysisResult) FreeAnalysis(slab);
		slabs.RemoveAt(pos);
	}
}
/*
Joint/Member manuplation
*/
void CmyDocument::OnJointRestraint() {
	JOINT tjoint;
	tjoint.restraint = UX | UY | UZ;
	CJointRestraintDia JointDia(AfxGetMainWnd(),&tjoint);
    if(JointDia.DoModal() == IDOK) {
		UBMP8 restraint = 0;
		for(int i = IUX;i <= IRZ;i++) {
			if(JointDia.res[i]) restraint |= (1 << i);
		}
		JOINT* joint;
		POSITION pos;
		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			if(joint->sel) {
				joint->restraint = restraint;
				joint->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::RefreshConstraints() {
	CONSTRAINT* pconst;
	JOINT* joint;
	POSITION pos,pos1;

	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->mconnect -= joint->cconnect;
		joint->cconnect = 0;
		joint->constraint = 0;
	}

	pos = constraints.GetHeadPosition();
	while(pos) {
         pconst = &constraints.GetNext(pos);

		 pos1 = pconst->jointlist.GetHeadPosition();
		 while(pos1) {
			 joint = pconst->jointlist.GetNext(pos1);
			 if(pconst->type == CONSTRAINT::DIAPHRAGM) {
				 if(pconst->constraint == IRX) joint->constraint |= UY | UZ | RX;
				 else if(pconst->constraint == IRY) joint->constraint |= UX | UZ | RY;
				 else if(pconst->constraint == IRZ) joint->constraint |= UX | UY | RZ; 
			 } else {
				 joint->constraint |= pconst->constraint;
			 }
			 joint->mconnect++;
			 joint->cconnect++;
		 }
	}
}
void CmyDocument::OnJointConstraint(UINT nID) {
	int type = nID - IDM_JOINT_CONSTRAINT - 1;
    CONSTRAINT* pbody;
	CDefineDia<CONSTRAINT> Dia(AfxGetMainWnd(),&constraints,"Constraint Database",&systems);
	if(Dia.DoModal() == IDOK && (pbody = (CONSTRAINT*)Dia.pT)) {
		POSITION pos,pos1;
		DOUBLE v;
		JOINT *joint;

		
		CList<DOUBLE,DOUBLE> zlist;
		CONSTRAINT constraint = *pbody;
		constraints.RemoveAt(constraints.Find(constraint));
		pbody = &constraint;


		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			if(!joint->sel) continue;

            if(type == IUX)  {
				if(zlist.Find(joint->p.x)) continue;
				zlist.AddTail(joint->p.x);
			} else if(type == IUY)  {
				if(zlist.Find(joint->p.y)) continue;
				zlist.AddTail(joint->p.y);
			} else if(type == IUZ)  {
				if(zlist.Find(joint->p.z)) continue;
				zlist.AddTail(joint->p.z);
			}
		}
		if(type == -1) {
			v = 0;
			zlist.AddTail(v);
		}

		pos1 = zlist.GetHeadPosition();
		while(pos1) {
			v = zlist.GetNext(pos1);
			
			pbody->jointlist.RemoveAll();
			
			pos = joints.GetHeadPosition();
			while(pos)	{
				joint = &joints.GetNext(pos);
				if(joint->sel) {
					if(type == -1){
					} else if(type == IUX) {
						if(v != joint->p.x) continue;
					} else if(type == IUY) {
						if(v != joint->p.y) continue;
					} else if(type == IUZ) {
						if(v != joint->p.z) continue;
					}
					joint->sel = false;
					pbody->jointlist.AddTail(joint);
				}
			}
			CString name = pbody->name;
			if(type != - 1) pbody->name.Format("%s_%.2f",name,v);
			constraints.AddTail(*pbody);
			pbody->name = name;
		}
	}
	RefreshConstraints();
	UpdateDrawing();
}
void CmyDocument::RemoveDesignCombos() {
	COMBINATION* pcomb;
    POSITION pos = combinations.GetHeadPosition(),ppos;
	while(pos) {
		ppos = pos;
		pcomb = &combinations.GetNext(pos);
		if(pcomb->designcombo) {
			combinations.RemoveAt(ppos);
		}
	}
}
void CmyDocument::AddDesignCombos() {
	const UINT NCOMB  = 5;
	BOOL hasWind = FALSE,hasEQ = FALSE;
	COMBINATION comb[NCOMB];
	LOADCASE* pcase;
	LOADCOMBO ldcombo;
	POSITION pos = loadcases.GetHeadPosition();
	while(pos) {
		pcase = &loadcases.GetNext(pos);
		ldcombo.loadcase = pcase;
        if(pcase->type == LOADCASE::DEAD) {
			ldcombo.FS = 1.3;
			comb[0].loadlist.AddTail(ldcombo);
			ldcombo.FS = 0.9 * 1.3;
			comb[1].loadlist.AddTail(ldcombo);
			comb[2].loadlist.AddTail(ldcombo);
			comb[3].loadlist.AddTail(ldcombo);
			comb[4].loadlist.AddTail(ldcombo);
		} else if(pcase->type == LOADCASE::LIVE) {
			ldcombo.FS = 1.6;
			comb[0].loadlist.AddTail(ldcombo);
			ldcombo.FS = 0.9 * 1.6;
			comb[1].loadlist.AddTail(ldcombo);
			comb[2].loadlist.AddTail(ldcombo);
			comb[3].loadlist.AddTail(ldcombo);
			comb[4].loadlist.AddTail(ldcombo);
		} else if(pcase->type == LOADCASE::WIND) {
			ldcombo.FS = 1.0;
			comb[1].loadlist.AddTail(ldcombo);
            ldcombo.FS = -1.0;
			comb[2].loadlist.AddTail(ldcombo);
			hasWind = TRUE;
		} else if(pcase->type == LOADCASE::QUAKE) {
			ldcombo.FS = 1.0;
			comb[3].loadlist.AddTail(ldcombo);
			ldcombo.FS = -1.0;
			comb[4].loadlist.AddTail(ldcombo);
			hasEQ = TRUE;
		}
	}

	/*add*/
	COMBINATION Envelope;
	Envelope.type = ENVELOPECOMBO;
	Envelope.name = "DENVELOPE";
	Envelope.designcombo = TRUE;

	CString name;
	for(int i = 0;i < NCOMB;i++) {
		if(!hasWind) {
			if(i == 1 || i == 2) continue;
		}
		if(!hasEQ) {
			if(i == 3 || i == 4) continue;
		}
		if(comb[i].loadlist.GetCount()) {
			name.Format("DCOMB%d",i + 1);
			comb[i].name = name;
			comb[i].designcombo = TRUE;
			combinations.AddTail(comb[i]);

			ldcombo.loadcase = &combinations.GetTail();
			ldcombo.FS = 1;
			Envelope.loadlist.AddTail(ldcombo);
		}
	}

	combinations.AddTail(Envelope);
}
void CmyDocument::OnAddDesignCombos() {
    RemoveDesignCombos();
	AddDesignCombos();
}
void CmyDocument::OnJointLoad(UINT nID) {
	int i;
	JLOAD load,*pload;
	JOINT tjoint;
	BOOL found;
	load.type = FORCE + (nID - IDM_JOINT_FORCE);
	CJointLoadDia JointDia(AfxGetMainWnd(),&loadcases,&systems,&tjoint,load.type,(LOADCASE*)c_AnalysisCase);
    if(JointDia.DoModal() == IDOK) {
		JOINT* joint;
		POSITION pos = joints.GetHeadPosition(),pos1,ppos;
		while(pos) {
			if((joint = &joints.GetNext(pos))->sel) {

				if(JointDia.mode >= F_REPLACE) {
                    pos1 = joint->load.GetHeadPosition();
					while(pos1) {
						ppos = pos1;
                        pload = &joint->load.GetNext(pos1);
						if(pload->type == JointDia.type 
							&& pload->loadcase == JointDia.ploadcase) {
							joint->load.RemoveAt(ppos);
						}
					}
				}
				found = FALSE;
				if(JointDia.mode != F_DELETE) {
					for(i = 0;i < 6;i++) {
						load.Q[i] = JointDia.Q[i];
						if(load.Q[i]) found = TRUE;
					}
					load.system = JointDia.system;
					load.loadcase = JointDia.ploadcase;
					if(found) {

						/*add*/
						found = FALSE;
                        pos1 = joint->load.GetHeadPosition();
						while(pos1) {
							pload = &joint->load.GetNext(pos1);
							if(pload->loadcase == load.loadcase 
								&& pload->type == load.type 
								&& pload->system == load.system) {
								for(i = 0;i < 6;i++) {
									pload->Q[i] += load.Q[i];
								}
								found = TRUE;
								break;
							}
						}
						if(!found)
						    joint->load.AddTail(load);
						/*end*/
					}
				}
				joint->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnJointAxis() {
	JOINT tjoint;
	CJointAxisDia JointDia(AfxGetMainWnd(),&tjoint);
    if(JointDia.DoModal() == IDOK) {
		JOINT* joint;
		POSITION pos = joints.GetHeadPosition();
		while(pos) {
			if((joint = &joints.GetNext(pos))->sel) {
				joint->rotation = JointDia.rotation;
				joint->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}

void CmyDocument::OnMemberAxis() {
	CAxisRotationDia Dia(AfxGetMainWnd(),0,"Member local axis");
    if(Dia.DoModal() == IDOK) {
		MEMBER* member;
		POSITION pos = members.GetHeadPosition();
		while(pos) {
			if((member = &members.GetNext(pos))->sel) {
				member->alpha = Dia.alpha;
				member->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnMemberDivisions() {
	CMemberDivisionsDia Dia(AfxGetMainWnd());
    if(Dia.DoModal() == IDOK) {
		MEMBER* member;
		POSITION pos = members.GetHeadPosition();
		while(pos) {
			if((member = &members.GetNext(pos))->sel) {
				member->nMinDiv = Dia.nMinDiv;
				member->nMinFrameDiv = Dia.nMinFrameDiv;
				member->DivAtInterim = Dia.DivAtInterim;
				member->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnSlabDivisions() {
	CSlabDivisionsDia Dia(AfxGetMainWnd());
    if(Dia.DoModal() == IDOK) {
		SLAB* slab;
		POSITION pos = slabs.GetHeadPosition();
		while(pos) {
			if((slab = &slabs.GetNext(pos))->sel) {
				slab->nDivx = Dia.nDivx;
				slab->nDivy = Dia.nDivy;
				slab->edge_constraint = Dia.constraint;
				slab->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnMemberPointLoad() {
	CMemberLoadDia MemberDia(AfxGetMainWnd(),&loadcases,&systems,FALSE,(LOADCASE*) c_AnalysisCase);
    if(MemberDia.DoModal() == IDOK) {
		MEMBER* member;
		LOAD load,*pload;
		BOOL found;
		POSITION pos = members.GetHeadPosition(),pos1,ppos;
		while(pos) {
			member = &members.GetNext(pos);
			if(member->sel) {
				if(MemberDia.mode >= F_REPLACE) {
					pos1 = member->load.GetHeadPosition();
					while(pos1) {
						ppos = pos1;
                        pload = &member->load.GetNext(pos1);
						if(pload->loadcase == MemberDia.ploadcase
							&& pload->type == CONCENTRATED)
							member->load.RemoveAt(ppos);
					}
				}
				if(MemberDia.mode != F_DELETE) {
					DOUBLE L = member->GetLength();
					for(int i = 0;i < 4;i++) {
						if(MemberDia.Q[i]) {
							load.type = CONCENTRATED;
							load.P = MemberDia.Q[i];
							load.x = MemberDia.d[i] * (MemberDia.use_relative ? L : 1);
							load.dir = MemberDia.dir;
							load.system = MemberDia.system;
							load.loadcase = MemberDia.ploadcase;

							/*add*/
							found = FALSE;
							pos1 = member->load.GetHeadPosition();
							while(pos1) {
								pload = &member->load.GetNext(pos1);
								if(pload->loadcase == load.loadcase 
									&& pload->type == load.type 
									&& pload->system == load.system
									&& pload->dir == load.dir
									&& EQUAL(pload->x,load.x)) {
									pload->P += load.P;
									found = TRUE;
									break;
								}
							}
							if(!found)
								member->load.AddTail(load);
							/*end*/
						}
					}
				}
				member->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnMemberDistLoad() {
	CMemberLoadDia MemberDia(AfxGetMainWnd(),&loadcases,&systems,TRUE,(LOADCASE*) c_AnalysisCase);
    if(MemberDia.DoModal() == IDOK) {
		MEMBER* member;
		LOAD load,*pload;
		BOOL found;
		POSITION pos = members.GetHeadPosition(),pos1,ppos;
		while(pos) {
			member = &members.GetNext(pos);
			if(member->sel) {
				if(MemberDia.mode >= F_REPLACE) {
					pos1 = member->load.GetHeadPosition();
					while(pos1) {
						ppos = pos1;
                        pload = &member->load.GetNext(pos1);
						if(pload->loadcase == MemberDia.ploadcase
							&& (pload->type == UNIFORM || pload->type == TRAPEZOIDAL))
							member->load.RemoveAt(ppos);
					}
				}
				if(MemberDia.mode != F_DELETE) {
					DOUBLE L = member->GetLength();
					for(int i = 1;i < 4;i++) {
						if((MemberDia.Q[i] || MemberDia.Q[i - 1]) &&
						   (MemberDia.d[i] > MemberDia.d[i - 1])
						   ) {
							load.type = TRAPEZOIDAL;
							load.P1 = MemberDia.Q[i];
							load.x1 = MemberDia.d[i] * (MemberDia.use_relative ? L : 1);
							load.P = MemberDia.Q[i - 1];
							load.x = MemberDia.d[i - 1] * (MemberDia.use_relative ? L : 1);
							load.dir = MemberDia.dir;
							load.system = MemberDia.system;
							load.loadcase = MemberDia.ploadcase;

							/*add*/
							found = FALSE;
							pos1 = member->load.GetHeadPosition();
							while(pos1) {
								pload = &member->load.GetNext(pos1);
								if(pload->loadcase == load.loadcase 
									&& pload->system == load.system
									&& pload->dir == load.dir
									) {
									if(pload->type == UNIFORM
										&& EQUAL(load.x,0) && EQUAL(load.x1,L)) {
										pload->type = TRAPEZOIDAL;
										pload->P = pload->P + load.P;
                                        pload->P1 = pload->P + load.P1;
										pload->x1 = load.x1;
										found = TRUE;
										break;
									} else if(pload->type == TRAPEZOIDAL 
										&& EQUAL(pload->x,load.x) && EQUAL(pload->x1,load.x1))  {
										pload->P += load.P;
										pload->P1 += load.P1;
										found = TRUE;
										break;
									}
								}
							}
							if(!found)
								member->load.AddTail(load);
							/*end*/
						}
					}
					if(MemberDia.Qu) {
						load.type = UNIFORM;
						load.P = MemberDia.Qu;
						load.dir = MemberDia.dir;
						load.system = MemberDia.system;
						load.loadcase = MemberDia.ploadcase;
						
						/*add*/
						found = FALSE;
						pos1 = member->load.GetHeadPosition();
						while(pos1) {
							pload = &member->load.GetNext(pos1);
							if(pload->loadcase == load.loadcase 
								&& pload->system == load.system
								&& pload->dir == load.dir
								) {
								if(pload->type == UNIFORM) {
									pload->P += load.P;
									found = TRUE;
									break;
								} else if(pload->type == TRAPEZOIDAL 
									&& EQUAL(pload->x,0) && EQUAL(pload->x1,L))  {
									pload->P += load.P;
									pload->P1 += load.P;
									found = TRUE;
									break;
								}
							}
						}
						if(!found)
							member->load.AddTail(load);
						/*end*/
					}
				}
				member->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnSlabAxis() {
	CAxisRotationDia Dia(AfxGetMainWnd(),0,"Slab local axis");
    if(Dia.DoModal() == IDOK) {
		SLAB* slab;
		POSITION pos = slabs.GetHeadPosition();
		while(pos) {
			if((slab = &slabs.GetNext(pos))->sel) {
				slab->alpha = Dia.alpha;
				slab->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnReverseSlabAxis() {
	SLAB* slab;
	POSITION pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		if(slab->sel) {
			UINT i,j;
			JOINT* temp;
			for(i = 0;i < slab->NJ;i++) {
				j = slab->NJ - 1 - i;
				if(i >= j) continue;
				temp = slab->jt[i];
				slab->jt[i] = slab->jt[j];
				slab->jt[j] = temp;
			}
			slab->sel = false;
		}
	}
	/*update*/
	SetModifiedFlag();
	UpdateDrawing();
}
void CmyDocument::OnSlabDistLoad() {
	CSlabDistLoadDia SlabDia(AfxGetMainWnd(),&loadcases,&systems,(LOADCASE*) c_AnalysisCase);
    if(SlabDia.DoModal() == IDOK) {
		BOOL found;
		SLAB* slab;
		LOAD load,*pload;
		POSITION pos,pos1;

		pos = slabs.GetHeadPosition();
		while(pos) {
			slab = &slabs.GetNext(pos);
			if(slab->sel) {
				if(SlabDia.mode >= F_REPLACE) {
					slab->load.RemoveAll();
				}
				if(SlabDia.mode != F_DELETE) {
					if(SlabDia.Qu) {
						load.type = UNIFORM;
						load.P = SlabDia.Qu;
						load.dir = SlabDia.dir;
						load.system = SlabDia.system;
						load.loadcase = SlabDia.ploadcase;

						/*add*/
						found = FALSE;
						pos1 = slab->load.GetHeadPosition();
						while(pos1) {
							pload = &slab->load.GetNext(pos1);
							if(pload->loadcase == load.loadcase 
								&& pload->system == load.system
								&& pload->dir == load.dir
								) {
								pload->P += load.P;
								found = TRUE;
								break;
							}
						}
						if(!found)
							slab->load.AddTail(load);
						/*end*/
					}
				}
				slab->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnStrain(UINT nID) {
	CStrainLoadDia myDia(AfxGetMainWnd(),&loadcases,(LOADCASE*)c_AnalysisCase,(nID == IDM_SLAB_STRAIN));
    if(myDia.DoModal() == IDOK) {
		BOOL found;
		LOAD load,*pload;
		POSITION pos,pos1;
		ELEMENT* pelem;
		PELEMENTLIST pelements;
		if(nID == IDM_MEMBER_STRAIN) pelements = (PELEMENTLIST)&members;
		else pelements = (PELEMENTLIST)&slabs;

		pos = pelements->GetHeadPosition();
		while(pos) {
			pelem = &pelements->GetNext(pos);
			if(pelem->sel) {
				
				if(myDia.mode >= F_REPLACE) {
					pelem->load.RemoveAll();
				}
				if(myDia.mode != F_DELETE) {
					if(myDia.Value) {
						if(myDia.type == myDia.TEMP_STRAIN) load.type = TEMPERATURE;
						else load.type = STRAIN;
						load.P = myDia.Value;
						load.dir = myDia.dir;
						load.loadcase = myDia.ploadcase;

						/*add*/
						found = FALSE;
						pos1 = pelem->load.GetHeadPosition();
						while(pos1) {
							pload = &pelem->load.GetNext(pos1);
							if(pload->loadcase == load.loadcase 
								&& pload->system == load.system
								&& pload->dir == load.dir
								) {
								pload->P += load.P;
								found = TRUE;
								break;
							}
						}
						if(!found)
							pelem->load.AddTail(load);
						/*end*/
					}
				}
				pelem->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
}
void CmyDocument::OnMemberStrain() {
	OnStrain(IDM_MEMBER_STRAIN);
}
void CmyDocument::OnSlabStrain() {
	OnStrain(IDM_SLAB_STRAIN);
}
void CmyDocument::OnAssignSection() {
	CDefineDia<SECTION> Dia(AfxGetMainWnd(),&sections,"Section Database",&materials,NULL,TRUE);
    if(Dia.DoModal() == IDOK) {
		MEMBER* mem;
		POSITION pos = members.GetHeadPosition();
		while(pos) {
			mem = &members.GetNext(pos); 
			if(mem->sel) {
				mem->section = (SECTION*)Dia.pT;
				mem->e_section = 0;
				mem->sel = false;
			}
		}
	}
	UpdateDrawing();
}
void CmyDocument::OnAssignSlabSection() {
	CDefineDia<ASECTION> Dia(AfxGetMainWnd(),&asections,"Area Section Database",&materials,NULL,TRUE);
    if(Dia.DoModal() == IDOK) {
		SLAB* slab;
		POSITION pos = slabs.GetHeadPosition();
		while(pos) {
			slab = &slabs.GetNext(pos); 
			if(slab->sel) {
				slab->section = (ASECTION*)Dia.pT;
				slab->sel = false;
			}
		}
	}
	UpdateDrawing();
}
void CmyDocument::OnAssignTaperedSection() {
	CTaperedSectionDia Dia(AfxGetMainWnd(),&sections,&materials);
    if(Dia.DoModal() == IDOK) {
		MEMBER* mem;
		POSITION pos = members.GetHeadPosition();
		while(pos) {
			mem = &members.GetNext(pos); 
			if(mem->sel) {
				mem->section = Dia.s_sec;
                mem->e_section = Dia.e_sec;
				mem->EIyy = Dia.EIyy;
				mem->EIzz = Dia.EIzz;
				mem->start_offset = Dia.soff;
				mem->end_offset = Dia.eoff;
				mem->sel = false;
			}
		}
	}
	UpdateDrawing();
}
void CmyDocument::OnAssignGroup() {
    GROUP* pgroup;
	CDefineDia<GROUP> Dia(AfxGetMainWnd(),&groups,"Group Database");
	if(Dia.DoModal() == IDOK) {
		POSITION pos;
		JOINT *joint;
		MEMBER* mem;
		SLAB* sla;
		pgroup = (GROUP*)Dia.pT;
		
		pgroup->jointlist.RemoveAll();
		pgroup->memberlist.RemoveAll();
		pgroup->slablist.RemoveAll();

		pos = joints.GetHeadPosition();
		while(pos)	{
			joint = &joints.GetNext(pos);
			if(joint->sel) {
				joint->sel = false;
				pgroup->jointlist.AddTail(joint);
			}
		}

		pos = members.GetHeadPosition();
		while(pos)	{
			mem = &members.GetNext(pos);
			if(mem->sel) {
				mem->sel = false;
				pgroup->memberlist.AddTail(mem);
			}
		}

		pos = slabs.GetHeadPosition();
		while(pos)	{
			sla = &slabs.GetNext(pos);
			if(sla->sel) {
				sla->sel = false;
				pgroup->slablist.AddTail(sla);
			}
		}
	}
	UpdateDrawing();
}
void CmyDocument::OnJoinMembers() {
	UINT count,total,i;
	MEMBERPLIST list,sel_list;
	MEMBER* member,*pmember = NULL;

	POSITION pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		if(member->sel) {
			sel_list.AddTail(member);
		}
	}

	count = 0;
	total = sel_list.GetCount();
	pos = members.GetHeadPosition();
	for(i = 0;i < total;i++) {
		member = &members.GetNext(pos);
		if(member->sel) {
			count++;
			if(count == 2) {
				if(member->j1 == pmember->j2 
					|| member->j1 == pmember->j1
					|| member->j2 == pmember->j1
					|| member->j2 == pmember->j2
					) {
					list.RemoveAll();
					if(member->j1 == pmember->j2 || member->j2 == pmember->j2) {
						if(member->j2 == pmember->j2) member->SwapEnds();
						list.AddTail(pmember);
						list.AddTail(member);
					} else {
						if(member->j1 == pmember->j1) member->SwapEnds();
						list.AddTail(member);
						list.AddTail(pmember);
					}
					
					MEMBER::Join(&list,&members);
					DeleteMember(member);
					DeleteMember(pmember);
					pmember = &members.GetTail();
				} else {
					pmember = member;
				}
				count--;
			} else {
				pmember = member;
			}
			member->sel = false;
		}
	}

	SetModifiedFlag();
	UpdateDrawing();
}
void CmyDocument::OnBreakMember() {
	CMemberBreakDia MemberDia(AfxGetMainWnd());
    if(MemberDia.DoModal() == IDOK) {

		UINT breakn = MemberDia.breakn,i;
		DOUBLE breakr = 1.0f / MemberDia.breakr,sum = 0,sum1;
		for(i = 0; i < breakn;i++) sum += pow(breakr,i);

		MEMBER* member;
		JOINT* pjoint;
		POSITION pos = members.GetHeadPosition(),ppos;
		while(pos) {
			ppos = pos;
			member = &members.GetNext(pos);
			if(member->sel) {
				member->sel = false;

				JOINTPLIST list;
				list.AddTail(member->j1);
				sum1 = 0;
				for(i = 0;i < breakn - 1;i++) {
					sum1 += pow(breakr,i);
					JOINT joint;
					joint.p = member->j1->p + (member->j2->p - member->j1->p) * (sum1 / sum);
					pjoint = AddJoint(&joint);
					list.AddTail(pjoint);
				}
				list.AddTail(member->j2);
				member->Break(&list,&members);
				DeleteMember(member);
			}
		}
	}
    SetModifiedFlag();
	UpdateDrawing();
}
void CmyDocument::OnMemberRelease() {
	MEMBER tmember;
	CMemberReleaseDia MemberDia(AfxGetMainWnd(),&tmember);
    if(MemberDia.DoModal() == IDOK) {

		UBMP8 nrelease = NOUR,frelease = NOUR;
		if(MemberDia.r[0]) nrelease |= UX;
		if(MemberDia.r[1]) nrelease |= UY;
		if(MemberDia.r[2]) nrelease |= UZ;
		if(MemberDia.r[3]) nrelease |= RX;
		if(MemberDia.r[4]) nrelease |= RY;
		if(MemberDia.r[5]) nrelease |= RZ;
		if(MemberDia.r[6]) frelease |= UX;
		if(MemberDia.r[7]) frelease |= UY;
		if(MemberDia.r[8]) frelease |= UZ;
		if(MemberDia.r[9]) frelease |= RX;
		if(MemberDia.r[10]) frelease |= RY;
		if(MemberDia.r[11]) frelease |= RZ;

		/*assign the members*/
		MEMBER* member;
		POSITION pos;
		pos = members.GetHeadPosition();
		while(pos) {
			member = &members.GetNext(pos);
			if(member->sel) {
				member->nrelease = nrelease;
				member->frelease = frelease;
				member->sel = false;
			}
		}
		/*update*/
		SetModifiedFlag();
		UpdateDrawing();
	}
	
}
void CmyDocument::OnJoinSlabs() {
}
void CmyDocument::OnBreakSlab() {
	CSlabDivisionsDia Dia(AfxGetMainWnd());
    if(Dia.DoModal() == IDOK) {
		POSITION pos;
		SLAB* slab;

        /*mesh selected slabs*/
		slab_elements.RemoveAll();
		pos = slabs.GetHeadPosition();
        while(pos) {
			slab = &slabs.GetNext(pos);
			if(slab->sel) {
                slab->sel = false;
				slab->nDivx = Dia.nDivx;
				slab->nDivy = Dia.nDivy;
				slab->edge_constraint = Dia.constraint;

				/*Mesh*/
				if(slab->Mesh(&slabs,&joints,&slab_elements))
					DeleteSlab(slab);
			}
		}
		/*delete slab objects*/
		SLABPLIST* pslist;
		pos = slab_elements.GetHeadPosition();
		while(pos) {
			pslist = slab_elements.GetNext(pos);
			pslist->RemoveAll();
			delete[] pslist;
		}
		slab_elements.RemoveAll();
	}
    SetModifiedFlag();
	UpdateDrawing();
}
void CmyDocument::OnDefineMaterial() {
	POSITION pos;
	SECTION* section;
	ASECTION* asection;

    pos = sections.GetHeadPosition();
	while(pos) {
		section = &sections.GetNext(pos);
		if(section->material->rank == LEVEL1)
            section->material->rank = LEVEL2;
	}

	pos = asections.GetHeadPosition();
	while(pos) {
		asection = &asections.GetNext(pos);
		if(asection->material->rank == LEVEL1)
            asection->material->rank = LEVEL2;
	}

	CDefineDia<MATERIAL> Dia(AfxGetMainWnd(),&materials,"Material Database",NULL,NULL,TRUE);
    Dia.DoModal();
}
void CmyDocument::OnDefineFrameSection() {
	POSITION pos;
	MEMBER* member;
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		if(member->section->rank == LEVEL1)
            member->section->rank = LEVEL2;
	}

	CDefineDia<SECTION> Dia(AfxGetMainWnd(),&sections,"Section Database",&materials,NULL,TRUE);
    Dia.DoModal();
}
void CmyDocument::OnDefineSlabSection() {
	POSITION pos;
	SLAB* slab;
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		if(slab->section->rank == LEVEL1)
            slab->section->rank = LEVEL2;
	}

	CDefineDia<ASECTION> Dia(AfxGetMainWnd(),&asections,"Section Database",&materials,NULL,TRUE);
    Dia.DoModal();
}
void CmyDocument::OnDefineSystem() {
	CDefineDia<SYSTEM> Dia(AfxGetMainWnd(),&systems,"Coordinate System Database");
	if(Dia.DoModal() == IDOK) {
		CmyView* pView;
		POSITION pos = GetFirstViewPosition();
		while(pos) {
			pView = (CmyView*) GetNextView(pos);
			if(pView->view != VIEW_3D) {
				GRIDLIST* pgrid = &global->grid[NORMGRID[pView->view]];
				if(!pgrid->IsEmpty()) {
					pView->position = pgrid->GetCount() - 1;
					pView->pvalue = pgrid->GetAt(pgrid->FindIndex(pView->position));
				}
			}
			pView->InitView();
		}
	}
}
void CmyDocument::CheckCaseInCombo(int type) {
	POSITION pos,pos1;
	COMBINATION* combo;
	NLCASE* nlcase;
	BUCKLINGCASE* bcase;

	LOADCOMBO* loadcombo;
	pos = combinations.GetHeadPosition();
	while(pos) {
		combo = &combinations.GetNext(pos);
		pos1 = combo->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &combo->loadlist.GetNext(pos1);
			if(loadcombo->loadcase->acase_type != type) continue;
			if(loadcombo->loadcase->rank == LEVEL1)
				loadcombo->loadcase->rank = LEVEL2;
		}
	}

	pos = nlcases.GetHeadPosition();
	while(pos) {
		nlcase = &nlcases.GetNext(pos);
		pos1 = nlcase->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &nlcase->loadlist.GetNext(pos1);
			if(loadcombo->loadcase->acase_type != type) continue;
			if(loadcombo->loadcase->rank == LEVEL1)
				loadcombo->loadcase->rank = LEVEL2;
		}
	}

	pos = bucklingcases.GetHeadPosition();
	while(pos) {
		bcase = &bucklingcases.GetNext(pos);
		pos1 = bcase->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &bcase->loadlist.GetNext(pos1);
			if(loadcombo->loadcase->acase_type != type) continue;
			if(loadcombo->loadcase->rank == LEVEL1)
				loadcombo->loadcase->rank = LEVEL2;
		}
	}
}
void CmyDocument::OnDefineLoadCase() {
	CheckCaseInCombo(ANALYSISCASE::LOAD_CASE);
	CDefineDia<LOADCASE> Dia(AfxGetMainWnd(),&loadcases,"Load Case Database");
    Dia.DoModal();
	RemoveDesignCombos();
	FillCases();
}
void CmyDocument::OnDefineLoadCombo() {
	CheckCaseInCombo(ANALYSISCASE::COMBO_CASE);
	CDefineDia<COMBINATION> Dia(AfxGetMainWnd(),&combinations,"Load Combination Database",
		&all_analysis_cases);
    Dia.DoModal();
}
void CmyDocument::OnDefineNonLinearCase() {
	CheckCaseInCombo(ANALYSISCASE::NL_CASE);
	CDefineDia<NLCASE> Dia(AfxGetMainWnd(),&nlcases,"Non linear analysis case database",
	     &all_analysis_cases);
    Dia.DoModal();
}
void CmyDocument::OnDefineBucklingCase() {
	CDefineDia<BUCKLINGCASE> Dia(AfxGetMainWnd(),&bucklingcases,"Buckling analysis case database",
	     &all_analysis_cases);
    Dia.DoModal();
}
void CmyDocument::OnDefineModalCase() {
	CDefineDia<MODALCASE> Dia(AfxGetMainWnd(),&modalcases,"Modal Case Database");
    Dia.DoModal();
}
void CmyDocument::OnDefineResponseHistory() {
	CheckCaseInCombo(ANALYSISCASE::RESPONSEH_CASE);
	CDefineDia<RESPONSEHIST> Dia(AfxGetMainWnd(),&responsecases,"Response Analysis Case Database",&modalcases,&rhfunctions);
    Dia.DoModal();
}
void CmyDocument::OnDefineResponseSpectrum() {
	CheckCaseInCombo(ANALYSISCASE::RESPONSES_CASE);
	CDefineDia<RESPONSESPEC> Dia(AfxGetMainWnd(),&responsespecs,"Response Spectrum Case Database",&modalcases,&rsfunctions);
    Dia.DoModal();
}
void CmyDocument::OnDefineConstraint() {
	CDefineDia<CONSTRAINT> Dia(AfxGetMainWnd(),&constraints,"Constraint Database",&systems);
	Dia.DoModal();
    RefreshConstraints();
	UpdateDrawing();
}
void CmyDocument::CheckCaseInResp() {
	POSITION pos,pos1;
	RESPONSEHIST* resp;
	SPECFUNC* func;
	pos = responsecases.GetHeadPosition();
	while(pos) {
		resp = &responsecases.GetNext(pos);
		pos1 = resp->funclist.GetHeadPosition();
		while(pos1) {
			func = &resp->funclist.GetNext(pos1);
			if(func->function->rank == LEVEL1)
				func->function->rank = LEVEL2;
		}
	}
}
void CmyDocument::OnDefineTimeHistoryFunction() {
	CheckCaseInResp();
	FUNCTION::defspectrumtype = FALSE;
	CDefineDia<FUNCTION> Dia(AfxGetMainWnd(),&rhfunctions,"Time History Function Database");
    Dia.DoModal();
}
void CmyDocument::OnDefineSpectrumFunction() {
	CheckCaseInResp();
	FUNCTION::defspectrumtype = TRUE;
	CDefineDia<FUNCTION> Dia(AfxGetMainWnd(),&rsfunctions,"Response Spectrum Function Database");
    Dia.DoModal();
}
void CmyDocument::OnDefineGroup() {
	CDefineDia<GROUP> Dia(AfxGetMainWnd(),&groups,"Group Database");
    Dia.DoModal();
}
void CmyDocument::OnSelectAllJoints() {
    POSITION pos = joints.GetHeadPosition();
	while(pos)	joints.GetNext(pos).sel = true;
	UpdateDrawing();
}
void CmyDocument::OnSelectAllMembers() {
	POSITION pos = members.GetHeadPosition();
	while(pos)	members.GetNext(pos).sel = true;
	UpdateDrawing();
}
void CmyDocument::OnSelectAllSlabs() {
	POSITION pos = slabs.GetHeadPosition();
	while(pos)	slabs.GetNext(pos).sel = true;
	UpdateDrawing();
}
void CmyDocument::OnSelectAll() {
	POSITION pos = joints.GetHeadPosition();
	while(pos)	joints.GetNext(pos).sel = true;
	pos = members.GetHeadPosition();
	while(pos)	members.GetNext(pos).sel = true;
	pos = slabs.GetHeadPosition();
	while(pos)	slabs.GetNext(pos).sel = true;
	UpdateDrawing();
}
void CmyDocument::OnDeselectAllJoints() {
    POSITION pos = joints.GetHeadPosition();
	while(pos)	joints.GetNext(pos).sel = false;
	UpdateDrawing();
}
void CmyDocument::OnDeselectAllMembers() {
	POSITION pos = members.GetHeadPosition();
	while(pos)	members.GetNext(pos).sel = false;
	UpdateDrawing();
}
void CmyDocument::OnDeselectAllSlabs() {
	POSITION pos = slabs.GetHeadPosition();
	while(pos)	slabs.GetNext(pos).sel = false;
	UpdateDrawing();
}
void CmyDocument::OnDeselectAll() {
	POSITION pos = joints.GetHeadPosition();
	while(pos)	joints.GetNext(pos).sel = false;
	pos = members.GetHeadPosition();
	while(pos)	members.GetNext(pos).sel = false;
	pos = slabs.GetHeadPosition();
	while(pos)	slabs.GetNext(pos).sel = false;
	UpdateDrawing();
}
void CmyDocument::OnSelectLabels(UINT nID) {
	BOOL select = (nID == IDM_SELECT_LABELS);
	CSelectLabelsDia Dia(AfxGetMainWnd(),select);
	if(Dia.DoModal() == IDOK)
		UpdateDrawing();
}
void CmyDocument::OnSelectSection(UINT nID) {
	POSITION pos,pos1;
	MEMBER* member;
	CString name;

	CSelectDia<SECTION> Dia(AfxGetMainWnd(),&sections,"Select Sections");
    if(Dia.DoModal() == IDOK) {
        pos = members.GetHeadPosition();
        while(pos) {
			member = &members.GetNext(pos);
			pos1 = Dia.selected_list.GetHeadPosition();
			while(pos1) {
				name = Dia.selected_list.GetNext(pos1);
				if(member->section->name == name) {
					member->sel = (nID == IDM_SELECT_SECTION);
				}
			}
		}
		UpdateDrawing();
	}
}
void CmyDocument::OnSelectSlabSection(UINT nID) {
	POSITION pos,pos1;
	SLAB* slab;
	CString name;

	CSelectDia<ASECTION> Dia(AfxGetMainWnd(),&asections,"Select Area Sections");
    if(Dia.DoModal() == IDOK) {
		pos = slabs.GetHeadPosition();
        while(pos) {
			slab = &slabs.GetNext(pos);
			pos1 = Dia.selected_list.GetHeadPosition();
			while(pos1) {
				name = Dia.selected_list.GetNext(pos1);
				if(slab->section->name == name) {
					slab->sel = (nID == IDM_SELECT_ASECTION);
				}
			}
		}
		UpdateDrawing();
	}
}
void CmyDocument::OnSelectConstraint(UINT nID) {
	POSITION pos,pos1;
	CONSTRAINT* cons;
	CONSTRAINT tcons;
	JOINT* joint;

	CSelectDia<CONSTRAINT> Dia(AfxGetMainWnd(),&constraints,"Select Constraints");
    if(Dia.DoModal() == IDOK) {
		pos1 = Dia.selected_list.GetHeadPosition();
		while(pos1) {
			tcons.name = Dia.selected_list.GetNext(pos1);

			if(pos = constraints.Find(tcons)) {
				cons = &constraints.GetAt(pos);
				pos = cons->jointlist.GetHeadPosition();
				while(pos) {
					joint = cons->jointlist.GetNext(pos);
					joint->sel = (nID == IDM_SELECT_CONSTRAINT);
				}
			}
		}
		UpdateDrawing();
	}
}
void CmyDocument::OnSelectGroup(UINT nID) {
	POSITION pos,pos1;
	GROUP* grp;
	GROUP tgrp;
	JOINT* joint;
	MEMBER* mem;
	SLAB* sla;

	CSelectDia<GROUP> Dia(AfxGetMainWnd(),&groups,"Select Groups");
    if(Dia.DoModal() == IDOK) {
		pos1 = Dia.selected_list.GetHeadPosition();
		while(pos1) {
			tgrp.name = Dia.selected_list.GetNext(pos1);

			if(pos = groups.Find(tgrp)) {
				grp = &groups.GetAt(pos);

				pos = grp->jointlist.GetHeadPosition();
				while(pos) {
					joint = grp->jointlist.GetNext(pos);
					joint->sel = (nID == IDM_SELECT_GROUP);
				}
				pos = grp->memberlist.GetHeadPosition();
				while(pos) {
					mem = grp->memberlist.GetNext(pos);
					mem->sel = (nID == IDM_SELECT_GROUP);
				}
				pos = grp->slablist.GetHeadPosition();
				while(pos) {
					sla = grp->slablist.GetNext(pos);
					sla->sel = (nID == IDM_SELECT_GROUP);
				}
			}
		}
		UpdateDrawing();
	}
}
void CmyDocument::OnCopy(UINT nID) {
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	POSITION pos;

	copyJoints.RemoveAll();
	copyMembers.RemoveAll();
	copySlabs.RemoveAll();

	if(nID == IDM_CUT)
		cutStatus = TRUE;
	else
		cutStatus = FALSE;

    pos = joints.GetHeadPosition();
	while(pos)	{
		joint = &joints.GetNext(pos);
		if(joint->sel) {
			copyJoints.AddTail(*joint);
			if(!cutStatus) joint->sel = false;
		}
	}
	
	pos = members.GetHeadPosition();
	while(pos)	{
		member = &members.GetNext(pos);
		if(member->sel) {
			copyMembers.AddTail(*member);
			if(!cutStatus) member->sel = false;
		}
	}

	pos = slabs.GetHeadPosition();
	while(pos)	{
		slab = &slabs.GetNext(pos);
		if(slab->sel) {
			copySlabs.AddTail(*slab);
			if(!cutStatus) slab->sel = false;
		}
	}
}
void CmyDocument::OnPaste() {
	JOINT* joint;
	MEMBER* member,*pmem;
	SLAB* slab;
	POSITION pos;

	CPasteDia PasteDia(AfxGetMainWnd());
    if(PasteDia.DoModal() == IDOK) {
		
		pos = copyJoints.GetHeadPosition();
		while(pos)	{
			joint = &copyJoints.GetNext(pos);
			JOINT j;
			j.p = joint->p + PasteDia.v;
			AddJoint(&j);

			if(cutStatus) DeleteJoint(joint);
		}
		
		pos = copyMembers.GetHeadPosition();
		while(pos)	{
			member = &copyMembers.GetNext(pos);
			JOINT j1,j2;
			j1.p = member->j1->p + PasteDia.v;
			j2.p = member->j2->p + PasteDia.v;
			pmem = AddMember(&j1,&j2);
			if(member->is_curved) {
                pmem->CopyCurveData(member);
			}

			if(cutStatus) DeleteMember(member);
		}

		pos = copySlabs.GetHeadPosition();
		while(pos)	{
			slab = &copySlabs.GetNext(pos);
			AddSlab(slab);
			if(cutStatus) DeleteSlab(slab);
		}
	}

	UpdateDrawing();
}
void CmyDocument::OnDelete() {
	POSITION pos;

	/*delete slab*/
	SLAB* slab;
    pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		if(slab->sel) {
			UndoRedo.Save(slab,FALSE);
			DeleteSlab(slab);
		}
	}

	/*delete member*/
	MEMBER* mem;
    pos = members.GetHeadPosition();
	while(pos) {
		mem = &members.GetNext(pos);
		if(mem->sel) {
			UndoRedo.Save(mem,FALSE);
			DeleteMember(mem);
		}
	}

	/*joints*/
	JOINT* joint;
    pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		if(joint->sel) {
			UndoRedo.Save(joint,FALSE);
			DeleteJoint(joint);
		}
	}

	SetModifiedFlag();
	UpdateDrawing();
}
void CmyDocument::OnInvertSelection() {
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
    POSITION pos = joints.GetHeadPosition();
	while(pos)	{
		joint = &joints.GetNext(pos);
		joint->sel = !joint->sel;
	}
	pos = members.GetHeadPosition();
	while(pos)	{
		member = &members.GetNext(pos);
		member->sel = !member->sel;
	}

	pos = slabs.GetHeadPosition();
	while(pos)	{
		slab = &slabs.GetNext(pos);
		slab->sel = !slab->sel;
	}
	UpdateDrawing();
}
void CmyDocument::OnSelectFromAll() {
	select_all = !select_all;
}
void CmyDocument::OnUpdateSelectFromAllUI(CCmdUI* pCmdUI) {
    pCmdUI->SetCheck(select_all);
}
void CmyDocument::OnUpdateDeselectUI(CCmdUI* pCmdUI) {
	
	BOOL enablem = FALSE,enablej = FALSE, enables = FALSE;

	POSITION pos = slabs.GetHeadPosition();
	while(pos)	{
		if(slabs.GetNext(pos).sel) {
			enables = TRUE;
			break;
		}
	}

    pos = members.GetHeadPosition();
	while(pos)	{
		if(members.GetNext(pos).sel) {
			enablem = TRUE;
			break;
		}
	}

    pos = joints.GetHeadPosition();
	while(pos)	{
		if(joints.GetNext(pos).sel) {
			enablej = TRUE;
			break;
		}
	}

	if(pCmdUI->m_nID == IDM_DESELECT_ALL_SLABS || pCmdUI->m_nID == IDM_DESELECT_ASECTION)
		pCmdUI->Enable(enables);
	else if(pCmdUI->m_nID == IDM_DESELECT_ALL_MEMBERS || pCmdUI->m_nID == IDM_DESELECT_SECTION)
		pCmdUI->Enable(enablem);
	else if(pCmdUI->m_nID == IDM_DESELECT_ALL_JOINTS || pCmdUI->m_nID == IDM_CONSTRAINT)
		pCmdUI->Enable(enablej);
	else if(pCmdUI->m_nID == IDM_DESELECT_CONSTRAINT)
		pCmdUI->Enable(enablej);
	else if(pCmdUI->m_nID == ID_DELETE || pCmdUI->m_nID == IDM_DELETE)
		pCmdUI->Enable((enables || enablem || enablej) && !Lock);
	else
		pCmdUI->Enable(enables || enablem || enablej);
}
void CmyDocument::OnUpdateAssignJointUI(CCmdUI* pCmdUI) {
	BOOL enable = FALSE;
    POSITION pos = joints.GetHeadPosition();
	while(pos)	{
		if(joints.GetNext(pos).sel) {
			enable = TRUE;
			break;
		}
	}
	pCmdUI->Enable(enable && !Lock);
}
void CmyDocument::OnUpdateAssignMemberUI(CCmdUI* pCmdUI) {
	BOOL enable = FALSE;
    POSITION pos = members.GetHeadPosition();
	while(pos)	{
		if(members.GetNext(pos).sel) {
			enable = TRUE;
			break;
		}
	}
	pCmdUI->Enable(enable && !Lock);
}
void CmyDocument::OnUpdateAssignSlabUI(CCmdUI* pCmdUI) {
	BOOL enable = FALSE;
    POSITION pos = slabs.GetHeadPosition();
	while(pos)	{
		if(slabs.GetNext(pos).sel) {
			enable = TRUE;
			break;
		}
	}
	pCmdUI->Enable(enable && !Lock);
}
void CmyDocument::OnSnapTo(UINT nID) {
	if(nID == ID_SNAP_TO_GRID) Snap_To_Grid = !Snap_To_Grid;
	else if(nID == ID_SNAP_TO_JOINT) Snap_To_Joint = !Snap_To_Joint;
    else if(nID == ID_SNAP_TO_MEMBERMID) Snap_To_MemberMid = !Snap_To_MemberMid;
	else if(nID == ID_SNAP_TO_MEMBER) Snap_To_Member = !Snap_To_Member;
}
void CmyDocument::OnUpdateSnapToUI(CCmdUI* pCmdUI) {
	UINT nID = pCmdUI->m_nID;
	if(nID == ID_SNAP_TO_GRID) pCmdUI->SetCheck(Snap_To_Grid);
	else if(nID == ID_SNAP_TO_JOINT) pCmdUI->SetCheck(Snap_To_Joint);
	else if(nID == ID_SNAP_TO_MEMBERMID) pCmdUI->SetCheck(Snap_To_MemberMid);
	else if(nID == ID_SNAP_TO_MEMBER) pCmdUI->SetCheck(Snap_To_Member);
}
void CmyDocument::FillCases() {
	BOOL first = TRUE;
	LOADCASE* loadcase;

	POSITION pos;
	LoadCaseListBox->ResetContent();
	pos = loadcases.GetHeadPosition();
	while(pos) {
		loadcase = &loadcases.GetNext(pos);
		if(first) {
			first = FALSE;
			c_AnalysisCase = loadcase;
		}
		LoadCaseListBox->AddString(loadcase->name);
	}
	LoadCaseListBox->SetCurSel(0);
}
void CmyDocument::FillAnalysisCases() {
	/*
	Add to list
	*/
	POSITION pos;
	PANALYSISCASELIST pcaselist;
	ANALYSISCASE* pcase;
	BOOL first = TRUE;

	LoadCaseListBox->ResetContent();
	for(UINT j = 0;j < CASETYPES;j++) {
		pcaselist = all_analysis_cases[j];
		
		pos = pcaselist->GetHeadPosition();
		while(pos) {
			pcase = &pcaselist->GetNext(pos);
			if(pcase->run && pcase->finished) {
				if(first) {
					first = FALSE;
                    c_AnalysisCase = pcase;
				}
				LoadCaseListBox->AddString(pcase->name);
			}
		}
	}
	/*reset*/
	LoadCaseListBox->SetCurSel(0);
	AssignIndex(c_AnalysisCase->index);
}
void CmyDocument::OnLock() { 
	if(Lock && AnalysisResult) {
		FillCases();
		
		int result = AfxMessageBox("Unlocking model will delete \nanalysis results",MB_OKCANCEL);
		if(result == IDCANCEL) 
			return;
		
		FreeAnalysis();
		
		AnalysisResult = FALSE;
		DesignResult = FALSE;

		RestoreObjects();

		CmyView* pView;
		POSITION pos = GetFirstViewPosition();
		while (pos != NULL) {
			pView = (CmyView*) GetNextView(pos);
			pView->SetTitle();
			pView->Invalidate();
		} 
	}
	Lock = !Lock;
}
void CmyDocument::FreeAnalysis(ENTITY* pentity) {
	POSITION pos1;
	PANALYSISCASELIST pcaselist;
	ANALYSISCASE* pcase;
	UINT i,j,repeats;

	
	for(j = 0;j < CASETYPES;j++) {
		pcaselist = all_analysis_cases[j];
		
		pos1 = pcaselist->GetHeadPosition();
		while(pos1) {
			pcase = &pcaselist->GetNext(pos1);
			if(pcase->finished) {
				switch(j) {
				case ANALYSISCASE::LOAD_CASE:
					repeats = 1;
					break;
				case ANALYSISCASE::MODAL_CASE:
					repeats = ((MODALCASE*)pcase)->runmodes;
					break;
				case ANALYSISCASE::RESPONSEH_CASE:
					repeats = ((RESPONSEHIST*)pcase)->N;
					break;
				case ANALYSISCASE::RESPONSES_CASE:
					repeats = 1;
					break;
				case ANALYSISCASE::COMBO_CASE:
					if(((COMBINATION*)pcase)->type == ENVELOPECOMBO) 
						repeats = 2;
					break;
				case ANALYSISCASE::NL_CASE:
					repeats = 1;
					break;
				case ANALYSISCASE::BUCKLING_CASE:
					repeats = ((BUCKLINGCASE*)pcase)->runmodes;
					break;
				}
				for(i = 0;i < repeats;i++) {
					pentity->FreeAnalysis(pcase->index + i);
				}
			}
		}
	}

	pentity->FreeVectors();
}

void CmyDocument::FreeAnalysis() { 
	
	POSITION pos;
	UINT ientity;
	ENTITY* pentity;
	PENTITYLIST pentities;

	for(ientity = 0;ientity < NENTITIES;ientity++) {
		pentities = entities[ientity];
		
		pos = pentities->GetHeadPosition();
		while(pos) {
			pentity = &pentities->GetNext(pos);
			FreeAnalysis(pentity);
		}
	}

	/*finish*/
	POSITION pos1;
	PANALYSISCASELIST pcaselist;
	ANALYSISCASE* pcase;
	UINT j;
	
	for(j = 0;j < CASETYPES;j++) {
		pcaselist = all_analysis_cases[j];
		
		pos1 = pcaselist->GetHeadPosition();
		while(pos1) {
			pcase = &pcaselist->GetNext(pos1);
			switch(j) {
			case ANALYSISCASE::MODAL_CASE:
				if(pcase->finished)
					vec_free(((MODALCASE*)pcase)->eigvalue);
				break;
			case ANALYSISCASE::BUCKLING_CASE:
				if(pcase->finished)
					vec_free(((BUCKLINGCASE*)pcase)->eigvalue);
				break;
			case ANALYSISCASE::NL_CASE:
				if(((NLCASE*)pcase)->end_stiffness)
					vec_free(((NLCASE*)pcase)->end_stiffness);
				break;
			}
			pcase->finished = FALSE;
		}
	}
}
void CmyDocument::OnUpdateLockUI(CCmdUI* pCmdUI) {
	pCmdUI->SetRadio(Lock);
}
void CmyDocument::OnUpdateLockMenuUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(!Lock);
}
void CmyDocument::ShowFile(int type) {
    CString sPath = "notepad.exe ";
	sPath += GetPathName();
	int len = sPath.GetLength();
	if(type == 0) {
		sPath.SetAt(len - 3,'A');
		sPath.SetAt(len - 2,'N');
		sPath.SetAt(len - 1,'A');
	} else {
		sPath.SetAt(len - 3,'D');
		sPath.SetAt(len - 2,'G');
		sPath.SetAt(len - 1,'N');
	}

	LPTSTR p = sPath.GetBuffer(0);
	STARTUPINFO si;
    PROCESS_INFORMATION pi;
    ZeroMemory( &si, sizeof(si));
    si.cb = sizeof(si);
    CreateProcess( NULL, p, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi);
}
void CmyDocument::OnShowResult() {
	ShowFile(0);
}
void CmyDocument::OnShowDesignResult() {
	ShowFile(1);
}
void CmyDocument::OnUpdateAnalysisResultUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(AnalysisResult);
}
void CmyDocument::OnUpdateDesignResultUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(DesignResult);
}
void CmyDocument::OnUpdateDesignUI(CCmdUI* pCmdUI) {
	BOOL enable = AnalysisResult;
	if(c_AnalysisCase->acase_type == ANALYSISCASE::MODAL_CASE)
		enable = FALSE;
	pCmdUI->Enable(enable);
}
void CmyDocument::AssignIndex(int index) {
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->forces = joint->forces_all[index];
		joint->disps = joint->disps_all[index]; 
	}
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->forces = member->forces_all[index];
		member->disps = member->disps_all[index];
	}
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->forces = slab->forces_all[index];
		slab->disps = slab->disps_all[index];
	}
}
void CmyDocument::OnSelectLoadCase() {
	LOADCASE ld;
	MODALCASE md;
	RESPONSEHIST rh;
	RESPONSESPEC rs;
	COMBINATION comb;
	NLCASE nlc;
	BUCKLINGCASE bc;
	POSITION pos;
	MODALCASE* pmodal;
	RESPONSEHIST* prh;
	COMBINATION* pcombo;
	BUCKLINGCASE* pbkl;
	int index;

    LoadCaseListBox->GetLBText(LoadCaseListBox->GetCurSel(), ld.name);
	md.name = ld.name;
	rh.name = ld.name;
    rs.name = ld.name;
	comb.name = ld.name;
	nlc.name = ld.name;
	bc.name = ld.name;

	if(pos = loadcases.Find(ld)) {
		c_AnalysisCase = &loadcases.GetAt(pos);
		index = c_AnalysisCase->index;
    } else if(pos = modalcases.Find(md)) {
		pmodal = &modalcases.GetAt(pos); 
		c_AnalysisCase = pmodal;
		CResultDia Dia(AfxGetMainWnd(),pmodal->runmodes,pmodal->selmode);
		if(Dia.DoModal() == IDOK) {
			pmodal->selmode = Dia.selected;
		}
		index = pmodal->index + pmodal->selmode;
    } else if(pos = responsecases.Find(rh)) {
		prh = &responsecases.GetAt(pos); 
		c_AnalysisCase = prh;
		CResultDia Dia(AfxGetMainWnd(),prh->N,prh->seltime);
		if(Dia.DoModal() == IDOK) {
			prh->seltime = Dia.selected;
		}
		index = prh->index + prh->seltime;
    } else if(pos = responsespecs.Find(rs)) {
		c_AnalysisCase = &responsespecs.GetAt(pos);
		index = c_AnalysisCase->index;
	} else if(pos = combinations.Find(comb)) {
        pcombo = &combinations.GetAt(pos);
		c_AnalysisCase = pcombo;
		if(pcombo->type == ENVELOPECOMBO) {
			CResultDia Dia(AfxGetMainWnd(),2,pcombo->selcombo,CResultDia::MINMAX);
			if(Dia.DoModal() == IDOK) {
				pcombo->selcombo = Dia.selected;
			}
		}
		index = pcombo->index + pcombo->selcombo;
	} else if(pos = nlcases.Find(nlc)) {
		c_AnalysisCase = &nlcases.GetAt(pos);
		index = c_AnalysisCase->index;
	} else if(pos = bucklingcases.Find(bc)) {
		pbkl = &bucklingcases.GetAt(pos); 
		c_AnalysisCase = pbkl;
		CResultDia Dia(AfxGetMainWnd(),pbkl->runmodes,pbkl->selmode);
		if(Dia.DoModal() == IDOK) {
			pbkl->selmode = Dia.selected;
		}
		index = pbkl->index + pbkl->selmode;
	} else {
		return;
	}
	
    if(AnalysisResult) {
		AssignIndex(index);
		DetermineScale();
		DesignResult = FALSE;
	}
	
	CmyView* pView;
	pos = GetFirstViewPosition();
	while (pos != NULL) {
		pView = (CmyView*) GetNextView(pos);
		pView->SetTitle();
		pView->Invalidate();
	} 
}
void CmyDocument::DetermineScale() {

	/*load scale*/
	DOUBLE maxload = 0;
	DOUBLE maxdisp = 0;
	MEMBER* member;
	SLAB* slab;
	LOAD* pload;
	UBMP32 i,j;

	
	for(j = IUX;j <= IRZ;j++) {
		force_scale[j] = 0;
	}
	
	POSITION pos = members.GetHeadPosition(),pos1;
	while(pos) {
		member = &members.GetNext(pos);

		/*load*/
	    pos1 = member->load.GetHeadPosition();
		while(pos1) {
           pload = &member->load.GetNext(pos1);
		   if(!pload->is_member_load())
			   continue;
		   if(fabs(pload->P) > maxload)  maxload = fabs(pload->P);
		   if(fabs(pload->P1) > maxload)  maxload = fabs(pload->P1);
		}
		/*force/displacment*/
		if(AnalysisResult) {
			DOUBLE H;
			for(i = 0;i < member->nDiv;i++) {
				for(j = IUX;j <= IRZ;j++) {
					H = fabs(member->forces[i][j]);
					if(H > force_scale[j]) force_scale[j] = H;
				}
			}
		}
		/*displacement*/
		if(AnalysisResult) {
			DOUBLE L = member->GetLength(),m;
			RPoint v;
			for(i = 0;i < member->nDiv;i++) {
				member->CalculateDeflection(member->station[i],v);
				m = v.magnitude();
				if(m > maxdisp) maxdisp = m;
			}
		}
	}

	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);

		/*load*/
	    pos1 = slab->load.GetHeadPosition();
		while(pos1) {
           pload = &slab->load.GetNext(pos1);
		   if(!pload->is_slab_load())
			   continue;
		   if(fabs(pload->P) > maxload)  maxload = fabs(pload->P);
		}
		
		/*displacement*/
		if(AnalysisResult) {
			DOUBLE m;
			for(j = 0;j < 6 * slab->NJ;j++) {
				m = fabs(slab->disps[j]);
				if(m > maxdisp) maxdisp = m;
			}
		}
	}

	/*joints*/
	JOINT* joint;
	JLOAD* pjload;
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		for(i = 0;i < 6;i++) {
			pos1 = joint->load.GetHeadPosition();
			while(pos1) {
                pjload = &joint->load.GetNext(pos1);
				if(!pjload->is_joint_load())
			         continue;
			    if(fabs(pjload->Q[i]) > maxload)  maxload = fabs(pjload->Q[i]);
			}
		}
	}
   
	if(maxload == 0) maxload  = 1;
	if(maxdisp == 0) maxdisp = 1;
	load_scale = 4 / maxload;
	disp_scale = 4 / maxdisp;

	for(j = IUX;j <= IRZ;j++) {
		if(force_scale[j]) force_scale[j] = 4 / force_scale[j];
		else force_scale[j] = 10;
	}
}

void CmyDocument::DetermineHistoryScale() {

	/*load scale*/
	DOUBLE maxdisp = 0;
	MEMBER* member;
	SLAB* slab;
	UBMP32 i,j,k;
	POSITION pos;
	RESPONSEHIST* pcase = (RESPONSEHIST*)c_AnalysisCase;

	
	for(j = IUX;j <= IRZ;j++) {
		force_scale[j] = 0;
	}
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		
		for(k = 0;k < pcase->N;k++) {
			/*force/displacment*/
			DOUBLE H;
			for(i = 0;i < member->nDiv;i++) {
				for(j = IUX;j <= IRZ;j++) {
					H = fabs(member->forces_all[pcase->index + k][i][j]);
					if(H > force_scale[j]) force_scale[j] = H;
				}
			}
			/*displacement*/
			DOUBLE L = member->GetLength(),m;
			RPoint v;
			for(i = 0;i < member->nDiv;i++) {
				member->CalculateDeflection(member->station[i],v,pcase->index + i);
				m = v.magnitude();
				if(m > maxdisp) maxdisp = m;
			}
		}
	}

	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		
		/*displacement*/
		if(AnalysisResult) {
			DOUBLE m;
			for(j = 0;j < 6 * slab->NJ;j++) {
				m = fabs(slab->disps[j]);
				if(m > maxdisp) maxdisp = m;
			}
		}
	}
   
	if(maxdisp == 0) maxdisp = 1;
	disp_scale = 3 / maxdisp;

	for(j = IUX;j <= IRZ;j++) {
		if(force_scale[j]) force_scale[j] = 6 / force_scale[j];
		else force_scale[j] = 10;
	}
}
ENTITY* CmyDocument::FindEntity(CString& name,int index) {
	ENTITY* entity;
	POSITION pos = entities[index]->GetHeadPosition();
	while(pos) {
		entity = &entities[index]->GetNext(pos);
		if(entity->name == name) {
			return entity;
		}
	}
	return NULL;

}
/*
combine loads
*/
void CmyDocument::CombineLoad(COMB_TYPE* combo) {
	register UINT index = combo->index,nindex,i,j;
	POSITION pos,pos1;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	LOADCOMBO* loadcombo;
	DOUBLE temp;
	int type = combo->type;
	if(combo->acase_type == ANALYSISCASE::BUCKLING_CASE)
		type = LINEARCOMBO;

	/*linear combination*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->AllocAnalysis(index);
		joint->forces = joint->forces_all[index];
		joint->disps = joint->disps_all[index]; 
		/*alloc one more*/
		if(type == ENVELOPECOMBO) {
			joint->AllocAnalysis(index + 1);
			for(i = 0;i < 6;i++) {
				joint->forces_all[index + 1][i] = 1e16;
				joint->disps_all[index + 1][i] = 1e16;
			}
		}

		pos1 = combo->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &combo->loadlist.GetNext(pos1);
			nindex = loadcombo->loadcase->index;

			for(i = 0;i < 6;i++) {
				if(type == LINEARCOMBO) {
					joint->forces[i] += (joint->forces_all[nindex][i] * loadcombo->FS);
					joint->disps[i] += (joint->disps_all[nindex][i] * loadcombo->FS);
				} else if(type == ABSSUMCOMBO) {
					joint->forces[i] += fabs(joint->forces_all[nindex][i] * loadcombo->FS);
					joint->disps[i] += fabs(joint->disps_all[nindex][i] * loadcombo->FS);
				} else if(type == SRSSCOMBO) {
					joint->forces[i] += pow(joint->forces_all[nindex][i] * loadcombo->FS , 2);
					joint->disps[i] += pow(joint->disps_all[nindex][i] * loadcombo->FS , 2);
				} else if(type == ENVELOPECOMBO) {
					temp = joint->forces_all[nindex][i] * loadcombo->FS;
					if(fabs(temp) > fabs(joint->forces[i])) 
						joint->forces[i] = temp;
					if(fabs(temp) < fabs(joint->forces_all[index + 1][i]))
						joint->forces_all[index + 1][i] = temp;

					temp = joint->disps_all[nindex][i] * loadcombo->FS;
					if(fabs(temp) > fabs(joint->disps[i]))
						joint->disps[i] = temp;
					if(fabs(temp) < fabs(joint->disps_all[index + 1][i]))
						joint->disps_all[index + 1][i] = temp;
				}
			}
		}
		/*SRSS*/
		if(type == SRSSCOMBO) {
			for(i = 0;i < 6;i++) {
				joint->forces[i] = sqrt(joint->forces[i]);
				joint->disps[i] = sqrt(joint->disps[i]);
			}
		}
	}
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->AllocAnalysis(index);
		member->forces = member->forces_all[index];
        member->disps = member->disps_all[index];

		/*alloc for others*/
		if(type == ENVELOPECOMBO)  {
			member->AllocAnalysis(index + 1);
			for(i = 0;i < member->nDiv; i++) {
				for(j = 0;j < 6;j++) {
                    member->forces_all[index + 1][i][j] = 1e16;
					member->disps_all[index + 1][i][j] = 1e16;
				}
			}
		}

		pos1 = combo->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &combo->loadlist.GetNext(pos1);
			nindex = loadcombo->loadcase->index;
			
			for(i = 0;i < member->nDiv; i++) {
				for(j = 0;j < 6;j++) {
					if(type == LINEARCOMBO) {
						member->forces[i][j] += (member->forces_all[nindex][i][j] * loadcombo->FS);
						member->disps[i][j] += (member->disps_all[nindex][i][j] * loadcombo->FS);
					} else if(type == ABSSUMCOMBO) {
						member->forces[i][j] += fabs(member->forces_all[nindex][i][j] * loadcombo->FS);
						member->disps[i][j] += fabs(member->disps_all[nindex][i][j] * loadcombo->FS);
					} else if(type == SRSSCOMBO) {
						member->forces[i][j] += pow(member->forces_all[nindex][i][j] * loadcombo->FS , 2);
						member->disps[i][j] += pow(member->disps_all[nindex][i][j] * loadcombo->FS , 2);
					} else if(type == ENVELOPECOMBO) {
						temp = member->forces_all[nindex][i][j] * loadcombo->FS ;
						if(fabs(temp) >  fabs(member->forces[i][j]))
							member->forces[i][j] = temp;
						if(fabs(temp) < fabs(member->forces_all[index + 1][i][j]))
							member->forces_all[index + 1][i][j] = temp;

						temp = member->disps_all[nindex][i][j] * loadcombo->FS ;
						if(fabs(temp) >  fabs(member->disps[i][j]))
							member->disps[i][j] = temp;
						if(fabs(temp) < fabs(member->disps_all[index + 1][i][j]))
							member->disps_all[index + 1][i][j] = temp;
					}
				}
			}
		}
		/*SRSS*/
		if(type == SRSSCOMBO) {
			for(i = 0;i < member->nDiv; i++) {
				for(j = 0;j < 6;j++) {
					member->forces[i][j] = sqrt(member->forces[i][j]);
					member->disps[i][j] = sqrt(member->disps[i][j]);
				}
			}
		}
	}
	
	
	UINT Total,Total1;
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->AllocAnalysis(index);
		slab->forces = slab->forces_all[index];
        slab->disps = slab->disps_all[index];
		
		Total = SLAB::NSTRESS * slab->NJ;
		Total1 = 6 * slab->NJ;
		
		/*alloc for others*/
		if(type == ENVELOPECOMBO)  {
			slab->AllocAnalysis(index + 1);
			for(j = 0;j < Total;j++) {
				slab->forces_all[index + 1][j] = 1e16;
			}
			for(j = 0;j < Total1;j++) {
				slab->disps_all[index + 1][j] = 1e16;
			}
		}
		
		pos1 = combo->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &combo->loadlist.GetNext(pos1);
			nindex = loadcombo->loadcase->index;
			
			for(j = 0;j < Total;j++) {
				if(type == LINEARCOMBO) {
					slab->forces[j] += (slab->forces_all[nindex][j] * loadcombo->FS);
				} else if(type == ABSSUMCOMBO) {
					slab->forces[j] += fabs(slab->forces_all[nindex][j] * loadcombo->FS);
				} else if(type == SRSSCOMBO) {
					slab->forces[j] += pow(slab->forces_all[nindex][j] * loadcombo->FS , 2);
				} else if(type == ENVELOPECOMBO) {
					temp = slab->forces_all[nindex][j] * loadcombo->FS ;
					if(fabs(temp) >  fabs(slab->forces[j]))
						slab->forces[j] = temp;
					if(fabs(temp) < fabs(slab->forces_all[index + 1][j]))
						slab->forces_all[index + 1][j] = temp;
				}
			}
			for(j = 0;j < Total1;j++) {
				if(type == LINEARCOMBO) {
					slab->disps[j] += (slab->disps_all[nindex][j] * loadcombo->FS);
				} else if(type == ABSSUMCOMBO) {
					slab->disps[j] += fabs(slab->disps_all[nindex][j] * loadcombo->FS);
				} else if(type == SRSSCOMBO) {
					slab->disps[j] += pow(slab->disps_all[nindex][j] * loadcombo->FS , 2);
				} else if(type == ENVELOPECOMBO) {
					temp = slab->disps_all[nindex][j] * loadcombo->FS ;
					if(fabs(temp) >  fabs(slab->disps[j]))
						slab->disps[j] = temp;
					if(fabs(temp) < fabs(slab->disps_all[index + 1][j]))
						slab->disps_all[index + 1][j] = temp;
				}
			}
		}
		/*SRSS*/
		if(type == SRSSCOMBO) {
			for(j = 0;j < Total;j++) {
				slab->forces[j] = sqrt(slab->forces[j]);
			}
			for(j = 0;j < Total1;j++) {
				slab->disps[j] = sqrt(slab->disps[j]);
			}
		}
	}
}
void CmyDocument::CombineModes(RESPONSESPEC* resp) {
	register UINT start,end,index;
	register UINT i,j,k,l,n,count;
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	DOUBLE bij,rij,ei = 0.05,ej = 0.05;

	/*Modal combination using ABSSUM/SRSS/CQC*/
	count = resp->funclist.GetCount();
	
	for(n = 0;n < count;n++) {
		
        index = resp->index + 1 + n * (resp->modalcase->runmodes + 1);
		start = index + 1;
		end = start + resp->modalcase->runmodes;
		
		/*
		JOINT
		*/
		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			joint->AllocAnalysis(index);
			joint->forces = joint->forces_all[index];
			joint->disps = joint->disps_all[index]; 
			
			if(resp->modal_comb == CQC) {
				for(j = start;j < end;j++) {
					for(i = start;i < end;i++) { 
						bij = sqrt(resp->modalcase->eigvalue[i - start]) / sqrt(resp->modalcase->eigvalue[j - start]);
						rij = 8 * sqrt(ei * ej) * (bij * ei + ej) * pow(bij,1.5);
						rij /= (pow(1 - bij * bij,2) + 4 * ei * ej * bij * (1 + bij * bij) + 4 * (ei * ei + ej * ej) * bij * bij);
						for(k = 0;k < 6;k++) {
							joint->forces[k] += (rij * joint->forces_all[i][k] * joint->forces_all[j][k]);
							joint->disps[k] += (rij * joint->disps_all[i][k] * joint->disps_all[j][k]);
						}
					}
				}
			} else {
				for(i = start;i < end;i++) {
					for(j = 0;j < 6;j++) {
						if(resp->modal_comb == ABSSUM) {
							joint->forces[j] += fabs(joint->forces_all[i][j]);
							joint->disps[j] += fabs(joint->disps_all[i][j]);
						} else if(resp->modal_comb == SRSS) {
							joint->forces[j] += pow(joint->forces_all[i][j], 2);
							joint->disps[j] += pow(joint->disps_all[i][j], 2);
						}
					}
				}
			}
			/*SRSS and CQC*/
			if(resp->modal_comb == SRSS || resp->modal_comb == CQC) {
				for(j = 0;j < 6;j++) {
					joint->forces[j] = sqrt(fabs(joint->forces[j]));
					joint->disps[j] = sqrt(fabs(joint->disps[j]));
				}
			}
		}
		/*
		MEMBER
		*/
		pos = members.GetHeadPosition();
		while(pos) {
			member = &members.GetNext(pos);
			member->AllocAnalysis(index);
			member->forces = member->forces_all[index];
			member->disps = member->disps_all[index];
						
			if(resp->modal_comb == CQC) {
                for(i = start;i < end;i++) { 
					for(j = start;j < end;j++) {
						bij = sqrt(resp->modalcase->eigvalue[i - start]) / sqrt(resp->modalcase->eigvalue[j - start]);
						rij = 8 * sqrt(ei * ej) * (bij * ei + ej) * pow(bij,1.5);
						rij /= (pow(1 - bij * bij,2) + 4 * ei * ej * bij * (1 + bij * bij) + 4 * (ei * ei + ej * ej) * bij * bij);
						
						for(k = 0;k < 6;k++) {
							for(l = 0;l < member->nDiv; l++) {
								member->forces[l][k] += (rij * member->forces_all[i][l][k] * member->forces_all[j][l][k]);
								member->disps[l][k] += (rij * member->disps_all[i][l][k] * member->disps_all[j][l][k]);
							}
						}
					}
				}
			} else {
				for(i = start;i < end;i++) {
					for(j = 0;j < 6;j++) {
						for(k = 0;k < member->nDiv; k++) {
							if(resp->modal_comb == ABSSUM) {
								member->forces[k][j] += fabs(member->forces_all[i][k][j]);
								member->disps[k][j] += fabs(member->disps_all[i][k][j]);
							} else if(resp->modal_comb == SRSS) {
								member->forces[k][j] += pow(member->forces_all[i][k][j],2);
								member->disps[k][j] += pow(member->disps_all[i][k][j],2);
							}
						}
					}
				}
			}
			/*SRSS and CQC*/
			if(resp->modal_comb == SRSS || resp->modal_comb == CQC) {
				for(j = 0;j < 6;j++) {
					for(k = 0;k < member->nDiv; k++) {
						member->forces[k][j] = sqrt(fabs(member->forces[k][j]));
						member->disps[k][j] = sqrt(fabs(member->disps[k][j]));
					}
				}
			}
		}
        /*
		SLAB
		*/
		UINT Total,Total1;
		pos = slabs.GetHeadPosition();
		while(pos) {
			slab = &slabs.GetNext(pos);
			slab->AllocAnalysis(index);
			slab->forces = slab->forces_all[index];
			slab->disps = slab->disps_all[index];
				
			Total = SLAB::NSTRESS * slab->NJ;
			Total1 = 6 * slab->NJ;

			if(resp->modal_comb == CQC) {
                for(i = start;i < end;i++) { 
					for(j = start;j < end;j++) {
						bij = sqrt(resp->modalcase->eigvalue[i - start]) / sqrt(resp->modalcase->eigvalue[j - start]);
						rij = 8 * sqrt(ei * ej) * (bij * ei + ej) * pow(bij,1.5);
						rij /= (pow(1 - bij * bij,2) + 4 * ei * ej * bij * (1 + bij * bij) + 4 * (ei * ei + ej * ej) * bij * bij);
						
						for(k = 0;k < Total;k++) {
							slab->forces[k] += (rij * slab->forces_all[i][k] * slab->forces_all[j][k]);
						}
						for(k = 0;k < Total1;k++) {
							slab->disps[k] += (rij * slab->disps_all[i][k] * slab->disps_all[j][k]);
						}
					}
				}
			} else {
				for(i = start;i < end;i++) {
					for(j = 0;j < Total;j++) {
						if(resp->modal_comb == ABSSUM) {
							slab->forces[j] += fabs(slab->forces_all[i][j]);
						} else if(resp->modal_comb == SRSS) {
							slab->forces[j] += pow(slab->forces_all[i][j],2);
						}
					}
					for(j = 0;j < Total1;j++) {
						if(resp->modal_comb == ABSSUM) {
							slab->disps[j] += fabs(slab->disps_all[i][j]);
						} else if(resp->modal_comb == SRSS) {
							slab->disps[j] += pow(slab->disps_all[i][j],2);
						}
					}
				}
			}
			/*SRSS and CQC*/
			if(resp->modal_comb == SRSS || resp->modal_comb == CQC) {
				for(j = 0;j < Total;j++) {
					slab->forces[j] = sqrt(fabs(slab->forces[j]));
				}
				for(j = 0;j < Total1;j++) {
					slab->disps[j] = sqrt(fabs(slab->disps[j]));
				}
			}
		}
	}

	
	/*directional combination*/
	count = resp->funclist.GetCount();
	index = resp->index;
	
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->AllocAnalysis(index);
		joint->forces = joint->forces_all[index];
		joint->disps = joint->disps_all[index]; 
		
		for(n = 0;n < count;n++) {
			i = resp->index + 1 + n * (resp->modalcase->runmodes + 1);
			for(j = 0;j < 6;j++) {
				if(resp->dir_comb == ABSSUM) {
					joint->forces[j] += fabs(joint->forces_all[i][j]);
					joint->disps[j] += fabs(joint->disps_all[i][j]);
				} else if(resp->dir_comb == SRSS) {
					joint->forces[j] += pow(joint->forces_all[i][j], 2);
					joint->disps[j] += pow(joint->disps_all[i][j], 2);
				}
			}
		}
		/*SRSS*/
		if(resp->dir_comb == SRSS) {
			for(j = 0;j < 6;j++) {
				joint->forces[j] = sqrt(joint->forces[j]);
				joint->disps[j] = sqrt(joint->disps[j]);
			}
		}
		/*free memory*/
		for(i = resp->index + 1;i < end;i++) {
			joint->FreeAnalysis(i);
		}
	}
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->AllocAnalysis(index);
		member->forces = member->forces_all[index];
		member->disps = member->disps_all[index];
			
		for(n = 0;n < count;n++) {
			i = resp->index + 1 + n * (resp->modalcase->runmodes + 1);
					
			for(j = 0;j < 6;j++) {
				for(k = 0;k < member->nDiv; k++) {
					if(resp->dir_comb == ABSSUM) {
						member->forces[k][j] += fabs(member->forces_all[i][k][j]);
						member->disps[k][j] += fabs(member->disps_all[i][k][j]);
					} else if(resp->dir_comb == SRSS) {
						member->forces[k][j] += pow(member->forces_all[i][k][j],2);
						member->disps[k][j] += pow(member->disps_all[i][k][j],2);
					}
				}
			}
		}
		/*SRSS*/
		if(resp->dir_comb == SRSS) {
			for(j = 0;j < 6;j++) {
				for(k = 0;k < member->nDiv; k++) {
					member->forces[k][j] = sqrt(member->forces[k][j]);
					member->disps[k][j] = sqrt(member->disps[k][j]);
				}
			}
		}
		/*free memory*/
		for(i = resp->index + 1;i < end;i++) {
			member->FreeAnalysis(i);
		}
	}

	UINT Total,Total1;
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->AllocAnalysis(index);
		slab->forces = slab->forces_all[index];
		slab->disps = slab->disps_all[index];

		Total = SLAB::NSTRESS * slab->NJ;
		Total1 = 6 * slab->NJ;
			
		for(n = 0;n < count;n++) {
			i = resp->index + 1 + n * (resp->modalcase->runmodes + 1);
			
			for(j = 0;j < Total;j++) {
				if(resp->dir_comb == ABSSUM) {
					slab->forces[j] += fabs(slab->forces_all[i][j]);
				} else if(resp->dir_comb == SRSS) {
					slab->forces[j] += pow(slab->forces_all[i][j],2);
				}
			}
			for(j = 0;j < Total1;j++) {
				if(resp->dir_comb == ABSSUM) {
					slab->disps[j] += fabs(slab->disps_all[i][j]);
				} else if(resp->dir_comb == SRSS) {
					slab->disps[j] += pow(slab->disps_all[i][j],2);
				}
			}
		}
		/*SRSS*/
		if(resp->dir_comb == SRSS) {
			for(j = 0;j < Total;j++) {
				slab->forces[j] = sqrt(slab->forces[j]);
			}
			for(j = 0;j < Total1;j++) {
				slab->disps[j] = sqrt(slab->disps[j]);
			}
		}
		/*free memory*/
		for(i = resp->index + 1;i < end;i++) {
			slab->FreeAnalysis(i);
		}
	}
}
void CmyDocument::LumpMasses() {
	DOUBLE m,em,L;
	POSITION pos,pos1;
	SLAB* slab;
	MEMBER* member;
	JOINT* joint;
	JLOAD* jload;
	UINT i;

	pos = joints.GetHeadPosition();
	while(pos) {
        joint = &joints.GetNext(pos);
		joint->mself = 0;
		for(i = IUX;i <= IRZ;i++) {
			joint->massembled[i] = 0;
		}
	}

    pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		L = member->GetLength();
		if(member->e_section) {
            DOUBLE d1,d2,A1,A2;
			d1 = member->section->material->density + (member->e_section->material->density - member->section->material->density) * member->start_offset;
            d2 = member->e_section->material->density + (member->section->material->density - member->e_section->material->density) * member->end_offset;;
            A1 = member->section->A + (member->e_section->A - member->section->A) * member->start_offset;
			A2 = member->e_section->A + (member->section->A - member->e_section->A) * member->end_offset;
			m  = ((3 * d1 + d2) / 4 ) * ((3 * A1 + A2) / 4 ) * (L / 2);
            em = ((3 * d2 + d1) / 4 ) * ((3 * A2 + A1) / 4 ) * (L / 2);
		} else {
			m = em = (member->section->material->density * member->section->A * (L / 2));
		}
		member->j1->mself += m;
		member->j2->mself += em;
	}

	DOUBLE Factors[9];
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->CalculateUniformFactors(Factors,SLAB::MASS_F);
		m = (slab->section->material->density * slab->section->h * slab->GetArea());
		for(i = 0;i < slab->NJ; i++) {
			slab->jt[i]->mself += (m * Factors[i]); 
		}
	}

	/*assemble masses
	TODO: global/local axis mass
	*/
	pos = joints.GetHeadPosition();
	while(pos) {
        joint = &joints.GetNext(pos);
		for(i = IUX;i <= IRZ;i++) {
			if(dynamic_dofs & (1 << i))
		       joint->massembled[i] += joint->mself;
		}
		
		pos1 = joint->load.GetHeadPosition();
		while(pos1) {
			jload = &joint->load.GetNext(pos1);
			if(jload->type == MASS) {
				for(i = IUX;i <= IRZ;i++) {
					joint->massembled[i] += jload->Q[i];
				}
			}
		} 
	}
}
void CmyDocument::AllocMemory(int index) {
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;

	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->AllocAnalysis(index);
		joint->forces = joint->forces_all[index];
		joint->disps = joint->disps_all[index]; 
	}
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->AllocAnalysis(index);
		member->forces = member->forces_all[index];
		member->disps = member->disps_all[index];
	}

	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->AllocAnalysis(index);
		slab->forces = slab->forces_all[index];
		slab->disps = slab->disps_all[index];
	}
}
void CmyDocument::FreeMemory(int index) {
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;

	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->FreeAnalysis(index);
	}
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->FreeAnalysis(index);
	}

	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->FreeAnalysis(index);
	}
}
void CmyDocument::PrintOutput() {
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	UINT i,j,k;

	/*joint displacement and reactions*/
	print("\n\t\t\t\t\t\tJOINT RESULT\n\t\t\t\t\t\t============\n");
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		print("\nJOINT %s\n=========\n",joint->name);
		for(int i = 0;i < 6;i++) {
			print("%3s %15.6g         %3s %15.6f\n",DirString[i % 6],joint->disps[i],ForceString[i % 6],joint->forces[i]);
		}
	}
	
	/*member displacement and forces*/
	print("\n\t\t\t\t\t\tMEMBER RESULT\n\t\t\t\t\t\t=============\n");
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		print("\n\t\t\tMember %s\n\t\t\t==========\n",member->name);
		print("\n%15s%15s%15s%15s%15s%15s%15s\n","X","VX","VY","VZ","MX","MY","MZ");
		for(i = 0;i < member->nDiv;i++) {
			print("%15.6f",member->station[i]);
			for(j = 0;j < 6;j++) {
				print("%15.6f",member->forces[i][j]);
			}
			print("\n");
		}
	}
	/*slab stresses*/
	print("\n\t\t\t\t\t\tSLAB RESULT\n\t\t\t\t\t\t=============\n");
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		print("\n\t\t\tSlab %s\n\t\t\t==========\n",slab->name);
		print("\n%15s","ID");
		for(k = 0;k < 23;k++) {
			print("%15s",StressString[k]);
		}
		print("\n");
		for(j = 0;j < slab->NJ;j++) {
			print("%15d",j);
			for(k = 0;k < 23;k++) {
				print("%15.6f",slab->forces[SLAB::NSTRESS * j + k]);
			}
			print("\n");
		}
		print("\n");
	}
}
void CmyDocument::CalculateForces(VECTOR Q,VECTOR D,UINT N,int index,BOOL add_to_prev) {
	POSITION pos;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;

	/*truncate*/
	if(!add_to_prev) {
		for(UINT i = 0;i < N;i++) {
			if(VEQUAL(Q[i],0)) Q[i] = 0;
			if(VEQUAL(D[i],0)) D[i] = 0;
		}
	}

	/*joint displacement and forces*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		for(int i = 0;i < 6;i++) {
			if(joint->number[i] != INVALID) {
				if(add_to_prev) {
					joint->forces[i] += Q[joint->number[i]];
					joint->disps[i] += D[joint->number[i]];
				} else {
					joint->forces[i] = Q[joint->number[i]];
					joint->disps[i] = D[joint->number[i]];
				}
			}
		}
	}

	/*member displacement and forces*/
	VECTOR q,d;
	vec_alloc(q,12);
	vec_alloc(d,12);

	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->CalcMF(Q,D,q,d);
		member->CalculateInternalForceAndDisplacement(q,d,add_to_prev);
	}

	vec_free(q);
	vec_free(d);

	/*slab stresses*/
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->CalculateStresses(Q,D);
	}
}
void CmyDocument::OnAnalysisOptions() {
	CAnalysisOptionsDia Dia(AfxGetMainWnd(),&static_dofs,&dynamic_dofs);
	if(Dia.DoModal() == IDOK)
		UpdateDrawing();
}
void CmyDocument::OnDesignOptions() {
	CDesignOptionsDia Dia(AfxGetMainWnd(),&detailing);
	Dia.DoModal();
}
void CmyDocument::OnDesignVerify() {
	POSITION pos;
	MEMBER* member;
	UINT mcount,i;

	mcount = 0;
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		for(i = 0;i < member->nDiv;i++) {
			if(member->design_data[i].flags) {
				mcount++;
				break;
			}
		}
	}

	if(mcount) {
		CString str = "Design of ";
		if(mcount) str.Format("%s%d members ",str,mcount);
		str += "failed. \nDo you want to select them?";

		if(AfxMessageBox(str,MB_OKCANCEL) == IDOK) {
			if(mcount) {
				pos = members.GetHeadPosition();
				while(pos) {
					member = &members.GetNext(pos);
					member->sel = 0;
					for(i = 0;i < member->nDiv;i++) {
						if(member->design_data[i].flags) {
							member->sel = true;
							break;
						}
					}
				}
			}
			UpdateDrawing();
		}
	} else {
		AfxMessageBox("All passed!");
	}
}
void CmyDocument::OnMiscOptions() {
	CMiscOptionsDia Dia(AfxGetMainWnd());
	if(Dia.DoModal() == IDOK)
		UpdateDrawing();
}
void CmyDocument::OnPlotFunctions() {
	CResultPlotDia Dia(AfxGetMainWnd());
	Dia.DoModal();
}
void CmyDocument::OnUpdatePlotFunctionsUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(AnalysisResult && (c_AnalysisCase->acase_type == ANALYSISCASE::RESPONSEH_CASE));
}
void CmyDocument::GetGridInfo(GRIDINFO* ginfo) {
	GRIDLIST* grids = global->grid;

	ginfo->Reset();

	if(!grids[IUX].IsEmpty()) {
		ginfo->xmin = grids[IUX].GetHead();
		ginfo->xmax = grids[IUX].GetTail();
	}
	if(!grids[IUY].IsEmpty()) {
		ginfo->ymin = grids[IUY].GetHead();
		ginfo->ymax = grids[IUY].GetTail();
	}
	if(!grids[IUZ].IsEmpty()) {
		ginfo->zmin = grids[IUZ].GetHead();
		ginfo->zmax = grids[IUZ].GetTail();
	}

}
void CmyDocument::CalculateEffectiveLength(MEMBER* target) {
	int index,i;
	DOUBLE c[2][2],b[2][2];
	DOUBLE Lt;
	DOUBLE restraint[12];
	MATRIX R;
	mat_alloc(R,3,3);
	target->GetRotMatrix(R,0);

	Lt = sqrt(pow(target->j2->p.x - target->j1->p.x,2) + 
	          pow(target->j2->p.y - target->j1->p.y,2) +
		      pow(target->j2->p.z - target->j1->p.z,2));
		
	c[0][0] = c[1][0] = (target->section->material->E * target->section->Iy) / Lt;
	c[0][1] = c[1][1] = (target->section->material->E * target->section->Iz) / Lt;
    b[0][0] = b[1][0] =	b[0][1] = b[1][1] = 0;


    MEMBER* member;
	JOINT *fjoint;
	POSITION pos;
	UINT INC;

	pos = members.GetHeadPosition();
	while(pos) {
        member = &members.GetNext(pos);
		if(member == target) continue;
		member->GetRestraint(restraint);

		if(target->j1 == member->j1 
			|| target->j1 == member->j2
			|| target->j2 == member->j1 
			|| target->j2 == member->j2) {

			if(target->j1 == member->j1) { fjoint = member->j2; index = 0;}
			else if(target->j1 == member->j2) { fjoint = member->j1; index = 0;}
			else if(target->j2 == member->j1) { fjoint = member->j2; index = 1;}
			else if(target->j2 == member->j2) { fjoint = member->j1; index = 1;}
			INC = (fjoint == member->j1) ? 0 : 6;

			VECTOR mem,vec;
			vec_alloc(mem,3);
			vec_alloc(vec,3);
			mem[0] = member->j2->p.x - member->j1->p.x;
			mem[1] = member->j2->p.y - member->j1->p.y;
			mem[2] = member->j2->p.z - member->j1->p.z;
			multiply(R,mem,vec,3,3);

			if(vec[0]) {
				c[index][0] += (member->section->material->E * member->section->Iy) / fabs(vec[0]);
				c[index][1] += (member->section->material->E * member->section->Iz) / fabs(vec[0]);
			} 
			if(vec[1]) {
				if((fjoint->mconnect > 1) || restraint[IRY + INC])
					b[index][0] += 1.0 * (member->section->material->E * member->section->Iy) / fabs(vec[1]);
				else if(restraint[IUZ + INC])
					b[index][0] += 0.5 * (member->section->material->E * member->section->Iy) / fabs(vec[1]);
			}
			if(vec[2]) {
				if((fjoint->mconnect > 1) || restraint[IRZ + INC])
					b[index][1] += 1.0 * (member->section->material->E * member->section->Iz) / fabs(vec[2]);
				else if(restraint[IUY + INC])
					b[index][1] += 0.5 * (member->section->material->E * member->section->Iz) / fabs(vec[2]);
			}

			vec_free(mem);
			vec_free(vec);
		}
	}

	DOUBLE Nsd,Ncr;
	DOUBLE a1,a2,z1,z2,am,g,tLe;
	BOOL swaymode,is_concrete;
	
	is_concrete = target->section->material->type == MATERIAL::CONCRETE;
	Nsd = max(fabs(target->forces[0][IUX]),fabs(target->forces[target->nDiv - 1][IUX]));

	/*
	Determine effecive length
	*/
#define MAXA 50

	target->GetRestraint(restraint);
		
	for(i = 0;i < 2;i++) {
		/*determine a1,a2,z1,z2*/
		if(b[0][i] == 0) a1 = MAXA;
		else {
			a1 = c[0][i] / b[0][i];
			if(a1 > MAXA) a1 = MAXA;
		}
		if(b[1][i] == 0) a2 = MAXA;
		else {
			a2 = c[1][i] / b[1][i];
			if(a2 > MAXA) a2 = MAXA;
		}

		swaymode = TRUE;
		if(i == 0) {
            if(restraint[IRY]) a1 = 0;
            if(restraint[6 + IRY]) a2 = 0;
			if(restraint[IUZ] && restraint[6 + IUZ]) swaymode = FALSE;
		} else {
			if(restraint[IRZ]) a1 = 0;
            if(restraint[6 + IRZ]) a2 = 0;
			if(restraint[IUY] && restraint[6 + IUY]) swaymode = FALSE;
		}
		
		z1 = a1 / (1 + a1);
		z2 = a2 / (1 + a2);

		if(swaymode) {
			if(i == 0) Ncr = PI * PI * target->section->material->E * target->section->Iy / pow(Lt, 2);
			else Ncr = PI * PI * target->section->material->E * target->section->Iz / pow(Lt, 2);
			if(is_concrete) Ncr *= 0.4;
			swaymode = ((Nsd / Ncr) > 0.1);
		}

		if(swaymode) {
			if((member->frame_type == DETAILING::NON_SWAY) ||
			   (member->frame_type == DETAILING::SPREF && detailing.frame_type == DETAILING::NON_SWAY) )
			     swaymode = FALSE;
		}
		
		if(!is_concrete) {
			if(swaymode) {
				g = sqrt((1 - 0.2 * (z1 + z2) - 0.12 * z1 * z2) / (1 - 0.8 * (z1 + z2) + 0.6 * z1 * z2));
				g = max(g, 1.15);
				tLe = g * Lt;
			} else {
				g = (1 + 0.145 * (z1 + z2) - 0.265 * z1 * z2) / (2 - 0.364 * (z1 + z2) - 0.247 * z1 * z2);
				g = max(g, 0.7);
				tLe = g * Lt;
			}
		} else {
			if(swaymode) {
				g = sqrt((7.5 + 4 * (a1 + a2) + 1.6 * a1 * a2) / (7.5 + a1 + a2));
				g = max(g, 1.15);
				tLe = g * Lt;
			} else {
				am = (a1 + a2) / 2;
				g = (am + 0.4) / (am + 0.8);
				g = max(g, 0.7);
				tLe = g * Lt;
			}
		}

        target->LeCalculated[i] = tLe;
	}
#undef MAXA

	vec_free(R);
}
/*
SLAB
*/
void CmyDocument::ApplyDeadLoad() {

	/*variables*/
	LOAD load,load1;
	MEMBER* member;
	SLAB* slab;
	DOUBLE L;
	POSITION pos,pos1;
	DOUBLE d1,d2,A1,A2;

	/*member self weight*/
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		L = member->GetLength();        
		if(member->e_section) {
			d1 = member->section->material->unitweight + (member->e_section->material->unitweight - member->section->material->unitweight) * member->start_offset;
            d2 = member->e_section->material->unitweight + (member->section->material->unitweight - member->e_section->material->unitweight) * member->end_offset;;
            A1 = member->section->A + (member->e_section->A - member->section->A) * member->start_offset;
			A2 = member->e_section->A + (member->section->A - member->e_section->A) * member->end_offset;
			load.type = TRAPEZOIDAL;
			load.x = 0;
		    load.P = -((3 * d1 + d2) / 4 ) * A1;
            load.x1 = L;
		    load.P1 = -((3 * d2 + d1) / 4 ) * A2;
		} else {
			load.type = UNIFORM;
		    load.P = -member->section->A * member->section->material->unitweight;
		}
		load.system = global;
		load.dir = IUZ;
		load.casetype = TEMPORARY;

		LOADCASE* pldcase;
		pos1 = loadcases.GetHeadPosition();
		while(pos1) {
			pldcase = &loadcases.GetNext(pos1);
			if(!pldcase->swm)
				continue;
			load1 = load;
			load1.P = pldcase->swm * load.P;
            load1.P1 = pldcase->swm * load.P1;
			load1.loadcase = pldcase;
			member->load.AddTail(load1);
		}
	}

	/*slab self weight*/
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
       
		load.type = UNIFORM;
        load.system = global;
		load.dir = IUZ;
		load.P = -slab->section->h * slab->section->material->unitweight;
		load.casetype = TEMPORARY;

		LOADCASE* pldcase;
		pos1 = loadcases.GetHeadPosition();
		while(pos1) {
			pldcase = &loadcases.GetNext(pos1);
			if(!pldcase->swm)
				continue;
			load1 = load;
			load1.P = pldcase->swm * load.P;
            load1.P1 = pldcase->swm * load.P1;
			load1.loadcase = pldcase;
			slab->load.AddTail(load1);
		}
	}
}
void CmyDocument::RemoveTemporaryLoads() {
	MEMBER* member;
	SLAB* slab;
	LOAD* pload;
	POSITION pos,pos1,post;

	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		pos1 = member->load.GetHeadPosition();
		while(pos1) {
			post = pos1;
			pload = &member->load.GetNext(pos1);
			if(pload->casetype == TEMPORARY)
				member->load.RemoveAt(post);
		}
	}

	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		pos1 = slab->load.GetHeadPosition();
		while(pos1) {
			post = pos1;
			pload = &slab->load.GetNext(pos1);
			if(pload->casetype == TEMPORARY)
				slab->load.RemoveAt(post);
		}
	}
}
/*
generate database and diplay it
*/
void CmyDocument::OnShowTables() {
	myTables.Clear();
	clear_tables();
	if(!Export(0,&myTables))
		return;
	CTableViewDia TableDia(AfxGetMainWnd(),&myTables);
	TableDia.DoModal();
}
/*
Display Options Dialog
*/
CDisplayOptionsDia::CDisplayOptionsDia(CWnd* parent,DISPLAY* pdisplay) : 
                    CDialog(IDD_DISPLAY,parent) {
	mydisplay = *pdisplay;
	ApplyToAll = FALSE;
}
void CDisplayOptionsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_DJOINT_LABEL, mydisplay.JointLabel);
	DDX_Check(pDX, IDC_DJOINT_RESTRAINT, mydisplay.JointRestraint);
	DDX_Check(pDX, IDC_DJOINT_SHOW, mydisplay.JointShow);
	DDX_Check(pDX, IDC_DJOINT_LOAD, mydisplay.JointLoad);
	DDX_Check(pDX, IDC_DJOINT_AXIS, mydisplay.JointAxis);
	DDX_Check(pDX, IDC_DMEMBER_LABEL, mydisplay.MemberLabel);
	DDX_Check(pDX, IDC_DMEMBER_SECTION, mydisplay.MemberSection);
	DDX_Check(pDX, IDC_DMEMBER_MATERIAL, mydisplay.MemberMaterial);
	DDX_Check(pDX, IDC_DMEMBER_LOAD, mydisplay.MemberLoad);
	DDX_Check(pDX, IDC_DMEMBER_AXIS, mydisplay.MemberAxis);
	DDX_Check(pDX, IDC_DMEMBER_RELEASE, mydisplay.MemberRelease);
	DDX_Check(pDX, IDC_DSLAB_LABEL, mydisplay.SlabLabel);
	DDX_Check(pDX, IDC_DSLAB_SECTION, mydisplay.SlabSection);
	DDX_Check(pDX, IDC_DSLAB_MATERIAL, mydisplay.SlabMaterial);
	DDX_Check(pDX, IDC_DSLAB_LOAD, mydisplay.SlabLoad);
	DDX_Check(pDX, IDC_DSLAB_AXIS, mydisplay.SlabAxis);
    DDX_Check(pDX, IDC_DORIGIN, mydisplay.Origin);
	DDX_Check(pDX, IDC_DGRID, mydisplay.Grid);
	DDX_Check(pDX, IDC_DCONSTRAINT, mydisplay.Constraints);
	DDX_Check(pDX, IDC_DNUMBERING, mydisplay.Numbering);
    DDX_Check(pDX, IDC_DEXTRUDE, mydisplay.Extrude);
    DDX_Check(pDX, IDC_APPLY_TO_ALL, ApplyToAll);
}
/*
Color Options Dialog
*/
CColorsDia::CColorsDia(CWnd* parent) : 
        CDialog(IDD_COLOR,parent) ,
		c_List1(CGridCtrl::TYPE3),
		c_List2(CGridCtrl::TYPE3),
		c_List3(CGridCtrl::TYPE3)
		{
}
void CColorsDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST1, c_List1);
    DDX_Control(pDX, IDC_LIST2, c_List2);
	DDX_Control(pDX, IDC_LIST3, c_List3);
}
BOOL CColorsDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	
	c_List1.SetTwoHeaders();
	c_List2.SetTwoHeaders();
	c_List3.SetTwoHeaders();

	CString str;
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::BackColor,"Background");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::JointColor,"Joint");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::MemberColor,"Member");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::SlabColor,"Slab");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::HighlightColor,"Highlight");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::SelectionColor,"Selection");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::RestraintColor,"Restraint");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::ConstraintColor,"Constraint");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::AxisColor,"Axis");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::ReleaseColor,"Release");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::TextColor,"Text");
	c_List1.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::GridColor,"Grid");
	for(UINT i = 0;i < CmyView::NCONTORS;i++) {
		str.Format("%d",i+1);
		c_List2.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::ContorColors[i],str);
	}
	c_List3.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::PositiveColor,"Positive");
	c_List3.InsertMyItem(GRIDCTRLDATA::COLOR,GRIDCTRLDATA::TINT,&CmyView::NegativeColor,"Negative");
	return TRUE;
}
void CColorsDia::OnOK() {
	c_List1.GetData();
	c_List2.GetData();
	c_List3.GetData();
	CDialog::OnOK();
}
/*
View class implementation
*/
IMPLEMENT_DYNCREATE(CmyView,CView)

BEGIN_MESSAGE_MAP(CmyView,CView)
    ON_WM_CREATE()
    ON_WM_SIZE()
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
    ON_WM_MOUSEWHEEL()
	ON_WM_SETCURSOR()
    ON_WM_SETFOCUS()
    ON_WM_KILLFOCUS()
	ON_WM_ERASEBKGND()
	ON_MESSAGE(WM_MOUSELEAVE, OnMouseLeave)
	ON_COMMAND(IDM_ZOOMIN,OnZoomIn)
    ON_COMMAND(IDM_ZOOMOUT,OnZoomOut)
	ON_COMMAND(IDM_RUBBERZOOM,OnRubberZoom)
    ON_COMMAND(IDM_RESTOREVIEW,OnRestoreView)
    ON_COMMAND(IDM_PAN,OnPan)
	ON_COMMAND(IDM_DISPLAY_OPTIONS,OnDisplayOptions)
	ON_COMMAND(IDM_SHOWGRIDS,OnShowGrids)
	ON_COMMAND(IDM_SHOWORIGIN,OnShowOrigin) 
	ON_COMMAND_RANGE(ID_VIEW_XY,ID_VIEW_3D,OnChangeView)
	ON_UPDATE_COMMAND_UI(ID_COORDINATE, OnCoordinate)
	ON_COMMAND(ID_ESCAPE,OnEscape)
	ON_COMMAND_RANGE(IDM_SELECT,IDM_LINE_SELECT,OnAdd)
    ON_UPDATE_COMMAND_UI_RANGE(IDM_SELECT,IDM_LINE_SELECT,OnUpdateAddUI)
	ON_COMMAND_RANGE(IDM_MOVEVIEWUP,IDM_MOVEVIEWDOWN,OnMoveView)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_MOVEVIEWUP,IDM_MOVEVIEWDOWN,OnUpdateMoveViewUI)
	ON_COMMAND(IDM_DIAGRAMS,OnDiagrams)
    ON_COMMAND(IDM_SLAB_DIAGRAMS,OnSlabDiagrams)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_DIAGRAMS,IDM_SHOW_LOADS,OnUpdateDiagramsUI)
    ON_UPDATE_COMMAND_UI(IDM_SLAB_DIAGRAMS,OnUpdateDiagramsUI)
	ON_COMMAND(IDM_DEFLECTED,OnDeflected)
	ON_COMMAND(IDM_DISPLACEMENT,OnDisplacement)
	ON_COMMAND(IDM_REACTIONS,OnReactions)
	ON_COMMAND(IDM_SHOW_LOADS,OnShowLoads)
	ON_COMMAND(IDM_PERSPECTIVE,OnPerspective)
	ON_COMMAND_RANGE(IDM_RIGHT,IDM_DOWN,OnRotate)
    ON_UPDATE_COMMAND_UI_RANGE(IDM_RIGHT,IDM_DOWN,OnUpdateRotateUI)
    ON_COMMAND(IDM_REFRESH,OnRefresh)
	ON_COMMAND(IDM_ANIMATE,OnAnimate)
	ON_UPDATE_COMMAND_UI(IDM_ANIMATE,OnUpdateAnimateUI)
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
END_MESSAGE_MAP()
/*
colors
*/
UINT CmyView::nViews = 0;
COLORREF CmyView::BackColor;
COLORREF CmyView::JointColor;
COLORREF CmyView::MemberColor;
COLORREF CmyView::SlabColor;
COLORREF CmyView::HighlightColor;
COLORREF CmyView::SelectionColor;
COLORREF CmyView::RestraintColor;
COLORREF CmyView::ConstraintColor;
COLORREF CmyView::AxisColor;
COLORREF CmyView::ReleaseColor;
COLORREF CmyView::TextColor;
COLORREF CmyView::GridColor;
COLORREF CmyView::PositiveColor;
COLORREF CmyView::NegativeColor;
UINT CmyView::NCONTORS = 15;
COLORREF CmyView::ContorColors[15];
/*
View construction
*/
CmyView::CmyView() { 
	view = CREATION[nViews]; 
	nViews++;
	pAnimationDia = NULL;
}
CmyView::~CmyView() { 
	nViews--; 
}
BOOL CmyView::PreCreateWindow(CREATESTRUCT &cs) {
	cs.lpszClass = AfxRegisterWndClass(NULL,(HCURSOR)::LoadCursor(NULL,IDC_ARROW));
	return CView::PreCreateWindow(cs);
}
int CmyView::OnCreate(LPCREATESTRUCT lpCreateStruct) {
	if(CView::OnCreate(lpCreateStruct) == -1)
		return -1;
    CmyDocument* pDoc = GetDocument();
	
	hasGridPoint = FALSE;
	hadGridPoint = FALSE;
    has_tracker = FALSE;
	had_tracker = FALSE;
    has_plane = FALSE;
	had_plane = FALSE;
	force_diagram = INONE;
	average_stress = TRUE;
	DrawingHasChanged = TRUE;
	perspective_toggle = view;
	if(view != VIEW_3D) {
		GRIDLIST* pgrid = &pDoc->global->grid[NORMGRID[view]];
		if(!pgrid->IsEmpty()) {
			if(view == VIEW_XY) position = pgrid->GetCount() - 1;
			else position = 0;
			pvalue = pgrid->GetAt(pgrid->FindIndex(position));
		} 
	}
	/*create tooltip*/
	m_toolTip.Create(this);
	m_ti.cbSize = sizeof(m_ti);
	m_ti.uFlags = TTF_IDISHWND|TTF_TRACK|TTF_ABSOLUTE|TTF_TRANSPARENT;
	m_ti.hwnd = m_hWnd;
	m_ti.hinst = NULL;
	m_ti.uId = (UINT)m_hWnd;
    m_ti.lpszText = NULL;
	m_ti.rect.left = 0;
    m_ti.rect.top = 0;
	m_ti.rect.right = 0;
	m_ti.rect.bottom = 0;
	m_toolTip.SendMessage(TTM_ADDTOOL, 0, (LPARAM)&m_ti);

	/*memory dc*/
	CDC* pDC = GetDC();
	memX = GetSystemMetrics(SM_CXSCREEN);
	memY = GetSystemMetrics(SM_CXSCREEN);
	memDC.CreateCompatibleDC(pDC);
	memBmp.CreateCompatibleBitmap(pDC,memX,memY);
	memDC.SelectObject(&memBmp);
	ReleaseDC(pDC);

	return 0;
}
/*
Default colors
*/
static COLORREF default_contor_colors[] = {
	RGB(255,0,  255),
    RGB(255,0,  128),
    RGB(255,0,    0),
	RGB(255, 64,  0),
	RGB(255,128,  0),
	RGB(255,192,  0),
	RGB(255,224,  0),
	RGB(255,255,  0),
    RGB(128,255,  0),
	RGB(  0,255,  0),
    RGB(  0,255,128),
	RGB(  0,255,255),
	RGB(  0,192,255),
	RGB(  0,128,255),
	RGB(  0,  0,255) 
};
void CmyView::ResetColors() {
	BackColor = RGB(0,0,0);
	JointColor = RGB(255,255,0);
	MemberColor = RGB(255,255,0);
	SlabColor = RGB(255,0,0);
	HighlightColor = RGB(255,0,0);
	SelectionColor = RGB(255,255,255);
	RestraintColor = RGB(0,255,0);
	ConstraintColor = RGB(0,255,255);
	AxisColor = RGB(100,255,255);
	ReleaseColor = RGB(0,255,0);
	TextColor = RGB(0,255,0);
	GridColor = RGB(100,100,100);
	PositiveColor = RGB(255,255,0);
	NegativeColor = RGB(255,0,0);
	for(UINT i = 0;i < NCONTORS;i++) {
		ContorColors[i] = default_contor_colors[i];
	}
}
/*
Animation Dialog
*/
UINT  AnimateProc(LPVOID pViewt) {
	CmyView* pView = (CmyView*) pViewt;
	pView->Animate();
    return 0;
}

BEGIN_MESSAGE_MAP(CAnimationDia , CDialog)
    ON_BN_CLICKED(IDC_SOUND,OnSound)
	ON_COMMAND_RANGE(IDC_FASTER,IDC_SLOWER,OnSpeed)
	ON_COMMAND(IDC_REDO,OnRedo)
END_MESSAGE_MAP()

CAnimationDia::CAnimationDia(CWnd* parent,CmyView* pViewt) :
           CDialog(IDD_ANIMATION,parent) {
	interval = 400;
	sound = FALSE;
	stop = FALSE;
	finished = FALSE;
	pView = pViewt;
}
void CAnimationDia::DoDataExchange(CDataExchange* pDX) {
    DDX_Control(pDX,IDC_TIME,c_Time);
    DDX_Control(pDX,IDC_REDO,c_Redo);
}
BOOL CAnimationDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
    c_Redo.EnableWindow(finished);
	return TRUE;
}
void CAnimationDia::OnOK() {
	stop = TRUE;
	while(!finished) Sleep(10);
	CDialog::OnOK();
}
void CAnimationDia::OnSound() {
	sound = !sound;
}
void CAnimationDia::OnSpeed(UINT nID) {
	if(nID == IDC_FASTER) {
		if(interval > 50) interval -= 50;
	} else interval += 50;
}
void CAnimationDia::OnRedo() {
	finished = FALSE;
	c_Redo.EnableWindow(finished);
	AfxBeginThread(AnimateProc,pView,THREAD_PRIORITY_NORMAL);
}
void CAnimationDia::Finish() {
	finished = TRUE;
	c_Redo.EnableWindow(finished);
}
/*
Animation
*/
void CmyView::Animate() {
	
    CmyDocument* pDoc = GetDocument();
	CFrameWnd* parent = GetParentFrame();
	CClientDC clientDC(this);
	CString str;

	pDoc->DetermineHistoryScale();

	/*animation*/
	CBrush brush(BackColor);
	CBrush* pOld = memDC.SelectObject(&brush);

    RESPONSEHIST* pcase = (RESPONSEHIST*) pDoc->c_AnalysisCase;
	for(UINT i = 0;i < pcase->N;i++) {
		if(pAnimationDia->stop)
			break;

		pDoc->AssignIndex(pcase->index + i);

		str.Format("Time %.2f",i * pcase->dt);
		pAnimationDia->c_Time.SetWindowText(str);
		
		/*Draw on memory DC and then BitBlt to avoid flicker*/
		memDC.PatBlt(0,0,memX,memY,PATCOPY);
		DrawView(&memDC);

		/*BitBlt on screen*/
		clientDC.BitBlt(0,0,memX,memY,&memDC,0,0,SRCCOPY);

		Sleep(pAnimationDia->interval);
		if(pAnimationDia->sound) 
			MessageBeep(0xFFFFFFFF);
	}
	
	memDC.SelectObject(pOld);
	/*end*/

    pAnimationDia->Finish();
}
void CmyView::OnAnimate() {
	CmyDocument* pDoc = GetDocument();

	CAnimationDia Dia(AfxGetMainWnd(),this);
	pAnimationDia = &Dia;
	AfxBeginThread(AnimateProc,this,THREAD_PRIORITY_NORMAL);
	Dia.DoModal();
	pAnimationDia = NULL;
	pDoc->AssignIndex(pDoc->c_AnalysisCase->index);
	InitView();
}
void CmyView::OnUpdateAnimateUI(CCmdUI* pCmdUI) {
	CmyDocument* pDoc = GetDocument();
	pCmdUI->Enable(pDoc->AnalysisResult && pDoc->c_AnalysisCase->acase_type == ANALYSISCASE::RESPONSEH_CASE
		&& force_diagram != INONE);
}
void CmyView::CalculateReye() {
	CmyDocument* pDoc = GetDocument();
	RPoint rp;

	if(pDoc->global->coordinate == CARTESIAN) rp = eye - gcenter;
	else rp = eye - RDtoCT(gcenter);	
	reye.x = sqrt(rp.x * rp.x + rp.y * rp.y + rp.z * rp.z);
	reye.y = acos(rp.z / reye.x);
	reye.z = myatan(rp.y , rp.x);
}
void CmyView::AdjustCenter() {
	CmyDocument* pDoc = GetDocument();
	GRIDINFO ginfo;
	pDoc->GetGridInfo(&ginfo);
	
	/*center of grid*/
	if(pDoc->global->coordinate == RADIAL && view != VIEW_XZ) {
		gcenter.x = (1 - (ginfo.ymin + ginfo.ymax) / 360) * ginfo.xmax;
	} else {
		gcenter.x = (ginfo.xmin + ginfo.xmax) / 2;
	}
	gcenter.y = (ginfo.ymin + ginfo.ymax) / 2;
	gcenter.z = (ginfo.zmin + ginfo.zmax) / 2;
	if(view == VIEW_XY) gcenter.z = pvalue; 
	else if(view == VIEW_XZ) gcenter.y = pvalue;
	else if(view == VIEW_YZ) gcenter.x = pvalue;
	
	/*Set Eye*/
	DOUBLE emax = 0,emin = 0;
	if(pDoc->global->coordinate == RADIAL) {
		if(6 * ginfo.xmax > emax) 
			emax = 6 * ginfo.xmax;
		if(!emin && emax) emin = emax;
	} else {
		if(ginfo.xmax - ginfo.xmin > emax)
			emax = ginfo.xmax - ginfo.xmin;
		if(!emin && emax) emin = emax;

		if(ginfo.ymax - ginfo.ymin > emax)
			emax = ginfo.ymax - ginfo.ymin;
		if(!emin && emax) emin = emax;
	}
	if(ginfo.zmax - ginfo.zmin > emax)
	   emax = ginfo.zmax - ginfo.zmin;
	if(!emin && emax) emin = emax;

	eye.x = 1.6 * emax;
	eye.y = 1.6 * emax;
	eye.z = -2 * emax;

	DOUBLE r = 0.6 * pow(emax / emin, 1/3.0);
	eye.x *= r;
	eye.y *= r;
	eye.z *= r;
	
	if(pDoc->global->coordinate == CARTESIAN)
		eye = eye + gcenter;

    CalculateReye();
}
void CmyView::AdjustScale() {
	CmyDocument* pDoc = GetDocument();

	CRect ClRect;
	GetClientRect(&ClRect); 
	Corner = ClRect.BottomRight();
	Center.x = ClRect.Width() / 2;
	Center.y = ClRect.Height() / 2;
	Scale = 1;

	RPoint myp;
	CPoint mycp,maxcp(-10000,-10000),mincp(10000,10000);
	SYSTEM* system;
	GRIDLIST* grids;
	POSITION pos,posx,posy,posz;

	pos = pDoc->systems.GetHeadPosition();
	while(pos) {
		system = &pDoc->systems.GetNext(pos);
		grids = system->grid;
		
		posx = grids[IUX].GetHeadPosition(); 
		while(posx) {
			myp.x = grids[IUX].GetNext(posx);
			if(view == VIEW_YZ && !EQUAL(myp.x , pvalue)) continue;
			
			posy = grids[IUY].GetHeadPosition(); 
			while(posy) {
				myp.y = grids[IUY].GetNext(posy);
				if(view == VIEW_XZ && !EQUAL(myp.y , pvalue)) continue;
				
				posz = grids[IUZ].GetHeadPosition(); 
				while(posz) {
					myp.z = grids[IUZ].GetNext(posz);
					if(view == VIEW_XY && !EQUAL(myp.z , pvalue)) continue;
					
					mycp = RPtoCP(myp,system,(system->coordinate == RADIAL));

					if(mycp.x > maxcp.x) maxcp.x = mycp.x;
					if(mycp.y > maxcp.y) maxcp.y = mycp.y;
					if(mycp.x < mincp.x) mincp.x = mycp.x;
					if(mycp.y < mincp.y) mincp.y = mycp.y;
				}
			}
		}
	}
	Scale = min(fabs((Corner.x) / DOUBLE(maxcp.x - mincp.x)),
		        fabs((Corner.y) / DOUBLE(maxcp.y - mincp.y)));
	
	Scale *= 0.8;
	if(view == VIEW_3D) Scale *= 0.9;

	SELSIZE = (1 / Scale) * 20;

	Center.x = Center.x - int(Scale * ((maxcp.x + mincp.x) / 2 - Center.x));
    Center.y = Center.y - int(Scale * ((maxcp.y + mincp.y) / 2 - Center.y));
}
void CmyView::AdjustPosition() {
	AdjustCenter();
	AdjustScale();
	DrawingHasChanged = TRUE;
}
void CmyView::OnSize(UINT nType,int cx,int cy) {
	CView::OnSize(nType,cx,cy);
	AdjustPosition();
}
void CmyView::OnRestoreView() {
	AdjustPosition();
	Invalidate();
}
void CmyView::OnZoomIn() {
	Scale *= 1.2;
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::OnZoomOut() {
	Scale /= 1.2;
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::OnRubberZoom() {
	CmyDocument* pDoc = GetDocument();
	pDoc->nClicks = 0;
	pDoc->state = RUBBER_ZOOM;
}
void CmyView::OnPan() {
	CmyDocument* pDoc = GetDocument();
	pDoc->nClicks = 0;
	pDoc->state = PAN;
}
BOOL CmyView::OnMouseWheel(UINT nFlags,short zDelta,CPoint pt) {
	DOUBLE factor = 1.05f;
	POINT myp = pt;
	ScreenToClient(&myp);
	CRect rect(0,0,Corner.x,Corner.y);
	if(!rect.PtInRect(myp))
		return FALSE;
	pt = myp;
	if(zDelta > 0) {
        Center.x = pt.x + int((Center.x - pt.x) * factor);
		Center.y = pt.y + int((Center.y - pt.y) * factor);
		Scale *= factor;
	} else {
        Center.x = pt.x + int((Center.x - pt.x) / factor);
		Center.y = pt.y + int((Center.y - pt.y) / factor);
		Scale /= factor;
	}
	DrawingHasChanged = TRUE;
	Invalidate();
	return TRUE;
}
void CmyView::OnDisplayOptions() {
	CmyDocument* pDoc = GetDocument();
	CDisplayOptionsDia DisplayDia(AfxGetMainWnd(),&mydisplay);
    if(DisplayDia.DoModal() == IDOK) {
		if(DisplayDia.ApplyToAll) {
			CmyView* pView;
			POSITION pos = pDoc->GetFirstViewPosition();
			while (pos != NULL) {
				pView = (CmyView*) pDoc->GetNextView(pos);
				pView->mydisplay = DisplayDia.mydisplay;
				pView->DrawingHasChanged = TRUE;
				pView->Invalidate();
			} 
		} else {
			mydisplay = DisplayDia.mydisplay;
			DrawingHasChanged = TRUE;
			Invalidate();
		}
	}
}
void CmyView::OnShowGrids() {
	mydisplay.Grid = !mydisplay.Grid;
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::OnShowOrigin() {
	mydisplay.Origin = !mydisplay.Origin;
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::OnRefresh() {
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::InitView() {
	SetTitle();
	AdjustPosition();
	Invalidate();
}
void CmyView::SetTitle() {
	CmyDocument* pDoc = GetDocument();
	CFrameWnd* parent = GetParentFrame();

	CString str,str1;

	if(pDoc->AnalysisResult && force_diagram != INONE) {
		if(force_diagram <= IRZ) {
			str = "Internal Force Diagram - " + ForceString[force_diagram];
		} else if(force_diagram >= ISX && force_diagram <= IVMAX) {
			if(force_diagram <= ISVMAX) str = "Slab Stress - " + StressString[force_diagram - ISX];
			else str = "Slab Force - " + StressString[force_diagram - ISX];
		} else {
			if(force_diagram == IDEFLECTED) str = "Deformed Shape";
			else if(force_diagram == IREACTIONS) str = "Joint Reactions";
            else if(force_diagram == IDISPLACEMENT) str = "Joint Displacements";

			if(pDoc->c_AnalysisCase->acase_type == ANALYSISCASE::MODAL_CASE) {
				MODALCASE* pcase = (MODALCASE*) pDoc->c_AnalysisCase;
                str1.Format(" - Mode %d - Period %.6f",pcase->selmode,(2 * PI) / sqrt(pcase->eigvalue[pcase->selmode]));
				str += str1;
			} else if(pDoc->c_AnalysisCase->acase_type == ANALYSISCASE::BUCKLING_CASE) {
				BUCKLINGCASE* pcase = (BUCKLINGCASE*) pDoc->c_AnalysisCase;
                str1.Format(" - Mode %d - Factor %.6f",pcase->selmode,pcase->eigvalue[pcase->selmode]);
				str += str1;
			}
		} 
	} else {
		if(view == VIEW_3D) {
			str = "3D-View";
		} else if(view == VIEW_XY) {
			str.Format("XY-View @ %.2f",pvalue);
		} else if(view == VIEW_XZ) {
			str.Format("XZ-View @ %.2f",pvalue);
		} else if(view == VIEW_YZ) {
			str.Format("YZ-View @ %.2f",pvalue);
		}
	}
	parent->SetWindowText(str);
	DrawingHasChanged = TRUE;
}
BOOL CmyView::InPlane(RPoint rp) {
	CmyDocument* pDoc = GetDocument();
	RPoint mrp = rp;
	if(pDoc->global->coordinate == RADIAL) {
		rp = CTtoRD(rp);
		if(view == VIEW_XZ) {
			if(EQUAL(rp.x,0) || 
			   EQUAL(rp.y,pvalue + 360) || 
			   EQUAL(pvalue,rp.y + 360)) 
				return TRUE;
		}
	}
	if(view == VIEW_XY) {
		if(!EQUAL(rp.z,pvalue)) return FALSE;
	} else if(view == VIEW_XZ) {
		if(!EQUAL(rp.y,pvalue)) return FALSE;
	} else if(view == VIEW_YZ) {
		if(!EQUAL(rp.x,pvalue)) return FALSE;
	}
	return TRUE;
}
RPoint CmyView::CPtoRP(CPoint cp) {
	CmyDocument* pDoc = GetDocument();
	RPoint myp;
	DOUBLE x = (cp.x - Center.x) / Scale;
	DOUBLE y = (Center.y - cp.y) / Scale;

	if(view == VIEW_3D) {
	
		RPoint v;
		DOUBLE d,q;
		myp.set(x,y,0);
		d = sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
		myp = aperspective(myp,-d);
		v.set(1,0,0);
		d = sqrt(eye.x * eye.x + eye.y * eye.y);
		q = (myatan(eye.z , d) + PI / 2);
		myp = rotate(myp,v,q);
		v.set(0,0,1);
		q = -myatan(eye.x , eye.y);
		myp = rotate(myp,v,q);

	} else if(pDoc->global->coordinate == CARTESIAN) {
		if(view == VIEW_XY) {
			myp.x = x; 
			myp.y = y; 
			myp.z = pvalue;
		} else if(view == VIEW_XZ) {
			myp.x = x; 
			myp.y = pvalue; 
			myp.z = y;
		} else if(view == VIEW_YZ) {
			myp.x = pvalue; 
			myp.y = x; 
			myp.z = y;
		}
	} else {
		if(view == VIEW_XY) {
			myp.x = sqrt(x * x + y * y); 
			myp.y = 180 * myatan(y , x) / PI;
			myp.z = pvalue;
		} else if(view == VIEW_XZ) {
			myp.x = x; 
			myp.y = pvalue; 
			myp.z = y;
		} else if(view == VIEW_YZ) {
			myp.x = pvalue; 
			myp.y = 180 * (x / myp.x) / PI; 
			myp.z = y;
		}
		myp = RDtoCT(myp);
	}
	
	if(pDoc->global->coordinate == CARTESIAN) {
		RPoint gc = gcenter;
		if(view == VIEW_XY) gc.z = 0;
		else if(view == VIEW_XZ) gc.y = 0;
		else if(view == VIEW_YZ) gc.x = 0;
		myp = myp + gc;
	}
	return myp;
}
RPoint CmyView::RPtoRP(RPoint rp,SYSTEM* system, BOOL radial) {
	CmyDocument* pDoc = GetDocument();
	RPoint myp;

	MATRIX R;
	mat_alloc(R,3,3);
	GetPointRotMatrix(system->rotation,R);

	/*cartesian*/
	if(radial)
		rp = RDtoCT(rp);
	rp = GetGlobal(rp,R);
	rp = rp + system->origin;
	if(pDoc->global->coordinate == CARTESIAN) 
		rp = rp - gcenter;
	
	if(view == VIEW_3D || pDoc->global->coordinate == CARTESIAN) {
		if(view == VIEW_3D) {
			RPoint v;
			DOUBLE d,q;		
			v.set(0,0,1);
			q = +(PI / 2 - myatan(eye.y , eye.x));
			myp = rotate(rp,v,q);
			v.set(1,0,0);
			d = sqrt(eye.x * eye.x + eye.y * eye.y);
			q = -(PI / 2 + myatan(eye.z , d));
			myp = rotate(myp,v,q);
			d = sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
			myp = perspective(myp,-d );
		} else if(view == VIEW_XY) {
			myp.x = rp.x; 
			myp.y = rp.y;
		} else if(view == VIEW_XZ) {
			myp.x = rp.x; 
			myp.y = rp.z;
		} else if(view == VIEW_YZ) {
			myp.x = rp.y; 
			myp.y = rp.z;
		} 
	} else {
		rp = CTtoRD(rp);
		
		if(view == VIEW_XY) {
			myp.x = rp.x * cos(PI * rp.y / 180); 
			myp.y = rp.x * sin(PI * rp.y / 180);
		} else if(view == VIEW_XZ) {
			myp.x = rp.x; 
			myp.y = rp.z;
		} else if(view == VIEW_YZ) {
			myp.x = (PI * rp.y / 180) * rp.x; 
			myp.y = rp.z;
		}
	}
	vec_free(R);
	return myp;
}

CPoint CmyView::RPtoCP(RPoint rp,SYSTEM* system,BOOL radial) {
	CmyDocument* pDoc = GetDocument();
	CPoint p;
	RPoint myp;
	if(!system) myp = RPtoRP(rp,pDoc->global,radial);
	else  myp = RPtoRP(rp,system,radial);
	p.x = Center.x + int(Scale * myp.x);
    p.y = Center.y - int(Scale * myp.y); 
	return p;
}

void CmyView::DrawGridPoint(RPoint& rp, BOOL found) {
	if(!mydisplay.Grid || (GetDocument()->AnalysisResult && force_diagram == IDEFLECTED))
		return;

	CClientDC clientDC(this);
	int mode = clientDC.SetROP2(R2_XORPEN );
	
	hasGridPoint = found;

	if(hadGridPoint) 
		DrawCircle(&clientDC,oldGridPoint,4,HighlightColor);
	if(hasGridPoint) {
		DrawCircle(&clientDC,rp,4,HighlightColor);
		oldGridPoint = rp;
	}
	
	hadGridPoint = found;

	clientDC.SetROP2(mode);
}
RPoint CmyView::SnapToGrid(CPoint& cp,int& found,JOINT** psel_joint) {
	RPoint rp = CPtoRP(cp),rp1;
	CmyDocument* pDoc = GetDocument();
	
    found = 0;
	
	if(pDoc->Snap_To_Joint) {
		JOINT* joint;
		if(joint = FindJoint(rp)) {
			found = SNAP_JOINT;
			rp = joint->p;
			if(psel_joint)
				*psel_joint = joint;
			return rp;
		}
	}
	
	if(pDoc->Snap_To_MemberMid || pDoc->Snap_To_Member) {
		MEMBER* member;
		RPoint pr;
		if(member = FindMember(rp,&pr)) {
			if(pDoc->Snap_To_MemberMid && !member->is_curved) {
				RPoint mid = (member->j1->p + member->j2->p) * 0.5;
				CPoint p1,p2;
				p1 = RPtoCP(rp);
				p2 = RPtoCP(mid);
				DOUBLE d = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
				if(d <= 5) {
					found = SNAP_MEMBERMID;
					return mid;
				}
			} else if(pDoc->Snap_To_Member) {
				found = SNAP_MEMBER;
				return pr;
			}
		}
	}
	
	if(pDoc->Snap_To_Grid) {
		SYSTEM* system;
		POSITION pos = pDoc->systems.GetHeadPosition();
		while(pos) {
			system = &pDoc->systems.GetNext(pos);
			rp1 = SnapToGridLocal(system,cp,found);
			if(found) {
				found = SNAP_GRID;
				return rp1;
			}
		}
	}

	return rp;
}
RPoint CmyView::SnapToGridLocal(SYSTEM* system,CPoint& cp,BOOL& found) {
	RPoint myp;
	CPoint mycp;
	GRIDLIST* grids = system->grid;
	POSITION posx,posy,posz;
	DOUBLE d;

	posx = grids[IUX].GetHeadPosition(); 
	while(posx) {
		myp.x = grids[IUX].GetNext(posx);
		if(view == VIEW_YZ && !EQUAL(myp.x , pvalue)) continue;
		
		posy = grids[IUY].GetHeadPosition(); 
		while(posy) {
			myp.y = grids[IUY].GetNext(posy);
			if(view == VIEW_XZ && !EQUAL(myp.y , pvalue)) continue;
         
			posz = grids[IUZ].GetHeadPosition(); 
			while(posz) {
				myp.z = grids[IUZ].GetNext(posz);
				if(view == VIEW_XY && !EQUAL(myp.z , pvalue)) continue;

				mycp = RPtoCP(myp,system,(system->coordinate == RADIAL));

				d = sqrt(pow(cp.x - mycp.x,2) + pow(cp.y - mycp.y,2));
				if(d <= 10) {
					MATRIX R;
					mat_alloc(R,3,3);

					GetPointRotMatrix(system->rotation,R);
					if(system->coordinate == RADIAL) 
						myp = RDtoCT(myp);
					myp = GetGlobal(myp,R);
					myp = myp + system->origin;

					vec_free(R);
					found = TRUE;
					return myp;
				}
			}
		}
	}
	
	found = FALSE;
	return myp;
}
static BOOL Intersect(CPoint& p1,CPoint& p2,CPoint& p3,CPoint& p4) {
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
			return FALSE;
		x = ((p3.y - m2 * p3.x) - (p1.y - m1 * p1.x)) / (m1 - m2);
		y = m1 * x + (p1.y - m1 * p1.x);
	}

    if( ((x >= p1.x && x <= p2.x) || (x >= p2.x && x <= p1.x)) &&
        ((y >= p1.y && y <= p2.y) || (y >= p2.y && y <= p1.y)) &&
	    ((x >= p3.x && x <= p4.x) || (x >= p4.x && x <= p3.x)) &&
        ((y >= p3.y && y <= p4.y) || (y >= p4.y && y <= p3.y)) )  
	    return TRUE;

	return FALSE;
}

void CmyView::DrawTracker(CPoint& point,int type) {
	CmyDocument* pDoc = GetDocument();
	CClientDC clientDC(this);
	COLORREF col = clientDC.SetBkColor(BackColor);
	int mode = clientDC.SetROP2(R2_XORPEN);
	CPen pen(PS_DOT,1,SelectionColor);
	CPen* old = clientDC.SelectObject(&pen);

	if(type == LINE_TRACK) {
		CPoint sel_point = RPtoCP(pDoc->selrp[0]);
		clientDC.MoveTo(sel_point);
		clientDC.LineTo(point);
		clientDC.MoveTo(sel_point);
		clientDC.LineTo(prev_mouse_cp);
	} else if(type == RECT_TRACK) {
		CPoint sel_point = RPtoCP(pDoc->selrp[0]);
		clientDC.MoveTo(sel_point);
		clientDC.LineTo(sel_point.x,point.y);
		clientDC.LineTo(point);
		clientDC.LineTo(point.x,sel_point.y);
		clientDC.LineTo(sel_point);
		
		clientDC.MoveTo(sel_point);
		clientDC.LineTo(sel_point.x,prev_mouse_cp.y);
		clientDC.LineTo(prev_mouse_cp);
		clientDC.LineTo(prev_mouse_cp.x,sel_point.y);
		clientDC.LineTo(sel_point);
	} else {
		CPoint p[NPOLY],wpoint;
		int found;
		UINT i,j;
		point = RPtoCP(SnapToGrid(point,found));
		CPoint sel_point = RPtoCP(pDoc->selrp[0]);
		
		for(i = 0;i < 2;i++) {
			if(i == 0) {
				if(had_tracker) wpoint = prev_mouse_cp;
				else continue;
			} else if(i == 1) {
				if(has_tracker) wpoint = point;
				else continue;
			}
			
			if(pDoc->state == ADD_QUAD_AREA) {

				for(j = 0; j < pDoc->nClicks;j++)
					p[j] = RPtoCP(pDoc->selrp[j]);
				p[j] = wpoint;

				clientDC.MoveTo(p[0]);
				clientDC.LineTo(p[pDoc->nClicks]);
				if(pDoc->nClicks > 1)
					clientDC.LineTo(p[pDoc->nClicks - 1]);
				if(pDoc->nClicks == 1) 
					clientDC.LineTo(p[0]);
			} else if(pDoc->state == ADD_RECT_AREA) {
				clientDC.MoveTo(sel_point);
				clientDC.LineTo(sel_point.x,wpoint.y);
				clientDC.LineTo(wpoint);
				clientDC.LineTo(wpoint.x,sel_point.y);
				clientDC.LineTo(sel_point);
			} else {
				clientDC.MoveTo(sel_point);
				clientDC.LineTo(wpoint);
			}
		}
		prev_mouse_cp = point;
	}

	clientDC.SetBkColor(col);
	clientDC.SetROP2(mode);
	clientDC.SelectObject(old);
}
void CmyDocument::OnEnter() {
	if(state == ADD_QUAD_AREA && nClicks) {
		if(nClicks >= NPOLY) {
			nClicks = 0;
			return;
		}
		
		SLAB* pslab;
		JOINT J[NPOLY] ;
		
		for(UINT i = 0;i < nClicks;i++)
			J[i].p = selrp[i];
		
		if(pslab = AddSlab(J,nClicks)) {
			pslab->nDivx = 10;
			pslab->nDivy = 10;
			SetModifiedFlag();
			UndoRedo.Save(pslab);
		}
		UpdateDrawing();
		nClicks = 0;
	}
}
void CmyView::OnLButtonDown(UINT nFlags,CPoint point) {
	
	CView::OnLButtonDown(nFlags,point);
    CmyDocument* pDoc = GetDocument();

	int found;
	RPoint rpoint = SnapToGrid(point,found);
	if(!found && (view == VIEW_3D) && (pDoc->state >= ADD_JOINT  && pDoc->state <= ADD_QUAD_AREA))
		return;

	if(pDoc->state == PAN) {
		pDoc->nClicks++;
		prev_mouse_cp = point;
		SetCapture();
	} else if(pDoc->state == RUBBER_ZOOM) {
		pDoc->nClicks++;
		prev_mouse_cp = point;
		pDoc->selrp[0] = rpoint;
		SetCapture();
	} else if(pDoc->state == SELECT) {
        if(SelectJoint(rpoint) ||
        SelectMember(rpoint) ||
		SelectSlab(rpoint)
		) {
			pDoc->UpdateDrawing();
			return;
		}
		pDoc->nClicks++;
		prev_mouse_cp = point;
		pDoc->selrp[0] = rpoint;
		SetCapture();
    } else if(pDoc->state == LINE_SELECT) {
		pDoc->nClicks++;
		prev_mouse_cp = point;
		pDoc->selrp[0] = rpoint;
		SetCapture();
	} else if(pDoc->state == ADD_JOINT) {
        JOINT joint;
		JOINT* pjoint;
		joint.p = rpoint;
	
		pjoint = pDoc->AddJoint(&joint);
		pDoc->SetModifiedFlag();
		pDoc->UpdateDrawing();
		pDoc->UndoRedo.Save(pjoint);
	} else if(pDoc->state == ADD_MEMBER ||
		pDoc->state == ADD_CURVED_MEMBER) {

		pDoc->nClicks++;
		if(pDoc->nClicks == 2) {

			had_tracker = TRUE;
			has_tracker = FALSE;
			DrawTracker(point,MEMBER_TRACK);

			if(pDoc->state == ADD_CURVED_MEMBER) {
				CCurvedMemberDia Dia(this);
				if(Dia.DoModal() == IDOK) {
					MEMBER member;
					member.is_curved = Dia.single;
					int N = Dia.nDiv,type = Dia.type;
					RPoint rp1,rp2,rp3,e;
					rp1 = pDoc->selrp[0];
					rp2 = Dia.rp;
					rp3 = rpoint;

					DOUBLE L = distance(rp3,rp1),Ye,Xe,X,Y,R,Q,q;
					e = ptprojection(rp2,rp1,rp3);
					Ye = distance(e,rp2);
					Xe = distance(rp1,e);
					e = rp2 - e;
					e = e.unitvector();

					JOINT j1,j2;
					j1.p = rp1;
					if(member.is_curved) {
						member.curved_list.AddTail(j1.p - rp1);
					}
					if(type == CCurvedMemberDia::CIRCULAR) {
						DOUBLE H = sqrt(Xe * Xe + Ye * Ye);
						q = atan(Ye / (L - Xe));
						R = H / (2 * sin(q));
						Q = q;
						Q += atan(Ye / Xe);
						Q *= 2;

						DOUBLE d1,d2;
                        for(int i = 0;i < N;i++) {
							q = Q * (float(i+1) / N);
                            d1 = 2 * R * sin (q / 2);
							d2 = 2 * R * sin ((Q - q) / 2);
							X = (L * L + d1 * d1 - d2 * d2) / (2 * L);
							Y = sqrt(d1 * d1 - X * X);
							j2.p = rp1 + (rp3 - rp1) * (X / L);
							j2.p = j2.p + e * Y;
							if(member.is_curved) {
								member.curved_list.AddTail(j2.p - rp1);
							} else {
								pDoc->AddMember(&j1,&j2);
							}
							j1 = j2;
						}
					} else {
						for(int i = 0;i < N;i++) {
							j2.p = rp1 + (rp3 - rp1) * (float(i+1) / N);
							X = distance(j2.p,rp1);
							if(type == CCurvedMemberDia::STRAIGHT) Y = 0;
							else if(type == CCurvedMemberDia::PARABOLIC) Y = X * (X - L) * Ye / (Xe * (Xe - L));
							j2.p = j2.p + e * Y;
							if(member.is_curved) {
								member.curved_list.AddTail(j2.p - rp1);
							} else {
								pDoc->AddMember(&j1,&j2);
							}
							j1 = j2;
						}
					}
					if(member.is_curved) {
						j1.p = member.curved_list.GetHead() + rp1;
                        j2.p = member.curved_list.GetTail() + rp1;
						MEMBER* pmem = pDoc->AddMember(&j1,&j2);
						pmem->CopyCurveData(&member);
					}

					pDoc->SetModifiedFlag();
					pDoc->UpdateDrawing();
				}
			} else {
				MEMBER* pmem;
				JOINT joint,joint1;
				joint.p = pDoc->selrp[0];
				joint1.p = rpoint;
				
				if(pmem = pDoc->AddMember(&joint,&joint1)) {
					pDoc->SetModifiedFlag();
					pDoc->UpdateDrawing();
					pDoc->UndoRedo.Save(pmem);
				}
			}
			pDoc->nClicks = 1;
		}
		pDoc->selrp[0] = rpoint;
	} else if(pDoc->state == ADD_RECT_AREA || pDoc->state == ADD_QUAD_AREA) {
		if(view == VIEW_3D)
			return;
		pDoc->nClicks++;

		if(pDoc->state == ADD_QUAD_AREA) {
			pDoc->selrp[pDoc->nClicks - 1] = rpoint;
			if(OnSameLine(pDoc->selrp[pDoc->nClicks - 3],
				          pDoc->selrp[pDoc->nClicks - 2],
						  pDoc->selrp[pDoc->nClicks - 1])) {
				pDoc->nClicks--;
				return;
			}

			/*Draw temporary object*/
			UINT i;
		    CPoint p[NPOLY];
			CClientDC clientDC(this);
	        CPen pen(PS_SOLID,2,TextColor),*pOld;
            pOld = clientDC.SelectObject(&pen);
			for(i = 0;i < pDoc->nClicks;i++)
				p[i] = RPtoCP(pDoc->selrp[i]);
			clientDC.MoveTo(p[0]);
			for(i = 0;i < pDoc->nClicks;i++)
				clientDC.LineTo(p[i]);
			clientDC.LineTo(p[0]);
			
            CPen pen1(PS_SOLID,3,BackColor);
			clientDC.SelectObject(&pen1);
			clientDC.MoveTo(p[0]);
			if(i > 3) 
				clientDC.LineTo(p[i - 2]);

			clientDC.SelectObject(pOld);
		} else {
			if(pDoc->nClicks == 2) {

				pDoc->nClicks = 0;
				
				SLAB* pslab;
				JOINT j1,j2,j3,j4;
				j1.p = pDoc->selrp[0];
				j3.p = rpoint;
				if(view == VIEW_XY) {
					j2.p.x = j3.p.x;
					j2.p.y = j1.p.y; 
					j4.p.x = j1.p.x;
					j4.p.y = j3.p.y; 
					j2.p.z = j1.p.z;
					j4.p.z = j1.p.z;
				} else if(view == VIEW_XZ) {
					j2.p.x = j3.p.x;
					j2.p.z = j1.p.z; 
					j4.p.x = j1.p.x;
					j4.p.z = j3.p.z;
					j2.p.y = j1.p.y;
					j4.p.y = j1.p.y;
				} else if(view == VIEW_YZ) {
					j2.p.y = j3.p.y;
					j2.p.z = j1.p.z; 
					j4.p.y = j1.p.y;
					j4.p.z = j3.p.z;
					j2.p.x = j1.p.x;
					j4.p.x = j1.p.x;
				}
				
				if(pslab = pDoc->AddSlab(&j1,&j2,&j3,&j4)) {
					pDoc->SetModifiedFlag();
					pDoc->UpdateDrawing();
					pDoc->UndoRedo.Save(pslab);
				}
			} else {
				prev_mouse_cp = point;
				pDoc->selrp[0] = rpoint;
			}
		}
	} else if(pDoc->state == QUICK_DRAW_MEMBER) {
		DOUBLE d,dmin1 = 1e8,dmin2 = 1e8;
		POSITION posx,posy,posz;
		RPoint myp,p1,p2;
		GRIDLIST* grids = pDoc->global->grid;
		CPoint mycp;
		UINT pcount = 0;
		
		/*select two nearest grid points*/
		posx = grids[IUX].GetHeadPosition(); 
		while(posx) {
			myp.x = grids[IUX].GetNext(posx);
			if(view == VIEW_YZ && !EQUAL(myp.x , pvalue)) continue;
			
			posy = grids[IUY].GetHeadPosition(); 
			while(posy) {
				myp.y = grids[IUY].GetNext(posy);
				if(view == VIEW_XZ && !EQUAL(myp.y , pvalue)) continue;
				
				posz = grids[IUZ].GetHeadPosition(); 
				while(posz) {
					myp.z = grids[IUZ].GetNext(posz);
					if(view == VIEW_XY && !EQUAL(myp.z , pvalue)) continue;
					
					mycp = RPtoCP(myp,pDoc->global,(pDoc->global->coordinate == RADIAL));
					d = sqrt(pow(point.x - mycp.x,2) + pow(point.y - mycp.y,2));
					if(pDoc->global->coordinate == RADIAL) myp = RDtoCT(myp);
					myp = myp + pDoc->global->origin;

					if(pcount == 0 || d <= dmin1) {
						dmin2 = dmin1;
						p2 = p1;
						p1 = myp;
						dmin1 = d;
						pcount++;
					} else if(pcount == 1 || d <= dmin2) {
						p2 = myp;
						dmin2 = d;
						pcount++;
					}
				}
			}
		}

        /*join them */
		if(pcount >= 2) {
			MEMBER* pmem;
			JOINT joint,joint1;
			joint.p = p1;
			joint1.p = p2;
			
			if(pmem = pDoc->AddMember(&joint,&joint1)) {
				pDoc->SetModifiedFlag();
				pDoc->UpdateDrawing();
			}
		}
	}
}
void CmyView::OnRButtonDown(UINT nFlags,CPoint point) {
	
	CView::OnRButtonDown(nFlags,point);
    CmyDocument* pDoc = GetDocument();
	JOINT* psel_joint = 0;

	int found;
	RPoint rpoint = SnapToGrid(point,found,&psel_joint);
		
    if(pDoc->state == SELECT) {
		if(pDoc->DesignResult) {
			MEMBER* mem;
			if(mem = FindMember(rpoint)) {
				CString name;
				pDoc->LoadCaseListBox->GetLBText(pDoc->LoadCaseListBox->GetCurSel(), name);
				CDesignResultDia myDia(AfxGetMainWnd(),mem,pDoc->c_AnalysisCase,this);
				myDia.DoModal();
				pDoc->UpdateDrawing();
				return;
			}
		} else if(pDoc->AnalysisResult) {
			MEMBER* mem;
			if(mem = FindMember(rpoint)) {
				CString name;
				pDoc->LoadCaseListBox->GetLBText(pDoc->LoadCaseListBox->GetCurSel(), name);
				CStressDiagramDia myDia(AfxGetMainWnd(),mem,pDoc->c_AnalysisCase,this);
				myDia.DoModal();
				pDoc->UpdateDrawing();
				return;
			}
			if(force_diagram == IDEFLECTED || 
				force_diagram == IDISPLACEMENT ||
				force_diagram == IREACTIONS) {
				if(psel_joint) {
					CJointReactionsDia myDia(AfxGetMainWnd(),psel_joint,force_diagram);			
					myDia.DoModal();
				    pDoc->UpdateDrawing();
				    return;
				}
			}
		} else {
			JOINT* joint;
			if(joint = FindJoint(rpoint)) {
				CJointInfoDia myDia(AfxGetMainWnd(),joint,&pDoc->loadcases);
				myDia.DoModal();
				pDoc->UpdateDrawing();
				return;
			}
			
			MEMBER* mem;
			if(mem = FindMember(rpoint)) {
				CMemberInfoDia myDia(AfxGetMainWnd(),mem,
					&pDoc->sections,&pDoc->materials,&pDoc->loadcases);
				myDia.DoModal();
				pDoc->UpdateDrawing();
				return;
			}

			SLAB* slab;
			if(slab = FindSlab(rpoint)) {
				CSlabInfoDia myDia(AfxGetMainWnd(),slab,
					&pDoc->asections,&pDoc->materials,&pDoc->loadcases);
				myDia.DoModal();
				pDoc->UpdateDrawing();
				return;
			}
		}
		
	} 

	/*popup menu*/
	CMenu menu;
	if (menu.CreatePopupMenu()) {
		menu.AppendMenu(MF_STRING, ID_ENTER, "Enter");
        menu.AppendMenu(MF_STRING, ID_ESCAPE, "Escape");
		menu.AppendMenu(MF_STRING, ID_DELETE, "Delete");
		menu.AppendMenu(MF_SEPARATOR);
		menu.AppendMenu(MF_STRING, IDM_SYSTEM, "Edit Grid");
		menu.AppendMenu(MF_STRING, IDM_REFRESH, "Refresh");
		
		ClientToScreen(&point);
		menu.TrackPopupMenu(TPM_LEFTALIGN,point.x, point.y,AfxGetMainWnd(),NULL);
		menu.DestroyMenu();
	}
}

void CmyView::OnLButtonUp(UINT nFlags,CPoint point) {
	
	CView::OnLButtonUp(nFlags,point);
	CmyDocument* pDoc = GetDocument();

	int found;
	RPoint rpoint = SnapToGrid(point,found);
	if(!found && (view == VIEW_3D) && (pDoc->state >= ADD_JOINT  && pDoc->state <= ADD_QUAD_AREA))
		return;
	
	if(pDoc->state == PAN)  {
		ReleaseCapture();
		pDoc->nClicks = 0;
		pDoc->state = SELECT;
   	} else if((pDoc->state == RUBBER_ZOOM || pDoc->state == SELECT || pDoc->state == LINE_SELECT) && pDoc->nClicks >= 1) {
		CPoint sel_point = RPtoCP(pDoc->selrp[0]);
		CRect rect(sel_point,prev_mouse_cp);
		rect.NormalizeRect();
		if(rect.Width() && rect.Height()) {
			if(pDoc->state == SELECT) {
				UINT i;
				CPoint pt,pt1;
				JOINT* joint;
				POSITION pos = pDoc->joints.GetHeadPosition();
				while(pos) {
                    joint = &pDoc->joints.GetNext(pos);
					if(!pDoc->select_all && !InPlane(joint->p))
						continue;
                    pt = RPtoCP(joint->p);
					if(rect.PtInRect(pt)) joint->sel = true;
				}
				MEMBER* member;
				pos = pDoc->members.GetHeadPosition();
				while(pos) {
                    member = &pDoc->members.GetNext(pos);
					if(!pDoc->select_all && !InPlane(member->j1->p))
						continue;
					if(!pDoc->select_all && !InPlane(member->j2->p))
						continue;
                    pt = RPtoCP(member->j1->p);
					pt1 = RPtoCP(member->j2->p);
					if(rect.PtInRect(pt) && rect.PtInRect(pt1) ) member->sel = true;
				}
				SLAB* slab;
				pos = pDoc->slabs.GetHeadPosition();
				while(pos) {
                    slab = &pDoc->slabs.GetNext(pos);
					
					for(i = 0;i < slab->NJ;i++) {
						if(!InPlane(slab->jt[i]->p))
							goto BOTTOM;
					}

					slab->sel = true;
					for(i = 0;i < slab->NJ;i++) {
						pt = RPtoCP(slab->jt[i]->p);
						if(rect.PtInRect(pt)) continue;
						slab->sel = false;
						break;
					}
BOTTOM:;
				}
            } else if(pDoc->state == LINE_SELECT) {
				CPoint pt,pt1;
				MEMBER* member;
				POSITION pos = pDoc->members.GetHeadPosition();
				while(pos) {
                    member = &pDoc->members.GetNext(pos);
					if(!pDoc->select_all && !InPlane(member->j1->p))
						continue;
					if(!pDoc->select_all && !InPlane(member->j2->p))
						continue;
					pt = RPtoCP(member->j1->p);
					pt1 = RPtoCP(member->j2->p);
					if(Intersect(sel_point,point,pt,pt1)) 
						member->sel = true;
				}
			} else if(pDoc->state == RUBBER_ZOOM) {
				DOUBLE tscale = min(Corner.x / rect.Width(),Corner.y / rect.Height());
				Center.x = long((Center.x - sel_point.x) * tscale);
				Center.y = long((Center.y - sel_point.y) * tscale);
				Scale *= tscale;
			}


			CmyView* pView;
			POSITION pos = pDoc->GetFirstViewPosition();
			while (pos != NULL) {
				pView = (CmyView*) pDoc->GetNextView(pos);
				pView->DrawingHasChanged = TRUE;
				pView->Invalidate();
			}  
		}
		pDoc->nClicks = 0;
		pDoc->state = SELECT;
		ReleaseCapture();
	}
}
void CmyView::OnMouseMove(UINT nFlags,CPoint point) {
	CView::OnMouseMove(nFlags,point);
	CmyDocument* pDoc = GetDocument();
	JOINT* psel_joint = 0;

	int found = FALSE;
	RPoint rpoint = SnapToGrid(point,found,&psel_joint);

	DrawGridPoint(rpoint,found);
	
	if(pDoc->state == PAN && pDoc->nClicks == 1) {
		if(abs(point.x - prev_mouse_cp.x) > 20 
			|| abs(point.y - prev_mouse_cp.y) > 20) {
		   Center += CSize(point - prev_mouse_cp);
		   prev_mouse_cp = point;
		   DrawingHasChanged = TRUE;
		   Invalidate();
		}
	} else if((pDoc->state == RUBBER_ZOOM || pDoc->state == SELECT) && pDoc->nClicks == 1) {
		DrawTracker(point,RECT_TRACK);
		prev_mouse_cp = point;
    } else if(pDoc->state == LINE_SELECT && pDoc->nClicks == 1) {
		DrawTracker(point,LINE_TRACK);
		prev_mouse_cp = point;
	} else if((pDoc->state >= ADD_MEMBER  && pDoc->state <= ADD_QUAD_AREA) && pDoc->nClicks >= 1) {
		if(!found && view == VIEW_3D)
			return;

		had_tracker = has_tracker;
		has_tracker = TRUE;
		DrawTracker(point,MEMBER_TRACK);
	}

	/*tooltip*/
	BOOL display_tip = FALSE;
	CString tip = "";

	if(pDoc->AnalysisResult && view != VIEW_3D) {
		MEMBER* member;
		SLAB* slab;
		RPoint pr;
		DOUBLE len;
		RPoint nonsnaprp = CPtoRP(point);
		if(force_diagram >= ISX && force_diagram <= IVMAX && (slab = FindSlab(rpoint))) {
			DOUBLE value = FindStress(slab,rpoint);
			tip.Format("%6.3f",value);
			display_tip = TRUE;
		} else if(force_diagram <= IRZ && (member = FindMember(rpoint,&pr,&len))) {
			DOUBLE value,x;
			if(member->is_curved) x = len;
			else x = distance(member->j1->p,rpoint);

			value = member->CalculateForce(x,force_diagram);
			tip.Format("%.3f at %.3f",value,x);
			display_tip = TRUE;
		}
	} else if(found) {
		if(found == SNAP_GRID) tip = "Grid Point";
		else if(found == SNAP_JOINT) tip = "Joint";
        else if(found == SNAP_MEMBERMID) tip = "Midpoint";
		else if(found == SNAP_MEMBER) tip = "Member";
		display_tip = TRUE;
	} 
	
	/*deflection values*/
	if(pDoc->AnalysisResult && force_diagram == IDEFLECTED && psel_joint) {
		tip = "";
		CString str;
		for(int i = IUX;i <= IRZ;i++) {
			str.Format("%s = %.6f ",DirString[i],psel_joint->disps[i]);
			tip += str;
		}
		display_tip = TRUE;
	}

	/*display tip*/
	if(display_tip) {

		//post a mouse leave message
		TRACKMOUSEEVENT tk;
		tk.cbSize = sizeof(tk);
		tk.dwFlags = TME_LEAVE;
		tk.dwHoverTime = 0;
		tk.hwndTrack = m_hWnd;
		_TrackMouseEvent(&tk);

		//update tip
		POINT pt = point;
		ClientToScreen(&pt);
		
        if(point.x + 100 >= Corner.x) pt.x -= 100;
		if(pDoc->AnalysisResult && force_diagram == IDEFLECTED) {
			if(point.y + 100 >= Corner.y) pt.y -= 100;
		} else {
		    if(point.y + 35 >= Corner.y) pt.y -= 35; 
		}
		pt.x += 16;
		pt.y += 16;
		m_toolTip.SendMessage(TTM_TRACKPOSITION, 0, (LPARAM)MAKELONG(pt.x, pt.y));

		m_ti.lpszText = tip.GetBuffer(0);
		m_toolTip.SendMessage(TTM_SETMAXTIPWIDTH, 0, 100);
		m_toolTip.SendMessage(TTM_UPDATETIPTEXT, 0, (LPARAM)&m_ti);	
		m_toolTip.SendMessage(TTM_TRACKACTIVATE, TRUE, (LPARAM)&m_ti);
	} else {
		m_toolTip.SendMessage(TTM_TRACKACTIVATE, FALSE, (LPARAM)&m_ti);
	}
}
LRESULT CmyView::OnMouseLeave(WPARAM wPawam, LPARAM lParam) {
     return m_toolTip.SendMessage(TTM_TRACKACTIVATE, FALSE, (LPARAM)&m_ti);
}
RPoint CmyView::GetJointPoint(JOINT* joint,CmyDocument* pDoc) {
	if(pDoc->AnalysisResult && force_diagram == IDEFLECTED) {
		MATRIX Rt;
		mat_alloc(Rt,3,3);
		joint->GetRotMatrix(Rt);
		transpose(Rt,3,3);

		RPoint p(joint->disps[IUX],joint->disps[IUY],joint->disps[IUZ]);
		DOUBLE disp_scale = (pDoc->disp_scale * pDoc->user_disp_scale * SELSIZE / 3);
		p = GetGlobal(p,Rt);
		p = joint->p + p * disp_scale;
	
		vec_free(Rt);
		return p;
	} else {
		return joint->p;
	}
}
JOINT* CmyView::FindJoint(RPoint rp,BOOL OriginalPos) {
	CmyDocument* pDoc = GetDocument();
	JOINT* joint;
	POSITION posi;
	DOUBLE d;
	CPoint p1,p2;

	posi = pDoc->joints.GetHeadPosition();
	while(posi) {
		joint = &pDoc->joints.GetNext(posi);
		if(!InPlane(joint->p)) 
			continue;

		p1 = RPtoCP(rp);
		if(OriginalPos) p2 = RPtoCP(joint->p);
		else p2 = RPtoCP(GetJointPoint(joint,pDoc));
		d = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
		if(d <= 5)
			return joint;
	}
	return NULL;
}
BOOL CmyView::SelectJoint(RPoint rp) {
	CmyDocument* pDoc = GetDocument();
	JOINT* joint;
	if(joint = FindJoint(rp)) {
		joint->sel = !joint->sel;
		return TRUE;
	}
	return FALSE;
}
BOOL CmyView::IsNear(RPoint rp,RPoint rp1,RPoint rp2) {
	DOUBLE d,d1,d2,n1,n2,q1,q2;
	CPoint pt,p1,p2;

	if(!InPlane(rp1)) 
		return FALSE;
	if(!InPlane(rp2)) 
		return FALSE;
	
	pt = RPtoCP(rp);
	p1 = RPtoCP(rp1);
	p2 = RPtoCP(rp2);
	
	n1 = (p2.x - p1.x) * (pt.x - p1.x) + (p2.y - p1.y) * (pt.y - p1.y);
	n2 = (p1.x - p2.x) * (pt.x - p2.x) + (p1.y - p2.y) * (pt.y - p2.y);
	d = sqrt(pow(p2.x - p1.x,2) + pow(p2.y - p1.y,2));
	d1 = sqrt(pow(pt.x - p1.x,2) + pow(pt.y - p1.y,2));
	d2 = sqrt(pow(pt.x - p2.x,2) + pow(pt.y - p2.y,2));
	q1 = acos((n1)/(d * d1));
	q2 = acos((n2)/(d * d2));
	
	if(q1 < PI / 2
		&& q2 < PI / 2
		&& d1 * sin(q1) <= 5
		) {
		return TRUE;
	}
	return FALSE;
}
MEMBER* CmyView::FindMember(RPoint rp,RPoint* tp,DOUBLE* len) {
	CmyDocument* pDoc = GetDocument();
	MEMBER* mem;
	POSITION posi;
	DOUBLE L;

	posi = pDoc->members.GetHeadPosition();
	while(posi) {
		mem = &pDoc->members.GetNext(posi);
		if(mem->is_curved) {
			RPoint p,p1;
            POSITION pos = mem->curved_list.GetHeadPosition();
			p1 = mem->curved_list.GetNext(pos) + mem->j1->p;
			L = 0;
			while(pos) {
				p = mem->curved_list.GetNext(pos) + mem->j1->p;
				if(IsNear(rp,p1,p)) {
					if(tp) {
						*tp = ptprojection(rp,p1,p);
						*len = L + distance(*tp,p1);
					}
					return mem;
				}
				L += distance(p,p1);
				p1 = p;
			}
		} else {
			RPoint p1 = GetJointPoint(mem->j1,pDoc);
            RPoint p2 = GetJointPoint(mem->j2,pDoc);
			if(IsNear(rp,p1,p2)) {
				if(tp) *tp = ptprojection(rp,p1,p2);
				return mem;
			}
		}
	}
    return NULL;
}
BOOL CmyView::SelectMember(RPoint rp) {
	CmyDocument* pDoc = GetDocument();
	MEMBER* mem;
	if(mem = FindMember(rp)) {
		mem->sel = !mem->sel;
		return TRUE;
	}
	return FALSE;
}

SLAB* CmyView::FindSlab(RPoint rp) {
	CmyDocument* pDoc = GetDocument();
	SLAB* slab;
	POSITION posi;
	CPoint pt,p[NPOLY];
	CRgn rgn;
	UINT j;
	BOOL result;

	posi = pDoc->slabs.GetHeadPosition();
	while(posi) {
		slab = &pDoc->slabs.GetNext(posi);

		for(j = 0;j < slab->NJ;j++) {
			if(!InPlane(slab->jt[j]->p))
				goto BOTTOM;
            p[j] = RPtoCP(slab->jt[j]->p);
		}

		rgn.CreatePolygonRgn(p,slab->NJ,ALTERNATE);
		pt = RPtoCP(rp);
        result = rgn.PtInRegion(pt);
		rgn.DeleteObject();

		if(result)
            return slab;
BOTTOM:;
	}
    return NULL;
}

BOOL CmyView::SelectSlab(RPoint rp) {
	CmyDocument* pDoc = GetDocument();
	SLAB* slab;
	if(slab = FindSlab(rp)) {
		slab->sel = !slab->sel;
		return TRUE;
	}
	return FALSE;
}
void CmyView::DrawArrow(CDC* pDC,RPoint& rp1,RPoint& rp2,int Arrow) {
    CPoint p1,p2;
	p1 = RPtoCP(rp1);
	p2 = RPtoCP(rp2);
	pDC->MoveTo(p1);
	pDC->LineTo(p2);

	if(Arrow >= SARROW) {
		CPoint p3,p4,p5;
		DOUBLE d,q,ARROW;
		d = sqrt(pow(p2.x - p1.x,2) + pow(p2.y - p1.y,2)); 
		if(p2.y - p1.y > 0) q = acos((p2.x - p1.x) / d);
		else q = -acos((p2.x - p1.x) / d);
		ARROW = d / 8; 
		
		if(Arrow >= SARROW) {
			p5.x = p2.x - (p2.x - p1.x) / 5;
			p5.y = p2.y - (p2.y - p1.y) / 5;
			p3.x = p5.x - int(ARROW * sin(q));
			p3.y = p5.y + int(ARROW * cos(q));
			p4.x = p5.x + int(ARROW * sin(q));
			p4.y = p5.y - int(ARROW * cos(q));
			pDC->MoveTo(p2);
			pDC->LineTo(p3);
			pDC->MoveTo(p2);
			pDC->LineTo(p4);
		}
		if(Arrow == DARROW) {
			p5.x = p1.x + (p2.x - p1.x) / 5;
			p5.y = p1.y + (p2.y - p1.y) / 5;
			p3.x = p5.x - int(ARROW * sin(q));
			p3.y = p5.y + int(ARROW * cos(q));
			p4.x = p5.x + int(ARROW * sin(q));
			p4.y = p5.y - int(ARROW * cos(q));
			pDC->MoveTo(p1);
			pDC->LineTo(p3);
			pDC->MoveTo(p1);
			pDC->LineTo(p4);
		}
	}
}
void CmyView::DrawRectangle(CDC* pDC,RPoint& rp1,RPoint& rp2,RPoint& rp3,RPoint& rp4) {
    CPoint pt[4];
	pt[0] = RPtoCP(rp1);
	pt[1] = RPtoCP(rp2);
    pt[2] = RPtoCP(rp3);
	pt[3] = RPtoCP(rp4);
	pDC->Polygon(pt,4);
}
void CmyView::DrawTriangle(CDC* pDC,RPoint& rp1,RPoint& rp2,RPoint& rp3) {
    CPoint pt[3];
	pt[0] = RPtoCP(rp1);
	pt[1] = RPtoCP(rp2);
    pt[2] = RPtoCP(rp3);
	pDC->Polygon(pt,3);
}
void CmyView::DrawText(CDC* pDC,CString str,CPoint* p,CPoint* p1) {
	LOGFONT logFont;
	LOGPEN logpen;
	COLORREF col;
	CPen* pen = pDC->GetCurrentPen();
	if(pen->GetLogPen(&logpen)) col = pDC->SetTextColor(logpen.lopnColor);
	else col = pDC->SetTextColor(TextColor);
    /*font*/
    
    logFont.lfHeight = 14;
    logFont.lfWidth = 0;
    logFont.lfEscapement = 0;
    logFont.lfOrientation = 0;
    logFont.lfWeight = FW_NORMAL;
    logFont.lfItalic = 0;
    logFont.lfUnderline = 0;
    logFont.lfStrikeOut = 0;
    logFont.lfCharSet = ANSI_CHARSET;
    logFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
    logFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
    logFont.lfQuality = PROOF_QUALITY;
    logFont.lfPitchAndFamily = VARIABLE_PITCH | FF_ROMAN;
    strcpy(logFont.lfFaceName, "Times New Roman");

	CPoint p3;
	if(!p1) {
		p3.x = p->x + 5;
		p3.y = p->y + 5;
	} else {
		DOUBLE d,q;
		p3.x = (p->x + p1->x) / 2;
		p3.y = (p->y + p1->y) / 2;
		
		d = sqrt(pow(p1->x - p->x,2) + pow(p1->y - p->y,2)); 
		q = acos((p1->x - p->x) / d);
		p3.x = p3.x - int(16 * cos(q));
		p3.y = p3.y + int(16 * sin(q));
		q = q * 1800 / PI;
		logFont.lfEscapement = LONG(q);
        logFont.lfOrientation = LONG(q);
	}

	CFont font;
	font.CreateFontIndirect(&logFont);
	CFont* oldFont = pDC->SelectObject(&font);
	int bkmode = pDC->SetBkMode(TRANSPARENT);
    pDC->TextOut(p3.x,p3.y,str);
	pDC->SetBkMode(bkmode);
	pDC->SelectObject(oldFont);
	/*end*/
	pDC->SetTextColor(col);
}
void CmyView::DrawOrigin(CDC* pDC) {
	CPen  pen(PS_SOLID,2,AxisColor);
	CPen* pOldPen = pDC->SelectObject(&pen);
	RPoint p1,C(0,0,0);
	DOUBLE ARROWSIZE = 2 * SELSIZE;
	CPoint cp1,cp2;

	if(view == VIEW_3D || view == VIEW_XY || view == VIEW_XZ) {
		p1.set(ARROWSIZE,0,0);
		DrawArrow(pDC,C,p1,SARROW);
		cp1 = RPtoCP(p1);
		p1 = p1 + RPoint(ARROWSIZE,0,0);
		cp2 = RPtoCP(p1);
		DrawText(pDC,"X",&cp1,&cp2);
	}
	if(view == VIEW_3D || view == VIEW_XY || view == VIEW_YZ) {
		p1.set(0,ARROWSIZE,0);
		DrawArrow(pDC,C,p1,SARROW);
		cp1 = RPtoCP(p1);
		p1 = p1 + RPoint(0,ARROWSIZE,0);
		cp2 = RPtoCP(p1);
		DrawText(pDC,"Y",&cp1,&cp2);
	}
	if(view == VIEW_3D || view == VIEW_XZ || view == VIEW_YZ) {
		p1.set(0,0,ARROWSIZE);
		DrawArrow(pDC,C,p1,SARROW);
		cp1 = RPtoCP(p1);
		p1 = p1 + RPoint(0,0,ARROWSIZE);
		cp2 = RPtoCP(p1);
		DrawText(pDC,"Z",&cp1,&cp2);
	}

    pDC->SelectObject(pOldPen);
}
void CmyView::DrawAxis(CDC* pDC,RPoint& p,MATRIX Rt,DOUBLE ARROWSIZE) {
    CPen* pPen,*pOldPen;
	RPoint p1,p2;

	p1 = p;

	pPen = new CPen(PS_SOLID,1,RGB(255,0,0));
	pOldPen = pDC->SelectObject(pPen);
	p2 = RPoint(ARROWSIZE,0,0);
    p2 = GetGlobal(p2,Rt) + p;
	DrawArrow(pDC,p1,p2,SARROW);
	
	delete pPen;
    pPen = new CPen(PS_SOLID,1,RGB(0,255,0));
	pDC->SelectObject(pPen);
	p2 = RPoint(0,ARROWSIZE,0);
    p2 = GetGlobal(p2,Rt) + p;
	DrawArrow(pDC,p1,p2,SARROW);
	
	delete pPen;
    pPen = new CPen(PS_SOLID,1,RGB(0,0,255));
	pDC->SelectObject(pPen);
	p2 = RPoint(0,0,ARROWSIZE);
    p2 = GetGlobal(p2,Rt) + p;
	DrawArrow(pDC,p1,p2,SARROW);
	
	delete pPen;
    pDC->SelectObject(pOldPen);
}

void CmyView::DrawRestraint(CDC* pDC,JOINT* joint,MATRIX Rt) {
	CmyDocument* pDoc = GetDocument();
	int res = joint->restraint;
	RPoint p2;
	DOUBLE RESSIZE = SELSIZE / 2;
	CString str;
	CPoint cp;

	RPoint p = GetJointPoint(joint,pDoc);
	CPen pen(PS_SOLID,2,RestraintColor),*pOldPen;
	pOldPen = pDC->SelectObject(&pen);
		
#define Convert(mp,op) (mp = GetGlobal(mp,Rt) + op)

	if((res & UX) && (res & UZ) && (res & RY)) {
		p2 = RPoint(-RESSIZE,0,0);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UX) && (res & UZ)) {
		p2 = RPoint(-RESSIZE,0,-2 * RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		pDC->LineTo(RPtoCP(p));
		p2 = RPoint(RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UZ)) {
		p2 = RPoint(-0.5 * RESSIZE,0,0);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0.5 * RESSIZE,0,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,-RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0.5 * RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-0.5 * RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,-RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-0.5 * RESSIZE,0,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,-2 * RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UX) || (res & UZ) || (res & RY)) {
		p2 = RPoint(-RESSIZE,0,-RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(RESSIZE,0,-RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(-RESSIZE,0,RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	}

	if((res & UY) && (res & UZ) && (res & RZ)) {
		p2 = RPoint(0,-RESSIZE,0);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UY) && (res & UZ)) {
		p2 = RPoint(0,-RESSIZE,-2 * RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		pDC->LineTo(RPtoCP(p));
		p2 = RPoint(0,RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UZ)) {
		p2 = RPoint(0,-0.5 * RESSIZE,0);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,0.5 * RESSIZE,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,-RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,0.5 * RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-0.5 * RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,-RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-0.5 * RESSIZE,0);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,-2 * RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,-2 * RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	} else if((res & UY) || (res & UZ) || (res & RZ)) {
		p2 = RPoint(0,-RESSIZE,-RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,RESSIZE,-RESSIZE);
		pDC->MoveTo(RPtoCP(Convert(p2,p)));
		p2 = RPoint(0,-RESSIZE,RESSIZE);
		pDC->LineTo(RPtoCP(Convert(p2,p)));
	}

#undef Convert

	if(force_diagram == IREACTIONS && pDoc->AnalysisResult) {
		cp = RPtoCP(p);
		for(int i = 0;i < 6;i++) {
			if(joint->number[i] == INVALID) continue;
			if(res & (1 << i)) {
				if(EQUAL(joint->forces[i],0)) continue;
				cp.y += 14;
				str.Format("%s = %.2f",ForceString[i],joint->forces[i]);
				DrawText(pDC,str,&cp);
			}
		}
	} else if(force_diagram == IDISPLACEMENT && pDoc->AnalysisResult) {
		cp = RPtoCP(p);
		for(int i = 0;i < 6;i++) {
			if(joint->number[i] == INVALID) continue;
			if(!(res & (1 << i))) {
				if(EQUAL(joint->disps[i],0)) continue;
				cp.y += 14;
				str.Format("%s = %.3g",DirString[i],joint->disps[i]);
				DrawText(pDC,str,&cp);
			}
		}
	} else if(mydisplay.Numbering) {
		cp = RPtoCP(p);
		for(int i = 0;i < 6;i++) {
			if(joint->number[i] == INVALID) continue;
			cp.y += 14;
			str.Format("%s = %d",DirString[i],joint->number[i]);
			DrawText(pDC,str,&cp);
		}
	}

    pDC->SelectObject(pOldPen);
}

void CmyView::DrawCircle(CDC* pDC,RPoint& rp,int JSIZE,COLORREF mycol) {
	COLORREF col = pDC->SetBkColor(BackColor);

	CPoint p = RPtoCP(rp);
	CBrush Brush(mycol);
	CBrush* pOldBrush = pDC->SelectObject(&Brush);
	CRect rect(p.x - JSIZE,p.y - JSIZE,p.x + JSIZE,p.y + JSIZE);
	pDC->Ellipse(&rect);
	pDC->SelectObject(pOldBrush);

	pDC->SetBkColor(col);
}

void CmyView::DrawJoint(JOINT& joint, BOOL constrained,CDC* pDC) {
	/*Draw those on the plane*/
	if(!InPlane(joint.p))
		return;

	/*Draw  Loads / restraints*/
	CmyDocument* pDoc = GetDocument(); 
	BOOL usedc = FALSE;
	if(!pDC) {
		pDC = GetDC();
		usedc = TRUE;
	}

	RPoint rp2;
	CPoint p1,p2,p3;
	CString str;
	DOUBLE H;

	/*rotation matrices*/
	MATRIX R,R1,R2;
	mat_alloc(R,3,3);
	mat_alloc(R1,3,3);
    mat_alloc(R2,3,3);
	joint.GetRotMatrix(R);
	transpose(R,3,3);

	/*loads*/
	if(mydisplay.JointLoad && !(force_diagram != INONE && pDoc->AnalysisResult)) {
		CPen Pen(PS_SOLID,1,RGB(255,255,0));
		CPen Pen1(PS_SOLID,1,RGB(255,0,255));
		CPen* pOldPen = pDC->SelectObject(&Pen);

		JLOAD* pload;
		POSITION pos = joint.load.GetHeadPosition();
		while(pos) {
			pload = &joint.load.GetNext(pos);

			if(pload->loadcase != pDoc->c_AnalysisCase)
				continue;

			if(!pload->is_joint_load())
				continue;
			
			/*rotate to local axis*/
			if(pload->system) {
				GetPointRotMatrix(pload->system->rotation,R1);
				multiply(R,R1,R2,3,3,3);
				equ(R1,R2,3,3);
			} else {
				equ(R1,R,3,3);
			}
			/*draw load arrow*/
			for(int i = IUX;i <= IRZ; i++) {
				if(pload->Q[i]) {
					if(i >= IRX) pDC->SelectObject(&Pen1);
					if(pload->Q[i] < 0) H = 3 * SELSIZE;
					else H = -3 * SELSIZE;

					if(i == IUX || i == IRX) rp2 = RPoint(H,0,0);
					else if(i == IUY || i == IRY) rp2 = RPoint(0,H,0);
					else if(i == IUZ || i == IRZ) rp2 = RPoint(0,0,H);

					rp2 = joint.p + GetGlobal(rp2,R1);

					DrawArrow(pDC,rp2,joint.p,SARROW);
					p1 = RPtoCP(rp2);
					if(pload->type == FORCE) str.Format("%.2f",pload->Q[i]);
					else str.Format("%.3g",pload->Q[i]);
					p1 = RPtoCP(rp2);
					p2 = RPtoCP(joint.p);
					DrawText(pDC,str,&p1,&p2);
				}
				
			}
		}
	    pDC->SelectObject(pOldPen);
	}

	/*restraints*/
	if(mydisplay.JointRestraint) {
		DrawRestraint(pDC,&joint,R);
	}

	/*Draw selection*/
	RPoint jointp = GetJointPoint(&joint,pDoc);
	if(joint.sel) {
		COLORREF col = pDC->SetBkColor(BackColor);
		CPen pen(PS_DOT,1,JointColor);
		CPen* pOldPen = pDC->SelectObject(&pen);
		p1 = p2 = RPtoCP(jointp);
		p1.x -= 10;
		p1.y -= 10;
		p2.x += 10;
		p2.y += 10;
		pDC->MoveTo(p1);
		pDC->LineTo(p2);
		p1 = p2 = RPtoCP(jointp);
		p1.x -= 10;
		p1.y += 10;
		p2.x += 10;
		p2.y -= 10;
		pDC->MoveTo(p1);
		pDC->LineTo(p2);
		pDC->SelectObject(pOldPen);
		pDC->SetBkColor(col);
	}

    /*Draw Ellipse*/
	if(mydisplay.JointShow) {
		int size = 4;
		if(constrained && mydisplay.Constraints)
			DrawCircle(pDC,jointp,size,ConstraintColor);
		else 
			DrawCircle(pDC,jointp,size,JointColor);
	}
	
	if(mydisplay.JointLabel) {
		CPen Pen(PS_SOLID,1,RGB(255,0,255));
		CPen* pOldPen = pDC->SelectObject(&Pen);

		CPoint p = RPtoCP(jointp);
		DrawText(pDC,joint.name,&p);

		pDC->SelectObject(pOldPen);
	}

	/*local axis*/
	if(mydisplay.JointAxis) {
		RPoint rpp = jointp;
		DOUBLE SIZE = 2 * SELSIZE;
	    DrawAxis(pDC,rpp,R,SIZE);
	}

	/*free memory*/
	vec_free(R);
	vec_free(R1);
	vec_free(R2);

	if(usedc) {
		ReleaseDC(pDC);
	}
	/*end*/
}
void CmyView::DrawMember(MEMBER& mem,CDC* pDC) {
	/*Draw those on the plane*/
	if(!InPlane(mem.j1->p))
		return;
	if(!InPlane(mem.j2->p))
		return;
	
    /*Draw*/  
	CmyDocument* pDoc = GetDocument();
	BOOL usedc = FALSE;
	if(!pDC) {
		pDC = GetDC();
		usedc = TRUE;
	}

	CPoint p1,p2;
	RPoint rp1,rp2,tp;
	CPen* pPen = new CPen(PS_SOLID,1,MemberColor);
    CPen pen1(PS_SOLID,1,RGB(255,0,255));
	CPen* oldPen = pDC->SelectObject(pPen);
	DOUBLE L = mem.GetLength();
	DOUBLE H;

	/*rotation matrices*/
	MATRIX R,R1,R2;
	mat_alloc(R,3,3);
	mat_alloc(R1,3,3);
    mat_alloc(R2,3,3);
	mem.GetRotMatrix(R,0);
	transpose(R,3,3);

	/*Loads*/
	if(mydisplay.MemberLoad && !(force_diagram != INONE && pDoc->AnalysisResult)) {
		LOAD* pload;
		POSITION pos = mem.load.GetHeadPosition();
		while(pos) {
			pload = &mem.load.GetNext(pos);

			if(pload->loadcase != pDoc->c_AnalysisCase)
				continue;

			if(!pload->is_member_load())
				continue;

			/*rotate to local axis*/
			if(pload->system) {
				GetPointRotMatrix(pload->system->rotation,R1);
				multiply(R,R1,R2,3,3,3);
				equ(R1,R2,3,3);
			} else {
				equ(R1,R,3,3);
			}

			/*plot after rotating*/
			DOUBLE P = pload->P;
			DOUBLE X = pload->x;
			DOUBLE P1 = pload->P1;
			DOUBLE X1 = pload->x1;
			DOUBLE S = (P1 - P) / (X1 - X);
			DOUBLE P0 = P - S * X;

			/*number of divisions*/
			int divcount;
			if(pload->type == CONCENTRATED) {
				divcount = 1;
			} else {
				if(mem.is_curved) divcount = mem.curved_list.GetCount();
				else divcount = 6;
			}

			/*draw*/
			for(int i = 0; i < divcount;i++) {
				
				if(mem.is_curved) {
                    int index = (i == divcount - 1) ? i - 1 : i;
					mem.GetRotMatrix(R,0,index,index + 1);
					transpose(R,3,3);
					
					if(pload->system) {
						multiply(R,R1,R2,3,3,3);
						equ(R1,R2,3,3);
					} else {
						equ(R1,R,3,3);
					}
				}

				if(pload->type == UNIFORM || pload->type == CONCENTRATED) {
					H = (-P / 2) * pDoc->load_scale * pDoc->user_load_scale * SELSIZE;
				} else {
					DOUBLE Px = P0 + S * (X + ((X1 - X) / (divcount - 1)) * i);
					H = (-Px / 2) * pDoc->load_scale * pDoc->user_load_scale * SELSIZE;
				}
				
				if(pload->type == CONCENTRATED) {
					rp1 = mem.FindPoint(X);
				} else if(pload->type == UNIFORM) {
					if(mem.is_curved) {
						rp1 = mem.curved_list.GetAt(mem.curved_list.FindIndex(i)) + mem.j1->p;
					} else {
						DOUBLE x = (L / (divcount - 1)) * i;
						rp1 = mem.FindPoint(x);
					}
				} else {
					DOUBLE x = X + (X1 - X) * (i / float(divcount - 1));
					rp1 = mem.FindPoint(x);
				}
				
				if(pload->dir == IUX || pload->dir == IRX) rp2 = RPoint(H,0,0);
				else if(pload->dir == IUY || pload->dir == IRY) rp2 = RPoint(0,H,0);
				else if(pload->dir == IUZ || pload->dir == IRZ) rp2 = RPoint(0,0,H);
				
				rp2 = rp1 + GetGlobal(rp2,R1);
				
				if(pload->dir <= IUZ) pDC->SelectObject(pPen);
				else pDC->SelectObject(&pen1);
				DrawArrow(pDC,rp2,rp1,SARROW);
				if(i) {
					DrawArrow(pDC,tp,rp2,NOARROW);
				} 
				if(!i || (pload->type == TRAPEZOIDAL && i == divcount - 1)){
					p1 = RPtoCP(rp1);
					p2 = RPtoCP(rp2);
					CString str;
					str.Format("%.2f",!i ? pload->P : pload->P1);
					DrawText(pDC,str,&p1,&p2);
				}
				tp = rp2;
			}
			/*end*/
		}
	}

    /*Draw Member*/
	if(!(pDoc->AnalysisResult && force_diagram == IDEFLECTED)) {
		delete pPen;
		if(mem.sel) pPen = new CPen(PS_DOT,1,MemberColor);
		else pPen = new CPen(PS_SOLID,1,MemberColor);

		pDC->SelectObject(pPen);

		COLORREF col = pDC->SetBkColor(BackColor);

		if(!mydisplay.Extrude || mem.is_curved) {
			rp1 = mem.j1->p;
		    rp2 = mem.j2->p;
			if(mydisplay.MemberRelease && mem.nrelease) {
				rp1 = rp1 + (rp2 - rp1) * (1 / 6.0);
				DrawCircle(pDC,rp1,3,ReleaseColor);
			}
			if(mydisplay.MemberRelease && mem.frelease) {
				rp2 = rp2 - (rp2 - rp1) * (1 / 6.0);
				DrawCircle(pDC,rp2,3,ReleaseColor);
			}
			if(mem.is_curved) {
				POSITION pos = mem.curved_list.GetHeadPosition();
				rp1 = mem.curved_list.GetNext(pos) + mem.j1->p;
				while(pos) {
					rp2 = mem.curved_list.GetNext(pos) + mem.j1->p;
					p1 = RPtoCP(rp1);
					p2 = RPtoCP(rp2);
					pDC->MoveTo(p1);
					pDC->LineTo(p2);
					rp1 = rp2;
				}
			} else {
				p1 = RPtoCP(rp1);
				p2 = RPtoCP(rp2);
				pDC->MoveTo(p1);
				pDC->LineTo(p2);
			}
		} else {
			int i,j,N,ni;
			DOUBLE r,r1,h,b,ratio;
			SECTION* sec[2],*psec;
			RPoint rp[2][2];
			sec[0] = mem.section;
			sec[1] = mem.e_section ? mem.e_section : sec[0];
			
			if(sec[0]->type == CIRCULAR || sec[1]->type == CIRCULAR) N = 24;
			else N = 4;

			CBrush brush;
			brush.CreateStockObject(NULL_BRUSH);
			CBrush* pOldBrush = pDC->SelectObject(&brush);

			for(i = 0;i < N;i ++) {
				ni = (i == N - 1) ? 0 : (i + 1);
				
				for(j = 0;j < 2; j++) {
					psec = sec[j];
					if(psec->type == CIRCULAR) {
						if(mem.e_section) {
							if(j == 0) {
								r = sec[0]->r + (sec[1]->r - sec[0]->r) * mem.start_offset;
							} else {
								r = sec[1]->r + (sec[0]->r - sec[1]->r) * mem.end_offset;
							}
							r1 = r;
						} else {
							r = psec->r;
						    r1 = psec->r;
						}
					} else {
						if(mem.e_section) {
							if(j == 0) {
								h = sec[0]->h + (sec[1]->h - sec[0]->h) * mem.start_offset;
								b = sec[0]->w + (sec[1]->w - sec[0]->w) * mem.start_offset;
							} else {
								h = sec[1]->h + (sec[0]->h - sec[1]->h) * mem.end_offset;
								b = sec[1]->w + (sec[0]->w - sec[1]->w) * mem.end_offset;
							}
							h /= 2;
							b /= 2;
						} else {
							h = psec->h / 2;
							b = psec->w / 2;
						}
						r = r1 = sqrt(h * h + b * b) / 2;
						ratio = (i % (N / 4)) / DOUBLE(N / 4);
						r = r * sqrt(pow(ratio,2) + pow(1 - ratio,2));
                        ratio = (ni % (N / 4)) / DOUBLE(N / 4);
						r1 = r1 * sqrt(pow(ratio,2) + pow(1 - ratio,2));
					}
					rp[j][0].set(j * L,r * cos ((2 * PI * i) / N + PI / 4),r * sin ((2 * PI * i) / N + PI / 4));
					rp[j][1].set(j * L,r1 * cos ((2 * PI * ni) / N + PI / 4),r1 * sin ((2 * PI * ni) / N + PI / 4)); 
					
					rp[j][0] = GetGlobal(rp[j][0],R);
					rp[j][0] = rp[j][0] + mem.j1->p;	
					rp[j][1] = GetGlobal(rp[j][1],R);
					rp[j][1] = rp[j][1] + mem.j1->p;	
				}
				DrawRectangle(pDC,rp[0][0],rp[0][1],rp[1][1],rp[1][0]);
			}
		    pDC->SelectObject(pOldBrush);
		}
		pDC->SetBkColor(col);
	}
	
	/*deflected shape*/
	if(pDoc->AnalysisResult && force_diagram == IDEFLECTED) {
		DrawDeflectedShape(pDC,mem);
	}
	/*local axis*/
	if(mydisplay.MemberAxis) {
		DOUBLE SIZE = distance(mem.j1->p,mem.j2->p) / 10;
		RPoint rpp = (mem.j1->p + mem.j2->p) * 0.5;
	    DrawAxis(pDC,rpp,R,SIZE);
	}

	/*force diagram*/
	if(force_diagram < INONE && pDoc->AnalysisResult) {
		UBMP32 i;
		RPoint trp1,trp2;
		DOUBLE tH;
		CBrush BrushPositive(PositiveColor);
		CBrush BrushNegative(NegativeColor);
		CPen PenPositive(PS_SOLID,1,PositiveColor);
		CPen PenNegative(PS_SOLID,1,NegativeColor);
		CBrush* pOldBrush;
		if(mydisplay.FillDiagram) pOldBrush = pDC->SelectObject(&BrushPositive);
		else pOldBrush = (CBrush*) pDC->SelectStockObject(NULL_BRUSH);

		/*find maximum and minminmum*/
		DOUBLE vmax = 1e-16,vmin = 1e16,imax,imin;
		for(i = 0;i < mem.nDiv;i++) {
			H = mem.forces[i][force_diagram];
			if(H > vmax) {
				vmax = H;
				imax = i;
			}
			if(H < vmin) {
				vmin = H;
				imin = i;
			}
		}

		DOUBLE X;
		DOUBLE scale = pDoc->force_scale[force_diagram] * pDoc->user_force_scale;
		for(i = 0;i < mem.nDiv;i++) {
			H = (mem.forces[i][force_diagram] / 2) * scale * SELSIZE;
			X = mem.station[i];
			if(mem.is_curved) {
				int index = mem.FindCurvedPoint(X);
				mem.GetRotMatrix(R,0,index,index + 1);
				transpose(R,3,3);
			}			
			if(view == VIEW_3D || view == VIEW_YZ) {
				if(force_diagram == IUX || force_diagram == IRX) rp2 = RPoint(0,H,0);
				else if(force_diagram == IUY || force_diagram == IRZ) rp2 = RPoint(0,H,0);
				else if(force_diagram == IUZ || force_diagram == IRY) rp2 = RPoint(0,0,H);
			} else if(view == VIEW_XY) {
				rp2 = RPoint(0,H,0);
			} else if(view == VIEW_XZ) {
				rp2 = RPoint(0,0,H);
			}

			rp1 = mem.FindPoint(X);
			rp2 = GetGlobal(rp2,R) + rp1;
		
			if(i) {
				if(H * tH >= 0) {
					if(H < 0 || (EQUAL(H,0) && tH < 0) ) {
						if(mydisplay.FillDiagram) pDC->SelectObject(&BrushPositive);
						pDC->SelectObject(&PenPositive);
					} else if(H > 0 || (EQUAL(H,0) && tH > 0)) {
						if(mydisplay.FillDiagram) pDC->SelectObject(&BrushNegative);
						pDC->SelectObject(&PenNegative);
					}
					DrawRectangle(pDC,trp1,trp2,rp2,rp1);
				} else {
					RPoint rp3 = trp2 + (rp2 - trp2) * (tH / (tH - H));
					
					if(mydisplay.FillDiagram) pDC->SelectObject(&BrushPositive);
					pDC->SelectObject(&PenPositive);
					if(H > 0) DrawTriangle(pDC,trp1,trp2,rp3);
					else  DrawTriangle(pDC,rp1,rp2,rp3);

					if(mydisplay.FillDiagram) pDC->SelectObject(&BrushNegative);
					pDC->SelectObject(&PenNegative);
					if(H < 0) DrawTriangle(pDC,trp1,trp2,rp3);
					else  DrawTriangle(pDC,rp1,rp2,rp3);
				}
			}
			trp1 = rp1;
            trp2 = rp2;
			tH = H;
			/*show minimum and maximum values*/
			if((i == imax || i == imin) && mydisplay.ShowValues && fabs(mem.forces[i][force_diagram]) >= 0.01) {
				p1 = RPtoCP(rp1);
				p2 = RPtoCP(rp2);
				CString str;
				str.Format("%.2f",mem.forces[i][force_diagram]);
				DrawText(pDC,str,&p2,&p1);
			}
		}

		pDC->SelectObject(pOldBrush);
	}
	/*labels*/
	if(!(pDoc->AnalysisResult && force_diagram == IDEFLECTED)) {
		p1 = RPtoCP(mem.j1->p);
		p2 = RPtoCP(mem.j2->p);
		if(mydisplay.MemberLabel) {
			DrawText(pDC,mem.name,&p1,&p2);
		}
		if(mydisplay.MemberSection) DrawText(pDC,LPCSTR(mem.section->name),&p1,&p2);
		if(mydisplay.MemberMaterial) DrawText(pDC,LPCSTR(mem.section->material->name),&p1,&p2);
	}
	/*clean*/
	delete pPen;
	pDC->SelectObject(oldPen);

	/*free memory*/
	vec_free(R);
	vec_free(R1);
	vec_free(R2);

	if(usedc) {
		ReleaseDC(pDC);
	}
	/*end*/
}
/*
contour colors
*/
void CmyView::InitMinMaxStress() {
	CmyDocument* pDoc = GetDocument();

	/*min and max*/
	POSITION pos;
	SLAB* slab;
	DOUBLE v;
    min_stress = 1e16;
	max_stress = -1e16;
	
	UINT j,findex;
	findex = force_diagram - ISX;
	pos = pDoc->slabs.GetHeadPosition();
	while(pos) {
		slab = &pDoc->slabs.GetNext(pos);
		for(j = 0;j < slab->NJ;j++) {
			v = slab->forces[SLAB::NSTRESS * j + findex];
			if(v > max_stress) max_stress = v;
			if(v < min_stress) min_stress = v;
		}
	}
}
void CmyView::GetEdgeStress(SLAB* slab,DOUBLE* V) {
	CmyDocument* pDoc = GetDocument();
	UINT i,j,findex = force_diagram - ISX;

	for(j = 0;j < slab->NJ;j++) {
		V[j] = slab->forces[SLAB::NSTRESS * j + findex];
		
		if(average_stress) {
			/*average at joints with adjacent slabs*/
			UINT count = 1;
			SLAB* pslab;
			POSITION pos = pDoc->slabs.GetHeadPosition();
			while(pos) {
				pslab = &pDoc->slabs.GetNext(pos);
				if(*slab == *pslab) continue;
				for(i = 0;i < pslab->NJ;i++) {
					if(slab->jt[j] == pslab->jt[i]) {
						count++;
						V[j] += pslab->forces[SLAB::NSTRESS * i + findex];
						break;
					}
				}
			}
			V[j] /= count;
			/*end*/
		} 
	}
}
DOUBLE CmyView::FindStress(SLAB* slab,RPoint rp) {
	/*get local joint coordinates*/
	RPoint p[NPOLY];
    slab->GetLocalJointCoordinates(p,ASECTION::Q4,&rp);

	/*find*/
	DOUBLE value,e,n,x[3],y[3],A,J1,J2,xp0,yp0,be,ce,bn,cn;

	x[0] = p[0].x - p[1].x + p[2].x - p[3].x;
	y[0] = p[0].y - p[1].y + p[2].y - p[3].y;
	x[1] = p[0].x + p[1].x - p[2].x - p[3].x;
	y[1] = p[0].y + p[1].y - p[2].y - p[3].y;
	x[2] = p[0].x - p[1].x - p[2].x + p[3].x;
	y[2] = p[0].y - p[1].y - p[2].y + p[3].y;

	A = 0.5 * ((p[2].x - p[0].x) * (p[3].y - p[1].y) - (p[3].x - p[1].x) * (p[2].y - p[0].y));
	J1 = (p[2].x - p[3].x) * (p[0].y - p[1].y) - (p[0].x - p[1].x) * (p[2].y - p[3].y);
    J2 = (p[1].x - p[2].x) * (p[0].y - p[3].y) - (p[0].x - p[3].x) * (p[1].y - p[2].y);
	xp0 = rp.x;
	yp0 = rp.y;
	be = A - xp0 * y[0] + yp0 * x[0];
    bn = -A - xp0 * y[0] + yp0 * x[0];
	ce = xp0 * y[1] - yp0 * x[1];
	cn = xp0 * y[2] - yp0 * x[2];

	e = 2 * ce / (-sqrt(be * be - 2 * J1 * ce) - be);
	n = 2 * cn / (sqrt(bn * bn + 2 * J2 * cn) - bn);

	/*calculate stress*/
	DOUBLE V[4],Nf[9];
    GetEdgeStress(slab,V);
	slab->GetShapeFunction(Nf,e,n,ASECTION::Q4);
	multiply(V,Nf,value,4);
	
	return value;
}
void CmyView::DrawColorTable(CDC* pDC) {
	CmyDocument* pDoc = GetDocument();

	/*Draw Contor Color table*/
	UINT i;
	DOUBLE v;
	CPoint p1,p2;
	p1.y = p2.y = Corner.y - 22;
	p1.x = 0;
	p2.x = Corner.x;

	CPen pen(PS_SOLID,1,RGB(255,255,255));
	CPen* pold = pDC->SelectObject(&pen);
	pDC->MoveTo(p1);
	pDC->LineTo(p2);
	pDC->SelectObject(pold);


	COLORREF col;
	for(i = 0;i < NCONTORS;i++) {
		col = ContorColors[i];
		
		CBrush brush1(col),*poldb;
		CPen pen1(PS_SOLID,1,col),*poldp;
		poldb = pDC->SelectObject(&brush1);
		poldp = pDC->SelectObject(&pen1);

		p1.x = i * Corner.x / NCONTORS;
		p1.y = Corner.y;
		p2.x = (i + 1) * Corner.x / NCONTORS;
		p2.y = Corner.y - 20;
		pDC->Rectangle(p1.x,p1.y,p2.x,p2.y);

		/*values*/
		pDC->SelectObject(GetStockObject(BLACK_PEN));
		v = min_stress + (max_stress - min_stress) * (DOUBLE(i + 0.5) / NCONTORS);
		p1.y = p2.y;
		CString str;
		str.Format("%d",int(v));
		DrawText(pDC,str,&p1);
		
		pDC->SelectObject(poldb);
		pDC->SelectObject(poldp);
	}
}
void CmyView::DrawSlabMiniBlockStress(CDC* pDC,DOUBLE* V,RPoint* bp) {
	COLORREF col,col1 = RGB(0,0,0);
	UINT i,j,k,t;
	DOUBLE v;
	CPoint p[8];
	RPoint rp;
	int count,next;

	if(fabs(min_stress - max_stress) <= 0.001) t = 1;
	else t = NCONTORS;

	for(k = 0;k < t;k++) {
		i = NCONTORS - k - 1;
		col = ContorColors[i];

		count = 0;
		v = min_stress + (max_stress - min_stress) * DOUBLE(i) / (t - 1);
		for(j = 0;j < 4;j++) {
			next = (j == 3) ? 0 : j + 1;
			if(V[j] < v || EQUAL(v,V[j])) p[count++] = RPtoCP(bp[j]);
			if((V[j] < v && V[next] > v) || (V[j] > v && V[next] < v)) {
				rp = bp[j] + (bp[next] - bp[j]) * ((v - V[j]) / (V[next] - V[j]));
				p[count++] = RPtoCP(rp);
			}			
		}
				
		CBrush brush1(col),*poldb;
		CPen pen1(PS_SOLID,1,col),*poldp;
		poldb = pDC->SelectObject(&brush1);
		poldp = pDC->SelectObject(&pen1);
		
		pDC->Polygon(p,count);
		
		pDC->SelectObject(poldb);
		pDC->SelectObject(poldp);
	}
}
void CmyView::DrawSlabBlockStress(CDC* pDC,SLAB& slab,RPoint* bp) {
	CmyDocument* pDoc = GetDocument();

	UINT i,j,it,jt,findex;
	int cindex[4];
	RPoint rp[4],sp[4];
	DOUBLE e,n,vstep;
	DOUBLE Nf[9],V[9];

	findex = force_diagram - ISX;

	GetEdgeStress(&slab,V);
	
	vstep = (max_stress - min_stress) / (NCONTORS - 1);
	for(i = 0;i < 4;i++) {
		cindex[i] = int(fabs(V[i] / vstep));
	}

	it = 1 + max(abs(cindex[0] - cindex[1]),abs(cindex[2] - cindex[3]));
	jt = 1 + max(abs(cindex[0] - cindex[3]),abs(cindex[1] - cindex[2]));
	it = max(it,4);
	jt = max(it,4);

	DOUBLE Vm[4],e1,n1;
	for(i = 0;i < it;i++) {
		e = -1 + i * 2.0 / it;
		sp[0] = bp[0] + (bp[1] - bp[0]) * (DOUBLE(i) / it);
		sp[1] = bp[3] + (bp[2] - bp[3]) * (DOUBLE(i) / it);
        sp[2] = bp[0] + (bp[1] - bp[0]) * (DOUBLE(i + 1) / it);
		sp[3] = bp[3] + (bp[2] - bp[3]) * (DOUBLE(i + 1) / it);
		for(j = 0;j < jt;j++) {
			n = -1 + j * 2.0 / jt;
			rp[0] = sp[0] + (sp[1] - sp[0]) * (DOUBLE(j) / jt);
			rp[1] = sp[2] + (sp[3] - sp[2]) * (DOUBLE(j) / jt);
            rp[2] = sp[2] + (sp[3] - sp[2]) * (DOUBLE(j + 1) / jt);
			rp[3] = sp[0] + (sp[1] - sp[0]) * (DOUBLE(j + 1) / jt);

			e1 = e;
			n1 = n;
			slab.GetShapeFunction(Nf,e1,n1,ASECTION::Q4);
			multiply(V,Nf,Vm[0],4);
			e1 = -1 + (i + 1) * 2.0 / it;
			n1 = -1 + j * 2.0 / jt;
			slab.GetShapeFunction(Nf,e1,n1,ASECTION::Q4);
			multiply(V,Nf,Vm[1],4);
			e1 = -1 + (i + 1) * 2.0 / it;
			n1 = -1 + (j + 1) * 2.0 / jt;
			slab.GetShapeFunction(Nf,e1,n1,ASECTION::Q4);
			multiply(V,Nf,Vm[2],4);
			e1 = -1 + (i) * 2.0 / it;
			n1 = -1 + (j + 1) * 2.0 / jt;
			slab.GetShapeFunction(Nf,e1,n1,ASECTION::Q4);
			multiply(V,Nf,Vm[3],4);

			DrawSlabMiniBlockStress(pDC,Vm,rp);
		}
	}
	/*draw border*/
	CPoint p[NPOLY];
	UINT NJ = slab.NJ,type = slab.section->type;
	if(type >= ASECTION::Q4 && type <= ASECTION::Q9)
		NJ = 4;

	CPen pen(PS_SOLID,1,RGB(0,0,0));
	CPen* pOld = pDC->SelectObject(&pen);

	for(i = 0;i < NJ;i++)
		p[i] = RPtoCP(bp[i]);
	pDC->MoveTo(p[0]);
	for(i = 1;i < NJ;i++)
		pDC->LineTo(p[i]);
	pDC->LineTo(p[0]);

	pDC->SelectObject(pOld);
}
void CmyView::DrawSlabStress(CDC* pDC,SLAB& slab,MATRIX Rt) {
	RPoint rp[NPOLY];
	for(UINT i = 0;i < slab.NJ;i++)
		rp[i] = slab.jt[i]->p;
	DrawSlabBlockStress(pDC,slab,rp);
}
void CmyView::DrawSlab(SLAB& slab,CDC* pDC) {
	UINT i,type = slab.section->type;
	RPoint rp0;

	/*Draw those on the plane*/
	for(i = 0;i < slab.NJ;i++) {
		if(!InPlane(slab.jt[i]->p))
			return;
	}

	/*centroid*/
	for(i = 0;i < slab.NJ;i++)
		rp0 = rp0 + slab.jt[i]->p;
	rp0 = rp0 * (1.0 / slab.NJ);

	/*rotation matrices*/
	MATRIX R,Rj;
	mat_alloc(R,3,3);
	mat_alloc(Rj,3,3);
	slab.GetRotMatrix(R,0);
	transpose(R,3,3);

	CmyDocument* pDoc = GetDocument();
	BOOL usedc = FALSE;
	if(!pDC) {
		pDC = GetDC();
		usedc = TRUE;
	}
	/*local axis*/
	if(mydisplay.SlabAxis) {
		DOUBLE d,maxd = 0;
		for(i = 0;i < slab.NJ;i++) {
			d = distance(slab.jt[i]->p,rp0);
			if(d > maxd) maxd = d;
		}
		DOUBLE SIZE = maxd / 4;
	    DrawAxis(pDC,rp0,R,SIZE);
	}
    
	/*gdi*/
	CPen *oldPen;
	CPen pen1(PS_SOLID,1,SlabColor);
	CPen pen2(PS_DOT,1,SlabColor);
	COLORREF col;
	
	col = pDC->SetBkColor(BackColor);
	oldPen = pDC->SelectObject(&pen1);

	/*Draw Slab*/
	CPoint p[NPOLY];
	BOOL stress_plot = (force_diagram >= ISX && force_diagram <= IVMAX && pDoc->AnalysisResult);
	
	if(stress_plot)
		DrawSlabStress(pDC,slab,R);

	if(pDoc->AnalysisResult && force_diagram == IDEFLECTED) {
		DOUBLE disp_scale = (pDoc->disp_scale * pDoc->user_disp_scale * SELSIZE / 3);
		RPoint rp[NPOLY],myrp;
		
		for(i = 0;i < slab.NJ;i++) {
			myrp.x = disp_scale * slab.disps[6 * i + IUX];
			myrp.y = disp_scale * slab.disps[6 * i + IUY];
			myrp.z = disp_scale * slab.disps[6 * i + IUZ];
			
			slab.jt[i]->GetRotMatrix(Rj);
			transpose(Rj,3,3);
			myrp = GetGlobal(myrp,Rj);
			
			rp[i] = slab.jt[i]->p + myrp;
			p[i] = RPtoCP(rp[i]);
		}
		
		if(type >= ASECTION::Q4 && type <= ASECTION::Q9) {
			pDC->MoveTo(p[0]);
			if(type >= ASECTION::Q5) pDC->LineTo(p[4]);
			pDC->LineTo(p[1]);
			if(type >= ASECTION::Q6) pDC->LineTo(p[5]);
			pDC->LineTo(p[2]);
			if(type >= ASECTION::Q7) pDC->LineTo(p[6]);
			pDC->LineTo(p[3]);
			if(type >= ASECTION::Q8) pDC->LineTo(p[7]);
			pDC->LineTo(p[0]);
		} else {
			pDC->MoveTo(p[0]);
			for(i = 1;i < slab.NJ;i++)
				pDC->LineTo(p[i]);
			pDC->LineTo(p[0]);
		}
	} else if(!stress_plot) {
		UINT NJ = slab.NJ;
		if(pDoc->AnalysisResult && type >= ASECTION::Q4 && type <= ASECTION::Q9)
			NJ = 4;
		
		pDC->SelectObject(&pen1);
		for(i = 0;i < NJ;i++)
			p[i] = RPtoCP(slab.jt[i]->p + (rp0 - slab.jt[i]->p) * 0.01);
		pDC->MoveTo(p[0]);
		for(i = 1;i < NJ;i++)
			pDC->LineTo(p[i]);
		pDC->LineTo(p[0]);
				
		if(slab.sel) {
			pDC->SelectObject(&pen2);
			for(i = 0;i < NJ;i++)
				p[i] = RPtoCP(slab.jt[i]->p + (rp0 - slab.jt[i]->p) * 0.1);
			pDC->MoveTo(p[0]);
			for(i = 1;i < NJ;i++)
				pDC->LineTo(p[i]);
			pDC->LineTo(p[0]);
		}

		pDC->SelectObject(&pen1);
				
		/*labels*/
		CPoint p1 = RPtoCP(rp0);
		if(mydisplay.SlabLabel) {
			DrawText(pDC,slab.name,&p1);
		}
		if(mydisplay.SlabSection) DrawText(pDC,LPCSTR(slab.section->name),&p1);
		if(mydisplay.SlabMaterial) DrawText(pDC,LPCSTR(slab.section->material->name),&p1);
	}
		
	/*delete*/
	pDC->SelectObject(oldPen);
	pDC->SetBkColor(col);

	/*free memory*/
	vec_free(R);
	vec_free(Rj);

	if(usedc) {
		ReleaseDC(pDC);
	}
	/*end*/
}
void CmyView::DrawDeflectedShape(CDC* pDC,MEMBER& mem) {
	CmyDocument* pDoc = GetDocument();
	CPoint p,p1;
	RPoint rp,rp1;
	CPen* pPen,*oldPen;
	DOUBLE disp_scale = (pDoc->disp_scale * pDoc->user_disp_scale * SELSIZE / 3);

	pPen = new CPen(PS_SOLID,1,MemberColor);
	oldPen = pDC->SelectObject(pPen);


	MATRIX Rt;
	mat_alloc(Rt,3,3);

	/*rotation matrices*/
	DOUBLE L = mem.GetLength(),X;
	if(pDoc->c_AnalysisCase->acase_type != ANALYSISCASE::RESPONSES_CASE &&
		pDoc->c_AnalysisCase->acase_type != ANALYSISCASE::COMBO_CASE
		) {

		mem.GetRotMatrix(Rt,0);
		transpose(Rt,3,3);
		for(UINT i = 0;i < mem.nDiv;i++) {
			RPoint v,rp;
			X = mem.station[i];

			if(mem.is_curved) {
				int index = mem.FindCurvedPoint(X);
				mem.GetRotMatrix(Rt,0,index,index + 1);
				transpose(Rt,3,3);
			}	
			
			mem.CalculateDeflection(mem.station[i],v);
			rp1 = v * disp_scale;
			rp = mem.FindPoint(X);
			rp1 = GetGlobal(rp1,Rt) + rp;

			p = RPtoCP(rp1);
			if(i) pDC->LineTo(p);
			else  pDC->MoveTo(p);
		}
	} else {
		rp = RPoint(mem.j1->disps[IUX],mem.j1->disps[IUY],mem.j1->disps[IUZ]) * disp_scale;
		rp1 = RPoint(mem.j2->disps[IUX],mem.j2->disps[IUY],mem.j2->disps[IUZ]) * disp_scale;

		mem.j1->GetRotMatrix(Rt);
		transpose(Rt,3,3);
		rp = GetGlobal(rp,Rt) + mem.j1->p;

		mem.j2->GetRotMatrix(Rt);
		transpose(Rt,3,3);
		rp1 = GetGlobal(rp1,Rt) + mem.j2->p;

		p = RPtoCP(rp);
		pDC->MoveTo(p);
		p = RPtoCP(rp1);
		pDC->LineTo(p);
	}

	delete pPen;
	pDC->SelectObject(oldPen);

	/*free memory*/
	vec_free(Rt);
	/*end*/
}
void CmyView::OnDeflected() {
    force_diagram = IDEFLECTED;
    SetTitle();
	Invalidate();
}
void CmyView::OnReactions() {
    force_diagram = IREACTIONS;
	SetTitle();
	Invalidate();
}
void CmyView::OnDisplacement() {
    force_diagram = IDISPLACEMENT;
	SetTitle();
	Invalidate();
}
void CmyView::OnShowLoads() {
	force_diagram = INONE;
	SetTitle();
	Invalidate();
}
void CmyView::DrawView(CDC* pDC) {
    CmyDocument* pDoc = GetDocument();
	POSITION pos;
	register int i;

	/*Draw grids*/
	DrawGrids(pDC);

	/*Draw Jonits*/
	JOINT* joint;
	pos = pDoc->joints.GetHeadPosition();
    for (i = 0;i < pDoc->joints.GetCount();i++) {
		joint = &pDoc->joints.GetNext(pos);
		DrawJoint(*joint,joint->cconnect,pDC);
	}

	/*Draw Slabs*/
	if(force_diagram >= ISX && force_diagram <= IVMAX && pDoc->AnalysisResult) {
		InitMinMaxStress();
	}
	pos = pDoc->slabs.GetHeadPosition();
    for (i = 0;i < pDoc->slabs.GetCount();i++)
		DrawSlab(pDoc->slabs.GetNext(pos),pDC);

	if(force_diagram >= ISX && force_diagram <= IVMAX && pDoc->AnalysisResult) {
		DrawColorTable(pDC);
	}

	/*Draw Members*/
	pos = pDoc->members.GetHeadPosition();
    for (i = 0;i < pDoc->members.GetCount();i++)
		DrawMember(pDoc->members.GetNext(pos),pDC);

	/*Draw Origin*/
	if(mydisplay.Origin) {
		DrawOrigin(pDC);
	}

	/*Draw Member Tracker*/
	if((pDoc->state >= ADD_MEMBER  && pDoc->state <= ADD_QUAD_AREA) && pDoc->nClicks >= 1) {
		DrawTracker(prev_mouse_cp,MEMBER_TRACK);
	}

	/*Draw Plane*/
	if(view == VIEW_3D && has_plane) {
		CmyView* pView = (CmyView*) (((CMDIFrameWnd*) AfxGetApp()->m_pMainWnd)->GetActiveFrame())->GetActiveView();
		if(pView->view != VIEW_3D) {
			had_plane = FALSE;
			DrawPlane(pView);
			had_plane = TRUE;
		}
	}
}
void CmyView::DrawPlaneGrid(CDC* pDC,SYSTEM* system,DOUBLE* value,int i1,int i2) {
	int first;
	CPoint p;
	RPoint rp;
	POSITION posx,posy;
	GRIDLIST* grids = system->grid;

	posx = grids[i1].GetHeadPosition(); 
	while(posx) {
		value[i1] = grids[i1].GetNext(posx);
		first = true;
		
		posy = grids[i2].GetHeadPosition(); 
		while(posy) {
			value[i2] = grids[i2].GetNext(posy);
			
			p = RPtoCP(RPoint(value[IUX],value[IUY],value[IUZ]),system,(system->coordinate == RADIAL));

			if(first) {
				pDC->MoveTo(p);
				first = false;
			} else pDC->LineTo(p);
		}
	}
}
void CmyView::DrawGridOnView(CDC* pDC,SYSTEM* system,int srcView,DOUBLE srcPvalue) {
	CPoint p;
	GRIDLIST* grids = system->grid;
	DOUBLE value[3];
	int i1,i2,i3;
	int count = 0;
	int index[12];
	POSITION pos;

	switch(srcView) {
	case VIEW_XY: 
		index[0] = IUX; index[1] = IUY; 
        index[2] = IUY; index[3] = IUX; 
		count = 4;
		break;
	case VIEW_XZ:
		index[0] = IUX; index[1] = IUZ; 
        index[2] = IUZ; index[3] = IUX; 
		count = 4;
		break;
	case VIEW_YZ: 
		index[0] = IUY; index[1] = IUZ; 
        index[2] = IUZ; index[3] = IUY; 
		count = 4;
		break;
	case VIEW_3D: 
		index[0] = IUX; index[1] = IUY; index[2] = IUZ; 
		index[3] = IUX; index[4] = IUZ; index[5] = IUY; 
        index[6] = IUY; index[7] = IUX; index[8] = IUZ;
        index[9] = IUY; index[10] = IUZ; index[11] = IUX;
		count = 12;
		break;
	}
	
	if(srcView == VIEW_3D) {
		for(int i = 0;i < count;i += 3) {
			i1 = index[i];
			i2 = index[i + 1];
			i3 = index[i + 2];
			pos = grids[i1].GetHeadPosition(); 
			while(pos) {
				value[i1] = grids[i1].GetNext(pos);
				DrawPlaneGrid(pDC,system,value,i2,i3);
			}
		}
	} else {
		for(int i = 0;i < count;i += 2) {
			i1 = index[i]; 
			i2 = index[i + 1];
			value[IUX] = srcPvalue;
			value[IUY] = srcPvalue;
			value[IUZ] = srcPvalue;
			DrawPlaneGrid(pDC,system,value,i1,i2);
		}
	}
}

void CmyView::DrawPlane(CmyView* pView) {
	if(pView->view == VIEW_3D)
		return;

	CmyDocument* pDoc = GetDocument();
	CClientDC clientDC(this);
	CPen pen(PS_SOLID,2,AxisColor);
	CPen* old = clientDC.SelectObject(&pen);

	COLORREF col = clientDC.SetBkColor(BackColor);
	int mode = clientDC.SetROP2(R2_XORPEN);

	if(had_plane) {
		DrawGridOnView(&clientDC,pDoc->global,oldview,oldpvalue);
	}
	if(has_plane) {
		DrawGridOnView(&clientDC,pDoc->global,pView->view,pView->pvalue);
		oldview = pView->view;
		oldpvalue = pView->pvalue;
	}

	clientDC.SetBkColor(col);
	clientDC.SetROP2(mode);
	clientDC.SelectObject(old);
}
void CmyView::DrawGrids(CDC* pDC) {
	if(!mydisplay.Grid || (GetDocument()->AnalysisResult && force_diagram == IDEFLECTED))
		return;

	CmyDocument* pDoc = GetDocument();
	CPen* pPen = new CPen(PS_SOLID,1,GridColor);
    CPen* oldPen = pDC->SelectObject(pPen);

	SYSTEM* system;
	POSITION pos = pDoc->systems.GetHeadPosition();
	while(pos) {
        system = &pDoc->systems.GetNext(pos);
		DrawGridOnView(pDC,system,view,pvalue);
	}
	
	/*finish*/
	delete pPen;
	pDC->SelectObject(oldPen);
}
BOOL CmyView::OnEraseBkgnd(CDC* pDC) {
	return FALSE;
}
BOOL CmyView::OnPreparePrinting(CPrintInfo* pInfo) {
	pInfo->SetMaxPage(1);
	return DoPreparePrinting(pInfo);
}
void CmyView::OnPrint(CDC* pDC,CPrintInfo* pInfo) {
	CmyDocument* pDoc = GetDocument();

	/*scale*/
	if(!pDoc->AnalysisResult) 
		pDoc->DetermineScale();

	/*Draw on memory DC and then BitBlt to avoid flicker*/
	CBrush brush(RGB(255,255,255));
	CBrush* pOld = memDC.SelectObject(&brush);
	
	memDC.PatBlt(0,0,memX,memY,PATCOPY);
	DrawView(&memDC);
	DrawingHasChanged = TRUE;
	
	memDC.SelectObject(pOld);
	/*set window/viewport extent*/
	CRect rect;
	GetClientRect (&rect);
	pDC->SetMapMode(MM_ISOTROPIC);
	pDC->SetWindowExt(rect.Width(), rect.Height());
	pDC->SetViewportExt(pInfo->m_rectDraw.Width(),pInfo->m_rectDraw.Height());
	pDC->BitBlt(0,0,memX,memY,&memDC,0,0,SRCCOPY);
	/*end*/
}
void CmyView::OnDraw(CDC* pDC) {
	CmyDocument* pDoc = GetDocument();
	
	/*scale*/
	if(!pDoc->AnalysisResult) 
		pDoc->DetermineScale();

	/*Draw on memory DC and then BitBlt to avoid flicker*/
	if(DrawingHasChanged) {
		CBrush brush(BackColor);
		CBrush* pOld = memDC.SelectObject(&brush);
		
		memDC.PatBlt(0,0,memX,memY,PATCOPY);
		DrawView(&memDC);
		memDC.SelectObject(pOld);

		DrawingHasChanged = FALSE;
	}
	/*BitBlt on screen*/
	pDC->BitBlt(0,0,memX,memY,&memDC,0,0,SRCCOPY);
}

BOOL CmyView::OnSetCursor(CWnd* pWnd,UINT nHitTest,UINT message) {
	CView::OnSetCursor(pWnd,nHitTest,message);
	CmyDocument* pDoc = GetDocument();
	if(pDoc->state == RUBBER_ZOOM) {
		SetCursor(AfxGetApp()->LoadCursor(IDC_RUBCURSOR));
	} else if(pDoc->state == PAN) {
		SetCursor(LoadCursor(NULL,IDC_SIZEALL));
	} else if(pDoc->state >= ADD_JOINT  && pDoc->state <= QUICK_DRAW_MEMBER) {
		SetCursor(LoadCursor(NULL,IDC_UPARROW));
		if(pDoc->nClicks == 1) {
			CFrameWnd* parent = GetParentFrame();
			parent->ActivateFrame();
		}
	} else {
		SetCursor(LoadCursor(NULL,IDC_ARROW));
	}
	return TRUE;
}
void CmyView::OnEscape() {
	CmyDocument* pDoc = GetDocument();
	if((pDoc->state >= ADD_MEMBER  && pDoc->state <= ADD_QUAD_AREA) && pDoc->nClicks >= 1) {
		has_tracker = FALSE;
		DrawTracker(prev_mouse_cp,MEMBER_TRACK);
		had_tracker = FALSE;
		pDoc->nClicks = 0;
	}
	pDoc->state = SELECT;
}
void CmyView::OnAdd(UINT nID) {
	CmyDocument* pDoc = GetDocument();
	if((pDoc->state >= ADD_MEMBER  && pDoc->state <= ADD_QUAD_AREA) && pDoc->nClicks >= 1) {
		has_tracker = FALSE;
		DrawTracker(prev_mouse_cp,MEMBER_TRACK);
		had_tracker = FALSE;
	}
    pDoc->state = STATE(nID - IDM_SELECT);
	pDoc->nClicks = 0;
}
void CmyView::OnSetFocus(CWnd* pOldWnd) {
	CView::OnSetFocus(pOldWnd);
	if(view == VIEW_3D)
		UpdatePlanes(FALSE);
	else
		UpdatePlanes(TRUE);
}
void CmyView::OnKillFocus(CWnd* pOldWnd) {
	CView::OnKillFocus(pOldWnd);
	CmyDocument* pDoc = GetDocument();
	if((pDoc->state >= ADD_MEMBER  && pDoc->state <= ADD_QUAD_AREA) && pDoc->nClicks >= 1) {
		has_tracker = FALSE;
		DrawTracker(prev_mouse_cp,MEMBER_TRACK);
		had_tracker = FALSE;
	}
	UpdatePlanes(FALSE);
}
void CmyView::OnUpdateAddUI(CCmdUI* pCmdUI) {
	CmyDocument* pDoc = GetDocument();
	pCmdUI->SetRadio(pCmdUI->m_nID == unsigned(pDoc->state) + IDM_SELECT);
	if(pCmdUI->m_nID != IDM_SELECT) pCmdUI->Enable(!pDoc->Lock);
}
void CmyView::OnPerspective() {
	CmyDocument* pDoc = GetDocument();
	if(view == VIEW_XY) {
		eye.x = 0;
		eye.y = 0.001;
		perspective_toggle = view;
		view = VIEW_3D;
	} else if(view == VIEW_XZ) {
		eye.x = 0;
		eye.z = 0;
		perspective_toggle = view;
		view = VIEW_3D;
	} else if(view == VIEW_YZ) {
		eye.y = 0;
		eye.z = 0;
		perspective_toggle = view;
		view = VIEW_3D;
	} else {
		if(perspective_toggle == VIEW_3D) {
			eye.x = -eye.x;
			eye.y = -eye.y;
			eye.z = -eye.z;
			perspective_toggle = view;
		    view = VIEW_3D;
		} else {
			view = perspective_toggle;
		}
	}

	CalculateReye();
	SetTitle();
	Invalidate();
}
void CmyView::OnRotate(UINT nID) {
	CmyDocument* pDoc = GetDocument();
	DOUBLE delta = PI / 24;
	if(nID == IDM_RIGHT) reye.z -= delta;
	else if(nID == IDM_LEFT) reye.z += delta;
	else if(nID == IDM_UP)	reye.y -= delta;
	else if(nID == IDM_DOWN) reye.y += delta;
	
	eye.x = reye.x * sin(reye.y) * cos(reye.z);
	eye.y = reye.x * sin(reye.y) * sin(reye.z);
	eye.z = reye.x * cos(reye.y);
	
	if(pDoc->global->coordinate == CARTESIAN)
		eye = eye + gcenter;
	DrawingHasChanged = TRUE;
	Invalidate();
}
void CmyView::OnUpdateRotateUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(view == VIEW_3D);
}

void CmyView::UpdatePlanes(BOOL value) {
	CmyDocument* pDoc = GetDocument();
	CmyView* pView;
	POSITION pos = pDoc->GetFirstViewPosition();
	while (pos != NULL) {
		pView = (CmyView*) pDoc->GetNextView(pos);
		if(pView->view == VIEW_3D) {
			pView->has_plane = value;
			pView->DrawPlane(this);
			pView->had_plane = value;
		}
	} 
}
void CmyView::OnMoveView(UINT nID) {
    CmyDocument* pDoc = GetDocument();
    if(nID == IDM_MOVEVIEWUP) {
		if(position < pDoc->global->grid[NORMGRID[view]].GetCount() - 1) position++;
		else position = 0;
	} else {
		if(position > 0) position--;
		else position = pDoc->global->grid[NORMGRID[view]].GetCount() - 1;
	}
	pvalue = pDoc->global->grid[NORMGRID[view]].GetAt(pDoc->global->grid[NORMGRID[view]].FindIndex(position));
    
	InitView();
    UpdatePlanes(TRUE);
}
void CmyView::OnUpdateMoveViewUI(CCmdUI* pCmdUI) {
	pCmdUI->Enable(view != VIEW_3D);
}
void CmyView::OnChangeView(UINT nID) {
	CmyDocument* pDoc = GetDocument();
	view = VIEW(nID - ID_VIEW_XY);
	perspective_toggle = view;
	if(view != VIEW_3D) {
		GRIDLIST* pgrid = &pDoc->global->grid[NORMGRID[view]];
		if(!pgrid->IsEmpty()) {
			if(view == VIEW_XY) position = pgrid->GetCount() - 1;
			else position = 0;
			pvalue = pgrid->GetAt(pgrid->FindIndex(position));
		}
	}
	InitView();
	UpdatePlanes(TRUE);
}
void CmyView::OnCoordinate(CCmdUI *pCmdUI) {
	CmyDocument* pDoc = GetDocument();
    POINT p;
	GetCursorPos(&p);
	ScreenToClient(&p);
	CRect rect(0,0,Corner.x,Corner.y);
	if(rect.PtInRect(p)) {
		CString str;
		int found;
		RPoint myp = SnapToGrid(CPoint(p),found);
		if(!found && (view == VIEW_3D)) {
			pCmdUI->Enable(); 
			pCmdUI->SetText("");
		} else {
			if(pDoc->global->coordinate == CARTESIAN)
				str.Format("X%.2f Y%.2f Z%.2f",myp.x,myp.y,myp.z);
			else {
				myp = CTtoRD(myp);
				str.Format("R%.2f Q%.2f Z%.2f",myp.x,myp.y,myp.z);
			}
			pCmdUI->Enable(); 
			pCmdUI->SetText(str);
		}
	} else {
		pCmdUI->Enable(); 
		pCmdUI->SetText("");
	}
}
void CmyView::OnDiagrams() {
	CForceDiagramDia ForceDia(AfxGetMainWnd());
	if(force_diagram < INONE)
		ForceDia.dir = force_diagram;
	else
		ForceDia.dir = IRZ;
	ForceDia.fill = mydisplay.FillDiagram;
	ForceDia.value = mydisplay.ShowValues;
    if(ForceDia.DoModal() == IDOK) {
		force_diagram = ForceDia.dir;
        mydisplay.FillDiagram = ForceDia.fill;
		mydisplay.ShowValues = ForceDia.value;
		Invalidate();
		SetTitle();
	}
}
void CmyView::OnSlabDiagrams() {
	CSlabForceDiagramDia ForceDia(AfxGetMainWnd());
	ForceDia.dir = force_diagram;
	ForceDia.average = average_stress;
    if(ForceDia.DoModal() == IDOK) {
		force_diagram = ForceDia.dir;
		average_stress = ForceDia.average;
		Invalidate();
		SetTitle();
	}
}
void CmyView::OnUpdateDiagramsUI(CCmdUI* pCmdUI) {
	CmyDocument* pDoc = GetDocument();
	if(pCmdUI->m_nID == IDM_DEFLECTED)
		pCmdUI->Enable(pDoc->AnalysisResult && pDoc->DisplacementResult);
	else 
        pCmdUI->Enable(pDoc->AnalysisResult);
}
/*
Interaction Chart
*/
BEGIN_MESSAGE_MAP(CChartDia , CDialog)
    ON_WM_PAINT()
    ON_WM_MOUSEMOVE()
	ON_COMMAND_RANGE(IDC_UP,IDC_BIAXIAL,OnPlot)
END_MESSAGE_MAP()

CChartDia::CChartDia(CWnd* parent,DESIGN* tdesign) : 
                          CDialog(IDD_CHART,parent){
	design = tdesign;
	eye.x = 10;
	eye.y = 10;
	eye.z = -7;
	isection = 0;
	type = SURFACE;
	WDIV = design->detail->nWDIV;
    QDIV = design->detail->nQDIV;

	rplist = new RPOINTLIST[WDIV + 1];
	allrplist = new RPOINTLIST*[WDIV + 1];
	for(UINT i = 0;i < WDIV + 1;i++)
	   allrplist[i] = new RPOINTLIST[QDIV + 1];

	Fill(design);
}
CChartDia::~CChartDia() {
	delete[] rplist;
    for(UINT i = 0;i < WDIV + 1;i++)
		delete[] allrplist[i];
	delete[] allrplist;
}
void CChartDia::DoDataExchange(CDataExchange* pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX,IDC_SURFACE,type);
	DDX_Control(pDX,IDC_CHART,c_Drawing);
	DDX_Control(pDX,IDC_INDEX,c_Index);
	DDX_Control(pDX,IDC_CURRENT,c_Current);
}
void CChartDia::Fill(DESIGN* design) {
	RPoint rp;
	DOUBLE Pu,Muz,Muy,x,w,q,inc;
	for(UINT j = 0;j <= WDIV;j++) {
		w = (j * 0.1);
		for(UINT i = 0;i <= QDIV;i++) {
			q = (i / (2.0 * QDIV)) * PI;
			allrplist[j][i].RemoveAll();
			x = 1e6;
			inc = 1;
			while(x > 0) {
				design->FindMatch(x,q,w,Pu,Muz,Muy);
				rp.x = Muz;
				rp.y = Muy;
				rp.z = Pu;
				
				if(rp.z < 0) {
					if(inc < 0.000001)  {
						allrplist[j][i].AddTail(rp);
						break;
					}
					x += inc;
					inc /= 2;
					x -= inc;
					continue;
				}
				if(inc < 0.01) {
					x -= inc;
					continue;
				}
				allrplist[j][i].AddTail(rp);
				if(x < 100)  {
					if(x > 50) inc = 10;
					else if(x > 20) inc = 5;
					else if(x > 10) inc = 1;
					else if(x > 5) inc = min(0.1,inc);
					else if(x > 2) inc = min(0.05,inc);
					else inc = min(0.01,inc);
					
					if(x > inc) 
						x -= inc;
					else  {
						x = 0.0001;
						design->FindMatch(x,q,w,Pu,Muz,Muy);
						rp.x = Muz + (rp.x - Muz) * (-Pu) / (rp.z - Pu);
						rp.y = Muy + (rp.y - Muy) * (-Pu) / (rp.z - Pu);
						rp.z = 0;
						allrplist[j][i].AddTail(rp);
						break;
					}
				} else x /= 2;
			}
		}
	}
}
void CChartDia::OnPlot(UINT nID) {
	if(nID == IDC_UP) {
		isection++;
	} else if(nID == IDC_DOWN) {
		if(isection > 0) isection--;
	} else {
		type = nID - IDC_SURFACE;
	}

	RPoint v,pv,t;
	DOUBLE target;
	UINT i,j;
	POSITION pos;

	for(j = 0;j <= WDIV;j++) {
		rplist[j].RemoveAll();
	}
	
	if(type == SURFACE) {
		if(isection > WDIV) isection = WDIV;
	} else if(type == UNIAXIAL) {
		if(isection >= QDIV) isection = QDIV;
		for(j = 0;j <= WDIV;j++) {
			pos = allrplist[j][isection].GetHeadPosition();
			while(pos) {
				v = allrplist[j][isection].GetNext(pos);
				rplist[j].AddTail(v);
			}
		}
	} else {
        if(isection > 7) isection = 7;
		target = isection * 0.2;
		for(j = 0;j <= WDIV;j++) {
			for(i = 0;i <= QDIV;i++) {
				pos = allrplist[j][i].GetHeadPosition();
				while(pos) {
					pv = v;
					v = allrplist[j][i].GetNext(pos);
					if(v.z <= target) {
						t.x = v.x + (pv.x - v.x) * (target - v.z) / (pv.z - v.z);
						t.y = v.y + (pv.y - v.y) * (target - v.z) / (pv.z - v.z);
						t.z = target;
						rplist[j].AddTail(t);
						break;
					}
				}
			}
		}
	}
	DrawChart();
}
BOOL CChartDia::OnInitDialog() {
	if(!CDialog::OnInitDialog())
		return FALSE;
	c_Drawing.parent = this;
	return TRUE;
}
void CChartDia::OnOK() {
	CDialog::OnOK();
}
void CChartDia::DrawChart(BOOL paint) {
	
	/*paint section drawing*/
	CRect rect;
	CDC* pDC;
	CPaintDC paintDC(&c_Drawing);
	if(paint) pDC = &paintDC;
	else pDC = c_Drawing.GetDC();
	CBrush wbrush(RGB(255,255,255));
	CBrush* oldBrush = pDC->SelectObject(&wbrush);
	c_Drawing.GetClientRect(&rect);
	

	if(type == SURFACE) {
		Center.x = rect.Width() / 2;
		Center.y = rect.Height() - 100;
	} else {
		Center.x = 1;
		Center.y = rect.Height() - 1;
	}
    Corner = rect.BottomRight();

	int H = rect.Height(),W = rect.Width();
	pDC->Rectangle(rect);
	CRgn rgn;
	rgn.CreateRectRgn(0,0,W,H);
	pDC->SelectClipRgn(&rgn);

	/*draw chart*/
	CPen pen(PS_SOLID,2,RGB(0,0,0));
	CPen* oldPen = pDC->SelectObject(&pen);
	POSITION pos;
	RPoint rp,rp1;
	UINT i,j,s;


#define MODIFY() {\
	switch(s) {\
	  case 0: break;\
	  case 1: rp.x = -rp.x; break;\
	  case 2: rp.x = -rp.x; rp.y = -rp.y; break;\
	  case 3: rp.y = -rp.y; break;\
	};\
};

	rp = allrplist[WDIV][0].GetHead();
	Scale = (H - 200) / rp.z;
	Scale = min(Scale , W / max(rp.x,rp.y));

	if(type == SURFACE) {
		for(s = 0;s < 4;s++) {
			for(i = 0;i <= QDIV;i++) {
				pos = allrplist[isection][i].GetHeadPosition();
				if(pos) {
					rp = allrplist[isection][i].GetNext(pos);
					MODIFY();
					pDC->MoveTo(RPtoCP(rp));
					while(pos) {
						rp = allrplist[isection][i].GetNext(pos);
						MODIFY();
						pDC->LineTo(RPtoCP(rp));
					}
				}
			}
			rp = allrplist[isection][0].GetTail();
			MODIFY();
			pDC->MoveTo(RPtoCP(rp));
			for(i = 1;i <= QDIV;i++) {
				rp = allrplist[isection][i].GetTail();
				MODIFY();
				pDC->LineTo(RPtoCP(rp));
			}
		}
	} else {
		if(type == UNIAXIAL) Scale = 1.65 * Scale;
		else Scale = 1.45 * Scale;
	
		/*grids*/
		CPen pen(PS_SOLID,1,RGB(200,200,200));
		CPen pen1(PS_SOLID,1,RGB(200,0,0));
		CPen* oldPen = pDC->SelectObject(&pen);
		for(i = 0;i < 800;i++) {
			DrawGrid(pDC,i);
		}
		pDC->SelectObject(&pen1);
		for(i = 0;i < 80;i++) {
			DrawGrid(pDC,i * 10);
		}
		pDC->SelectObject(oldPen);

        /*chart*/
		for(j = 0;j <= WDIV;j++) {
			pos = rplist[j].GetHeadPosition();
			if(pos) {
				rp = rplist[j].GetNext(pos);
				pDC->MoveTo(RPtoCP(rp));
				while(pos) {
					rp = rplist[j].GetNext(pos);
					pDC->LineTo(RPtoCP(rp));
				}
			}
		}
	}
	
	/*index*/
	CString str;
	if(type == SURFACE) {
		str.Format("w = %2.1f",isection * 0.1);
	} else if(type == UNIAXIAL) {
		str.Format("q = %.2f PI",(isection * 0.5) / QDIV);
	} else {
		str.Format("v = %2.1f",isection * 0.2);
	}
	c_Index.SetWindowText(str);

	pDC->SelectObject(oldPen);
    /*Draw origin*/
	DrawOrigin(pDC);

	/*release dc*/
	c_Drawing.ReleaseDC(pDC);
	rgn.DeleteObject();
}
CPoint CChartDia::RPtoCP(RPoint rp) {
	CPoint p;
	RPoint myp;

	rp.x = 4 * rp.x;
	rp.y = 4 * rp.y;

	if(type == SURFACE) {
		RPoint v;
		DOUBLE d,q;		
		v.set(0,0,1);
		q = +(PI / 2 - myatan(eye.y , eye.x));
		myp = rotate(rp,v,q);
		v.set(1,0,0);
		d = sqrt(eye.x * eye.x + eye.y * eye.y);
		q = -(PI / 2 + myatan(eye.z , d));
		myp = rotate(myp,v,q);
		d = sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
		myp = perspective(myp,-d);
	} else if(type == UNIAXIAL) {
		myp.x = sqrt(rp.x * rp.x + rp.y * rp.y);
		myp.y = rp.z;
	} else {
        myp.x = rp.x;
		myp.y = rp.y;
	}
	
	p.x = Center.x + int(Scale * myp.x);
    p.y = Center.y - int(Scale * myp.y); 
	return p;
}
RPoint CChartDia::CPtoRP(CPoint cp) {
	RPoint myp;
	DOUBLE x = (cp.x - Center.x) / Scale;
	DOUBLE y = (Center.y - cp.y) / Scale;
	DOUBLE q;

	if(type == SURFACE) {
	
	} else if(type == UNIAXIAL) {
        q = (isection * 0.5 * PI) / QDIV;
		myp.x = x * cos(q); 
		myp.y = x * sin(q); 
		myp.z = y;
	} else if(type == BIAXIAL) {
		myp.x = x; 
		myp.y = y; 
		myp.z = isection * 0.2;
	}
	myp.x = myp.x / 4;
	myp.y = myp.y / 4;
	return myp;
}
void CChartDia::DrawArrow(CDC* pDC,RPoint& rp1,RPoint& rp2) {
    CPoint p1,p2;
	p1 = RPtoCP(rp1);
	p2 = RPtoCP(rp2);
	pDC->MoveTo(p1);
	pDC->LineTo(p2);
}
void CChartDia::DrawOrigin(CDC* pDC) {
	CPen* pPen,*pOldPen;
	RPoint C(0,0,0);
	RPoint p1;
	DOUBLE ARROWSIZE = 4;

	pPen = new CPen(PS_SOLID,1,RGB(255,0,0));
	pOldPen = pDC->SelectObject(pPen);
	p1.set(ARROWSIZE,0,0);
	DrawArrow(pDC,C,p1);

	delete pPen;
    pPen = new CPen(PS_SOLID,1,RGB(0,255,0));
	pDC->SelectObject(pPen);
	p1.set(0,ARROWSIZE,0);
	DrawArrow(pDC,C,p1);

	delete pPen;
    pPen = new CPen(PS_SOLID,1,RGB(0,0,255));
	pDC->SelectObject(pPen);
	p1.set(0,0,ARROWSIZE);
	DrawArrow(pDC,C,p1);

	delete pPen;
    pDC->SelectObject(pOldPen);
}

void CChartDia::DrawGrid(CDC* pDC,int i) {
	DOUBLE SIZE = 0.005;
	RPoint rp,rp1;
	if(type == BIAXIAL) {
		rp = RPoint(i * SIZE,0,0);
		rp1 = RPoint(i * SIZE,10,0);
		pDC->MoveTo(RPtoCP(rp));
		pDC->LineTo(RPtoCP(rp1));
		
		rp = RPoint(0,i * SIZE,0);
		rp1 = RPoint(10,i * SIZE,0);
		pDC->MoveTo(RPtoCP(rp));
		pDC->LineTo(RPtoCP(rp1));
	} else {
		rp = RPoint(i * SIZE,0,0);
		rp1 = RPoint(i * SIZE,0,10);
		pDC->MoveTo(RPtoCP(rp));
		pDC->LineTo(RPtoCP(rp1));
		
		rp = RPoint(0,0,i * SIZE);
		rp1 = RPoint(0,10,i * SIZE);
		pDC->MoveTo(RPtoCP(rp));
		pDC->LineTo(RPtoCP(rp1));
	}
}
void CChartDia::OnCoordinate(CPoint& cp) {
	RPoint rp = CPtoRP(cp);
	CString str;
	if(type == SURFACE) {
		str.Format("");
	} else {
		DOUBLE q,w;
		design->vt = rp.z;
		design->ubt = rp.x;
		design->uht = rp.y;
		if(rp.y == 0) q = PI / 2;
		else if(rp.x == 0) q = 0;
		else q = myatan(rp.x , rp.y);
        w = design->FindW(q);
		str.Format("v = %.2f uh = %.2f  ub = %.2f w = %.2f",rp.z,rp.x,rp.y,w);
	}
	c_Current.SetWindowText(str);

}
void CChartDia::OnPaint() {
	CDialog::OnPaint();
	DrawChart(TRUE);
}
BEGIN_MESSAGE_MAP(CChart , CStatic)
	ON_WM_MOUSEMOVE()
END_MESSAGE_MAP()

void CChart::OnMouseMove(UINT nFlags,CPoint point) {
	CStatic::OnMouseMove(nFlags,point);
	parent->OnCoordinate(point);
}

/*
Child Frame
*/
IMPLEMENT_DYNCREATE(CChildFrame, CMDIChildWnd)

BEGIN_MESSAGE_MAP(CChildFrame,CMDIChildWnd)
    ON_WM_CLOSE()
END_MESSAGE_MAP()

BOOL CChildFrame::PreCreateWindow(CREATESTRUCT& cs) {
   cs.style &= ~(LONG)FWS_ADDTOTITLE;
   return CMDIChildWnd::PreCreateWindow(cs);
}
void CChildFrame::OnClose() {
	CMDIChildWnd::OnClose();
	CMDIFrameWnd* pMDIFrameWnd=(CMDIFrameWnd*) AfxGetMainWnd();
	pMDIFrameWnd->MDITile(MDITILE_VERTICAL);
}
/*
Custom tool bar
*/
#define PLACEHOLDER 24

BEGIN_MESSAGE_MAP(CMyToolBar, CToolBar)
    ON_WM_CREATE()
END_MESSAGE_MAP()

int CMyToolBar::OnCreate(LPCREATESTRUCT lpCreateStruct) {
    CToolBar::OnCreate(lpCreateStruct);

    LoadToolBar(IDR_MAINFRAME);
    SetButtonInfo(PLACEHOLDER, IDM_SELECT_LOADCASE, TBBS_SEPARATOR, 100);
    
    CRect rect;
    GetItemRect(PLACEHOLDER, &rect);
    rect.left += 5;
    rect.right -= 5;
    rect.bottom = rect.top + 100;

    c_List.Create(WS_CHILD | WS_VISIBLE | WS_VSCROLL | 
        CBS_DROPDOWNLIST, rect, this, IDM_SELECT_LOADCASE);

    return 0;
}
/*
Main Frame
*/
IMPLEMENT_DYNCREATE(CMainFrame,CMDIFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame , CMDIFrameWnd)
    ON_WM_CREATE()
	ON_WM_SIZE()
	ON_COMMAND_RANGE(IDM_MAINTOOL,IDM_ASSIGNTOOL,OnToolbar)
	ON_UPDATE_COMMAND_UI_RANGE(IDM_MAINTOOL,IDM_ASSIGNTOOL,OnUpdateToolbarUI)
	ON_COMMAND(IDM_RESET_TOOLBARS,ResetToolBars)
	ON_UPDATE_COMMAND_UI(IDM_SELECT_LOADCASE,OnUpdateLoadcaseUI)
END_MESSAGE_MAP()

CToolBar* CMainFrame::GetTool(UINT nID) {
	CToolBar* pTool;
	switch(nID) {
	case IDM_MAINTOOL:pTool = &myToolBar;break;
	case IDM_ANALYSISTOOL:pTool = &myAnalysisToolBar;break;
	case IDM_DESIGNTOOL:pTool = &myDesignToolBar;break;
	case IDM_DRAWTOOL:pTool = &myDrawToolBar;break;
	case IDM_SELECTTOOL:pTool = &mySelectToolBar;break;
	case IDM_VIEWTOOL:pTool = &myViewToolBar;break;
	case IDM_SNAPTOOL:pTool = &mySnapToolBar;break;
	case IDM_EDITTOOL:pTool = &myEditToolBar;break;
	case IDM_ASSIGNTOOL:pTool = &myAssignToolBar;break;
	}
    return pTool;
}
void CMainFrame::OnToolbar(UINT nID) {
	CToolBar* pTool = GetTool(nID);
	ShowControlBar(pTool,!(pTool->IsVisible()),FALSE);
}
void CMainFrame::OnUpdateToolbarUI(CCmdUI* pCmdUI) {
	CToolBar* pTool = GetTool(pCmdUI->m_nID);
	pCmdUI->SetCheck(pTool->IsVisible());
}
void CMainFrame::DockControlBarBottomOf(CToolBar* Bar,CToolBar* LeftOf) {
	CRect rect;
	DWORD dw;
	UINT n;

	RecalcLayout();
	LeftOf->GetWindowRect(&rect);
	rect.OffsetRect(1,1);
	dw = LeftOf->GetBarStyle();
	n = 0;
	n = (dw & CBRS_ALIGN_TOP) ? AFX_IDW_DOCKBAR_TOP : n;
	n = (dw & CBRS_ALIGN_BOTTOM && n == 0) ? AFX_IDW_DOCKBAR_BOTTOM : n;
	n = (dw & CBRS_ALIGN_LEFT && n == 0) ? AFX_IDW_DOCKBAR_LEFT : n;
	n = (dw & CBRS_ALIGN_RIGHT && n == 0) ? AFX_IDW_DOCKBAR_RIGHT : n;

	DockControlBar(Bar,n,&rect);
}

void CMainFrame::ResetToolBars() {

	ShowControlBar(&myToolBar,TRUE,FALSE);
	ShowControlBar(&myAnalysisToolBar,TRUE,FALSE);
	ShowControlBar(&myDesignToolBar,TRUE,FALSE);
	ShowControlBar(&myAssignToolBar,TRUE,FALSE);
	ShowControlBar(&myDrawToolBar,TRUE,FALSE);
	ShowControlBar(&mySelectToolBar,TRUE,FALSE);
	ShowControlBar(&myViewToolBar,TRUE,FALSE);
	ShowControlBar(&mySnapToolBar,TRUE,FALSE);
	ShowControlBar(&myEditToolBar,TRUE,FALSE);

	DockControlBar(&myToolBar,AFX_IDW_DOCKBAR_TOP);
	DockControlBarBottomOf(&myAnalysisToolBar,&myToolBar);
	DockControlBarBottomOf(&myDesignToolBar,&myAnalysisToolBar);

	DockControlBar(&myAssignToolBar,AFX_IDW_DOCKBAR_RIGHT);

	DockControlBar(&myDrawToolBar,AFX_IDW_DOCKBAR_LEFT);
	DockControlBarBottomOf(&mySelectToolBar,&myDrawToolBar);
	DockControlBarBottomOf(&myViewToolBar,&mySelectToolBar);
    DockControlBarBottomOf(&mySnapToolBar,&myViewToolBar);
	DockControlBarBottomOf(&myEditToolBar,&mySnapToolBar);
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct) {

	CMDIFrameWnd::OnCreate(lpCreateStruct);

	/*Enable docking*/
	EnableDocking(CBRS_ALIGN_ANY);

	/*status bar*/
    myStatusBar.Create(this);
    const UINT indicators[] = { 
	   ID_SEPARATOR,
	   ID_COORDINATE
	};
    myStatusBar.SetIndicators(indicators,sizeof(indicators) / sizeof(UINT));

	/*Main toolbar*/
	myToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myToolBar.EnableDocking(CBRS_ALIGN_TOP);
	
	/*Analysis toolbar*/
	myAnalysisToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myAnalysisToolBar.LoadToolBar(IDR_ANALYSISTOOL);
	myAnalysisToolBar.EnableDocking(CBRS_ALIGN_ANY);
	
	/*Design toolbar*/
	myDesignToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myDesignToolBar.LoadToolBar(IDR_DESIGNTOOL);
	myDesignToolBar.EnableDocking(CBRS_ALIGN_ANY);

	/*Assign toolbar*/
	myAssignToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myAssignToolBar.LoadToolBar(IDR_ASSIGNTOOL);
	myAssignToolBar.EnableDocking(CBRS_ALIGN_ANY);

	/*Draw toolbar*/
	myDrawToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_LEFT
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myDrawToolBar.LoadToolBar(IDR_DRAWTOOL);
	myDrawToolBar.EnableDocking(CBRS_ALIGN_ANY);
	
	/*Select toolbar*/
	mySelectToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_LEFT
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	mySelectToolBar.LoadToolBar(IDR_SELECTTOOL);
	mySelectToolBar.EnableDocking(CBRS_ALIGN_ANY);
	
	/*3D view toolbar*/
	myViewToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_LEFT
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myViewToolBar.LoadToolBar(IDR_VIEWTOOL);
	myViewToolBar.EnableDocking(CBRS_ALIGN_ANY);
	
	/*Snap toolbar*/
	mySnapToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_LEFT
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	mySnapToolBar.LoadToolBar(IDR_SNAPTOOL);
	mySnapToolBar.EnableDocking(CBRS_ALIGN_ANY);

	/*Edit toolbar*/
	myEditToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_LEFT
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);
	myEditToolBar.LoadToolBar(IDR_EDITTOOL);
	myEditToolBar.EnableDocking(CBRS_ALIGN_ANY);

	ResetToolBars();
	return 0;
}
void CMainFrame::OnUpdateLoadcaseUI(CCmdUI* pCmdUI) {
	if(MDIGetActive()) pCmdUI->Enable(TRUE);
	else pCmdUI->Enable(FALSE);
}
void CMainFrame::OnSize(UINT nType,int cx,int cy) {
	ResetToolBars();
	CMDIFrameWnd::OnSize(nType,cx,cy);
	MDITile(MDITILE_VERTICAL);
}
/*
Application
*/
CmyApplication myApplication;

BEGIN_MESSAGE_MAP(CmyApplication, CWinApp)
    ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
    ON_COMMAND(ID_FILE_NEW, OnFileNew) 
	ON_COMMAND(ID_FILE_OPEN, OnFileOpen)
	ON_COMMAND_RANGE(IDM_ONE_WINDOWS,IDM_FOUR_WINDOWS,OnWindows)
	ON_COMMAND_RANGE(IDM_IMPORT_EXCEL,IDM_IMPORT_AUTOCAD,OnImport)
	ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()

void CmyApplication::CreateMyView() {

	POSITION pos = GetFirstDocTemplatePosition();
	CDocTemplate* pDocTemplate = GetNextDocTemplate(pos);
	CMDIChildWnd* pNewFrameWnd = ( CMDIChildWnd* ) ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
	CDocument* pDocument = pNewFrameWnd->GetActiveDocument();
	CmyView* pmyView;

	pNewFrameWnd = ( CMDIChildWnd* )(pDocTemplate->CreateNewFrame( pDocument, NULL));
    pDocTemplate->InitialUpdateFrame(pNewFrameWnd, pDocument, TRUE );
	pmyView = (CmyView*) pNewFrameWnd->GetActiveView();
    pmyView->InitView();
}

void CmyApplication::DestroyMyView() {
    CFrameWnd* pNewFrameWnd = ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
	pNewFrameWnd->DestroyWindow();
}
void CmyApplication::OpenViews(UINT count) {
	if(CmyView::nViews < count) {
		while(CmyView::nViews < count) {
			CreateMyView();
		}
	} else if(CmyView::nViews > count) {
		while(CmyView::nViews > count) {
	        DestroyMyView();		
		}
	}

	((CMDIFrameWnd*)m_pMainWnd)->MDITile(MDITILE_VERTICAL);

	UINT id,i = 0;
	CFrameWnd* pNewFrameWnd = ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
	CmyView* pmyView = (CmyView*) pNewFrameWnd->GetActiveView();
	CmyDocument* pDoc = pmyView->GetDocument();
	POSITION pos = pDoc->GetFirstViewPosition();
	while(pos) {
		pmyView = (CmyView*) pDoc->GetNextView(pos);
		id = ID_VIEW_XY + CmyView::CREATION[count - 1 - i];
		pmyView->OnChangeView(id);
		i++;
	}
}
void CmyApplication::OnWindows(UINT nID) {
	UINT count;
	count = nID - IDM_ONE_WINDOWS + 1;
	OpenViews(count);
}
void CmyApplication::OnImport(UINT nID) {
	skip_model = TRUE;
	OnFileNew();
    CFrameWnd* pNewFrameWnd = ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
	CmyDocument* pDoc = (CmyDocument*) pNewFrameWnd->GetActiveDocument();
	pDoc->OnImport(nID);
	skip_model = FALSE;
}
void CmyApplication::OnFileNew() {
	CWinApp::CloseAllDocuments(FALSE);
	CWinApp::OnFileNew();
	CMDIFrameWnd* myframe = ((CMDIFrameWnd*)m_pMainWnd);
	if(myframe->MDIGetActive()) {
		CFrameWnd* pNewFrameWnd = ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
		CmyView* pmyView = (CmyView*) pNewFrameWnd->GetActiveView();
		pmyView->InitView();
		OpenViews(2);
	}
}
void CmyApplication::OnFileOpen() {
 	char strFilter[] = { "Structure files (*.STR)|*.STR|" };
	CFileDialog FileDlg(TRUE, ".STR", NULL, 0, strFilter);
	if( FileDlg.DoModal() != IDOK )	return;

	CString fileName = FileDlg.GetFileName();
	CloseAllDocuments(FALSE);
	OpenDocumentFile(fileName);
	OpenViews(CmyDocument::nViews);
	
	CMDIChildWnd* pNewFrameWnd = ( CMDIChildWnd* ) ((CFrameWnd*) m_pMainWnd)->GetActiveFrame();
	CmyDocument* pDocument = (CmyDocument*) pNewFrameWnd->GetActiveDocument();
	pDocument->InitViews(LPCSTR(fileName));
}
BOOL CmyApplication::InitInstance() {
   	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	AfxEnableControlContainer();

	SetRegistryKey(_T("Structural Analysis"));
	LoadStdProfileSettings(4);

	CMultiDocTemplate* pDocTemplate;
	pDocTemplate = new CMultiDocTemplate(
		IDR_CHILDFRAME,
		RUNTIME_CLASS(CmyDocument),
		RUNTIME_CLASS(CChildFrame),
		RUNTIME_CLASS(CmyView));
	AddDocTemplate(pDocTemplate);

	CMainFrame* pMainFrame = new CMainFrame;
	if (!pMainFrame->LoadFrame(IDR_MAINFRAME))
		return FALSE;
	m_pMainWnd = pMainFrame;

	pMainFrame->ShowWindow(SW_SHOWMAXIMIZED);
	pMainFrame->UpdateWindow();

	GetCurrentDirectory(256,app_path);
	open_file("log.txt",FALSE,TRUE);
    return TRUE;
}
void CmyApplication::OnAppAbout() {
	CDialog AboutDia(IDD_ABOUT,AfxGetMainWnd());
	AboutDia.DoModal();
}
