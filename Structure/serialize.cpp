#include "common.h"
#include "GUI.h"

static const CString dir_str[] = {
	"UX","UY","UZ","RX","RY","RZ","0"
};

static const CString stress_str[] = {
	"SX","SY","SXY","SMAX","SMIN","SVM","SYZ","SXZ","SVMAX",
	"FX","FY","FXY","FMAX","FMIN","FVM",
    "MX","MY","MXY","MMAX","MMIN",
    "VYZ","VXZ","VMAX"
};
/*
Table
*/
int TABLE::current_view;

#define D(x) TABLE::x

TABLE all_tables[NTABLES] = {
	{"General Settings",true,NULL,FALSE,D(MODEL),D(MSettings)},
	{"Color Settings",true,NULL,FALSE,D(MODEL),D(MSettings)},
	{"Materials",true,NULL,FALSE,D(MODEL),D(MProperties)},
	{"Frame sections",true,NULL,FALSE,D(MODEL),D(MProperties)},
	{"Area sections",true,NULL,FALSE,D(MODEL),D(MProperties)},
	{"Systems",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Grid lines",true,NULL,FALSE,D(MODEL),D(GSystems)},
	{"Load cases",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
	{"Modal cases",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
	{"Response history functions",true,NULL,FALSE,D(MODEL),D(MFunctions)},
	{"Response spectrum functions",true,NULL,FALSE,D(MODEL),D(MFunctions)},
	{"Response history cases",true,NULL,FALSE,D(MODEL),D(ARH)},
	{"Response history damping",true,NULL,FALSE,D(MODEL),D(ARH)},
	{"Response history loading",true,NULL,FALSE,D(MODEL),D(ARH)},
	{"Response spectrum cases",true,NULL,FALSE,D(MODEL),D(ARS)},
	{"Response spectrum damping",true,NULL,FALSE,D(MODEL),D(ARS)},
	{"Response spectrum loading",true,NULL,FALSE,D(MODEL),D(ARS)},
    {"Combinations",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
	{"Non-linear cases",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
	{"Buckling cases",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
	{"Joint General",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Joint Coordinates",true,NULL,FALSE,D(MODEL),D(JGeneral)},
	{"Joint Restraints",true,NULL,FALSE,D(MODEL),D(JGeneral)},
	{"Joint Local Axis",true,NULL,FALSE,D(MODEL),D(JGeneral)},
	{"Joint Loads",true,NULL,FALSE,D(MODEL),D(JGeneral)},
	{"Member General",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Member Connectivity",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Member Section",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Member Local Axis",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Member Loads",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Member Releases",true,NULL,FALSE,D(MODEL),D(MGeneral)},
    {"Member Divisions",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Member Curved Points",true,NULL,FALSE,D(MODEL),D(MGeneral)},
	{"Slab General",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Slab Connectivity",true,NULL,FALSE,D(MODEL),D(SGeneral)},
	{"Slab Section",true,NULL,FALSE,D(MODEL),D(SGeneral)},
	{"Slab Local Axis",true,NULL,FALSE,D(MODEL),D(SGeneral)},
	{"Slab Loads",true,NULL,FALSE,D(MODEL),D(SGeneral)},
	{"Constraints",true,NULL,FALSE,D(MODEL),D(Model)}, 
    {"Constrained Joints",true,NULL,FALSE,D(MODEL),D(GConstraints)}, 
    {"Groups",true,NULL,FALSE,D(MODEL),D(Model)}, 
    {"Grouped Joints",true,NULL,FALSE,D(MODEL),D(GGroups)}, 
    {"Grouped Members",true,NULL,FALSE,D(MODEL),D(GGroups)}, 
	{"Grouped Slabs",true,NULL,FALSE,D(MODEL),D(GGroups)}, 
	{"Detailing",true,NULL,FALSE,D(MODEL),D(Model)}, 
	{"Shear Bars",true,NULL,FALSE,D(MODEL),D(GDetailing)}, 
	{"Beam  Bars",true,NULL,FALSE,D(MODEL),D(GDetailing)}, 
	{"Column Bars",true,NULL,FALSE,D(MODEL),D(GDetailing)}, 
	{"Slab Bars",true,NULL,FALSE,D(MODEL),D(GDetailing)}, 

	{"Joint Forces",true,NULL,FALSE,D(ANALYSIS),D(AJoint)},
	{"Joint Displacements",true,NULL,FALSE,D(ANALYSIS),D(AJoint)},
	{"Member Forces",true,NULL,FALSE,D(ANALYSIS),D(AMember)},
	{"Member Displacements",true,NULL,FALSE,D(ANALYSIS),D(AMember)},
	{"Slab Forces",true,NULL,FALSE,D(ANALYSIS),D(ASlab)},
	{"Slab Displacements",true,NULL,FALSE,D(ANALYSIS),D(ASlab)},
	{"Modal shapes",true,NULL,FALSE,D(ANALYSIS),D(AModes)},
	{"Buckling shapes",true,NULL,FALSE,D(ANALYSIS),D(AModes)},
	
	{"Model",true,NULL,FALSE,D(MODEL),-1},
    {"Settings",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Property_Definitions",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Function_Definitions",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Analysis_Cases",true,NULL,FALSE,D(MODEL),D(Model)},
	{"Response History",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},
    {"Response Spectrum",true,NULL,FALSE,D(MODEL),D(MAnalysis_Cases)},

	{"Analysis",true,NULL,FALSE,D(ANALYSIS),-1},
	{"Joint",true,NULL,FALSE,D(ANALYSIS),D(Analysis)},
	{"Member",true,NULL,FALSE,D(ANALYSIS),D(Analysis)},
	{"Slab",true,NULL,FALSE,D(ANALYSIS),D(Analysis)},
    {"Modes",true,NULL,FALSE,D(ANALYSIS),D(Analysis)},

	{"Design",true,NULL,FALSE,D(DESIGN),-1},
	{"Member Design Result",true,NULL,FALSE,D(DESIGN),D(Design)},
};

#undef D

void clear_tables() {
	for (int i = 0; i < NTABLES; i++) {
		all_tables[i].first_record = true;
		all_tables[i].prs = NULL;
		all_tables[i].save = FALSE;
	}
}
/*
Export / Import
*/
static const CString DbNames[] = {"Excel","Access","AutoCad","0"};
static const CString DbExtension[] = {"xls","mdb","dxf","0"};

static CString GetDriver(CString test) {
    char szBuf[2001];
    WORD cbBufMax = 2000;
    WORD cbBufOut;
    char *pszBuf = szBuf;
    CString sDriver;
	
    if (!SQLGetInstalledDrivers(szBuf, cbBufMax, &cbBufOut))
        return "";
    
    do {
        if (strstr(pszBuf, test) != 0) {
            sDriver = CString(pszBuf);
            break;
        }
        pszBuf = strchr(pszBuf, '\0') + 1;
    } while (pszBuf[1] != '\0');
	
    return sDriver;
}

void CmyDocument::OnExport(UINT nID) {
	Export(nID);
}
BOOL CmyDocument::Export(UINT nID,CmyDatabase* mydatabase) {
	int type;
	CDatabase database;
	CString sDriver,sSql,sPath;
	
	if(!mydatabase) {
		type = nID - IDM_EXPORT_EXCEL;
		sPath = GetPathName();
		if(sPath.IsEmpty()) {
			AfxMessageBox("Save file before continuing");
			return FALSE;
		}
		
		int len = sPath.GetLength();
		sPath.SetAt(len - 3,DbExtension[type][0]);
		sPath.SetAt(len - 2,DbExtension[type][1]);
		sPath.SetAt(len - 1,DbExtension[type][2]);
		
		
		if(nID == IDM_EXPORT_AUTOCAD) {
			ExportToAutocad(sPath);
			return TRUE;
		}
		
		sDriver = GetDriver(DbNames[type]);
		if (sDriver.IsEmpty()) {
			AfxMessageBox("Can't find odbc driver.");
			return FALSE;
		}
	}

	/*select tables*/
	if(DesignResult) TABLE::current_view = TABLE::DESIGN;
	else if(AnalysisResult) TABLE::current_view = TABLE::ANALYSIS;
	else TABLE::current_view = TABLE::MODEL;
    CSelectTableDia Dia(AfxGetMainWnd(),all_analysis_cases);
	if(Dia.DoModal() != IDOK)
		return FALSE;
		
	TRY 
	{
		if(!mydatabase) {
			sSql.Format("CREATE_DB=\"%s\" General\0",sPath);
			SQLConfigDataSource(NULL,ODBC_ADD_DSN,sDriver,sSql);
			sSql.Format("DRIVER={%s};DSN='';FIRSTROWHASNAMES=1;READONLY=FALSE;DBQ=%s",sDriver, sPath);
		}
		if(mydatabase ||  database.OpenEx(sSql,CDatabase::noOdbcDialog) ) {
			CSerializer sr;
			if(mydatabase) sr.Set(*mydatabase,TRUE);
			else sr.Set(database,TRUE);
			Serialize(sr);
			return TRUE;
		}
	} 
	CATCH(CDBException, e)
    {
        AfxMessageBox("Database error: " + e->m_strError);
    }
    END_CATCH;
	return FALSE;
}
void CmyDocument::OnImport(UINT nID) {
	Import(nID);
}
BOOL CmyDocument::Import(UINT nID,CmyDatabase* mydatabase) {
	int type;
	CDatabase database;
	CString sDriver,sSql,sPath,sFilter;
	
	if(nID == IDM_IMPORT_EXCEL) sFilter = "Excel files (*.xls)|*.xls|All Files (*.*)|*.*||";
	else if(nID == IDM_IMPORT_ACCESS) sFilter = "Access files (*.mdb)|*.mdb|All Files (*.*)|*.*||";
    else sFilter = "Autocad files (*.dxf)|*.dxf|All Files (*.*)|*.*||";
	
	CFileDialog FileDlg(TRUE,".xls",NULL,0,sFilter);
	
	if(FileDlg.DoModal() == IDOK) {
		sPath = FileDlg.GetPathName();
	} else {
		return FALSE;
	}
	
	type = nID - IDM_IMPORT_EXCEL;
	if(nID == IDM_IMPORT_AUTOCAD) {
		ImportFromAutocad(sPath);
		return TRUE;
	}
	
	sDriver = GetDriver(DbNames[type]);
	if (sDriver.IsEmpty()) {
		AfxMessageBox("Can't find odbc driver.");
		return FALSE;
	}
	
	TRY 
	{
		if(!mydatabase) {
			sSql.Format("DRIVER={%s};DSN='';FIRSTROWHASNAMES=1;READONLY=FALSE;DBQ=%s",sDriver, sPath);
		}
		if(mydatabase || database.OpenEx(sSql,CDatabase::noOdbcDialog) ) {
			CSerializer sr;
			if(mydatabase) sr.Set(*mydatabase,FALSE);
			else sr.Set(database,FALSE);
			Serialize(sr);
			return TRUE;
		}
	}
	CATCH(CDBException, e)
    {
        AfxMessageBox("Database error: " + e->m_strError);
    }
    END_CATCH;
	return FALSE;
}

/*
my Database
*/
void CmyData::GetFieldsInString(CString str,CString& sFields,int& nCount,int start,BOOL skip) {
	BOOL copy = FALSE,wasWhite = FALSE;
	int i,len;
	unsigned char c;
	CString temp;
	len = str.GetLength();
	temp = " ";
	sFields = "";
	nCount = 0;
	for(i = start;i <= len;i++) {
		c = str.GetAt(i);
		if(c == ' ' || c == ',' || c == ')') {
			if(!wasWhite) {
				copy = !copy;
				if(copy || !skip) {
					sFields += temp;
					nCount++;
				}
			}
			temp = " ";
			wasWhite = TRUE;
			if(c == ')') break;
			continue;
		}
		temp += c;
		wasWhite = FALSE;
	}
}
int CmyData::tokenize(CString str, CString* tokens, const CString str2) {
	CString temp = "";
	BOOL wasWhite = TRUE;
	int i,nCount = 0,len = str.GetLength();
	unsigned char c;
	for(i = 0;i <= len;i++) {
		c = str.GetAt(i);
		if(str2.Find(c) != -1 || i == len) {
			if(!wasWhite) {
				if(i == len) temp += c;
				tokens[nCount++] = temp;
			}
			wasWhite = TRUE;
			temp = "";
			if(i == len) break;
			continue;
		} 
		temp += c;
		wasWhite = FALSE;
	}
	return nCount;
}
void CmyDatabase::Create(CString str) {
	CString temp;
	int i,count;
	CmyData data;

	count = str.Find("\"",2);
	for(i = 0;i <= count;i++)
		temp += str.GetAt(i);
	data.name = temp;
	
	count = str.Find("(",0);
	data.GetFieldsInString(str,data.field_names,data.nFields,count + 1,TRUE);
	database.AddTail(data);
}
void CmyDatabase::Insert(CString str) {
	CString temp;
	int i,count;
	CmyData data,*pdata;
	
	count = str.Find("\"",2);
	for(i = 0;i <= count;i++)
		temp += str.GetAt(i);
	data.name = temp;
    pdata = &database.GetAt(database.Find(data));

	CString str1,str2;
	int nCount;
	count = str.Find("(",0);
	data.GetFieldsInString(str,str1,nCount,count + 1,FALSE);
	count = str.Find("(",count + 1);
	data.GetFieldsInString(str,str2,nCount,count + 1,FALSE);
		
	CString fields[50],newfields[50],newvalues[50];
	CmyData::tokenize(pdata->field_names,fields," ");
	CmyData::tokenize(str1,newfields," ");
    CmyData::tokenize(str2,newvalues," ");

	CString newstr = "";
	BOOL found;
	for(int j = 0;j < pdata->nFields;j++) {
		found = FALSE;
		for(i = 0;i < nCount;i++) {
			if(newfields[i] == fields[j]) {
				newstr += newvalues[i];
				found = TRUE;
                break;
			}
		}
		if(!found) newstr += " ";
		newstr += "\n";
	}
	pdata->data.AddTail(newstr);
}
void CmyDatabase::Clear() {
	database.RemoveAll();
}
/*
Export to excel/access/database
*/
static BOOL OpenRecord(CRecordset* prs,CString sSql) {
		TRY {
			prs->Open(CRecordset::forwardOnly,sSql);
			return TRUE;
		}
		CATCH(CDBException, e)
		{
			return FALSE;
		}
		END_CATCH;
}

void CSerializer::CreateTable(CString str) {
	if(type == MYDB) {
		pmydb->Create(str);
	} else if(type == DB) {
		TRY {
			CString sSql = "CREATE TABLE " + str;
			pdb->ExecuteSQL(sSql);
		}
		CATCH(CDBException, e)
		{
			CString sSql;
			int count = str.Find("\"",2);
			for(int i = 0;i <= count;i++)
				sSql += str.GetAt(i);
			sSql= "DROP TABLE " + sSql;
			pdb->ExecuteSQL(sSql);
			sSql = "CREATE TABLE " + str;
			pdb->ExecuteSQL(sSql);
		}
		END_CATCH;
	}
}
void CSerializer::SetStart() {
	if(!is_db()) 
		return;
	TABLE* table = &all_tables[table_index];
	//if(!table->save)
	//	return;
	
	if(!is_storing) {
		CString sSql;
		if(table->first_record) {
			sSql.Format("SELECT * FROM \"%s\"",table->name);
			table->prs = new CRecordset(pdb);
			OpenRecord(table->prs,sSql);
			table->first_record = false;
		}
	} else {
		str_values = "";
		str_name = "";
		str_name_type = "";
	}
	moving_index = 0;
}
void CSerializer::SetEnd() {
	if(!is_db()) 
		return;
	TABLE* table = &all_tables[table_index];
	//if(!table->save)
	//	return;

	if(is_storing) {
		CString sSql;
		if(table->first_record) {
			sSql.Format("\"%s\" (%s)",table->name,str_name_type);
			CreateTable(sSql);
			table->first_record = false;
		}
		sSql.Format("\"%s\" (%s) VALUES (%s)",table->name,str_name,str_values);
		if(type == MYDB) {
			pmydb->Insert(sSql);
		} else if(type == DB) {
			CString str = "INSERT INTO " + sSql;
			pdb->ExecuteSQL(str);
		}
	} else {
		table->prs->MoveNext();
	}
}
/*
start/end writing
*/
#define StartRecord() if(!ar.skip) ar.SetStart()
#define EndRecord() if(!ar.skip) ar.SetEnd()

#define StartSkipRecord() if(!ar.skip) {ar.SetStart();ar.skip = true;}
#define EndSkipRecord() if(ar.skip) {ar.SetEnd();ar.skip = false;}

/*
save
*/
template<class T>
void FormatSql(T& x,CString& str1,CString& str2) {
	str1.Format("%d",x);
	str2 = "INTEGER";
}
void FormatSql(CString& x,CString& str1,CString& str2) {
	str1.Format("'%s'",x);
	str2 = "VARCHAR";
}
void FormatSql(DOUBLE& x,CString& str1,CString& str2) {
	str1.Format("%.4g",x);
	str2 = "DOUBLE";
}
void FormatBool(BOOL* x,CString& str1,CString& str2) {
	str1.Format("%d",*x);
	str2 = "BIT";
}
/*
load
*/
template<class T>
void FormatSql(T& x,CDBVariant& varValue) {
	x = varValue.m_iVal;
}
void FormatSql(CString& x,CDBVariant& varValue) {
	x = *varValue.m_pstring;
}
void FormatSql(DOUBLE& x,CDBVariant& varValue) {
	x = varValue.m_dblVal;
}
void FormatBool(BOOL& x,CDBVariant& varValue) {
	x = varValue.m_boolVal;
}

template<class T>
void _Save(CSerializer& sr, T& x,CString name,int type = -1) {
	if(sr.is_storing) {
		switch(sr.type) {
		case CSerializer::ARCHIVE :
			(*sr.par) << x;
			break;
		case CSerializer::DB :
        case CSerializer::MYDB :
			CString str1,str2;
			if(type == 0) {
				FormatBool((BOOL*)&x,str1,str2);
			} else {
				FormatSql(x,str1,str2);
			}
			if(sr.str_name == "") {
				sr.str_values += (str1);
				sr.str_name += ("'" + name + "' ");
				sr.str_name_type += ("'" + name + "' " + str2);
			} else {
				sr.str_values += ("," + str1);
				sr.str_name += (",'" + name + "' ");
				sr.str_name_type += (",'" + name + "' " + str2);
			}
			break;
		}
	} else {
		switch(sr.type) {
		case CSerializer::ARCHIVE :
			(*sr.par) >> x;
			break;
		case CSerializer::DB :
        case CSerializer::MYDB :
			CDBVariant varValue;
			TABLE* table = &all_tables[sr.table_index];
            table->prs->GetFieldValue(sr.moving_index++,varValue);
			FormatSql(x,varValue);
			break;
		}
	}
}
void _Save(CSerializer& sr, RPoint& x,CString name) {
	_Save(sr,x.x,name + "x");
	_Save(sr,x.y,name + "y");
	_Save(sr,x.z,name + "z");
}
#define Save(sr,x) {                   \
   _Save(sr,x,#x);                     \
};
#define SaveBool(sr,x) {               \
    _Save(sr,x,#x,0);                  \
};
#define SaveN(sr,x,name) {             \
   _Save(sr,x,name);                   \
};
/*
Find item
*/
#define FindItemIndex(index,item,data) {                    \
	index = -1;                                             \
	int found = false;                                      \
    POSITION mpos = (data).GetHeadPosition();               \
	while(mpos) {                                           \
	    index++;                                            \
		if(item == &(data).GetNext(mpos)) {                 \
		    found = true;                                   \
	        break;                                          \
		}                                                   \
	}                                                       \
	if(!found) index = -1;                                  \
};
/*
save pointer
*/
template <class T,class TL>
void _SavePtr(CSerializer& sr,T& item,TL& data,CString name) {                                
	int i1;                                                 
	if(sr.is_storing) {                                     
		FindItemIndex(i1,item,data);                    
		SaveN(sr,i1,name);                                    
	} else {                                                
		SaveN(sr,i1,name);                                        
		if(i1 == -1) item = 0;                              
		else item = &(data).GetAt((data).FindIndex(i1));    
    }                                                       
};
#define SavePtr(sr,item,data) {                             \
   _SavePtr(sr,item,data,#item);                            \
};

/*
save list
*/
#define _SaveList(data,MYDATATYPE,ser_code,header,footer) { \
	    POSITION pos;                                       \
	    UINT i,count;                                       \
		if(ar.is_storing) {                                 \
             count = (data).GetCount();                     \
		     Save(ar,count);                                \
		     pos = (data).GetHeadPosition();                \
			 for(i = 0; i < count;i++) {                    \
			     header;                                    \
			     MYDATATYPE mydat;                          \
                 mydat  = (data).GetNext(pos);              \
			     ser_code;                                  \
				 footer;                                    \
			 }                                              \
		} else {                                            \
		     Save(ar,count);                                \
			 for(i = 0;i < count;i++) {                     \
				 header;                                    \
			     MYDATATYPE mydat;                          \
			     ser_code;                                  \
			     (data).AddTail(mydat);                     \
				 footer;                                    \
			 }                                              \
        }                                                   \
};


#define SaveList0(data,MYDATATYPE) {                               \
	    _SaveList(data,MYDATATYPE,Save(ar,mydat),                  \
        StartRecord(),EndRecord())                                 \
};
#define SaveList1(data,MYDATATYPE) {                               \
	    _SaveList(data,MYDATATYPE,mydat.Serialize(ar),             \
         0,0)                                                      \
};
#define SaveList2(data,MYDATATYPE,param) {                         \
	    _SaveList(data,MYDATATYPE,mydat.Serialize(ar,param),       \
         0,0)													   \
};
#define SavePtrList(data,data_src,MYDATATYPE) {                    \
        _SaveList(data,MYDATATYPE*,SavePtr(ar,mydat,data_src),     \
        StartRecord(),EndRecord())                                 \
};
/*
Serialize units
*/

IMPLEMENT_SERIAL(DISPLAY,CObject,1)
void DISPLAY::Serialize(CSerializer& ar) {
	SaveBool(ar,JointLabel ); 
	SaveBool(ar,JointRestraint ); 
	SaveBool(ar,JointShow ); 
	SaveBool(ar,JointLoad ); 
	SaveBool(ar,JointAxis ); 
	SaveBool(ar,MemberLabel ); 
	SaveBool(ar,MemberSection ); 
	SaveBool(ar,MemberMaterial ); 
	SaveBool(ar,MemberLoad ); 
	SaveBool(ar,MemberAxis ); 
	SaveBool(ar,MemberRelease ); 
	SaveBool(ar,SlabLabel ); 
	SaveBool(ar,SlabSection ); 
	SaveBool(ar,SlabMaterial ); 
	SaveBool(ar,SlabLoad ); 
	SaveBool(ar,SlabAxis ); 
	SaveBool(ar,Origin ); 
	SaveBool(ar,Grid ); 
	SaveBool(ar,FillDiagram ); 
	SaveBool(ar,ShowValues ); 
	SaveBool(ar,Constraints ); 
	SaveBool(ar,Numbering ); 
	SaveBool(ar,Extrude );
}
void DESIGNDATA::Serialize(CSerializer& ar) {
	Save(ar,flags); 
	Save(ar,As); 
	Save(ar,sd); 
	Save(ar,s); 
	_Save(ar,Area[0],"AreaBot"); 
	_Save(ar,Area[1],"AreaTop");
	_Save(ar,diam[0][0],"DiamTop0"); 
	_Save(ar,diam[0][1],"DiamTop1"); 
	_Save(ar,diam[1][0],"DiamBot0"); 
	_Save(ar,diam[1][1],"DiamBot1");
	_Save(ar,count[0][0],"CountTop0"); 
	_Save(ar,count[0][1],"CountTopt1"); 
	_Save(ar,count[1][0],"CountBot0"); 
	_Save(ar,count[1][1],"CountBot1");
	_Save(ar,total[0],"TotalTop"); 
	_Save(ar,total[1],"TotalBot"); 
	_Save(ar,ratios[0],"Ratio0"); 
	_Save(ar,ratios[1],"Ratio1"); 
	_Save(ar,ratios[2],"Ratio2");
}
void DETAILING::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::GDetailing;

	StartRecord();

    SaveBool(ar,single_shear); 
	SaveBool(ar,cut_at_L3); 
	Save(ar,min_spacing); 
	Save(ar,max_layers); 
	Save(ar,le_option); 
	Save(ar,frame_type);
    SaveBool(ar,accidental_eccentricity); 
	SaveBool(ar,individual_buckling); 
	Save(ar,nWDIV); 
	Save(ar,nQDIV);

	EndRecord();

	ar.table_index = TABLE::GDShear;
	SaveList0(shearbarlist,DOUBLE);

	ar.table_index = TABLE::GDBeam;
	SaveList0(beambarlist,DOUBLE);

	ar.table_index = TABLE::GDColumn;
	SaveList0(columnbarlist,DOUBLE);

	ar.table_index = TABLE::GDSlab;
	SaveList0(slabbarlist,DOUBLE);
}
void REBAR::Serialize(CSerializer& ar) {
	REBAR* rb = this;
    Save(ar,rb->design); 
	Save(ar,rb->type); 
	Save(ar,rb->nz); 
	Save(ar,rb->ny); 
	Save(ar,rb->nt); 
	Save(ar,rb->barsize); 
	Save(ar,rb->stirrup_barsize); 
	Save(ar,rb->cover);
}
void MATERIAL::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::PMaterial;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank); 
	Save(ar,type); 
	Save(ar,E); 
	Save(ar,nu); 
	Save(ar,unitweight); 
	Save(ar,density); 
	Save(ar,alphac); 
	Save(ar,G); 
	Save(ar,fck); 
	Save(ar,fctk); 
	Save(ar,fyk); 
	Save(ar,fyks); 
	Save(ar,ftk); 
	Save(ar,ftks); 
	Save(ar,Es);

	EndRecord();
}

void SECTION::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();

	ar.table_index = TABLE::PMSection;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank);
	Save(ar,type); 
	Save(ar,w); 
	Save(ar,h); 
	Save(ar,r); 
	Save(ar,tw); 
	Save(ar,tf); 
	Save(ar,bd); 
	Save(ar,A); 
	Save(ar,Ay); 
	Save(ar,Az); 
	Save(ar,Ix); 
	Save(ar,Iy); 
	Save(ar,Iz); 
	Save(ar,ry); 
	Save(ar,rz); 
	Save(ar,Zy); 
	Save(ar,Zz); 
	Save(ar,Sy); 
	Save(ar,Sz); 
	Save(ar,fA); 
	Save(ar,fAy); 
	Save(ar,fAz); 
	Save(ar,fIx); 
	Save(ar,fIy); 
	Save(ar,fIz);
	rebar.Serialize(ar);
	SavePtr(ar,material,pDoc->materials);

	EndRecord();
}
void ASECTION::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();

	ar.table_index = TABLE::PASection;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank); 
	Save(ar,type); 
	Save(ar,state_type);
	Save(ar,h); 
	rebar.Serialize(ar);
	SavePtr(ar,material,pDoc->materials);

	EndRecord();
}
void ANALYSISCASE::Serialize(CSerializer& ar) {
    Save(ar,name); 
	Save(ar,rank); 
	Save(ar,index); 
	SaveBool(ar,finished); 
	SaveBool(ar,run); 
	Save(ar,acase_type); 
	Save(ar,start_case);
}
void LOADCASE::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::ALoad_cases;

	StartRecord();

	ANALYSISCASE::Serialize(ar);
    Save(ar,swm); 
	Save(ar,type);

	EndRecord();
}
void MODALCASE::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::AModal_cases;

	StartRecord();

	ANALYSISCASE::Serialize(ar);
    Save(ar,minm); 
	Save(ar,maxm); 
	Save(ar,selmode); 
	Save(ar,runmodes); 
	Save(ar,shift); 
	Save(ar,tolerance);

	EndRecord();
}
void SPECFUNC::Serialize(CSerializer& ar,PFUNCTIONLIST flist) {
	StartRecord();

	SavePtr(ar,function,*flist);
	Save(ar,dir); 
	Save(ar,scale);

	EndRecord();
}
void DAMPING::Serialize(CSerializer& ar) {
	StartRecord();

	Save(ar,w1); 
	Save(ar,w2); 
	Save(ar,e1); 
	Save(ar,e2); 
	Save(ar,fm); 
	Save(ar,fk); 
	Save(ar,cdamping); 
	Save(ar,type);
	RPoint mydat;
	Save(ar,mydat);

	EndRecord();

	SaveList0(dvalues.points,RPoint);
}
void RESPONSE::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();
	ANALYSISCASE::Serialize(ar);
	SavePtr(ar,modalcase,pDoc->modalcases);
}
void RESPONSE::Serialize(CSerializer& ar,PFUNCTIONLIST flist) {
	if(acase_type == RESPONSEH_CASE) ar.table_index = TABLE::RH_damping;
	else ar.table_index = TABLE::RS_damping;
	damping.Serialize(ar);

	if(acase_type == RESPONSEH_CASE) ar.table_index = TABLE::RH_loading;
	else ar.table_index = TABLE::RS_loading;
	SaveList2(funclist,SPECFUNC,flist);
}
void RESPONSEHIST::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::RH_cases;

	StartRecord();

	CmyDocument* pDoc = GetMyDocument();
	RESPONSE::Serialize(ar);
	Save(ar,dt); 
	Save(ar,N); 
	Save(ar,ana_type); 
	Save(ar,type); 
	Save(ar,seltime);

	EndRecord();

	RESPONSE::Serialize(ar,&pDoc->rhfunctions);
}
void RESPONSESPEC::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::RS_cases;

	StartRecord();

	CmyDocument* pDoc = GetMyDocument();
	RESPONSE::Serialize(ar);
    Save(ar,modal_comb); 
	Save(ar,dir_comb);

	EndRecord();

	RESPONSE::Serialize(ar,&pDoc->rsfunctions);
}
void BASE_LOAD::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();
	Save(ar,type); 
	Save(ar,casetype);
	SavePtr(ar,loadcase,pDoc->loadcases);
	SavePtr(ar,system,pDoc->systems);
}
void JLOAD::Serialize(CSerializer& ar) {
	StartRecord();

	BASE_LOAD::Serialize(ar);
	for(int i = IUX;i <= IRZ;i++) 
		_Save(ar,Q[i],dir_str[i]);

	EndRecord();
}
void LOAD::Serialize(CSerializer& ar) {
	StartRecord();

	BASE_LOAD::Serialize(ar);
    Save(ar,dir); 
	Save(ar,P); 
	Save(ar,x); 
	Save(ar,P1); 
	Save(ar,x1);

	EndRecord();
}
void LOADCOMBO::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();
	PANALYSISCASELIST* list = pDoc->all_analysis_cases;
	
	POSITION pos;
	int type,case_index;
    ANALYSISCASE* acase;

	StartRecord();

	Save(ar,FS);
	
	if(ar.is_storing) {
        for(type = 0;type < CASETYPES;type++) {
			case_index = 0;
			pos = list[type]->GetHeadPosition();
			while(pos) {
                acase = &list[type]->GetNext(pos);
				if(acase->name == loadcase->name)
					goto END;
				case_index++;
			}
		}
END:
		Save(ar,type); 
		Save(ar,case_index);
	} else {
		Save(ar,type); 
		Save(ar,case_index);
		loadcase = &list[type]->GetAt(list[type]->FindIndex(case_index));
	}

	EndRecord();
}
/*
Joint
*/
static void SaveJointCoor(CSerializer& ar,JOINT* jt) {
	StartRecord();

	Save(ar,jt->name);
	Save(ar,jt->p);

	EndRecord();
}
static void SaveJointRes(CSerializer& ar,JOINT* jt) {
	StartRecord();

	BOOL res[6];
    UINT i;

	Save(ar,jt->name);

	if(ar.is_storing) {
		for(i = IUX;i <= IRZ;i++)
			res[i] = ((jt->restraint >> i) & 1);
	}

	for(i = IUX;i <= IRZ;i++)
		_Save(ar,res[i],dir_str[i],0);

	if(!ar.is_storing) {
		jt->restraint = 0;
		for(i = IUX;i <= IRZ;i++)
			if(res[i]) jt->restraint |= (1 << i);
	}

	EndRecord();
}
static void SaveJointAxis(CSerializer& ar,JOINT* jt) {
	StartRecord();

	Save(ar,jt->name);
	Save(ar,jt->rotation);

	EndRecord();
}

void JOINT::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();


	ar.table_index = TABLE::JGeneral;
	StartRecord();

	Save(ar,name);
	SaveBool(ar,sel);
	Save(ar,constraint); 
	Save(ar,mconnect); 
	Save(ar,cconnect); 
	Save(ar,sconnect); 
#define S(i) _Save(ar,number[i],"numb"#i);
	S(IUX);S(IUY);S(IUZ);S(IRX);S(IRY);S(IRZ);
#undef S

	EndRecord();

	ar.table_index = TABLE::JCoordinates;
	SaveJointCoor(ar,this);

	ar.table_index = TABLE::JRestraints;
	SaveJointRes(ar,this);

	ar.table_index = TABLE::JLocal_axis;
	SaveJointAxis(ar,this);

	ar.table_index = TABLE::JLoads;
	StartSkipRecord();
    Save(ar,name);
	JLOAD ld;
	ld.Serialize(ar);
	EndSkipRecord();

	SaveList1(load,JLOAD);
}
/*
Member
*/
static void SaveMemberConn(CSerializer& ar,MEMBER* mem) {
	CmyDocument* pDoc = GetMyDocument();
	StartRecord();

	Save(ar,mem->name);
	SavePtr(ar,mem->j1,pDoc->joints);
	SavePtr(ar,mem->j2,pDoc->joints);

	EndRecord();
}
static void SaveMemberSec(CSerializer& ar,MEMBER* mem) {
	CmyDocument* pDoc = GetMyDocument();
	StartRecord();

	Save(ar,mem->name);
	SavePtr(ar,mem->section,pDoc->sections);
	SavePtr(ar,mem->e_section,pDoc->sections);

	EndRecord();
}
static void SaveMemberAxis(CSerializer& ar,MEMBER* mem) {
	StartRecord();

	Save(ar,mem->name);
	Save(ar,mem->alpha);

	EndRecord();
}
static void SaveMemberRel(CSerializer& ar,MEMBER* mem) {
	StartRecord();

	Save(ar,mem->name);
	Save(ar,mem->nrelease);
	Save(ar,mem->frelease);

	EndRecord();
}
static void SaveMemberDiv(CSerializer& ar,MEMBER* mem) {
	StartRecord();

	Save(ar,mem->name);
	Save(ar,mem->nDiv);
	Save(ar,mem->nMinDiv);
	Save(ar,mem->nMinFrameDiv);
	SaveBool(ar,mem->DivAtInterim);

	EndRecord();
}
void MEMBER::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();

	ar.table_index = TABLE::MGeneral;
	StartRecord();

	Save(ar,name); 
	SaveBool(ar,sel); 
	Save(ar,start_offset); 
	Save(ar,end_offset); 
	SaveBool(ar,is_curved);
	Save(ar,redistribution_factor); 
	Save(ar,frame_type); 
	Save(ar,le_option); 
#define S(i) _Save(ar,Le[i],"Le"#i); 
	S(0);S(1);
#undef S
#define S(i) _Save(ar,LeCalculated[i],"LeCalculated"#i); 
	S(0);S(1);
#undef S
	Save(ar,EIyy); 
	Save(ar,EIzz);

	EndRecord();

	ar.table_index = TABLE::MConnectivity;
	SaveMemberConn(ar,this);

	ar.table_index = TABLE::MSection;
	SaveMemberSec(ar,this);

	ar.table_index = TABLE::MLocal_axis;
	SaveMemberAxis(ar,this);
	
	ar.table_index = TABLE::MLoads;
	StartSkipRecord();
	Save(ar,name);
	LOAD ld;
	ld.Serialize(ar);
	EndSkipRecord();

	SaveList1(load,LOAD);

	ar.table_index = TABLE::MReleases;
	SaveMemberRel(ar,this);

	ar.table_index = TABLE::MDivisions;
	SaveMemberDiv(ar,this);

	ar.table_index = TABLE::MCurvedMember;
	StartSkipRecord();
	Save(ar,name);
	RPoint mydat;
	Save(ar,mydat);
	EndSkipRecord();

	SaveList0(curved_list,RPoint);
}
/*
Slab
*/
static void SaveSlabConn(CSerializer& ar,SLAB* sla) {
	CmyDocument* pDoc = GetMyDocument();
	for(UINT i = 0;i < sla->NJ;i++) {
		StartRecord();

		Save(ar,sla->name);
		_SavePtr(ar,sla->jt[i],pDoc->joints,"Joint");

        EndRecord();
	}
}
static void SaveSlabSec(CSerializer& ar,SLAB* sla) {
	CmyDocument* pDoc = GetMyDocument();
	StartRecord();

	Save(ar,sla->name);
	SavePtr(ar,sla->section,pDoc->asections);

	EndRecord();
}
static void SaveSlabAxis(CSerializer& ar,SLAB* sla) {
	StartRecord();

	Save(ar,sla->name);
	Save(ar,sla->alpha);

	EndRecord();
}

void SLAB::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();

	ar.table_index = TABLE::SGeneral;
	StartRecord();

	Save(ar,name); 
	SaveBool(ar,sel); 
	Save(ar,NJ); 
	Save(ar,nDivx); 
	Save(ar,nDivy); 
	SaveBool(ar,edge_constraint);

	EndRecord();

	ar.table_index = TABLE::SConnectivity;
	SaveSlabConn(ar,this);

	ar.table_index = TABLE::SSection;
	SaveSlabSec(ar,this);

	ar.table_index = TABLE::SLocal_axis;
	SaveSlabAxis(ar,this);

	ar.table_index = TABLE::SLoads;
	StartSkipRecord();
	Save(ar,name);
	LOAD ld;
	ld.Serialize(ar);
	EndSkipRecord();

	SaveList1(load,LOAD);
}

#define AddHeader() {             \
	StartRecord();               \
	Save(ar,name);                \
	Save(ar,mydat);               \
	EndRecord();                 \
}
void CONSTRAINT::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();
	int mydat = -1;

	ar.table_index = TABLE::GConstraints;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank); 
	Save(ar,type); 
	Save(ar,constraint); 
	SaveBool(ar,isauto); 
	Save(ar,tolerance);
	SavePtr(ar,axis,pDoc->systems);

	EndRecord();

	ar.table_index = TABLE::GC_Joints;
    AddHeader();
	SavePtrList(jointlist,pDoc->joints,JOINT);
}
void GROUP::Serialize(CSerializer& ar) {
	CmyDocument* pDoc = GetMyDocument();
	int mydat = -1;

	ar.table_index = TABLE::GGroups;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank);
	Save(ar,color); 

	EndRecord();

	ar.table_index = TABLE::GR_Joints;
	AddHeader();
	SavePtrList(jointlist,pDoc->joints,JOINT);

	ar.table_index = TABLE::GR_Members;
	AddHeader();
	SavePtrList(memberlist,pDoc->members,MEMBER);

	ar.table_index = TABLE::GR_Slabs;
	AddHeader();
	SavePtrList(slablist,pDoc->slabs,SLAB);
}
void COMB_TYPE::Serialize(CSerializer& ar) {
	SaveList1(loadlist,LOADCOMBO);
}
#undef AddHeader

static void AddHeader(CSerializer& ar) {
	DOUBLE FS = 0;
	int type = 0,case_index = 0;
	Save(ar,FS);
	Save(ar,type);
	Save(ar,case_index);
}
void COMBINATION::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::ACombinations;

	StartRecord();

	AddHeader(ar);
	ANALYSISCASE::Serialize(ar);
	Save(ar,designcombo); 
	Save(ar,selcombo);

	EndRecord();

	COMB_TYPE::Serialize(ar);
}
void NLCASE::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::ANonLinear_cases;

	StartRecord();

	AddHeader(ar);
	ANALYSISCASE::Serialize(ar);
	Save(ar,nSteps); 
	Save(ar,nIteration); 
	Save(ar,tolerance); 
	Save(ar,save_stiffness);

	EndRecord();

	COMB_TYPE::Serialize(ar);
}
void BUCKLINGCASE::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::ABuckling_cases;

	StartRecord();

	AddHeader(ar);
	ANALYSISCASE::Serialize(ar);
	Save(ar,nmodes); 
	Save(ar,runmodes); 
	Save(ar,tolerance);

	EndRecord();

	COMB_TYPE::Serialize(ar);
}
void SYSTEM::Serialize(CSerializer& ar) {
	ar.table_index = TABLE::GSystems;

	StartRecord();

	Save(ar,name); 
	Save(ar,rank); 
	Save(ar,coordinate); 
	Save(ar,origin); 
	Save(ar,rotation); 

	EndRecord();


	ar.table_index = TABLE::GGridLines;

	DOUBLE position;
	POSITION pos;
	for(int dir = IUX;dir <= IUZ;dir++) {
		if(!ar.is_db()) {
			SaveList0(grid[dir],DOUBLE);
		} else {
			pos = grid[dir].GetHeadPosition();
			while(pos)  {
				position = grid[dir].GetNext(pos);
				
				StartRecord();
				Save(ar,name);
				Save(ar,dir);
				Save(ar,position);
				EndRecord();
			}
		}
	}
}
void FUNCTION::Serialize(CSerializer& ar) {
	StartRecord();

	Save(ar,name); 
	Save(ar,rank); 
	Save(ar,period); 
	Save(ar,amplitude); 
	Save(ar,cycles); 
	Save(ar,divisions); 
	Save(ar,type); 
	Save(ar,Ag); 
	Save(ar,soil_type);
	RPoint mydat;
	Save(ar,mydat);

	EndRecord();

	SaveList0(points,RPoint);
}
/*
Serialize colors
*/
void CmyView::SerializeColors(CSerializer& ar) {
	ar.table_index = TABLE::GColors;

	StartRecord();

	Save(ar,BackColor);
	Save(ar,JointColor);
	Save(ar,MemberColor);
	Save(ar,SlabColor);
	Save(ar,HighlightColor);
	Save(ar,SelectionColor);
	Save(ar,RestraintColor);
	Save(ar,ConstraintColor);
	Save(ar,AxisColor);
	Save(ar,ReleaseColor);
	Save(ar,TextColor);
	Save(ar,GridColor);
	Save(ar,PositiveColor);
	Save(ar,NegativeColor);

#define S(i) _Save(ar,ContorColors[i],"Contor"#i);

	S(0);S(1);S(2);S(3);S(4);
	S(5);S(6);S(7);S(8);S(9);
	S(10);S(11);S(12);S(13);S(14);

#undef S
	
	EndRecord();
};
void CmyDocument::SerializeGeneral(CSerializer& ar) {
	ar.table_index = TABLE::GSettings;

	StartRecord();

	if(ar.is_storing) {
		UINT& numViews = CmyView::nViews; 
		Save(ar,numViews);
	} else {
		UINT& numViews = CmyDocument::nViews;
		Save(ar,numViews);
	}
	SaveBool(ar,AnalysisResult); 
	SaveBool(ar,DesignResult); 
	SaveBool(ar,Lock);
	Save(ar,SLAB::TotalNumber); 
	Save(ar,MEMBER::TotalNumber); 
	Save(ar,JOINT::TotalNumber);
	Save(ar,ELEMENT::TotalLoadCases);
	Save(ar,static_dofs); 
	Save(ar,dynamic_dofs);
	Save(ar,user_force_scale); 
	Save(ar,user_disp_scale); 
	Save(ar,user_load_scale);

	EndRecord();
}
/*
Serialize whole document
*/
void CmyDocument::Serialize(CArchive& ar) {
	CSerializer sr;
	if(ar.IsStoring()) sr.Set(ar,TRUE);
	else sr.Set(ar,FALSE);
	Serialize(sr);
}
void CmyDocument::Serialize(CSerializer& ar) {
	SerializeGeneral(ar);
	DisplacementResult = AnalysisResult;
    CmyView::SerializeColors(ar);
	SaveList1(materials,MATERIAL);
	SaveList1(sections,SECTION);
	SECTION::defmaterial = &materials.GetHead();
	SaveList1(asections,ASECTION);
	ASECTION::defmaterial = &materials.GetHead();
	SaveList1(systems,SYSTEM);
	global = &systems.GetHead();
	SaveList1(loadcases,LOADCASE);
	SaveList1(modalcases,MODALCASE);
	RESPONSE::defmodalcase = &modalcases.GetHead();
	ar.table_index = TABLE::RH_function;
	SaveList1(rhfunctions,FUNCTION);
	RESPONSE::deffunction = &rhfunctions.GetHead();
	ar.table_index = TABLE::RS_function;
	SaveList1(rsfunctions,FUNCTION);
	FUNCTION::defspectrumtype = FALSE;
	SaveList1(responsecases,RESPONSEHIST);
	SaveList1(responsespecs,RESPONSESPEC);
	SaveList1(combinations,COMBINATION);
	SaveList1(nlcases,NLCASE);
	SaveList1(bucklingcases,BUCKLINGCASE);
    SaveList1(joints,JOINT);
	SaveList1(members,MEMBER);
	MEMBER::defsection = &sections.GetHead();
    SaveList1(slabs,SLAB);
	SLAB::defsection = &asections.GetHead();
    SaveList1(constraints,CONSTRAINT);
	SaveList1(groups,GROUP);
	detailing.Serialize(ar);
	if(AnalysisResult) {
		SerializeAnalysisResult(ar);
	}
	/*fill cases*/
	if(!ar.is_storing) {
		if(AnalysisResult) {
			FillAnalysisCases();
			DetermineScale();
		} else {
			FillCases();
		}
	}
	/*close record sets*/
	if(ar.is_db()) {
		TABLE* table;
		for(int i = 0;i < NTABLES;i++) {
			table = &all_tables[i];
			if(table->prs  && table->save) {
				table->prs->Close();
				delete table->prs;
				table->prs = NULL;
				table->first_record = false;
			}
		}
	}
}
/*
Serialize analysis result
*/
static void AddHeader(CSerializer& ar,ANALYSISCASE* pcase,ENTITY* pent) {
	if(!ar.is_db()) return;
	_Save(ar,pent->name,"Name");
	_Save(ar,pcase->name,"Case");
	_Save(ar,pcase->acase_type,"Type");
	_Save(ar,pcase->selected,"Selected");
}

void JOINT::SerializeAnalysis(CSerializer& ar,int index,ANALYSISCASE* pcase) {
	UINT i;
	if(!ar.is_storing)
		AllocAnalysis(index);

	ar.table_index = TABLE::JForces;

	StartRecord();
	AddHeader(ar,pcase,this);
	for(i = 0;i < 6; i++) {
		_Save(ar,forces_all[index][i],dir_str[i]); 
	}
	EndRecord();

	ar.table_index = TABLE::JDisplacements;
	StartRecord();
	AddHeader(ar,pcase,this);
	for(i = 0;i < 6; i++) {
		_Save(ar,disps_all[index][i],dir_str[i]);
	}
	EndRecord();
}

void MEMBER::SerializeAnalysis(CSerializer& ar,int index,ANALYSISCASE* pcase) {
	CmyDocument* pDoc = GetMyDocument();
	UINT i,j;
	if(!ar.is_storing)
		AllocAnalysis(index);

	ar.table_index = TABLE::MForces;
	for(i = 0;i < nDiv; i++) {
		StartRecord();

		AddHeader(ar,pcase,this);
		_Save(ar,station[i],"Station");
		for(j = 0;j < 6; j++) {
			_Save(ar,forces_all[index][i][j],dir_str[j]); 
		}

		EndRecord();
	}

	ar.table_index = TABLE::MDisplacements;
	for(i = 0;i < nDiv; i++) {
		StartRecord();

		AddHeader(ar,pcase,this);
		_Save(ar,station[i],"Station");
		for(j = 0;j < 6; j++) {
			_Save(ar,disps_all[index][i][j],dir_str[j]);
		}

		EndRecord();
	}


	if(pDoc->DesignResult) {
		ar.table_index = TABLE::MDesign;
		for(i = 0;i < nDiv; i++) {
			StartRecord();
			
			AddHeader(ar,pcase,this);
			_Save(ar,station[i],"Station");
			design_data[i].Serialize(ar);
			
			EndRecord();
		}
	}
}
void SLAB::SerializeAnalysis(CSerializer& ar,int index,ANALYSISCASE* pcase) {
	UINT i,j;
	if(!ar.is_storing)
		AllocAnalysis(index);

	ar.table_index = TABLE::SForces;
	for(j = 0;j < NJ;j++) {
		StartRecord();

		AddHeader(ar,pcase,this);
		_Save(ar,j,"Joint");
		for(i = 0;i < NSTRESS; i++) {
			_Save(ar,forces_all[index][j * NSTRESS + i],stress_str[i]);
		}

		EndRecord();
	}

	ar.table_index = TABLE::SDisplacements;
	for(j = 0;j < NJ;j++) {
		StartRecord();

		AddHeader(ar,pcase,this);
		_Save(ar,j,"Joint");
		for(i = 0;i < 6; i++) {
			_Save(ar,disps_all[index][j * 6 + i],dir_str[i]);
		}

		EndRecord();
	}
}
void CmyDocument::SerializeAnalysisResult(CSerializer& ar) { 
	if(!AnalysisResult)
		return;
	
	POSITION pos,pos1;
	PANALYSISCASELIST pcaselist;
	ANALYSISCASE* pcase;
	UINT i,j,ientity,repeats;
	ENTITY* pentity;
	PENTITYLIST pentities;

	for(ientity = 0;ientity < NENTITIES;ientity++) {
		pentities = entities[ientity];
		
		pos = pentities->GetHeadPosition();
		while(pos) {
			pentity = &pentities->GetNext(pos);
			
			if(!ar.is_storing)
				pentity->AllocVectors();
			
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
							pcase->selected = i;
							pentity->SerializeAnalysis(ar,pcase->index + i,pcase);
						}
					}
				}
			}
		}
	}
	
	/*
	eigen values
	*/
#define SaveModes(data,MYDATATYPE) {                           \
	UINT i;                                                    \
	MYDATATYPE* pcase;                                         \
	pos = (data).GetHeadPosition();                            \
	while(pos) {                                               \
		pcase = &(data).GetNext(pos);                          \
		if(pcase->finished) {                                  \
		    if(!ar.is_storing)  {                              \
				vec_alloc(pcase->eigvalue,pcase->runmodes);    \
			}                                                  \
			for(i = 0;i < pcase->runmodes;i++) {               \
			    StartRecord();                                 \
			    _Save(ar,pcase->name,"Name");                  \
				_Save(ar,i,"Mode");                            \
				_Save(ar,pcase->eigvalue[i],"EigVal");         \
				EndRecord();                                   \
			}                                                  \
		}                                                      \
	}                                                          \
};

	ar.table_index = TABLE::MModal;
	SaveModes(modalcases,MODALCASE);

	ar.table_index = TABLE::MBuckling;
    SaveModes(bucklingcases,BUCKLINGCASE);
}
/*
Serialize view
*/
void CmyView::Serialize(CSerializer& ar) {
	mydisplay.Serialize(ar);
	
	int myview = view;
	Save(ar,myview); 
	view = VIEW(myview);
	
	Save(ar,position); 
	Save(ar,pvalue); 
	Save(ar,Scale); 
	Save(ar,SELSIZE); 
	Save(ar,force_diagram); 
	Save(ar,average_stress);
	
    if(!ar.is_storing)
        InitView();
}