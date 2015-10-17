#include "common.h"
#include "GUI.h"

#define FindJoint(x)   (JOINT*)FindEntity((x),0)
#define FindMember(x)  (MEMBER*)FindEntity((x),1)
#define FindSlab(x)    (SLAB*)FindEntity((x),2)


static const CString DbNames[] = {"Excel","Access","AutoCad","0"};
static const CString DbExtension[] = {"xls","mdb","dxf","0"};
static const CString sec_type_str[] = {"RECTANGULAR", "CIRCULAR", "WIDEFLANGE","CHANNEL","DOUBLE_CHANNEL","TEE","ANGLE","DOUBLE_ANGLE","BOX","PIPE","GENERAL","0"};
static const CString mat_type_str[] = {"CONCRETE", "STEEL","FILLER","0"};
static const CString design_str[] = {"COULMN", "BEAM", "SLAB","0"};
static const CString rank_str[] = {"LEVEL1", "LEVEL2", "LEVEL3","LEVEL4","0"};
static const CString loadtype_str[] = {"CONCENTRATED","UNIFORM","TRAPEZOIDAL","FORCE","DISPLACEMENT","MASS","TEMPERATURE","STRAIN","0"};
static const CString dir_str[] = {"UX","UY","UZ","RX","RY","RZ","0"};
static const CString design_case_str[] = {"DEAD", "LIVE","WIND","QUAKE","0"};
static const CString analysis_case_str[] = {"STATIC","MODAL","RESPONSE HISTORY","RESPONSES SPECTRUM","COMBINATION","NONLINEAR_STATIC","BUCKLING","0"};
static const CString step_type_str[] = {"","MODAL","TIME","","","0"};
static const CString history_str[] = {"PERIODIC","TRANSIENT","0"};
static const CString history_ana_str[] = {"MODAL","DIRECT_INTEGRATION","0"};
static const CString comb_str[] = { "ABSSUM", "SRSS", "CQC","0"};
static const CString comb_type_str[] = {"LINEAR", "ENVELOPE", "ABSSUM", "SRSS","0"};
static const CString buckle_type_str[] = {"LINEAR","0"};
static const CString nlcase_type_str[] = {"PDELA_linearized","PDELA_consistent","0"};
static const CString func_str[] = {"EBCS_DESIGN","EBCS_ELASTIC","USER","SINE","COSINE","TRIANGULAR","PERIODIC","0"};
static const CString soil_str[] = {"SOILA", "SOILB", "SOILC","0"};
static const CString damp_str[] = {"RAYLEIGH", "CONSTANT","NONE","0"};
static const CString cons_str[] = {"BODY","DIAPHRAGM","EQUAL","PLATE","BEAM","ROD","LINE","LOCAL","WELD","0"};
static const CString system_str[] = {"CARTESIAN","RADIAL","0"};

/*
my Database
*/
void CmyData::GetFieldsInString(CString str,CString& sFields,int& nCount,int start,BOOL skip) {
	BOOL copy = FALSE,wasWhite = FALSE;
	int i,c,len;
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
	int i,c,nCount = 0,len = str.GetLength();
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
static int FindString(CDBVariant& varValue,CString* list) {
	int i = 0;
	while(true) {
		if(*varValue.m_pstring == list[i]) {
			return i;
		}
		if(*varValue.m_pstring == "0")
		    break;
		i++;
	}
	return INVALID;
}

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
static void CreateTable(CDatabase& database,CmyDatabase* mydatabase, CString str) {
	
	if(mydatabase)
		mydatabase->Create(str);
	else {
		TRY {
			CString sSql = "CREATE TABLE " + str;
			database.ExecuteSQL(sSql);
			
		}
		CATCH(CDBException, e)
		{
			CString sSql;
			int count = str.Find("\"",2);
			for(int i = 0;i <= count;i++)
				sSql += str.GetAt(i);
			sSql= "DROP TABLE " + sSql;
			database.ExecuteSQL(sSql);
			sSql = "CREATE TABLE " + str;
			database.ExecuteSQL(sSql);
		}
		END_CATCH;
	}
}

static void InsertInTable(CDatabase& database,CmyDatabase* mydatabase, CString str) {
   	if(mydatabase)
		mydatabase->Insert(str);
	else {
		CString sSql = "INSERT INTO " + str;
		database.ExecuteSQL(sSql);
	}
}

void CmyDocument::OnExport(UINT nID) {
	Export(nID);
}
BOOL CmyDocument::Export(UINT nID,CmyDatabase* mydatabase) {
	int type,i;
	UINT maxn;
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

	/*select tabels*/
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

			POSITION pos,pos1;
			JOINT* joint;
			MEMBER* member;
			SLAB* slab;
			MATERIAL* pmat;
			SECTION* psec;
			FUNCTION* pfunc;
			UBMP8 res;
			RPoint v;
			BOOL first;
			
			/*
			General Settings
			*/
#define InsertText(x,y) {\
	sSql.Format("\"General Settings\" (item,tvalue) VALUES ('%s','%s')",x,y);\
	InsertInTable(database,mydatabase,sSql);\
			};
#define InsertInt(x,y) {\
	sSql.Format("\"General Settings\" (item,nvalue) VALUES ('%s',%d)",x,y);\
	InsertInTable(database,mydatabase,sSql);\
			};
#define InsertDouble(x,y) {\
	sSql.Format("\"General Settings\" (item,nvalue) VALUES ('%s',%f)",x,y);\
	InsertInTable(database,mydatabase,sSql);\
			};
			if(Dia.state[Dia.GSettings]) {
				sSql = "\"General Settings\" (item TEXT,nvalue NUMBER,tvalue TEXT)";
				CreateTable(database,mydatabase,sSql);
				InsertInt("Number_of_Views",CmyView::nViews);
				InsertText("Analysis_Result",AnalysisResult ? "YES" : "NO");
				InsertText("Design_Result",DesignResult ? "YES" : "NO");
				InsertText("Locked",Lock ? "YES" : "NO");
				InsertDouble("Force_Scale",user_force_scale);
				InsertDouble("Displacement_scale",user_disp_scale);
				InsertDouble("Load_Scale",user_load_scale);
				InsertInt("Number_of_joints",JOINT::TotalNumber);
				InsertInt("Number_of_members",MEMBER::TotalNumber);
				InsertInt("Number_of_slabs",SLAB::TotalNumber);
				InsertInt("Number_of_loadcases",JOINT::TotalLoadCases);
				InsertInt("Static_Dofs",static_dofs);
				InsertInt("Dynamic_Dofs",dynamic_dofs); 
			}
			
#undef InsertText
#undef InsertInt
#undef InsertDouble
			
			/*Systems*/
			if(Dia.state[Dia.GCoordinate_systems]) {
				SYSTEM* psys;
				sSql = "\"Coordinate Systems\" (name TEXT,type TEXT,X NUMBER,Y NUMBER,Z NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
				CreateTable(database,mydatabase,sSql);
				sSql = "\"Grid lines\" (name TEXT,dir TEXT,coordinate NUMBER)";
				CreateTable(database,mydatabase,sSql);

				pos = systems.GetHeadPosition();
				while(pos) {
					psys = &systems.GetNext(pos);
					sSql.Format("\"Coordinate Systems\" (name,type,X,Y,Z,RX,RY,RZ)"
						"VALUES ('%s','%s',%f,%f,%f,%f,%f,%f)",
						psys->name,system_str[psys->coordinate],psys->origin.x,psys->origin.y,psys->origin.z,
						psys->rotation.x,psys->rotation.y,psys->rotation.z);
					InsertInTable(database,mydatabase,sSql);

					DOUBLE v;
					for(i = 0;i < 3;i++) {
					    pos1 = psys->grid[i].GetHeadPosition();
						while(pos1) {
							v = psys->grid[i].GetNext(pos1);
							sSql.Format("\"Grid lines\" (name,dir,coordinate)"
								"VALUES ('%s','%s',%f)",psys->name,dir_str[i],v);
							InsertInTable(database,mydatabase,sSql);
						}
					}
				}
			}

			/*Joints*/
			JLOAD* pjload;
			
			first = TRUE;
			pos = joints.GetHeadPosition();
			while(pos) {
				joint = &joints.GetNext(pos);
				
				if(first) {
					if(Dia.state[Dia.JCoordinates]) {
						sSql = "\"Joint Coordinates\" (Joint TEXT, X NUMBER, Y NUMBER, Z NUMBER,nmembers NUMBER,nslabs NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.JRestraints]) {
						sSql = "\"Joint Restraint Assignments\" (Joint TEXT, UX TEXT,UY TEXT,UZ TEXT,RX TEXT,RY TEXT,RZ TEXT)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.JLoads]) {
						sSql = "\"Joint Loads\" (Joint TEXT,Case TEXT,dirtype TEXT,loadtype TEXT,UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}

					if(Dia.state[Dia.JAddedMass]) {
						sSql = "\"Joint Added Mass\" (Joint TEXT,UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.JLocal_axis]) {
						sSql = "\"Joint Local Axis Assignments\" (Joint TEXT, qx NUMBER, qy NUMBER, qz NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
				}
				
				/*Joint coordinate table*/
				if(Dia.state[Dia.JCoordinates]) {
					sSql.Format("\"Joint Coordinates\" (Joint,X,Y,Z,nmembers,nslabs) VALUES ('%s',%.2f,%.2f,%.2f,%d,%d)",
						joint->name,joint->p.x,joint->p.y,joint->p.z,joint->mconnect,joint->sconnect);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Joint restraint table*/
				if(Dia.state[Dia.JRestraints]) {
					if(joint->restraint) {
#define TX(x) ((res & x) ? "YES" : "NO")
						res = joint->restraint;
						sSql.Format("\"Joint Restraint Assignments\" (Joint,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s','%s','%s','%s')",
							joint->name,TX(UX),TX(UY),TX(UZ),TX(RX),TX(RY),TX(RZ));
						InsertInTable(database,mydatabase,sSql);
#undef TX
					}
				}
                /*Joint load table*/
				if(Dia.state[Dia.JLoads]) {
					pos1 = joint->load.GetHeadPosition();
					while(pos1) {
						pjload = &joint->load.GetNext(pos1);
						if(pjload->type == MASS) continue;
						sSql.Format("\"Joint Loads\" (Joint,Case,dirtype,loadtype,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s',%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)",
							joint->name,pjload->loadcase->name,pjload->system ? pjload->system->name : "LOCAL",loadtype_str[pjload->type],
							pjload->Q[0],pjload->Q[1],pjload->Q[2],pjload->Q[3],pjload->Q[4],pjload->Q[5]);
						InsertInTable(database,mydatabase,sSql);
					}
				}
				/*Joint added mass table*/
				if(Dia.state[Dia.JAddedMass]) {
					pos1 = joint->load.GetHeadPosition();
					while(pos1) {
						pjload = &joint->load.GetNext(pos1);
						if(pjload->type != MASS) continue;
						sSql.Format("\"Joint Added Mass\" (Joint,UX,UY,UZ,RX,RY,RZ) VALUES ('%s',%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)",
							joint->name,pjload->Q[0],pjload->Q[1],pjload->Q[2],pjload->Q[3],pjload->Q[4],pjload->Q[5]);
						InsertInTable(database,mydatabase,sSql);
					}
				}
				/*Joint Local Axis*/
				if(Dia.state[Dia.JLocal_axis]) {
					RPoint rotation = joint->rotation * (180 / PI);
					sSql.Format("\"Joint Local Axis Assignments\" (Joint,qx,qy,qz) VALUES ('%s',%f,%f,%f)",
						joint->name,rotation.x,rotation.y,rotation.z);
					InsertInTable(database,mydatabase,sSql);
				}
				
				first = FALSE;
			}
			/*Member*/
			LOAD* pload;
			first = TRUE;
			pos = members.GetHeadPosition();
			while(pos) {
				member = &members.GetNext(pos);
				if(first) {
					if(Dia.state[Dia.MConnectivity]) {
						sSql = "\"Member Connectivity\" (Member TEXT, JointI NUMBER, JointJ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MSection]) {
						sSql = "\"Member Section Assignments\" (Member TEXT, sSection TEXT, eSection TEXT,sFactor NUMBER,eFactor NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MLoads]) {
						sSql = "\"Member Loads\" (Member TEXT,Case TEXT,dirtype TEXT,loadtype TEXT,dir TEXT,X NUMBER,P NUMBER,X1 NUMBER,P1 NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MLocal_axis]) {
						sSql = "\"Member Local Axis Assignments\" (Member TEXT, Angle NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MRelease]) {
						sSql = "\"Member Release Assignments\" (Member TEXT,nUX TEXT,nUY TEXT,nUZ TEXT,nRX TEXT,nRY TEXT,nRZ TEXT,fUX TEXT,fUY TEXT,fUZ TEXT,fRX TEXT,fRY TEXT,fRZ TEXT)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MDivisions]) {
						sSql = "\"Member Divisions\" (Member TEXT, nMinDiv NUMBER, nDiv NUMBER, nOutputDiv NUMBER ,DivAtInterim TEXT)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.MStrain_loads]) {
						sSql = "\"Member Strain Loads\" (Member TEXT, temperature NUMBER, faberror NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
				}
				
				
				/*Member connectivity table*/
				if(Dia.state[Dia.MConnectivity]) {
					sSql.Format("\"Member Connectivity\" (Member,JointI,JointJ) VALUES ('%s','%s','%s')",
						member->name,member->j1->name,member->j2->name);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Member Section Assignments table*/
				if(Dia.state[Dia.MSection]) {
					sSql.Format("\"Member Section Assignments\" (Member,sSection,eSection,sFactor,eFactor) VALUES ('%s','%s','%s',%f,%f)",
						member->name,member->section->name,member->e_section ? member->e_section->name:"",member->start_offset,member->end_offset);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Member Load Assignments*/
				if(Dia.state[Dia.MLoads]) {
					pos1 = member->load.GetHeadPosition();
					while(pos1) {
						pload = &member->load.GetNext(pos1);
												
						sSql.Format("\"Member Loads\" (Member,Case,dirtype,loadtype,dir,X,P,X1,P1) VALUES ('%s','%s','%s','%s','%s',%.2f,%.2f,%.2f,%.2f)",
							member->name,pload->loadcase->name,pload->system ? pload->system->name : "LOCAL",loadtype_str[pload->type],dir_str[pload->dir],
							pload->x,pload->P,pload->x1,pload->P1);
						InsertInTable(database,mydatabase,sSql);
					}
				}
				/*Member Local Axis*/
				if(Dia.state[Dia.MLocal_axis]) {
					sSql.Format("\"Member Local Axis Assignments\" (Member,Angle) VALUES ('%s',%f)",
						member->name,(member->alpha * 180) / PI);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Member Releases*/
				if(Dia.state[Dia.MRelease]) {
					if(member->nrelease || member->frelease) {
#define TN(x) ((nrel & x) ? "YES" : "NO")
#define TF(x) ((frel & x) ? "YES" : "NO")
						UBMP8 nrel = member->nrelease,frel = member->frelease;
						sSql.Format("\"Member Release Assignments\" (Member,nUX,nUY,nUZ,nRX,nRY,nRZ,fUX,fUY,fUZ,fRX,fRY,fRZ) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')",
							member->name,TN(UX),TN(UY),TN(UZ),TN(RX),TN(RY),TN(RZ),TF(UX),TF(UY),TF(UZ),TF(RX),TF(RY),TF(RZ));
						InsertInTable(database,mydatabase,sSql);
#undef TF
#undef TN
					}
				}
				/*Member Divisions*/
				if(Dia.state[Dia.MDivisions]) {
					sSql.Format("\"Member Divisions\" (Member,nMinDiv,nDiv,nOutputDiv,DivAtInterim) VALUES ('%s',%d,%d,%d,'%s')",
						member->name,member->nMinFrameDiv,member->nDiv,member->nMinDiv,member->DivAtInterim ? "YES" : "NO");
					InsertInTable(database,mydatabase,sSql);
				}
				
				first = FALSE;
			}
			
			/*Slabs*/
			first = TRUE;
			pos = slabs.GetHeadPosition();
			while(pos) {
				slab = &slabs.GetNext(pos);
				if(first) {
					if(Dia.state[Dia.SConnectivity]) {
						sSql = "\"Slab Connectivity\" (Slab TEXT, Joint1 NUMBER, Joint2 NUMBER,Joint3 NUMBER,Joint4 NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.SSection]) {
						sSql = "\"Slab Section Assignments\" (Slab TEXT, Section TEXT)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.SLoads]) {
						sSql = "\"Slab Loads\" (Slab TEXT,Case TEXT,dirtype TEXT,loadtype TEXT,dir TEXT,P NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					if(Dia.state[Dia.SLocal_axis]) {
						sSql = "\"Slab Local Axis Assignments\" (Slab TEXT, Angle NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
				}
				
				/*Slab connectivity table*/
				if(Dia.state[Dia.SConnectivity]) {
					sSql.Format("\"Slab Connectivity\" (Slab,Joint1,Joint2,Joint3,Joint4) VALUES ('%s','%s','%s','%s','%s')",
						slab->name,slab->jt[0]->name,slab->jt[1]->name,slab->jt[2]->name,slab->jt[3]->name);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Slab Section Assignments table*/
				if(Dia.state[Dia.SSection]) {
					sSql.Format("\"Slab Section Assignments\" (Slab,Section) VALUES ('%s','%s')",
						slab->name,slab->section->name);
					InsertInTable(database,mydatabase,sSql);
				}
				
				/*Slab Load Assignments*/
				if(Dia.state[Dia.SLoads]) {
					pos1 = slab->load.GetHeadPosition();
					while(pos1) {
						pload = &slab->load.GetNext(pos1);
						
						sSql.Format("\"Slab Loads\" (Slab,Case,dirtype,loadtype,dir,P) VALUES ('%s','%s','%s','%s','%s',%.2f)",
							slab->name,pload->loadcase->name,pload->system ? pload->system->name : "LOCAL",loadtype_str[pload->type],dir_str[pload->dir],
							pload->P);
						InsertInTable(database,mydatabase,sSql);
					}
				}
				
				/*Slab Local Axis*/
				if(Dia.state[Dia.SLocal_axis]) {
					sSql.Format("\"Slab Local Axis Assignments\" (Slab,Angle) VALUES ('%s',%f)",
						slab->name,(slab->alpha * 180) / PI);
					InsertInTable(database,mydatabase,sSql);
				}
				
				first = FALSE;
			}
			/*Constraints*/
			if(Dia.state[Dia.JConstraints]) {
#define TX(x) ((con & x) ? "YES" : "NO")
				CONSTRAINT* pcons;
				JOINT* pjoint;
				UBMP8 con;
				sSql = "\"Joint Constraints\" (name TEXT,joints NUMBER,type TEXT,rank TEXT,UX TEXT,UY TEXT,UZ TEXT,RX TEXT,RY TEXT,RZ TEXT)";
				CreateTable(database,mydatabase,sSql);
				pos = constraints.GetHeadPosition();
				while(pos) {
					pcons = &constraints.GetNext(pos);
					
					con = pcons->constraint;
					sSql.Format("\"Joint Constraints\" (name,type,rank,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s')",
						pcons->name,cons_str[pcons->type],rank_str[pcons->rank],TX(UX),TX(UY),TX(UZ),TX(RX),TX(RY),TX(RZ));
					InsertInTable(database,mydatabase,sSql);
					
					pos1 = pcons->jointlist.GetHeadPosition();
					while(pos1) {
						pjoint = pcons->jointlist.GetNext(pos1);
						
						sSql.Format("\"Joint Constraints\" (name,joints) VALUES ('%s',%d)",pcons->name,pjoint->name);
						InsertInTable(database,mydatabase,sSql);
					}
				}
#undef TX
			}
			/*Material Properties table*/
			if(Dia.state[Dia.DMaterial]) {
				sSql = "\"Material Properties\" (name TEXT,type TEXT,E NUMBER,nu NUMBER,unitweight NUMBER,density NUMBER, alphac NUMBER, G NUMBER,rank TEXT,fck NUMBER,fctk NUMBER,fyk NUMBER,fyks NUMBER,ftk NUMBER,ftks NUMBER,Es NUMBER)";
				CreateTable(database,mydatabase,sSql);
				pos = materials.GetHeadPosition();
				while(pos) {
					pmat = &materials.GetNext(pos);
					sSql.Format("\"Material Properties\" (name,type,E,nu,unitweight,density,alphac,G,rank,fck,fctk,fyk,fyks,ftk,ftks,Es)"
						"VALUES ('%s','%s',%g,%g,%g,%g,%g,%g,'%s',%g,%g,%g,%g,%g,%g,%g)",
						pmat->name,mat_type_str[pmat->type],pmat->E,pmat->nu,pmat->unitweight,pmat->density,pmat->alphac,pmat->G,rank_str[pmat->rank],
						pmat->fck,pmat->fctk,pmat->fyk,pmat->fyks,pmat->ftk,pmat->ftks,pmat->Es);
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*Section Properties table*/
			if(Dia.state[Dia.DSection]) {
				sSql = "\"Section Properties\" (name TEXT,material TEXT,type TEXT,w NUMBER,h NUMBER,r NUMBER, tw NUMBER, tf NUMBER, bd NUMBER,A NUMBER,Ay NUMBER,Az NUMBER, Ix NUMBER, Iy NUMBER,Iz NUMBER,ry NUMBER,rz NUMBER,Zy NUMBER,Zz NUMBER,Sy NUMBER,Sz NUMBER,fA NUMBER,fAy NUMBER,fAz NUMBER,fIx NUMBER,fIy NUMBER,fIz NUMBER,rank TEXT,design TEXT,rtype TEXT,nz NUMBER,ny NUMBER,nt NUMBER,cover NUMBER,barsize NUMBER,sbarsize NUMBER)";
				CreateTable(database,mydatabase,sSql);
				pos = sections.GetHeadPosition();
				while(pos) {
					psec = &sections.GetNext(pos);
					sSql.Format("\"Section Properties\" (name,material,type,w,h,r,tw,tf,bd,A,Ay,Az,Ix,Iy,Iz,ry,rz,Zy,Zz,Sy,Sz,fA,fAy,fAz,fIx,fIy,fIz,rank,design,rtype,nz,ny,nt,cover,barsize,sbarsize)"
						"VALUES ('%s','%s','%s',%f,%f,%f,%f,%f,%f,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%f,%f,%f,%f,%f,%f,'%s','%s','%s',%d,%d,%d,%f,%f,%f)",
						psec->name,psec->material->name,sec_type_str[psec->type],psec->w,psec->h,psec->r,psec->tw,psec->tf,psec->bd,
						psec->A,psec->Ay,psec->Az,psec->Ix,psec->Iy,psec->Iz,psec->ry,psec->rz,psec->Zy,psec->Zz,psec->Sy,psec->Sz,
						psec->fA,psec->fAy,psec->fAz,psec->fIx,psec->fIy,psec->fIz,
						rank_str[psec->rank],design_str[psec->rebar.design],sec_type_str[psec->rebar.type],psec->rebar.nz,psec->rebar.ny,psec->rebar.nt,psec->rebar.cover,psec->rebar.barsize,psec->rebar.stirrup_barsize);
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*response history functions*/
            if(Dia.state[Dia.LR_history]) {
				sSql = "\"Response History Functions\" (name TEXT,ctime NUMBER,accel NUMBER,type TEXT,rank TEXT)";
				CreateTable(database,mydatabase,sSql);
				pos = rhfunctions.GetHeadPosition();
				while(pos) {
					pfunc = &rhfunctions.GetNext(pos);
					
					sSql.Format("\"Response History Functions\" (name,type,rank) VALUES ('%s','%s','%s')",
						pfunc->name,func_str[pfunc->type],rank_str[pfunc->rank]);
					InsertInTable(database,mydatabase,sSql);
					
					pos1 = pfunc->points.GetHeadPosition();
					while(pos1) {
						v = pfunc->points.GetNext(pos1);
						sSql.Format("\"Response History Functions\" (name,ctime,accel) VALUES ('%s',%f,%f)",
							pfunc->name,v.x,v.y);
						InsertInTable(database,mydatabase,sSql);
					}
				}
			}
			/*response spectrum functions*/
			if(Dia.state[Dia.LR_spectrum]) {
				sSql = "\"Response Spectrum Functions\" (name TEXT,period NUMBER,accel NUMBER,type TEXT,rank TEXT,Ag NUMBER,soiltype TEXT)";
				CreateTable(database,mydatabase,sSql);
				pos = rsfunctions.GetHeadPosition();
				while(pos) {
					pfunc = &rsfunctions.GetNext(pos);
					
					sSql.Format("\"Response Spectrum Functions\" (name,type,rank,Ag,soiltype) VALUES ('%s','%s','%s',%f,'%s')",
						pfunc->name,func_str[pfunc->type],rank_str[pfunc->rank],pfunc->Ag,soil_str[pfunc->soil_type]);
					InsertInTable(database,mydatabase,sSql);
					
					pos1 = pfunc->points.GetHeadPosition();
					while(pos1) {
						v = pfunc->points.GetNext(pos1);
						sSql.Format("\"Response Spectrum Functions\" (name,period,accel) VALUES ('%s',%f,%f)",
							pfunc->name,v.x,v.y);
						InsertInTable(database,mydatabase,sSql);
					}
				}
			}
			/*Analysis cases*/
			ANALYSISCASE* pcase;
			if(Dia.state[Dia.ADefinitions]) {
				PANALYSISCASELIST pcaselist;
				sSql = "\"Analysis Case Definitions\" (name TEXT,type TEXT,run TEXT,finished TEXT,cindex NUMBER)";
				CreateTable(database,mydatabase,sSql);
				for(i = 0;i < CASETYPES;i++) {
					pcaselist = all_analysis_cases[i];
					if(pcaselist) {
						pos1 = pcaselist->GetHeadPosition();
						while(pos1) {
							pcase = &pcaselist->GetNext(pos1);
							sSql.Format("\"Analysis Case Definitions\" (name,type,run,finished,cindex) VALUES ('%s','%s','%s','%s',%d)",
								pcase->name,analysis_case_str[pcase->acase_type],pcase->run ? "YES" : "NO",pcase->finished ? "YES" : "NO",pcase->index);
							InsertInTable(database,mydatabase,sSql);
						}
					}
				}
			}
			/*Load cases*/
			if(Dia.state[Dia.ALoad_cases]) {
				LOADCASE* ploadcase;
				sSql = "\"Load Case Definitions\" (name TEXT,type TEXT,selfwtmutl NUMBER)";
				CreateTable(database,mydatabase,sSql);
				pos = loadcases.GetHeadPosition();
				while(pos) {
					ploadcase = &loadcases.GetNext(pos);
					sSql.Format("\"Load Case Definitions\" (name,type,selfwtmutl) VALUES ('%s','%s',%f)",
						ploadcase->name,design_case_str[ploadcase->type],ploadcase->swm);
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*Non linear static cases*/
			if(Dia.state[Dia.ANonLinear_cases]) {
				NLCASE* pnlcase;
				sSql = "\"Non Linear Case Definitions\" (name TEXT,nSteps NUMBER,nIteration NUMBER,tolerance NUMBER,savestiffness TEXT)";
				CreateTable(database,mydatabase,sSql);
				
				pos = nlcases.GetHeadPosition();
				while(pos) {
					pnlcase = &nlcases.GetNext(pos);
					sSql.Format("\"Non Linear Case Definitions\" (name,nSteps,nIteration,tolerance,savestiffness) VALUES ('%s',%d,%d,%g,'%s')",
						pnlcase->name,pnlcase->nSteps,pnlcase->nIteration,pnlcase->tolerance,pnlcase->save_stiffness ? "YES":"NO");
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*Modal cases*/
			if(Dia.state[Dia.AModal_cases]) {
				MODALCASE* pmodalcase;
				sSql = "\"Modal Case Definitions\" (name TEXT,minmodes NUMBER,maxmodes NUMBER,runmodes NUMBER,shift NUMBER,tolerance NUMBER)";
				CreateTable(database,mydatabase,sSql);
				
				pos = modalcases.GetHeadPosition();
				while(pos) {
					pmodalcase = &modalcases.GetNext(pos);
					sSql.Format("\"Modal Case Definitions\" (name,minmodes,maxmodes,runmodes,shift,tolerance) VALUES ('%s',%d,%d,%d,%g,%g)",
						pmodalcase->name,pmodalcase->minm,pmodalcase->maxm,pmodalcase->runmodes,pmodalcase->shift,pmodalcase->tolerance);
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*Buckling cases*/
			if(Dia.state[Dia.ABuckling_cases]) {
				BUCKLINGCASE* pbklcase;
				sSql = "\"Buckling Case Definitions\" (name TEXT,nmodes NUMBER,runmodes NUMBER,tolerance NUMBER)";
				CreateTable(database,mydatabase,sSql);
				
				pos = bucklingcases.GetHeadPosition();
				while(pos) {
					pbklcase = &bucklingcases.GetNext(pos);
					sSql.Format("\"Buckling Case Definitions\" (name,nmodes,runmodes,tolerance) VALUES ('%s',%d,%d,%g)",
						pbklcase->name,pbklcase->nmodes,pbklcase->runmodes,pbklcase->tolerance);
					InsertInTable(database,mydatabase,sSql);
				}
			}
			/*Response history*/
			SPECFUNC* psfunc;
			DAMPING* pdamp;
			
			if(Dia.state[Dia.ARH_cases]) {
				RESPONSEHIST* prhcase;
				first = TRUE;
				sSql = "\"Response History Case Definitions\" (name TEXT,modes TEXT,type TEXT,analysistype TEXT,steps NUMBER,stepsize NUMBER)";
				CreateTable(database,mydatabase,sSql);
				pos = responsecases.GetHeadPosition();
				while(pos) {
					prhcase = &responsecases.GetNext(pos);
					
					sSql.Format("\"Response History Case Definitions\" (name,modes,type,analysistype,steps,stepsize) VALUES ('%s','%s','%s','%s',%d,%f)",
						prhcase->name,prhcase->modalcase->name,history_str[prhcase->type],history_ana_str[prhcase->ana_type],prhcase->N,prhcase->dt);
					InsertInTable(database,mydatabase,sSql);
					
					/*loads*/
					if(first) {
						sSql = "\"Response History Loading\" (name TEXT,loading TEXT,dir TEXT,scale NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					pos1 = prhcase->funclist.GetHeadPosition();
					while(pos1) {
						psfunc = &prhcase->funclist.GetNext(pos1);
						sSql.Format("\"Response History Loading\" (name,loading,dir,scale) VALUES ('%s','%s','%s',%f)",
							prhcase->name,psfunc->function->name,dir_str[psfunc->dir],psfunc->scale);
						InsertInTable(database,mydatabase,sSql);
					}
					
					/*damping*/
					if(first) {
						sSql = "\"Response History Damping\" (name TEXT,frequency NUMBER,damping NUMBER,type TEXT,w1 NUMBER,w2 NUMBER,e1 NUMBER,e2 NUMBER,fm NUMBER,fk NUMBER,cdamping NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					pdamp = &prhcase->damping;
					sSql.Format("\"Response History Damping\" (name,type,w1,w2,e1,e2,fm,fk,cdamping) VALUES ('%s','%s',%f,%f,%f,%f,%f,%f,%f)",
						prhcase->name,damp_str[pdamp->type],pdamp->w1,pdamp->w2,pdamp->e1,pdamp->e2,pdamp->fm,pdamp->fk,pdamp->cdamping);
					InsertInTable(database,mydatabase,sSql);
					
					pos1 = pdamp->dvalues.points.GetHeadPosition();
					while(pos1) {
						v = pdamp->dvalues.points.GetNext(pos1);
						sSql.Format("\"Response History Damping\" (name,frequency,damping) VALUES ('%s',%f,%f)",
							prhcase->name,v.x,v.y);
						InsertInTable(database,mydatabase,sSql);
					}
					
					first = FALSE;
				}
			}
			/*Response spectrum*/
			if(Dia.state[Dia.ARS_cases]) {
				RESPONSESPEC* prscase;
				
				first = TRUE;
				sSql = "\"Response Spectrum Case Definitions\" (name TEXT,modes TEXT,modalcomb TEXT,dircomb TEXT)";
				CreateTable(database,mydatabase,sSql);
				pos = responsespecs.GetHeadPosition();
				while(pos) {
					prscase = &responsespecs.GetNext(pos);
					
					sSql.Format("\"Response Spectrum Case Definitions\" (name,modes,modalcomb,dircomb) VALUES ('%s','%s','%s','%s')",
						prscase->name,prscase->modalcase->name,comb_str[prscase->modal_comb],comb_str[prscase->dir_comb]);
					InsertInTable(database,mydatabase,sSql);
					
					/*loading*/
					if(first) {
						sSql = "\"Response Spectrum Loading\" (name TEXT,loading TEXT,dir TEXT,scale NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					pos1 = prscase->funclist.GetHeadPosition();
					while(pos1) {
						psfunc = &prscase->funclist.GetNext(pos1);
						sSql.Format("\"Response Spectrum Loading\" (name,loading,dir,scale) VALUES ('%s','%s','%s',%f)",
							prscase->name,psfunc->function->name,dir_str[psfunc->dir],psfunc->scale);
						InsertInTable(database,mydatabase,sSql);
					}
					
					/*damping*/
					if(first) {
						sSql = "\"Response Spectrum Damping\" (name TEXT,frequency NUMBER,damping NUMBER,type TEXT,w1 NUMBER,w2 NUMBER,e1 NUMBER,e2 NUMBER,fm NUMBER,fk NUMBER,cdamping NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					pdamp = &prscase->damping;
					sSql.Format("\"Response Spectrum Damping\" (name,type,w1,w2,e1,e2,fm,fk,cdamping) VALUES ('%s','%s',%f,%f,%f,%f,%f,%f,%f)",
						prscase->name,damp_str[pdamp->type],pdamp->w1,pdamp->w2,pdamp->e1,pdamp->e2,pdamp->fm,pdamp->fk,pdamp->cdamping);
					InsertInTable(database,mydatabase,sSql);
					
					pos1 = pdamp->dvalues.points.GetHeadPosition();
					while(pos1) {
						v = pdamp->dvalues.points.GetNext(pos1);
						sSql.Format("\"Response Spectrum Damping\" (name,frequency,damping) VALUES ('%s',%f,%f)",
							prscase->name,v.x,v.y);
						InsertInTable(database,mydatabase,sSql);
					}
					
					first = FALSE;
				}
			}
			/*Combination*/
			if(Dia.state[Dia.ACombinations]) {
				COMB_TYPE* pcombo;
				LOADCOMBO* pldcombo;
				sSql = "\"Combination Definitions\" (name TEXT,type TEXT,casetype TEXT,casename TEXT,FS NUMBER)";
				CreateTable(database,mydatabase,sSql);

				CList<COMB_TYPE,COMB_TYPE&>* pcombolist;
				const CString* type_str;
				for(int i = 0;i < 3;i++) {
					if(i == 0) {
						pcombolist = (CList<COMB_TYPE,COMB_TYPE&>*) &combinations;
						type_str = comb_type_str;
					} else if(i == 1) {
						pcombolist = (CList<COMB_TYPE,COMB_TYPE&>*) &bucklingcases;
						type_str = buckle_type_str;
					} else if(i == 2) {
						pcombolist = (CList<COMB_TYPE,COMB_TYPE&>*) &nlcases;
						type_str = nlcase_type_str;
					}

					pos = pcombolist->GetHeadPosition();
					while(pos) {
						pcombo = &pcombolist->GetNext(pos);
						
						pos1 = pcombo->loadlist.GetHeadPosition();
						while(pos1) {
							pldcombo = &pcombo->loadlist.GetNext(pos1);
							pcase = pldcombo->loadcase;
							
							sSql.Format("\"Combination Definitions\" (name,type,casetype,casename,FS) VALUES ('%s','%s','%s','%s',%f)",
								pcombo->name,type_str[pcombo->type],analysis_case_str[pcase->acase_type],pcase->name,pldcombo->FS);
							InsertInTable(database,mydatabase,sSql);
						}
					}
				}
			}
			/*
			Analysis result
			*/
			if(AnalysisResult && Dia.state[Dia.Analysis]) {
				int index;
				PANALYSISCASELIST pcaselist;
				ANALYSISCASE* pcase;
				/*
				Modal Eigen Values
				*/
				if(Dia.state[Dia.GModes]) {
					DOUBLE p;
					MODALCASE* pmodalcase;
					sSql = "\"Modal Eigen Values\" (name TEXT,mode NUMBER,eigvalue NUMBER,frequency NUMBER,period NUMBER)";
					CreateTable(database,mydatabase,sSql);
					pos = modalcases.GetHeadPosition();
					while(pos) {
						pmodalcase = &modalcases.GetNext(pos);
						if(!pmodalcase->finished) continue;
						for(UINT j = 0;j < pmodalcase->runmodes;j++) {
                            p = pmodalcase->eigvalue[j];
							sSql.Format("\"Modal Eigen Values\" (name,mode,eigvalue,frequency,period) VALUES ('%s',%d,%f,%f,%f)",
								pmodalcase->name,j,p,sqrt(p),2 * PI / sqrt(p));
							InsertInTable(database,mydatabase,sSql);
						}
					}
				}
				/*
				Buckling modes
				*/
				if(Dia.state[Dia.GBucklingModes]) {
					DOUBLE p;
					BUCKLINGCASE* pbklcase;
					sSql = "\"Buckling modes\" (name TEXT,mode NUMBER,factor NUMBER)";
					CreateTable(database,mydatabase,sSql);
					pos = bucklingcases.GetHeadPosition();
					while(pos) {
						pbklcase = &bucklingcases.GetNext(pos);
						if(!pbklcase->finished) continue;
						for(UINT j = 0;j < pbklcase->runmodes;j++) {
                            p = pbklcase->eigvalue[j];
							sSql.Format("\"Buckling modes\" (name,mode,factor) VALUES ('%s',%d,%f)",
								pbklcase->name,j,p);
							InsertInTable(database,mydatabase,sSql);
						}
					}
				}
				/*
				Joint
				*/
				if(Dia.state[Dia.AJoint]) {
					if(Dia.state[Dia.JReactions]) {
						sSql = "\"Joint Reactions\" (Joint TEXT, Case TEXT, CaseType TEXT,StepType TEXT,StepNum NUMBER,UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					if(Dia.state[Dia.JDisplacements]) {
						sSql = "\"Joint Displacements\" (Joint TEXT, Case TEXT, CaseType TEXT,StepType TEXT,StepNum NUMBER,UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					if(Dia.state[Dia.JMasses]) {
						sSql = "\"Joint Masses\" (Joint TEXT, UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
						pos = joints.GetHeadPosition();
						while(pos) {
							joint = &joints.GetNext(pos);
							sSql.Format("\"Joint Masses\" (Joint,UX,UY,UZ,RX,RY,RZ) VALUES ('%s',%f,%f,%f,%f,%f,%f)",
								joint->name,joint->massembled[IUX],joint->massembled[IUY],joint->massembled[IUZ],
								joint->massembled[IRX],joint->massembled[IRY],joint->massembled[IRZ]);
							InsertInTable(database,mydatabase,sSql);
						}
					}
					if(Dia.state[Dia.JReactions] || Dia.state[Dia.JDisplacements]) {
						pos = joints.GetHeadPosition();
						while(pos) {
							joint = &joints.GetNext(pos);

							for(i = 0;i < CASETYPES;i++) {
								pcaselist = all_analysis_cases[i];
								if(pcaselist) {
									pos1 = pcaselist->GetHeadPosition();
									while(pos1) {
										pcase = &pcaselist->GetNext(pos1);
										if(pcase->finished && pcase->selected) {
											
											if(pcase->acase_type == ANALYSISCASE::MODAL_CASE ||
                                                pcase->acase_type == ANALYSISCASE::BUCKLING_CASE ||
												pcase->acase_type == ANALYSISCASE::RESPONSEH_CASE
												) {
												if(pcase->acase_type == ANALYSISCASE::MODAL_CASE) maxn = ((MODALCASE*) pcase)->runmodes;
												else if(pcase->acase_type == ANALYSISCASE::BUCKLING_CASE) maxn = ((BUCKLINGCASE*) pcase)->runmodes;
												else maxn = ((RESPONSEHIST*) pcase)->N;
												for(UINT j = 0;j < maxn;j++) {
													index = pcase->index + j;
													/*reactions*/
													if(Dia.state[Dia.JReactions]) {
														if(joint->restraint) {
															sSql.Format("\"Joint Reactions\" (Joint,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
																joint->name,pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],j,
																joint->forces_all[index][0],joint->forces_all[index][1],joint->forces_all[index][2],joint->forces_all[index][3],joint->disps_all[index][4],joint->disps_all[index][5]);
															InsertInTable(database,mydatabase,sSql);
														}
													}
													/*displacement*/
													if(Dia.state[Dia.JDisplacements]) {
														sSql.Format("\"Joint Displacements\" (Joint,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
															joint->name,pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],j,
															joint->disps_all[index][0],joint->disps_all[index][1],joint->disps_all[index][2],joint->disps_all[index][3],joint->disps_all[index][4],joint->disps_all[index][5]);
														InsertInTable(database,mydatabase,sSql);
													}
												}
											} else {
												index = pcase->index;
												
												/*reactions*/
												if(Dia.state[Dia.JReactions]) {
													if(joint->restraint) {
														sSql.Format("\"Joint Reactions\" (Joint,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
															joint->name,pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],0,
															joint->forces_all[index][0],joint->forces_all[index][1],joint->forces_all[index][2],joint->forces_all[index][3],joint->disps_all[index][4],joint->disps_all[index][5]);
														InsertInTable(database,mydatabase,sSql);
													}
												}
												/*displacement*/
												if(Dia.state[Dia.JDisplacements]) {
													sSql.Format("\"Joint Displacements\" (Joint,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s','%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
														joint->name,pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],0,
														joint->disps_all[index][0],joint->disps_all[index][1],joint->disps_all[index][2],joint->disps_all[index][3],joint->disps_all[index][4],joint->disps_all[index][5]);
													InsertInTable(database,mydatabase,sSql);
												}
											}
										}
									}
								}
							}
						} 
					}
				}
				/*
				Member output
				*/
				if(Dia.state[Dia.AMember]) {
					if(Dia.state[Dia.MInternal_forces]) {
						sSql = "\"Member Internal Forces\" (Member TEXT, Station NUMBER,Case TEXT, CaseType TEXT,StepType TEXT,StepNum NUMBER,UX NUMBER,UY NUMBER,UZ NUMBER,RX NUMBER,RY NUMBER,RZ NUMBER)";
						CreateTable(database,mydatabase,sSql);
					}
					
					pos = members.GetHeadPosition();
					while(pos) {
						member = &members.GetNext(pos);
						
						for(i = 0;i < CASETYPES;i++) {
							pcaselist = all_analysis_cases[i];
							if(pcaselist) {
								pos1 = pcaselist->GetHeadPosition();
								while(pos1) {
									pcase = &pcaselist->GetNext(pos1);
									if(pcase->finished && pcase->selected) {
										
										if(pcase->acase_type == ANALYSISCASE::MODAL_CASE ||
											pcase->acase_type == ANALYSISCASE::BUCKLING_CASE ||
											pcase->acase_type == ANALYSISCASE::RESPONSEH_CASE
											) {
											if(pcase->acase_type == ANALYSISCASE::MODAL_CASE) maxn = ((MODALCASE*) pcase)->runmodes;
											else if(pcase->acase_type == ANALYSISCASE::BUCKLING_CASE) maxn = ((BUCKLINGCASE*) pcase)->runmodes;
											else maxn = ((RESPONSEHIST*) pcase)->N;
											
											for(UINT k = 0;k < maxn;k++) {
												index = pcase->index + k;
												if(Dia.state[Dia.MInternal_forces]) {
													for(UINT j = 0;j < member->nDiv;j++) {
														sSql.Format("\"Member Internal Forces\" (Member,Station,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s',%f,'%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
															member->name,member->station[j],pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],k,
															member->forces_all[index][j][0],member->forces_all[index][j][1],member->forces_all[index][j][2],member->forces_all[index][j][3],member->forces_all[index][j][4],member->forces_all[index][j][5]);
														InsertInTable(database,mydatabase,sSql);
													}
												}
											}
										} else {
											index = pcase->index;
											if(Dia.state[Dia.MInternal_forces]) {
												for(UINT j = 0;j < member->nDiv;j++) {
													sSql.Format("\"Member Internal Forces\" (Member,Station,Case,CaseType,StepType,StepNum,UX,UY,UZ,RX,RY,RZ) VALUES ('%s',%f,'%s','%s','%s',%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g)",
														member->name,member->station[j],pcase->name,analysis_case_str[pcase->acase_type],step_type_str[pcase->acase_type],0,
														member->forces_all[index][j][0],member->forces_all[index][j][1],member->forces_all[index][j][2],member->forces_all[index][j][3],member->forces_all[index][j][4],member->forces_all[index][j][5]);
													InsertInTable(database,mydatabase,sSql);
												}
											}
										}
									}
								}
							}
						}
					} 
				}
				/*
				end
				*/
			}
			/*
			Design Result
			*/
            if(DesignResult && Dia.state[Dia.Design]) {
				DESIGNDATA* pdata;
				
				if(Dia.state[Dia.DResult]) {
					sSql = "\"Concrete Design\" (Member TEXT, Station NUMBER,Status TEXT,TopArea NUMBER,BotArea NUMBER,ShearArea NUMBER,"
						"TopBarTypeCount NUMBER,TopCount1 NUMBER,TopBar1 NUMBER,TopCount2 NUMBER,TopBar2 NUMBER," 
						"BotBarTypeCount NUMBER,BotCount1 NUMBER,BotBar1 NUMBER,BotCount2 NUMBER,BotBar2 NUMBER,"
						"ShearSpacing NUMBER,ShearBar NUMBER" 
						")";
					CreateTable(database,mydatabase,sSql);

					sSql = "\"Steel Design\" (Member NUMBER,Status TEXT,"
						"RatioAxial NUMBER,RatioMinor NUMBER,RatioMajor NUMBER" 
						")";
					CreateTable(database,mydatabase,sSql);
					
					pos = members.GetHeadPosition();
					while(pos) {
						member = &members.GetNext(pos);
						if(member->section->material->type == MATERIAL::STEEL) {
							pdata = &member->design_data[0];
							sSql.Format("\"Steel Design\" (Member,Status,RatioAxial,RatioMinor,RatioMajor)"
								"VALUES ('%s','%s',%f,%f,%f)",
								member->name,pdata->flags ? "FAIL" : "OK",
								pdata->ratios[0],pdata->ratios[1],pdata->ratios[2]
								);
							InsertInTable(database,mydatabase,sSql);
						} else {
							for(UINT i = 0;i < member->nDiv;i++) {
								pdata = &member->design_data[i];
								sSql.Format("\"Concrete Design\" (Member,Station,Status,TopArea,BotArea,ShearArea,TopBarTypeCount,TopCount1,TopBar1,TopCount2,TopBar2,BotBarTypeCount,BotCount1,BotBar1,BotCount2,BotBar2,ShearSpacing,ShearBar)"
									"VALUES ('%s',%f,'%s',%f,%f,%f,%d,%d,%f,%d,%f,%d,%d,%f,%d,%f,%d,%f)",
									member->name,member->station[i],pdata->flags ? "FAIL" : "OK",pdata->Area[0],pdata->Area[1],pdata->As,
									pdata->total[0],pdata->count[0][0],pdata->diam[0][0],pdata->count[0][1],pdata->diam[0][1],
									pdata->total[1],pdata->count[1][0],pdata->diam[1][0],pdata->count[1][1],pdata->diam[1][1],
									pdata->s,pdata->sd
									);
								InsertInTable(database,mydatabase,sSql);
							}
						}
					} 
				}
			}
		}      
	} 
    CATCH(CDBException, e)
    {
        AfxMessageBox("Database error: " + e->m_strError);
		return FALSE;
    }
    END_CATCH;
	
	
    database.Close();
	return TRUE;
}
/*
Import
*/
static BOOL OpenRecord(CRecordset& rs,CString sSql) {
	TRY {
		rs.Open(CRecordset::forwardOnly,sSql);
		return TRUE;
	}
	CATCH(CDBException, e)
    {
		return FALSE;
    }
    END_CATCH;
}
void CmyDocument::OnImport(UINT nID) {
	int type,i;
	CDatabase database;
	CString sDriver,sSql,sPath,sFilter;
	
	if(nID == IDM_IMPORT_EXCEL) sFilter = "Excel files (*.xls)|*.xls|All Files (*.*)|*.*||";
	else if(nID == IDM_IMPORT_ACCESS) sFilter = "Access files (*.mdb)|*.mdb|All Files (*.*)|*.*||";
    else sFilter = "Autocad files (*.dxf)|*.dxf|All Files (*.*)|*.*||";
	
	CFileDialog FileDlg(TRUE,".xls",NULL,0,sFilter);
	
	if(FileDlg.DoModal() == IDOK) {
		sPath = FileDlg.GetPathName();
	} else {
		return;
	}
	
	type = nID - IDM_IMPORT_EXCEL;
	if(nID == IDM_IMPORT_AUTOCAD) {
		ImportFromAutocad(sPath);
		return;
	}
	
	sDriver = GetDriver(DbNames[type]);
	if (sDriver.IsEmpty()) {
		AfxMessageBox("Can't find odbc driver.");
		return;
	}
	
	TRY 
	{
		sSql.Format("DRIVER={%s};DSN='';FIRSTROWHASNAMES=1;READONLY=FALSE;DBQ=%s",sDriver, sPath);
		
		if( database.OpenEx(sSql,CDatabase::noOdbcDialog) ) {
			
			CRecordset rs( &database );
			CDBVariant varValue;
			short index;
			CString name;
			POSITION pos;
			RPoint v;
			
			/*
			General Settings
			*/
            if(OpenRecord(rs,_T( "SELECT * FROM \"General Settings\"" ) ) ) {
				index = 1;
				rs.GetFieldValue(index, varValue);
				CmyView::nViews = (int)varValue.m_dblVal;
				rs.MoveNext( );
				index = 2;
				rs.GetFieldValue(index, varValue);
				AnalysisResult = (*varValue.m_pstring == "YES");
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				DesignResult = (*varValue.m_pstring == "YES");
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				Lock = (*varValue.m_pstring == "YES");
				rs.MoveNext( );
				index = 1;
				rs.GetFieldValue(index, varValue);
				user_force_scale = varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				user_disp_scale = varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				user_load_scale = varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				JOINT::TotalNumber = (UINT)varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				MEMBER::TotalNumber = (UINT)varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				SLAB::TotalNumber = (UINT)varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				SLAB::TotalLoadCases = MEMBER::TotalLoadCases = JOINT::TotalLoadCases = (UINT)varValue.m_dblVal;
				rs.MoveNext( );
				rs.GetFieldValue(index, varValue);
				static_dofs = (UBMP8)varValue.m_dblVal;
				rs.MoveNext( ); 
				rs.GetFieldValue(index, varValue);
				dynamic_dofs = (UBMP8)varValue.m_dblVal;
				rs.MoveNext( ); 
				rs.Close();
			}
			/*Systems*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Coordinate Systems\"" ) ) ) {
				SYSTEM sys;
				BOOL first = TRUE;
				while( !rs.IsEOF( ) ) {
					index = 0;
					rs.GetFieldValue(index++, varValue);
					sys.name = *varValue.m_pstring;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)system_str);
					if(i != INVALID) sys.coordinate = i;
					rs.GetFieldValue(index++, varValue);
				    sys.origin.x = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
				    sys.origin.y = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
				    sys.origin.z = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
				    sys.rotation.x = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
				    sys.rotation.y = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
				    sys.rotation.z = varValue.m_dblVal;

					if(first) *global = sys;
					else systems.AddTail(sys);
					first = FALSE;					

					rs.MoveNext( );
				}
				rs.Close();
			}

			if(OpenRecord(rs,_T( "SELECT * FROM \"Grid lines\"" ) ) ) {
				SYSTEM sys,*psys;
				DOUBLE v;
				while( !rs.IsEOF( ) ) {
					index = 0;
					rs.GetFieldValue(index++, varValue);
					sys.name = *varValue.m_pstring;
					if(pos = systems.Find(sys)) {
						psys = &systems.GetAt(pos);
					}
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)dir_str);
					rs.GetFieldValue(index++, varValue);
				    v = varValue.m_dblVal;
					psys->grid[i].AddTail(v);

					rs.MoveNext( );
				}
				rs.Close();
			}

			/*materials*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Material Properties\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					MATERIAL mat,*pmat;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					mat.name = *varValue.m_pstring;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)mat_type_str);
					if(i != INVALID) mat.type = i;
					rs.GetFieldValue(index++, varValue);
					mat.E = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.nu = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.unitweight = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.density = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.alphac = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.G = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)rank_str);
					if(i != INVALID) mat.rank = i;
					
					rs.GetFieldValue(index++, varValue);
					mat.fck = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.fctk = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.fyk = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.fyks = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.ftk = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.ftks = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					mat.Es = varValue.m_dblVal;
					
					pos = materials.Find(mat);
					if(pos) {
						pmat = &materials.GetAt(pos);
						*pmat = mat;
					} else {
						materials.AddTail(mat);
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			/*sections*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Section Properties\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					SECTION sec,*psec;
					MATERIAL mat;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					sec.name = *varValue.m_pstring;
					
					rs.GetFieldValue(index++, varValue);
					mat.name = *varValue.m_pstring;
					pos = materials.Find(mat);
					sec.material = &materials.GetAt(pos);
					
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)sec_type_str);
					if(i != INVALID) sec.type = i;
					rs.GetFieldValue(index++, varValue);
					sec.w = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.h = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.r = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.tw = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.tf = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.bd = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.A = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Ay = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Az = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Ix = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Iy = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Iz = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.ry = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rz = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Zy = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Zz = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Sy = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.Sz = varValue.m_dblVal;
                    rs.GetFieldValue(index++, varValue);
					sec.fA = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.fAy = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.fAz = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.fIx = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.fIy = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.fIz = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)rank_str);
					if(i != INVALID) sec.rank = i;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)design_str);
					if(i != INVALID) sec.rebar.design = i;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)sec_type_str);
					if(i != INVALID) sec.rebar.type = i;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.nz = (UBMP8)varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.ny = (UBMP8)varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.nt = (UBMP8)varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.cover = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.barsize = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					sec.rebar.stirrup_barsize = varValue.m_dblVal;
					
					pos = sections.Find(sec);
					if(pos) {
						psec = &sections.GetAt(pos);
						*psec = sec;
					} else {
						sections.AddTail(sec);
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			if(JOINT::TotalNumber) {
				/*Joint coordinate table*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Coordinates\" ORDER BY Joint" ) ) ) {
					while( !rs.IsEOF( ) ) {
						JOINT joint;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						joint.name = *varValue.m_pstring;
						rs.GetFieldValue(index++, varValue);
						joint.p.x = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						joint.p.y = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						joint.p.z = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						joint.mconnect = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						joint.sconnect = (UINT)varValue.m_dblVal;
						
						joints.AddTail(joint);
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Joint restraint table*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Restraint Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						JOINT* joint = FindJoint(*varValue.m_pstring);
						if(joint) {
							for(i = IUX;i <= IRZ;i++) {
								rs.GetFieldValue(index++, varValue);
								if(*varValue.m_pstring == "YES") joint->restraint |= (1 << i);
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Joint Local Axis*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Local Axis Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						JOINT* joint = FindJoint(*varValue.m_pstring);
						if(joint) {
							RPoint rotation;
							rs.GetFieldValue(index++, varValue);
							rotation.x = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							rotation.y = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							rotation.z = varValue.m_dblVal;
							joint->rotation = rotation * (PI / 180);
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			if(MEMBER::TotalNumber) {
				/*Member connectivity table*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Connectivity\" ORDER BY Member" ) ) ) {
					while( !rs.IsEOF( ) ) {
						MEMBER member;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						member.name = *varValue.m_pstring;
						for(int i = 0;i < 2;i++) {
							rs.GetFieldValue(index++, varValue);
							if(i == 0) member.j1 = FindJoint(*varValue.m_pstring);
							else member.j2 = FindJoint(*varValue.m_pstring);
						}
						members.AddTail(member);
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Member section Assignment*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Section Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						
						MEMBER* member = FindMember(*varValue.m_pstring);
						if(member) {
							SECTION sec;
							for(i = 0;i < 2;i++) {
								rs.GetFieldValue(index++, varValue);
								sec.name = *varValue.m_pstring;
								pos = sections.Find(sec);
								if(pos) {
									if(i == 0) member->section = &sections.GetAt(pos);
									else member->e_section = &sections.GetAt(pos);
								}
							}
							rs.GetFieldValue(index++, varValue);
							member->start_offset = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							member->end_offset = varValue.m_dblVal;
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Member Local Axis*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Local Axis Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						MEMBER* member = FindMember(*varValue.m_pstring);
						if(member) {
							rs.GetFieldValue(index++, varValue);
							member->alpha = (PI * varValue.m_dblVal) / 180;
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Member Release Assignments*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Release Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						MEMBER* member = FindMember(*varValue.m_pstring);
						if(member) {
							for(i = 0;i < 12;i++) {
								rs.GetFieldValue(index++, varValue);
								if(*varValue.m_pstring == "YES") {
									if(i < 6) member->nrelease |= (1 << i);
									else member->frelease |= (1 << (i - 6));
								}
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
				/*Member Divisions*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Divisions\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						MEMBER* member = FindMember(*varValue.m_pstring);
						if(member) {
							rs.GetFieldValue(index++, varValue);
							member->nMinFrameDiv = (UINT) varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							member->nDiv = (UINT) varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							member->nMinDiv = (UINT) varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							if(*varValue.m_pstring == "YES") member->DivAtInterim = TRUE;
							else member->DivAtInterim = FALSE;
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
            } 
			if(SLAB::TotalNumber) {
				/*Slab connectivity table*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Slab Connectivity\" ORDER BY Slab" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SLAB slab;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						slab.name = *varValue.m_pstring;
						rs.GetFieldValue(index++, varValue);
						for(int i = 0;i < 4;i++) {
							rs.GetFieldValue(index++, varValue);
							slab.jt[i] = FindJoint(*varValue.m_pstring);
						}
						slabs.AddTail(slab);
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Slab section Assignment*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Slab Section Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						SLAB* slab = FindSlab(*varValue.m_pstring);
						if(slab) {
							SECTION sec;
							rs.GetFieldValue(index++, varValue);
							sec.name = *varValue.m_pstring;
							pos = sections.Find(sec);
							if(pos) {
								slab->section = &asections.GetAt(pos);
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Slab Local Axis*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Slab Local Axis Assignments\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						SLAB* slab = FindSlab(*varValue.m_pstring);
						if(slab) {
							rs.GetFieldValue(index++, varValue);
							slab->alpha = (PI * varValue.m_dblVal) / 180;
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			/*response history functions*/
			FUNCTION func,*pfunc = NULL;
			if(OpenRecord(rs,_T( "SELECT * FROM \"Response History Functions\" ORDER BY name,ctime" ) ) ) {
				while( !rs.IsEOF( ) ) {
					
					index = 0;
					rs.GetFieldValue(index++, varValue);
					func.name = *varValue.m_pstring;
					
					if(!pfunc || pfunc->name != func.name) {
						pos = rhfunctions.Find(func);
						if(pos) {
							pfunc = &rhfunctions.GetAt(pos);
							pfunc->points.RemoveAll();
						} else {
							rhfunctions.AddTail(func);
							pfunc = &rhfunctions.GetTail();
						}
						index+=2;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)func_str);
						if(i != INVALID) pfunc->type = i;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)rank_str);
						if(i != INVALID) pfunc->rank = i;
						rs.MoveNext();
						continue;
					}
					
					rs.GetFieldValue(index++, varValue);
					v.x = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					v.y = varValue.m_dblVal;
					pfunc->points.AddSorted(v);
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			/*response spectrum functions*/
			pfunc = NULL;
			if(OpenRecord(rs,_T( "SELECT * FROM \"Response Spectrum Functions\"ORDER BY name,period" ) ) ) {
				while( !rs.IsEOF( ) ) {
					
					index = 0;
					rs.GetFieldValue(index++, varValue);
					func.name = *varValue.m_pstring;
					
					if(!pfunc || pfunc->name != func.name) {
						pos = rsfunctions.Find(func);
						if(pos) {
							pfunc = &rsfunctions.GetAt(pos);
							pfunc->points.RemoveAll();
						} else {
							rsfunctions.AddTail(func);
							pfunc = &rsfunctions.GetTail();
						}
						index+=2;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)func_str);
						if(i != INVALID) pfunc->type = i;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)rank_str);
						if(i != INVALID) pfunc->rank = i;
						rs.GetFieldValue(index++, varValue);
						pfunc->Ag = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)soil_str);
						if(i != INVALID) pfunc->soil_type = i;
						
						rs.MoveNext();
						continue;
					}
					
					rs.GetFieldValue(index++, varValue);
					v.x = varValue.m_dblVal;
					rs.GetFieldValue(index++, varValue);
					v.y = varValue.m_dblVal;
					pfunc->points.AddSorted(v);
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			/*Analysis cases*/
			ANALYSISCASE acase,*pcase;
			if(OpenRecord(rs,_T( "SELECT * FROM \"Analysis Case Definitions\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					index = 0;
					rs.GetFieldValue(index++, varValue);
					acase.name = *varValue.m_pstring;
					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)analysis_case_str);
					if(i != INVALID) {
						acase.acase_type = i;
						rs.GetFieldValue(index++, varValue);
						acase.run = (*varValue.m_pstring == "YES");
						rs.GetFieldValue(index++, varValue);
						acase.finished = (*varValue.m_pstring == "YES");
						rs.GetFieldValue(index++, varValue);
						acase.index = (int) varValue.m_dblVal;
						
						if(pos = all_analysis_cases[i]->Find(acase)) {
							pcase = &all_analysis_cases[i]->GetAt(pos);
							*pcase = acase;
						} else {
							if(i == ANALYSISCASE::LOAD_CASE) {
								LOADCASE ldcase;
								ldcase.name = acase.name;
								ldcase.run = acase.run;
								ldcase.finished = acase.finished;
								ldcase.index = acase.index;
								((PLOADCASELIST)all_analysis_cases[i])->AddTail(ldcase);
							} else if(i == ANALYSISCASE::MODAL_CASE) {
								MODALCASE mcase;
								mcase.name = acase.name;
								mcase.run = acase.run;
								mcase.finished = acase.finished;
								mcase.index = acase.index;
								((PMODALCASELIST)all_analysis_cases[i])->AddTail(mcase);
							} else if(i == ANALYSISCASE::BUCKLING_CASE) {
								BUCKLINGCASE bcase;
								bcase.name = acase.name;
								bcase.run = acase.run;
								bcase.finished = acase.finished;
								bcase.index = acase.index;
								((PBUCKLINGCASELIST)all_analysis_cases[i])->AddTail(bcase);
							} else if(i == ANALYSISCASE::NL_CASE) {
								NLCASE ncase;
								ncase.name = acase.name;
								ncase.run = acase.run;
								ncase.finished = acase.finished;
								ncase.index = acase.index;
								((PNLCASELIST)all_analysis_cases[i])->AddTail(ncase);
							} else if(i == ANALYSISCASE::RESPONSEH_CASE) {
								RESPONSEHIST rhcase;
								rhcase.name = acase.name;
								rhcase.run = acase.run;
								rhcase.finished = acase.finished;
								rhcase.index = acase.index;
								((PRESPONSECASELIST)all_analysis_cases[i])->AddTail(rhcase);
							} else if(i == ANALYSISCASE::RESPONSES_CASE) {
								RESPONSESPEC rscase;
								rscase.name = acase.name;
								rscase.run = acase.run;
								rscase.finished = acase.finished;
								rscase.index = acase.index;
								((PRESPONSESPECLIST)all_analysis_cases[i])->AddTail(rscase);
							} else if(i == ANALYSISCASE::COMBO_CASE) {
								COMBINATION combo;
								combo.name = acase.name;
								combo.run = acase.run;
								combo.finished = acase.finished;
								combo.index = acase.index;
								((PCOMBOLIST)all_analysis_cases[i])->AddTail(combo);
							}
						}
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			/*Load cases*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Load Case Definitions\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					LOADCASE ldcase,*pcase;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					ldcase.name = *varValue.m_pstring;
					pos = loadcases.Find(ldcase);
					if(pos) {
						pcase = &loadcases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)design_case_str);
						if(i != INVALID) pcase->type = i;
						rs.GetFieldValue(index++, varValue);
						pcase->swm = varValue.m_dblVal;
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			/*Modal cases*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Modal Case Definitions\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					MODALCASE mcase,*pcase;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					mcase.name = *varValue.m_pstring;
					pos = modalcases.Find(mcase);
					if(pos) {
						pcase = &modalcases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						pcase->minm = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->maxm = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->runmodes = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->shift = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->tolerance = varValue.m_dblVal;
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			/*Buckling cases*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Non Linear Case Definitions\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					NLCASE ncase,*pcase;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					ncase.name = *varValue.m_pstring;
					pos = nlcases.Find(ncase);
					if(pos) {
						pcase = &nlcases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						pcase->nSteps = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->nIteration = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->tolerance = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->save_stiffness = (*varValue.m_pstring == "YES");
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			/*Buckling cases*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Buckling Case Definitions\"" ) ) ) {
				while( !rs.IsEOF( ) ) {
					BUCKLINGCASE bcase,*pcase;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					bcase.name = *varValue.m_pstring;
					pos = bucklingcases.Find(bcase);
					if(pos) {
						pcase = &bucklingcases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						pcase->nmodes = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->runmodes = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						pcase->tolerance = varValue.m_dblVal;
					}
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			DAMPING *pdamp;
			if(responsecases.GetCount()) {
				/*Response history*/
				RESPONSEHIST rhcase,*prhcase;
				MODALCASE mcase;
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response History Case Definitions\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rhcase.name = *varValue.m_pstring;
						pos = responsecases.Find(rhcase);
						prhcase = &responsecases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						mcase.name = *varValue.m_pstring;
						prhcase->modalcase = &modalcases.GetAt(modalcases.Find(mcase));
						
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)history_str);
						if(i != INVALID) prhcase->type = i;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)history_ana_str);
						if(i != INVALID) prhcase->ana_type = i;
						rs.GetFieldValue(index++, varValue);
						prhcase->N = (UINT)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						prhcase->dt = varValue.m_dblVal;
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Response history Damping functions*/
				prhcase = NULL;
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response History Damping\" ORDER BY name,frequency" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rhcase.name = *varValue.m_pstring;
						
						if(!prhcase || prhcase->name != rhcase.name) {
							pos = responsecases.Find(rhcase);
							prhcase = &responsecases.GetAt(pos);
							pdamp = &prhcase->damping;
							pdamp->dvalues.points.RemoveAll();
							
							index+=2;
							rs.GetFieldValue(index++, varValue);
							i = FindString(varValue,(CString*)damp_str);
							if(i != INVALID) pdamp->type = i;
							
							rs.GetFieldValue(index++, varValue);
							pdamp->w1 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->w2 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->e1 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->e2 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->fm = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->fk = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->cdamping = varValue.m_dblVal;
							rs.MoveNext();
							continue;
						}
						
						rs.GetFieldValue(index++, varValue);
						v.x = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						v.y = varValue.m_dblVal;
						pdamp->dvalues.points.AddSorted(v);
						
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Response history loading*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response History Loading\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SPECFUNC sfunc;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rhcase.name = *varValue.m_pstring;
						pos = responsecases.Find(rhcase);
						prhcase = &responsecases.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						func.name = *varValue.m_pstring;
						sfunc.function = &rhfunctions.GetAt(rhfunctions.Find(func));
						rs.GetFieldValue(index++, varValue);
						sfunc.dir = (int)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						sfunc.scale = varValue.m_dblVal;
						prhcase->funclist.AddTail(sfunc);
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			
			if(responsespecs.GetCount()) {
				
				/*Response spectrum*/
				MODALCASE mcase;
				RESPONSESPEC rscase,*prscase;
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response Spectrum Case Definitions\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rscase.name = *varValue.m_pstring;
						pos = responsespecs.Find(rscase);
						prscase = &responsespecs.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						mcase.name = *varValue.m_pstring;
						prscase->modalcase = &modalcases.GetAt(modalcases.Find(mcase));
						
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)comb_str);
						if(i != INVALID) prscase->modal_comb = i;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)comb_str);
						if(i != INVALID) prscase->dir_comb = i;
						
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Response Spectrum Damping functions*/
				prscase = NULL;
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response Spectrum Damping\" ORDER BY name,frequency" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rscase.name = *varValue.m_pstring;
						
						if(!prscase || prscase->name != rscase.name) {
							pos = responsespecs.Find(rscase);
							prscase = &responsespecs.GetAt(pos);
							pdamp = &prscase->damping;
							pdamp->dvalues.points.RemoveAll();
							
							index+=2;
							rs.GetFieldValue(index++, varValue);
							i = FindString(varValue,(CString*)damp_str);
							if(i != INVALID) pdamp->type = i;
							
							rs.GetFieldValue(index++, varValue);
							pdamp->w1 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->w2 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->e1 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->e2 = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->fm = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->fk = varValue.m_dblVal;
							rs.GetFieldValue(index++, varValue);
							pdamp->cdamping = varValue.m_dblVal;
							rs.MoveNext();
							continue;
						}
						
						rs.GetFieldValue(index++, varValue);
						v.x = varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						v.y = varValue.m_dblVal;
						pdamp->dvalues.points.AddSorted(v);
						
						rs.MoveNext( );
					}
					rs.Close();
				}
				
				/*Response spectrum loading*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Response Spectrum Loading\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SPECFUNC sfunc;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						rscase.name = *varValue.m_pstring;
						pos = responsespecs.Find(rscase);
						prscase = &responsespecs.GetAt(pos);
						rs.GetFieldValue(index++, varValue);
						func.name = *varValue.m_pstring;
						sfunc.function = &rsfunctions.GetAt(rsfunctions.Find(func));
						rs.GetFieldValue(index++, varValue);
						sfunc.dir = (int)varValue.m_dblVal;
						rs.GetFieldValue(index++, varValue);
						sfunc.scale = varValue.m_dblVal;
						prscase->funclist.AddTail(sfunc);
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			/*Joint Constraints*/
			CONSTRAINT constraint,*pcons = NULL;
			JOINT* pjoint;
			if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Constraints\" ORDER BY name,joints" ) ) ) {
				while( !rs.IsEOF( ) ) {
					index = 0;
					rs.GetFieldValue(index++, varValue);
					constraint.name = *varValue.m_pstring;
					
					if(!pcons || pcons->name != constraint.name) {
						pos = constraints.Find(constraint);
						if(pos) {
							pcons = &constraints.GetAt(pos);
						} else {
							constraints.AddTail(constraint);
							pcons = &constraints.GetTail();
						}
						index++;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)cons_str);
						if(i != INVALID) pcons->type = i;
						rs.GetFieldValue(index++, varValue);
						i = FindString(varValue,(CString*)rank_str);
						if(i != INVALID) pcons->rank = i;
						
						pcons->constraint = 0;
						for(i = IUX;i <= IRZ;i++) {
							rs.GetFieldValue(index++, varValue);
							if(*varValue.m_pstring == "YES") pcons->constraint |= (1 << i);
						}
						
						rs.MoveNext();
						continue;
					}
					
					rs.GetFieldValue(index++, varValue);
					pjoint = FindJoint(*varValue.m_pstring);
					pcons->jointlist.AddTail(pjoint);
					
					rs.MoveNext( );
				}
				rs.Close();
			}
			
			if(JOINT::TotalNumber) {
				/*Joint Load Assignments*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Loads\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SYSTEM sys;
						JLOAD load;
						LOADCASE ldcase;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						JOINT* joint = FindJoint(*varValue.m_pstring);
						if(joint) {
							rs.GetFieldValue(index++, varValue);
							ldcase.name = *varValue.m_pstring;
							pos = loadcases.Find(ldcase);
							if(pos) {
								load.loadcase = &loadcases.GetAt(pos);
								
								rs.GetFieldValue(index++, varValue);
								if(*varValue.m_pstring == "LOCAL") load.system = 0;
								else  {
									sys.name = *varValue.m_pstring;
									load.system = &systems.GetAt(systems.Find(sys));
								}

								rs.GetFieldValue(index++, varValue);
								i = FindString(varValue,(CString*)loadtype_str);
								if(i != INVALID) load.type = i;
								for(i = IUX;i <= IRZ;i++) {
									rs.GetFieldValue(index++, varValue);
									load.Q[i] = varValue.m_dblVal;
								}
								joint->load.AddTail(load);
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			if(MEMBER::TotalNumber) {
				/*Member Load Assignments*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Member Loads\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SYSTEM sys;
						LOAD load;
						LOADCASE ldcase;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						MEMBER* member = FindMember(*varValue.m_pstring);
						
						if(member) {
							rs.GetFieldValue(index++, varValue);
							ldcase.name = *varValue.m_pstring;
							pos = loadcases.Find(ldcase);
							if(pos) {
								load.loadcase = &loadcases.GetAt(pos);
								
								rs.GetFieldValue(index++, varValue);
								if(*varValue.m_pstring == "LOCAL") load.system = 0;
								else  {
									sys.name = *varValue.m_pstring;
									load.system = &systems.GetAt(systems.Find(sys));
								}

								rs.GetFieldValue(index++, varValue);
								i = FindString(varValue,(CString*)loadtype_str);
								if(i != INVALID) load.type = i;
								rs.GetFieldValue(index++, varValue);
								i = FindString(varValue,(CString*)dir_str);
								if(i != INVALID) load.dir = i;
								rs.GetFieldValue(index++, varValue);
								load.x = varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								load.P = varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								load.x1 = varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								load.P1 = varValue.m_dblVal;
								
								member->load.AddTail(load);
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			if(SLAB::TotalNumber) {
				/*Slab Load Assignments*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Slab Loads\"" ) ) ) {
					while( !rs.IsEOF( ) ) {
						SYSTEM sys;
						LOAD load;
						LOADCASE ldcase;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						SLAB* slab = FindSlab(*varValue.m_pstring);
						
						if(slab) {
							rs.GetFieldValue(index++, varValue);
							ldcase.name = *varValue.m_pstring;
							pos = loadcases.Find(ldcase);
							if(pos) {
								load.loadcase = &loadcases.GetAt(pos);
								
								rs.GetFieldValue(index++, varValue);
								if(*varValue.m_pstring == "LOCAL") load.system = 0;
								else  {
									sys.name = *varValue.m_pstring;
									load.system = &systems.GetAt(systems.Find(sys));
								}

								rs.GetFieldValue(index++, varValue);
								i = FindString(varValue,(CString*)loadtype_str);
								if(i != INVALID) load.type = i;
								rs.GetFieldValue(index++, varValue);
								i = FindString(varValue,(CString*)dir_str);
								if(i != INVALID) load.dir = i;
								rs.GetFieldValue(index++, varValue);
								load.P = varValue.m_dblVal;
								slab->load.AddTail(load);
							}
						}
						rs.MoveNext( );
					}
					rs.Close();
				}
			}
			/*Combinations*/
			if(OpenRecord(rs,_T( "SELECT * FROM \"Combination Definitions\"" ) ) ) {
				ANALYSISCASE acase;
				COMB_TYPE combo,*pcombo;
				CList<COMB_TYPE,COMB_TYPE&>* pcombolist[3];
				pcombolist[0] = (CList<COMB_TYPE,COMB_TYPE&>*) &combinations;
				pcombolist[1] = (CList<COMB_TYPE,COMB_TYPE&>*) &bucklingcases;
				pcombolist[2] = (CList<COMB_TYPE,COMB_TYPE&>*) &nlcases;

				while( !rs.IsEOF( ) ) {
					LOADCOMBO ldcombo;
					index = 0;
					rs.GetFieldValue(index++, varValue);
					combo.name = *varValue.m_pstring;
					for(i = 0;i < 3;i++) {
						if(pos = pcombolist[i]->Find(combo)) {
                            pcombo = &pcombolist[i]->GetAt(pos);
							break;
						}
					}

					rs.GetFieldValue(index++, varValue);
					if(pcombo->acase_type == ANALYSISCASE::COMBO_CASE) {
						i = FindString(varValue,(CString*)comb_type_str);
					} else if(pcombo->acase_type == ANALYSISCASE::NL_CASE) {
						i = FindString(varValue,(CString*)nlcase_type_str);
					} else {
						i = FindString(varValue,(CString*)buckle_type_str);
					}
					if(i != INVALID) pcombo->type = i;

					rs.GetFieldValue(index++, varValue);
					i = FindString(varValue,(CString*)analysis_case_str);
					rs.GetFieldValue(index++, varValue);
					acase.name = *varValue.m_pstring;
					if(pos = all_analysis_cases[i]->Find(acase)) {
						ldcombo.loadcase = &all_analysis_cases[i]->GetAt(pos);
					}
					rs.GetFieldValue(index++, varValue);
					ldcombo.FS = varValue.m_dblVal;
					
					pcombo->loadlist.AddTail(ldcombo);
					rs.MoveNext( );
				}
				rs.Close();
			}
			/*
			Analysis result
			*/
			if(AnalysisResult) {
				ANALYSISCASE acase,*pcase;
				JOINT* joint;
				MEMBER* member;
				SLAB* slab;
				int aindex;
				
				pos = joints.GetHeadPosition();
				while(pos) {
					joint = &joints.GetNext(pos);
					joint->AllocVectors();
				}
				pos = members.GetHeadPosition();
				while(pos) {
					member = &members.GetNext(pos);
					member->AllocVectors();
				}
				pos = slabs.GetHeadPosition();
				while(pos) {
					slab = &slabs.GetNext(pos);
					slab->AllocVectors();
				}
				
				/*Modal Eigen Values*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Modal Eigen Values\" ORDER BY name,mode" ) ) ) {
					while( !rs.IsEOF( ) ) {
						MODALCASE mcase,*pcase;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						mcase.name = *varValue.m_pstring;
						pos = modalcases.Find(mcase);
						if(pos) {
							index++;
							pcase = &modalcases.GetAt(pos);
							vec_alloc(pcase->eigvalue,pcase->runmodes);
							index = 2;
							for(UINT j = 0;j < pcase->runmodes;j++) {
								rs.GetFieldValue(index, varValue);
								pcase->eigvalue[j] = (UINT)varValue.m_dblVal;
								rs.MoveNext( );
							}
						} else {
							break;
						}
					}
					rs.Close();
				}
				/*Buckling modes*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Buckling modes\" ORDER BY name,mode" ) ) ) {
					while( !rs.IsEOF( ) ) {
						BUCKLINGCASE bcase,*pcase;
						index = 0;
						rs.GetFieldValue(index++, varValue);
						bcase.name = *varValue.m_pstring;
						pos = bucklingcases.Find(bcase);
						if(pos) {
							index++;
							pcase = &bucklingcases.GetAt(pos);
							vec_alloc(pcase->eigvalue,pcase->runmodes);
							index = 2;
							for(UINT j = 0;j < pcase->runmodes;j++) {
								rs.GetFieldValue(index, varValue);
								pcase->eigvalue[j] = (UINT)varValue.m_dblVal;
								rs.MoveNext( );
							}
						} else {
							break;
						}
					}
					rs.Close();
				}
				
				/*
				JOINT
				*/
				if(JOINT::TotalNumber) {
					/*joint displacements*/
					if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Displacements\" ORDER BY Joint,Case,CaseType,StepType,StepNum" ) ) ) {
						while( !rs.IsEOF( ) ) {
							index = 0;
							rs.GetFieldValue(index++, varValue);
							joint = FindJoint(*varValue.m_pstring);
							
							if(joint) {
								rs.GetFieldValue(++index, varValue);
								i = FindString(varValue,(CString*)analysis_case_str);
								if(i != INVALID) acase.acase_type = i;
								rs.GetFieldValue(--index, varValue);
								acase.name = *varValue.m_pstring;
								pcase = &all_analysis_cases[i]->GetAt(all_analysis_cases[i]->Find(acase));
								index+=3;
								aindex = pcase->index;
								rs.GetFieldValue(index++, varValue);
								aindex += (int)varValue.m_dblVal;
								
								joint->AllocAnalysis(aindex);
								for(i = IUX;i <= IRZ;i++) {
									rs.GetFieldValue(index++, varValue);
									joint->disps_all[aindex][i] = varValue.m_dblVal;
								}
							}
							rs.MoveNext( );
						}
						rs.Close();
					}
					
					/*joint reactions*/
					if(OpenRecord(rs,_T( "SELECT * FROM \"Joint Reactions\" ORDER BY Joint,Case,Case,CaseType,StepNum" ) ) ) {
						while( !rs.IsEOF( ) ) {
							index = 0;
							rs.GetFieldValue(index++, varValue);
							joint = FindJoint(*varValue.m_pstring);
							
							if(joint) {
								rs.GetFieldValue(++index, varValue);
								i = FindString(varValue,(CString*)analysis_case_str);
								if(i != INVALID) acase.acase_type = i;
								rs.GetFieldValue(--index, varValue);
								acase.name = *varValue.m_pstring;
								pcase = &all_analysis_cases[i]->GetAt(all_analysis_cases[i]->Find(acase));
								index+=3;
								aindex = pcase->index;
								rs.GetFieldValue(index++, varValue);
								aindex += (int)varValue.m_dblVal;
								
								for(i = IUX;i <= IRZ;i++) {
									rs.GetFieldValue(index++, varValue);
									joint->forces_all[aindex][i] = varValue.m_dblVal;
								}
							}
							rs.MoveNext( );
						}
						rs.Close();
					}
				}
				if(MEMBER::TotalNumber) {
					/*member internal forces*/
					if(OpenRecord(rs,_T( "SELECT * FROM \"Member Internal Forces\" ORDER BY Member,Case,Case,CaseType,StepNum,Station" ) ) ) {
						while( !rs.IsEOF( ) ) {
							index = 0;
							rs.GetFieldValue(index++, varValue);
							member = FindMember(*varValue.m_pstring);
							
							if(member) {
								index++;
								rs.GetFieldValue(++index, varValue);
								i = FindString(varValue,(CString*)analysis_case_str);
								if(i != INVALID) acase.acase_type = i;
								rs.GetFieldValue(--index, varValue);
								acase.name = *varValue.m_pstring;
								pcase = &all_analysis_cases[i]->GetAt(all_analysis_cases[i]->Find(acase));
								index+=3;
								aindex = pcase->index;
								rs.GetFieldValue(index++, varValue);
								aindex += (int)varValue.m_dblVal;
								member->AllocAnalysis(aindex);
								
								for(UINT k = 0;k < member->nDiv;k++) {
									index = 1;
									rs.GetFieldValue(index, varValue);
									member->station[k] = varValue.m_dblVal;
									index = 6;
									for(i = IUX;i <= IRZ;i++) {
										rs.GetFieldValue(index++, varValue);
										member->forces_all[aindex][k][i] = varValue.m_dblVal;
									}
									rs.MoveNext();
								}
							}
						}
						rs.Close();
					}
                }
            }
            
			/*
			Design Result
			*/

            if(DesignResult) {
				DESIGNDATA* pdata;
				MEMBER* member;

				/*Concrete design*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Concrete Design\" ORDER BY Member,Station" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						if(member = FindMember(*varValue.m_pstring)) {
							for(UINT k = 0;k < member->nDiv;k++) {
								pdata = &member->design_data[k];
								index = 1;
								rs.GetFieldValue(index++, varValue);
								member->station[k] = varValue.m_dblVal;
                                rs.GetFieldValue(index++, varValue);
								pdata->flags = (*varValue.m_pstring == "FAIL");
								rs.GetFieldValue(index++, varValue);
								pdata->Area[0] = varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								pdata->Area[1] = varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								pdata->As = varValue.m_dblVal;
								for(i = 0;i < 2;i++) {
									rs.GetFieldValue(index++, varValue);
									pdata->total[i] = (int)varValue.m_dblVal;
									rs.GetFieldValue(index++, varValue);
									pdata->count[i][0] = (int)varValue.m_dblVal;
									rs.GetFieldValue(index++, varValue);
									pdata->diam[i][0] = varValue.m_dblVal;
									rs.GetFieldValue(index++, varValue);
									pdata->count[i][1] = (int)varValue.m_dblVal;
									rs.GetFieldValue(index++, varValue);
									pdata->diam[i][1] = varValue.m_dblVal;
								}
								rs.GetFieldValue(index++, varValue);
								pdata->s = (int)varValue.m_dblVal;
								rs.GetFieldValue(index++, varValue);
								pdata->sd = varValue.m_dblVal;
								rs.MoveNext();
							}
						}
					}
					rs.Close();
				}
				/*Steel design*/
				if(OpenRecord(rs,_T( "SELECT * FROM \"Steel Design\" ORDER BY Member" ) ) ) {
					while( !rs.IsEOF( ) ) {
						index = 0;
						rs.GetFieldValue(index++, varValue);
						if(member = FindMember(*varValue.m_pstring)) {
							pdata = &member->design_data[0];
							rs.GetFieldValue(index++, varValue);
							pdata->flags = (*varValue.m_pstring == "FAIL");
							for(i = 0;i < 3;i++) {
								rs.GetFieldValue(index++, varValue);
								pdata->ratios[i] = varValue.m_dblVal;
							}
						}
						rs.MoveNext();
					}
					rs.Close();
				}
			}
			/*
			end
			*/
		}      
	} 
    CATCH(CDBException, e)
    {
        AfxMessageBox("Database error: " + e->m_strError);
    }
    END_CATCH;
	
	/*close*/
    database.Close();

	DisplacementResult = FALSE;

	/*fill cases*/
	if(AnalysisResult) {
		FillAnalysisCases();
		DetermineScale();
	} else {
		FillCases();
	}
	
	/*views*/
	CmyView* pView;
	POSITION pos = GetFirstViewPosition();
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
}