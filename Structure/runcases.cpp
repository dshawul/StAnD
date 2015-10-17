#include "common.h"
#include "GUI.h"

enum {
	STATIC_DOF = 1,DYNAMIC_DOF = 2,RESTRAINED_DOF = 4,CONSTRAINED_DOF = 8
};

void RESPONSE::ConstructDampingMatrix(MATRAN C,UINT SJ) {
	DOUBLE en,wn;
	for(UINT i = 0;i < SJ;i++) {
		wn = sqrt(modalcase->eigvalue[i]);
		en = damping.dvalues.Get(wn,TRANSIENT,1,TRUE);
		C[i] = 2 * wn * en;
	}
}

UINT  ThreadProc(LPVOID pDoct) {
	CmyDocument* pDoc = (CmyDocument*)pDoct;
	WaitForSingleObject(pDoc->threadStart.m_hObject, INFINITE);
	pDoc->RunAnalysis();
    return 0;
}

void CmyDocument::PrintProgress(CString str,BOOL onsameline) {
	int index;
	if(pAnalysisProgressDia) {
		if(onsameline) {
			index = pAnalysisProgressDia->c_List.GetCurSel();
			pAnalysisProgressDia->c_List.DeleteString(index);
			index = pAnalysisProgressDia->c_List.AddString(str);
			pAnalysisProgressDia->c_List.SetCurSel(index);
		} else {
			index = pAnalysisProgressDia->c_List.AddString(str);
            pAnalysisProgressDia->c_List.SetCurSel(index);
		}
	}
}

void CmyDocument::PrintStylish(CString str,BOOL year) {
	CTime endTime;
	endTime = CTime::GetCurrentTime();
	if(year) str += endTime.Format( "%Y/%m/%d  %H:%M:%S");
	else str += endTime.Format( "%H:%M:%S");
	PrintProgress("");
	PrintProgress(str);
	PrintProgress("");
}

void CmyDocument::NumberJoint(JOINT* joint,UINT& N,UINT& SUN,UINT& DUN,UINT& RUN,INTEGER* tdof_flags) {

	if(joint->mconnect + joint->sconnect + joint->cconnect == 0)
		return;

	int i,count,jnumber;
	count = 0;
	for(i = 0;i < 6;i++) {
		
		if(!(static_dofs & (1 << i))) 
			continue;

		joint->number[i] = jnumber = N;
		N++;

		tdof_flags[jnumber] = 0;

		if(joint->restraint & (1 << i)) {
			tdof_flags[jnumber] |= RESTRAINED_DOF;
			RUN++;
			continue;
		}
	
		tdof_flags[jnumber] |= STATIC_DOF;
		SUN++;

		if(!EQUAL(joint->massembled[i],0) && (dynamic_dofs & (1 << i))) {
			tdof_flags[jnumber] |= DYNAMIC_DOF;
			DUN++;
		}
	}
}

void CmyDocument::NumberDofs(UINT& N,UINT& SUN,UINT& DUN,UINT& RUN,INTEGER* tdof_flags) {
	POSITION pos,pos1,pos2;
	MEMBER* member;
	SLAB* slab;
	JOINT* joint,*pjoint;
	CONSTRAINT* pconst;
	
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->NodeNumber = INVALID;
	}

	JOINTPLIST parents,neighbors,save;
	int JN = 0;

	while(true) {
		BOOL found = FALSE;
		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			if(joint->NodeNumber == INVALID) {
				found = TRUE;
				break;
			}
		}
		if(!found)
			break;

		joint->NodeNumber = JN++;
		NumberJoint(joint,N,SUN,DUN,RUN,tdof_flags);
		parents.RemoveAll();
		parents.AddTail(joint);
		
		while(parents.GetCount()) {
		/*
		Iterate for each parent
			*/
			save.RemoveAll();
			pos = parents.GetHeadPosition();
			while(pos) {
				joint = parents.GetNext(pos);
				
				neighbors.RemoveAll();

				/*neighbors by member*/
				pos1 = members.GetHeadPosition();
				while(pos1) {
					member = &members.GetNext(pos1);
					if(joint == member->j1) {
						if(member->j2->NodeNumber == INVALID) {
							neighbors.AddTail(member->j2);
						}
					} else if(joint == member->j2) {
						if(member->j1->NodeNumber == INVALID) {
							neighbors.AddTail(member->j1);
						}
					}
				}

				/*neighbors by slab*/
				pos1 = slabs.GetHeadPosition();
				while(pos1) {
					slab = &slabs.GetNext(pos1);
					for(UINT i = 0;i < slab->NJ;i++) {
						if(joint == slab->jt[i]) {
							for(UINT j = 0;j < slab->NJ;j++) {
								if(i == j) continue;
								if(slab->jt[j]->NodeNumber == INVALID) {
									if(!neighbors.Find(slab->jt[j])) neighbors.AddTail(slab->jt[j]);
								}
							}
						}
					}
				}

				/*neighbors by use of constraints as well!*/
				pos1 = constraints.GetHeadPosition();
				while(pos1) {
					pconst = &constraints.GetNext(pos1);
					if(joint == pconst->jointlist.GetHead()) {
						pos2 = pconst->jointlist.GetHeadPosition();
						while(pos2) {
							pjoint = pconst->jointlist.GetNext(pos2);
							if(joint != pjoint && pjoint->NodeNumber == INVALID) {
								if(!neighbors.Find(pjoint)) neighbors.AddTail(pjoint);
							}
						}
					} else if(pconst->jointlist.Find(joint)) {
						pjoint = pconst->jointlist.GetHead();
						if(joint != pjoint && pjoint->NodeNumber == INVALID) {
							if(!neighbors.Find(pjoint)) neighbors.AddTail(pjoint);
						}
					}
				}
				/*sort*/
				JOINT** pj1,**pj2,*temp;
				pos1 = neighbors.GetHeadPosition();
				while(pos1) {
					pj1 = &neighbors.GetNext(pos1);
					pos2 = pos1;
					while(pos2) {
						pj2 = &neighbors.GetNext(pos2);
						if(((*pj2)->mconnect + (*pj2)->sconnect + (*pj2)->cconnect) > 
						   ((*pj1)->mconnect + (*pj1)->sconnect + (*pj1)->cconnect)) {
							temp = *pj1;
							*pj1 = *pj2;
							*pj2 = temp;
						}
					}
				}
				/*number*/
				pos1 = neighbors.GetHeadPosition();
				while(pos1) {
					temp = neighbors.GetNext(pos1);
					save.AddTail(temp);
					temp->NodeNumber = JN++;
					NumberJoint(temp,N,SUN,DUN,RUN,tdof_flags);	
				}
			}
			
			/*fill new parents*/
			parents.RemoveAll();
			pos = save.GetHeadPosition();
			while(pos) {
				joint = save.GetNext(pos);
				parents.AddTail(joint);
			}
		}
	}
}

void CmyDocument::DetermineBandwidth(UINT N,UINT& NB) {
	POSITION pos;
	MEMBER* member;
	SLAB* slab;
	UINT i,j,k,delta,maxd = 0;
	CString str;

	/*members*/
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		for(i = 0;i < 6;i++) {
			if(member->j1->number[i] == INVALID) continue;
			delta = abs(member->j1->number[i] - member->j2->number[i]);
			if(delta > maxd) maxd = delta;
			break;
		}
	}

	/*slabs*/
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		for(i = 0;i < 6;i++) {
			if(slab->jt[0]->number[i] == INVALID) continue;
			break;
		}
		for(k = 0;k < slab->NJ;k++) {
			for(j = 0;j < slab->NJ;j++) {
				if(k == j) continue;
				delta = abs(slab->jt[k]->number[i] - slab->jt[j]->number[i]);
				if(delta > maxd) maxd = delta;
			}
		}
	}

	/*constraints*/
	CONSTRAINT* pconst;
	POSITION pos2;
	JOINT*joint,*pjoint;
	pos = constraints.GetHeadPosition();
	while(pos) {
		pconst = &constraints.GetNext(pos);

		pos2 = pconst->jointlist.GetHeadPosition();
		joint = pconst->jointlist.GetNext(pos2);
		while(pos2) {
			pjoint = pconst->jointlist.GetNext(pos2);
			
			for(i = 0;i < 6;i++) {
				if(joint->number[i] == INVALID) continue;
				delta = abs(joint->number[i] - pjoint->number[i]);
				if(delta > maxd) maxd = delta;
				break;
			}
		}
	}

    /*adjust NB*/
	NB = maxd + 1;
	for(i = 0;i < 6;i++) {
		if(static_dofs & (1 << i)) NB++;
	}
	if(NB > N) NB = N;

	str.Format("Bandwidth = %12d",NB);
	PrintProgress(str);
}
/*
Construct Non-Linear stiffness matrix
*/
void CmyDocument::FormK(MATRAN K,UINT N,UINT NB,BOOL formpenalty) {
	POSITION pos;
	MEMBER* member;
	SLAB* slab;
	CONSTRAINT* body;
	
	clear(K,NB * N);
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->AddToK(K,NB);
	}
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->AddToK(K,NB);
	}
	/*determine penalty weight*/
	if(formpenalty) {
		DOUBLE maxK = 0;
		for(UINT i = 0;i < N;i++) {
			if(fabs(K[BANDI(i,i)]) > maxK) maxK = fabs(K[BANDI(i,i)]); 
		}
		CONSTRAINT::PenaltyWeight = pow(10,log10(maxK) + 16 / 2);
        CONSTRAINT::PenaltyWeight = min(CONSTRAINT::PenaltyWeight,1e19);
	}
	/*Edge constraints*/
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		if(slab->edge_constraint)
			slab->AddEdgeConstraint(&joints,global,K,NB);
	}
	/*Other constraints*/
	pos = constraints.GetHeadPosition();
	while(pos) {
        body = &constraints.GetNext(pos);
		body->AddToK(K,NB);
	}
}
void CmyDocument::RestrainK(MATRAN K,UINT N,UINT NB,INTEGER* dof_flags,DOUBLE v1,DOUBLE v2) {
	UINT i,j,diag;

	for(i = 0;i < N;i++) {
		diag = BANDI(i,i);
		if(dof_flags[i] & RESTRAINED_DOF) {
			for(j = 0; j < NB;j++) {
				if(i + j < N) 
					K[NB - 1 - j + (i + j) * NB] = 0;
				K[NB - 1 - j + i * NB] = 0;
			}
			K[diag] = v1;
		} 
		if(EQUAL(K[diag], 0)) {
			K[diag] = v2;
		}
	}
}
BOOL CmyDocument::FactorK(MATRAN K,MATRAN FAC,UINT N,UINT NB,INTEGER* dof_flags,INTEGER* IPVT) {
	int code;
	
	/*Modify K*/
	equ(FAC,K,NB * N);
	RestrainK(FAC,N,NB,dof_flags);

	/*factor*/
	if(IPVT) {
		MATRAN Copy;
		vec_alloc(Copy,NB * N);
		equ(Copy,FAC,NB * N);
		convert_symm_full(Copy,FAC,N,NB,3 * NB - 2);
		vec_free(Copy);
		
		code = factor_banded_full(FAC,N,NB,IPVT);
		if(code) {
			PrintProgress("");
			PrintProgress("Factoring failed: Input matrix is singular!!");
			PrintProgress("");
			return FALSE;
		}
	} else {
		code = factor_banded(FAC,N,NB);
		if(code) {
			PrintProgress("");
			PrintProgress("Factoring failed: Matrix is not positive definite!!");
			PrintProgress("");
			return FALSE;
		}
	}
	return TRUE;
}
void CmyDocument::RestoreObjects() {
	POSITION pos,pos1;
	MEMBERPLIST* pmlist;
    SLABPLIST* pslist;
	MEMBER* member;
	SLAB* slab;

	pos = removed_members.GetHeadPosition();
	while(pos) {
		member = &removed_members.GetNext(pos);
		member->j1->mconnect++;
        member->j2->mconnect++;
		members.AddTail(*member);
	}

	pos = removed_slabs.GetHeadPosition();
	while(pos) {
		slab = &removed_slabs.GetNext(pos);
		for(UINT i = 0;i < slab->NJ;i++)
			slab->jt[i]->sconnect++;
		slabs.AddTail(*slab);
	}

	/*delete member objects*/
	pos = member_elements.GetHeadPosition();
	while(pos) {
		pmlist = member_elements.GetNext(pos);
        
		pos1 = pmlist->GetHeadPosition();
		while(pos1) {
			member = pmlist->GetNext(pos1);
			DeleteMember(member);
		}
		pmlist->RemoveAll();
		delete[] pmlist;
	}
	member_elements.RemoveAll();

	/*delete slab objects*/
	pos = slab_elements.GetHeadPosition();
	while(pos) {
		pslist = slab_elements.GetNext(pos);
        
		pos1 = pslist->GetHeadPosition();
		while(pos1) {
			slab = pslist->GetNext(pos1);
			DeleteSlab(slab);
		}
		pslist->RemoveAll();
		delete[] pslist;
	}
	slab_elements.RemoveAll();
}
void CmyDocument::OnRunAnalysis() {
	POSITION pos;
	JOINTPLIST list;
	MEMBER* member;
	SLAB* slab;
	UINT i;

	/*open analyis file*/
    CString sPath;
	sPath = GetPathName();
	if(sPath.IsEmpty()) {
		AfxMessageBox("Save file before continuing");
		return;
	}

	int len = sPath.GetLength();
	sPath.SetAt(len - 3,'A');
	sPath.SetAt(len - 2,'N');
	sPath.SetAt(len - 1,'A');
	open_file((LPCSTR)sPath , AnalysisResult ? TRUE : FALSE);

	/*Analyis dialog*/
	CAnalysisDia AnalysisDia(AfxGetMainWnd(),all_analysis_cases);
	if(AnalysisDia.DoModal() != IDC_RUN)
		return;

	/*Total Number*/
	UINT jTotalNumber = JOINT::TotalNumber;
	UINT mTotalNumber = MEMBER::TotalNumber;
    UINT sTotalNumber = SLAB::TotalNumber;
	UINT count;

	/*1.Mesh slabs =>mesh slabs first so that members will be meshed
    with joints generated by slab meshing*/
	slab_elements.RemoveAll();
	removed_slabs.RemoveAll();

	count = slabs.GetCount();
	pos = slabs.GetHeadPosition();
	for(i = 0;i < count;i++) {
		slab = &slabs.GetNext(pos);
		if(slab->Mesh(&slabs,&joints,&slab_elements)) {
			removed_slabs.AddTail(*slab);
			DeleteSlab(slab);
		}
	}

	/*2.Mesh members*/
	member_elements.RemoveAll();
	removed_members.RemoveAll();

	count = members.GetCount();
	pos = members.GetHeadPosition();
	for(i = 0;i < count;i++) {
		member = &members.GetNext(pos);
		if(member->is_curved) {
			if(member->MeshCurved(&members,&joints,&member_elements)) {
				removed_members.AddTail(*member);
				DeleteMember(member);
			}
		}
	}
	
	count = members.GetCount();
	pos = members.GetHeadPosition();
	for(i = 0;i < count;i++) {
		member = &members.GetNext(pos);
		if(member->Mesh(&members,&joints,&member_elements)) {
			removed_members.AddTail(*member);
			DeleteMember(member);
		}
	}
	
	/*start*/
	AfxBeginThread(ThreadProc,this,THREAD_PRIORITY_NORMAL);

	CAnalysisProgressDia Dia(AfxGetMainWnd(),&threadStart);
	pAnalysisProgressDia = &Dia;
    Dia.DoModal();
	pAnalysisProgressDia = NULL;
	AnalysisResult = TRUE;

	/*Total Number*/
	JOINT::TotalNumber = jTotalNumber;
	MEMBER::TotalNumber = mTotalNumber;
    SLAB::TotalNumber = sTotalNumber;

	/*display stress diagram*/
	Lock = TRUE;
	DisplacementResult = TRUE;
	state = SELECT;
	DetermineScale();

	CmyView* pView;
	pos = GetFirstViewPosition();
	while (pos != NULL) {
		pView = (CmyView*) GetNextView(pos);
		pView->SetTitle();
		pView->Invalidate();
	} 
}

/*
Run analysis
*/
void CmyDocument::RunAnalysis() {
	if(AnalysisResult)
		FreeAnalysis();

	PrintStylish("B  E  G  I  N     A  N  A  L  Y  S  I  S"
		"                                                                                           ",TRUE);

	/*variables*/
	CString str;
	UINT i,j,count;
	POSITION pos,pos1;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	UBMP32 N,DUN,SUN,RUN,NB;
	MATRAN K,FAC;
	INTEGER *dof_flags;
	INTEGER index;
	BOOL rundynamic;

	/*Total number of load cases*/
	ELEMENT::TotalLoadCases = loadcases.GetCount() + combinations.GetCount() + nlcases.GetCount() + responsespecs.GetCount();
	
	MODALCASE* modalcase;
	pos1 = modalcases.GetHeadPosition();
	while(pos1) {
		modalcase = &modalcases.GetNext(pos1);
		ELEMENT::TotalLoadCases += modalcase->maxm;
	}

	BUCKLINGCASE* bcase;
	pos1 = bucklingcases.GetHeadPosition();
	while(pos1) {
		bcase = &bucklingcases.GetNext(pos1);
		ELEMENT::TotalLoadCases += bcase->nmodes;
	}
	
	RESPONSEHIST* responsecase;
	pos1 = responsecases.GetHeadPosition();
	while(pos1) {
		responsecase = &responsecases.GetNext(pos1);
		ELEMENT::TotalLoadCases += responsecase->N;
	}
	
	ELEMENT::TotalLoadCases += 124;
	
	/*Number of divisions*/
	LOAD* pload;
	pos1 = members.GetHeadPosition();
	while(pos1) {
		member = &members.GetNext(pos1);
		count = 0;
		pos = member->load.GetHeadPosition();
		while(pos) {
			pload = &member->load.GetNext(pos);
			if(pload->type == CONCENTRATED) count += 2;
			else if(pload->type == TRAPEZOIDAL) count += 4;
		}
		member->nDiv = member->nMinDiv + count;
	}
	/*Allocate vectors*/
	pos = joints.GetHeadPosition();
	while(pos) {
		joint = &joints.GetNext(pos);
		joint->AllocVectors();
	}
	
	pos = members.GetHeadPosition();
	while(pos) {
		member = &members.GetNext(pos);
		member->AllocVectors();
	    member->AssignSections();
	}
	
	pos = slabs.GetHeadPosition();
	while(pos) {
		slab = &slabs.GetNext(pos);
		slab->AllocVectors();
	}
	
	str.Format("Total Number of Joints    =  %d\n", joints.GetCount());
	PrintProgress(str);
	str.Format("Total Number of Members   =  %d\n", members.GetCount());
	PrintProgress(str);
	str.Format("Total Number of Slabs     =  %d\n", slabs.GetCount());
	PrintProgress(str);

	/*Lump masses at joints*/
	LumpMasses();

	/*
	Number degree of freedoms
	*/
	PrintStylish("F  O  R  M  I  N  G    A  N  D    F  A  C  T  O  R  I  N  G    M  A  T  R  I  X"
		"                                                       ");

	/*number dofs*/
	INTEGER* tdof_flags = new INTEGER[6 * joints.GetCount()];

	pos = joints.GetHeadPosition();
	while(pos) {
        joint = &joints.GetNext(pos);
		for(int i = IUX;i <= IRZ;i++) {
		    joint->number[i] = INVALID;
		}
	}

	N = 0;
	DUN = 0;
	SUN = 0;
	RUN = 0;

	NumberDofs(N,SUN,DUN,RUN,tdof_flags);

	/*print survey result*/
	str.Format("Total number of equations  = %12d",N);
	PrintProgress(str);
	str.Format("Dynamic degrees of freedom = %12d",DUN);
	PrintProgress(str);
	str.Format("Static degrees of freedom  = %12d",SUN);
	PrintProgress(str);
	/*
	Determine band-width
	*/
	DetermineBandwidth(N,NB);

	/*
	Form stiffness matrix
	*/
	dof_flags = new INTEGER[N];
	for(i = 0;i < N;i++)
	    dof_flags[i] = tdof_flags[i];
	delete[] tdof_flags;


	vec_alloc(K,N * NB);
	vec_alloc(FAC,N * NB);

	FormK(K,N,NB,TRUE);
	FactorK(K,FAC,N,NB,dof_flags);
	   	
  	/*
	Apply self weight
	*/
	ApplyDeadLoad();

	index = 0;
	/*
	run non-linear static load case
	*/
   	SolveStaticNonLinearAnalysisCases(N,NB,dof_flags,index);
	/*
	run static load cases
	*/
   	SolveStaticAnalysisCases(K,FAC,N,NB,dof_flags,index);
	/*
	run buckling cases
	*/
   	SolveBucklingAnalysisCases(K,N,DUN,NB,dof_flags,index);
	/*
	Remove self-weight loads
	*/
	RemoveTemporaryLoads();
    /*
	Dyanmic analysis
	*/
	rundynamic = FALSE;
	ANALYSISCASE* pcase;
	pos1 = modalcases.GetHeadPosition();
	while(pos1) {
        pcase = &modalcases.GetNext(pos1);
		if(pcase->run) {
			rundynamic = TRUE;
			break;
		}
	}

	if(rundynamic && !(DUN == 0 || N == DUN)) {
		/*fill mass matrix*/
		VECTOR M;
		INTEGER* Dir;
		vec_alloc(M,N);
		Dir = new INTEGER[N];
		
		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			for(j = IUX;j <= IRZ;j++) {
				count = joint->number[j];
				if(count == INVALID) continue;
				if(dof_flags[count] & DYNAMIC_DOF) {
					M[count] = joint->massembled[j];
					Dir[count] = j;
				} else {
					M[count] = 0;
					Dir[count] = INVALID;
				}
			}
		}
		/*
		run modal analysis cases
		*/
		SolveModalAnalysisCases(K,M,N,DUN,NB,dof_flags,index);
		/*
		run time history analysis cases
		*/
		SolveTimeHistoryAnalysisCases(K,M,Dir,N,NB,dof_flags,index);
		/*
		Response spectrum
		*/
		SolveResponseSpectrumAnalysisCases(K,M,FAC,Dir,N,NB,dof_flags,index);
		
		vec_free(M);
		delete[] Dir;
    }
	/*
	Load combination
	*/
	SolveCombination(index);
	
	/*fill cases*/
	FillAnalysisCases();
    
	/*free eigenvectors*/
	pos = modalcases.GetHeadPosition();
	while(pos) {
		modalcase = &modalcases.GetNext(pos);
		if(modalcase->finished) {
			vec_free(modalcase->eigvector);
		}
	}

	pos = bucklingcases.GetHeadPosition();
	while(pos) {
		bcase = &bucklingcases.GetNext(pos);
		if(bcase->finished) {
			vec_free(bcase->eigvector);
		}
	}

	/*free memory*/
	vec_free(K);
	vec_free(FAC);
	delete[] dof_flags;

	/*close analysis file*/
	close_file();

	/*end*/
    PrintStylish("A  N  A  L  Y  S  I  S     C  O  M  P  L  E  T  E"
		"                                                                            ",TRUE);

	if(pAnalysisProgressDia)
		pAnalysisProgressDia->AnalysisEnd();
}
/*
Run design
*/
void CmyDocument::OnDesign(UINT nID) {

	BOOL design_selected = (nID == IDM_DESIGN_SELECTED);

	/*open design file*/
    CString sPath;
	sPath = GetPathName();
	if(sPath.IsEmpty()) {
		AfxMessageBox("Save file before continuing");
		return;
	}

	int len = sPath.GetLength();
	sPath.SetAt(len - 3,'D');
	sPath.SetAt(len - 2,'G');
	sPath.SetAt(len - 1,'N');
	open_file((LPCSTR)sPath , DesignResult ? TRUE : FALSE);

	/*start design*/
	print("\tDesign of beams and columns using Ethiopian Building Code Standards\n");
	print("\t===================================================================\n\n");
	print("\t                   DESIGN CASE: %s\n",c_AnalysisCase->name);
    print("\t                   ============\n\n",c_AnalysisCase->name);


	MEMBER* member;
	POSITION pos;
	int le_option;

	pos = members.GetHeadPosition();
	while(pos) {
        member = &members.GetNext(pos);
		if(design_selected && !member->sel) {
			continue;
		}
		print("\n\t======================\n\t| Member %s: %s |\n\t======================\n\n",
			member->name,(member->section->rebar.design == COULMN) ? "COULMN" : "  BEAM");

		if(detailing.le_option == DETAILING::LCALCULATED) {
			CalculateEffectiveLength(member);
		} else if(detailing.le_option == DETAILING::LORIGINAL){
			member->LeCalculated[0] = member->LeCalculated[1] = member->GetLength();
		}

		le_option = member->le_option;
		if(le_option == DETAILING::LPREF)
			le_option = detailing.le_option;

		if(le_option == DETAILING::LCALCULATED) {
			member->Le[0] = member->LeCalculated[0];
			member->Le[1] = member->LeCalculated[1];
		} else if(le_option == DETAILING::LORIGINAL) {
			member->Le[0] = member->Le[1] = member->GetLength();
		} else if(le_option == DETAILING::LASSIGNED) {
		}
		member->Design(&detailing);
	}

	DesignResult = TRUE;
	close_file();
}

/*
solve static analysis cases
*/
void CmyDocument::SolveStaticAnalysisCases(MATRAN KT,MATRAN FACT,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index) {

	PrintStylish("L  I  N  E  A  R     S  T  A  T  I  C     C  A  S  E  S"
		"                                                                                            ");
	MATRAN K,FAC;
	VECTOR Q,D,Qcopy;
	POSITION pos,pos1;
	CString str;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	LOADCOMBO ldcombo;
	LOADCOMBOLIST current_load_list;

	vec_alloc(Q,N);
	vec_alloc(D,N);
	vec_alloc(Qcopy,N);

	UINT i,code;
	LOADCASE* loadcase;
	pos1 = loadcases.GetHeadPosition();
	while(pos1) {
		if(!pAnalysisProgressDia) break;

        loadcase = &loadcases.GetNext(pos1);
		if(!loadcase->run || loadcase->finished) {
			index++;
			continue;
		}

		loadcase->index = index;

		/*stiffness from non-linear case*/
		if(loadcase->start_case == "") {
			FAC = FACT;
			K = KT;
		} else {
			NLCASE ncase,*pnl_case;
			ncase.name = loadcase->start_case;
			pnl_case = &nlcases.GetAt(nlcases.Find(ncase));

			K = pnl_case->end_stiffness;
			vec_alloc(FAC,NB * N);
			FactorK(K,FAC,N,NB,dof_flags);
		}
		/*construct load case list*/
		ldcombo.loadcase = loadcase;
		ldcombo.FS = 1;
		current_load_list.RemoveAll();
		current_load_list.AddTail(ldcombo);
		BASE_LOAD::cloadcombolist = &current_load_list;

		/*print start*/
		str.Format("case : %s",loadcase->name);
		PrintProgress(str);

		print("\n\t\t\t\t\t==========================\n\t\t\t\t\tLOAD CASE  %15s\n\t\t\t\t\t==========================\n\n",
			loadcase->name);

		/*apply load*/
		clear(Q,N);
		clear(D,N);

		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			joint->ApplyLoad(Q,D);
		}
		
		pos = members.GetHeadPosition();
		while(pos) {
			member = &members.GetNext(pos);
			member->ApplyLoad(Q);
		}

		pos = slabs.GetHeadPosition();
		while(pos) {
			slab = &slabs.GetNext(pos);
			slab->ApplyLoad(Q);
		}

		equ(Qcopy,Q,N);

		/*modify load vector of restrained dofs*/
		VECTOR Qf;
		vec_alloc(Qf,N);
		multiply_banded(K,D,Qf,N,NB);
		for(i = 0;i < N;i++) {
			if(dof_flags[i] & RESTRAINED_DOF) Q[i] = D[i];
			else if(dof_flags[i] & STATIC_DOF) Q[i] -= Qf[i];
		}
		clear(D,N);
		vec_free(Qf);	

  		/*
		SOLVE
		*/
		code = solve_from_factor_banded(FAC,N,NB,Q,D);

		if(code == 1) {
			PrintProgress("");
		    PrintProgress("Ill conditioned matrix!!");
			PrintProgress("");
		} else if(code == 2) {
			PrintProgress("");
		    PrintProgress("Unstable structure!!");
			PrintProgress("");
		}
		
		multiply_banded(K,D,Q,N,NB);
		for(i = 0;i < N;i++) {
			if(dof_flags[i] & RESTRAINED_DOF) {
				Q[i] -= Qcopy[i];
			}
		}
		
		/*end*/

		AllocMemory(index);
		CalculateForces(Q,D,N,index);
		PrintOutput();

		index++;

		/*stiffness from non-linear case*/
		if(loadcase->start_case != "") {
			vec_free(FAC);
		}
		/*finish*/
		loadcase->finished = TRUE;
	}

	vec_free(Q);
	vec_free(D);
	vec_free(Qcopy);
}
/*
solve static non-linear analysis cases
*/
void CmyDocument::SolveStaticNonLinearAnalysisCases(UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index) {
	PrintStylish("N O N - L  I  N  E  A  R     S  T  A  T  I  C     C  A  S  E  S"
		"                                                                               ");

	UINT step,i,code,LDFAC = 3 * NB - 2;
	CString str;
	VECTOR Q,D,Qcopy;
	POSITION pos,pos1;
	JOINT* joint;
	MEMBER* member;
	SLAB* slab;
	MATRAN Knl,FACnl;
	INTEGER* IPVT;

	vec_alloc(Q,N);
	vec_alloc(D,N);
	vec_alloc(Qcopy,N);
	vec_alloc(Knl,NB * N);
	vec_alloc(FACnl,LDFAC * N);
	IPVT = new INTEGER[N];

	NLCASE* nlcase;
	pos1 = nlcases.GetHeadPosition();
	while(pos1) {
		
		if(!pAnalysisProgressDia) break;

        nlcase = &nlcases.GetNext(pos1);
		if(!nlcase->run || nlcase->finished) {
			index++;
			continue;
		}

		nlcase->index = index;
		BASE_LOAD::cloadcombolist = &nlcase->loadlist;

		str.Format("case : %s",nlcase->name);
		PrintProgress(str);

		print("\n\t\t\t\t\t================================\n\t\t\t\t\tNON-LINEAR CASE  %15s\n\t\t\t\t\t================================\n\n",
			nlcase->name);

        /*apply load*/
		clear(Q,N);
		clear(D,N);

		pos = joints.GetHeadPosition();
		while(pos) {
			joint = &joints.GetNext(pos);
			joint->ApplyLoad(Q,D);
		}
		
		pos = members.GetHeadPosition();
		while(pos) {
			member = &members.GetNext(pos);
			member->ApplyLoad(Q);
		}

		pos = slabs.GetHeadPosition();
		while(pos) {
			slab = &slabs.GetNext(pos);
			slab->ApplyLoad(Q);
		}

		equ(Qcopy,Q,N);
				
		/*
		Incrementation
		*/
		VECTOR Qinc,Dinc,Qsum,Dsum,Qt,Dt;
		vec_alloc(Qinc,N);
        vec_alloc(Dinc,N);
		vec_alloc(Dsum,N);
		vec_alloc(Qsum,N);
		vec_alloc(Qt,N);
		vec_alloc(Dt,N);

		AllocMemory(index);

		/*modify load vector of restrained dofs*/
		multiply_banded(Knl,D,Qinc,N,NB);
		for(i = 0;i < N;i++) {
			if(dof_flags[i] & RESTRAINED_DOF) Q[i] = D[i];
			else if(dof_flags[i] & STATIC_DOF) Q[i] -= Qinc[i];
		}
		clear(D,N);

        /*form stiffness matrix*/		
		FormK(Knl,N,NB);
		FactorK(Knl,FACnl,N,NB,dof_flags,IPVT);
		
		/*limits*/
		const DOUBLE tolerance = nlcase->tolerance;
		const UINT NSTEPS = nlcase->nSteps;
		const UINT NITER = nlcase->nIteration;
		
        /*
		For each load step
		*/
		BOOL converged;
		DOUBLE factor,total_factor = 0;
		MEMBER::add_geometric_stiffness = TRUE;
		MEMBER::p_delta = nlcase->type;
		
		for(step = 1;step <= NSTEPS;step++) {

			/*print iteration info*/
			str.Format("Increment = %d",step);
			PrintProgress(str,(step == 1) ? FALSE : TRUE);

			while(total_factor < 1) {

				/*save*/
				converged = FALSE;
				equ(Dt,Dsum,N);

				/*increase load and save info*/
				factor = 1.0 / (1 << (step - 1));
				total_factor += factor;
				BASE_LOAD::load_step_factor = total_factor;

				multiply(Q,factor,Qinc,N);
				multiply(Q,total_factor,Qt,N);
								
				DOUBLE norm = 0;
				for(i = 0;i < N;i++) norm += Qt[i] * Qt[i];
				norm = sqrt(norm);

				/*
                newton-raphson
				*/
				for(UINT iter = 0;iter < NITER;iter++) {

					/*solve*/
					code = solve_from_factor_banded_full(FACnl,N,NB,IPVT,Qinc,Dinc);
					if(code == 1) {
						PrintProgress("");
						PrintProgress("Ill conditioned matrix!!");
						PrintProgress("");
					} else if(code == 2) {
						PrintProgress("");
						PrintProgress("Unstable structure!!");
						PrintProgress("");
					}
					multiply_banded(Knl,Dinc,Qinc,N,NB);
					CalculateForces(Qinc,Dinc,N,index,TRUE);
					
					add(Dsum,Dinc,N);
					
					/*form stiffness matrix*/
					FormK(Knl,N,NB);
					FactorK(Knl,FACnl,N,NB,dof_flags,IPVT);
					
					/*form residual force vector*/
					multiply_banded(Knl,Dsum,Qsum,N,NB);
					for(i = 0;i < N;i++) {
						if(dof_flags[i] & RESTRAINED_DOF) {
							Qinc[i] = 0;
						} else {
							Qinc[i] = Qt[i] - Qsum[i];
						}
					}
					
					/*converged*/
					for(i = 0;i < N;i++) {
						if(fabs(Qinc[i]) >= tolerance * norm)
							break;
					}
					if(i == N) {
						converged = TRUE;
						break;
					}
				}
                /*
				No convergence : retrieve previous state
				*/
				if(converged) {
					multiply_banded(Knl,Dsum,Qsum,N,NB);
					CalculateForces(Qsum,Dsum,N,index,FALSE);
				} else {
					equ(Dsum,Dt,N);
					total_factor -= factor;
					BASE_LOAD::load_step_factor = total_factor;

					multiply_banded(Knl,Dsum,Qsum,N,NB);
					CalculateForces(Qsum,Dsum,N,index,FALSE);
					break;
				}
			}
			if(converged)
				break;
		}
		/*
		end
		*/
		BASE_LOAD::load_step_factor = 1;
		multiply_banded(Knl,Dsum,Qsum,N,NB);
		for(i = 0;i < N;i++) {
			if(dof_flags[i] & RESTRAINED_DOF) {
				Qsum[i] -= Qcopy[i];
			}
		}
		CalculateForces(Qsum,Dsum,N,index,FALSE);
		MEMBER::add_geometric_stiffness = FALSE;

		/*print*/
		PrintOutput();
		index++;

		/*free*/
		vec_free(Qinc);
		vec_free(Dinc);
		vec_free(Qsum);
		vec_free(Dsum);
		vec_free(Qt);
		vec_free(Dt);

		/*save end stiffness*/
		if(nlcase->save_stiffness) {
			vec_alloc(nlcase->end_stiffness,NB * N);
			equ(nlcase->end_stiffness,Knl,NB * N);
		}
		/*finish*/
		nlcase->finished = TRUE;
	}

	delete[] IPVT;
	vec_free(Knl);
	vec_free(FACnl);
	vec_free(Q);
	vec_free(D);
	vec_free(Qcopy);
}
/*
solve modal analysis cases
*/
void CmyDocument::SolveModalAnalysisCases(MATRAN KT,MATRAN M,UINT N,UINT DUN,UINT NB,INTEGER* dof_flags,INTEGER& index) {
	PrintStylish("E  I  G  E  N     M  O  D  A  L     A  N  A  L  Y  S  I  S     C  A  S  E  S"
		"                                                              ");

	MATRAN K,Km;
	POSITION pos1;
	CString str;
	UINT i;
	VECTOR Q,D;

	vec_alloc(Q,N);
	vec_alloc(D,N);
	vec_alloc(Km,NB * N);


	MODALCASE* modalcase;
	pos1 = modalcases.GetHeadPosition();
	while(pos1) {
		if(!pAnalysisProgressDia) break;

        modalcase = &modalcases.GetNext(pos1);
		if(!modalcase->run || modalcase->finished) {
			index += modalcase->maxm;
			continue;
		}
		
		modalcase->index = index;
		BASE_LOAD::cloadcombolist = NULL;

		str.Format("case : %s",modalcase->name);
		PrintProgress(str);

		PrintProgress("");
		str.Format("Maximum number of eigen modes sought  =  %12d",modalcase->maxm);
		PrintProgress(str);
        str.Format("Minimum number of eigen modes sought  =  %12d",modalcase->minm);
		PrintProgress(str);

		/*stiffness from non-linear case*/
		if(modalcase->start_case == "") {
			K = KT;
		} else {
			NLCASE ncase,*pnl_case;
			ncase.name = modalcase->start_case;
			pnl_case = &nlcases.GetAt(nlcases.Find(ncase));
			K = pnl_case->end_stiffness;
		}
		equ(Km,K,N*NB);
		RestrainK(Km,N,NB,dof_flags);

		for(i = 0;i < N;i++) {
			if(M[i] == 0 || dof_flags[i] & RESTRAINED_DOF) {
				M[i] = EXTREMESMALLNUMBER;
			}
		}

		/*find mode shapes*/
		VECTOR EigVec;
		VECTOR EigVal;
		UINT ncv,emax = min(modalcase->maxm,DUN);
		ncv = emax + 1;

		vec_alloc(modalcase->eigvector,emax * N);
		vec_alloc(modalcase->eigvalue,emax);
		vec_alloc(EigVec,emax * N);
		vec_alloc(EigVal,emax);

		int code = find_eigen(Km,M,N,NB,emax,ncv,modalcase->shift,modalcase->tolerance,EigVal,EigVec);
		if(code) {
			PrintProgress("");
			PrintProgress("Problem with Eigen Solution!!");
			PrintProgress("");
		}
		
		/*sort mode shapes*/
		UINT current = 0;
		for(UINT n = 0;n < emax;n++) {
			if(EigVal[n] > 0) {
				modalcase->eigvalue[current] = EigVal[n]; 
				equ(&modalcase->eigvector[current * N],&EigVec[n * N],N); 
				current++;
				if(current >= modalcase->maxm)
					break;
			}
		}
		modalcase->runmodes = current;

		/*calcualte internal forces*/
		PrintProgress("");

		for(i = 0;i < modalcase->runmodes;i++) {
			clear(Q,N);
			clear(D,N);

			equ(D,&modalcase->eigvector[i * N],N);
			multiply_banded(K,D,Q,N,NB);

			str.Format("Found mode %12d  of  %12d, Eigen value = %.6e, Period = %.6f",
				i + 1,modalcase->runmodes,modalcase->eigvalue[i],(2 * PI) / sqrt(modalcase->eigvalue[i]));
			PrintProgress(str);

			print("\n\t\t\t\t\t==========================\n\t\t\t\t\tMODE SHAPE %15d\n\t\t\t\t\t==========================\n\n",i + 1);
			print("period\n======\n%.6f sec\n\n",(2 * PI) / sqrt(modalcase->eigvalue[i]));
			print("eigenvector\n===========\n");
			print_vector(D,N);
		
			AllocMemory(index);
			CalculateForces(Q,D,N,index);
			
			index++;
		}

		/*free memory*/
		vec_free(EigVec);
		vec_free(EigVal);

		/*finish*/
		modalcase->finished = TRUE;
	}

	vec_free(Km);
	vec_free(Q);
	vec_free(D);
}
/*
Buckling analysis cases
*/
void CmyDocument::SolveBucklingAnalysisCases(MATRAN KT,UINT N,UINT DUN,UINT NB,INTEGER* dof_flags,INTEGER& index) {
	PrintStylish("B  U  C  K  L  I  N  G     A  N  A  L  Y  S  I  S     C  A  S  E  S"
		"                                                                          ");

	MATRAN K,Km,G;
	POSITION pos1;
	CString str;
	UINT i;
	VECTOR Q,D;

	vec_alloc(Q,N);
	vec_alloc(D,N);
	vec_alloc(Km,NB * N);
	vec_alloc(G,NB * N);


	BUCKLINGCASE* bcase;
	pos1 = bucklingcases.GetHeadPosition();
	while(pos1) {
		if(!pAnalysisProgressDia) break;

        bcase = &bucklingcases.GetNext(pos1);
		if(!bcase->run || bcase->finished) {
			index += bcase->nmodes;
			continue;
		}
		
		bcase->index = index;
		BASE_LOAD::cloadcombolist = NULL;

		str.Format("case : %s",bcase->name);
		PrintProgress(str);

		PrintProgress("");
		str.Format("Number of buckling modes sought  =  %12d",bcase->nmodes);
		PrintProgress(str);

		/*stiffness from non-linear case*/
		if(bcase->start_case == "") {
			K = KT;
		} else {
			NLCASE ncase,*pnl_case;
			ncase.name = bcase->start_case;
			pnl_case = &nlcases.GetAt(nlcases.Find(ncase));
			K = pnl_case->end_stiffness;
		}
		equ(Km,K,N*NB);
		RestrainK(Km,N,NB,dof_flags);
		print_vector(Km,N,NB);

		/*stiffness from given loadcase*/
		CombineLoad(bcase);
		
		MEMBER::add_geometric_stiffness = TRUE;
		MEMBER::only_geometric_stiffness = TRUE;

		MEMBER::p_delta = MEMBER::CONSISTENT;
		FormK(G,N,NB);
		multiply(G,-1,N * NB);
		RestrainK(G,N,NB,dof_flags,0,EXTREMESMALLNUMBER);

		MEMBER::add_geometric_stiffness = FALSE;
		MEMBER::only_geometric_stiffness = FALSE;

		FreeMemory(bcase->index);

		/*find mode shapes*/
		VECTOR EigVec;
		VECTOR EigVal;
		UINT ncv,emax = min(bcase->nmodes,N - 1);
		ncv = emax + 1;

		vec_alloc(bcase->eigvector,emax * N);
		vec_alloc(bcase->eigvalue,emax);
		vec_alloc(EigVec,emax * N);
		vec_alloc(EigVal,emax);

		int code = find_buckling(Km,G,N,NB,emax,ncv,0,bcase->tolerance,EigVal,EigVec);
		if(code) {
			PrintProgress("");
			PrintProgress("Problem with Buckling Solution!!");
			PrintProgress("");
		}

		/*sort mode shapes*/
		UINT current = 0;
		for(UINT n = 0;n < emax;n++) {
			if(!EQUAL(EigVal[n],0) && EigVal[n] < 1 / SMALLNUMBER) {
				bcase->eigvalue[current] = EigVal[n]; 
				equ(&bcase->eigvector[current * N],&EigVec[n * N],N); 
				current++;
				if(current >= bcase->nmodes)
					break;
			}
		}
		bcase->runmodes = current;

		/*calcualte internal forces*/
		PrintProgress("");

		for(i = 0;i < bcase->runmodes;i++) {
			clear(Q,N);
			clear(D,N);

			equ(D,&bcase->eigvector[i * N],N);
			multiply_banded(K,D,Q,N,NB);

			str.Format("Buckling mode %12d  of  %12d, Eigen value = %.6e",i + 1,bcase->runmodes,bcase->eigvalue[i]);
			PrintProgress(str);

			print("\n\t\t\t\t\t==========================\n\t\t\t\t\tBUCKLING MODE %12d\n\t\t\t\t\t==========================\n\n",i + 1);
			print("factor\n======\n%.6g\n\n",bcase->eigvalue[i]);
			print("eigenvector\n===========\n");
			print_vector(D,N);
			AllocMemory(index);
			CalculateForces(Q,D,N,index);
			
			index++;
		}

		/*free memory*/
		vec_free(EigVec);
		vec_free(EigVal);

		/*finish*/
		bcase->finished = TRUE;
	}

	vec_free(Km);
	vec_free(G);
	vec_free(Q);
	vec_free(D);
}
/*
solve time-history analysis cases
*/
void CmyDocument::GetExcitationLoad(RESPONSEHIST* responsecase,VECTOR P,VECTOR M,INTEGER* Dir,UINT N,UINT SJ,VECTOR E,DOUBLE step) {
	POSITION pos;
	DOUBLE accel,sum;
	VECTOR tv;
	SPECFUNC* pfunc;
	UINT j,l,count;

	vec_alloc(tv,N);

	count = 0;
	pos = responsecase->funclist.GetHeadPosition();
	while(pos) {
		pfunc = &responsecase->funclist.GetNext(pos);
		accel = pfunc->function->Get(step,responsecase->type,pfunc->scale);
		for(j = 0;j < N;j++) {
			if(Dir[j] == pfunc->dir) 
				tv[j] += (-M[j] * accel);
		}
		count++;
	}
	
	for(l = 0;l < SJ;l++) {
		sum = 0;
		for(j = 0;j < N;j++) {
			sum += E[l * N + j] * tv[j];
		}
		P[l] = sum;
	}

	vec_free(tv);
}
void CmyDocument::SolveTimeHistoryAnalysisCases(MATRAN KT,MATRAN M,INTEGER* Dir,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index) {
	PrintStylish("T  I  M  E     H  I  S  T  O  R  Y     A  N  A  L  Y  S  I  S     C  A  S  E  S"
		"                                                           ");

	MATRAN K;
	POSITION pos1;
	CString str;
	UINT i,j;
	VECTOR Q,D;

	vec_alloc(Q,N);
	vec_alloc(D,N);

	RESPONSEHIST* responsecase;
	MODALCASE* modalcase;

	pos1 = responsecases.GetHeadPosition();
	while(pos1) {
		if(!pAnalysisProgressDia) break;

        responsecase = &responsecases.GetNext(pos1);
		modalcase = responsecase->modalcase;
		if(!responsecase->run || 
			responsecase->finished ||
			!modalcase->finished) {
			index += responsecase->N;
			continue;
		}
		
		responsecase->index = index;
		BASE_LOAD::cloadcombolist = NULL;

		str.Format("case : %s",responsecase->name);
		PrintProgress(str);
		PrintProgress("");

		str.Format("Using modes from case:   %s",responsecase->modalcase->name);
		PrintProgress(str);
		str.Format("Time step size              =  %12.6f",responsecase->dt);
		PrintProgress(str);
        str.Format("Total number of time steps  =  %12d",responsecase->N);
		PrintProgress(str);
		str.Format("Total time length           =  %12.6f",responsecase->dt * responsecase->N);
		PrintProgress(str);

		print("\n\t\t\t\t\t==========================\n\t\t\t\t\t   Time History Analysis\n\t\t\t\t\t==========================\n\n");

		/*stiffness from non-linear case*/
		if(responsecase->start_case == "") {
			K = KT;
		} else {
			NLCASE ncase,*pnl_case;
			ncase.name = responsecase->start_case;
			pnl_case = &nlcases.GetAt(nlcases.Find(ncase));
			K = pnl_case->end_stiffness;
		}
		/*
		number of modes to consider
		*/
		UINT SJ = responsecase->modalcase->runmodes;
		MATRAN E = modalcase->eigvector;

		/*
		form modal mass/stiffness/damping matrix
		*/
        MATRAN Md,Kd,Cd;
		vec_alloc(Md,SJ);
		vec_alloc(Kd,SJ);
		vec_alloc(Cd,SJ);

		for(i = 0;i < SJ;i++) {
			Md[i] = 1.0;
			Kd[i] = modalcase->eigvalue[i];
		}
		responsecase->ConstructDampingMatrix(Cd,SJ);

		/*
		determine step size
		*/
		DOUBLE dt;
		UINT dN,Ninc;
		
		dt = 2.0 / sqrt(responsecase->modalcase->eigvalue[SJ - 1]);
		dt /= 2;
		Ninc = int(ceil(responsecase->dt / dt));
		dN = int(responsecase->N * Ninc) + 1;
		dt = responsecase->dt / (Ninc);

		print("Increment = %.5f\n\n",dt);
		
		/*
		Modal Superposition => Central-Difference Method
		*/
		MATRIX disp;
		MATRAN Kp,a,b;
		MATRIX u,tu,q,tq;
		VECTOR du,dq,ddq,P;
		DOUBLE sum;
		
		vec_alloc(Kp,SJ);
		vec_alloc(a,SJ);
		vec_alloc(b,SJ);
		
		mat_alloc(tu,dN + 1,N);
		u = &tu[1];
		mat_alloc(disp,dN,N);
		vec_alloc(du,N);
		mat_alloc(tq,dN + 1,SJ);
		q = &tq[1];
		

		vec_alloc(dq,SJ);
		vec_alloc(ddq,SJ); 
		vec_alloc(P,SJ);
		
		//initial values
		for(i = 0;i < N;i++) {
			u[0][i] = 0;
			du[i] = 0;
		}
		
		//modal displacement
		for(i = 0;i < SJ;i++) {
			sum = 0;
			for(j = 0;j < N;j++) {
				sum += E[i * N + j] * M[j] * u[0][j];
			}
			q[0][i] = sum;
		}
		
		//modal velocity
		for(i = 0;i < SJ;i++) {
			sum = 0;
			for(j = 0;j < N;j++) {
				sum += E[i * N + j] * M[j] * du[j];
			}
			dq[i] = sum;
		}
		
		//get initial excitation load
		GetExcitationLoad(responsecase,P,M,Dir,N,SJ,E,0);
		
		//modal acceleration
		for(i = 0;i < SJ;i++) {
			ddq[i] = P[i] - Cd[i] * dq[i] - Kd[i] * q[0][i];
		}
		
		//qminus1
		for(i = 0;i < SJ;i++) {
			q[-1][i] = q[0][i] - dt * dq[i] - (dt * dt / 2) * ddq[i];
		}
				
		//Kp
		for(i = 0;i < SJ;i++) {
			Kp[i] = 1 / (dt * dt) + Cd[i] / (2 * dt);
		}
		
		//a
		for(i = 0;i < SJ;i++) {
			a[i] = 1 / (dt * dt) - Cd[i] / (2 * dt);
		}
		
		//b
		for(i = 0;i < SJ;i++) {
			b[i] = Kd[i] -  2 / (dt * dt);
		}
		
		/*
	    Iterate over each time step
		*/
		for(UINT step = 0;step < dN;step++) {

			//save output displacment
			if(!(step % Ninc)) {
				print("t = %.5f\n",step * dt);
				print_vector(u[step],N);
				equ(disp[step],u[step],N);
			}
			
			if(step == dN - 1) 
				continue;
						
			//excitation load at current time
			GetExcitationLoad(responsecase,P,M,Dir,N,SJ,E,step * dt);
			
			//determine displacement in modal coordinates
			for(i = 0;i < SJ;i++) {
				q[step + 1][i] = (P[i] - a[i] * q[step - 1][i] - b[i] * q[step][i]) / Kp[i];
			}

			//actual displacment
			for(j = 0;j < N;j++) {
				sum = 0;
				for(i = 0;i < SJ;i++) {
					sum += E[i * N + j] * q[step + 1][i];
				}
				u[step + 1][j] = sum;
			}
		}
		
		/*free memory*/
		vec_free(tu);
		vec_free(tq);
		vec_free(du);
		vec_free(dq);
		vec_free(ddq);
		vec_free(P);
		vec_free(Kp);
		vec_free(a);
		vec_free(b);
		
		/*Calculate internal forces*/
		for(i = 0;i < dN;i += Ninc) {
			clear(Q,N);
			clear(D,N);
			
			equ(D,disp[i],N);
			multiply_banded(K,D,Q,N,NB);
			
			print("\n\t\t\t==========================\n\t\t\tTime Step %15d\n\t\t\t==========================\n\n",(i / Ninc));
			print("Joint Displacemnt\n=================\n");
			print_vector(D,N);
			
			AllocMemory(index);
			CalculateForces(Q,D,N,index);
			
			index++;
		}
		
        /*end*/
		vec_free(Md);
		vec_free(Kd);
		vec_free(Cd);
		vec_free(disp);
        
		/*finish*/
		responsecase->finished = TRUE;
	}

	vec_free(Q);
	vec_free(D);
}

/*
solve response spectrum analysis cases
*/
void CmyDocument::SolveResponseSpectrumAnalysisCases(MATRAN K,MATRAN M,MATRAN FAC,INTEGER* Dir,UINT N,UINT NB,INTEGER* dof_flags,INTEGER& index) {
	PrintStylish("R  E  S  P  O  N  S  E     S  P  E  C  T  R  U  M     A  N  A  L  Y  S  I  S     C  A  S  E  S"
		"                                  ");

	POSITION pos,pos1;
	CString str;
	UINT i,j,code;
	VECTOR Q,D;
	vec_alloc(Q,N);
	vec_alloc(D,N);

	RESPONSESPEC* responsespec;
	MODALCASE* modalcase;
	pos1 = responsespecs.GetHeadPosition();
	while(pos1) {
		if(!pAnalysisProgressDia) break;

        responsespec = &responsespecs.GetNext(pos1);
		modalcase = responsespec->modalcase;
		if(!responsespec->run || 
			responsespec->finished ||
			!modalcase->finished) {
			index++;
			continue;
		}
		
		responsespec->index = index;
		BASE_LOAD::cloadcombolist = NULL;

		str.Format("case : %s",responsespec->name);
		PrintProgress(str);
		PrintProgress("");

		str.Format("Using modes from case:   %s",responsespec->modalcase->name);
		PrintProgress(str);
		str.Format("Number of modes to be used  =  %12d",responsespec->modalcase->maxm);
		PrintProgress(str);

		print("\n\t\t\t\t\t==========================\n\t\t\t\t\tResponse Spectrum Analysis\n\t\t\t\t\t==========================\n\n");

		/*Construct Modal matrices*/
		UINT SJ = responsespec->modalcase->runmodes;
		MATRAN E = modalcase->eigvector;

		/*for modal analysis*/
		VECTOR s;
		VECTOR sn;
		VECTOR tn;

		vec_alloc(sn,SJ * N);
		vec_alloc(tn,N);
		vec_alloc(s,N);
				
		/*
		RSA
		*/
		index++;

		DOUBLE accel,period;
		SPECFUNC* pfunc;
		pos = responsespec->funclist.GetHeadPosition();
		UINT cc = 0;
		while(pos) {
			pfunc = &responsespec->funclist.GetNext(pos);

			index++;

			cc++;
			print("Spectrum Loading Case %3d\n=========================\n",cc);

			//influence vector
			for(i = 0;i < N;i++) {
				if(Dir[i] == pfunc->dir) 
					s[i] = -M[i];
				else s[i] = 0;
			}
		
			//modal contribution factor
			for(i = 0;i < SJ;i++) {
				tn[i] = 0;
				for(j = 0;j < N;j++) {
					tn[i] += E[i * N + j] * s[j];
				}
			}
		
			//modal expansion of excitation vector
			for(i = 0;i < SJ;i++) {
				for(j = 0;j < N;j++) {
					sn[i * N + j] = E[i * N + j] * M[j] * tn[i];
				}
			}
			
			/*for each mode*/
			for(UINT mode = 0;mode < SJ;mode++) {
				
				period = 2.0 * PI / sqrt(responsespec->modalcase->eigvalue[mode]);
				accel = pfunc->function->Get(period,TRANSIENT,pfunc->scale);
				
				/*clear*/
				clear(Q,N);
				clear(D,N);

				/*apply equivalent static load*/
				for(i = 0;i < N;i++) {
					Q[i] = sn[mode * N + i] * accel;
				}

				/*print*/
				print("\nEquivalent Mode %d Static Load\n============================\n",mode + 1);
				for(i = 0;i < N;i++) {
					print("%15.3f\n",Q[i]);
				}
				/*
				SOLVE
				*/
				code = solve_from_factor_banded(FAC,N,NB,Q,D);
				if(code == 1) {
					PrintProgress("");
					PrintProgress("Ill conditioned matrix!!");
					PrintProgress("");
				} else if(code == 2) {
					PrintProgress("");
					PrintProgress("Unstable structure!!");
					PrintProgress("");
				}
				multiply_banded(K,D,Q,N,NB);
				/*end*/

				AllocMemory(index);
				CalculateForces(Q,D,N,index);
				
				index++;
				/*end*/
			}
		}
		/*Modal combination*/
		CombineModes(responsespec);
		index = responsespec->index + 1;

		/*free*/
		vec_free(s);
		vec_free(sn);
		vec_free(tn);

		/*finish*/
		responsespec->finished = TRUE;
	}

	vec_free(Q);
	vec_free(D);
}

/*
solve combinatin cases
*/
void CmyDocument::SolveCombination(INTEGER& index) {
	PrintStylish("L  O  A  D     C  O  M  B  I  N  A  T  I  O  N"
		"                                                                                                     ");

	POSITION pos,pos1;
	CString str;
    COMBINATION* combo;
	LOADCOMBO* loadcombo;
	BOOL combine;

	pos = combinations.GetHeadPosition();
	while(pos) {
		combo = &combinations.GetNext(pos);
		if(!pAnalysisProgressDia) break;

		/*check if all load cases are run*/
		combine = TRUE;
        pos1 = combo->loadlist.GetHeadPosition();
		while(pos1) {
			loadcombo = &combo->loadlist.GetNext(pos1);
			if(!loadcombo->loadcase->finished) {
				combine = FALSE;
				break;
			}
		}

		if(!combine || combo->finished) {
			index++;
			continue;
		}

		str.Format("case : %s",combo->name);
		PrintProgress(str);
		PrintProgress("");

		combo->index = index;

		CombineLoad(combo);
        
		index++;
        if(combo->type == ENVELOPECOMBO) 
			index++;
		/*finish*/
		combo->finished = TRUE;
	}
}
