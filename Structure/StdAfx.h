// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__0ACDC08E_BE38_11D5_8598_9B1670516327__INCLUDED_)
#define AFX_STDAFX_H__0ACDC08E_BE38_11D5_8598_9B1670516327__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#undef _WIN32_WINNT
#define _WIN32_WINNT    0x0500

#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers

#include <afxwin.h>         // MFC core and standard components
#include <afxext.h>         // MFC extensions
#include <afxdisp.h>        // MFC Automation classes
#include <afxdtctl.h>		// MFC support for Internet Explorer 4 Common Controls
#include <afxtempl.h>
#include <afxmt.h>
#include <afxdb.h>
#include <odbcinst.h>


#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			// MFC support for Windows Common Controls
#endif // _AFX_NO_AFXCMN_SUPPORT


#if _MSC_VER > 1200

#else

#ifdef _AFXDLL
#define BEGIN_TEMPLATE_MESSAGE_MAP1(theTemplate, theClass, baseClass) \
  template theTemplate const AFX_MSGMAP* PASCAL theClass::_GetBaseMessageMap() \
      { return &baseClass::messageMap; } \
  template theTemplate const AFX_MSGMAP* theClass::GetMessageMap() const \
      { return &theClass::messageMap; } \
  template theTemplate AFX_DATADEF const AFX_MSGMAP theClass::messageMap = \
  { &theClass::_GetBaseMessageMap, &theClass::_messageEntries[0] }; \
  template theTemplate const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
  { \

#else
#define BEGIN_TEMPLATE_MESSAGE_MAP1(theTemplate, theClass, baseClass) \
  template theTemplate const AFX_MSGMAP* theClass::GetMessageMap() const \
      { return &theClass::messageMap; } \
  template theTemplate AFX_DATADEF const AFX_MSGMAP theClass::messageMap = \
  { &baseClass::messageMap, &theClass::_messageEntries[0] }; \
  template theTemplate const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
  { \

#endif

#define BEGIN_TEMPLATE_MESSAGE_MAP(CDefineDia, T , CDialog) BEGIN_TEMPLATE_MESSAGE_MAP1(<class T>,CDefineDia<T>,CDialog)

#endif


//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__0ACDC08E_BE38_11D5_8598_9B1670516327__INCLUDED_)
