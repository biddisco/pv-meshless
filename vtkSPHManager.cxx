/*=========================================================================

  Project                 : vtkCSCS
  Module                  : vtkSPHManager.h
  Revision of last commit : $Rev: 1584 $
  Author of last commit   : $Author: soumagne $
  Date of last commit     : $Date:: 2010-02-26 18:09:23 +0100 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing 
  1) This copyright notice appears on all copies of source code 
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it 
  must not be reformatted such that the indentation, bracketing or 
  overall style is modified significantly. 

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
#include "vtkSPHManager.h"
//
#include "vtkObjectFactory.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtkstd/vector>
//
#include "vtkSmartPointer.h"

//----------------------------------------------------------------------------
// #ifdef JB_DEBUG__
//----------------------------------------------------------------------------
#ifdef NO_WIN32
  #define OUTPUTTEXT(a) vtkOutputWindowDisplayText(a);
#else
  #define OUTPUTTEXT(a) std::cout << (a) << std::endl;
#endif

#undef vtkDebugMacro
#define vtkDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg a << endl; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }

#undef vtkErrorMacro
#define vtkErrorMacro(a) vtkDebugMacro(a)
// #endif
//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkSPHManager, "$Revision: 1584 $");
vtkStandardNewMacro(vtkSPHManager);
//----------------------------------------------------------------------------
vtkSPHManager::vtkSPHManager() 
{
}
//----------------------------------------------------------------------------
vtkSPHManager::~vtkSPHManager()
{ 
}
//----------------------------------------------------------------------------
void vtkSPHManager::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
