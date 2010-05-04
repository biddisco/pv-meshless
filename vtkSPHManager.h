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
// .NAME vtkSPHManager - Create/Expose an Xdmf DSM to Paraview
// .SECTION Description
// Create/Expose an Xdmf DSM to Paraview

#ifndef __vtkSPHManager_h
#define __vtkSPHManager_h

#include "vtkToolkits.h"     // For VTK_USE_MPI 
#include "vtkObject.h"       // Base class

class VTK_EXPORT vtkSPHManager : public vtkObject
{
public:
  static vtkSPHManager *New();
  vtkTypeRevisionMacro(vtkSPHManager,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);   


protected:
     vtkSPHManager();
    ~vtkSPHManager();

private:
    vtkSPHManager(const vtkSPHManager&);  // Not implemented.
    void operator=(const vtkSPHManager&);  // Not implemented.
};

#endif
