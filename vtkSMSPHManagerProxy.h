/*=========================================================================

  Program:   ParaView
  Module:    vtkSMProxy.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSMSPHManagerProxy - proxy for a VTK object(s) on a server
// .SECTION Description
// Singleton for SPH manager

#ifndef __vtkSMSPHManagerProxy_h
#define __vtkSMSPHManagerProxy_h

#include "vtkSMProxy.h"

class vtkCallbackCommand;
class vtkSPHManager;

class VTK_EXPORT vtkSMSPHManagerProxy : public vtkSMProxy
{
public:
  static vtkSMSPHManagerProxy* New();
  vtkTypeMacro(vtkSMSPHManagerProxy, vtkSMProxy);

  static void SPHModifiedCallback(vtkObject* caller, unsigned long eid,
                                        void* clientdata, void* calldata);
  
  virtual void CreateVTKObjects();

protected:
   vtkSMSPHManagerProxy();
  ~vtkSMSPHManagerProxy();

  vtkCallbackCommand *SPHObserver;
  vtkSPHManager      *SPHSingleton;
  static vtkSMSPHManagerProxy *ReferenceProxy;
};

#endif
