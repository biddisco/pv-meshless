/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkSMCustomBoxRepresentationProxy.h,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSMCustomBoxRepresentationProxy - proxy for vtkCustomBoxRepresentation
// .SECTION Description
// vtkSMCustomBoxRepresentationProxy is a proxy for vtkCustomBoxRepresentation. A
// specialization is needed to set the tranform on the vtkCustomBoxRepresentation.

#ifndef __vtkSMCustomBoxRepresentationProxy_h
#define __vtkSMCustomBoxRepresentationProxy_h

#include "vtkSMWidgetRepresentationProxy.h"

class VTK_EXPORT vtkSMCustomBoxRepresentationProxy : public vtkSMWidgetRepresentationProxy
{
public:
  static vtkSMCustomBoxRepresentationProxy* New();
  vtkTypeRevisionMacro(vtkSMCustomBoxRepresentationProxy, vtkSMWidgetRepresentationProxy);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void UpdateVTKObjects()
    {this->Superclass::UpdateVTKObjects();}
  virtual void UpdatePropertyInformation();
  virtual void UpdatePropertyInformation(vtkSMProperty* prop)
    { this->Superclass::UpdatePropertyInformation(prop); }

//BTX
protected:
  vtkSMCustomBoxRepresentationProxy();
  ~vtkSMCustomBoxRepresentationProxy();

  // This method is overridden to set the transform on the vtkWidgetRepresentation.
  virtual void CreateVTKObjects();
  virtual void UpdateVTKObjects(vtkClientServerStream&);

private:
  vtkSMCustomBoxRepresentationProxy(const vtkSMCustomBoxRepresentationProxy&); // Not implemented
  void operator=(const vtkSMCustomBoxRepresentationProxy&); // Not implemented
//ETX
};

#endif

