/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkSMCustomBoxRepresentationProxy.cxx,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMCustomBoxRepresentationProxy.h"

#include "vtkCustomBoxRepresentation.h"
#include "vtkClientServerStream.h"
#include "vtkObjectFactory.h"
#include "vtkProcessModule.h"
#include "vtkTransform.h"

vtkStandardNewMacro(vtkSMCustomBoxRepresentationProxy);
vtkCxxRevisionMacro(vtkSMCustomBoxRepresentationProxy, "$Revision: 1.2 $");
//----------------------------------------------------------------------------
vtkSMCustomBoxRepresentationProxy::vtkSMCustomBoxRepresentationProxy()
{
}

//----------------------------------------------------------------------------
vtkSMCustomBoxRepresentationProxy::~vtkSMCustomBoxRepresentationProxy()
{
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::CreateVTKObjects()
{
  if (this->ObjectsCreated)
    {
    return;
    }

  this->Superclass::CreateVTKObjects();

  vtkClientServerStream stream;
  stream << vtkClientServerStream::Invoke
         << VTKOBJECT(this)
         << "SetTransform"
         << VTKOBJECT(this->GetSubProxy("Transform"))
         << vtkClientServerStream::End;
  this->ExecuteStream(stream);
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::UpdateVTKObjects(vtkClientServerStream& stream)
{
  if (this->InUpdateVTKObjects)
    {
    return;
    }

  int something_changed = this->ArePropertiesModified();

  this->Superclass::UpdateVTKObjects();

  if (something_changed)
    {
    vtkClientServerStream stream;
    stream << vtkClientServerStream::Invoke
      << VTKOBJECT(this)
      << "SetTransform"
      << VTKOBJECT(this->GetSubProxy("Transform"))
      << vtkClientServerStream::End;
    this->ExecuteStream(stream);
    }
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::UpdatePropertyInformation()
{
  vtkCustomBoxRepresentation* repr = vtkCustomBoxRepresentation::SafeDownCast(
    this->GetClientSideObject());
  vtkTransform* transform = vtkTransform::SafeDownCast(
      this->GetSubProxy("Transform")->GetClientSideObject());
  repr->GetTransform(transform);

  this->Superclass::UpdatePropertyInformation();
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}


