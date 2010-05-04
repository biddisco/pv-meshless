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
  stream  << vtkClientServerStream::Invoke
          << this->GetID()
          << "SetTransform"
          << this->GetSubProxy("Transform")->GetID()
          << vtkClientServerStream::End;
  vtkProcessModule::GetProcessModule()->SendStream(
    this->GetConnectionID(),
    this->GetServers(), 
    stream);
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::UpdateVTKObjects(vtkClientServerStream& stream)
{
  if (this->InUpdateVTKObjects)
    {
    return;
    }

  int something_changed = this->ArePropertiesModified();

  this->Superclass::UpdateVTKObjects(stream);

  if (something_changed)
    {
    stream  << vtkClientServerStream::Invoke
            << this->GetID()
            << "SetTransform"
            << this->GetSubProxy("Transform")->GetID()
            << vtkClientServerStream::End;
    }
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::UpdatePropertyInformation()
{
  vtkCustomBoxRepresentation* repr = vtkCustomBoxRepresentation::SafeDownCast(
    this->GetClientSideObject());
  vtkTransform* transform = vtkTransform::SafeDownCast(
    vtkProcessModule::GetProcessModule()->GetObjectFromID(
      this->GetSubProxy("Transform")->GetID()));
  repr->GetTransform(transform);

  this->Superclass::UpdatePropertyInformation();
}

//----------------------------------------------------------------------------
void vtkSMCustomBoxRepresentationProxy::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}


