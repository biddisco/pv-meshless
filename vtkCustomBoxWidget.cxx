/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCustomBoxWidget.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCustomBoxWidget.h"
#include "vtkCustomBoxRepresentation.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkCustomBoxWidget, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkCustomBoxWidget);
//----------------------------------------------------------------------------
vtkCustomBoxWidget::vtkCustomBoxWidget()
{
  this->RotationEnabled = 0;
}
//----------------------------------------------------------------------------
vtkCustomBoxWidget::~vtkCustomBoxWidget()
{  
}

//----------------------------------------------------------------------
void vtkCustomBoxWidget::CreateDefaultRepresentation()
{
  if ( ! this->WidgetRep )
    {
    this->WidgetRep = vtkCustomBoxRepresentation::New();
    }
}
//----------------------------------------------------------------------------


