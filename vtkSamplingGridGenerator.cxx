/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkSamplingGridGenerator.cpp,v $
  Date:      $Date: 2002/12/30 22:27:27 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
#include "vtkSamplingGridGenerator.h"
#include "vtkObjectFactory.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataSet.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkImageData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkBoundingBox.h"
#include "vtkVariantArray.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkOBBTree.h"
#include "vtkImplicitFunction.h"
#include "vtkPlane.h"
#include "vtkBox.h"
#include "vtkImageData.h"
#include "vtkCutter.h"
#include "vtkTransform.h"
#include "vtkHomogeneousTransform.h"
//
vtkCxxRevisionMacro(vtkSamplingGridGenerator, "$Revision: 1.4 $");
vtkStandardNewMacro(vtkSamplingGridGenerator);
vtkCxxSetObjectMacro(vtkSamplingGridGenerator, CutFunction, vtkImplicitFunction);
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkSG_Copy(vtkDataSet *d) {
  vtkSmartPointer<vtkDataSet> result;
  result.TakeReference(d->NewInstance());
  result->ShallowCopy(d);
  return result;
}
// --------------------------------------------------------------------------------------
vtkSamplingGridGenerator::vtkSamplingGridGenerator(void) 
{
  this->CutFunction = NULL;
  this->Box = vtkImageData::New();
  this->Box->SetDimensions(2,2,2);
  this->Cutter = vtkCutter::New();
  this->Cutter->SetInput(this->Box);
}
// --------------------------------------------------------------------------------------
vtkSamplingGridGenerator::~vtkSamplingGridGenerator(void) 
{
  this->SetCutFunction(NULL);
  this->Box->Delete();
  this->Cutter->Delete();
}
//----------------------------------------------------------------------------
int vtkSamplingGridGenerator::RequiredDataType()
{
  if (this->GetCutFunction() && vtkPlane::SafeDownCast(this->GetCutFunction())) {
//    return VTK_POLY_DATA;
  }
  return vtkRegularGridSource::RequiredDataType();
}
//----------------------------------------------------------------------------
int vtkSamplingGridGenerator::ComputeInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  
  vtkDataSet *inData = inInfo ? vtkDataSet::SafeDownCast
    (inInfo->Get(vtkDataObject::DATA_OBJECT())) : NULL;

  double bounds[6], lengths[3], zerovec[3] = {0.0, 0.0, 0.0};
  vtkImplicitFunction *cf = this->GetCutFunction();
  if (cf && cf->IsA("vtkPlane")) {
    vtkPlane *plane = vtkPlane::SafeDownCast(cf);
    //
    inData->GetBounds(bounds);
    vtkBoundingBox box(bounds);
    box.Inflate(this->Delta);
    this->Box->SetOrigin(bounds[0], bounds[2], bounds[4]);
    this->Box->SetSpacing(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]);
    this->Cutter->SetCutFunction(plane);
    this->Cutter->SetInput(this->Box);
    this->Cutter->Update();
    vtkSmartPointer<vtkPolyData>  polys = this->Cutter->GetOutput();
    vtkSmartPointer<vtkPoints> pts = polys->GetPoints();
    double p1[3], p2[3];
    if (pts->GetNumberOfPoints()>0) pts->GetPoint(0,this->origin);
    else this->Box->GetOrigin(this->origin);
    if (pts->GetNumberOfPoints()>1) pts->GetPoint(1,p1);
    else memcpy(p1, this->origin, 3*sizeof(double));
    if (pts->GetNumberOfPoints()>2) pts->GetPoint(2,p2);
    else memcpy(p2, this->origin, 3*sizeof(double));
    //
    vtkMath::Subtract(p1, this->origin, this->axesvectors[0]);
    vtkMath::Subtract(p2, this->origin, this->axesvectors[1]);
    vtkMath::Subtract(this->origin, this->origin, this->axesvectors[2]);
    for (int i=0; i<3; i++) {
      lengths[i] = sqrt(vtkMath::Distance2BetweenPoints(zerovec, axesvectors[i]));
    }
  }
  else if (cf && cf->IsA("vtkBox")) {
    // Get origin, p1, p2 from implicit box
    vtkBox *iBox = vtkBox::SafeDownCast(cf);
    iBox->GetBounds(bounds);
    inData->GetBounds(bounds);
    vtkBoundingBox box(bounds);
//    box.Inflate(this->Delta);
    vtkAbstractTransform *trans = iBox->GetTransform();
    vtkAbstractTransform *itrans = trans->GetInverse();

    this->origin[0] = bounds[0];
    this->origin[1] = bounds[2];
    this->origin[2] = bounds[4];
    double p0[3] = { bounds[1], bounds[2], bounds[4]};
    double p1[3] = { bounds[0], bounds[3], bounds[4]};
    double p2[3] = { bounds[0], bounds[2], bounds[5]};
    vtkMath::Subtract(p0, this->origin, axesvectors[0]); 
    vtkMath::Subtract(p1, this->origin, axesvectors[1]); 
    vtkMath::Subtract(p2, this->origin, axesvectors[2]); 
//    double p0[3] = { 1,0,0};
//    double p1[3] = { 0,1,0};
//    double p2[3] = { 0,0,1};
    //
    itrans->TransformPoint(this->origin,this->origin);
    itrans->TransformVectorAtPoint(p0,axesvectors[0],axesvectors[0]);
    itrans->TransformVectorAtPoint(p1,axesvectors[1],axesvectors[1]);
    itrans->TransformVectorAtPoint(p2,axesvectors[2],axesvectors[2]);
    lengths[0] = vtkMath::Norm(axesvectors[0]);
    lengths[1] = vtkMath::Norm(axesvectors[1]);
    lengths[2] = vtkMath::Norm(axesvectors[2]);
//    trans->TransformPoint(lengths,lengths);
    //
  }
  //
  // Define sampling box...
  //
  double o2[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; i++) {
    if (this->Resolution[i]>1) {
      for (int j=0; j<3; j++) o2[j] += axesvectors[i][j];
    }
    if (this->Spacing[i]<=0.0 && this->Resolution[i]>1) {
      this->scaling[i] = 1.0/(this->Resolution[i]-1);
      this->spacing[i] = lengths[i]*this->scaling[i];
    }
    else{
      this->spacing[i] = this->Spacing[i];
      if (lengths[i]>0.0) {
        this->scaling[i] = this->Spacing[i]/lengths[i];
      }
      else {
        this->scaling[i] = 0.0;
      }
    }
    if (this->scaling[i]>0.0 && this->spacing[i]>1E-8) {
      this->Dimension[i] = vtkMath::Round(1.0/this->scaling[i] + 0.5);
    }
    else {
      this->Dimension[i] = 1;
    }
  }
  for (int i=0; i<3; i++) {
    this->centre[i] = this->origin[i] + o2[i]/2.0;
  }
  //
  this->NumPoints = Dimension[0]*Dimension[1]*Dimension[2];
  //
  return 1;
}
//----------------------------------------------------------------------------

