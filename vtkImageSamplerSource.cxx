/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkImageSamplerSource.cpp,v $
  Date:      $Date: 2002/12/30 22:27:27 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
#include "vtkImageSamplerSource.h"
//
#include "vtkObjectFactory.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataSet.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
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
#include "vtkExtentTranslator.h"
//
#include "vtkBoundsExtentTranslator.h"
//
#include "vtkDummyController.h"
//
#include <set>
#include <algorithm>
#include <functional>
//
vtkStandardNewMacro(vtkImageSamplerSource);
//----------------------------------------------------------------------------
#define REGULARGRID_SAXPY(a,x,y,z) \
  z[0] = a*x[0] + y[0]; \
  z[1] = a*x[1] + y[1]; \
  z[2] = a*x[2] + y[2]; 

#define REGULARGRID_ISAXPY(a,s,x2,y,z) \
  a[0] = (s[0]*x2[0])>0 ? (z[0] - y[0])/(s[0]*x2[0]) : 0; \
  a[1] = (s[1]*x2[1])>0 ? (z[1] - y[1])/(s[1]*x2[1]) : 0; \
  a[2] = (s[2]*x2[2])>0 ? (z[2] - y[2])/(s[2]*x2[2]) : 0;
// --------------------------------------------------------------------------------------
vtkImageSamplerSource::vtkImageSamplerSource(void) {
  this->GlobalOrigin[0]         = 0.0;
  this->GlobalOrigin[1]         = 0.0;
  this->GlobalOrigin[2]         = 0.0;
  this->LocalOrigin[0]          = 0.0;
  this->LocalOrigin[1]          = 0.0;
  this->LocalOrigin[2]          = 0.0;
  this->Spacing[0]              = 0.0;
  this->Spacing[1]              = 0.0;
  this->Spacing[2]              = 0.0;
  this->Delta                   = 0.0;
  this->Resolution[0]           = 0;
  this->Resolution[1]           = 0;
  this->Resolution[2]           = 0;
  this->WholeDimension[0]       = 1;
  this->WholeDimension[1]       = 1;
  this->WholeDimension[2]       = 1;
}

//----------------------------------------------------------------------------
int vtkImageSamplerSource::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}
//----------------------------------------------------------------------------
int vtkImageSamplerSource::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageSamplerSource::RequiredDataType()
{
  return VTK_IMAGE_DATA;
}
//----------------------------------------------------------------------------
int vtkImageSamplerSource::RequestDataObject(
  vtkInformation *, 
  vtkInformationVector  **vtkNotUsed(inputVector), 
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  bool ok = (output!=NULL);
  //
  ok = (ok && output->GetDataObjectType()==this->RequiredDataType());
  //
  vtkDataSet *newOutput = NULL;
  if (!ok) {
    switch (this->RequiredDataType()) {
      case VTK_IMAGE_DATA:
        newOutput = vtkImageData::New();
        break;
    }
    outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
    newOutput->Delete();
    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageSamplerSource::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  int piece, numPieces, ghostLevels;
  //
  piece       = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces   = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevels = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  //
  // Pass the piece request through
  //
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  std::cout << "Imagesampler (" << piece << ") set num pieces to " << numPieces << std::endl;
  
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageSamplerSource::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  this->ComputeInformation(request, inputVector, outputVector);
  //
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //
	int maxpieces = -1;//outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 
    0, this->WholeDimension[0]-1, 
    0, this->WholeDimension[1]-1, 
    0, this->WholeDimension[2]-1 );
//	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),maxpieces);
  outInfo->Set(vtkDataObject::ORIGIN(), this->GlobalOrigin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->spacing, 3);

  //
  //std::cout << "RI WHOLE_EXTENT {";
  //for (int i=0; i<3; i++) std::cout << WholeDimension[i] << (i<2 ? "," : "}");
  //std::cout << std::endl;

  return 1;
}
//----------------------------------------------------------------------------
void vtkImageSamplerSource::ComputeAxesFromBounds(vtkDataSet *inputData, double lengths[3], bool inflate)
{
  //
  // Define box...
  //
  vtkBoundingBox box;
  double bounds[6];
  inputData->GetBounds(bounds);
  box.SetBounds(bounds);
  //
  double bmin[3], bmn[3] = {bounds[0], bounds[2], bounds[4]};
  double bmax[3], bmx[3] = {bounds[1], bounds[3], bounds[5]};
  vtkSmartPointer<vtkMultiProcessController> Controller = vtkMultiProcessController::GetGlobalController();
  if (Controller == NULL) {
    Controller = vtkSmartPointer<vtkDummyController>::New();
  }
  Controller->AllReduce(bmn, bmin, 3, vtkCommunicator::MIN_OP);
  Controller->AllReduce(bmx, bmax, 3, vtkCommunicator::MAX_OP);
  box.SetMinPoint(bmin);
  box.SetMaxPoint(bmax);
  if (inflate) {
    box.Inflate(this->Delta);
  }
  box.GetMinPoint(this->GlobalOrigin[0], this->GlobalOrigin[1], this->GlobalOrigin[2]);
  box.GetLengths(lengths);
}
//----------------------------------------------------------------------------
int vtkImageSamplerSource::ComputeInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet  *inputData = inInfo ? vtkDataSet::SafeDownCast
    (inInfo->Get(vtkDataObject::DATA_OBJECT())) : NULL;
  //
  double lengths[3];
  //
  if (!inputData) {
    return 0;
  }
  this->ComputeAxesFromBounds(inputData, lengths, true);

  //
  // Define sampling box...
  // This represents the complete box which may be distributed over N pieces
  // Hence we compute WholeDimension (=WholeExtent+1)
  //
  double o2[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; i++) {
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
      this->WholeDimension[i] = vtkMath::Round(1.0/this->scaling[i] + 0.5);
    }
    else {
      this->WholeDimension[i] = 1;
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
// only valid for axis aligned grids
void vtkImageSamplerSource::BoundsToExtent(double *bounds, int *extent, int updatePiece) 
{
  double axesvec[3] = {1,1,1};
  vtkBoundingBox processRegion(bounds);
  vtkBoundingBox dataRegion;
  dataRegion.SetMinPoint(this->GlobalOrigin);
  dataRegion.SetMaxPoint(this->GlobalOrigin[0]+axesvec[0],
                         this->GlobalOrigin[1]+axesvec[1],
                         this->GlobalOrigin[2]+axesvec[2]);
  //
  if (!dataRegion.IntersectBox(processRegion)) {
    for (int i=0; i<3; i++) { extent[i*2  ] = 0; }
    for (int i=0; i<3; i++) { extent[i*2+1] = -1; }
    return;
  }
  const double *minvec = dataRegion.GetMinPoint();
  const double *maxvec = dataRegion.GetMaxPoint();
  //
  double updateExtentLo[3],updateExtentHi[3];
  REGULARGRID_ISAXPY(updateExtentLo, this->scaling, axesvec, this->GlobalOrigin, minvec);
  REGULARGRID_ISAXPY(updateExtentHi, this->scaling, axesvec, this->GlobalOrigin, maxvec);
  for (int i=0; i<3; i++) { extent[i*2  ] = static_cast<int>(updateExtentLo[i]+0.5); }
  for (int i=0; i<3; i++) { extent[i*2+1] = static_cast<int>(updateExtentHi[i]+0.5); }
  //
  std::cout << updatePiece << "Setting Extent to {";
  for (int i=0; i<6; i++) std::cout << extent[i] << (i<5 ? "," : "}");
  std::cout << std::endl;
}
//----------------------------------------------------------------------------
int vtkImageSamplerSource::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkImageData *outGrid = NULL;
  //
  switch (this->RequiredDataType()) {
    case VTK_IMAGE_DATA:
      outGrid = this->GetOutput();
      break;
  }
  //
  vtkInformation *inInfo = inputVector[0] ? inputVector[0]->GetInformationObject(0) : NULL;
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //
  int outUpdateExt[6];
  int outWholeExt[6] = {
    0, this->WholeDimension[0]-1, 
    0, this->WholeDimension[1]-1, 
    0, this->WholeDimension[2]-1 };
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  vtkBoundsExtentTranslator *bet = vtkBoundsExtentTranslator::SafeDownCast(translator);
  if (bet) {
    int updatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    double *bounds = bet->GetBoundsForPiece(updatePiece);
    // which one is right!
//    this->BoundsToExtent(bounds,outUpdateExt,updatePiece);
    bet->BoundsToExtentThreadSafe(bounds, outWholeExt,outUpdateExt); 
    std::cout << "Image sampler " << updatePiece << " Setting Extent to {";
    for (int i=0; i<6; i++) std::cout << outUpdateExt[i] << (i<5 ? "," : "}");
    std::cout << std::endl;
  }
  else {
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outUpdateExt);
  }
  //
  int dims[3]= {1+outUpdateExt[1]-outUpdateExt[0],
                1+outUpdateExt[3]-outUpdateExt[2], 
                1+outUpdateExt[5]-outUpdateExt[4]};
  //
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outWholeExt, 6); 
  outGrid->SetExtent(outUpdateExt);
  outGrid->SetOrigin(this->GlobalOrigin);
  outGrid->SetSpacing(this->spacing);
  return 1;
}
