/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkSPHImageResampler.cpp,v $
  Date:      $Date: 2002/12/30 22:27:27 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
#include "vtkSPHImageResampler.h"
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
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"
//
#include "vtkBoundsExtentTranslator.h"
#include "vtkSPHProbeFilter.h"
//
#include "vtkDummyController.h"
//
#include <set>
#include <algorithm>
#include <functional>
#include <sstream>
//
vtkStandardNewMacro(vtkSPHImageResampler);
vtkCxxSetObjectMacro(vtkSPHImageResampler, Controller, vtkMultiProcessController);
vtkCxxSetObjectMacro(vtkSPHImageResampler, SPHManager, vtkSPHManager);
//
#include "vtkZoltanV1PartitionFilter.h"
//----------------------------------------------------------------------------
#if 0

  #define OUTPUTTEXT(a) std::cout <<(a); std::cout.flush();

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    if (this->UpdatePiece>=0) { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg << this->UpdatePiece << " : " a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    } \
  }

  #undef  vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
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
//----------------------------------------------------------------------------
vtkSPHImageResampler::vtkSPHImageResampler(void) 
{
  this->UpdatePiece               = 0;
  this->UpdateNumPieces           = 1;
  this->GlobalResamplingOrigin[0] = 0.0;
  this->GlobalResamplingOrigin[1] = 0.0;
  this->GlobalResamplingOrigin[2] = 0.0;
  this->Spacing[0]                = 0.0;
  this->Spacing[1]                = 0.0;
  this->Spacing[2]                = 0.0;
  this->Delta                     = 0.0;
  this->Resolution[0]             = 0;
  this->Resolution[1]             = 0;
  this->Resolution[2]             = 0;
  this->WholeDimensionWithoutDelta[0] = this->WholeDimensionWithDelta[0] = 1;
  this->WholeDimensionWithoutDelta[0] = this->WholeDimensionWithDelta[1] = 1;
  this->WholeDimensionWithoutDelta[0] = this->WholeDimensionWithDelta[2] = 1;
  //
  this->DensityScalars              = NULL;
  this->MassScalars                 = NULL;
  this->VolumeScalars               = NULL;
  this->HScalars                    = NULL;
  this->ModifiedNumber              = 0;
  this->ComputeDensityFromNeighbourVolume = 0;
  this->BoundsInitialized           = false;
  //
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
  //
  this->SPHManager    = vtkSPHManager::New();
  this->SPHProbe      = vtkSmartPointer<vtkSPHProbeFilter>::New();
  this->ProbeProgress = vtkSmartPointer<vtkSPHProbeProgress>::New();
  this->ProbeProgress->Self = this;
  this->ProbeProgress->Offset = 0.0;
  this->ProbeProgress->Scale  = 1.0;
  this->SPHProbe->AddObserver(vtkCommand::ProgressEvent, this->ProbeProgress);
  this->ExtentTranslator = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
}
//----------------------------------------------------------------------------
vtkSPHImageResampler::~vtkSPHImageResampler(void) 
{
  this->SetController(NULL);
}
//----------------------------------------------------------------------------
unsigned long vtkSPHImageResampler::GetMTime()
{
  unsigned long mtime = this->Superclass::GetMTime();
  unsigned long sm_mtime = this->SPHManager->GetMTime();
  mtime = mtime > sm_mtime? mtime : sm_mtime;
  return mtime;
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkSPHImageResampler::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation  *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  this->UpdatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
/*
  //
  int updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), updateExtent);
  std::stringstream temp;
  temp << "Image sampler " << this->UpdatePiece << " UPDATE_EXTENT {";
  for (int i=0; i<3; i++) temp << updateExtent[i*2] << "," << updateExtent[i*2+1] << (i<2 ? ":" : "}");
  vtkDebugMacro( << temp.str() );
*/
  return 1;
}

//----------------------------------------------------------------------------
int vtkSPHImageResampler::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation  *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet       *input = inInfo ? vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT())) : NULL;
  //
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

  // Get the extent translator from upstream, assumes data partitioned
  vtkBoundsExtentTranslator *bet = inInfo ? vtkBoundsExtentTranslator::SafeDownCast(
    inInfo->Get(vtkBoundsExtentTranslator::META_DATA())) : NULL;

  // copy the extent translator to the output for futher use downstream
  outInfo->Set(vtkBoundsExtentTranslator::META_DATA(), this->ExtentTranslator);

  // Work out how big our output will be
  this->ComputeGlobalInformation(input, this->WholeDimensionWithoutDelta, false);
  this->ComputeGlobalInformation(input, this->WholeDimensionWithDelta, true);

  //
  //int maxpieces = -1;//outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 
    0, this->WholeDimensionWithDelta[0]-1, 
    0, this->WholeDimensionWithDelta[1]-1, 
    0, this->WholeDimensionWithDelta[2]-1 );

  outInfo->Set(vtkDataObject::ORIGIN(), this->GlobalResamplingOrigin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->spacing, 3);

  inInfo->Set(vtkZoltanV1PartitionFilter::ZOLTAN_SAMPLE_ORIGIN(), 0.0, 0.0, 0.0);
  inInfo->Append(vtkExecutive::KEYS_TO_COPY(),vtkZoltanV1PartitionFilter::ZOLTAN_SAMPLE_ORIGIN());

  std::cout << "Setting Zoltan origin request " << std::endl;

  return 1;
}
//----------------------------------------------------------------------------
void vtkSPHImageResampler::ComputeAxesFromBounds(vtkDataSet *inputData, double samplinglengths[3], double datalengths[3], bool inflate)
{
  //
  // Define box...
  //
  double bounds[6];
  inputData->GetBounds(bounds);
  this->BoundsInitialized = vtkMath::AreBoundsInitialized(bounds);
  this->GlobalResamplingBounds.SetBounds(bounds);
  //
  double bmin[3], bmn[3] = {bounds[0], bounds[2], bounds[4]};
  double bmax[3], bmx[3] = {bounds[1], bounds[3], bounds[5]};
  Controller->AllReduce(bmn, bmin, 3, vtkCommunicator::MIN_OP);
  Controller->AllReduce(bmx, bmax, 3, vtkCommunicator::MAX_OP);
  this->GlobalResamplingBounds.SetMinPoint(bmin);
  this->GlobalResamplingBounds.SetMaxPoint(bmax);
  this->GlobalDataBounds = this->GlobalResamplingBounds;
  if (inflate) {
    this->GlobalResamplingBounds.Inflate(this->Delta);
  }
  this->GlobalResamplingBounds.GetMinPoint(this->GlobalResamplingOrigin[0], this->GlobalResamplingOrigin[1], this->GlobalResamplingOrigin[2]);
  this->GlobalResamplingBounds.GetLengths(samplinglengths);
  this->GlobalDataBounds.GetLengths(datalengths);
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::ComputeGlobalInformation(vtkDataSet *input, int dimensions[3], bool inflate)
{
  double samplinglengths[3], datalengths[3];
  this->ComputeAxesFromBounds(input, samplinglengths, datalengths, inflate);

  //
  // The data will sit inside some bounding box, but the resampling operation might exist in a larger space (delta>0)
  // so the extent we compute will not map directly to the extent of the data.
  //
  if (this->BoundsInitialized) {
    //
    // Define sampling box...
    // This represents the complete box which may be distributed over N pieces
    // Hence we compute WholeDimension (=WholeExtent+1)
    //
    for (int i=0; i<3; i++) {
      if (this->Spacing[i]<=0.0 && this->Resolution[i]>1) {
        this->scaling[i] = 1.0/(this->Resolution[i]-1);
        this->spacing[i] = samplinglengths[i]*this->scaling[i];
        samplinglengths[i] = vtkMath::Round(samplinglengths[i]/this->spacing[i] + 0.5)*this->spacing[i];
      }
      else{
        this->spacing[i] = this->Spacing[i];
        if (samplinglengths[i]>0.0) {
          this->scaling[i] = this->Spacing[i]/samplinglengths[i];
        }
        else {
          this->scaling[i] = 0.0;
        }
      }
      if (this->scaling[i]>0.0 && this->spacing[i]>1E-8) {
        dimensions[i] = vtkMath::Round(1.0/this->scaling[i] + 0.5);
        if (this->GlobalResamplingOrigin[i]>0.0) {
          this->GlobalResamplingOrigin[i] = vtkMath::Round(this->GlobalResamplingOrigin[i]/this->spacing[i] + 0.5)*this->spacing[i];
        }
        else if (this->GlobalResamplingOrigin[i]<0.0) {
          this->GlobalResamplingOrigin[i] = -vtkMath::Round(-this->GlobalResamplingOrigin[i]/this->spacing[i] + 0.5)*this->spacing[i];
        }
      }
      else {
        dimensions[i] = 1;
      }
    }
  }
  else {
    for (int i=0; i<3; i++) {
      dimensions[i] = this->Resolution[i];
    }
  }
  std::stringstream temp;
  temp << "Image sampler " << this->UpdatePiece << " set dimensions to {";
  for (int i=0; i<3; i++) temp << dimensions[i] << (i<2 ? "," : "}");
  vtkDebugMacro( << temp.str() );
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::ComputeLocalInformation(vtkBoundsExtentTranslator *bet, int localExtents[6])
{
  int WholeExt[6] = {
    0, this->WholeDimensionWithDelta[0]-1, 
    0, this->WholeDimensionWithDelta[1]-1, 
    0, this->WholeDimensionWithDelta[2]-1 };
  int WholeExtNoDelta[6] = {
    0, this->WholeDimensionWithoutDelta[0]-1, 
    0, this->WholeDimensionWithoutDelta[1]-1, 
    0, this->WholeDimensionWithoutDelta[2]-1 };
  //
  // The Local Extent must be calculated carefully by taking the extent we would have
  // if no delta had been added, then extending each piece outwards to the bounds
  // but without extending the edges which are not exterior faces of the domain
  //
  // 1) get the local extent based on delta=0
  int     localExt[6]; 
  double *localbounds = NULL;
  if (bet) {
    localbounds = bet->GetBoundsForPiece(this->UpdatePiece);
    bet->BoundsToExtentThreadSafe(localbounds, WholeExtNoDelta, localExt); 
  }

  // 2) compare local extent for zero delta with the global extent for zero delta
  //    any bound which is a local max/min in a dimension must be pushed to the global max/min
  for (int i=0; i<3; i++) {
    this->ExtentOffset[i] = (GlobalDataBounds.GetMinPoint()[i]-GlobalResamplingBounds.GetMinPoint()[i])/this->spacing[i];
    std::cout << "Extent offset " << i << " " << this->ExtentOffset[i] << std::endl;
    // min bound
    if (localExt[i*2] == WholeExtNoDelta[i*2]) {
      localExtents[i*2] = WholeExt[i*2];
      localbounds[i*2]  = this->GlobalResamplingBounds.GetBound(i*2);
    }
    else localExtents[i*2] = localExt[i*2] + this->ExtentOffset[i];
    // max bound
    if (localExt[i*2+1] == WholeExtNoDelta[i*2+1]) {
      localExtents[i*2+1] = WholeExt[i*2+1];
      localbounds[i*2+1]  = this->GlobalResamplingBounds.GetBound(i*2+1);
    }
    else localExtents[i*2+1] = localExt[i*2+1] + this->ExtentOffset[i];
  }
  //
  if (bet) this->ExtentTranslator->ShallowCopy(bet);
  if (this->Delta>0.0) {
    int BETExtents[6];
    this->ExtentTranslator->SetUserBoundsEnabled(1);
    this->ExtentTranslator->SetBoundsForPiece(this->UpdatePiece,localbounds);
    this->ExtentTranslator->ExchangeBoundsForAllProcesses(this->Controller, localbounds);
    this->ExtentTranslator->BoundsToExtentThreadSafe(localbounds, WholeExt, BETExtents);
    //
    std::stringstream temp;
    temp << "Image sampler changing local extents " << this->UpdatePiece << " from {";
    for (int i=0; i<6; i++) temp << localExt[i] << (i<5 ? "," : "} to {");
    for (int i=0; i<6; i++) temp << localExtents[i] << (i<5 ? "," : "} but BET produced {");
    for (int i=0; i<6; i++) temp << BETExtents[i] << (i<5 ? "," : "}");
    vtkDebugMacro( << temp.str() );
  }

  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation  *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet       *input = vtkDataSet::SafeDownCast(this->GetInput());
  vtkImageData  *outImage = this->GetOutput();
  //
  vtkBoundsExtentTranslator *bet = inInfo ? vtkBoundsExtentTranslator::SafeDownCast(
    inInfo->Get(vtkBoundsExtentTranslator::META_DATA())) : NULL;
  //
  if (!this->BoundsInitialized) {
    this->ComputeGlobalInformation(input, this->WholeDimensionWithoutDelta, false);
    this->ComputeGlobalInformation(input, this->WholeDimensionWithDelta, true);
  }
  this->ComputeLocalInformation(bet, this->LocalExtent);
  //
  int outWholeExt[6] = {
    0, this->WholeDimensionWithDelta[0]-1, 
    0, this->WholeDimensionWithDelta[1]-1, 
    0, this->WholeDimensionWithDelta[2]-1 };

  //
  if (!bet) {
    int updatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    int updateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    vtkSmartPointer<vtkExtentTranslator> translator = vtkSmartPointer<vtkExtentTranslator>::New();
    outInfo->Set(vtkBoundsExtentTranslator::META_DATA(), translator);
    translator->SetWholeExtent(outWholeExt);
    translator->SetPiece(updatePiece);
    translator->SetNumberOfPieces(updateNumPieces);
    translator->SetGhostLevel(0);
    translator->PieceToExtent();
    translator->GetExtent(this->LocalExtent);
  }

  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outWholeExt, 6);
  outImage->SetExtent(this->LocalExtent);
  outImage->SetOrigin(this->GlobalResamplingOrigin);
  outImage->SetSpacing(this->spacing);

  //
  // Now resample the input data onto the generated grid
  //
  this->SPHProbe->SetSPHManager(this->SPHManager);
  this->SPHProbe->SetController(this->Controller);
  this->SPHProbe->SetProgressFrequency(100);
  this->SPHProbe->SetDensityScalars(this->DensityScalars);
  this->SPHProbe->SetMassScalars(this->MassScalars);
  this->SPHProbe->SetVolumeScalars(this->VolumeScalars);
  this->SPHProbe->SetHScalars(this->HScalars);
  this->SPHProbe->SetComputeDensityFromNeighbourVolume(this->ComputeDensityFromNeighbourVolume);
  this->SPHProbe->SetAbortLongCalculations(this->GetAbortLongCalculations());
  this->SPHProbe->InitializeVariables(input);
  this->SPHProbe->InitializeKernelCoefficients();
  this->SPHProbe->ProbeMeshless(input, outImage, outImage);

  //
  return 1;
}
