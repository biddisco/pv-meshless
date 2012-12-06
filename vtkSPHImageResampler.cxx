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
#include "vtkBoundingBox.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"
//
#include "vtkPVExtentTranslator.h"
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
//----------------------------------------------------------------------------
#if 0
  #define OUTPUTTEXT(a) std::cout << (a);

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << /* this->UpdatePiece << " : " */ a << endl; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }
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
vtkSPHImageResampler::vtkSPHImageResampler(void) {
  this->GlobalOrigin[0]         = 0.0;
  this->GlobalOrigin[1]         = 0.0;
  this->GlobalOrigin[2]         = 0.0;
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
  vtkDebugMacro( "Imagesampler (" << piece << ") set num pieces to " << numPieces );
  //  
  return 1;
}

//----------------------------------------------------------------------------
int vtkSPHImageResampler::RequestInformation(
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
  //vtkDebugMacro( "RI WHOLE_EXTENT {";
  //for (int i=0; i<3; i++) std::cout << WholeDimension[i] << (i<2 ? "," : "}");
  //std::cout << std::endl;

  return 1;
}
//----------------------------------------------------------------------------
void vtkSPHImageResampler::ComputeAxesFromBounds(vtkDataSet *inputData, double lengths[3], bool inflate)
{
  //
  // Define box...
  //
  vtkBoundingBox box;
  double bounds[6];
  inputData->GetBounds(bounds);
  this->BoundsInitialized = vtkMath::AreBoundsInitialized(bounds);
  box.SetBounds(bounds);
  //
  double bmin[3], bmn[3] = {bounds[0], bounds[2], bounds[4]};
  double bmax[3], bmx[3] = {bounds[1], bounds[3], bounds[5]};
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
int vtkSPHImageResampler::ComputeInformation(
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
  vtkDebugMacro( "vtkSPHImageResampler::ComputeInformation BoundsInitialized " << BoundsInitialized );

  if (this->BoundsInitialized) {
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
  }
  else {
    for (int i=0; i<3; i++) {
      this->WholeDimension[i] = this->Resolution[i];
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHImageResampler::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkImageData *outImage = this->GetOutput();
  vtkDataSet      *input = vtkDataSet::SafeDownCast(this->GetInput());
  //
  if (!this->BoundsInitialized) {
    this->ComputeInformation(request, inputVector, outputVector);
  }
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
    bet->BoundsToExtentThreadSafe(bounds, outWholeExt, outUpdateExt); 
    std::stringstream temp;
    temp << "Image sampler " << updatePiece << " Setting Extent to {";
    for (int i=0; i<6; i++) temp << outUpdateExt[i] << (i<5 ? "," : "}");
//    vtkDebugMacro( temp.str() );
  }
  else {
//    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outUpdateExt);
    int updatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    int updateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    vtkSmartPointer<vtkPVExtentTranslator> translator = vtkSmartPointer<vtkPVExtentTranslator>::New();
    outInfo->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), translator);
    translator->SetWholeExtent(0,this->WholeDimension[0]-1,0,WholeDimension[1]-1,0,WholeDimension[2]-1);
    translator->SetPiece(updatePiece);
    translator->SetNumberOfPieces(updateNumPieces);
    translator->SetGhostLevel(0);
    translator->PieceToExtent();
    translator->GetExtent(outUpdateExt);
  }
  //
  int dims[3]= {1+outUpdateExt[1]-outUpdateExt[0],
                1+outUpdateExt[3]-outUpdateExt[2], 
                1+outUpdateExt[5]-outUpdateExt[4]};
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outWholeExt, 6);
//  outImage->SetWholeExtent(outWholeExt);
  outImage->SetExtent(outUpdateExt);
  outImage->SetOrigin(this->GlobalOrigin);
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
  this->SPHProbe->InitializeVariables(input);
  this->SPHProbe->InitializeKernelCoefficients();
  this->SPHProbe->ProbeMeshless(input, outImage, outImage);
  //
  return 1;
}
