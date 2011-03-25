/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSPHProbeFilter3.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSPHProbeFilter3.h"
#include "vtkSPHManager.h"

#include "vtkSmartPointer.h"
#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkCharArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointLocator.h"
#include "vtkPointSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkPointSet.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkTimerLog.h"
#include "vtkVariantArray.h"
#include "vtkImplicitFunction.h"
#include "vtkCutter.h"
#include "vtkExtractGeometry.h"
#include "vtkProbeFilter.h"
#include "vtkDataObjectTypes.h"
#include "vtkBoundingBox.h"
#include "vtkSamplingGridGenerator.h"
#include "vtkPlane.h"
//
#include "KernelGaussian.h"
#include "KernelWendland.h"
#include "KernelQuadratic.h"
#include "KernelSpline3rdOrder.h"
#include "KernelSpline5thOrder.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSPHProbeFilter3);
vtkCxxSetObjectMacro(vtkSPHProbeFilter3, CutFunction, vtkImplicitFunction);
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkSPH3_Copy(vtkDataSet *d) {
  vtkSmartPointer<vtkDataSet> result;
  result.TakeReference(d->NewInstance());
  result->ShallowCopy(d);
  return result;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkSPHProbeFilter3::vtkSPHProbeFilter3()
{
  this->CutFunction = NULL;
  this->Resolution[0] = 0;
  this->Resolution[1] = 0;
  this->Resolution[2] = 0;
  this->GenerateConnectedCells = 1;
  this->UseKernelDistanceForSampling = 0;
  this->KernelDistanceMultiplier = 1.0;
  //
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkSPHProbeFilter3::~vtkSPHProbeFilter3()
{
}
//----------------------------------------------------------------------------
unsigned long vtkSPHProbeFilter3::GetMTime()
{
  unsigned long mTime=this->Superclass::GetMTime();
  unsigned long time;

  if ( this->CutFunction != NULL )
    {
    time = this->CutFunction->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }
  if ( this->Locator != NULL )
    {
    time = this->Locator->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  return mTime;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> vtkSPHProbeFilter3::GenerateProbePts(vtkDataSet *data)
{
  vtkSmartPointer<vtkSamplingGridGenerator> grid = vtkSmartPointer<vtkSamplingGridGenerator>::New();
  // always copy data when adding to internal filters to prevent spurious updates
  vtkSmartPointer<vtkDataSet> copy = vtkSPH3_Copy(data);
  grid->SetInput(copy);
  grid->SetCutFunction(this->GetCutFunction());    
  //
  double delta = 0.0;
  if (this->UseKernelDistanceForSampling && this->KernelFunction) {
    delta = this->GetMaxKernelCutoffDistance();
    double md = this->KernelDistanceMultiplier*delta/2.0;
    grid->SetSpacing(md, md, md);
    grid->SetResolution(0, 0, 0);
  }
  else {
    if (this->GetCutFunction() && vtkPlane::SafeDownCast(this->GetCutFunction())) {
      grid->SetResolution(this->Resolution[0], this->Resolution[1], 1);
    }
    else {
      grid->SetResolution(this->Resolution[0], this->Resolution[1], this->Resolution[2]);
    }
  }
  grid->SetGenerateConnectedCells(this->GenerateConnectedCells);
  grid->SetDelta(delta);
  grid->Update();
  return vtkPointSet::SafeDownCast(grid->GetOutput());
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter3::RequiredDataType()
{
  if (!this->GenerateConnectedCells) return VTK_POLY_DATA;
  //
  if (this->GetCutFunction() && vtkPlane::SafeDownCast(this->GetCutFunction())) {
//    return VTK_POLY_DATA;
  }
  return VTK_STRUCTURED_GRID;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter3::RequestDataObject(
  vtkInformation *, 
  vtkInformationVector  **vtkNotUsed(inputVector), 
  vtkInformationVector *outputVector)
{
  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));
  bool ok = (output!=NULL);
  //
  ok = (ok && output->GetDataObjectType()==this->RequiredDataType());
  //
  vtkDataSet *newOutput = NULL;
  if (!ok) {
    switch (this->RequiredDataType()) {
      case VTK_POLY_DATA:
        newOutput = vtkPolyData::New();
        break;
      case VTK_STRUCTURED_GRID:
        newOutput = vtkStructuredGrid::New();
        break;
    }
    newOutput->SetPipelineInformation(info);
    newOutput->Delete();
    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter3::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);

  // get the input and output, check for multiblock probe/output geometry
  vtkPointSet *particles = vtkPointSet::SafeDownCast(
    dataInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  //
  // Make sure scalar array pointers are setup
  //
  this->InitializeVariables(particles);

  //
  // Precompute any kernel coefficients that we already know 
  // * * before * * generating probe pts as we need Kernel for delta
  //
  if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
    this->InitializeKernelCoefficients();
  }
  //
  vtkSmartPointer<vtkPointSet> probepts = this->GenerateProbePts(particles);

  //
  this->UpdateProgress(0.0);
  //
  this->Timer->StartTimer();
  if (probepts) {
    this->ProbeMeshless(particles, probepts, output);
    vtkStructuredGrid *grid = vtkStructuredGrid::SafeDownCast(probepts);
    if (grid) {
      int dims[3];
      grid->GetDimensions(dims);
      output->SetWholeExtent(0, dims[0]-1, 0,dims[1]-1, 0,dims[2]-1);
    }
  }
  this->Timer->StopTimer();
  double time = this->Timer->GetElapsedTime();
//  vtkErrorMacro(<<"Probe time is " << time);
  this->UpdateProgress(1.0);
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter3::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);
/*
  I don't remember why I wanted this in, but it stops the probe re-execting
  when time increments - so leave it out
  outInfo->CopyEntry(dataInfo, 
                     vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->CopyEntry(dataInfo, 
                     vtkStreamingDemandDrivenPipeline::TIME_RANGE());
*/
  // be careful, our output extent is taken from the probe pts which is the 
  // second input
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
//               probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
//               6);
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
//               probePtsInfo->Get(
//                 vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()));

  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter3::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // we must request the whole extent from probe pts
//  probePtsInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
//              probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
//              6);
    
  return 1;
}
