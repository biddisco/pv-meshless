/*=========================================================================

  Project:   vtkCSCS
  Module:    vtkTemporalPlotValueFilter.cxx

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
#include "vtkTemporalPlotValueFilter.h"
#include "vtkPolyData.h"
#include "vtkPointSet.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
//
#include <vector>
#include <list>
#include <map>
#include <string>
#include <stdexcept>
#include <cmath>
//---------------------------------------------------------------------------
vtkStandardNewMacro(vtkTemporalPlotValueFilter);
//----------------------------------------------------------------------------
vtkTemporalPlotValueFilter::vtkTemporalPlotValueFilter()
{
  this->NumberOfTimeSteps    = 0;
  this->LatestTime           = 01E10;
  this->Vertices             = vtkSmartPointer<vtkCellArray>::New();
  this->Points               = vtkSmartPointer<vtkPoints>::New();
  this->Values               = vtkSmartPointer<vtkPointData>::New();
  this->TimeData             = vtkSmartPointer<vtkDoubleArray>::New();
  this->StorageSize          = 1000;
  this->TimeData->SetName("TimeData");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1); // Lines and points
}
//----------------------------------------------------------------------------
vtkTemporalPlotValueFilter::~vtkTemporalPlotValueFilter()
{
}
//----------------------------------------------------------------------------
int vtkTemporalPlotValueFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkTemporalPlotValueFilter::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//----------------------------------------------------------------------------
int vtkTemporalPlotValueFilter::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  outInfo0->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  //vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  //if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
  //  this->NumberOfTimeSteps = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  //}
  return 1;
}
//---------------------------------------------------------------------------
int vtkTemporalPlotValueFilter::RequestData(
  vtkInformation *vtkNotUsed(information),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkDataSet* input = vtkDataSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);

  vtkInformation* inputInfo = input? input->GetInformation() : 0;
  vtkInformation* outputInfo = outputVector->GetInformationObject(0);

  double currenttime = 0;
  if (inputInfo && inputInfo->Has(vtkDataObject::DATA_TIME_STEP()))
    {
    currenttime = inputInfo->Get(vtkDataObject::DATA_TIME_STEP());
    }
  else if (outputInfo && outputInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    currenttime = outputInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    }
  //
  // if the user rewound the animation, the time will be wrong, reset
  //
  if (currenttime<this->LatestTime) {
    this->Flush();
  }
  this->LatestTime = currenttime;
  //
  // Add the latest time to the time values array
  //
  TimeData->InsertNextTuple1(currenttime);
  vtkIdType numT = TimeData->GetNumberOfTuples();

  //
  // At start we must initialize the field arrays
  //
  if (numT==1) {
    Values->CopyAllocate(input->GetPointData());  
  }

  // copy point field data from input
  Values->CopyData(input->GetPointData(), 0, numT-1);
  //
  // create a dummy 3D 'point' to display data in normal plot mode
  //
  double pt[3] = { currenttime, 0.0, 0.0 };
  vtkIdType Id = Points->InsertNextPoint(pt);
  Vertices->InsertNextCell(1,&Id);

  // to stop annoying can't plot with one point error messages
  // we will always do the first point twice
  if (numT==1) {
    TimeData->InsertNextTuple1(currenttime);
    Values->CopyData(input->GetPointData(), 0, numT);
    Vertices->InsertNextCell(1,&Id);
  }

  if (this->StorageSize>0 && this->Values->GetNumberOfTuples()>this->StorageSize) {
    for (int i=0; i<this->StorageSize; i++) {
      this->Values->CopyData(this->Values, i+1, i);
      this->TimeData->SetValue(i, this->TimeData->GetValue(i+1));
    }
    for (int i=0; i<this->Values->GetNumberOfArrays(); i++) {
      this->Values->GetArray(i)->SetNumberOfTuples(this->StorageSize);
    }
    this->TimeData->SetNumberOfTuples(this->StorageSize);
  }
  //
  output->GetPointData()->ShallowCopy(Values);
  output->GetPointData()->AddArray(TimeData);
  return 1;
}
//---------------------------------------------------------------------------
void vtkTemporalPlotValueFilter::Flush()
{
  this->Vertices->Initialize();
  this->Points->Initialize();
  this->Values->Initialize();
  this->TimeData->Initialize();
}
//---------------------------------------------------------------------------
void vtkTemporalPlotValueFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//-----------------------------------------------------------------------------
