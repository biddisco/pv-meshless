/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiBlockExtractDataOverTime.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMultiBlockExtractDataOverTime.h"

#include "vtkSmartPointer.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(vtkMultiBlockExtractDataOverTime);

//----------------------------------------------------------------------------
vtkMultiBlockExtractDataOverTime::vtkMultiBlockExtractDataOverTime()
{
  this->NumberOfTimeSteps = 0;
  this->CurrentTimeIndex = 0;
  this->PointIndex = 0; 
}

//----------------------------------------------------------------------------
void vtkMultiBlockExtractDataOverTime::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Point Index: " << this->PointIndex << endl;
  os << indent << "NumberOfTimeSteps: " << this->NumberOfTimeSteps << endl;
}


//----------------------------------------------------------------------------
int vtkMultiBlockExtractDataOverTime::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  if ( inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) )
    {
    this->NumberOfTimeSteps = 
      inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
    }
  else
    {
    this->NumberOfTimeSteps = 0;
    }
  // The output of this filter does not contain a specific time, rather 
  // it contains a collection of time steps. Also, this filter does not
  // respond to time requests. Therefore, we remove all time information
  // from the output.
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_RANGE()))
    {
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
    }

  return 1;
}


//----------------------------------------------------------------------------
int vtkMultiBlockExtractDataOverTime::ProcessRequest(vtkInformation* request,
                                           vtkInformationVector** inputVector,
                                           vtkInformationVector* outputVector)
{
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return this->RequestInformation(request, inputVector, outputVector);
    }
  else if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    // get the requested update extent
    double *inTimes = inputVector[0]->GetInformationObject(0)
      ->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    if (inTimes)
      {
      double timeReq[1];
      timeReq[0] = inTimes[this->CurrentTimeIndex];
      inputVector[0]->GetInformationObject(0)->Set
        ( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
      }
    return 1;
    }
  
  // generate the data
  else if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
  {
    if (this->NumberOfTimeSteps == 0)
      {
      vtkErrorMacro("No Time steps in input time data!");
      return 0;
      }

    // get the output data object
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *output = 
      vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // and input data object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet *input = 
      vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // is this the first request
    if (!this->CurrentTimeIndex)
      {
      // Tell the pipeline to start looping.
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      this->AllocateOutputData(input, output);
      }

    int blocknum = 0;
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(input->NewIterator());
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
      vtkPointSet* indata = vtkPointSet::SafeDownCast(iter->GetCurrentDataObject());
      vtkPointSet* outdata = vtkPointSet::SafeDownCast(output->GetBlock(blocknum));
      if (indata && outdata) {
        // extract the actual data
        outdata->GetPoints()->SetPoint( this->CurrentTimeIndex, indata->GetPoints()->GetPoint(this->PointIndex) );
        outdata->GetPointData()->CopyData(indata->GetPointData(), this->PointIndex, this->CurrentTimeIndex);
        if (indata->GetPointData()->GetArray("Time")) {
          outdata->GetPointData()->GetArray("TimeData")->SetTuple1
            (this->CurrentTimeIndex, indata->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP()));
        }
        else {
          outdata->GetPointData()->GetArray("Time")->SetTuple1
            (this->CurrentTimeIndex, indata->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP()));
        }
      }
      blocknum++;
    }

    // increment the time index
    this->CurrentTimeIndex++;
    if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
    {
      // Tell the pipeline to stop looping.
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
    }  
    return 1;
  }
  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}
//----------------------------------------------------------------------------
int vtkMultiBlockExtractDataOverTime::AllocateOutputData(vtkMultiBlockDataSet *input, vtkMultiBlockDataSet *output)
{
  // by default vtkMultiBlockDataSet::RequestDataObject already 
  // created an output of the same type as the input
  if (!output)
    {
    vtkErrorMacro("Output not created as expected!");
    return 0;
    }

  int blocknum = 0;
  vtkSmartPointer<vtkCompositeDataIterator> iter;
  iter.TakeReference(input->NewIterator());
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
  {
    vtkPointSet* dObj = vtkPointSet::SafeDownCast(iter->GetCurrentDataObject());
    if (dObj) {
      // create an output
      vtkSmartPointer<vtkPolyData> out = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      out->SetPoints( points );
      points->SetNumberOfPoints( this->NumberOfTimeSteps );
      // init the point data
      out->GetPointData()->CopyAllocate(dObj->GetPointData(), this->NumberOfTimeSteps);
      // and add an array to hold the time at each step
      vtkSmartPointer<vtkDoubleArray> timeArray = vtkSmartPointer<vtkDoubleArray>::New();
      timeArray->SetNumberOfComponents(1);
      timeArray->SetNumberOfTuples(this->NumberOfTimeSteps);
      if (dObj->GetPointData()->GetArray("Time")) {
        timeArray->SetName("TimeData");
      }
      else {
        timeArray->SetName("Time");
      }
      out->GetPointData()->AddArray(timeArray);
      output->SetBlock(blocknum,out);

      const char *dname = input->GetMetaData(blocknum)->Get(vtkCompositeDataSet::NAME());
      if (dname) {
        output->GetMetaData(blocknum)->Set(vtkCompositeDataSet::NAME(), dname);
      }
    }
    else {
      output->SetBlock(blocknum,NULL);
    }
    blocknum++;
  }
  return 1;
}


