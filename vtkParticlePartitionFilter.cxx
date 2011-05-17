/*=========================================================================

  Project                 : vtkCSCS
  Module                  : vtkParticlePartitionFilter.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
#include "vtkParticlePartitionFilter.h"
#include "vtkPolyData.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkMultiProcessController.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticlePartitionFilter);
vtkCxxSetObjectMacro(vtkParticlePartitionFilter,Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::vtkParticlePartitionFilter()
{
  this->UpdatePiece         = 0;
  this->UpdateNumPieces     = 0;
  this->NumberOfLocalPoints = 0;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkParticlePartitionFilter::~vtkParticlePartitionFilter()
{
   this->SetController(0);
}

//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);

  // Setup the output
  vtkPolyData* output = vtkPolyData::GetData(outputVector);

  //
  vtkSmartPointer<vtkPoints>   newPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  output->SetPoints(newPoints);
  output->SetVerts(vertices);

/*
  // Also saving it as a data array for easy csv export
  vtkSmartPointer<vtkFloatArray> c_of_m = vtkSmartPointer<vtkFloatArray>::New();
  c_of_m->SetName("CentreOfMass");
  c_of_m->SetNumberOfComponents(3);
  c_of_m->SetNumberOfTuples(1);
  c_of_m->SetTuple(0, dbCenterOfMass);
  output->GetPointData()->AddArray(c_of_m);
*/
  //
  timer->StopTimer();
  if (this->UpdatePiece==0) { 
//    vtkErrorMacro(<< "Particle partitioning : " << timer->GetElapsedTime() << " seconds\n");
  }
  return 1;
}
