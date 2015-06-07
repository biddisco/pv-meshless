/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkParticleBoxTreeRepresentation.cpp
  Revision of last commit : $Rev: 155 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2006-07-13 10:23:31 +0200 #$

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This code is derived from an earlier work and is distributed
  with permission from, and thanks to

  ------------------------------------------
  Copyright (C) 2000-2004 John Biddiscombe
  Skipping Mouse Software Ltd,
  Blewbury, England
  ------------------------------------------

=========================================================================*/
#include "vtkParticleBoxTreeRepresentation.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkParticleBoxTree.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include <cmath>
//---------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticleBoxTreeRepresentation);
//---------------------------------------------------------------------------
vtkParticleBoxTreeRepresentation::vtkParticleBoxTreeRepresentation(void)
{
  this->Level           = -1;
  this->ParticleSize    = 0.001;
  this->MaxDepth        = 32;
  this->MaxCellsPerNode = 32;
  this->BSPTree = vtkSmartPointer<vtkParticleBoxTree>::New();
}
//---------------------------------------------------------------------------
vtkParticleBoxTreeRepresentation::~vtkParticleBoxTreeRepresentation(void) {
}
//----------------------------------------------------------------------------
int vtkParticleBoxTreeRepresentation::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticleBoxTreeRepresentation::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticleBoxTreeRepresentation::RequestData(vtkInformation *request,
                                       vtkInformationVector** inputVector,
                                       vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  vtkPointSet *input  = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData  *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  if (input->GetNumberOfPoints()==0 || input->GetPointData()==NULL) return 0;
  //
  BSPTree->SetDataSet(input);
  BSPTree->SetParticleSize(this->ParticleSize);
  BSPTree->SetNumberOfCellsPerNode(this->MaxCellsPerNode);
  BSPTree->SetMaxLevel(this->MaxDepth);
  BSPTree->SetCacheCellBounds(1);
  BSPTree->BuildLocator();
  //
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  output->SetPoints(pts);
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  output->SetLines(lines);
  vtkSmartPointer<vtkIntArray> levels = vtkSmartPointer<vtkIntArray>::New();
  levels->SetName("Level");
  output->GetPointData()->AddArray(levels);
  //
  BSPTree->GenerateRepresentation(this->Level,output);
  return 1;
}
//---------------------------------------------------------------------------


