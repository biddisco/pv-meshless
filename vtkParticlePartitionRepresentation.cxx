/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkParticlePartitionRepresentation.cpp
  Revision of last commit : $Rev: 155 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2006-07-13 10:23:31 +0200 #$

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
#include "vtkParticlePartitionRepresentation.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkBoundingBox.h"
//
#include "vtkBoundsExtentTranslator.h"
#include "vtkAppendPolyData.h"
#include "vtkOutlineSource.h"
//
#include <cmath>
//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkParticlePartitionRepresentation, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkParticlePartitionRepresentation);
//---------------------------------------------------------------------------
vtkParticlePartitionRepresentation::vtkParticlePartitionRepresentation(void)
{
  this->AllBoxesOnAllProcesses = 0;
  this->InflateFactor = 1.0;
}
//---------------------------------------------------------------------------
vtkParticlePartitionRepresentation::~vtkParticlePartitionRepresentation(void) {
}
//----------------------------------------------------------------------------
int vtkParticlePartitionRepresentation::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionRepresentation::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionRepresentation::RequestData(vtkInformation *request,
                                       vtkInformationVector** inputVector,
                                       vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  vtkPointSet *input  = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  vtkBoundsExtentTranslator *bet = vtkBoundsExtentTranslator::SafeDownCast(translator);
  //
  vtkSmartPointer<vtkAppendPolyData> polys = vtkSmartPointer<vtkAppendPolyData>::New();
  size_t boxes = 1;
  if (bet) {
    boxes = bet->GetNumberOfPieces();  
  }
  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  vtkSmartPointer<vtkIntArray> processIds = vtkSmartPointer<vtkIntArray>::New();
  processIds->SetName("ProcessId");
  //
  double bounds[6];
  vtkBoundingBox box;
  for (int i=0; i<boxes; i++)
  {
    bool add = false;
    if (this->AllBoxesOnAllProcesses && bet) {
      box.SetBounds(bet->GetBoundsForPiece(i));
      add = true;
    }
    else if (i==piece) {
      if (bet) box.SetBounds(bet->GetBoundsForPiece(i));
      else box.SetBounds(input->GetBounds());
      add = true;
    }
    if (add) {
      double p1[3],p2[3];
      box.GetMaxPoint(p1[0], p1[1], p1[2]);
      box.GetMinPoint(p2[0], p2[1], p2[2]);
      for (int j=0; j<3; j++) {
        double l = box.GetLength(j);
        double d = (l-l*this->InflateFactor);
        p1[j] -= d/2.0; 
        p2[j] += d/2.0; 
      }
      vtkSmartPointer<vtkOutlineSource> cube = vtkSmartPointer<vtkOutlineSource>::New();
      cube->SetBounds(p1[0],p2[0],p1[1],p2[1],p1[2],p2[2]);
      cube->Update();
      polys->AddInput(cube->GetOutput());
      processIds->InsertNextValue(i);
    }
  }
  polys->Update();
  output->SetPoints(polys->GetOutput()->GetPoints());
  output->SetLines(polys->GetOutput()->GetLines());
  output->GetPointData()->AddArray(processIds);

  return 1;
}
//---------------------------------------------------------------------------


