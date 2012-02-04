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
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
//
#include "vtkBoundsExtentTranslator.h"
#include "vtkAppendPolyData.h"
#include "vtkCubeSource.h"
//
#include <cmath>
//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkParticlePartitionRepresentation, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkParticlePartitionRepresentation);
//---------------------------------------------------------------------------
vtkParticlePartitionRepresentation::vtkParticlePartitionRepresentation(void)
{
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
  double *bounds = NULL;
  if (bet) {
    int updatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    bounds = bet->GetBoundsForPiece(updatePiece);
  }
  else {
    bounds = input->GetBounds();
  }
  //
  // Ok, now create cube(oid)s and stuff'em into a polydata thingy
  vtkAppendPolyData *polys = vtkAppendPolyData::New();
  size_t s = 1;
  for (size_t i=0; i<s; i++)
    {
    vtkCubeSource *cube = vtkCubeSource::New();
    cube->SetBounds( bounds );
    cube->Update();
    polys->AddInput(cube->GetOutput());
    cube->Delete();
    }
  polys->Update();
  output->SetPoints(polys->GetOutput()->GetPoints());
  output->SetPolys(polys->GetOutput()->GetPolys());
  polys->Delete();

  return 1;
}
//---------------------------------------------------------------------------


