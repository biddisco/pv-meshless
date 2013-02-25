/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkRegularGridSource.cpp,v $
  Date:      $Date: 2002/12/30 22:27:27 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
//
#include "vtkRegularGridSource.h"
//
#include "vtkObjectFactory.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataSet.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
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
vtkStandardNewMacro(vtkRegularGridSource);
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
vtkRegularGridSource::vtkRegularGridSource(void) {
  this->Spacing[0]              = 0.0;
  this->Spacing[1]              = 0.0;
  this->Spacing[2]              = 0.0;
  this->Origin[0]               = 0.0;
  this->Origin[1]               = 0.0;
  this->Origin[2]               = 0.0;
  this->Point1[0]               = 1.0;
  this->Point1[1]               = 0.0;
  this->Point1[2]               = 0.0;
  this->Point2[0]               = 0.0;
  this->Point2[1]               = 1.0;
  this->Point2[2]               = 0.0;
  this->Point3[0]               = 0.0;
  this->Point3[1]               = 0.0;
  this->Point3[2]               = 1.0;
  this->Delta                   = 0.0;
  this->Resolution[0]           = 0;
  this->Resolution[1]           = 0;
  this->Resolution[2]           = 0;
  this->GenerateConnectedCells  = 0;
  this->UseAutoPlacement        = 0;
  this->WholeDimension[0]       = 1;
  this->WholeDimension[1]       = 1;
  this->WholeDimension[2]       = 1;
}
//----------------------------------------------------------------------------
void vtkRegularGridSource::TagDataSet(vtkDataSet *output, char *name)
{
  // tag the dataset with the data so we can use an implicit function 
  // if needed when probing mesh based data
  //
  double norm[3];
  vtkMath::Cross(this->axesvectors[0], this->axesvectors[1], norm);
  double sx = vtkMath::Norm(this->axesvectors[0]);
  double sy = vtkMath::Norm(this->axesvectors[1]);

  // array format is #Label X Y Z Nx Ny Nz X_size Y_size
  vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
  row->SetName(name);
  row->InsertNextValue(vtkVariant("Gui Generated Probe"));
//  row->InsertNextValue(vtkVariant(this->centre[0]));
//  row->InsertNextValue(vtkVariant(this->centre[1]));
//  row->InsertNextValue(vtkVariant(this->centre[2]));
  row->InsertNextValue(vtkVariant(norm[0]));
  row->InsertNextValue(vtkVariant(norm[1]));
  row->InsertNextValue(vtkVariant(norm[2]));
  row->InsertNextValue(vtkVariant(sx));
  row->InsertNextValue(vtkVariant(sy));
  output->GetFieldData()->AddArray(row);
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}
//----------------------------------------------------------------------------
double *vtkRegularGridSource::GetNormal()
{
  double v0v1[3] = {this->Point1[0]-this->Origin[0], 
                    this->Point1[1]-this->Origin[1], 
                    this->Point1[2]-this->Origin[2] };
  double v0v2[3] = {this->Point2[0]-this->Origin[0], 
                    this->Point2[1]-this->Origin[1], 
                    this->Point2[2]-this->Origin[2] };
  vtkMath::Cross(v0v1, v0v2, this->normal);
  vtkMath::Normalize(this->normal);
  return this->normal;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::RequiredDataType()
{
  if (!this->GenerateConnectedCells) return VTK_POLY_DATA;
  //
  return VTK_STRUCTURED_GRID;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::RequestDataObject(
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
      case VTK_POLY_DATA:
        newOutput = vtkPolyData::New();
        break;
      case VTK_STRUCTURED_GRID:
        newOutput = vtkStructuredGrid::New();
        break;
    }
    outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
    newOutput->FastDelete();
    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkRegularGridSource::RequestUpdateExtent(
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
  if (inInfo) {
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  } 
  return 1;
}

//----------------------------------------------------------------------------
int vtkRegularGridSource::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  this->ComputeInformation(request, inputVector, outputVector);
  //
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //
	int maxpieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 
    0, this->WholeDimension[0]-1, 
    0, this->WholeDimension[1]-1, 
    0, this->WholeDimension[2]-1 );
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),maxpieces);
  outInfo->Set(vtkDataObject::ORIGIN(), this->origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->spacing, 3);

  //
  //std::cout << "RI WHOLE_EXTENT {";
  //for (int i=0; i<3; i++) std::cout << WholeDimension[i] << (i<2 ? "," : "}");
  //std::cout << std::endl;

  return 1;
}
//----------------------------------------------------------------------------
int dominantAxis(double vec[3]) {
  if (vec[0]>vec[1] && vec[0]>vec[2]) return 0;
  if (vec[1]>vec[0] && vec[1]>vec[2]) return 1;
  if (vec[2]>vec[0] && vec[2]>vec[1]) return 2;
  return -1;
}
//----------------------------------------------------------------------------
int missingAxis(int axis[3]) {
  std::set<int> present;
  for (int i=0; i<3; i++) present.insert(axis[i]);
  for (int i=0; i<3; i++) if (present.find(i)==present.end()) return i;
  return -1;
}
//----------------------------------------------------------------------------
// @TODO Check this, delta not included
void vtkRegularGridSource::ComputeAxesFromPoints(double lengths[3], bool inflate) 
{
  double tempaxes[3][3], templengths[3];
  int axes[3];
  for (int i=0; i<3; i++) { 
    this->origin[i]  = this->Origin[i];
    tempaxes[0][i] = this->Point1[i] - this->Origin[i];
    tempaxes[1][i] = this->Point2[i] - this->Origin[i];
    tempaxes[2][i] = this->Point3[i] - this->Origin[i];
  }
  // we sort the axes to be as close as we can to {X,Y,Z} because the Extents 
  // are always ordered this way. (Not valid when non axis-aligned).
  for (int i=0; i<3; i++) {
    axes[i] = dominantAxis(tempaxes[i]);
    templengths[i] = sqrt(vtkMath::Norm(tempaxes[i]));
  }
  for (int i=0; i<3; i++) {
    if (axes[i]!=-1) {
      this->axesvectors[axes[i]][0] = tempaxes[i][0];
      this->axesvectors[axes[i]][1] = tempaxes[i][1];
      this->axesvectors[axes[i]][2] = tempaxes[i][2];
      lengths[axes[i]] = templengths[i];
    }
    else {
      int missing = missingAxis(axes);
      this->axesvectors[missing][0] = 0.0;
      this->axesvectors[missing][1] = 0.0;
      this->axesvectors[missing][2] = 0.0;
      lengths[missing] = 0.0;
    }
  }
}
//----------------------------------------------------------------------------
void vtkRegularGridSource::ComputeAxesFromBounds(vtkDataSet *inputData, double lengths[3], bool inflate)
{
  //
  // Define box...
  //
  vtkBoundingBox box;
  double bounds[6];
  inputData->GetBounds(bounds);
  box.SetBounds(bounds);
      
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
  box.GetMinPoint(this->origin[0], this->origin[1], this->origin[2]);
  box.GetLengths(lengths);
  //
  this->axesvectors[0][0] = this->axesvectors[1][0] = this->axesvectors[2][0] = 0.0;
  this->axesvectors[0][1] = this->axesvectors[1][1] = this->axesvectors[2][1] = 0.0;
  this->axesvectors[0][2] = this->axesvectors[1][2] = this->axesvectors[2][2] = 0.0;
  this->axesvectors[0][0] += lengths[0];
  this->axesvectors[1][1] += lengths[1];
  this->axesvectors[2][2] += lengths[2];
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::ComputeInformation(
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
  if (inputData && this->UseAutoPlacement) {
    this->ComputeAxesFromBounds(inputData, lengths, true);
  }
  else {
    this->ComputeAxesFromPoints(lengths, true);
  }

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
void vtkRegularGridSource::BoundsToExtent(double *bounds, int *extent, int updatePiece) 
{
  double axesvec[3] = {this->axesvectors[0][0], this->axesvectors[1][1], this->axesvectors[2][2]};
  vtkBoundingBox processRegion(bounds);
  vtkBoundingBox dataRegion;
  dataRegion.SetMinPoint(this->origin);
  dataRegion.SetMaxPoint(this->origin[0]+axesvec[0],
                         this->origin[1]+axesvec[1],
                         this->origin[2]+axesvec[2]);
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
  REGULARGRID_ISAXPY(updateExtentLo, this->scaling, axesvec, this->origin, minvec);
  REGULARGRID_ISAXPY(updateExtentHi, this->scaling, axesvec, this->origin, maxvec);
  for (int i=0; i<3; i++) { extent[i*2  ] = static_cast<int>(updateExtentLo[i]+0.5); }
  for (int i=0; i<3; i++) { extent[i*2+1] = static_cast<int>(updateExtentHi[i]+0.5); }
  //
  std::cout << updatePiece << "Setting Extent to {";
  for (int i=0; i<6; i++) std::cout << extent[i] << (i<5 ? "," : "}");
  std::cout << std::endl;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkDataSet         *output = NULL;
  vtkPolyData       *outPoly = NULL;
  vtkStructuredGrid *outGrid = NULL;
  //
  switch (this->RequiredDataType()) {
    case VTK_POLY_DATA:
      outPoly = this->GetPolyDataOutput();
      output = outPoly;
      break;
    case VTK_STRUCTURED_GRID:
      outGrid = this->GetStructuredGridOutput();
      output = outGrid;
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
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  vtkBoundsExtentTranslator *bet = vtkBoundsExtentTranslator::SafeDownCast(translator);
  if (bet) {
    int updatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    double *bounds = bet->GetBoundsForPiece(updatePiece);
//    this->BoundsToExtent(bounds,outUpdateExt,updatePiece);
    bet->BoundsToExtentThreadSafe(bounds, outWholeExt,outUpdateExt); 
  }
  else {
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outUpdateExt);
  }
  //
  int dims[3]= {1+outUpdateExt[1]-outUpdateExt[0],
                1+outUpdateExt[3]-outUpdateExt[2], 
                1+outUpdateExt[5]-outUpdateExt[4]};
  vtkIdType NumPoints = dims[0]*dims[1]*dims[2];
  //

  vtkSmartPointer<vtkPoints>     newpoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> dataarray = vtkFloatArray::SafeDownCast(newpoints->GetData());
  float *pointdata = dataarray->WritePointer(0,NumPoints*3);
  //
  vtkIdType totalpoints = 0;
  double pos1[3], pos2[3];

  //
  // Walk bounds and create points
  //
  for (double k=outUpdateExt[4]; k<=outUpdateExt[5]; k++) {
    // increment along Z axis by k
    REGULARGRID_SAXPY(k*this->scaling[2], this->axesvectors[2], this->origin, pos2);
    for (double j=outUpdateExt[2]; j<=outUpdateExt[3]; j++) {
      // increment along Y axis by j
      REGULARGRID_SAXPY(j*this->scaling[1], this->axesvectors[1], pos2, pos1);
      for (double i=outUpdateExt[0]; i<=outUpdateExt[1]; i++) {
        // write final value directly to the points array
        float *position = &pointdata[totalpoints*3];
        // increment start point along Y axis by j
        REGULARGRID_SAXPY(i*this->scaling[0], this->axesvectors[0], pos1, position);
        totalpoints++;
      }
    }
  }
  newpoints->SetNumberOfPoints(totalpoints);
  //
  if (outPoly) {
    vtkCellArray *outVerts   = vtkCellArray::New();
    vtkIdType *arraydata = outVerts->WritePointer(totalpoints, 2*totalpoints);
    //
    for (int i=0; i<totalpoints; i++) {
      arraydata[i*2]   = 1;
      arraydata[i*2+1] = i;
    }
    outPoly->SetVerts(outVerts);
    outVerts->FastDelete();
    outPoly->SetPoints(newpoints);
  }
  else if (outGrid) {
    outGrid->SetExtent(outUpdateExt);
    outGrid->SetPoints(newpoints);
  }
  if (this->WholeDimension[2]==1) {
    this->TagDataSet(output, "ImplicitPlaneData");
    //
    vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(totalpoints);
    normals->SetName("Normals");
    double *normal = this->GetNormal();
    for (int i=0; i<totalpoints; i++) {
      normals->SetTuple(i, normal);
    }
    output->GetPointData()->SetNormals(normals);
  }
  return 1;
}
//----------------------------------------------------------------------------
void vtkRegularGridSource::ComputeOBB()
{
  vtkPointSet *input = vtkPointSet::SafeDownCast(this->GetInput());
  if (!input) return;
  vtkPoints *pts = input->GetPoints();
  //
  double corner[3], max[3], mid[3], min[3], size[3];
  vtkSmartPointer<vtkOBBTree> obb = vtkSmartPointer<vtkOBBTree>::New();
  obb->ComputeOBB(pts, corner, max, mid, min, size);
  //
  for (int i=0; i<3; i++) {
    this->Origin[i] = corner[i];
    this->Point1[i] = max[i]+corner[i];
    this->Point2[i] = mid[i]+corner[i];
    this->Point3[i] = min[i]+corner[i];
  }
  //
  // the OBB returned by vtkOBBTree may have the z axis defined in the opposite
  // sense to what we are accustomed to (left/right handed axes)
  //
  double vz[3] = {
    this->Point3[0]-this->Origin[0], 
    this->Point3[1]-this->Origin[1], 
    this->Point3[2]-this->Origin[2]};

  REGULARGRID_SAXPY( 1.0, vz, this->Origin, this->Origin);
  REGULARGRID_SAXPY( 1.0, vz, this->Point1, this->Point1);
  REGULARGRID_SAXPY( 1.0, vz, this->Point2, this->Point2);
  REGULARGRID_SAXPY(-1.0, vz, this->Point3, this->Point3);

  this->Modified();
}
//----------------------------------------------------------------------------
