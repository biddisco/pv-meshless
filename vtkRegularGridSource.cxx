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
#include "vtkRegularGridSource.h"
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
//
vtkCxxRevisionMacro(vtkRegularGridSource, "$Revision: 1.4 $");
vtkStandardNewMacro(vtkRegularGridSource);
//----------------------------------------------------------------------------
#define REGULARGRID_SAXPY(a,x,y,z) \
  z[0] = a*x[0] + y[0]; \
  z[1] = a*x[1] + y[1]; \
  z[2] = a*x[2] + y[2]; 
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
  this->Dimension[0]            = 1;
  this->Dimension[1]            = 1;
  this->Dimension[2]            = 1;
}
//----------------------------------------------------------------------------
void vtkRegularGridSource::TagDataSet(vtkDataSet *output, char *name)
{
  // tag the dataset with the data so we can use an implicit function 
  // if needed when probing mesh based data
  //
  double norm[3];
  vtkMath::Cross(axesvectors[0], axesvectors[1], norm);
  double sx = vtkMath::Norm(axesvectors[0]);
  double sy = vtkMath::Norm(axesvectors[1]);

  // array format is #Label X Y Z Nx Ny Nz X_size Y_size
  vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
  row->SetName(name);
  row->InsertNextValue(vtkVariant("Gui Generated Probe"));
  row->InsertNextValue(vtkVariant(this->centre[0]));
  row->InsertNextValue(vtkVariant(this->centre[1]));
  row->InsertNextValue(vtkVariant(this->centre[2]));
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
void vtkRegularGridSource::GetNormal(double &n1, double &n2, double &n3)
{
  this->GetNormal();
  n1 = this->normal[0];
  n2 = this->normal[1];
  n3 = this->normal[2];
}
//----------------------------------------------------------------------------
void vtkRegularGridSource::GetNormal(double n[3])
{
  this->GetNormal();
  n[0] = this->normal[0];
  n[1] = this->normal[1];
  n[2] = this->normal[2];
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::RequiredDataType()
{
  if (!this->GenerateConnectedCells) return VTK_POLY_DATA;
  //
  if (this->Dimension[2]==1) {
//    return VTK_POLY_DATA;
  }
  return VTK_STRUCTURED_GRID;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::RequestDataObject(
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
int vtkRegularGridSource::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
//  if (!this->StructuredOutput) {
    int numInputPorts = this->GetNumberOfInputPorts();
    for (int i=0; i<numInputPorts; i++)
      {
      int numInputConnections = this->GetNumberOfInputConnections(i);
      for (int j=0; j<numInputConnections; j++)
        {
        vtkInformation* inputInfo = inputVector[i]->GetInformationObject(j);
        inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
        }
      }
//  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkRegularGridSource::ComputeInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  
  vtkDataSet *inData = inInfo ? vtkDataSet::SafeDownCast
    (inInfo->Get(vtkDataObject::DATA_OBJECT())) : NULL;

  //
  // Define box...
  //
  double bounds[6], lengths[3], zerovec[3] = {0.0, 0.0, 0.0};
  if (inData && this->UseAutoPlacement) {
    inData->GetBounds(bounds);
    vtkBoundingBox box(bounds);
    box.Inflate(this->Delta);
    box.GetLengths(lengths);
    box.GetBounds(bounds);
    this->origin[0] = bounds[0];
    this->origin[1] = bounds[2];
    this->origin[2] = bounds[4];
    axesvectors[0][0] = axesvectors[1][0] = axesvectors[2][0] = 0.0;
    axesvectors[0][1] = axesvectors[1][1] = axesvectors[2][1] = 0.0;
    axesvectors[0][2] = axesvectors[1][2] = axesvectors[2][2] = 0.0;
    axesvectors[0][0] += lengths[0];
    axesvectors[1][1] += lengths[1];
    axesvectors[2][2] += lengths[2];
  }
  else {
    for (int i=0; i<3; i++) { // @TODO Check this, delta not included
      this->origin[i]  = this->Origin[i];
      axesvectors[0][i] = this->Point1[i] - this->Origin[i];
      axesvectors[1][i] = this->Point2[i] - this->Origin[i];
      axesvectors[2][i] = this->Point3[i] - this->Origin[i];
      lengths[i] = sqrt(vtkMath::Distance2BetweenPoints(zerovec, axesvectors[i]));
    }
  }
  //
  // Define sampling box...
  //
  double o2[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; i++) {
    if (this->Resolution[i]>1) {
      for (int j=0; j<3; j++) o2[j] += axesvectors[i][j];
    }
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
      this->Dimension[i] = vtkMath::Round(1.0/this->scaling[i] + 0.5);
    }
    else {
      this->Dimension[i] = 1;
    }
  }
  for (int i=0; i<3; i++) {
    this->centre[i] = this->origin[i] + o2[i]/2.0;
  }
  //
  this->NumPoints = Dimension[0]*Dimension[1]*Dimension[2];
  //
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
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 
    0, this->Dimension[0]-1, 
    0, this->Dimension[1]-1, 
    0, this->Dimension[2]-1 );
  // Make sure these are correctly set
  outInfo->Set(vtkDataObject::ORIGIN(), this->origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->spacing, 3);
  return 1;
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
  this->ComputeInformation(request, inputVector, outputVector);
  //
  vtkSmartPointer<vtkPoints>     newpoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> dataarray = vtkFloatArray::SafeDownCast(newpoints->GetData());
  float *pointdata = dataarray->WritePointer(0,NumPoints*3);
  //
  vtkIdType totalpoints = 0;
  double    pos1[3], pos2[3];

  //
  // Walk bounds and create points
  //
  for (double k=0; k<Dimension[2]; k++) {
    // increment along Z axis by k
    REGULARGRID_SAXPY(k*this->scaling[2], this->axesvectors[2], this->origin, pos2);
    for (double j=0; j<Dimension[1]; j++) {
      // increment along Y axis by j
      REGULARGRID_SAXPY(j*this->scaling[1], this->axesvectors[1], pos2, pos1);
      for (double i=0; i<Dimension[0]; i++) {
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
    outVerts->Delete();
    outPoly->SetPoints(newpoints);
  }
  else if (outGrid) {
    outGrid->SetDimensions(this->Dimension);
    outGrid->SetWholeExtent(0, Dimension[0]-1, 
                            0, Dimension[1]-1, 
                            0, Dimension[2]-1);
    outGrid->SetPoints(newpoints);
  }
  if (this->Dimension[2]==1) {
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
