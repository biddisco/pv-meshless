/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkExtractValueFilter.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkExtractValueFilter.h"

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
#include "vtkCellArray.h"
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
#include "vtkContourFilter.h"
//
#ifdef VTK_USE_MPI
#include "vtkMPI.h"
#include "vtkMultiProcessController.h"
#include "vtkMPICommunicator.h"
#endif
//
vtkStandardNewMacro(vtkExtractValueFilter);
//
#ifdef VTK_USE_MPI
vtkCxxSetObjectMacro(vtkExtractValueFilter, Controller, vtkMultiProcessController);
#endif
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkExtractValueFilter::vtkExtractValueFilter()
{
  this->ExtractByMaximum = 1;
  this->ExtractionScalars = NULL;
  this->ExtractByCoordinate = 0;
  this->Component = 2;
#ifdef VTK_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkExtractValueFilter::~vtkExtractValueFilter()
{
  delete []this->ExtractionScalars;
}
//----------------------------------------------------------------------------
int vtkExtractValueFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);

  // get the input and output, check for multiblock probe/output geometry
  vtkDataSet *data = vtkDataSet::SafeDownCast(
    dataInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  //
  this->UpdateProgress(0.0);
  //
  vtkDataArray *DensityArray = this->ExtractionScalars ?
    data->GetPointData()->GetArray(this->ExtractionScalars) : NULL;
  //
  this->FindMaximum(data, output, DensityArray);
  //
  this->UpdateProgress(1.0);
  return 1;
}
//----------------------------------------------------------------------------
void vtkExtractValueFilter::FindMaximum(vtkDataSet *input, vtkPolyData *output, vtkDataArray *scalars)
{
  double vmax = VTK_DOUBLE_MIN;
  double vmin = VTK_DOUBLE_MAX;
  vtkIdType i, imax = -1, imin = VTK_INT_MAX;
  
  double val;
  double *p = NULL;
  for (vtkIdType i=0; i<input->GetNumberOfPoints(); i++) {
    if (this->ExtractByCoordinate) {
      p = input->GetPoint(i);
      if (p[this->Component]>vmax) {
        vmax = p[this->Component];
        imax = i;
      }
      if (p[this->Component]<vmin) {
        vmin = p[this->Component];
        imin = i;
      }
    }
    else {
      val = scalars->GetTuple1(i);
      if (val>vmax) {
        vmax = val;
        imax = i;
      }
      if (val<vmin) {
        vmin = val;
        imin = i;
      }
    }

  }
  //
  // Value has been found, do parallel reduction
  //
  bool validresult = true;
#ifdef VTK_USE_MPI
  if (this->Controller) {
    double result;
    vtkIdType rank, lowRank;
    rank = this->Controller->GetLocalProcessId();
    if (this->ExtractByMaximum) {
      this->Controller->AllReduce(&vmax, &result/*(double*)MPI_IN_PLACE*/, 1, vtkCommunicator::MAX_OP);
      if (result > vmax) {
        validresult = false;
      }
    }
    else {
      this->Controller->AllReduce(&vmin, &result/*(double*)MPI_IN_PLACE*/, 1, vtkCommunicator::MIN_OP);
      if (result < vmin) {
        validresult = false;
      }
    }
    // If many processes has the same peak value (like zero!), we choose by rank
    if (!validresult) {
      rank = VTK_INT_MAX;
    }
    this->Controller->AllReduce(&rank, &lowRank, 1, vtkCommunicator::MIN_OP);
    if (validresult) {
      if (rank!=lowRank) {
        validresult = false;
      }
    }
  }
#endif

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
  if (validresult && (imax != -1) && (imin != VTK_INT_MAX)) {
    vtkIdType pt = 0;
    if (this->ExtractByMaximum) {
      i = imax;
    }
    else {
      i = imin;
    }
    points->InsertNextPoint(input->GetPoint(i));
    verts->InsertNextCell(1, &pt);
    output->GetPointData()->CopyAllocate(input->GetPointData());
    output->GetPointData()->CopyData(input->GetPointData(), i, 0);
  }
  output->SetPoints(points);
  output->SetVerts(verts);
}
//----------------------------------------------------------------------------
void vtkExtractValueFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
//  os << indent << "Probe: " << probepts << "\n";
//  os << indent << "SpatialMatch: " << ( this->SpatialMatch ? "On" : "Off" ) << "\n";
}
