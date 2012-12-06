/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTemporalParticleDataInterpolator.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*

  The Custom particle interpolator currently only works with the modified
  time pipeline which is a work in progress in the pv-meshless branch of paraview.

  This filter will not produce any output when compiled against the standard paraview

*/

#include "vtkTemporalParticleDataInterpolator.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

#include <algorithm>
#include <vector>

vtkStandardNewMacro(vtkTemporalParticleDataInterpolator);
//----------------------------------------------------------------------------
template <class T>
void vtkTemporalParticleDataInterpolatorFillIndexMap(
  vtkDataArray *indices,
  vtkIdType numTuple, vtkIdType offset,
  std::vector<vtkIdType> &IdToIndexMap, T *)
{
  T *data = static_cast<T*>(indices->GetVoidPointer(0));
  for (vtkIdType i=0; i<numTuple; i++)
  {
	  vtkIdType index = static_cast<vtkIdType>(static_cast<T>(data[i]) - offset);
	  IdToIndexMap[index] = i;
  }
}
//----------------------------------------------------------------------------
// In practice, Id array will be vtkInt or vtkIdType so static_cast will be
// fine as we don't need to round values.
class IndexMap {
	public:
		IndexMap(vtkDataArray *indices, vtkIdType minindex, vtkIdType maxindex)
		{
      // if the first and last particles in the list are not highest/lowest indices, we're stuffed
		  this->MinIndex = minindex;
		  this->MaxIndex = maxindex;
      vtkIdType arraySize = maxindex - minindex + 1;
			this->IdToIndexMap.assign(arraySize, -1);
      switch (indices->GetDataType())
        {
        vtkTemplateMacro(
          vtkTemporalParticleDataInterpolatorFillIndexMap(
          indices, indices->GetNumberOfTuples(), this->MinIndex, 
          this->IdToIndexMap, static_cast<VTK_TT *>(0)));
        }
		}

		// @TODO bounds checking
		vtkIdType GetIndexOfParticle(vtkIdType id) const
		{
			vtkIdType index = id - this->MinIndex;
      if (index >= static_cast<vtkIdType>(this->IdToIndexMap.size()) || index<0) {
				return -1;
      }
      else {
				return this->IdToIndexMap[index];
      }
		}

		vtkIdType operator[](vtkIdType id) const
		{
			return GetIndexOfParticle(id);
		}

    vtkIdType GetMinIndex() const { return this->MinIndex; }
    vtkIdType GetMaxIndex() const { return this->MaxIndex; }

	private:
    std::vector<vtkIdType> IdToIndexMap;
		vtkIdType                 MinIndex;
		vtkIdType                 MaxIndex;
};
//----------------------------------------------------------------------------
vtkTemporalParticleDataInterpolator::vtkTemporalParticleDataInterpolator()
{
  this->NumberOfDataSets        = 0;
  this->ParticleIdArray         = NULL;
  this->ParticleIdMaps          = NULL;
  this->PolynomialInterpolation = 0;
  this->VelocityArray           = NULL;
}
//----------------------------------------------------------------------------
vtkTemporalParticleDataInterpolator::~vtkTemporalParticleDataInterpolator()
{
  this->DeleteParticleIdMaps();
  delete []this->ParticleIdArray;
  delete []this->VelocityArray;
}
//----------------------------------------------------------------------------
void vtkTemporalParticleDataInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "ParticleIdArray: " << this->ParticleIdArray << "\n";
  os << indent << "VelocityArray: " << this->VelocityArray << "\n";
  os << indent << "PolynomialInterpolation: " << this->PolynomialInterpolation << "\n";
}
//----------------------------------------------------------------------------
// @todo : Check ModifiedTime to see if we really need to rebuild maps
void vtkTemporalParticleDataInterpolator::DeleteParticleIdMaps()
{
  for (int i=0; this->ParticleIdMaps!=NULL && i<this->NumberOfDataSets; i++) {
    if (this->ParticleIdMaps[i]) delete this->ParticleIdMaps[i];
  }
  delete []this->ParticleIdMaps;
  this->ParticleIdMaps = NULL;
}
//----------------------------------------------------------------------------
// @todo : Check ModifiedTime to see if we really need to rebuild maps
void vtkTemporalParticleDataInterpolator::GenerateParticleIdMaps(vtkDataSet **datasets, int N)
{
  if (!this->ParticleIdArray) return;
  // 
  this->DeleteParticleIdMaps();
  this->NumberOfDataSets = N;
  //
  this->ParticleIdMaps = new IndexMap*[this->NumberOfDataSets];
	for (vtkIdType i=0; i<this->NumberOfDataSets; i++)
	{
    vtkDataArray  *Ids = datasets[i]->GetPointData()->GetArray(this->ParticleIdArray);
    vtkIdType minindex, maxindex;
    this->GetMinMaxIndices(Ids, minindex, maxindex);
		this->ParticleIdMaps[i] = new IndexMap(Ids, minindex, maxindex);
	}
  this->MapBuildTimes.Modified();
}
//----------------------------------------------------------------------------
template <class T>
void vtkTemporalParticleDataInterpolatorMinMax(vtkTemporalParticleDataInterpolator *,
                                    vtkIdType &minindex,
                                    vtkIdType &maxindex,
                                    vtkDataArray *indices,
                                    vtkIdType numTuple,
                                    T *)
{
  minindex = VTK_INT_MAX;
  maxindex = -1;
  T *data = static_cast<T*>(indices->GetVoidPointer(0));
  for (vtkIdType idx=0; idx<numTuple; ++idx) {
    minindex = std::min(minindex, static_cast<vtkIdType>(data[idx]));
    maxindex = std::max(maxindex, static_cast<vtkIdType>(data[idx]));
  }
}
//----------------------------------------------------------------------------
void vtkTemporalParticleDataInterpolator::GetMinMaxIndices(
  vtkDataArray *indices, 
  vtkIdType &minindex,
  vtkIdType &maxindex)
{
  // now do the interpolation
  switch (indices->GetDataType())
    {
    vtkTemplateMacro(vtkTemporalParticleDataInterpolatorMinMax
      (this, minindex, maxindex, indices, indices->GetNumberOfTuples(),
                      static_cast<VTK_TT *>(0)));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }
}
//----------------------------------------------------------------------------
bool vtkTemporalParticleDataInterpolator
::VerifyArrays(vtkDataArray **vtkNotUsed(arrays), int vtkNotUsed(N))
{
  return true;
}
//----------------------------------------------------------------------------
// templated linear interpolation function
template <class T>
void vtkTemporalParticleDataInterpolatorLinear(vtkTemporalParticleDataInterpolator *,
                                    double ratio,
                                    vtkDataArray *output,
                                    vtkDataArray **arrays,
                                    IndexMap **ParticleIdMaps,
                                    vtkIdType numComp,                                    
                                    T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  T *inData0 = static_cast<T*>(arrays[0]->GetVoidPointer(0));
  T *inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
  //
  vtkIdType MinIndex = std::min(ParticleIdMaps[0]->GetMinIndex(), ParticleIdMaps[1]->GetMinIndex());
	vtkIdType MaxIndex = std::max(ParticleIdMaps[0]->GetMaxIndex(), ParticleIdMaps[1]->GetMaxIndex());
  //
  double oneMinusRatio = 1.0 - ratio;
  vtkIdType missingIds = 0;
  for (vtkIdType t=MinIndex; t<=MaxIndex; ++t)
  {
    vtkIdType idx0 = ParticleIdMaps[0]->GetIndexOfParticle(t);
    vtkIdType idx1 = ParticleIdMaps[1]->GetIndexOfParticle(t);
    if (idx0!=-1 && idx1!=-1) {
      T *x0 = &inData0[idx0*numComp];
      T *x1 = &inData1[idx1*numComp];
      for (int c=0; c<numComp; ++c)
        *outData++ = static_cast<T>(x0[c]*oneMinusRatio + x1[c]*ratio);
    }
    else if (idx0!=-1) {
      T *x0 = &inData0[idx0*numComp];
      for (int c=0; c<numComp; ++c) *outData++ = static_cast<T>(x0[c]);
    }
    else if (idx1!=-1) {
      T *x1 = &inData1[idx1*numComp];
      for (int c=0; c<numComp; ++c) *outData++ = static_cast<T>(x1[c]);
    }
    else {
      missingIds++;
    }
  }
  output->SetNumberOfTuples(MaxIndex - MinIndex + 1 - missingIds);
}
//----------------------------------------------------------------------------
// templated polynomial interpolation function
template <class T>
void vtkTemporalParticleDataInterpolatorPolynomial(vtkTemporalParticleDataInterpolator *,
                                    double s, double h,
                                    vtkDataArray *output,
                                    vtkDataArray **positions,
                                    vtkDataArray **velocities,
                                    IndexMap **ParticleIdMaps,
                                    vtkIdType numComp,                                    
                                    T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  T *inData0 = static_cast<T*>(positions[0]->GetVoidPointer(0));
  T *inData1 = static_cast<T*>(positions[1]->GetVoidPointer(0));
  T *inVel0  = static_cast<T*>(velocities[0]->GetVoidPointer(0));
  T *inVel1  = static_cast<T*>(velocities[1]->GetVoidPointer(0));
  //
  vtkIdType MinIndex = std::min(ParticleIdMaps[0]->GetMinIndex(), ParticleIdMaps[1]->GetMinIndex());
	vtkIdType MaxIndex = std::max(ParticleIdMaps[0]->GetMaxIndex(), ParticleIdMaps[1]->GetMaxIndex());
  //
  double h2 = h*h, h3 = h*h*h;
  double s2 = s*s, s3 = s*s*s;
  double sh3 = s3/h3;
  double sh2 = s2/h2;
  //
  vtkIdType missingIds = 0;
  for (vtkIdType t=MinIndex; t<=MaxIndex; ++t)
  {
    vtkIdType idx0 = ParticleIdMaps[0]->GetIndexOfParticle(t);
    vtkIdType idx1 = ParticleIdMaps[1]->GetIndexOfParticle(t);
    if (idx0!=-1 && idx1!=-1) {
      T *x0 = &inData0[idx0*numComp];
      T *x1 = &inData1[idx1*numComp];
      T *v0 = &inVel0[idx0*numComp];
      T *v1 = &inVel1[idx1*numComp];
      for (int c=0; c<numComp; ++c)
      {
        double dx2  = 2.0*(x1[c]-x0[c]);
        double dx3  = 3.0*(x1[c]-x0[c]);
        *outData++ = static_cast<T>(
         (h*(v0[c] + v1[c]) - dx2)*sh3 - (h*(2.0*v0[c] + v1[c]) - dx3)*sh2 + v0[c]*s + x0[c] 
        );
      }
    }
    else if (idx0!=-1) {
      T *x0 = &inData0[idx0*numComp];
      for (int c=0; c<numComp; ++c) *outData++ = static_cast<T>(x0[c]);
    }
    else if (idx1!=-1) {
      T *x1 = &inData1[idx1*numComp];
      for (int c=0; c<numComp; ++c) *outData++ = static_cast<T>(x1[c]);
    }
    else {
      missingIds++;
    }
  }
  output->SetNumberOfTuples(MaxIndex - MinIndex + 1 - missingIds);
}
//----------------------------------------------------------------------------
vtkDataArray *vtkTemporalParticleDataInterpolator
::InterpolateDataArrayPolynomial(double s, double h, 
    vtkDataArray **positions, vtkDataArray **velocities, vtkIdType N)
{
  //
  // Create the output
  //
  vtkAbstractArray *aa = positions[0]->CreateArray(positions[0]->GetDataType());
  vtkDataArray *output = vtkDataArray::SafeDownCast(aa);
  
  int Nc = positions[0]->GetNumberOfComponents();

  if (positions[0]->GetDataType()!=velocities[0]->GetDataType()) {
    vtkErrorMacro(<< "Currently, the velocity and coordinate arrays must all be float or double, mixed types are not supported");
    return 0;
  }
  //
  // initialize the output
  //
  output->SetNumberOfComponents(Nc);
  output->SetNumberOfTuples(N);
  output->SetName(positions[0]->GetName());

  // now do the interpolation
  switch (positions[0]->GetDataType())
    {
    vtkTemplateMacro(vtkTemporalParticleDataInterpolatorPolynomial
      (this, s, h, output, positions, velocities, this->ParticleIdMaps, Nc,
        static_cast<VTK_TT *>(0)));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return 0;
    }

  return output;
}
//----------------------------------------------------------------------------
vtkDataArray *vtkTemporalParticleDataInterpolator
::InterpolateDataArray(double ratio, vtkDataArray **arrays, vtkIdType N)
{
  //
  // Create the output
  //
  vtkAbstractArray *aa = arrays[0]->CreateArray(arrays[0]->GetDataType());
  vtkDataArray *output = vtkDataArray::SafeDownCast(aa);
  
  int Nc = arrays[0]->GetNumberOfComponents();

  //
  // initialize the output
  //
  output->SetNumberOfComponents(Nc);
  output->SetNumberOfTuples(N);
  output->SetName(arrays[0]->GetName());

  // now do the interpolation
  switch (arrays[0]->GetDataType())
    {
    vtkTemplateMacro(vtkTemporalParticleDataInterpolatorLinear
      (this, ratio, output, arrays, this->ParticleIdMaps, Nc,
        static_cast<VTK_TT *>(0)));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return 0;
    }

  return output;
}
//----------------------------------------------------------------------------
vtkDataSet *vtkTemporalParticleDataInterpolator
::InterpolateDataSet(vtkDataSet *in1, vtkDataSet *in2, double ratio)
{
  vtkDataSet *input[2];
  input[0] = in1;
  input[1] = in2;

  vtkDataArray *velocity[2] = {NULL,NULL};
  if (this->PolynomialInterpolation) {
    if (!this->VelocityArray) {
      vtkErrorMacro(<<"Polynomial Interpolation requires a velocity array");
      return NULL;
    }
    for (int i=0; i<2; i++) {
      velocity[i] = input[i]->GetPointData()->GetArray(this->VelocityArray);
      if (!velocity[i]) {
        vtkErrorMacro(<<"Velocity array " << VelocityArray 
          << " not found for Polynomial Interpolation");
        return NULL;
      }
    }
  }

  //
  if (this->GetMTime()>this->MapBuildTimes ||
      input[0]->GetMTime()>this->MapBuildTimes ||
      input[1]->GetMTime()>this->MapBuildTimes) 
  {
    this->GenerateParticleIdMaps(input, 2);
  }

  //
  vtkDataSet *output = input[0]->NewInstance();
  output->CopyStructure(input[0]);
  //
  // Interpolate points if the dataset is a vtkPointSet
  //
  vtkPointSet *inPointSet1 = vtkPointSet::SafeDownCast(input[0]);
  vtkPointSet *inPointSet2 = vtkPointSet::SafeDownCast(input[1]);
  vtkPointSet *outPointSet = vtkPointSet::SafeDownCast(output);
  if (inPointSet1 && inPointSet2) 
    {
    vtkDataArray *outarray = NULL;
    vtkPoints *outpoints;
    vtkIdType Nt;
    if (inPointSet1->GetNumberOfPoints()>0 && inPointSet2->GetNumberOfPoints()>0)
      {
      vtkDataArray *arrays[2];
      arrays[0] = inPointSet1->GetPoints()->GetData();
      arrays[1] = inPointSet2->GetPoints()->GetData();
      vtkIdType nt0 = arrays[0]->GetNumberOfTuples();
      vtkIdType nt1 = arrays[1]->GetNumberOfTuples();
      Nt = std::max(nt0,nt1);
      if (this->PolynomialInterpolation) {
        outarray = this->InterpolateDataArrayPolynomial(this->Tfrac, this->DeltaT, arrays, velocity, Nt);
      }
      else {
        outarray = this->InterpolateDataArray(this->Ratio, arrays, Nt);
      }
      // Do not shallow copy points from either input, because otherwise when
      // we set the actual point coordinate data we overwrite the original
      // we must instantiate a new points object 
      // (ie we override the copystrucure above)
      vtkPoints *inpoints = inPointSet1->GetPoints();
      outpoints = inpoints->NewInstance();
      outPointSet->SetPoints(outpoints);
      }
    else
      {
      // not much we can do really
      outpoints = vtkPoints::New();
      outPointSet->SetPoints(outpoints);
      Nt = 0;
      }
    outpoints->SetDataType(outarray->GetDataType());
    outpoints->SetNumberOfPoints(outarray->GetNumberOfTuples());
    outpoints->SetData(outarray);
    outpoints->Delete();
    if (outarray)
      {
      outarray->Delete();
      }

    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
    verts->Allocate(outpoints->GetNumberOfPoints());
    for (vtkIdType c=0; c<outpoints->GetNumberOfPoints(); c++) {
      verts->InsertNextCell(1, &c);
    }

    vtkPolyData *polydata = vtkPolyData::SafeDownCast(outPointSet);
    if (polydata) 
      {
      polydata->SetVerts(verts);
      polydata->SetLines(NULL);
      polydata->SetPolys(NULL);
      polydata->SetStrips(NULL);
      }
    vtkUnstructuredGrid *grid = vtkUnstructuredGrid::SafeDownCast(outPointSet);
    if (grid) 
      {
      grid->SetCells(VTK_VERTEX, verts);
      }
    }
  //
  // Interpolate pointdata if present
  //
  output->GetPointData()->ShallowCopy(input[0]->GetPointData());
  for (int s=0; s < input[0]->GetPointData()->GetNumberOfArrays(); ++s) 
    {
    std::vector<vtkDataArray*> arrays;
    char *scalarname = NULL;
    for (int i=0; i<2; ++i) 
      {
      //
      // On some data, the scalar arrays are consistent but ordered
      // differently on each time step, so we will fetch them by name if
      // possible.
      //
      if (i==0 || (scalarname==NULL)) 
        {
        vtkDataArray *dataarray = input[i]->GetPointData()->GetArray(s);
        scalarname = dataarray->GetName();
        arrays.push_back(dataarray);
        }
      else 
        {
        vtkDataArray *dataarray = 
          input[i]->GetPointData()->GetArray(scalarname);
        arrays.push_back(dataarray);
        }
      }
    vtkIdType nt0 = arrays[0]->GetNumberOfTuples();
    vtkIdType nt1 = arrays[1]->GetNumberOfTuples();
    vtkIdType Nt = std::max(nt0,nt1);
    vtkDataArray *outarray = this->InterpolateDataArray(ratio, &arrays[0], Nt);
    output->GetPointData()->AddArray(outarray);
    outarray->Delete();
    }
  //
  // Interpolate celldata if present
  //
  output->GetCellData()->ShallowCopy(input[0]->GetCellData());
  for (int s=0; s<input[0]->GetCellData()->GetNumberOfArrays(); ++s) 
    {
    // copy the structure
    std::vector<vtkDataArray*> arrays;
    char *scalarname = NULL;
    for (int i=0; i<2; ++i) 
      {
      //
      // On some data, the scalar arrays are consistent but ordered
      // differently on each time step, so we will fetch them by name if
      // possible.
      //
      if (i==0 || (scalarname==NULL)) 
        {
        vtkDataArray *dataarray = input[i]->GetCellData()->GetArray(s);
        scalarname = dataarray->GetName();
        arrays.push_back(dataarray);
        }
      else 
        {
        vtkDataArray *dataarray = 
          input[i]->GetCellData()->GetArray(scalarname);
        arrays.push_back(dataarray);
        }
      }
    vtkIdType nt0 = arrays[0]->GetNumberOfTuples();
    vtkIdType nt1 = arrays[1]->GetNumberOfTuples();
    vtkIdType Nt = std::max(nt0,nt1);
    vtkDataArray *outarray = this->InterpolateDataArray(ratio, &arrays[0], Nt);
    output->GetCellData()->AddArray(outarray);
    outarray->Delete();
    }
  if (in1->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()) &&
      in2->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
    {
    output->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
    }
  return output;
}
