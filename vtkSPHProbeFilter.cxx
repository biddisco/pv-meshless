/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSPHProbeFilter.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSPHProbeFilter.h"

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
//
#include "KernelGaussian.h"
#include "KernelQuadratic.h"
#include "KernelSpline3rdOrder.h"
#include "KernelSpline5thOrder.h"

vtkCxxRevisionMacro(vtkSPHProbeFilter, "$Revision: 1.91 $");
vtkStandardNewMacro(vtkSPHProbeFilter);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int CheckIfPointInsideNeighbors(vtkDataSet *input, int Mesh,
														double *Pt, vtkIdType numNeighbors, vtkIdList *NeighborArray, double range)
{
	double avg[3] = {0.0, 0.0, 0.0};
	double angle;
	double max_angle = 0.0;
	double u[3];
	double v[3];
	double xx[3];
	double lu, lv;
	
	int cnt = 0;

	if(Mesh)
		return 1;
	
#if 1

	for(int j = 0; j < numNeighbors; j++)
	{
		input->GetPoint(NeighborArray->GetId(j), xx);
		if(!Mesh)
		{
			if(sqrt(vtkMath::Distance2BetweenPoints(xx, Pt)) < range)
			{
				avg[0] += xx[0];	
				avg[1] += xx[1];	
				avg[2] += xx[2];
				cnt++;
			}
		}
	}

  if (!Mesh && cnt==0) return 0;

	avg[0] = avg[0] / (double)cnt;
	avg[1] = avg[1] / (double)cnt;
	avg[2] = avg[2] / (double)cnt;

#else

	double den[3] = {0.0, 0.0, 0.0};
	double nom = 0.0;
	double tt = 0.0;
	double h = 0.01;

	for(int j = 0; j < numNeighbors; j++)
	{
		input->GetPoint(NeighborArray->GetId(j), xx);
		tt = exp((-1.0) * vtkMath::Distance2BetweenPoints(xx, Pt) / (h*h));
		den[0] += tt * xx[0];
		den[1] += tt * xx[1];
		den[2] += tt * xx[2];
		nom += tt;
	}


	avg[0] = den[0] / nom;
	avg[1] = den[1] / nom;
	avg[2] = den[2] / nom;
	
	//printf("%f %f %f\n", avg[0], avg[1], avg[2]);

#endif
	
	for(int j = 0; j < numNeighbors; j++)
	{
		input->GetPoint(NeighborArray->GetId(j), xx);
		if(!Mesh)
		{
			if(sqrt(vtkMath::Distance2BetweenPoints(xx, Pt)) < range)
			{
				u[0] = xx[0] - Pt[0];
				u[1] = xx[1] - Pt[1];
				u[2] = xx[2] - Pt[2];
				v[0] = avg[0] - Pt[0];
				v[1] = avg[1] - Pt[1];
				v[2] = avg[2] - Pt[2];
				lu = vtkMath::Normalize(u);
				lv = vtkMath::Normalize(v);
				if(lu != 0.0 && lv != 0.0)
					angle = acos(vtkMath::Dot(u, v));
				else
					//angle = (PI);
					angle = 0.0;

				if(fabs(angle) > max_angle)
					max_angle = fabs(angle);
			}
		}
	}
	//printf("\n");

	//if(Mesh)
	//	return 1;

  if(max_angle > (vtkMath::Pi()) * (1.0 / 2.0))
		return 1;
	else
		return 0;
}
//----------------------------------------------------------------------------
float Avg_Min_Distance_Between_Pts(vtkDataSet *input, double *Pt, 
												int numNeighbors, vtkIdList *NeighborArray)
{
	double xx[3];
	double yy[3];
	double t1 = 1.42383234E20;
	double t2 = 0.0;
	double tt;
	for(int j = 0; j < numNeighbors; j++)
	{
		input->GetPoint(NeighborArray->GetId(j), xx);
		t1 = 1.42383234E20;
		for(int k = 0; k < numNeighbors; k++)
		{	
			if(j != k)
			{
				input->GetPoint(NeighborArray->GetId(k), yy);
				tt = sqrt(vtkMath::Distance2BetweenPoints(xx, yy));
				if(tt != 0.0)
					if(t1 > tt)
						t1 = tt;
			}
		}
		t2 += t1;
	}
	//printf("avg distance = %f\n", t2 / (double) numNeighbors);

	return (3.0 * t2 / (double)numNeighbors);
}
//----------------------------------------------------------------------------
int Cal_Weights_ShepardMethod(double x[3], vtkDataSet *data, vtkIdList *NearestPoints, double *weights)
{
	int i;
	int num_neighbors = NearestPoints->GetNumberOfIds();

	int Inside = 0;
	int Mesh = 0;

	double range = Avg_Min_Distance_Between_Pts(data, x, num_neighbors, NearestPoints);
	Inside = CheckIfPointInsideNeighbors(data, Mesh, x, num_neighbors, NearestPoints, range);

	double total_w = 0.0;

	if(num_neighbors>0)
	{
		for(i = 0; i < num_neighbors; i++)
		{
			vtkIdType index = NearestPoints->GetId(i);
			double *point = data->GetPoint(index);

			if(sqrt(vtkMath::Distance2BetweenPoints(point, x)) < range)
			{
				weights[i] = 1.0 / vtkMath::Distance2BetweenPoints(point, x);
				total_w += weights[i];
			}
		}
		for(i = 0; i < num_neighbors; i++)
		{
			if(total_w != 0.0)
				weights[i] = weights[i] / total_w;
			else
				weights[i] = 0.0;
		}
	}
	else
	{
		for(i = 0; i < num_neighbors; i++)
			weights[i] = 0.0;
	}

	//
	return Inside;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkSPHProbeFilter::vtkSPHProbeFilter()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
  this->SpatialMatch                = 0; // not used yet
  // Point based interpolation
  this->Locator                     = vtkSmartPointer<vtkPointLocator>::New();
  this->Timer                       = vtkSmartPointer<vtkTimerLog>::New();
  this->InterpolationMethod         = POINT_INTERPOLATION_KERNEL;
  // Shepard method
  this->LimitSearchByNeighbourCount = 1;
  this->MaximumNeighbours           = 128;
  this->MaximumRadius               = 0.001;
  //
  this->DefaultParticleSideLength   = 0.0183333333333333333;
  this->HCoefficient                = 1.5;
  this->KernelType                  = 0;
  this->KernelDimension             = 3;
  this->KernelFunction              = NULL;
  this->DefaultParticleVolume       = 1.0;
  this->DensityScalars              = NULL;
  this->MassScalars                 = NULL;
  this->VolumeScalars               = NULL;
  this->HScalars                    = NULL;
}

//----------------------------------------------------------------------------
vtkSPHProbeFilter::~vtkSPHProbeFilter()
{
  if (this->DensityScalars)       delete []this->DensityScalars;
  if (this->MassScalars)          delete []this->MassScalars;
  if (this->VolumeScalars)        delete []this->VolumeScalars;
  if (this->HScalars)             delete []this->HScalars;
  delete KernelFunction;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  }
  else if (port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // we might be multiblock, we might be dataset
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
  return 1;
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::SetProbeConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::SetProbe(vtkDataObject *probe)
{
  this->SetInput(1, probe);
}
//----------------------------------------------------------------------------
vtkDataObject *vtkSPHProbeFilter::GetProbe()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  
  return this->GetExecutive()->GetInputData(1, 0);
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::InitializeKernelCoefficients()
{
  double dpsl = this->DefaultParticleSideLength>0 ? this->DefaultParticleSideLength : 1E-6;
  this->DefaultParticleVolume  = pow(this->DefaultParticleSideLength, this->KernelDimension);
  this->DefaultParticleMass    = this->DefaultDensity*this->DefaultParticleVolume;
  this->ScaleCoefficient       = this->DefaultParticleMass/this->DefaultDensity;
  double H                     = this->HCoefficient*this->DefaultParticleSideLength;
  //
  if (this->KernelFunction) {
    delete this->KernelFunction;
  }
  switch (this->KernelType) {
    case SPH_KERNEL_GAUSSIAN:
      this->KernelFunction = new KernelGaussian(this->KernelDimension, H);
      break;
    case SPH_KERNEL_QUADRATIC:
      this->KernelFunction = new KernelQuadratic(this->KernelDimension, H);
      break;
    case SPH_KERNEL_SPLINE_3RD:
      this->KernelFunction = new KernelSpline3rdOrder(this->KernelDimension, H);
      break;
    case SPH_KERNEL_SPLINE_5TH:
      this->KernelFunction = new KernelSpline5thOrder(this->KernelDimension, H);
      break;
   default :
      this->KernelFunction = new KernelSpline3rdOrder(this->KernelDimension, H);
      break;
  };
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::GetImplicitFunctionStatus(vtkDataSet *probepts)
{
  vtkVariantArray *implicitfn = vtkVariantArray::SafeDownCast(probepts->GetFieldData()->GetAbstractArray("ImplicitPlaneData"));
  return (implicitfn!=NULL);
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::OutputType(vtkDataSet *probepts) 
{
  int outputType = probepts->GetDataObjectType();
  int implicit   = this->GetImplicitFunctionStatus(probepts);
  //
  return probepts->GetDataObjectType();
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkSPHProbeFilter::NewOutput(vtkDataSet *probepts) 
{
  int outputType = this->OutputType(probepts);
  vtkSmartPointer<vtkDataSet> newoutput = vtkDataSet::SafeDownCast(vtkDataObjectTypes::NewDataObject(outputType));
  newoutput->Delete(); // dec ref count
  return newoutput;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestDataObject(
  vtkInformation *, 
  vtkInformationVector  **inputVector, 
  vtkInformationVector *outputVector)
{
  vtkInformation *probePtsInfo = inputVector[1]->GetInformationObject(0);
  // get the input probe geometry
  vtkPointSet *probepts = vtkPointSet::SafeDownCast(
    probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet *mbprobepts = vtkMultiBlockDataSet::SafeDownCast(
    probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));
  //
  // do we have an implicit function desciption?
  // if so the output data type will be UnstructuredGrid, because that's what
  // the cutter exports (if an implicit function is being used for cutting).
  int outputType = -1; 
  if (probepts) {
    outputType = this->OutputType(probepts);
  }
  else if (mbprobepts) outputType = VTK_MULTIBLOCK_DATA_SET;
  //
  bool ok = (output!=NULL);
  ok = (ok && output->IsA(vtkDataObjectTypes::GetClassNameFromTypeId(outputType))!=0);
  //
  vtkDataObject *newOutput = NULL;
  if (!ok) {
    newOutput = vtkDataObjectTypes::NewDataObject(outputType);
    newOutput->SetPipelineInformation(info);
    newOutput->Delete();
    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *probePtsInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);

  // get the input and output, check for multiblock probe/output geometry
  vtkDataSet *data = vtkDataSet::SafeDownCast(
    dataInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet *particles = vtkPointSet::SafeDownCast(
    dataInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet *probepts = vtkPointSet::SafeDownCast(
    probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet *mbprobepts = vtkMultiBlockDataSet::SafeDownCast(
    probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet *mboutput = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!probepts && !mbprobepts)
    {
    return 0;
    }

  //
  this->UpdateProgress(0.0);
  //
  this->Timer->StartTimer();
  if (probepts) {
    this->ProbeMeshless(particles, probepts, output);
    vtkStructuredGrid *grid = vtkStructuredGrid::SafeDownCast(probepts);
    if (grid) {
      int dims[3];
      grid->GetDimensions(dims);
      output->SetWholeExtent(0, dims[0]-1, 0,dims[1]-1, 0,dims[2]-1);
    }
  }
  else if (mbprobepts) {
    int block = 0;
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(mbprobepts->NewIterator());
    // count total number so we can do a progress bar
    int numDataSets = 0, count = 0;
    while (!iter->IsDoneWithTraversal()) {
      numDataSets++;
      iter->GoToNextItem();
    }
    //
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
      vtkPointSet *probepts = vtkPointSet::SafeDownCast(iter->GetCurrentDataObject());
      vtkStdString name = mbprobepts->GetMetaData(iter)->Get(vtkCompositeDataSet::NAME());
      //
      if (probepts) {
        vtkSmartPointer<vtkDataSet> newOutput = this->NewOutput(probepts);
        bool nonempty;
        nonempty = this->ProbeMeshless(particles, probepts, newOutput);
        if (nonempty) {
          mboutput->SetBlock(block, newOutput);
          mboutput->GetMetaData(block++)->Set(vtkCompositeDataSet::NAME(), name.c_str());
        }
      }
      count++;
      this->UpdateProgress((double)count/numDataSets);
    }
  }
  this->Timer->StopTimer();
  double time = this->Timer->GetElapsedTime();
//  vtkErrorMacro(<<"Probe time is " << time);
  this->UpdateProgress(1.0);
  return 1;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
double vtkSPHProbeFilter::GetMaxKernelCutoffDistance()
{
  double volume, mass, rho, h;
  double cutoff = 0.0;
  double kappa  = this->KernelFunction->getDilationFactor();
  double dpower = 1.0/this->KernelDimension;
  //
  if (this->HData || this->MassData || this->VolumeData) {
    for (vtkIdType index=0; index<this->NumInputParticles; index++) {
      //
      mass    = this->MassData    ? this->MassData[index]    : this->DefaultParticleMass;
      rho     = this->DensityData ? this->DensityData[index] : this->DefaultDensity;
      //
      if (this->HData) {
        h      = this->HData[index];
        cutoff = vtkstd::max(cutoff, h*kappa);
      }
      else if (this->MassData) {
        volume = mass/rho;
        h      = vtkstd::pow(volume, dpower)*this->HCoefficient;
        cutoff = vtkstd::max(cutoff, h*kappa);
      }
      else if (this->VolumeData) {
        volume = this->VolumeData[index];
        h      = vtkstd::pow(volume, dpower)*this->HCoefficient;
        cutoff = vtkstd::max(cutoff, h*kappa);
      }
    }
  }
  else {
    cutoff = this->KernelFunction->maxDistance();
  }
  return cutoff;
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::KernelCompute(
  double x[3], vtkPointSet *data, vtkIdList *NearestPoints, double *gradW)
{
  int N = NearestPoints->GetNumberOfIds();
  double *point;
  double  shepard_coeff = 0;
  double  weight, volume, mass, rho, d, dr, h=0;
  double  dpower = 1.0/this->KernelDimension;
  Vector _gradW(0.0);
  //
  for (int i=0; i<N; i++) {
    vtkIdType index = NearestPoints->GetId(i);
    point = data->GetPoint(index);
    d = sqrt(vtkMath::Distance2BetweenPoints(x, point));
    //
    mass    = this->MassData    ? this->MassData[index]    : this->DefaultParticleMass;
    rho     = this->DensityData ? this->DensityData[index] : this->DefaultDensity;
    //
    if (this->HData) {
      h      = this->HData[index];
      dr     = h/this->HCoefficient;
      volume = dr*dr*dr;
      weight = this->KernelFunction->w(h, d);
    }
    else if (this->MassData) {
      volume = mass/rho;
      h      = vtkstd::pow(volume, dpower)*this->HCoefficient;
      weight = this->KernelFunction->w(h, d);
    }
    else if (this->VolumeData) {
      volume = this->VolumeData[index];
      h      = vtkstd::pow(volume, dpower)*this->HCoefficient;
      weight = this->KernelFunction->w(h, d);
    }
    else {
      volume = mass/rho;
      weight = this->KernelFunction->w(d);
    }

    // 
    // Weight and shepard  summation
    //
    this->weights[i] = weight*volume;
    shepard_coeff   += this->weights[i];

    //
    // Gradient vector
    //
    Vector distanceVector(x[0]-point[0], x[1]-point[1], x[2]-point[2]);
    if (d>0) { distanceVector = distanceVector/d; }
    if (h>0) {
      _gradW += this->KernelFunction->gradW(h, d, distanceVector);
    }
    else {
      _gradW += this->KernelFunction->gradW(d, distanceVector);
    }
  }
  gradW[0] = _gradW[0];
  gradW[1] = _gradW[1];
  gradW[2] = _gradW[2];
  this->ScaleCoefficient  = shepard_coeff;
  this->GradientMagnitude = _gradW.magnitude();
}
//----------------------------------------------------------------------------
bool vtkSPHProbeFilter::ProbeMeshless(vtkPointSet *data, vtkPointSet *probepts, vtkDataSet *output)
{
  vtkIdType ptId, N;
  double x[3];
  double bounds[6], cutoff; 
  vtkPointData *pd, *outPD;

  vtkDebugMacro(<<"Probing data");

  Vector::dim = this->KernelDimension;
  this->NumInputParticles = data->GetNumberOfPoints();
  // 
  // Find the arrays to be used for Mass/Density
  // if not present, we will use default values based on particle size
  //
  this->DensityData = NULL;
  this->MassData    = NULL;
  this->VolumeData  = NULL;
  this->HData       = NULL;
  //
  this->MassArray    = this->MassScalars ? 
    data->GetPointData()->GetArray(this->MassScalars) : NULL;
  this->DensityArray = this->DensityScalars ?
    data->GetPointData()->GetArray(this->DensityScalars) : NULL;
  this->VolumeArray  = this->VolumeScalars ?
    data->GetPointData()->GetArray(this->VolumeScalars) : NULL;
  this->HArray       = this->HScalars ?
    data->GetPointData()->GetArray(this->HScalars ) : NULL;
  //
  if (MassArray && vtkFloatArray::SafeDownCast(MassArray)) {
    this->MassData = vtkFloatArray::SafeDownCast(MassArray)->GetPointer(0);
  }
  if (DensityArray && vtkFloatArray::SafeDownCast(DensityArray)) {
    this->DensityData = vtkFloatArray::SafeDownCast(DensityArray)->GetPointer(0);
  }
  if (VolumeArray && vtkFloatArray::SafeDownCast(VolumeArray)) {
    this->VolumeData = vtkFloatArray::SafeDownCast(VolumeArray)->GetPointer(0);
  }
  if (HArray && vtkFloatArray::SafeDownCast(HArray)) {
    this->HData = vtkFloatArray::SafeDownCast(HArray)->GetPointer(0);
  }

  //
  // Precompute any coefficients that we know already
  //
  if (this->InterpolationMethod==POINT_INTERPOLATION_KERNEL) {
    this->InitializeKernelCoefficients();
  }

  //
  // Locator optimization
  //
  double bins[3] = {50, 50, 50};
  if (this->InterpolationMethod==POINT_INTERPOLATION_KERNEL) {
    cutoff = this->GetMaxKernelCutoffDistance(); 
    data->GetBounds(bounds);
    bins[0] = (bounds[1]-bounds[0])/cutoff;
    bins[1] = (bounds[3]-bounds[2])/cutoff;
    bins[2] = (bounds[5]-bounds[4])/cutoff;
  }
  else if (this->InterpolationMethod==POINT_INTERPOLATION_SHEPARD &&
    !this->LimitSearchByNeighbourCount) 
  {
    //data->GetBounds(bounds);
    //bins[3] = { (bounds[1]-bounds[0])/this->MaximumRadius,
    //            (bounds[3]-bounds[2])/this->MaximumRadius,
    //            (bounds[5]-bounds[4])/this->MaximumRadius };
  }
  double total=bins[0]*bins[1]*bins[2];
  for (int i=0; i<3; i++) { 
    if (total>50E6 && bins[i]>100) bins[i]=100;
  }
  vtkDebugMacro(<<"Bins are  " << (int)bins[0] << " " << (int)bins[1] << " " << (int)bins[2]);

  //
  // initialize locator
  //
  this->Locator->SetDataSet(data);
  this->Locator->SetDivisions((int)bins[0],(int)bins[1],(int)bins[2]);
  this->Locator->BuildLocator();

  // 
  // setup output dataset
  // 
  pd = data->GetPointData();
  this->NumOutputPoints = probepts->GetNumberOfPoints();
  // Copy the probe structure to the output
  output->CopyStructure( probepts );
  output->CopyInformation( probepts );
  // Allocate storage for output PointData
  outPD = output->GetPointData();
  outPD->InterpolateAllocate(pd, this->NumOutputPoints, this->NumOutputPoints);

  //
  // we will add some new arrays to the point data as well as the interpolated ones.
  //
  vtkSmartPointer<vtkCharArray> InsideArray = vtkSmartPointer<vtkCharArray>::New();
  InsideArray->SetName("Inside");
  InsideArray->SetNumberOfTuples(this->NumOutputPoints);
  //
  vtkSmartPointer<vtkFloatArray> GradArray = vtkSmartPointer<vtkFloatArray>::New();
  GradArray->SetName("GradW");
  GradArray->SetNumberOfTuples(this->NumOutputPoints);
  //
  vtkSmartPointer<vtkFloatArray> ShepardArray = vtkSmartPointer<vtkFloatArray>::New();
  ShepardArray->SetName("ShepardCoeff");
  ShepardArray->SetNumberOfTuples(this->NumOutputPoints);

  // Nearest Neighbours list setup
  vtkSmartPointer<vtkIdList> NearestPoints = vtkSmartPointer<vtkIdList>::New();

  //
  // Loop over all probe points, interpolating particle data
  //
  int abort=0;
  vtkIdType progressInterval=this->NumOutputPoints/20 + 1;
  for (ptId=0; ptId<this->NumOutputPoints && !abort; ptId++) {
    
    if ( !(ptId % progressInterval) ) {
      this->UpdateProgress((double)ptId/this->NumOutputPoints);
      abort = GetAbortExecute();
    }

    // Get the xyz coordinate of the point in the probe dataset
    probepts->GetPoint(ptId, x);

    // get neighbours of this point
    if (this->InterpolationMethod==POINT_INTERPOLATION_KERNEL) {
      this->Locator->FindPointsWithinRadius(
        cutoff, x, NearestPoints);
    }
    else if (this->InterpolationMethod==POINT_INTERPOLATION_SHEPARD) {
      if (this->LimitSearchByNeighbourCount) {
        this->Locator->FindClosestNPoints(this->MaximumNeighbours, x, NearestPoints);
      }
      else {
        this->Locator->FindPointsWithinRadius(this->MaximumRadius, x, NearestPoints);
      }
    }
    N = NearestPoints->GetNumberOfIds();
    if (N>PROBE_MAX_NEIGHBOURS) {
      NearestPoints->SetNumberOfIds(PROBE_MAX_NEIGHBOURS);
      N = PROBE_MAX_NEIGHBOURS;
    }
    
    // compute the interpolated scalar value(s)
    if (N>0) {
      if (this->InterpolationMethod==POINT_INTERPOLATION_KERNEL) {
        // Compute the weights (equivalent) for each nearest point
        double grad[3];
        this->KernelCompute(x, data, NearestPoints, grad);
        double gradmag = vtkMath::Norm(grad);
        GradArray->SetValue(ptId, gradmag/this->ScaleCoefficient);
        // Interpolate the point data
        outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
        ShepardArray->SetValue(ptId, this->ScaleCoefficient);
      }
      else if (this->InterpolationMethod==POINT_INTERPOLATION_SHEPARD) {
        // if Cal_Interpolation_ShepardMethod_All is called, 
        // then please comment out outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
        // because Cal_Interpolation_ShepardMethod_All calcualtes interpolations.       
        //int inside = Cal_Interpolation_ShepardMethod_All(data, false, x, N, NearestPoints, ptId, outPD);
        int inside = Cal_Weights_ShepardMethod(x, data, NearestPoints, this->weights);
        ShepardArray->SetValue(ptId, 0.0);
        GradArray->SetValue(ptId, 0.0);
        InsideArray->SetValue(ptId, static_cast<char>(inside));
        // need to call this interpolation function 
        // becasue Cal_Weights_ShepardMethod only calculates weights.
        outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
      }
    }
    else {
      outPD->NullPoint(ptId);
      ShepardArray->SetValue(ptId, 0.0);
      GradArray->SetValue(ptId, 0.0);
      InsideArray->SetValue(ptId, static_cast<char>(0));
    }
  }

  // Add Grad and Shepard arrays
  if (this->InterpolationMethod==POINT_INTERPOLATION_KERNEL) {
    outPD->AddArray(GradArray);
    outPD->AddArray(ShepardArray);
  }
  // Add inside/outside array to point data
  if (this->InterpolationMethod==POINT_INTERPOLATION_SHEPARD) {
    outPD->AddArray(InsideArray);
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *probePtsInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);
/*
  I don't remember why I wanted this in, but it stopf the probe re-execting
  when time increments - so leave it out
  outInfo->CopyEntry(dataInfo, 
                     vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->CopyEntry(dataInfo, 
                     vtkStreamingDemandDrivenPipeline::TIME_RANGE());
*/
  // be careful, our output extent is taken from the probe pts which is the 
  // second input
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
               probePtsInfo->Get(
                 vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()));

  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *probePtsInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // we must request the whole extent from probe pts
  probePtsInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
              probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
              6);
    
  return 1;
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkDataObject *probepts = this->GetProbe();

  this->Superclass::PrintSelf(os,indent);
  os << indent << "Probe: " << probepts << "\n";
  os << indent << "SpatialMatch: " << ( this->SpatialMatch ? "On" : "Off" ) << "\n";
}
