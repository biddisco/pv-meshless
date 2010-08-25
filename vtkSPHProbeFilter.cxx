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
#include "vtkSPHManager.h"

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

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSPHProbeFilter);
vtkCxxSetObjectMacro(vtkSPHProbeFilter, SPHManager, vtkSPHManager);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
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
void Cal_Weights_ShepardMethod(double x[3], vtkDataSet *data, vtkIdList *NearestPoints, double *weights)
{
	int num_neighbors = NearestPoints->GetNumberOfIds();
	double range = Avg_Min_Distance_Between_Pts(data, x, num_neighbors, NearestPoints);
	double total_w = 0.0;

	if (num_neighbors>0) {
		for (int i = 0; i < num_neighbors; i++) {
			vtkIdType index = NearestPoints->GetId(i);
			double *point = data->GetPoint(index);

			if (sqrt(vtkMath::Distance2BetweenPoints(point, x)) < range) {
				weights[i] = 1.0 / vtkMath::Distance2BetweenPoints(point, x);
				total_w += weights[i];
			}
		}
		for (int i = 0; i < num_neighbors; i++) {
			if(total_w != 0.0)
				weights[i] = weights[i] / total_w;
			else
				weights[i] = 0.0;
		}
	}
	else {
		for (int i = 0; i < num_neighbors; i++)
			weights[i] = 0.0;
	}
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkSPHProbeFilter::vtkSPHProbeFilter()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
  // Point based interpolation
  this->Locator                     = vtkSmartPointer<vtkPointLocator>::New();
  this->Timer                       = vtkSmartPointer<vtkTimerLog>::New();
  //
  this->KernelFunction              = NULL;
  this->DefaultParticleVolume       = 1.0;
  this->DefaultParticleMass         = 1.0;
  this->DensityScalars              = NULL;
  this->MassScalars                 = NULL;
  this->VolumeScalars               = NULL;
  this->HScalars                    = NULL;
  //
  this->SPHManager                  = vtkSPHManager::New();
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
  //
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
    case vtkSPHManager::SPH_KERNEL_GAUSSIAN:
      this->KernelFunction = new KernelGaussian(this->KernelDimension, H);
      break;
    case vtkSPHManager::SPH_KERNEL_QUADRATIC:
      this->KernelFunction = new KernelQuadratic(this->KernelDimension, H);
      break;
    case vtkSPHManager::SPH_KERNEL_SPLINE_3RD:
      this->KernelFunction = new KernelSpline3rdOrder(this->KernelDimension, H);
      break;
    case vtkSPHManager::SPH_KERNEL_SPLINE_5TH:
      this->KernelFunction = new KernelSpline5thOrder(this->KernelDimension, H);
      break;
   default :
      this->KernelFunction = new KernelSpline3rdOrder(this->KernelDimension, H);
      break;
  };
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::OutputType(vtkDataSet *probepts) 
{
  int outputType = probepts->GetDataObjectType();
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
    this->InitializeVariables(particles);
    if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
      this->InitializeKernelCoefficients();
    }
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
        this->InitializeVariables(particles);
        if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
          this->InitializeKernelCoefficients();
        }
        bool nonempty = this->ProbeMeshless(particles, probepts, newOutput);
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
bool vtkSPHProbeFilter::InitializeVariables(vtkDataSet *data)
{
  this->InterpolationMethod         = this->SPHManager->GetInterpolationMethod();
  this->DefaultParticleSideLength   = this->SPHManager->GetDefaultParticleSideLength();
  this->DefaultDensity              = this->SPHManager->GetDefaultDensity();
  this->HCoefficient                = this->SPHManager->GetHCoefficient();
  this->KernelType                  = this->SPHManager->GetKernelType();
  this->KernelDimension             = this->SPHManager->GetKernelDimension();
  this->LimitSearchByNeighbourCount = this->SPHManager->GetLimitSearchByNeighbourCount();
  this->MaximumNeighbours           = this->SPHManager->GetMaximumNeighbours();
  this->MaximumRadius               = this->SPHManager->GetMaximumRadius();
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
  return true;
}
//----------------------------------------------------------------------------
bool vtkSPHProbeFilter::ProbeMeshless(vtkPointSet *data, vtkPointSet *probepts, vtkDataSet *output)
{
  vtkIdType ptId, N;
  double x[3];
  double bounds[6], cutoff; 
  vtkPointData *pd, *outPD;

  vtkDebugMacro(<<"Probing data");

  //
  // Locator optimization
  //
  double bins[3] = {50, 50, 50};
  if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
    cutoff = this->GetMaxKernelCutoffDistance(); 
    data->GetBounds(bounds);
    bins[0] = (bounds[1]-bounds[0])/cutoff;
    bins[1] = (bounds[3]-bounds[2])/cutoff;
    bins[2] = (bounds[5]-bounds[4])/cutoff;
  }
  else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD &&
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
    if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
      this->Locator->FindPointsWithinRadius(
        cutoff, x, NearestPoints);
    }
    else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
      if (this->LimitSearchByNeighbourCount) {
        this->Locator->FindClosestNPoints(this->MaximumNeighbours, x, NearestPoints);
      }
      else {
        this->Locator->FindPointsWithinRadius(this->MaximumRadius, x, NearestPoints);
      }
    }
    N = NearestPoints->GetNumberOfIds();
    if (N>KERNEL_MAX_NEIGHBOURS) {
      NearestPoints->SetNumberOfIds(KERNEL_MAX_NEIGHBOURS);
      N = KERNEL_MAX_NEIGHBOURS;
    }
    
    // compute the interpolated scalar value(s)
    if (N>0) {
      if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
        // Compute the weights (equivalent) for each nearest point
        double grad[3];
        this->KernelCompute(x, data, NearestPoints, grad);
        double gradmag = vtkMath::Norm(grad);
        GradArray->SetValue(ptId, gradmag/this->ScaleCoefficient);
        // Interpolate the point data
        outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
        ShepardArray->SetValue(ptId, this->ScaleCoefficient);
      }
      else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
        // if Cal_Interpolation_ShepardMethod_All is called, 
        // then please comment out outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
        // because Cal_Interpolation_ShepardMethod_All calcualtes interpolations.       
        //int inside = Cal_Interpolation_ShepardMethod_All(data, false, x, N, NearestPoints, ptId, outPD);
        Cal_Weights_ShepardMethod(x, data, NearestPoints, this->weights);
        ShepardArray->SetValue(ptId, 0.0);
        GradArray->SetValue(ptId, 0.0);
        // need to call this interpolation function 
        // becasue Cal_Weights_ShepardMethod only calculates weights.
        outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
      }
    }
    else {
      outPD->NullPoint(ptId);
      ShepardArray->SetValue(ptId, 0.0);
      GradArray->SetValue(ptId, 0.0);
    }
  }

  // Add Grad and Shepard arrays
  if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
    outPD->AddArray(GradArray);
    outPD->AddArray(ShepardArray);
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
  I don't remember why I wanted this in, but it stops the probe re-execting
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
