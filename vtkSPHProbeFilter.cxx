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
// For PARAVIEW_USE_MPI 
#include "vtkPVConfig.h"     
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkDummyController.h"
//
#include "vtkSPHProbeFilter.h"
#include "vtkSPHManager.h"
//
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
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkDataSet.h"
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
#include "vtkCellArray.h"
#include "vtkBitArray.h"
//
#include "vtkParticleBoxTree.h"
//

#include "KernelGaussian.h"
#include "KernelWendland.h"
#include "KernelQuadratic.h"
#include "KernelSpline3rdOrder.h"
#include "KernelSpline5thOrder.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSPHProbeFilter);
vtkCxxSetObjectMacro(vtkSPHProbeFilter, SPHManager, vtkSPHManager);
vtkCxxSetObjectMacro(vtkSPHProbeFilter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
#if 0

  #define OUTPUTTEXT(a) std::cout <<(a); std::cout.flush();

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    if (this->UpdatePiece>=0) { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg << this->UpdatePiece << " : " a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    } \
  }

  #undef  vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void FloatOrDoubleArrayPointer(vtkDataArray *dataarray, float *&F, double *&D) {
  if (dataarray && vtkFloatArray::SafeDownCast(dataarray)) {
    F = vtkFloatArray::SafeDownCast(dataarray)->GetPointer(0);
    D = NULL;
  }
  if (dataarray && vtkDoubleArray::SafeDownCast(dataarray)) {
    D = vtkDoubleArray::SafeDownCast(dataarray)->GetPointer(0);
    F = NULL;
  }
  //
  if (dataarray && !F && !D) {
    vtkGenericWarningMacro(<< dataarray->GetName() << "must be float or double");
  }
}
//----------------------------------------------------------------------------
#define FloatOrDouble(F, D, index) F ? F[index] : D[index]
#define FloatOrDoubleorDefault(F, D, def, index) F ? F[index] : (D ? D[index] : def)
#define FloatOrDoubleSet(F, D) ((F!=NULL) || (D!=NULL))
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkSPHProbeFilter::vtkSPHProbeFilter()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
  // Point based interpolation
  this->Locator                           = vtkSmartPointer<vtkPointLocator>::New();
  this->Timer                             = vtkSmartPointer<vtkTimerLog>::New();
  //
  this->KernelFunction                    = NULL;
  this->DefaultParticleVolume             = 1.0;
  this->DefaultParticleMass               = 1.0;
  this->DensityScalars                    = NULL;
  this->MassScalars                       = NULL;
  this->VolumeScalars                     = NULL;
  this->HScalars                          = NULL;
  this->ComputeDensityFromNeighbourVolume = 0;
  this->ApplyShepardNormalization         = 1;
  this->PassScalars                       = 0;
  this->MaximumSearchRadius               = 0.0;
  this->ParticleTree                      = NULL;
  this->UseParticleTree                   = 0;
  this->ProgressFrequency                 = 20;
  //
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
  //
  this->SPHManager                  = vtkSPHManager::New();
  this->ModifiedNumber              = 0;
  this->UpdatePiece                 = 0;
  this->UpdateNumPieces             = 1;
  this->TraversalAlgorithm          = vtkSPHProbeFilter::LINEAR_TRAVERSAL;
//  this->TraversalAlgorithm          = vtkSPHProbeFilter::NEIGHBOURHOOD_TILED_TRAVERSAL;
}

//----------------------------------------------------------------------------
vtkSPHProbeFilter::~vtkSPHProbeFilter()
{
  this->SetController(NULL);
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
  this->SetInputData(1, probe);
}
//----------------------------------------------------------------------------
vtkDataObject *vtkSPHProbeFilter::GetProbe()
{
  if (this->GetNumberOfInputConnections(1) < 1) {
    return NULL;
  }
  
  return this->GetExecutive()->GetInputData(1, 0);
}
//----------------------------------------------------------------------------
unsigned long vtkSPHProbeFilter::GetMTime()
{
  unsigned long mtime = this->Superclass::GetMTime();
  unsigned long sm_mtime = this->SPHManager->GetMTime();
  mtime = mtime > sm_mtime? mtime : sm_mtime;
  return mtime;
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
    case vtkSPHManager::SPH_KERNEL_WENDLAND:
      this->KernelFunction = new KernelWendland(this->KernelDimension, H);
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
  vtkSmartPointer<vtkDataSet> newoutput = 
    vtkDataSet::SafeDownCast(vtkDataObjectTypes::NewDataObject(outputType));
  newoutput->FastDelete(); // dec ref count
  return newoutput;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestDataObject(
  vtkInformation *, 
  vtkInformationVector  **inputVector, 
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  vtkInformation     *probePtsInfo = NULL;
  vtkDataSet            *probepts = NULL;
  vtkMultiBlockDataSet *mbprobepts = NULL;
  if (this->GetNumberOfInputPorts()>1 && this->GetNumberOfInputConnections(1)>0) {
    probePtsInfo = inputVector[1]->GetInformationObject(0);
    probepts = vtkDataSet::SafeDownCast(
      probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
    mbprobepts = vtkMultiBlockDataSet::SafeDownCast(
      probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
  }
  //
  int outputType = -1; 
  if (probepts) {
    outputType = this->OutputType(probepts);
  }
  else if (mbprobepts) {
    outputType = VTK_MULTIBLOCK_DATA_SET;
  }
  //
  bool ok = (output!=NULL);
  ok = (ok && output->IsA(vtkDataObjectTypes::GetClassNameFromTypeId(outputType))!=0);
  //
  vtkDataObject *newOutput = NULL;
  if (!ok) {
    newOutput = vtkDataObjectTypes::NewDataObject(outputType);
    outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
    newOutput->FastDelete();
    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), newOutput->GetExtentType());
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
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);
  vtkInformation *probePtsInfo = NULL;
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  if (this->GetNumberOfInputPorts()>1 && this->GetNumberOfInputConnections(1)>0) {
    probePtsInfo = inputVector[1]->GetInformationObject(0);
    //
    // Whole extent is taken from 2nd input (probe points)
    //
    int *wh = probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               probePtsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
               probePtsInfo->Get(
               vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()));
    outInfo->Set(vtkDataObject::ORIGIN(), probePtsInfo->Get(vtkDataObject::ORIGIN()),3);
    outInfo->Set(vtkDataObject::SPACING(), probePtsInfo->Get(vtkDataObject::SPACING()),3);
  }
  else {
    outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestUpdateExtent(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *dataInfo     = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);
  vtkInformation *probePtsInfo = NULL;

  int piece, numPieces, ghostLevels;

  piece       = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces   = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevels = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  
  if (this->GetNumberOfInputPorts()>1 && this->GetNumberOfInputConnections(1)>0) {
    probePtsInfo = inputVector[1]->GetInformationObject(0);
  }
  //
  // Pass the piece request through to both inputs (we operate in parallel like this)
  //
  dataInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  dataInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  //
  if (probePtsInfo) {
    probePtsInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
    probePtsInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkSPHProbeFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *dataInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo  = outputVector->GetInformationObject(0);
  //
  vtkDataSet *particles = vtkDataSet::SafeDownCast(
    dataInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet *mboutput = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  vtkInformation     *probePtsInfo = NULL;
  vtkDataSet            *probepts = NULL;
  vtkMultiBlockDataSet *mbprobepts = NULL;
  if (this->GetNumberOfInputPorts()>1 && this->GetNumberOfInputConnections(1)>0) {
    probePtsInfo = inputVector[1]->GetInformationObject(0);
    probepts = vtkDataSet::SafeDownCast(
      probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
    mbprobepts = vtkMultiBlockDataSet::SafeDownCast(
      probePtsInfo->Get(vtkDataObject::DATA_OBJECT()));
  }

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
    this->InitOutput(particles, probepts, output);
    vtkStructuredGrid *grid = vtkStructuredGrid::SafeDownCast(probepts);
    if (grid) {
      int dims[3];
      grid->GetDimensions(dims);
//      output->SetWholeExtent(0, dims[0]-1, 0,dims[1]-1, 0,dims[2]-1);
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
      vtkDataSet *probepts = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
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
          this->InitOutput(particles, probepts, newOutput);
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
  bool usingH = FloatOrDoubleSet(this->HDataF,this->HDataD);
  bool usingV = FloatOrDoubleSet(this->VolumeDataF,this->VolumeDataD);
  bool usingM = FloatOrDoubleSet(this->MassDataF,this->MassDataD);
  bool usingR = FloatOrDoubleSet(this->DensityDataF,this->DensityDataD);
  if (usingH || usingV || (usingM && usingR)) {
    for (vtkIdType index=0; index<this->NumInputParticles; index++) {
      if (usingH) {
        h = FloatOrDouble(this->HDataF, this->HDataD, index);
      }
      else if (usingV) {
        volume = FloatOrDouble(this->VolumeDataF, this->VolumeDataD, index);
        h      = std::pow(volume, dpower)*this->HCoefficient;
      }
      else if (usingM && usingR) {
        mass   = FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
        rho    = FloatOrDoubleorDefault(this->DensityDataF, this->DensityDataD, this->DefaultDensity, index);
        volume = mass/rho;
        h      = std::pow(volume, dpower)*this->HCoefficient;
      }
      cutoff = std::max(cutoff, h*kappa);
    }
  }
  else {
    cutoff = this->KernelFunction->maxDistance();
  }
  if (this->MaximumSearchRadius>0.0 && this->MaximumSearchRadius<cutoff) {
    vtkWarningMacro(<<"Kernel Cutoff radius of " << cutoff << " is being overridden by user setting of MaximumSearchRadius = " << this->MaximumSearchRadius);
    cutoff = this->MaximumSearchRadius;
  }
  return cutoff;
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::KernelCompute(
  double x[3], vtkDataSet *data, vtkIdList *TestPoints, vtkIdList *NearestPoints, 
  double *gradW, double &totalmass, double &maxDistance)
{
  int N = TestPoints->GetNumberOfIds();
  double *point;
  double  shepard_coeff = 0.0;
  double  weight, volume, mass, rho, d, dr, h=0;
  double  dpower = 1.0/this->KernelDimension;
  Vector _gradW(0.0);
  maxDistance = 0.0;
  totalmass   = 0.0;
  int validpoints = 0;
  //
  bool usingH = FloatOrDoubleSet(this->HDataF,this->HDataD);
  bool usingV = FloatOrDoubleSet(this->VolumeDataF,this->VolumeDataD);
  bool usingM = FloatOrDoubleSet(this->MassDataF,this->MassDataD);
  bool usingR = FloatOrDoubleSet(this->DensityDataF,this->DensityDataD);
  //
  double Dc = this->KernelFunction->maxDistance();
  for (int i=0; i<N; i++) {
    vtkIdType index = TestPoints->GetId(i);
    point = data->GetPoint(index);
    d = sqrt(vtkMath::Distance2BetweenPoints(x, point));
    if (d>Dc) {
      continue;
    }
    if (d>maxDistance) maxDistance=d;
    //
    weight = this->KernelFunction->w(d);
    //
    // each case must set a value for mass, volume, weight
    //
    if (usingH) {
      h      = h = FloatOrDouble(this->HDataF, this->HDataD, index);
      dr     = h/this->HCoefficient;
      volume = dr*dr*dr;
      mass   = FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
    }
    else if (usingM && usingR) {
      mass   = FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
      rho    = FloatOrDoubleorDefault(this->DensityDataF, this->DensityDataD, this->DefaultDensity, index);
      volume = mass/rho;
      h      = std::pow(volume, dpower)*this->HCoefficient;
    }
    else if (usingV) {
      volume = FloatOrDouble(this->VolumeDataF, this->VolumeDataD, index);
      h      = std::pow(volume, dpower)*this->HCoefficient;
      mass   = FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
    }
    else {
      mass   = FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
      rho    = FloatOrDoubleorDefault(this->DensityDataF, this->DensityDataD, this->DefaultDensity, index);
      volume = mass/rho;
    }
    //
    // Mass sum for local density 
    //
    totalmass += mass;
    //
    // Weight and shepard summation
    //
    this->weights[validpoints] = weight*volume;
    shepard_coeff             += this->weights[validpoints];
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
    //
    // Stop if we are using too many points
    // 
    NearestPoints->SetId(validpoints, index); 
    if ((++validpoints)>=KERNEL_MAX_NEIGHBOURS) {
      i=N;
    }
  }
  NearestPoints->SetNumberOfIds(validpoints); 
  gradW[0] = _gradW[0];
  gradW[1] = _gradW[1];
  gradW[2] = _gradW[2];
  if (this->ApplyShepardNormalization) {
    for (int i=0; i<validpoints; i++) {
      this->weights[i] /= shepard_coeff;
    }
  }
  this->ScaleCoefficient  = shepard_coeff;
  this->GradientMagnitude = _gradW.magnitude();
}
//----------------------------------------------------------------------------
void vtkSPHProbeFilter::ShepardCompute(
  double x[3], vtkDataSet *data, vtkIdList *NearestPoints, 
  double &totalmass, double &maxDistance)
{
  int num_neighbors = NearestPoints->GetNumberOfIds();
  double total_w = 0.0;
  maxDistance = 0.0;
  totalmass   = 0.0;
  //
  // compute weights using inverse distance squared
  //
  for (int i=0; i<num_neighbors; i++) {
    vtkIdType index = NearestPoints->GetId(i);
    double *point = data->GetPoint(index);
    double d = vtkMath::Distance2BetweenPoints(point, x);
    if (d!=0.0) {
      this->weights[i] = 1.0 / d;
    }
    else {
      this->weights[i] = 1E6;
    }
    total_w += this->weights[i];
    // mass sum for local density (astro)
    totalmass += FloatOrDoubleorDefault(this->MassDataF, this->MassDataD, this->DefaultParticleMass, index);
    // distance to furthest point for volume calculation
    if (d>maxDistance) maxDistance=d;
  }
  //
  // Normalize
  //
  for (int i=0; i<num_neighbors; i++) {
    if (total_w != 0.0)
      this->weights[i] = this->weights[i] / total_w;
    else
      this->weights[i] = 0.0;
  }
  maxDistance = sqrt(maxDistance);
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
  this->MaximumNeighbours           = this->SPHManager->GetMaximumNeighbours();
  this->MaximumSearchRadius               = this->SPHManager->GetMaximumSearchRadius();
  if (this->MaximumNeighbours>KERNEL_MAX_NEIGHBOURS) {
    this->MaximumNeighbours = KERNEL_MAX_NEIGHBOURS;
  }
  Vector::dim = this->KernelDimension;

  this->NumInputParticles = data->GetNumberOfPoints();
  // 
  // Find the arrays to be used for Mass/Density
  // if not present, we will use default values based on particle size
  //
  this->DensityDataF = NULL;
  this->DensityDataD = NULL;
  this->MassDataF    = NULL;
  this->MassDataD    = NULL;
  this->VolumeDataF  = NULL;
  this->VolumeDataD  = NULL;
  this->HDataF       = NULL;
  this->HDataD       = NULL;
  //
  vtkDataArray *MassArray = this->MassScalars ? 
    data->GetPointData()->GetArray(this->MassScalars) : NULL;
  FloatOrDoubleArrayPointer(MassArray, this->MassDataF, this->MassDataD);
  //
  vtkDataArray *DensityArray = this->DensityScalars ?
    data->GetPointData()->GetArray(this->DensityScalars) : NULL;
  FloatOrDoubleArrayPointer(DensityArray, this->DensityDataF, this->DensityDataD);
  //
  vtkDataArray *VolumeArray  = this->VolumeScalars ?
    data->GetPointData()->GetArray(this->VolumeScalars) : NULL;
  FloatOrDoubleArrayPointer(VolumeArray, this->VolumeDataF, this->VolumeDataD);
  //
  vtkDataArray *HArray = this->HScalars ?
    data->GetPointData()->GetArray(this->HScalars ) : NULL;
  FloatOrDoubleArrayPointer(HArray, this->HDataF, this->HDataD);
  //
  return true;
}
//----------------------------------------------------------------------------
bool vtkSPHProbeFilter::ProbeMeshless(vtkDataSet *data, vtkDataSet *probepts, vtkDataSet *output)
{
  double x[3];
  double bounds[6], cutoff; 
  vtkPointData *pd, *outPD;

  vtkDebugMacro(<<"Probing data");
  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

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
  double total=bins[0]*bins[1]*bins[2];
  for (int i=0; i<3; i++) { 
    if (total>50E6 && bins[i]>100) bins[i]=100;
  }
  vtkDebugMacro(<<"Bins are  " << (int)bins[0] << " " << (int)bins[1] << " " << (int)bins[2]);

#ifdef PARAVIEW_USE_MPI
  //
  // for debug in parallel
  //
  vtkMPICommunicator *com = NULL;
  vtkMPIController   *con = vtkMPIController::SafeDownCast(this->Controller);
  if (con) {
    com = vtkMPICommunicator::SafeDownCast(con->GetCommunicator());
    this->UpdatePiece     = com ? com->GetLocalProcessId() : 0;
    this->UpdateNumPieces = com ? com->GetNumberOfProcesses() : 1;
  }
#endif
  
  //
  // ghost cell stuff
  //
  unsigned char updateLevel = static_cast<unsigned char>
    (probepts->GetInformation()->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  unsigned char *ghostdata = 0;
  vtkDataArray* ghostarray = 0;
  vtkIdType G = 0, nonGhost = 0;
  if (probepts->GetPointData()) {
    ghostarray = probepts->GetPointData()->GetArray("vtkGhostLevels");
  }
  // we don't need this check, we'll take any datatype for ghost flags
  if ( (!ghostarray) || (ghostarray->GetDataType() != VTK_UNSIGNED_CHAR) || (ghostarray->GetNumberOfComponents() != 1)) {
    vtkDebugMacro("No appropriate ghost levels field available.");
  }
  else {
    G = ghostarray->GetNumberOfTuples();
    ghostdata=static_cast<vtkUnsignedCharArray *>(ghostarray)->GetPointer(0);
  }
  for (vtkIdType i=0; ghostdata!=NULL && i<G; i++) {
    if (ghostdata[i]==0) nonGhost++;
  }
#ifdef PARAVIEW_USE_MPI
  vtkDebugMacro(<< "Non ghost N = " << nonGhost << " from total " << G);
#endif

  // 
  // setup output dataset
  // 
  this->NumOutputPoints = (nonGhost==0) ? probepts->GetNumberOfPoints() : nonGhost;
  vtkIdType numInputPoints = probepts->GetNumberOfPoints();
  //
  // Allocate storage for output PointData
  //
  pd = data->GetPointData();
  outPD = output->GetPointData();
  outPD->InterpolateAllocate(pd, this->NumOutputPoints, this->NumOutputPoints);

  bool passdata = false;
  bool computesmootheddensity = false;
  bool computesmoothedradius  = false;

  //
  // We may add some new arrays to the point data as well as the interpolated ones.
  //
  vtkSmartPointer<vtkFloatArray> GradArray, ShepardArray;
  vtkSmartPointer<vtkFloatArray> SmoothedDensity,SmoothedRadius;
  if (this->ComputeDensityFromNeighbourVolume && FloatOrDoubleSet(this->MassDataF,this->MassDataD)) {
    SmoothedDensity = vtkSmartPointer<vtkFloatArray>::New();
    SmoothedDensity->SetName("SmoothedDensity");
    SmoothedDensity->SetNumberOfTuples(this->NumOutputPoints);
    computesmootheddensity = true;
    // if we are resampling data onto itself, we can use the actual mass for a particle
    // and go from smoothed density to a smoothed radius
    if (data==probepts) {
      SmoothedRadius = vtkSmartPointer<vtkFloatArray>::New();
      SmoothedRadius->SetName("SmoothedRadius");
      SmoothedRadius->SetNumberOfTuples(this->NumOutputPoints);
      computesmoothedradius = true;
    }
    passdata = (this->PassScalars!=0);
  }
  if (!passdata) {
    GradArray = vtkSmartPointer<vtkFloatArray>::New();
    GradArray->SetName("GradW");
    GradArray->SetNumberOfTuples(this->NumOutputPoints);
    //
    ShepardArray = vtkSmartPointer<vtkFloatArray>::New();
    ShepardArray->SetName("ShepardCoeff");
    ShepardArray->SetNumberOfTuples(this->NumOutputPoints);
  }
  //

  // Nearest Neighbours list setup
  vtkSmartPointer<vtkIdList> TestPoints = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> NearestPoints = vtkSmartPointer<vtkIdList>::New();

  TestPoints->Allocate(1000);
  NearestPoints->Allocate(1000);

  // Nearest Neighbours list setup
//  vtkSmartPointer<vtkBitArray> Visited = vtkSmartPointer<vtkBitArray>::New();
//  Visited->SetNumberOfTuples(numInputPoints);
//  for (vtkIdType ptId=0; ptId<numInputPoints; ptId++) {
//    Visited->SetValue(ptId, 0);
//  }

  //
  // Loop over all probe points, interpolating particle data
  //
  int abort=0, abortdecision=0;;
  vtkIdType progressInterval=this->NumOutputPoints/this->ProgressFrequency + 1;
  vtkIdType outId = 0;
  vtkIdType NeighbourCount= 0;
  double grad[3] = {0,0,0}, totalmass, maxDistance;

  //
  // initialize locator
  //
  if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
    this->Locator->SetDataSet(data);
    this->Locator->SetDivisions((int)bins[0],(int)bins[1],(int)bins[2]);
    this->Locator->BuildLocator();
  }
  else {
/*
    this->ParticleTree = vtkSmartPointer<vtkParticleBoxTree>::New();
    this->ParticleTree->SetDataSet(data);
    this->ParticleTree->SetCacheCellBounds(1);
    this->ParticleTree->SetNumberOfCellsPerNode(8);
    this->ParticleTree->SetParticleSize(cutoff*2.0);
    this->ParticleTree->BuildLocator();
*/
    this->Locator->SetDataSet(data);
    this->Locator->SetDivisions((int)bins[0],(int)bins[1],(int)bins[2]);
    this->Locator->BuildLocator();
  }
  for (vtkIdType ptId=0; ptId<numInputPoints && !abort; ptId++) {
    if (ghostdata && ghostdata[ptId]>0) {
      continue;
    }
    if ( !(ptId % progressInterval) ) {
      double progress = (double)(ptId)/(double)(numInputPoints);
/*
      std::cout << std::setw(5) << this->UpdatePiece << " Sending Progress of "
        << std::setw(6) << ptId << " " << std::setw(6) << numInputPoints << " " 
        << progress << std::endl;
*/
      this->UpdateProgress(progress);
      abort = GetAbortExecute();
    }

    // Get the xyz coordinate of the point in the probe dataset
    probepts->GetPoint(ptId, x);

    TestPoints->Reset();
    NearestPoints->Reset();
    vtkIdType N;
    // get neighbours of this point
    if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
      this->Locator->FindPointsWithinRadius(cutoff, x, TestPoints);
//      this->ParticleTree->FindCellsFast(x, TestPoints);
      N = TestPoints->GetNumberOfIds();
      if (N>KERNEL_MAX_NEIGHBOURS*2) {
        if (++abortdecision >1000) {
          vtkWarningMacro(<< "Too many large neighbour searches were found and calculation has been aborted to save time");
          abort = 1;
        }
      }
    }
    else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
      this->Locator->FindClosestNPoints(this->MaximumNeighbours, x, NearestPoints);
      N = NearestPoints->GetNumberOfIds();
    }
    
    //
    // compute weights from neighbours
    //
    if (N>0) {
      // For local density calculation we need totalmass, maxDistance
      // For gradients we need grad
      if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
        // KernelCompute will reject any testpoints that were outside the kernel radius
        // the final list is placed in NearestPoints and should be used from here on
        this->KernelCompute(x, data, TestPoints, NearestPoints, grad, totalmass, maxDistance);
        // prevent weights array being mangled
        if (NearestPoints->GetNumberOfIds()>KERNEL_MAX_NEIGHBOURS) {
          NearestPoints->SetNumberOfIds(KERNEL_MAX_NEIGHBOURS);
        }
      }
      else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
        this->ShepardCompute(x, data, NearestPoints, totalmass, maxDistance);
      }
    }
    //
    // if weights were calculated ok, do the interpolation
    //
    N=NearestPoints->GetNumberOfIds();
    if (N>0) {
      NeighbourCount += N;
      // If we are only computing smoothed density and passing input scalars straight to output
      if (passdata) {
        outPD->CopyData(pd, ptId, outId);
      }
      else {
        double gradmag = vtkMath::Norm(grad);
        // Interpolate the point scalar/field data
        outPD->InterpolatePoint(pd, outId, NearestPoints, weights);
        // set our extra computed values
        GradArray->SetValue(outId, gradmag/this->ScaleCoefficient);
//        ShepardArray->SetValue(outId, this->ScaleCoefficient);
        ShepardArray->SetValue(outId, (1.0+rand()%100)/100.0);
      }
      // for astrophysics plots we might be computing a smoothed density
      if (computesmootheddensity) {
        double volume = (4.0/3.0)*M_PI*maxDistance*maxDistance*maxDistance;
        if (volume>0.0) {
          // smoothed density, using total mass in whole neighbourhood
          double smootheddensity = totalmass/volume;
          SmoothedDensity->SetValue(outId, smootheddensity);
          // computed radius based on smoothed density and actual mass
          if (computesmoothedradius) {
            // actual mass of this particle from input data
            double mass = FloatOrDouble(this->MassDataF,this->MassDataD, ptId);
            double smoothedradius = pow(0.75*mass/smootheddensity,0.33333333);
            SmoothedRadius->SetValue(outId, smoothedradius);
          }
        }
        else {
          SmoothedDensity->SetValue(outId, 0.0);
          if (computesmoothedradius) {
            SmoothedRadius->SetValue(outId, 0.0);
          }
        }
      }
    }
    else {
      outPD->NullPoint(outId);
        ShepardArray->SetValue(outId, (1.0+rand()%100)/100.0);
      if (!passdata) {
        ShepardArray->SetValue(outId, 0.0);
        GradArray->SetValue(outId, 0.0);
      }
      if (computesmootheddensity) {
        SmoothedDensity->SetValue(outId, 0.0);
        if (computesmoothedradius) {
          SmoothedRadius->SetValue(outId, 0.0);
        }
      }
    }
    outId++;
  }

/*
  if (this->TraversalAlgorithm==NEIGHBOURHOOD_TILED_TRAVERSAL) {
    this->ParticleTree = vtkSmartPointer<vtkParticleBoxTree>::New();
    this->ParticleTree->SetDataSet(data);
    this->ParticleTree->SetCacheCellBounds(1);
    this->ParticleTree->SetNumberOfCellsPerNode(8);
    this->ParticleTree->SetParticleSize(cutoff);
    this->ParticleTree->BuildLocator();

    //
    // initialize locator
    //
    this->Locator->SetDataSet(data);
    this->Locator->SetDivisions((int)bins[0],(int)bins[1],(int)bins[2]);
    this->Locator->BuildLocator();

    for (vtkIdType ptId=0; ptId<numInputPoints && !abort; ptId++) {
      if (ghostdata && ghostdata[ptId]>0) {
        continue;
      }
      if ( !(ptId % progressInterval) ) {
        this->UpdateProgress((double)ptId/numInputPoints);
        abort = GetAbortExecute();
      }

      // Get the xyz coordinate of the point in the probe dataset
      probepts->GetPoint(ptId, x);
      
      //
      // NB. Inflate extends box on both sides, so we use 1/2
      //
//      double *bounds = this->ParticleTree->GetCellBounds(ptId);
//      vtkBoundingBox box(bounds);
//      box.Inflate(cutoff*0.5); 
      SecondPoints->Reset();
      NearestPoints->Reset();
//      this->ParticleTree->FindCellsWithinBounds(box, NearestPoints); 
      this->ParticleTree->FindCellsFast(x, NearestPoints);
      this->Locator->FindPointsWithinRadius(cutoff, x, SecondPoints);

//      vtkIdType N = 0;
      vtkIdType N = SecondPoints->GetNumberOfIds();
      //
      if (N>KERNEL_MAX_NEIGHBOURS) {
        SecondPoints->SetNumberOfIds(KERNEL_MAX_NEIGHBOURS);
        N = KERNEL_MAX_NEIGHBOURS;
      }

      if (N>0) {
        if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
          // Compute the weights (equivalent) for each nearest point
          // also sum masses and find radius of neighbourhood
          double grad[3], totalmass, maxDistance;
          this->KernelCompute(x, data, SecondPoints, grad, totalmass, maxDistance, N);
        
          // If we are passing input scalars straight to output
          if (passdata) {
            outPD->CopyData(pd, ptId, outId);
            maxDistance = 0.01;
            totalmass = 1E8;
          }
          else {
            double gradmag = vtkMath::Norm(grad);
            // Interpolate the point scalar/field data
            outPD->InterpolatePoint(pd, outId, NearestPoints, weights);
            // set our extra computed values
            GradArray->SetValue(outId, gradmag/this->ScaleCoefficient);
            ShepardArray->SetValue(outId, this->ScaleCoefficient);
          }
          // for astrophysics plots we might be computing a smoothed density
          if (computesmootheddensity) {
            double volume = (4.0/3.0)*M_PI*maxDistance*maxDistance*maxDistance;
            if (volume>0.0) {
              // smoothed density, using total mass in whole neighbourhood
              double smootheddensity = totalmass/volume;
              // actual mass of this particle from input data
              double mass = FloatOrDouble(this->MassDataF,this->MassDataD, ptId);
              // computed radius based on smoothed density and actual mass
              double smoothedradius = pow(0.75*mass/smootheddensity,0.33333333);
              SmoothedDensity->SetValue(outId, smootheddensity);
              SmoothedRadius->SetValue(outId, smoothedradius);
            }
            else {
              SmoothedDensity->SetValue(outId, 0.0);
              SmoothedRadius->SetValue(outId, 0.0);
            }
          }
        }
        else if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_SHEPARD) {
          // if Cal_Interpolation_ShepardMethod_All is called, 
          // then please comment out outPD->InterpolatePoint(pd, ptId, NearestPoints, weights);
          // because Cal_Interpolation_ShepardMethod_All calcualtes interpolations.       
          //int inside = Cal_Interpolation_ShepardMethod_All(data, false, x, N, NearestPoints, ptId, outPD);
          Cal_Weights_ShepardMethod(x, data, NearestPoints, this->weights);
          ShepardArray->SetValue(outId, 0.0);
          GradArray->SetValue(outId, 0.0);
          // need to call this interpolation function 
          // because Cal_Weights_ShepardMethod only calculates weights.
          outPD->InterpolatePoint(pd, outId, NearestPoints, weights);
          if (computesmootheddensity) {
            SmoothedDensity->SetValue(outId, 0.0);
            SmoothedRadius->SetValue(outId, 0.0);
          }
        }
      }
      else {
        outPD->NullPoint(outId);
        ShepardArray->SetValue(outId, 0.0);
        GradArray->SetValue(outId, 0.0);
        if (computesmootheddensity) {
          SmoothedDensity->SetValue(outId, 0.0);
          SmoothedRadius->SetValue(outId, 0.0);
        }
      }
      outId++;
    }
*/

  double averageneighbours = NeighbourCount/outId;
  vtkDebugMacro(<< "Average Neighbour count is " << averageneighbours );

  //
  // Make sure arrays are valid if we aborted mid calculation
  //
  if (abort) {
    for (vtkIdType ptId=outId; ptId<numInputPoints; ptId++) {
      if (ghostdata && ghostdata[ptId]>0) {
        continue;
      }
      outPD->NullPoint(outId);
      if (!passdata) {
        ShepardArray->SetValue(outId, 0.0);
        GradArray->SetValue(outId, 0.0);
      }
      if (computesmootheddensity) {
        SmoothedDensity->SetValue(outId, 0.0);
        if (computesmoothedradius) {
          SmoothedRadius->SetValue(outId, 0.0);
        }
      }
      outId++;
    }
  }

  // Add Grad and Shepard arrays
  if (this->InterpolationMethod==vtkSPHManager::POINT_INTERPOLATION_KERNEL) {
    outPD->AddArray(GradArray);
    outPD->AddArray(ShepardArray);
  }
  if (computesmootheddensity) {
    outPD->AddArray(SmoothedDensity);
    if (computesmoothedradius) {
      outPD->AddArray(SmoothedRadius);
    }
  }

  timer->StopTimer();
  double elapsed = timer->GetElapsedTime();
#ifdef PARAVIEW_USE_MPI
  if (com) com->Barrier();
#endif
  vtkDebugMacro(<< "Probe filter complete : exported N = " << outId << "\t Time : " << elapsed << " seconds" );
  return 1;
}

//----------------------------------------------------------------------------
bool vtkSPHProbeFilter::InitOutput(vtkDataSet *data, vtkDataSet *probepts, vtkDataSet *output)
{
  vtkIdType numInputPoints = probepts->GetNumberOfPoints();

  unsigned char *ghostdata = 0;
  vtkDataArray* ghostarray = 0;
  vtkIdType G = 0, nonGhost = 0;
  if (probepts->GetPointData()) {
    ghostarray = probepts->GetPointData()->GetArray("vtkGhostLevels");
  }
  // we don't need this check, we'll take any datatype for ghost flags
  if ( (!ghostarray) || (ghostarray->GetDataType() != VTK_UNSIGNED_CHAR) || (ghostarray->GetNumberOfComponents() != 1)) {
    vtkDebugMacro("No appropriate ghost levels field available.");
  }
  else {
    G = ghostarray->GetNumberOfTuples();
    ghostdata=static_cast<vtkUnsignedCharArray *>(ghostarray)->GetPointer(0);
  }
  for (vtkIdType i=0; ghostdata!=NULL && i<G; i++) {
    if (ghostdata[i]==0) nonGhost++;
  }

  if (this->NumOutputPoints==numInputPoints) {
    // Copy the probe structure to the output
    output->CopyStructure( probepts );
//    output->CopyInformation( probepts );
  }
  else if (vtkPointSet::SafeDownCast(output)) {
    vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
    newPts->SetNumberOfPoints(this->NumOutputPoints);

    vtkPointSet::SafeDownCast(output)->SetPoints(newPts);
    vtkIdType outId = 0;
    for (vtkIdType i=0; i<numInputPoints; i++) {
      if (ghostdata[i]==0) {
        newPts->SetPoint(outId++, probepts->GetPoint(i));
      }
    }
    //
    // Create Vertices for each point (assumption of point probing)
    //
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
    verts->Allocate(this->NumOutputPoints);
    for (vtkIdType c=0; c<this->NumOutputPoints; c++) {
      verts->InsertNextCell(1, &c);
    }
    vtkPolyData *polydata = vtkPolyData::SafeDownCast(output);
    if (polydata) {
      polydata->SetVerts(verts);
      polydata->SetLines(NULL);
      polydata->SetPolys(NULL);
      polydata->SetStrips(NULL);
    }
    vtkUnstructuredGrid *grid = vtkUnstructuredGrid::SafeDownCast(output);
    if (grid) {
      grid->SetCells(VTK_VERTEX, verts);
    }
  }
  else {
    vtkErrorMacro(<<"Unsupported data type with Ghost Cells found");
    return false;
  }
  return true;
}
