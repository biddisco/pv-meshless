/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkCenterOfMassFilter.cxx,v $
=========================================================================*/
#include "vtkCenterOfMassFilter.h"
#include "vtkPolyData.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
//
#include "vtkDummyController.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCenterOfMassFilter);
vtkCxxSetObjectMacro(vtkCenterOfMassFilter,Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
vtkCenterOfMassFilter::vtkCenterOfMassFilter()
{
  this->UpdatePiece      = 0;
  this->UpdateNumPieces  = 0;
  //this->SetInputDataArrayToProcess(
  //  0,
  //  0,
  //  0,
  //  vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
  //  vtkDataSetAttributes::SCALARS);
  this->MassArray  = NULL;
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::~vtkCenterOfMassFilter()
{
   this->SetController(0);
}

//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::ComputeCenterOfMassFinal(
  double &totalMass, double totalWeightedMass[], double *result)
{
  // calculating the result
  if(totalMass!=0) 
    {
    for(int i = 0; i < 3; ++i)
      {
      result[i] = totalWeightedMass[i]/totalMass;
      }
    }
  else
    {
    vtkErrorMacro("total mass is zero, cannot calculate center of mass, setting center to 0,0,0");
    for(int i = 0; i < 3; ++i)
      {
      result[i] = 0;  
      }
    }
}
//----------------------------------------------------------------------------
template <typename T>
void WeightedSum(T* data, double mass, double *weightedsum) {
  weightedsum[0] += data[0]*mass;
  weightedsum[1] += data[1]*mass;
  weightedsum[2] += data[2]*mass;
}
//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::UpdateCenterOfMassVariables(  
  double &totalMass, double totalWeightedMass[])
{
  vtkIdType index = 0;
  for (vtkIdType id=0; id<this->NumberOfPoints; ++id)
  {
    double mass=this->MassArray->GetTuple1(id);
    if (this->dPointData) {
      WeightedSum<double>(&this->dPointData[index], mass, totalWeightedMass);
    }
    else {
      WeightedSum<float>(&this->fPointData[index], mass, totalWeightedMass);
    }
    totalMass += mass;
    index += 3;
  }
}
//----------------------------------------------------------------------------  
double* vtkCenterOfMassFilter::ComputeWeightedMass(double& mass,double* point)
{
  double* weightedMass = new double[3];
  for(int i = 0; i < 3; ++i)
  {
  weightedMass[i]=mass*point[i];
  }
  return weightedMass;
}

//----------------------------------------------------------------------------
bool vtkCenterOfMassFilter::ComputeCenterOfMass(
  vtkPoints *points, vtkDataArray *mass, double COM[3])
{
  //
  // Setup input data pointers
  //
  this->fPointData = NULL;
  this->dPointData = NULL;
  this->MassArray  = mass;
  //
  if (vtkFloatArray::SafeDownCast(points->GetData())) {
    this->fPointData = vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0);
  }
  else {
    this->dPointData = vtkDoubleArray::SafeDownCast(points->GetData())->GetPointer(0);
  }
  this->NumberOfPoints = points->GetNumberOfPoints();

  //
  // Check parallel operation
  //
  if (this->Controller) {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
  else {
    this->UpdateNumPieces = 1;
    this->UpdatePiece = 0;
  }

  // Allocating data arrays and setting to zero
  double totalMass = 0.0;
  double totalWeightedMass[3] = {0,0,0};
  // testing to make sure we can get to work with D3
  if (this->UpdateNumPieces>1)
    {
    if (this->UpdatePiece!=0)
      {
      // We are at non-root process so simply update and move on
      // Private variables to aid computation of COM
      this->UpdateCenterOfMassVariables(totalMass, totalWeightedMass);
      // Sending to root
      this->Controller->Send(&totalMass,1,0,TOTAL_MASS);
      this->Controller->Send(&totalWeightedMass[0],3,0,TOTAL_WEIGHTED_MASS);
      this->MassArray = NULL;
      return false;
      }
    else
      {
      // We are at root process so update results from root process 
      this->UpdateCenterOfMassVariables(totalMass, totalWeightedMass);

      // Now gather results from each process other than this one
      for (int proc = 1; proc < this->UpdateNumPieces; ++proc)
        {
        double recTotalMass;
        double recTotalWeightedMass[3];

        // Receiving
        this->Controller->Receive(&recTotalMass, 1, proc, TOTAL_MASS);
        this->Controller->Receive(&recTotalWeightedMass[0], 3, proc, TOTAL_WEIGHTED_MASS);

        // Updating
        totalMass += recTotalMass;
        for(int i = 0; i < 3; ++i)
          {
          totalWeightedMass[i] += recTotalWeightedMass[i];         
          }
        }
      this->ComputeCenterOfMassFinal(totalMass, totalWeightedMass, COM);
      }
    }
  else
    {
    // we aren't using MPI or have only one process
    this->UpdateCenterOfMassVariables(totalMass, totalWeightedMass);
    this->ComputeCenterOfMassFinal(totalMass, totalWeightedMass, COM);
    }
  this->MassArray = NULL;
  return true;
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);

  // Get name of data array containing mass
  this->MassArray = this->GetInputArrayToProcess(0, inputVector);
  if (!this->MassArray)
    {
    vtkErrorMacro("Failed to locate mass array");
    return 0;
    }

  // Setup the output
  vtkPolyData* output = vtkPolyData::GetData(outputVector);

  //
  vtkSmartPointer<vtkPoints>   newPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  output->SetPoints(newPoints);
  output->SetVerts(vertices);

  // Compute centre of Mass on local points
  double dbCenterOfMass[3] = {0,0,0};
  bool ok = this->ComputeCenterOfMass(input->GetPoints(), this->MassArray, dbCenterOfMass);
  if (ok)
    {
    // we are in serial or at process 0
    newPoints->SetNumberOfPoints(1);
    newPoints->SetPoint(0, dbCenterOfMass);
    vtkIdType *cells = vertices->WritePointer(1, 2);
    cells[0] = 1;
    cells[1] = 0;
    }

  // Also saving it as a data array for easy csv export
  vtkSmartPointer<vtkFloatArray> c_of_m = vtkSmartPointer<vtkFloatArray>::New();
  c_of_m->SetName("CentreOfMass");
  c_of_m->SetNumberOfComponents(3);
  c_of_m->SetNumberOfTuples(1);
  c_of_m->SetTuple(0, dbCenterOfMass);
  output->GetPointData()->AddArray(c_of_m);
  //
  timer->StopTimer();
  if (this->UpdatePiece==0) { 
//    vtkErrorMacro(<< "Centre Of Mass Calculation : " << timer->GetElapsedTime() << " seconds\n");
  }
  return 1;
}
