/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkMomentsOfInertiaFilter.cxx,v $
=========================================================================*/
#include "vtkMomentsOfInertiaFilter.h"
#include "AstroVizHelpers.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkCellArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkMultiProcessController.h"
#include "vtkCenterOfMassFilter.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkUnsignedCharArray.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkMomentsOfInertiaFilter);
vtkCxxSetObjectMacro(vtkMomentsOfInertiaFilter,Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::vtkMomentsOfInertiaFilter()
{
  this->UpdatePiece      = 0;
  this->UpdateNumPieces  = 0;
	this->SetInputArrayToProcess(
    0,
    0,
    0,
    vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
    vtkDataSetAttributes::SCALARS);
	this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::~vtkMomentsOfInertiaFilter()
{
 	this->SetController(0);
}

//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkMomentsOfInertiaFilter::FillInputPortInformation(int, 
	vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkMomentsOfInertiaFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::ComputeInertiaTensor(vtkPointSet* input,
	vtkstd::string massArrayName, double* centerPoint, 
	double inertiaTensor[3][3])
{
	this->UpdateInertiaTensor(input,massArrayName,centerPoint,inertiaTensor);
	this->UpdateInertiaTensorFinal(input,centerPoint,inertiaTensor);
}

//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::UpdateInertiaTensor(vtkPointSet* input, 
	vtkstd::string massArrayName, double* centerPoint, 
	double inertiaTensor[3][3])
{
         for (int i=0; i<3; i++) 
           {
           for (int j=0; j<3; j++) 
             {
             inertiaTensor[i][j] = 0.0;
             }
           }
	for(unsigned long nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		// GetTuple function which operates with float
		double* mass=GetDataValue(input,massArrayName.c_str(),nextPointId);
		// get distance from nextPoint to center point
		double* radius = PointVectorDifference(nextPoint,centerPoint);
		// update the components of the inertia tensor
		inertiaTensor[0][0]+=mass[0]*(pow(nextPoint[1],2)+pow(nextPoint[2],2));
		inertiaTensor[1][1]+=mass[0]*(pow(nextPoint[0],2)+pow(nextPoint[2],2));
		inertiaTensor[2][2]+=mass[0]*(pow(nextPoint[0],2)+pow(nextPoint[1],2));
		inertiaTensor[0][1]+=mass[0]*nextPoint[0]*nextPoint[1];		
		inertiaTensor[0][2]+=mass[0]*nextPoint[0]*nextPoint[2];		
		inertiaTensor[1][2]+=mass[0]*nextPoint[1]*nextPoint[2];		
		// Finally, some memory management
		delete [] radius;
		delete [] mass;
		delete [] nextPoint;
		}
}
//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::UpdateInertiaTensorFinal(vtkPointSet* input, 
	double* centerPoint, double inertiaTensor[3][3])
{
	// Update the signs of off diagonal elements
	inertiaTensor[0][1]*=-1;		
	inertiaTensor[0][2]*=-1;		
	inertiaTensor[1][2]*=-1;
	// We didn't compute these components as we know the tensor is symmetric
	// so symmetrizing based on the components we computed
	inertiaTensor[1][0]=inertiaTensor[0][1];
	inertiaTensor[2][0]=inertiaTensor[0][2];		
	inertiaTensor[2][1]=inertiaTensor[1][2];
}
//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::DisplayVectorsAsLines(vtkPointSet* input,
 	vtkPolyData* output, double vectors[3][3], double* centerPoint)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	//setup the colors array
  vtkSmartPointer<vtkUnsignedCharArray> momentNumber = \
 		vtkSmartPointer<vtkUnsignedCharArray>::New();
		momentNumber->SetNumberOfComponents(1);
		momentNumber->SetNumberOfValues(3);
		momentNumber->SetName("moment number");
	// setting origin
	points->InsertNextPoint(centerPoint);
	double scale=ComputeMaxR(input,centerPoint);
	for(int i = 0; i < 3; ++i)
		{
		VecMultConstant(vectors[i],scale);
		points->InsertNextPoint(vectors[i]);
		// creating the lines
		vtkSmartPointer<vtkLine> nextLine = vtkSmartPointer<vtkLine>::New();
			// setting the first point of the line to be the origin
			nextLine->GetPointIds()->SetId(0,0); 
			// setting the second point of the line to be the scaled vector
			nextLine->GetPointIds()->SetId(1,i+1); // i+1 as origin is 0
		// adding the line to the cell array
		lines->InsertNextCell(nextLine);
		// saving the moment number so the lines can be colored by this
		momentNumber->SetValue(i,i+1); // i+1 as want to index array by 1
		}
	// ready to update the output
	output->SetPoints(points);
	output->SetLines(lines);
	output->GetCellData()->AddArray(momentNumber);
}

//----------------------------------------------------------------------------
int vtkMomentsOfInertiaFilter::RequestData(vtkInformation*,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get input and output data
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// Get name of data array containing mass
	vtkDataArray* massArray = this->GetInputArrayToProcess(0, inputVector);
  if (!massArray)
    {
    vtkErrorMacro("Failed to locate mass array");
    return 0;
    }
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	output->Initialize();

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

	// computing the center of mass, works in parallel if necessary
	vtkSmartPointer<vtkCenterOfMassFilter> centerOfMassFilter = \
		vtkSmartPointer<vtkCenterOfMassFilter>::New();
	centerOfMassFilter->SetController(this->Controller);
	// will be != null only for root process or serial
  double calcCenterOfMass[3];
  centerOfMassFilter->ComputeCenterOfMass(input->GetPoints(), massArray, calcCenterOfMass); 
	// finally calculation
	double inertiaTensor[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];
	if(this->UpdateNumPieces<=1)
		{
		// computing the moment of inertia tensor 3x3 matrix, and its
		// eigenvalues and eigenvectors
		this->ComputeInertiaTensor(input,massArray->GetName(),
			calcCenterOfMass,inertiaTensor);
		// finally perform final computation
		vtkMath::Diagonalize3x3(inertiaTensor,eigenvalues,eigenvectors);
		// displaying eigenvectors
		this->DisplayVectorsAsLines(input,output,eigenvectors,calcCenterOfMass);
		return 1;
		}
	else
		{
		double syncedCenterOfMass[3];
		for(int i = 0; i < 3; ++i)
			{
			syncedCenterOfMass[i]=calcCenterOfMass[i];
			}
		// syncs the value of centerOfMass from root to rest of all processes
		this->Controller->Broadcast(syncedCenterOfMass,3,0);
		// computing the moment of inertia tensor 3x3 matrix
		this->UpdateInertiaTensor(input,massArray->GetName(), 
			syncedCenterOfMass,inertiaTensor);
		// TODO: can send this as one array instead of 3
		if(this->UpdatePiece!=0)
			{
			// send result to root
			this->Controller->Send(inertiaTensor[0],3,0,INERTIA_TENSOR_COLUMN_ZERO);
			this->Controller->Send(inertiaTensor[1],3,0,INERTIA_TENSOR_COLUMN_ONE);
			this->Controller->Send(inertiaTensor[2],3,0,INERTIA_TENSOR_COLUMN_TWO);
			return 1;	
			}
		else
			{
			// we are at proc 0, the last proc, we recieve data from all procs > 0
			for (int proc = 1; proc < this->UpdateNumPieces; ++proc)
				{
				double recInertiaTensorColumnZero[3];
				double recInertiaTensorColumnOne[3];
				double recInertiaTensorColumnTwo[3];
				// Receiving
				this->Controller->Receive(recInertiaTensorColumnZero,
					3,proc,INERTIA_TENSOR_COLUMN_ZERO);
				this->Controller->Receive(recInertiaTensorColumnOne,
					3,proc,INERTIA_TENSOR_COLUMN_ONE);
				this->Controller->Receive(recInertiaTensorColumnTwo,
					3,proc,INERTIA_TENSOR_COLUMN_TWO);
				// Updating inertia tensor
				for(int i = 0; i < 3; ++i)
					{
					inertiaTensor[0][i]+=recInertiaTensorColumnZero[i];
					inertiaTensor[1][i]+=recInertiaTensorColumnOne[i];
					inertiaTensor[2][i]+=recInertiaTensorColumnTwo[i];
					}
				}
			// finally perform final computation
			this->UpdateInertiaTensorFinal(input,syncedCenterOfMass,inertiaTensor);
			vtkMath::Diagonalize3x3(inertiaTensor,eigenvalues,eigenvectors);
			// displaying eigenvectors
			this->DisplayVectorsAsLines(input,output,
				eigenvectors,syncedCenterOfMass);
			return 1;
			}
		}
}
