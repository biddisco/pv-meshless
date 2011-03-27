/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: AstroVizHelpers.cxx,v $
=========================================================================*/
#define _USE_MATH_DEFINES
#include "AstroVizHelpers.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformationVector.h"
#include "vtkDataSetAttributes.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkTable.h"
#include "vtkPointLocator.h"
#include "vtkSphereSource.h"
#include "vtkSmartPointer.h"
#include "vtkMultiProcessController.h"
#include "vtkMath.h"
#include <cmath>
/*----------------------------------------------------------------------------
*
* Work with vtkDataArray
*
*---------------------------------------------------------------------------*/
void InitializeDataArray(vtkDataArray* dataArray, const char* arrayName,
	int numComponents, unsigned long numTuples)
{
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetNumberOfTuples(numTuples);
	dataArray->SetName(arrayName);
	//initializes everything to zero
	for(int i=0; i < numComponents; ++i)
		{
		dataArray->FillComponent(i, 0.0);
		}
}

/*----------------------------------------------------------------------------
*
* Work with vtkTable
*
*---------------------------------------------------------------------------*/
void AllocateDataArray(vtkTable* output, const char* arrayName,
 			int numComponents, unsigned long numTuples)
{
	vtkFloatArray* dataArray=vtkFloatArray::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);		
  output->AddColumn(dataArray);
}
/*----------------------------------------------------------------------------
*
* Work with VtkPolyData (a derived class from vtkPointSet)
*
*---------------------------------------------------------------------------*/

//----------------------------------------------------------------------------
vtkIdType SetPointValue(vtkPolyData* output,float pos[])
{
	vtkIdType id=output->GetPoints()->InsertNextPoint(pos);
	output->GetVerts()->InsertNextCell(1, &id);
	return id;
}

//----------------------------------------------------------------------------
float* DoublePointToFloat(double point[])
{
	float* floatPoint = new float[3];
	for(int i = 0; i < 3; ++i)
	{
		floatPoint[i]=static_cast<float>(point[i]);
	}
	return floatPoint;
}

//----------------------------------------------------------------------------
void CreateSphere(vtkPolyData* output,double radius,double center[])
{
	vtkSmartPointer<vtkSphereSource> sphere = \
	 															vtkSmartPointer<vtkSphereSource>::New();
	sphere->SetRadius(radius);
	sphere->SetCenter(center);
	sphere->Update();
	//Setting the points in the output to be those of the sphere
	output->SetPoints(sphere->GetOutput()->GetPoints());
	output->SetVerts(sphere->GetOutput()->GetVerts());
	output->SetPolys(sphere->GetOutput()->GetPolys());
}

/*----------------------------------------------------------------------------
*
* Work with VtkPointSet
*
*---------------------------------------------------------------------------*/
//----------------------------------------------------------------------------
vtkIdType SetPointValue(vtkPointSet* output,float pos[])
{
	vtkIdType id=output->GetPoints()->InsertNextPoint(pos);
	return id;
}

//----------------------------------------------------------------------------
void AllocateDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, unsigned long numTuples)
{
	vtkSmartPointer<vtkFloatArray> dataArray=\
		vtkSmartPointer<vtkFloatArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,
 	int numComponents, unsigned long numTuples)
{
	vtkSmartPointer<vtkDoubleArray> dataArray=\
		vtkSmartPointer<vtkDoubleArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
void AllocateIntDataArray(vtkPointSet* output, const char* arrayName,
 	int numComponents, unsigned long numTuples)
{
	vtkSmartPointer<vtkIntArray> dataArray = \
		vtkSmartPointer<vtkIntArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
void AllocateIdTypeDataArray(vtkPointSet* output, const char* arrayName,
 	int numComponents, unsigned long numTuples)
{
	vtkSmartPointer<vtkIdTypeArray> idArray = \
		vtkSmartPointer<vtkIdTypeArray>::New();
	InitializeDataArray(idArray,arrayName,1,numTuples);
	output->GetPointData()->AddArray(idArray);
}

//----------------------------------------------------------------------------
double* GetPoint(vtkPointSet* output,vtkIdType id)
{
	double* nextPoint=new double[3]; 
	output->GetPoints()->GetPoint(id,nextPoint);
	return nextPoint;
}

//----------------------------------------------------------------------------
void SetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id,float data[])
{
	output->GetPointData()->GetArray(arrayName)->SetTuple(id,data);
}

//----------------------------------------------------------------------------
void SetIdTypeValue(vtkPointSet* output, const char* arrayName,
	const vtkIdType indexId,const vtkIdType globalId)
{
	vtkIdTypeArray::SafeDownCast(output->GetPointData()->GetArray(
		arrayName))->SetValue(indexId,globalId);
}

//----------------------------------------------------------------------------
void SetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id,double data[])
{
	output->GetPointData()->GetArray(arrayName)->SetTuple(id,data);
}


//----------------------------------------------------------------------------
double* GetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id)
{
	double* data=new double[3];
	output->GetPointData()->GetArray(arrayName)->GetTuple(id,data);
	return data;
}

/*----------------------------------------------------------------------------
*
* Work with vtkInformationVector and vtkInformation objets
*
*---------------------------------------------------------------------------*/
//----------------------------------------------------------------------------
vtkInformationVector** DeepCopyInputVector(vtkInformationVector** inputVector,
	unsigned long inputVectorSize)
{
	vtkInformationVector** newInputVector = \
		new vtkInformationVector *[inputVectorSize];
	for(unsigned long i = 0; i < inputVectorSize; ++i)
		{
		vtkInformationVector* newInput = vtkInformationVector::New();
		// performas a deep copy of inputVector
		newInput->Copy(inputVector[i],1);
		newInputVector[i]=newInput;
		}
	return newInputVector; // note caller must manage this memory
												// see helper function DeleteDeepCopyInput below
}

/*----------------------------------------------------------------------------
*
* Work with vtkInformationVector and vtkInformation objets
*
*---------------------------------------------------------------------------*/
//----------------------------------------------------------------------------
bool RunInParallel(vtkMultiProcessController* controller)
{
	return (controller != NULL && controller->GetNumberOfProcesses() > 1);
}

//----------------------------------------------------------------------------
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,\
											double r,double s,double xacc,double yacc,\
											int *pnIter) 
{
	// This code copied from Doug Potter's and Joachim Stadel's
	// pkdgrav, class master.c, method illinois
	// Changed to return -1 if there are problems finding root
  const int maxIter = 100;
  double t,fr,fs,ft,phis,phir,gamma;
  int i;

  fr = func(r,ctx);
  fs = func(s,ctx);
	if(fr*fs>0)
		{
			return -1;
		}
  t = (s*fr - r*fs)/(fr - fs);

  for(i=0; i<maxIter && fabs(t-s) > xacc; ++i) 
		{
		ft = func(t,ctx);
		if (fabs(ft)<=yacc)
		 {
		 break;
		 }
		if(ft*fs<0) 
			{
	   	/*
	   	** Unmodified step.
	   	*/
	   	r = s;
	   	s = t;
	   	fr = fs;
	   	fs = ft;
			}
		else 
			{
	   	/*
	   	** Modified step to make sure we do not retain the 
	   	** endpoint r indefinitely.
	   	*/
			phis = ft/fs;
	    phir = ft/fr;
	    gamma = 1 - (phis/(1-phir));  /* method 3, for true illinois gamma=0.5*/
	    if(gamma < 0) 
				{
				gamma = 0.5;
	   		}
	   	fr *= gamma;
	   	s = t;
	   	fs = ft;
			}
		t = (s*fr - r*fs)/(fr - fs);
		}	
		
  if(pnIter)
		{
		*pnIter = i;
		}
		
  return(t);
}

//----------------------------------------------------------------------------
double ComputeMaxRadiusInParallel(
	vtkMultiProcessController* controller,vtkPointSet* input, double point[])
{
	double maxR=ComputeMaxR(input,point);
	if(RunInParallel(controller))
		{
		int procId=controller->GetLocalProcessId();
		int numProc=controller->GetNumberOfProcesses();
		if(procId==0)
			{
			// collecting and updating maxR from other processors			
			for(int proc= 1; proc < numProc; ++proc)
				{
				double recMaxR=0;
				controller->Receive(&recMaxR,1,proc,MAX_R);
				maxR=vtkstd::max(maxR,recMaxR);
				}
			// Syncronizing global maxR results
			controller->Broadcast(&maxR,1,0);
			}
		else
			{
			// sending to process 0, which will compare all results and compute
			// global maximum
			controller->Send(&maxR,1,0,MAX_R);
			// syncronizing global maxR results
			controller->Broadcast(&maxR,1,0);
			}
		}
	return maxR;
}
//----------------------------------------------------------------------------

double ComputeMaxR(vtkPointSet* input,double point[])
{
	double bounds[6]; //xmin,xmax,ymin,ymax,zmin,zmax
	input->GetPoints()->ComputeBounds();
	input->GetPoints()->GetBounds(bounds);
	double maxR=0;
	double testR=0;
	// for each of the 8 corners of the bounding box, compute the 
	// distance to the point. maxR is the max distance.
	for(int x = 0; x < 2; ++x)
		{
		for(int y = 2; y <4; ++y)
			{
			for(int z = 4; z < 6; ++z)
				{
				double testCorner[3] = {bounds[x],bounds[y],bounds[z]};
				testR = sqrt(vtkMath::Distance2BetweenPoints(testCorner,point));
				// only if our test R is greater than the current max do we update
				maxR=vtkstd::max(maxR,testR);
				}
			}
		}
	return maxR;
}

//----------------------------------------------------------------------------
double OverDensityInSphere(double r,void* inputVirialRadiusInfo)
{
	VirialRadiusInfo* virialRadiusInfo = \
		static_cast<VirialRadiusInfo*>(inputVirialRadiusInfo);
	vtkIdList* pointsInRadius = \
		FindPointsWithinRadius(r,virialRadiusInfo->center,
		virialRadiusInfo->locator);
	// calculating the average mass, dividing this by the volume of the sphere
	// to get the density
	double totalMass=0;
	vtkPointSet* dataSet=\
		vtkPointSet::SafeDownCast(
		virialRadiusInfo->locator->GetDataSet());
	for(unsigned long pointLocalId = 0; 
			pointLocalId < pointsInRadius->GetNumberOfIds(); 
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		double* nextPoint=GetPoint(dataSet,pointGlobalId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have float method
		double* mass=GetDataValue(dataSet,virialRadiusInfo->massArrayName.c_str(),
			pointGlobalId);
		totalMass+=mass[0];
		// Finally, some memory management
		delete [] mass;
		delete [] nextPoint;
		}
	// If we are running in parallel, update result based on that of other 
	// processors
	if(RunInParallel(virialRadiusInfo->controller))
		{
		int procId=virialRadiusInfo->controller->GetLocalProcessId();
		int numProc=virialRadiusInfo->controller->GetNumberOfProcesses();
		if(procId!=0)
			{
			// Sending to root
			virialRadiusInfo->controller->Send(&totalMass,1,0,TOTAL_MASS_IN_SPHERE);
			// Receiving back the final totalMass from root
			virialRadiusInfo->controller->Broadcast(&totalMass,1,0);
			}
		else
			{
			// Now gather results from each process other than this one
			for(int proc = 1; proc < numProc; ++proc)
				{
				double recTotalMass[1];
				// Receiving
				virialRadiusInfo->controller->Receive(recTotalMass,
					1,proc,TOTAL_MASS_IN_SPHERE);
				// Updating
				totalMass+=recTotalMass[0];
				}
			virialRadiusInfo->controller->Broadcast(&totalMass,1,0);
			}
		}
	// Returning the density minus the critical density. Density is defined
	// as zero if the number points within the radius is zero
	double density = totalMass/(4./3*M_PI*pow(r,3));
	double overdensity = density - 	virialRadiusInfo->criticalValue;
	// managing memory as instructed
	pointsInRadius->Delete();
	return overdensity;
}

//----------------------------------------------------------------------------
double OverNumberInSphere(double r,void* inputVirialRadiusInfo)
{
	VirialRadiusInfo* virialRadiusInfo = \
		static_cast<VirialRadiusInfo*>(inputVirialRadiusInfo);
	vtkIdList* pointsInRadius = \
		FindPointsWithinRadius(r,virialRadiusInfo->center,
		virialRadiusInfo->locator);
	// If we are running in parallel, update result based on that of other 
	// processors
	unsigned long totalNumberInSphere = pointsInRadius->GetNumberOfIds();
	if(RunInParallel(virialRadiusInfo->controller))
		{
		int procId=virialRadiusInfo->controller->GetLocalProcessId();
		int numProc=virialRadiusInfo->controller->GetNumberOfProcesses();
		if(procId!=0)
			{
			// Sending to root
			virialRadiusInfo->controller->Send(&totalNumberInSphere,1,0,
				TOTAL_NUMBER_IN_SPHERE);
			// Making sure we are done, so that the function must finish
			// on all processes before it is attempted to be called again
			virialRadiusInfo->controller->Broadcast(&totalNumberInSphere,1,0);
			}
		else
			{
			// Now gather results from each process other than this one
			for(int proc = 1; proc < numProc; ++proc)
				{
				unsigned long recTotalNumber[1];
				// Receiving
				virialRadiusInfo->controller->Receive(recTotalNumber,
					1,proc,TOTAL_NUMBER_IN_SPHERE);
				// Updating
				totalNumberInSphere+=recTotalNumber[0];
				}
			virialRadiusInfo->controller->Broadcast(&totalNumberInSphere,1,0);
			}
		}
	// Returning the number minus the critical number
	double overNumberInSphere = totalNumberInSphere - \
	 	virialRadiusInfo->criticalValue;
	// Managing memory first, as instructed by FindPoints..
	pointsInRadius->Delete();
	return overNumberInSphere;
}

//----------------------------------------------------------------------------
vtkIdList* FindPointsWithinRadius(double r, 
	double* center, vtkPointLocator* locatorOfThisProcess)
{
	// find points within r, all will need this
	// THIS MEMORY MUST BE MANAGED BY ROOT/CALLER
	vtkIdList* pointsInRadius = vtkIdList::New();
	pointsInRadius->Initialize();
	locatorOfThisProcess->FindPointsWithinRadius(r,
		center,
		pointsInRadius);
	// serial, simply return result
	return pointsInRadius;
}

//----------------------------------------------------------------------------
double* CalculateCenter(vtkDataSet* source)
{
	double* center = new double[3];
	center = source->GetCenter();
	return center;
}

//----------------------------------------------------------------------------
VirialRadiusInfo ComputeVirialRadius(
	vtkMultiProcessController* controller, vtkPointLocator* locator,
	vtkstd::string massArrayName, double softening,double overdensity,
	double maxR,double center[])
{
		// Building the struct to use as argument to root finder and density
		// functions. Contains locator, center, softening info and stores virial
		// radius info for output
		VirialRadiusInfo virialRadiusInfo;
		virialRadiusInfo.locator=locator;
		virialRadiusInfo.controller=controller;
		for(int i = 0; i < 3; ++i)
		{
			virialRadiusInfo.center[i]=center[i];
		}
		virialRadiusInfo.softening=softening;
		virialRadiusInfo.virialRadius = -1; // if stays -1 means not found
		virialRadiusInfo.massArrayName = massArrayName;
		virialRadiusInfo.criticalValue=maxR;
		// but IllinoisRootFinder takes in a void pointer
		void* pntrVirialRadiusInfo = &virialRadiusInfo;
		// 3. Define necessary variables to find virial radius, then search for 
		// it
		int numIter=0; // don't ever use this info, but root finder needs it
	 // keeps track of our guesses and their associated overdensities
	 /// initial guess is the softening
		double guessR[3]={softening,softening,softening};
		double denGuessR[3]={0,0,0}; // keeps track of the density within each R
		int fib[2]={1,1};
		while(guessR[2]<maxR)
			{
			// if our last three guesses have been monotonically decreasing
			// in density, then try to calculate the root
			if(denGuessR[0]>denGuessR[1]>denGuessR[2])
				{
				virialRadiusInfo.criticalValue=overdensity;
				virialRadiusInfo.virialRadius = \
					IllinoisRootFinder(OverDensityInSphere,
					pntrVirialRadiusInfo,
					guessR[0],guessR[2],
					softening,softening,
					&numIter);
				// we are done trying to find the root if the virial radius found is 
				// greater than zero, as rootfinder returns -1 if there were problems 
				if(virialRadiusInfo.virialRadius>0)
					{
					break;
					}
				}
			int nextFib=fib[0]+fib[1];
			// updating the fibonacci sequence
			shiftLeftUpdate(fib,2,nextFib);
			// Updating guessR
			double nextR=nextFib*virialRadiusInfo.softening;
			shiftLeftUpdate(guessR,3,nextR);
			// Updating density estimates
			// Means that OverDensityInSphere will just return DensityInSphere
			virialRadiusInfo.criticalValue=0; 
			shiftLeftUpdate(denGuessR,3,OverDensityInSphere(nextR,
				pntrVirialRadiusInfo));
		}
		return virialRadiusInfo;
}

//----------------------------------------------------------------------------
template <class T> void shiftLeftUpdate(T* array,int size, T updateValue)
{
	// for everything but the last, value is equal to item one to right
	for(int i = 0; i < size-1; ++i)
		{
		array[i]=array[i+1];
		}
	// for last item, value is equal to updateValue
	array[size-1]=updateValue;
}
//----------------------------------------------------------------------------
vtkPointSet* CopyPointsAndData(vtkPointSet* dataSet, vtkIdList*
 	pointsInRadius)
{
	// TODO: I was using CopyCells method of vtkPolyData
	// but this wasn't working so I decided to do manually
	// go back to finding the way using the VTK api to do this
	unsigned long numNewPoints=pointsInRadius->GetNumberOfIds();
	// Initilizing
	vtkPointSet* newDataSet = vtkPolyData::New(); // this memory must be managed
		// Initializing points and verts
	  newDataSet->SetPoints(vtkSmartPointer<vtkPoints>::New());
	  // Initializing data
	vtkSmartPointer<vtkDataArray> nextArray;
	for(int i = 0; i < dataSet->GetPointData()->GetNumberOfArrays(); ++i)
		{
		nextArray = dataSet->GetPointData()->GetArray(i);
		AllocateDataArray(newDataSet,
			nextArray->GetName(),
			nextArray->GetNumberOfComponents(),
			numNewPoints);
		}
	// Copying
	for(unsigned long pointLocalId = 0; 
			pointLocalId < pointsInRadius->GetNumberOfIds(); 
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		double* dbNextPoint=GetPoint(dataSet,pointGlobalId);
		float* nextPoint=DoublePointToFloat(dbNextPoint);
		vtkIdType newId =SetPointValue(newDataSet,nextPoint);
		// adding this to the newDataSet
		// adding this point's data to the newDataSet, for each data array
			for(int i = 0; i < dataSet->GetPointData()->GetNumberOfArrays(); ++i)
				{
				nextArray = dataSet->GetPointData()->GetArray(i);
				double* nextData = GetDataValue(dataSet,
					nextArray->GetName(),pointGlobalId);
				SetDataValue(newDataSet,nextArray->GetName(),newId,nextData);
				delete [] nextData;
				}
			delete [] nextPoint;
			delete [] dbNextPoint;
			}
	return newDataSet;
}



//----------------------------------------------------------------------------
vtkPointSet* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo)
{

	vtkIdList* pointsInRadius = \
		FindPointsWithinRadius(virialRadiusInfo.virialRadius,
		virialRadiusInfo.center, virialRadiusInfo.locator);
  vtkPointSet* dataSet = \
		vtkPointSet::SafeDownCast(virialRadiusInfo.locator->GetDataSet());	
	// Creating a new dataset
	// first allocating
	vtkPointSet* newDataSet = \
		CopyPointsAndData(dataSet,pointsInRadius);
	// Managing memory
	pointsInRadius->Delete();
	return newDataSet;
}


//----------------------------------------------------------------------------
double* ComputeRadialVelocity(double v[],double r[])
{
	return ComputeProjection(v,r);
}

//----------------------------------------------------------------------------
double* ComputeTangentialVelocity(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
  double* vTan=PointVectorDifference(v,vRad);
	delete [] vRad;
	return vTan;	
}
//----------------------------------------------------------------------------
double* ComputeAngularMomentum(double v[], double r[])
{
	double* angularMomentum = new double[3];
	vtkMath::Cross(v,r,angularMomentum);
	return angularMomentum;
}

//----------------------------------------------------------------------------
double* ComputeVelocitySquared(double v[],double r[])
{
	double* velocitySquared = new double[1];
	velocitySquared[0]=vtkMath::Dot(v,v);
	return velocitySquared;
}

//----------------------------------------------------------------------------
double* ComputeRadialVelocitySquared(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
	double* vRadSquared=new double[1];
	vRadSquared[0]=vtkMath::Dot(vRad,vRad);
	delete [] vRad;
	return vRadSquared;
}

//----------------------------------------------------------------------------
double* ComputeTangentialVelocitySquared(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
	double* vTan=ComputeTangentialVelocity(v,r);
	double* vTanSquared=new double[1];
	vTanSquared[0]=vtkMath::Dot(vTan,vTan);
	delete [] vRad;
	delete [] vTan;
	return vTanSquared;	
}

//----------------------------------------------------------------------------
double* ComputeVelocityDispersion(vtkVariant vSquaredAve, vtkVariant vAve)
{
	// vSquared ave required to be a variant which holds a double,
	// vAve required to be a variant which holds a double array with 3 
	// components
	double* velocityDispersion = new double[3];
	for(int comp = 0; comp < 3; ++comp)
		{
		velocityDispersion[comp] = sqrt(fabs(vSquaredAve.ToDouble() -
			pow(vAve.ToArray()->GetVariantValue(comp).ToDouble(),2)));
		}
	return velocityDispersion;
}
//----------------------------------------------------------------------------
double* ComputeCircularVelocity(vtkVariant cumulativeMass, 
	vtkVariant binRadius)
{
	double* circularVelocity = new double[1];
	circularVelocity[0]=cumulativeMass.ToDouble()/binRadius.ToDouble();
	return circularVelocity;
}

//----------------------------------------------------------------------------
double* ComputeDensity(vtkVariant cumulativeMass, 
	vtkVariant binRadius)
{
	double* density = new double[1];
	density[0] = cumulativeMass.ToDouble()/(4./3*vtkMath::Pi()*pow(
		binRadius.ToDouble(),3));
	return density;
}

//----------------------------------------------------------------------------
double* ComputeProjection(double  vectorOne[],double vectorTwo[])
{
	double normVectorTwo = vtkMath::Norm(vectorTwo);
	double projectionMagnitude = \
		vtkMath::Dot(vectorOne,vectorTwo)/normVectorTwo;
	double* projectionVector = new double[3];
	for(int i = 0; i < 3; ++i)
	{
		projectionVector[i] = projectionMagnitude * vectorTwo[i] /normVectorTwo;
	}
	return projectionVector;
}

//----------------------------------------------------------------------------
double* PointVectorDifference(double vectorOne[], double vectorTwo[])
{
	double* pointVectorDifference = new double[3];
	for(int i = 0; i < 3; ++i)
	{
		pointVectorDifference[i] = vectorOne[i] - vectorTwo[i];
	}
	return pointVectorDifference;
}

//----------------------------------------------------------------------------
double* ComputeMidpoint(double pointOne[], double pointTwo[])
{
	double* midpoint = new double[3];
	for(int i = 0; i < 3; ++i)
	{
	midpoint[i] = (pointOne[i] + pointTwo[i])/2;
	}
	return midpoint;
}

//----------------------------------------------------------------------------
void VecMultConstant(double vector[],double constant)
{
	for(int i = 0; i < 3; ++i)
	{
		vector[i] *= constant;
	}
}


//----------------------------------------------------------------------------
	
	


