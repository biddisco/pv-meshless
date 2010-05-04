/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkCenterOfMassFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCenterOfMassFilter 
// .SECTION Description
// vtkCenterOfMassFilter 
// Finds the center of mass of a collection of particles. Either of all marked
// particles or of all particles. Fully parallel.

#ifndef __vtkCenterOfMassFilter_h
#define __vtkCenterOfMassFilter_h

#include "vtkPointSetAlgorithm.h" // superclass

class vtkMultiProcessController;
class vtkPoints;

enum CenterOfMassMPIData 
{
	TOTAL_MASS,
	TOTAL_WEIGHTED_MASS
};
class VTK_EXPORT vtkCenterOfMassFilter : public vtkPointSetAlgorithm
{
public:
  static vtkCenterOfMassFilter *New();
  vtkTypeMacro(vtkCenterOfMassFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // By defualt this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

	// Description:
	// Computes the center of mass of the vtkPointSet input
	// Functions in parallel if a controller is set. 
  // Returns false if in parallel and process id != 0
	// (result isn't ready until process 0). So check for this.
	//BTX
  bool ComputeCenterOfMass(vtkPoints *points, vtkDataArray *mass, double COM[3]);

protected:
  vtkCenterOfMassFilter();
  ~vtkCenterOfMassFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

	// Override to specify different type of output
	virtual int FillOutputPortInformation(int vtkNotUsed(port), 
		vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  vtkMultiProcessController *Controller;
  //
  int           UpdatePiece;
  int           UpdateNumPieces;
  vtkDataArray *MassArray;
  float        *fPointData;
  double       *dPointData;
  float        *fMassData;
  double       *dMassData;
  vtkIdType     NumberOfPoints;
  //
private:
  vtkCenterOfMassFilter(const vtkCenterOfMassFilter&);  // Not implemented.
  void operator=(const vtkCenterOfMassFilter&);  // Not implemented.
	// Description
	// helper function called by each process to increment totalMass and
	// totalWeighted mass based on the dataset of that process
	// compute COM is called with these variables as input
	// at the very last stage
	void UpdateCenterOfMassVariables(
		double& totalMass, 
		double totalWeightedMass[]);

	// Description:
	// helper function to compute center of mass of point set, must be called
	// at last stage, once update COM vars has been called on each process
  void ComputeCenterOfMassFinal(double &totalMass, 
    double totalWeightedMass[], double *result);

  // Description:
	// ComputeCOM helper function to calculate [m*x,m*y,m*z]
	double* ComputeWeightedMass(double& mass,double* point);
//ETX
};

#endif
