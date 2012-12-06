/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkMomentsOfInertiaFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMomentsOfInertiaFilter 
// vtkMomentsOfInertiaFilter 
// Finds the moment of inertia tensor of a collection of particles, then
// displays graphically the principle moments of inertia. Fully parallel.

#ifndef __vtkMomentsOfInertiaFilter_h
#define __vtkMomentsOfInertiaFilter_h

#include "vtkPointSetAlgorithm.h"
#include <string> // some functions use as argument

enum MomentsMPIData 
{
	INERTIA_TENSOR_COLUMN_ZERO,
	INERTIA_TENSOR_COLUMN_ONE,
	INERTIA_TENSOR_COLUMN_TWO
};

class vtkMultiProcessController;
class VTK_EXPORT vtkMomentsOfInertiaFilter : public vtkPointSetAlgorithm
{
public:
  static vtkMomentsOfInertiaFilter *New();
  vtkTypeMacro(vtkMomentsOfInertiaFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // By defualt this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

//BTX

  // Description:
  // Helper function, for serial applications calls UpdateInertiaTensor
  // then UpdateInertiaTensorFinal
  void ComputeInertiaTensor(vtkPointSet* input, std::string massArrayName,
        double* centerPoint,double inertiaTensor[3][3]);

protected:
  vtkMomentsOfInertiaFilter();
  ~vtkMomentsOfInertiaFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

	// Override to specify different type of output
	virtual int FillOutputPortInformation(int vtkNotUsed(port), 
		vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
   	vtkInformationVector**, vtkInformationVector*);

  vtkMultiProcessController *Controller;
  //
  int           UpdatePiece;
  int           UpdateNumPieces;

private:
  vtkMomentsOfInertiaFilter(const vtkMomentsOfInertiaFilter&);  // Not implemented.
  void operator=(const vtkMomentsOfInertiaFilter&);  // Not implemented.

	// Description:
	// helper function to compute the moment of inertia tensor, returns
	// result in inertiaTensor, incremental update, UpdateInertiaTensorFinal
	// must be called to symmetrize, then update signs
	// I00=sum i=1 to n: m_i(y_i^2+z_i^2)
	// I11=sum i=1 to n: m_i(x_i^2+z_i^2)
	// I22=sum i=1 to n: m_i(x_i^2+y_i^2)
	// I01=sum i=1 to n: m_i*x_i*y_i
	// I02=sum i=1 to n: m_i*x_i*y_i
	// I12=sum i=1 to n: m_i*x_i*y_i
	void UpdateInertiaTensor(vtkPointSet* input,
	 	std::string massArrayName, double* centerPoint,
		double inertiaTensor[3][3]);

	// Description
	// I final update which changes the signs of the components as appropriate
	// and symmetrizes
	// I=[[I00,-I01,-I02],[-I10,I11,-I12],[-I20,-I21,I22]]
	// And I is symmetric so
	// I10=I01
	// I20=I02
	// I21=I12
	void UpdateInertiaTensorFinal(vtkPointSet* input, 
		double* centerPoint, double inertiaTensor[3][3]);

	// Description:
	// Create three lines to display in the output, one for each vector
	// extending from the center point in the direction of the vector
	// until the bounds of the data set
	void DisplayVectorsAsLines(vtkPointSet* input, vtkPolyData* output,
		double vectors[3][3], double* centerPoint);
//ETX
};
#endif
