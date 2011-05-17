/*=========================================================================

  Project                 : vtkCSCS
  Module                  : vtkParticlePartitionFilter.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkParticlePartitionFilter 
// .SECTION Description
// vtkParticlePartitionFilter 
// Finds the center of mass of a collection of particles. Either of all marked
// particles or of all particles. Fully parallel.

#ifndef __vtkParticlePartitionFilter_h
#define __vtkParticlePartitionFilter_h

#include "vtkPointSetAlgorithm.h" // superclass

class vtkMultiProcessController;
class vtkPoints;

enum CenterOfMassMPIData 
{
	TOTAL_MASS,
	TOTAL_WEIGHTED_MASS
};
class VTK_EXPORT vtkParticlePartitionFilter : public vtkPointSetAlgorithm
{
public:
  static vtkParticlePartitionFilter *New();
  vtkTypeMacro(vtkParticlePartitionFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // By default this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

protected:
  vtkParticlePartitionFilter();
  ~vtkParticlePartitionFilter();

  // Override to specify support for vtkPointSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

	// Override to specify different type of output
	virtual int FillOutputPortInformation(int vtkNotUsed(port), 
		vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  //
  vtkMultiProcessController *Controller;
  //
  int           UpdatePiece;
  int           UpdateNumPieces;
  vtkIdType     NumberOfLocalPoints;
  //
private:
  vtkParticlePartitionFilter(const vtkParticlePartitionFilter&);  // Not implemented.
  void operator=(const vtkParticlePartitionFilter&);  // Not implemented.
};

#endif
