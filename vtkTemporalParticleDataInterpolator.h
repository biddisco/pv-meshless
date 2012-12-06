/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTemporalParticleDataInterpolator.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTemporalParticleDataInterpolator - interpolate temporal datasets
// .SECTION Description
// vtkTemporalParticleDataInterpolator interpolates between two time steps to
// produce new particle data for an arbitrary T.
// vtkTemporalParticleDataInterpolator 

#ifndef __vtkTemporalParticleDataInterpolator_h
#define __vtkTemporalParticleDataInterpolator_h

#include "vtkTemporalInterpolator.h"

class vtkDataSet;
//BTX
class IndexMap;
//ETX

class VTK_EXPORT vtkTemporalParticleDataInterpolator : public vtkTemporalInterpolator
{
public:
  static vtkTemporalParticleDataInterpolator *New();
  vtkTypeMacro(vtkTemporalParticleDataInterpolator, vtkTemporalInterpolator);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Particles whose position in the point list changes need to be referenced
  // by a unique ID. This array must be specified here and is mandatory.
  // If no ID array exists, then the particles are assumed to be in a fixed order
  // and the class vtkTemporalInterpolator should be used instead.
  vtkSetStringMacro(ParticleIdArray);
  vtkGetStringMacro(ParticleIdArray);

  // Description:
  // By default, all parameters are interpolated using linear interpolation.
  // If linear interpolation does not produce good interpolations of particle
  // positions, then polynomial interpolation should be enabled. 
  // When enabled, the velocity vector must be supplied for each particle
  // and this is used to constrain the derivative of a polynomial interpolation
  // based on a cubic kernel. Particle positions and derivaties at {T0,T1}
  // are used to generate a path at intermediate time steps.
  vtkSetMacro(PolynomialInterpolation, int);
  vtkGetMacro(PolynomialInterpolation, int);
  vtkBooleanMacro(PolynomialInterpolation, int);

  // Description:
  // Specify the name of the Velocity vectors array here. 
  // This array is mandatory when using polynomial interpolation
  vtkSetStringMacro(VelocityArray);
  vtkGetStringMacro(VelocityArray);

protected:
  vtkTemporalParticleDataInterpolator();
  ~vtkTemporalParticleDataInterpolator();
  
  void GetMinMaxIndices(vtkDataArray *indices, 
    vtkIdType &minindex,
    vtkIdType &maxindex);

  // Description:
  // Root level interpolation for a concrete dataset object.
  // Point/Cell data and points are interpolated.
  // Needs improving if connectivity is to be handled
  virtual vtkDataSet *InterpolateDataSet(vtkDataSet *in1, 
                                         vtkDataSet *in2,
                                         double ratio);

  // Description:
  // Interpolate the point using polynomial interpolation
  virtual vtkDataArray *InterpolateDataArrayPolynomial(
      double s, double h,
      vtkDataArray **positions, 
      vtkDataArray **velocities, vtkIdType N);

  // Description:
  // Interpolate a single vtkDataArray. Called from the Interpolation routine
  // on the points and pointdata/celldata
  virtual vtkDataArray *InterpolateDataArray(
      double ratio,
      vtkDataArray **arrays, vtkIdType N);

  // Description:
  // Called just before interpolation to ensure each data arrayhas the same 
  // number of tuples
  bool VerifyArrays(vtkDataArray **arrays, int N);

  // Description:
  // Erase existing maps
  void DeleteParticleIdMaps();

  // Description:
  // Generate a Map of Ids to point indexes 
  void GenerateParticleIdMaps(vtkDataSet **arrays, int N);

  char *ParticleIdArray;
  char *VelocityArray;
  int   PolynomialInterpolation;
  int   NumberOfDataSets;
  vtkTimeStamp MapBuildTimes;
//BTX
  IndexMap **ParticleIdMaps;
//ETX

private:
  vtkTemporalParticleDataInterpolator(const vtkTemporalParticleDataInterpolator&);  // Not implemented.
  void operator=(const vtkTemporalParticleDataInterpolator&);  // Not implemented.
};



#endif



