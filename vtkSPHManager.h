/*=========================================================================

  Project                 : vtkCSCSMeshless
  Module                  : vtkSPHManager.h

  Authors:
     John Biddiscombe     Jerome Soumagne
     biddisco@cscs.ch     soumagne@cscs.ch

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
// .NAME vtkSPHProbeFilter2 - SPH Convenience Manager Class
// .SECTION Description

#ifndef __vtkSPHManager_h
#define __vtkSPHManager_h

#include "vtkObject.h"
#include "vtkSmartPointer.h" // for smartpointers

//BTX
class Kernel;
//ETX

#define KERNEL_MAX_NEIGHBOURS 256

class VTK_EXPORT vtkSPHManager : public vtkObject
{
public:
  static vtkSPHManager *New();
  vtkTypeMacro(vtkSPHManager,vtkObject);

  //BTX
  enum {
    POINT_INTERPOLATION_KERNEL=0,
    POINT_INTERPOLATION_SHEPARD,
  };
  //ETX

  // Description:
  // Use either SPH Kernel interpolation or Linear interpolation
  vtkSetMacro(InterpolationMethod, int);
  vtkGetMacro(InterpolationMethod, int);
  void SetInterpolationMethodToKernel(int k) { 
    this->SetKernelType(POINT_INTERPOLATION_KERNEL); }
  void SetInterpolationMethodToLinear(int k) { 
    this->SetKernelType(POINT_INTERPOLATION_SHEPARD); }

  // Description:
  // if VolumeScalars are provided (per particle), then the kernel
  // summation term m/rho uses the array provided and the parameters
  // (DensityScalars, DefaultDensity, MassScalars) are ignored.
  vtkSetStringMacro(VolumeScalarsRegex);
  vtkGetStringMacro(VolumeScalarsRegex);

  // Description:
  // For a variable-h simulation, an h-array must be supplied
  // this determines the kernel cutoff on a per-particle basis.
  vtkSetStringMacro(HScalarsRegex);
  vtkGetStringMacro(HScalarsRegex);

  // Description:
  // The Kernel H Factor with a Default value of 1.5.
  // It is also referred to as the LengthRatio and converts the 
  // side length of a cubic particle of volume v into the h value
  // <verbatim>
  // h = HCoefficient * DefaultParticleSideLength
  // </verbatim>
  // the Kernel cutoff is usually 2h (3rd order spline) 
  // or 3h (5th order spline)
  vtkSetMacro(HCoefficient, double);
  vtkGetMacro(HCoefficient, double);
  
  // Description:
  // The DefaultParticleSideLength is the length (in metres) of 1 side of a cube 
  // holding the particle at the start of the simulation. 
  // When using a 3D kernel the volume of a particle will
  // be the cube of DefaultParticleSideLength, and this is used to normalize
  // the summation of kernel contributions given the density.
  // For EDF Test Case 2, DefaultParticleSideLength = 0.01*(55.0/30.0)m
  // Where the 0.01 factor converts cm to m. 
  // The units used must be consistent with the position coordinates
  // of the particles (ie metres, cm, etc)
  // Default value = 0.0183333333333333333
  vtkSetMacro(DefaultParticleSideLength, double);
  vtkGetMacro(DefaultParticleSideLength, double);

  // Description:
  // if each particle has a density value, specify the array name here.
  // if not specified then the value of DefaultDensity will be used
  vtkSetStringMacro(DensityScalarsRegex);
  vtkGetStringMacro(DensityScalarsRegex);

  // Description:
  // Material density is 1000 (Kg/m^3) by default, but can be changed if
  // the simulation requires it, overridden if DensityScalars are supplied.
  vtkSetMacro(DefaultDensity, double);
  vtkGetMacro(DefaultDensity, double);

  // Description:
  // if each particle has a mass value, specify the array name here.
  // if not specified then the value of mass will be computed from particle size
  // and density etc.
  vtkSetStringMacro(MassScalarsRegex);
  vtkGetStringMacro(MassScalarsRegex);

  //BTX
  enum {
    SPH_KERNEL_GAUSSIAN=0,
    SPH_KERNEL_QUADRATIC,
    SPH_KERNEL_SPLINE_3RD,
    SPH_KERNEL_SPLINE_5TH,
    SPH_KERNEL_CUSP,
    SPH_KERNEL_BOX,
    SPH_KERNEL_WENDLAND,
  };
  //ETX

  // Description:
  // Set which kernel is being used, currently the following are supported
  // CubicSPline3D, CubicSPline2D, Cusp3D
  vtkSetMacro(KernelType, int);
  vtkGetMacro(KernelType, int);
  void SetKernelTypeToGaussian(int k) { 
    this->SetKernelType(SPH_KERNEL_GAUSSIAN); }
  void SetKernelTypeToWendland(int k) { 
    this->SetKernelType(SPH_KERNEL_WENDLAND); }
  void SetKernelTypeToQuadratic3D(int k) { 
    this->SetKernelType(SPH_KERNEL_QUADRATIC); }
  void SetKernelTypeToCubicSpline(int k) { 
    this->SetKernelType(SPH_KERNEL_SPLINE_3RD); }
  void SetKernelTypeToQuinticSpline(int k) { 
    this->SetKernelType(SPH_KERNEL_SPLINE_5TH); }
  void SetKernelTypeToCusp(int k) { 
    this->SetKernelType(SPH_KERNEL_CUSP); }

  // Description:
  // Set the kernel dimension, 2D or 3D are permitted
  vtkSetMacro(KernelDimension, int);
  vtkGetMacro(KernelDimension, int);
  
  // Description:
  // Set to true to limit by MaximumNeighbours count,
  // false to limit by radius
  vtkSetMacro(LimitSearchByNeighbourCount, int);
  vtkGetMacro(LimitSearchByNeighbourCount, int);
  vtkBooleanMacro(LimitSearchByNeighbourCount, int);

  // Description:
  // Set the Maximum number of Neighbours to use in SHEPARD mode
  vtkSetMacro(MaximumNeighbours, int);
  vtkGetMacro(MaximumNeighbours, int);
  
  // Description:
  // Set the Maximum radius to use in SHEPARD mode
  vtkSetMacro(MaximumRadius, double);
  vtkGetMacro(MaximumRadius, double);
  
protected:
   vtkSPHManager();
  ~vtkSPHManager();

  // switch shepard/sph
  int    InterpolationMethod;

  // Shepard Mode
  int    LimitSearchByNeighbourCount;
  int    MaximumNeighbours;
  double MaximumRadius;

  // SPH Mode
  int                KernelType;
  int                KernelDimension;
  char              *HScalarsRegex;
  char              *VolumeScalarsRegex;
  char              *MassScalarsRegex;
  char              *DensityScalarsRegex;

  // SPH Variables
  double             DefaultParticleSideLength;
  double             DefaultDensity;
  double             HCoefficient;
  //

  static vtkSPHManager *SPHManagerSingleton;
private:
  vtkSPHManager(const vtkSPHManager&);  // Not implemented.
  void operator=(const vtkSPHManager&);  // Not implemented.
};

#endif
