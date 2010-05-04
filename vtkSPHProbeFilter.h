/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSPHProbeFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSPHProbeFilter - Probe Mesless or Mesh based data with a common interface
// .SECTION Description
// vtkSPHProbeFilter is a filter that takes two inputs, one is a dataset
// which will be probed, the other is geometry representing probe locations.
// The filter operates on particle data, or on mesh-based data and detects
// which type of data is being operated upon. 
// If meshed data is supplied, and the probe is a planar slice,
// the filter will attempt to use clip/cut filters to preserve the original
// mesh so that the interpolation is exact (providing the original mesh
// consisted of linear cells). To activate this smart probing, the user
// should use the vtkRegularGridSource filter to generate the input plane.
// vtkRegularGridSource tags the output data with plane parameters so that
// an implicit plane and box can be created for clip/cut operations.
// If meshed data is supplied and the probe type is some arbitrary geometry
// then the filter internally uses a standard vtkProbe filter to create the
// required output.
// When meshless data is probed, the filter uses custom interpolation
// routines based on SPH kernel or Shepard interpolations to generate
// output values.

#ifndef __vtkSPHProbeFilter_h
#define __vtkSPHProbeFilter_h

#include "vtkDataSetAlgorithm.h"
#include "vtkSmartPointer.h" // for smartpointers

class vtkIdTypeArray;
class vtkPointLocator;
class vtkDataSet;
class vtkPointSet;
class vtkIdList;
class vtkTimerLog;
//BTX
class Kernel;
//ETX

#define PROBE_MAX_NEIGHBOURS 1024

class VTK_EXPORT vtkSPHProbeFilter : public vtkDataSetAlgorithm
{
public:
  static vtkSPHProbeFilter *New();
  vtkTypeRevisionMacro(vtkSPHProbeFilter,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the point locations used to probe input. Any geometry
  // can be used. Old style. Do not use unless for backwards compatibility.
  void SetProbe(vtkDataObject *probe);
  vtkDataObject *GetProbe();

  // Description:
  // Specify the point locations used to probe input. Any geometry
  // can be used.
  void SetProbeConnection(vtkAlgorithmOutput* algOutput);

  //BTX
  enum {
    POINT_INTERPOLATION_KERNEL=0,
    POINT_INTERPOLATION_SHEPARD,
  };
  //ETX

  // Description:
  // An internal method used to test if an implicit function exists 
  // on the probe input
  int GetImplicitFunctionStatus(vtkDataSet *probepts);

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
  vtkSetStringMacro(VolumeScalars);
  vtkGetStringMacro(VolumeScalars);

  // Description:
  // For a variable-h simulation, an h-array must be supplied
  // this determines the kernel cutoff on a per-particle basis.
  vtkSetStringMacro(HScalars);
  vtkGetStringMacro(HScalars);

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
  vtkSetStringMacro(DensityScalars);
  vtkGetStringMacro(DensityScalars);

  // Description:
  // Material density is 1000 (Kg/m^3) by default, but can be changed if
  // the simulation requires it, overridden if DensityScalars are supplied.
  vtkSetMacro(DefaultDensity, double);
  vtkGetMacro(DefaultDensity, double);

  // Description:
  // if each particle has a mass value, specify the array name here.
  // if not specified then the value of mass will be computed from particle size
  // and density etc.
  vtkSetStringMacro(MassScalars);
  vtkGetStringMacro(MassScalars);

  //BTX
  enum {
    SPH_KERNEL_GAUSSIAN=0,
    SPH_KERNEL_QUADRATIC,
    SPH_KERNEL_SPLINE_3RD,
    SPH_KERNEL_SPLINE_5TH,
    SPH_KERNEL_CUSP,
    SPH_KERNEL_BOX,
  };
  //ETX

  // Description:
  // Set which kernel is being used, currently the following are supported
  // CubicSPline3D, CubicSPline2D, Cusp3D
  vtkSetMacro(KernelType, int);
  vtkGetMacro(KernelType, int);
  void SetKernelTypeToGaussian(int k) { 
    this->SetKernelType(SPH_KERNEL_GAUSSIAN); }
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
   vtkSPHProbeFilter();
  ~vtkSPHProbeFilter();

  int SpatialMatch;

  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//BTX
  // create the correct output type when looping over blocks
  // with multiblock probes
  vtkSmartPointer<vtkDataSet> NewOutput(vtkDataSet *probepts);
  virtual int OutputType(vtkDataSet *probepts);
//ETX
  // main execution of loop over probe points
  bool ProbeMeshless(vtkPointSet *data, vtkPointSet *probepts, vtkDataSet *output);

  // kernel specific functions
  void   InitializeKernelCoefficients();
  double GetMaxKernelCutoffDistance();
  void   KernelCompute(double x[3], vtkPointSet *source, vtkIdList *NearestPoints, double *gradW);

  //
  // Variables
  //
//BTX
  vtkSmartPointer<vtkPointLocator>   Locator;
  vtkSmartPointer<vtkTimerLog>       Timer;
//ETX

  // switch shepard/sph
  int    InterpolationMethod;

  // Shepard Mode
  int    LimitSearchByNeighbourCount;
  int    MaximumNeighbours;
  double MaximumRadius;

  // SPH Mode
  int                KernelType;
  int                KernelDimension;
  char              *HScalars;
  char              *VolumeScalars;
  char              *MassScalars;
  char              *DensityScalars;
  vtkDataArray      *MassArray;
  vtkDataArray      *DensityArray;
  vtkDataArray      *VolumeArray;
  vtkDataArray      *HArray;

  // SPH Variables
  double             DefaultParticleSideLength;
  double             DefaultParticleVolume;
  double             DefaultParticleMass;
  double             DefaultDensity;
  double             HCoefficient;
  double             weights[PROBE_MAX_NEIGHBOURS];
  float             *HData;
  float             *MassData;
  float             *DensityData;
  float             *VolumeData;
  //
  double             ScaleCoefficient;
  double             GradientMagnitude;
  vtkIdType          NumInputParticles;
  vtkIdType          NumOutputPoints;
//BTX
  Kernel            *KernelFunction;
//ETX

private:
  vtkSPHProbeFilter(const vtkSPHProbeFilter&);  // Not implemented.
  void operator=(const vtkSPHProbeFilter&);  // Not implemented.
};

#endif
