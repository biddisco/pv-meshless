/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSPHProbeFilter3.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSPHProbeFilter3 - Probe Meshless data 
// .SECTION Description
// vtkSPHProbeFilter3 is a filter that takes two inputs, one is a dataset
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

#ifndef __vtkSPHProbeFilter3_h
#define __vtkSPHProbeFilter3_h

#include "vtkDataSetAlgorithm.h"
#include "vtkSmartPointer.h" // for smartpointers
#include "vtkSPHManager.h"   // for vtkSPHManager

class vtkIdTypeArray;
class vtkPointLocator;
class vtkDataSet;
class vtkPointSet;
class vtkIdList;
class vtkTimerLog;
class vtkImplicitFunction;
//BTX
class Kernel;
//ETX

class VTK_EXPORT vtkSPHProbeFilter3 : public vtkDataSetAlgorithm
{
public:
  static vtkSPHProbeFilter3 *New();
  vtkTypeMacro(vtkSPHProbeFilter3,vtkDataSetAlgorithm);

  // Description
  // Specify the implicit function to perform the cutting.
  virtual void SetCutFunction(vtkImplicitFunction*);
  vtkGetObjectMacro(CutFunction,vtkImplicitFunction);

  // Description:
  // Set/Get the SPH manager that will look after Kernel parameters
  void SetSPHManager(vtkSPHManager *SPHManager);
  vtkGetObjectMacro(SPHManager, vtkSPHManager);

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
  // if each particle has a density value, specify the array name here.
  // if not specified then the value of DefaultDensity will be used
  vtkSetStringMacro(DensityScalars);
  vtkGetStringMacro(DensityScalars);

  // Description:
  // if each particle has a mass value, specify the array name here.
  // if not specified then the value of mass will be computed from particle size
  // and density etc.
  vtkSetStringMacro(MassScalars);
  vtkGetStringMacro(MassScalars);

    // Description
    // If Resolution[X/Y/Z]are all Non-zero, then
    // the spacing is ignored and the box defined by the points is 
    // sampled using the specified resolutions
    vtkSetVector3Macro(Resolution, int);
    vtkGetVector3Macro(Resolution, int);

    // Description:
    // By default, the filter will create an array of points which fill
    // the grid - the output will be a PolyData object with vertices at 
    // each point. When GenerateConnectedCells is set, the output will
    // be a PolyData Plane in plane mode, or a StructuredGrid of 3D cells
    // in Box mode.
    vtkSetMacro(GenerateConnectedCells,int);
    vtkGetMacro(GenerateConnectedCells,int);
    vtkBooleanMacro(GenerateConnectedCells,int);

    // Description:
    // If UseKernelDistanceForSampling is enabled, the sample distance
    // will be taken from Kernel->MaxCutoffRadius/2
    vtkSetMacro(UseKernelDistanceForSampling,int);
    vtkGetMacro(UseKernelDistanceForSampling,int);
    vtkBooleanMacro(UseKernelDistanceForSampling,int);

    // Description:
    // To reduce the sample spacing, when UseKernelDistanceForSampling is
    // enabled - increase the Multiplier by some amount
    vtkSetMacro(KernelDistanceMultiplier, double);
    vtkGetMacro(KernelDistanceMultiplier, double);

    unsigned long GetMTime();

protected:
   vtkSPHProbeFilter3();
  ~vtkSPHProbeFilter3();

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

//BTX
  vtkSmartPointer<vtkPointSet> GenerateProbePts(vtkDataSet *data);
//ETX
  int RequiredDataType();

  //
  // Variables
  //
//BTX
  vtkSmartPointer<vtkPointLocator>   Locator;
  vtkSmartPointer<vtkTimerLog>       Timer;
//ETX

  //  
  // SPHManager
  //
  vtkSPHManager *SPHManager;

  //  
  // Implicit Function for cutting
  //
  vtkImplicitFunction *CutFunction;
  int Resolution[3];
  int GenerateConnectedCells;
  int UseKernelDistanceForSampling;
  double KernelDistanceMultiplier;

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
  double             weights[KERNEL_MAX_NEIGHBOURS];
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
  vtkSPHProbeFilter3(const vtkSPHProbeFilter3&);  // Not implemented.
  void operator=(const vtkSPHProbeFilter3&);  // Not implemented.
};

#endif
