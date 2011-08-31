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
// .NAME vtkSPHProbeFilter - vtkSPHProbeFilter is a abstract base class for SPH probe filters.
// .SECTION Description
// vtkSPHProbeFilter is a abstract base class for SPH probe filters.
// It contains the main logic, neighbour searching and interpolation
// routines that are used by all SPH probe filters.
// Concrete sub-classes of this may use geometry/implicit functions
// file readers etc to generate the probe locations.
// The Probe filter reads its kernel parameters from the SPH manager class.

#ifndef __vtkSPHProbeFilter_h
#define __vtkSPHProbeFilter_h

#include "vtkDataSetAlgorithm.h"
#include "vtkSmartPointer.h" // for smartpointers
#include "vtkSPHManager.h"   // for vtkSPHManager

class vtkIdTypeArray;
class vtkPointLocator;
class vtkDataSet;
class vtkDataSet;
class vtkIdList;
class vtkTimerLog;
//BTX
class Kernel;
//ETX

class VTK_EXPORT vtkSPHProbeFilter : public vtkDataSetAlgorithm
{
public:
  static vtkSPHProbeFilter *New();
  vtkTypeMacro(vtkSPHProbeFilter,vtkDataSetAlgorithm);

  // Description:
  // Set/Get the geometry used for probing/sampling
  void SetProbeConnection(vtkAlgorithmOutput* algOutput);

  // Description:
  // Set/Get the geometry used for probing/sampling
  void SetProbe(vtkDataObject *probe);
  vtkDataObject *GetProbe();

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

  // Description:
  // When enabled (and mass is present), the density for each sample point is 
  // computed by smoothing the mass and using the volume from the sphere 
  // which surrounds all N neighbours
  vtkSetMacro(ComputeDensityFromNeighbourVolume,int);
  vtkGetMacro(ComputeDensityFromNeighbourVolume,int);
  vtkBooleanMacro(ComputeDensityFromNeighbourVolume,int);

  // overridden to handle SPHManager changes
  unsigned long GetMTime();

  vtkSetMacro(ModifiedNumber,int);
  vtkGetMacro(ModifiedNumber,int);

protected:
   vtkSPHProbeFilter();
  ~vtkSPHProbeFilter();


  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

//BTX
  // create the correct output type when looping over blocks
  // with multiblock probes
  vtkSmartPointer<vtkDataSet> NewOutput(vtkDataSet *probepts);
  virtual int OutputType(vtkDataSet *probepts);
//ETX

  // Setup of variables prior to main probing routine
  bool InitializeVariables(vtkDataSet *data);

  // main execution of loop over probe points
  bool ProbeMeshless(vtkDataSet *data, vtkDataSet *probepts, vtkDataSet *output);

  // kernel specific functions
  void   InitializeKernelCoefficients();
  double GetMaxKernelCutoffDistance();
  void   KernelCompute(double x[3], vtkDataSet *source, 
    vtkIdList *NearestPoints, double *gradW, double &totalmass, double &maxDistance);

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

  // switch shepard/sph
  int    InterpolationMethod;

  // Compute Density
  int    ComputeDensityFromNeighbourVolume;

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
  float             *MassDataF;
  double            *MassDataD;
  float             *DensityData;
  float             *VolumeData;

  //
  double             ScaleCoefficient;
  double             GradientMagnitude;
  vtkIdType          NumInputParticles;
  vtkIdType          NumOutputPoints;

  int                ModifiedNumber;

  // Parallel support
  int                UpdatePiece;
  int                UpdateNumPieces;

//BTX
  Kernel            *KernelFunction;
//ETX

private:
  vtkSPHProbeFilter(const vtkSPHProbeFilter&);  // Not implemented.
  void operator=(const vtkSPHProbeFilter&);  // Not implemented.
};

#endif
