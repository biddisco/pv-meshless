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

#include "vtkSPHProbeFilter.h"

class vtkImplicitFunction;

class VTK_EXPORT vtkSPHProbeFilter3 : public vtkSPHProbeFilter
{
  public:
    static vtkSPHProbeFilter3 *New();
    vtkTypeMacro(vtkSPHProbeFilter3,vtkSPHProbeFilter);

    // Description
    // Specify the implicit function to perform the sampling.
    virtual void SetCutFunction(vtkImplicitFunction*);
    vtkGetObjectMacro(CutFunction,vtkImplicitFunction);

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

    virtual int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  //BTX
    vtkSmartPointer<vtkPointSet> GenerateProbePts(vtkDataSet *data);
  //ETX

    int RequiredDataType();

    //
    // Variables
    //

    //  
    // Implicit Function for cutting
    //
    vtkImplicitFunction *CutFunction;
    int Resolution[3];
    int GenerateConnectedCells;
    int UseKernelDistanceForSampling;
    double KernelDistanceMultiplier;

  private:
    vtkSPHProbeFilter3(const vtkSPHProbeFilter3&);  // Not implemented.
    void operator=(const vtkSPHProbeFilter3&);  // Not implemented.
};

#endif
