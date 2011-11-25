/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkSPHImageResampler.h,v $
  Date:      $Date: 2004/02/06 07:31:49 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
// .NAME vtkSPHImageResampler - Generate points in a box-like volume
//
// .SECTION Description
// vtkSPHImageResampler creates points along a regular grid
// in3 dimensions. The grid is defined by the X/Y/Z Spacing fields
// and a supplied bounding volume. In general the bounding volume
// is taken from a supplied input dataset.
//
// .SECTION See Also
// vtkGeneratePointsFilter

#ifndef _vtkSPHImageResampler_h
#define _vtkSPHImageResampler_h

#include "vtkToolkits.h"     // For VTK_USE_MPI 
#include "vtkImageAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkSPHManager;
class vtkSPHProbeFilter;
class vtkMultiProcessController;


class VTK_EXPORT vtkSPHImageResampler : public vtkImageAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    static vtkSPHImageResampler *New();
    vtkTypeRevisionMacro(vtkSPHImageResampler, vtkImageAlgorithm);

    // Description:
    // Specify the point spacing on the X/Y/Z axis
    vtkSetVector3Macro(Spacing, double);
    vtkGetVector3Macro(Spacing, double);

    // Description
    // If Resolution[X/Y/Z]are all Non-zero, then
    // the spacing is ignored and the box defined by the points is 
    // sampled using the specified resolutions
    vtkSetVector3Macro(Resolution, int);
    vtkGetVector3Macro(Resolution, int);

    // Description
    // Add a delta to each axis to expand the box
    vtkSetMacro(Delta, double);
    vtkGetMacro(Delta, double);

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

//BTX
    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);
//ETX

  protected:
    vtkSPHImageResampler();

    // Pipeline mechanism
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    // Description:
    // This is called within ProcessRequest when a request asks for 
    // Information. Typically an algorithm provides whatever lightweight 
    // information about its output that it can here without doing any 
    // lengthy computations. This happens in the first pass of the pipeline
    // execution.
    virtual int RequestInformation(vtkInformation*, 
                                   vtkInformationVector**, 
                                   vtkInformationVector*);
    
    //// Description:
    virtual int RequestUpdateExtent(vtkInformation*,
                                    vtkInformationVector**,
                                    vtkInformationVector*);
    
    // Description:
    // This is called within ProcessRequest when a request asks the algorithm
    // to do its work. This is the method you should override to do whatever the
    // algorithm is designed to do. This happens during the fourth pass in the
    // pipeline execution process.
    virtual int RequestData(vtkInformation*, 
                            vtkInformationVector**, 
                            vtkInformationVector*);

    virtual int ComputeInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector *);

    virtual void ComputeAxesFromBounds(vtkDataSet *inputData, double lengths[3], bool inflate);

    // properties
    double   GlobalOrigin[3];
    double   Spacing[3];
    int      Resolution[3];
    double   Delta;

    // internal values which may differ from above ones
    int      WholeDimension[3];
    double   scaling[3];
    double   spacing[3];

    //  
    // SPHManager
    //
    vtkSPHManager *SPHManager;
    vtkSmartPointer<vtkSPHProbeFilter> SPHProbe;
    vtkMultiProcessController         *Controller;
    //      
    // SPH Mode
    //
    char              *HScalars;
    char              *VolumeScalars;
    char              *MassScalars;
    char              *DensityScalars;
    int                ModifiedNumber;
    int                ComputeDensityFromNeighbourVolume;
    bool               BoundsInitialized;


private:
  vtkSPHImageResampler(const vtkSPHImageResampler&);  // Not implemented.
  void operator=(const vtkSPHImageResampler&);  // Not implemented.
};

#endif
