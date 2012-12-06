/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkImageSamplerSource.h,v $
  Date:      $Date: 2004/02/06 07:31:49 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
// .NAME vtkImageSamplerSource - Generate points in a box-like volume
//
// .SECTION Description
// vtkImageSamplerSource creates points along a regular grid
// in3 dimensions. The grid is defined by the X/Y/Z Spacing fields
// and a supplied bounding volume. In general the bounding volume
// is taken from a supplied input dataset.
//
// .SECTION See Also
// vtkGeneratePointsFilter

#ifndef _vtkImageSamplerSource_h
#define _vtkImageSamplerSource_h

#include "vtkImageAlgorithm.h"

class VTK_EXPORT vtkImageSamplerSource : public vtkImageAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    static vtkImageSamplerSource *New();
    vtkTypeMacro(vtkImageSamplerSource, vtkImageAlgorithm);

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


  protected:
    vtkImageSamplerSource();

    // Pipeline mechanism
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    // Description:
    // A utility function to help us decide what type of data to export
    virtual int RequiredDataType();

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
    // This is called within ProcessRequest to when a request asks the
    // algorithm to create empty output data objects. This typically happens
    // early on in the execution of the pipeline. The default behavior is to 
    // create an output DataSet of the same type as the input for each 
    // output port. This method can be overridden to change the output 
    // data type of an algorithm. This happens in the third pass of the 
    // pipeline execution.
    virtual int RequestDataObject(vtkInformation* request, 
                                  vtkInformationVector** inputVector, 
                                  vtkInformationVector* outputVector);
    
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
    virtual void BoundsToExtent(double *bounds, int *extent, int updatePiece);

    // properties
    double   GlobalOrigin[3];
    double   LocalOrigin[3];
    double   Spacing[3];
    int      Resolution[3];
    double   Delta;

    // internal values which may differ from above ones
    int      WholeDimension[3];
    double   scaling[3];
    double   spacing[3];

private:
  vtkImageSamplerSource(const vtkImageSamplerSource&);  // Not implemented.
  void operator=(const vtkImageSamplerSource&);  // Not implemented.
};

#endif
