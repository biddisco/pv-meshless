/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkSamplingGridGenerator.h,v $
  Date:      $Date: 2004/02/06 07:31:49 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
// .NAME vtkSamplingGridGenerator - Generate points in a box-like volume
//
// .SECTION Description
// vtkSamplingGridGenerator creates points along a regular grid
// in3 dimensions. The grid is defined by the X/Y/Z Spacing fields
// and a supplied bounding volume. In general the bounding volume
// is taken from a supplied input dataset.
//
// .SECTION See Also
// vtkGeneratePointsFilter

#ifndef _vtkSamplingGridGenerator_h
#define _vtkSamplingGridGenerator_h

#include "vtkRegularGridSource.h"

class vtkImplicitFunction;
class vtkImageData;
class vtkCutter;

class VTK_EXPORT vtkSamplingGridGenerator : public vtkRegularGridSource {
  public:
    // Description:
    // Standard Type-Macro
    static vtkSamplingGridGenerator *New();
    vtkTypeRevisionMacro(vtkSamplingGridGenerator, vtkDataSetAlgorithm);

    // Description
    // Specify the implicit function to perform the cutting.
    virtual void SetCutFunction(vtkImplicitFunction*);
    vtkGetObjectMacro(CutFunction,vtkImplicitFunction);

  protected:
     vtkSamplingGridGenerator();
    ~vtkSamplingGridGenerator(); 

    virtual int RequiredDataType();

    virtual int ComputeInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector *);

    //  
    // Implicit Function for cutting
    //
    vtkImplicitFunction *CutFunction;
    vtkImageData        *Box;
    vtkCutter           *Cutter;

private:
  vtkSamplingGridGenerator(const vtkSamplingGridGenerator&);  // Not implemented.
  void operator=(const vtkSamplingGridGenerator&);  // Not implemented.
};

#endif
