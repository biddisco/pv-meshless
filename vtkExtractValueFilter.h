/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkExtractValueFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkExtractValueFilter - Extract the max/min value from a dataset 
// .SECTION Description
// vtkExtractValueFilter can be used to extract the maximum or minimum value from
// a dataset cordinates, or scalars. It is particularly useful for further processing 
// the output of a probe filter for example to extract a value such surface height
// for ploitting over time. It generates a single point as output with a scalar value
// set to the extracted value.

#ifndef __vtkExtractValueFilter_h
#define __vtkExtractValueFilter_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h" // for smartpointers
//
class vtkMultiProcessController;
//

class VTK_EXPORT vtkExtractValueFilter : public vtkPolyDataAlgorithm
{
  public:
    static vtkExtractValueFilter *New();
    vtkTypeMacro(vtkExtractValueFilter,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Set the extraction algorithm to operate on the coordinates
    // instead of the scalars.
    // Default is true.
    vtkSetMacro(ExtractByMaximum,int);
    vtkGetMacro(ExtractByMaximum,int);
    vtkBooleanMacro(ExtractByMaximum,int);

    // Description:
    // Set/Get the scalar value to use for testing and extraction
    vtkSetStringMacro(ExtractionScalars);
    vtkGetStringMacro(ExtractionScalars);

    // Description:
    // Set the extraction algorithm to operate on the coordinates
    // instead of the scalars. Default is true.
    vtkSetMacro(ExtractByCoordinate,int);
    vtkGetMacro(ExtractByCoordinate,int);
    vtkBooleanMacro(ExtractByCoordinate,int);

    // Description:
    // Set the component of the array to search for, the default value is
    // 2 (so that extraction by maximum Z is the default operation {x=0,y=1,z=2})
    vtkSetMacro(Component,int);
    vtkGetMacro(Component,int);

    // Description:
    // Set/Get the controller used for coordinating parallel writing
    // (set to the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);

  protected:
     vtkExtractValueFilter();
    ~vtkExtractValueFilter();

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  //BTX
    void FindMaximum(vtkDataSet *input, vtkPolyData *output, vtkDataArray *scalars);
  //ETX

    int ExtractByMaximum;
    int ExtractByCoordinate;
    int Component;
    char *ExtractionScalars;

    vtkMultiProcessController* Controller;

  private:
    vtkExtractValueFilter(const vtkExtractValueFilter&);  // Not implemented.
    void operator=(const vtkExtractValueFilter&);  // Not implemented.
};

#endif
