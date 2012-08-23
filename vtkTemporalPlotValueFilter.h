/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTemporalPlotValueFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTemporalPlotValueFilter - Prepare data that is changing over time for an XY-plot
//
// .SECTION Description
// 
// .SECTION See Also
//
// .SECTION Thanks
// John Bidiscombe of 
// CSCS - Swiss National Supercomputing Centre
// for creating and contributing this class.

#ifndef _vtkTemporalPlotValueFilter_h
#define _vtkTemporalPlotValueFilter_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h" // for memory safety

class vtkPoints;
class vtkCellArray;
class vtkDoubleArray;
class vtkPointData;

class VTK_EXPORT vtkTemporalPlotValueFilter : public vtkPolyDataAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    static vtkTemporalPlotValueFilter *New();
    vtkTypeMacro(vtkTemporalPlotValueFilter,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
  
    // Description:
    // Flush will wipe any existing data so that traces can be restarted from
    // whatever time step is next supplied.
    void Flush();

  protected:
     vtkTemporalPlotValueFilter();
    ~vtkTemporalPlotValueFilter();

    //
    // Make sure the pipeline knows what type we expect as input
    //
    virtual int FillInputPortInformation (int port, vtkInformation* info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    // Description:
    // The necessary parts of the standard pipeline update mechanism
    virtual int RequestInformation (vtkInformation *,
                                    vtkInformationVector **,
                                    vtkInformationVector *);
    //
    virtual int RequestData(vtkInformation *request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

    // internal data variables
    int           NumberOfTimeSteps;
    int           FirstTime;
    double        LatestTime;
    //
//BTX
    vtkSmartPointer<vtkCellArray>   Vertices;
    vtkSmartPointer<vtkPoints>      Points;
    vtkSmartPointer<vtkPointData>   Values;
    vtkSmartPointer<vtkDoubleArray> TimeData;
//ETX
    //
  private:
    vtkTemporalPlotValueFilter(const vtkTemporalPlotValueFilter&);  // Not implemented.
    void operator=(const vtkTemporalPlotValueFilter&);  // Not implemented.
};

#endif
