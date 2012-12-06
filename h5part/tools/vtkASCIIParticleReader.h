/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkASCIIParticleReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkASCIIParticleReader - 
//

#ifndef _vtkASCIIParticleReader_h
#define _vtkASCIIParticleReader_h

#include "vtkCSVBaseReader.h"
#include <utility>
#include <vector>
#include <string>

//BTX
#ifdef CSCS_PARAVIEW_INTERNAL
  #define VTK_vtkCSCSCommon_EXPORT CSCS_EXPORT
#else
  #include "vtkCSCSCommonConfigure.h"
#endif
//ETX

#include <vector>
#include "vtkSmartPointer.h"
class vtkPolyData;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkIdList;

typedef std::vector< vtkSmartPointer<vtkFloatArray> > ScalarList;

class vtkASCIIParticleReader : public vtkCSVBaseReader {
  public:
    static vtkASCIIParticleReader *New();

    // Description:
    // Standard Type-Macro
    vtkTypeMacro(vtkASCIIParticleReader,vtkCSVBaseReader);

    // Description:
    // Set/Get the Time step to read
    vtkSetMacro(TimeStep, int);
    vtkGetMacro(TimeStep, int);

    // Description:
    // Pass in a string of substrings to tell the reader to combine different 
    // scalars into vector fields. 
    // For example to tell the reader to use
    // columns 0 1 2 as position x,y,z : 3,4,5 as velocity x y z and 6,7,8 as vorticity
    // use reader->SetFieldIndices("0 1 2, 3 4 5, 6 7 8") - 
    // if there is no Z coordinate use -1 as the index
    // space delimits columns and comma delimits fields.
    vtkSetStringMacro(FieldIndices);
    vtkGetStringMacro(FieldIndices);

    // Description:
    // Pass in a string of substrings to supply names for arrays
    // Example fields 
    // 0,1,2 position (x,y,z): 
    // 3,4,5 velocity (x,y,z):
    // 6     mass 
    // 7     density
    // 8,9,10 vorticity
    // ->SetFieldIndices("0 1 2, 3 4 5, 6, 7, 8 9 10")
    // ->SetFieldNames("Position, Velocity, mass, density, vorticity")
    vtkSetStringMacro(FieldNames);
    vtkGetStringMacro(FieldNames);
    
    // Description:
    vtkSetMacro(DFormat,int);
    vtkGetMacro(DFormat,int);
    vtkBooleanMacro(DFormat,int);

    // Description:
    vtkSetMacro(MultiFileScalarMode,int);
    vtkGetMacro(MultiFileScalarMode,int);
    vtkBooleanMacro(MultiFileScalarMode,int);    

    // Description:
    vtkSetMacro(MultiFileCollectNode,int);
    vtkGetMacro(MultiFileCollectNode,int);
    vtkBooleanMacro(MultiFileCollectNode,int);    

    // Description:
    vtkSetMacro(TimeValue,double);
    vtkGetMacro(TimeValue,double);

    // Description:
    vtkSetStringMacro(TimeExpression);
    vtkGetStringMacro(TimeExpression);

    // Description:
    vtkSetStringMacro(IgnoreExpression);
    vtkGetStringMacro(IgnoreExpression);

    // Description:
    // Save the range of valid timestep index values. This can be used by the ParaView GUI
    // The TimeStepRange will only be meaningful if the InterpolationTimeStep
    // is set to a valid number.
    // For example, if the input has 3 time steps at 0,1,2 setting the 
    // NumberOfSubDivisions to 10 will produce an output with 21 values corresponding
    // to 0.0, 0.1, 0.2,.....1.8, 1.9, 1.0 
    // These values can be accessed as TimeStep 0 to TimeStep 21
    vtkGetVector2Macro(TimeStepRange, int);
    
    void GetTimeStepValues(std::vector<double> &values);

    vtkPolyData *GetPolyDataOutput();

    void CollectMultiFileScalars(ScalarList &scalars);

  protected:
    std::vector<long int>  FileMarkers;
    std::vector<int>       RowsPerStep;
    std::vector<double>    TimeStepValues;
    int                    DFormat;
    int                    MultiFileScalarMode;
    int                    MultiFileCollectNode;
    double                 TimeValue;
    int                    InvalidTime;
    std::string         LastFileRead;
    char                  *FieldIndices;
    char                  *FieldNames;
    char                  *TimeExpression;
    char                  *IgnoreExpression;
    vtkPoints             *Points;
    vtkCellArray          *Cells;
    vtkFloatArray        **Scalars;
    typedef std::pair<std::string, std::vector<int> > ScalarPair;
    std::vector< ScalarPair > ScalarNamesList;

    //
    virtual bool ParseFields(vtkInformationVector* outputVector) ;
    virtual bool AllocateOutputData(vtkInformationVector* outputVector, int numrows, int numcols, int &actualcolumns) ;
    virtual bool ProcessOneLine(vtkInformationVector* outputVector, vtkIdType rownum, const char *line);
    virtual void ProcessOneRow(vtkInformationVector* outputVector, vtkIdType rownum, double *onerow);
    //
    virtual void CleanUp(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector);
    bool ReadMetaInformation();
    bool WriteMetaInformation();

    //
    // Fill output information 
    //
    virtual int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info);
    //
    // Get information about the data
    //
    virtual int RequestInformation(vtkInformation* request,
                                   vtkInformationVector** inputVector,
                                   vtkInformationVector* outputVector);
    //
    // Read the data
    //
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);
  protected:
     vtkASCIIParticleReader();
    ~vtkASCIIParticleReader();
    //
    friend class SumSizes;

private:
  vtkASCIIParticleReader(const vtkASCIIParticleReader&);  // Not implemented.
  void operator=(const vtkASCIIParticleReader&);  // Not implemented.
};

#endif
