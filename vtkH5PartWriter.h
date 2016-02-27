/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5PartWriter.h
  Revision of last commit : $Rev: 153 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2006-07-12 10:09:37 +0200 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkH5PartWriter - Write H5Part (HDF5) Particle files
// .SECTION Description
// vtkH5PartWriter writes compatible with H5Part : documented here
// http://amas.web.psi.ch/docs/H5Part-doc/h5part.html

#ifndef __vtkH5PartWriter_h
#define __vtkH5PartWriter_h

#include "vtkSmartPointer.h" // For vtkSmartPointer
#include <string>     // for strings
#include <vector>     // for vectors
#include "vtkAbstractParticleWriter.h"
//
class vtkMultiProcessController;
//

struct H5PartFile;
class vtkPointSet;
class vtkDataArray;
class vtkPointData;

class VTK_EXPORT vtkH5PartWriter : public vtkAbstractParticleWriter
{
public:
  static vtkH5PartWriter *New();
  vtkTypeMacro(vtkH5PartWriter,vtkAbstractParticleWriter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the input to this writer.
  vtkPointSet* GetInput();
  vtkPointSet* GetInput(int port);

  // Description:
  // Quesry the file to see if a timestep was written previously
  // this should really be part of vtkHDF5Reader, but when the
  // Writer is used as a cache file, we need a function to check
  // if the timestep is already present
  bool IsTimeStepPresent(int timestep);
  void DeleteTimeStep(int timestep);

  // Description:
  // Make this public so that files can be closed between time steps and
  // the file might survive an application crash.
  void CloseFile();

  // Description:
  // Usually, we want to open the file in WRITE mode, but when being used
  // as a file cache, we want to write and read so use READWRITE mode
  vtkSetMacro(FileMode,int);
  vtkGetMacro(FileMode,int);
  void SetFileModeToWrite();
  void SetFileModeToReadWrite();

  // Description:
  // If VectorsWithStridedWrite is true, performance may be impacted
  // but no copy of {X,Y,Z] out of the triple vector into a single
  // flat array takes place. Test show that using a strided write
  // is terrible in parallel. Hope to fix this one day.
  vtkSetMacro(VectorsWithStridedWrite,int);
  vtkGetMacro(VectorsWithStridedWrite,int);
  vtkBooleanMacro(VectorsWithStridedWrite,int);

  // Description:
  // Set/Get the controller used for coordinating parallel writing
  // (set to the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  // Description:
  // When writing in parallel, all processes must send data to HDF...
  // if one node has no particles, it does nothing - but this causes the
  // others to hang. We therefore gather information and write an
  // empty array of the correct type
  vtkSetMacro(DisableInformationGather,int);
  vtkGetMacro(DisableInformationGather,int);
  vtkBooleanMacro(DisableInformationGather,int);

  // Description:
  // Set/Get the name used for each time step (Usually Step#0, Step#1 etc)
  vtkGetStringMacro(StepName);
  vtkSetStringMacro(StepName);

protected:
   vtkH5PartWriter();
  ~vtkH5PartWriter();
  //
  int   OpenFile();

  // Override superclass' Write method
  virtual void WriteData();

  void CopyFromVector(int offset, vtkDataArray *source, vtkDataArray *dest);
  void WriteDataArray(int i, vtkDataArray *array);

  // Overide information to only permit PolyData as input
  virtual int FillInputPortInformation(int, vtkInformation *info);
  virtual int FillOutputPortInformation(int, vtkInformation* info);

  // Description:
  virtual int RequestInformation(vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector);

  // Description:
  // If a certain process has zero particles, the dataarrays for
  // points and scalars will be effectivle empty, or even NULL
  // in this case, the collective parallel IO write may fail because
  // the zero data process does not know what datatype to 'write'
  // or dataset names to create. We therefore provide a gather call
  // before writing to ensure that all processes 'agree' on what they are writing.
//BTX
  bool GatherDataArrayInfo(vtkDataArray *data, int &datatype,
    std::string &dataname, int &numComponents);
  bool GatherScalarInfo(vtkPointData *pd, int N, int &numScalar);
//ETX

  //
  // Internal Variables
  //
  int           NumberOfTimeSteps;
  long long     NumberOfParticles;
  int           FileMode;
  int           VectorsWithStridedWrite;
  H5PartFile   *H5FileId;
  char         *StepName;
  //BTX
  std::vector<double>  InputTimeValues;
  //ETX
  int           ActualTimeStep;

  // Used for Parallel write
  int     UpdatePiece;
  int     UpdateNumPieces;
  int     DisableInformationGather;

  vtkMultiProcessController* Controller;

private:
  vtkH5PartWriter(const vtkH5PartWriter&);  // Not implemented.
  void operator=(const vtkH5PartWriter&);  // Not implemented.
};

#endif
