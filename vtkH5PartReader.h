/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5PartReader.h
  Revision of last commit : $Rev: 754 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2009-01-09 13:40:38 +0100 #$

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
// .NAME vtkH5PartReader - Write H5Part (HDF5) Particle files
// .SECTION Description
// vtkH5PartReader reads compatible with H5Part : documented here
// http://amas.web.psi.ch/docs/H5Part-doc/h5part.html 

#ifndef __vtkH5PartReader_h
#define __vtkH5PartReader_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkBoundingBox.h"
#include <vtkstd/string>
#include <vtkstd/vector>

class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkBoundsExtentTranslator;

struct H5PartFile;

class VTK_EXPORT vtkH5PartReader : public vtkPolyDataAlgorithm
{
public:
  static vtkH5PartReader *New();
  vtkTypeMacro(vtkH5PartReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);   

  // Description:
  // Specify file name.
  void SetFileName(char *filename);
  vtkGetStringMacro(FileName);

  // Description:
  // Set/Get the array that will be used for the X coordinates
  vtkGetStringMacro(Xarray);
  vtkSetStringMacro(Xarray);

  // Description:
  // Set/Get the array that will be used for the Y coordinates
  vtkGetStringMacro(Yarray);
  vtkSetStringMacro(Yarray);

  // Description:
  // Set/Get the array that will be used for the Z coordinates
  vtkGetStringMacro(Zarray);
  vtkSetStringMacro(Zarray);

  // Description:
  // Set/Get the name used for each time step (Usually Step#0, Step#1 etc)
  vtkGetStringMacro(StepName);
  vtkSetStringMacro(StepName);

  // Description:
  // Set/Get the timestep to be read
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);
  
  // Description:
  // Export time values as 0,1...N-1 regardless of real time values in file
  vtkSetMacro(IntegerTimeStepValues,int);
  vtkGetMacro(IntegerTimeStepValues,int);
  vtkBooleanMacro(IntegerTimeStepValues,int);
  
  // Description:
  // Get the number of timesteps in the file
  vtkGetMacro(NumberOfTimeSteps,int);

  // Description:
  // When set (default no), the reader will generate a vertex cell
  // for each point/particle read. When using the points directly
  // this is unnecessary and time can be saved by omitting cell generation
  // vtkPointSpriteMapper does not require them.
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable
  vtkSetMacro(GenerateVertexCells, int);
  vtkGetMacro(GenerateVertexCells, int);
  vtkBooleanMacro(GenerateVertexCells, int);

  // Description:
  // If set, the reader will not attempt to read or use Bounding Boxes
  // for partitions, display of them will also be disabled
  vtkSetMacro(IgnorePartitionBoxes,int);
  vtkGetMacro(IgnorePartitionBoxes,int);
  vtkBooleanMacro(IgnorePartitionBoxes,int);
  
  // Description:
  // If set and present in the file, bounding boxes of each parallel
  // partition will be displayed (as lines) along with particles.
  // Note that the Partition Boxes are the ones present in the file
  // and the Piece Boxes are the ones we export from the reader
  // which will consist of 1 or more partitions joined together.
  vtkSetMacro(DisplayPartitionBoxes,int);
  vtkGetMacro(DisplayPartitionBoxes,int);
  vtkBooleanMacro(DisplayPartitionBoxes,int);
  
  // Description:
  // If set and present in the file, bounding boxes of each parallel
  // piece will be exported (as lines) along with particles.
  // Note that the Partition Boxes are the ones present in the file
  // and the Piece Boxes are the ones we export from the reader
  // which will consist of 1 or more partitions joined together.
  vtkSetMacro(DisplayPieceBoxes,int);
  vtkGetMacro(DisplayPieceBoxes,int);
  vtkBooleanMacro(DisplayPieceBoxes,int);
  
  // Description:
  // For testing, randomize extents when partitioning by pieces.
  // This is intended to permute slightly the start and end particle indices
  // or each partition so that some randomness is introduced into reads.
  // For debugging only.
  vtkSetMacro(RandomizePartitionExtents,int);
  vtkGetMacro(RandomizePartitionExtents,int);
  vtkBooleanMacro(RandomizePartitionExtents,int);

  // Description:
  // When this option is set, scalar fields with names which form a pattern
  // of the form scalar_0, scalar_1, scalar_2 will be combined into a single
  // vector field with N components
  vtkSetMacro(CombineVectorComponents, int);
  vtkGetMacro(CombineVectorComponents, int);
  vtkBooleanMacro(CombineVectorComponents, int);

  // Description:
  // When 3 separate arrays such as {x, y,,z} are loaded into a coordinate/vector array, they
  // can also be added to the field data for visualization. Default is off as it uses extra memory.
  vtkSetMacro(MultiComponentArraysAsFieldData, int);
  vtkGetMacro(MultiComponentArraysAsFieldData, int);
  vtkBooleanMacro(MultiComponentArraysAsFieldData, int);

  // Description:
  // When on UseStridedMultiComponentRead tells the reader to use  a hyperslab
  // with strides to read an array such as "X" into a coordinate array.
  // This can be very slow, but save some memory and is off by default.
  // When off, arrays such as "X, "Y", "Z" are read independently, then
  // copied into the vector array one by one.
  vtkSetMacro(UseStridedMultiComponentRead, int);
  vtkGetMacro(UseStridedMultiComponentRead, int);
  vtkBooleanMacro(UseStridedMultiComponentRead, int);

  // Description:
  // Normally, a request for data at time t=x, where x is either before the start of
  // time for the data, or after the end, will result in the first or last
  // timestep of data to be retrieved (time is clamped to max/min values).
  // Forsome applications/animations, it may be desirable to not display data
  // for invalid times. When MaskOutOfTimeRangeOutput is set to ON, the reader
  // will return an empty dataset for out of range requests. This helps
  // avoid corruption of animations.
  vtkSetMacro(MaskOutOfTimeRangeOutput, int);
  vtkGetMacro(MaskOutOfTimeRangeOutput, int);
  vtkBooleanMacro(MaskOutOfTimeRangeOutput, int);

  bool HasStep(int Step);

  // Description:
  // An H5Part file may contain multiple arrays
  // a GUI (eg Paraview) can provide a mechanism for selecting which data arrays
  // are to be read from the file. The PointArray variables and members can
  // be used to query the names and number of arrays available
  // and set the status (on/off) for each array, thereby controlling which
  // should be read from the file. Paraview queries these point arrays after
  // the (update) information part of the pipeline has been updated, and before the
  // (update) data part is updated.
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAll();
  void        EnableAll();
  void        Disable(const char* name);
  void        Enable(const char* name);
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

  int         GetNumberOfCoordinateArrays() { return GetNumberOfPointArrays(); }
  const char* GetCoordinateArrayName(int index) { return GetPointArrayName(index); }
  int         GetCoordinateArrayStatus(const char* name);
  void        SetCoordinateArrayStatus(const char* name, int status);

  // Description:
  // Set/Get the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  void SetFileModified();

protected:
   vtkH5PartReader();
  ~vtkH5PartReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // read bboxes information into box/partition arrays
  vtkIdType ReadBoundingBoxes();

  // create polygons/lines from boxes for visual display
  // Note : we only colour particles actually read on this process 
  // i.e. those between extent0 and extent1
  vtkIdType DisplayBoundingBoxes(vtkDataArray *coords, vtkPolyData *output, vtkIdType extent0, vtkIdType extent1);

  // Use simple 1...N, divide particle list to partition data for reading
  int PartitionByExtents(vtkIdType N, std::vector<vtkIdType> &startend);
  int PartitionByExtentsRandomized(vtkIdType N, std::vector<vtkIdType> &startend);
  
  // Use the BoundingBoxes read to partition data for reading
  int PartitionByBoundingBoxes(
    std::vector<vtkIdType> &minIds, 
    std::vector<vtkIdType> &maxIds,
    std::vector<vtkBoundingBox> &PieceBounds,
    std::vector<vtkBoundingBox> &PieceHaloBounds);

  //
  virtual int  OpenFile();
  virtual void CloseFile();
  
  // Under normal circumstances, when loading data and animating though timesteps
  // one does not want to close the file between steps (calls to ExecuteInfo/Data)
  // but subclasses (e.g. dsm) do need to close the file for real. We therefore
  // call CloseFileIntermediate when we can leave it open and subclasses can decide
  // to act on it, or do nothing. By default, do nothing.
  virtual void CloseFileIntermediate();

  // Copy a scalar array into one component field of a vector dataset
  // example X using offset 0, into {X,-,-}, etc etc
  template <class T1>
  void CopyIntoVector(int offset, vtkDataArray *source, vtkDataArray *dest);

  // returns 0 if no, returns 1,2,3,45 etc for the first, second...
  // example : if CombineVectorComponents is true, then 
  // velocity_0 returns 1, velocity_1 returns 2 etc
  // if CombineVectorComponents is false, then 
  // velocity_0 returns 0, velocity_1 returns 0 etc
  int             IndexOfVectorComponent(const char *name);
//BTX
  vtkstd::string  NameOfVectorComponent(const char *name);
//ETX

  //
  // Internal Variables
  //
  char         *FileName;
  int           NumberOfTimeSteps;
  int           TimeStep;
  int           ActualTimeStep;
  double        TimeStepTolerance;
  int           CombineVectorComponents;
  int           MultiComponentArraysAsFieldData;
  int           UseStridedMultiComponentRead;
  int           GenerateVertexCells;
  H5PartFile   *H5FileId;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           MaskOutOfTimeRangeOutput;
  int           TimeOutOfRange;
  int           IntegerTimeStepValues;
  int           IgnorePartitionBoxes;
  int           DisplayPartitionBoxes;
  int           DisplayPieceBoxes;
  int           UseLinearBoxPartitioning;
  int           RandomizePartitionExtents;
  //
  char         *Xarray;
  char         *Yarray;
  char         *Zarray;
  char         *StepName;
//BTX
  vtkstd::vector<double>                  TimeStepValues;
  typedef vtkstd::vector<vtkstd::string>  stringlist;
  vtkstd::vector<stringlist>              FieldArrays;
  // For Bounding boxes if present
  std::vector<vtkIdType>      PartitionCount;
  std::vector<vtkIdType>      PartitionOffset;
  std::vector<vtkIdType>      PieceId;
  std::vector<double>         PartitionBoundsTable;
  std::vector<double>         PartitionBoundsTableHalo;
  std::vector<vtkBoundingBox> PieceBounds;
  std::vector<vtkBoundingBox> PieceBoundsHalo;
  vtkBoundsExtentTranslator  *ExtentTranslator;
//ETX

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* CoordinateSelection;

  vtkMultiProcessController* Controller;

private:
  vtkH5PartReader(const vtkH5PartReader&);  // Not implemented.
  void operator=(const vtkH5PartReader&);  // Not implemented.
};

#endif
