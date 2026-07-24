// SPDX-FileCopyrightText: Copyright (C) CSCS - Swiss National Supercomputing
// Centre SPDX-License-Identifier: LicenseRef-CSCS

#ifndef vtkH5PartReaderV2_h
#define vtkH5PartReaderV2_h

#include "vtkBoundingBox.h"
#include "vtkPolyDataAlgorithm.h"

#include <string>
#include <vector>

class vtkBoundsExtentTranslator;
class vtkDataArraySelection;
class vtkMultiProcessController;

class VTK_EXPORT vtkH5PartReaderV2 : public vtkPolyDataAlgorithm {
public:
  static vtkH5PartReaderV2 *New();
  vtkTypeMacro(vtkH5PartReaderV2, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  ///@{
  /**
   * Specify file name.
   */
  void SetFileName(VTK_FILEPATH char *filename);
  vtkGetFilePathMacro(FileName);
  ///@}

  ///@{
  /**
   * Set/Get the array that will be used for the X coordinates
   */
  vtkGetStringMacro(Xarray);
  vtkSetStringMacro(Xarray);
  ///@}

  ///@{
  /**
   * Set/Get the array that will be used for the Y coordinates
   */
  vtkGetStringMacro(Yarray);
  vtkSetStringMacro(Yarray);
  ///@}

  ///@{
  /**
   * Set/Get the array that will be used for the Z coordinates
   */
  vtkGetStringMacro(Zarray);
  vtkSetStringMacro(Zarray);
  ///@}

  ///@{
  /**
   * Set/Get the name used for each time step (Usually Step#0, Step#1 etc)
   */
  vtkGetStringMacro(StepName);
  vtkSetStringMacro(StepName);
  ///@}

  ///@{
  /**
   * Set/Get the timestep to be read
   */
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
  ///@}

  ///@{
  /**
   * Get the number of timesteps in the file
   */
  vtkGetMacro(NumberOfTimeSteps, int);
  ///@}

  ///@{
  /**
   * Export time values as 0,1...N-1 regardless of real time values in file
   */
  vtkSetMacro(IntegerTimeStepValues, int);
  vtkGetMacro(IntegerTimeStepValues, int);
  vtkBooleanMacro(IntegerTimeStepValues, int);
  ///@}

  ///@{
  /**
   * When set (default no), the reader will generate a vertex cell
   * for each point/particle read.
   */
  vtkSetMacro(GenerateVertexCells, int);
  vtkGetMacro(GenerateVertexCells, int);
  vtkBooleanMacro(GenerateVertexCells, int);
  ///@}

  ///@{
  /**
   * If set, the reader will not attempt to read or use Bounding Boxes
   * for partitions, display of them will also be disabled
   */
  vtkSetMacro(IgnorePartitionBoxes, int);
  vtkGetMacro(IgnorePartitionBoxes, int);
  vtkBooleanMacro(IgnorePartitionBoxes, int);
  ///@}

  ///@{
  /**
   * If set and present in the file, bounding boxes of each parallel
   * partition will be displayed (as lines) along with particles.
   */
  vtkSetMacro(DisplayPartitionBoxes, int);
  vtkGetMacro(DisplayPartitionBoxes, int);
  vtkBooleanMacro(DisplayPartitionBoxes, int);
  ///@}

  ///@{
  /**
   * If set and present in the file, bounding boxes of each parallel
   * piece will be exported (as lines) along with particles.
   */
  vtkSetMacro(DisplayPieceBoxes, int);
  vtkGetMacro(DisplayPieceBoxes, int);
  vtkBooleanMacro(DisplayPieceBoxes, int);
  ///@}

  ///@{
  /**
   * For testing, randomize extents when partitioning by pieces.
   */
  vtkSetMacro(RandomizePartitionExtents, int);
  vtkGetMacro(RandomizePartitionExtents, int);
  vtkBooleanMacro(RandomizePartitionExtents, int);
  ///@}

  ///@{
  /**
   * When set, scalar fields with names scalar_0, scalar_1, scalar_2
   * will be combined into a single vector field.
   */
  vtkSetMacro(CombineVectorComponents, int);
  vtkGetMacro(CombineVectorComponents, int);
  vtkBooleanMacro(CombineVectorComponents, int);
  ///@}

  ///@{
  /**
   * When 3 separate arrays such as {x, y, z} are loaded into a
   * coordinate/vector array, they can also be added to the field data.
   */
  vtkSetMacro(MultiComponentArraysAsFieldData, int);
  vtkGetMacro(MultiComponentArraysAsFieldData, int);
  vtkBooleanMacro(MultiComponentArraysAsFieldData, int);
  ///@}

  ///@{
  /**
   * When on, use a strided hyperslab to read an array such as "X" into
   * a coordinate array.  Off by default.
   */
  vtkSetMacro(UseStridedMultiComponentRead, int);
  vtkGetMacro(UseStridedMultiComponentRead, int);
  vtkBooleanMacro(UseStridedMultiComponentRead, int);
  ///@}

  ///@{
  /**
   * Maximum number of particles to read per rank.  0 means unlimited
   * (read the full partition).  When > 0, each rank reads only the
   * first MaxParticlesPerRank particles from its partition, providing
   * a spatially distributed sample of the full dataset for interactive
   * exploration of very large files.
   */
  vtkSetMacro(MaxParticlesPerRank, vtkIdType);
  vtkGetMacro(MaxParticlesPerRank, vtkIdType);
  ///@}

  ///@{
  /**
   * Return true if the requested step exists in the file.
   */
  bool HasStep(int step);
  ///@}

  ///@{
  /**
   * When this option is set, a request for data outside the file time range
   * returns an empty dataset instead of clamping to the nearest valid time.
   */
  vtkSetMacro(MaskOutOfTimeRangeOutput, int);
  vtkGetMacro(MaskOutOfTimeRangeOutput, int);
  vtkBooleanMacro(MaskOutOfTimeRangeOutput, int);
  ///@}

  ///@{
  /**
   * Point array selection API used by ParaView.
   */
  int GetNumberOfPointArrays();
  const char *GetPointArrayName(int index);
  int GetPointArrayStatus(const char *name);
  void SetPointArrayStatus(const char *name, int status);
  void DisableAll();
  void EnableAll();
  void Disable(const char *name);
  void Enable(const char *name);

  int GetNumberOfPointArrayStatusArrays() {
    return this->GetNumberOfPointArrays();
  }
  const char *GetPointArrayStatusArrayName(int index) {
    return this->GetPointArrayName(index);
  }
  int GetPointArrayStatusArrayStatus(const char *name) {
    return this->GetPointArrayStatus(name);
  }
  void SetPointArrayStatusArrayStatus(const char *name, int status) {
    this->SetPointArrayStatus(name, status);
  }
  ///@}

  ///@{
  /**
   * Coordinate array selection API used by ParaView.
   */
  int GetNumberOfCoordinateArrays() { return this->GetNumberOfPointArrays(); }
  const char *GetCoordinateArrayName(int index) {
    return this->GetPointArrayName(index);
  }
  int GetCoordinateArrayStatus(const char *name);
  void SetCoordinateArrayStatus(const char *name, int status);
  ///@}

  ///@{
  /**
   * Set/Get the controller use in parallel operations.
   */
  virtual void SetController(vtkMultiProcessController *controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  ///@}

  /**
   * Mark the file as modified. Used internally when the filename changes.
   */
  void SetFileModified();

protected:
  vtkH5PartReaderV2();
  ~vtkH5PartReaderV2() override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int OpenFile();
  void CloseFile();
  void CloseFileIntermediate();

  int IndexOfVectorComponent(const char *name);
  std::string NameOfVectorComponent(const char *name);

  vtkIdType ReadBoundingBoxes();
  vtkIdType DisplayBoundingBoxes(vtkDataArray *coords, vtkPolyData *output,
                                 vtkIdType extent0, vtkIdType extent1);
  int PartitionByBoundingBoxes(int numPieces, std::vector<vtkIdType> &minIds,
                               std::vector<vtkIdType> &maxIds,
                               std::vector<vtkBoundingBox> &pieceBounds,
                               std::vector<vtkBoundingBox> &pieceHaloBounds);
  int PartitionByExtents(vtkIdType n, int piece, int numPieces,
                         std::vector<vtkIdType> &startend);
  int PartitionByExtentsRandomized(vtkIdType n, int piece, int numPieces,
                                   std::vector<vtkIdType> &startend);
  int SplitExtent(int piece, int numPieces, vtkIdType *ext);

  template <class T1>
  void CopyIntoVector(int offset, vtkDataArray *source, vtkDataArray *dest);

  char *FileName;
  int NumberOfTimeSteps;
  int TimeStep;
  int ActualTimeStep;
  double TimeStepTolerance;
  int CombineVectorComponents;
  int MultiComponentArraysAsFieldData;
  int UseStridedMultiComponentRead;
  vtkIdType MaxParticlesPerRank;
  int GenerateVertexCells;
  uintptr_t H5FileId;
  vtkTimeStamp FileModifiedTime;
  vtkTimeStamp FileOpenedTime;
  int MaskOutOfTimeRangeOutput;
  int TimeOutOfRange;
  int IntegerTimeStepValues;
  int IgnorePartitionBoxes;
  int DisplayPartitionBoxes;
  int DisplayPieceBoxes;
  int UseLinearBoxPartitioning;
  int RandomizePartitionExtents;

  char *Xarray;
  char *Yarray;
  char *Zarray;
  char *StepName;

  std::vector<double> TimeStepValues;

  // Bounding-box tables read from the file (if present).
  std::vector<vtkIdType> PartitionCount;
  std::vector<vtkIdType> PartitionOffset;
  std::vector<vtkIdType> PieceId;
  std::vector<double> PartitionBoundsTable;
  std::vector<double> PartitionBoundsTableHalo;
  std::vector<vtkBoundingBox> PieceBounds;
  std::vector<vtkBoundingBox> PieceBoundsHalo;
  vtkBoundsExtentTranslator *ExtentTranslator;

  vtkDataArraySelection *PointDataArraySelection;

  vtkMultiProcessController *Controller;

private:
  vtkH5PartReaderV2(const vtkH5PartReaderV2 &) = delete;
  void operator=(const vtkH5PartReaderV2 &) = delete;
};

#endif
