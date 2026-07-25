// SPDX-FileCopyrightText: Copyright (C) CSCS - Swiss National Supercomputing
// Centre SPDX-License-Identifier: LicenseRef-CSCS

#include "vtkH5PartReaderV2.h"

#include "vtkAppendPolyData.h"
#include "vtkBoundingBox.h"
#include "vtkCellArray.h"
#include "vtkCharArray.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkMathUtilities.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineSource.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSOATypeFloat32Array.h"
#include "vtkSOATypeFloat64Array.h"
#include "vtkSOATypeInt32Array.h"
#include "vtkSOATypeInt64Array.h"
#include "vtkShortArray.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringScanner.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"

#include "vtkBoundsExtentTranslator.h"
#include "vtkDummyController.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"

#include "vtkH5hutHelper.h"

#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkH5PartReaderV2);

//----------------------------------------------------------------------------
static void vtkPickArray(char *&arrayPtr,
                         const std::initializer_list<const char *> &values,
                         vtkDataArraySelection *selection) {
  if (arrayPtr != nullptr && arrayPtr[0] != '\0') {
    return;
  }

  for (int cc = 0, max = selection->GetNumberOfArrays(); cc < max; ++cc) {
    const char *aname = selection->GetArrayName(cc);
    for (const char *value : values) {
      if (vtksys::SystemTools::Strucmp(aname, value) == 0) {
        arrayPtr = vtksys::SystemTools::DuplicateString(aname);
        return;
      }
    }
  }
}

//----------------------------------------------------------------------------
static int GetVTKDataType(h5_int64_t h5hut_datatype) {
  switch (h5hut_datatype) {
  case H5_FLOAT32_T:
    return VTK_TYPE_FLOAT32;
  case H5_FLOAT64_T:
    return VTK_TYPE_FLOAT64;
  case H5_INT32_T:
    return VTK_TYPE_INT32;
  case H5_INT64_T:
    return VTK_TYPE_INT64;
  default:
    return VTK_VOID;
  }
}

//----------------------------------------------------------------------------
vtkH5PartReaderV2::vtkH5PartReaderV2() {
  this->SetNumberOfInputPorts(0);

  this->NumberOfTimeSteps = 0;
  this->TimeStep = 0;
  this->ActualTimeStep = 0;
  this->TimeStepTolerance = 1E-6;
  this->CombineVectorComponents = 1;
  this->MultiComponentArraysAsFieldData = 0;
  this->UseStridedMultiComponentRead = 0;
  this->MaxParticlesPerRank = 0;
  this->UseStridedMaxParticlesPerRank = 0;
  this->GenerateVertexCells = 0;
  this->FileName = nullptr;
  this->H5FileId = 0;
  this->Xarray = nullptr;
  this->Yarray = nullptr;
  this->Zarray = nullptr;
  this->StepName = nullptr;
  this->TimeOutOfRange = 0;
  this->MaskOutOfTimeRangeOutput = 0;
  this->IntegerTimeStepValues = 0;
  this->IgnorePartitionBoxes = 0;
  this->DisplayPartitionBoxes = 0;
  this->DisplayPieceBoxes = 0;
  this->UseLinearBoxPartitioning = 1;
  this->RandomizePartitionExtents = 0;
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->ExtentTranslator = vtkBoundsExtentTranslator::New();
  this->SetXarray("Coords_0");
  this->SetYarray("Coords_1");
  this->SetZarray("Coords_2");
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == nullptr) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}

//----------------------------------------------------------------------------
vtkH5PartReaderV2::~vtkH5PartReaderV2() {
  this->SetController(nullptr);
  this->CloseFile();

  delete[] this->FileName;
  this->FileName = nullptr;
  delete[] this->Xarray;
  this->Xarray = nullptr;
  delete[] this->Yarray;
  this->Yarray = nullptr;
  delete[] this->Zarray;
  this->Zarray = nullptr;
  delete[] this->StepName;
  this->StepName = nullptr;

  if (this->PointDataArraySelection) {
    this->PointDataArraySelection->Delete();
    this->PointDataArraySelection = nullptr;
  }
  if (this->ExtentTranslator) {
    this->ExtentTranslator->Delete();
    this->ExtentTranslator = nullptr;
  }
}

//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkH5PartReaderV2, Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::SetFileName(char *filename) {
  if (this->FileName == nullptr && filename == nullptr) {
    return;
  }
  if (this->FileName && filename && (!strcmp(this->FileName, filename))) {
    return;
  }
  delete[] this->FileName;
  this->FileName = nullptr;

  if (filename) {
    this->FileName = vtksys::SystemTools::DuplicateString(filename);
    this->SetFileModified();
  }
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::SetFileModified() {
  this->FileModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::CloseFile() {
  if (this->H5FileId != 0) {
    H5CloseFile(this->H5FileId);
    this->H5FileId = 0;
  }
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::CloseFileIntermediate() {}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::OpenFile() {
  if (!this->FileName) {
    vtkErrorMacro(<< "FileName must be specified.");
    return 0;
  }

  if (this->FileModifiedTime > this->FileOpenedTime) {
    this->CloseFile();
  }

  if (!this->H5FileId) {
    this->H5FileId = H5OpenFile(this->FileName, H5_O_RDONLY, H5_PROP_DEFAULT);
    this->FileOpenedTime.Modified();
  }

  if (!this->H5FileId) {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::IndexOfVectorComponent(const char *name) {
  if (!this->CombineVectorComponents) {
    return 0;
  }
  vtksys::RegularExpression re1(".*_([0-9]+)");
  if (re1.find(name)) {
    int index = std::stoi(re1.match(1));
    return index + 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
std::string vtkH5PartReaderV2::NameOfVectorComponent(const char *name) {
  if (!this->CombineVectorComponents) {
    return name;
  }
  vtksys::RegularExpression re1("(.*)_[0-9]+");
  if (re1.find(name)) {
    return re1.match(1);
  }
  return name;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector) {
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

  if (!this->OpenFile()) {
    return 0;
  }

  this->NumberOfTimeSteps = static_cast<int>(H5GetNumSteps(this->H5FileId));
  H5SetStep(this->H5FileId, 0);
  int nds = static_cast<int>(H5PartGetNumDatasets(this->H5FileId));
  char name[512];
  for (int i = 0; i < nds; ++i) {
    H5PartGetDatasetName(this->H5FileId, i, name, 512);
    this->PointDataArraySelection->AddArray(name);
  }

  // If the coordinate arrays are not present under the default names,
  // try to auto-detect common alternatives.
  vtkPickArray(this->Xarray, {"x", "X", "Coords_0"},
               this->PointDataArraySelection);
  vtkPickArray(this->Yarray, {"y", "Y", "Coords_1"},
               this->PointDataArraySelection);
  vtkPickArray(this->Zarray, {"z", "Z", "Coords_2"},
               this->PointDataArraySelection);

  this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);
  int validTimes = 0;
  for (int i = 0; i < this->NumberOfTimeSteps; ++i) {
    H5SetStep(this->H5FileId, i);
    h5_int64_t numAttribs = H5GetNumStepAttribs(this->H5FileId);
    if (numAttribs > 0) {
      char attribName[128];
      h5_int64_t attribType = 0;
      h5_size_t attribNelem = 0;
      for (h5_int64_t a = 0; a < numAttribs; ++a) {
        h5_err_t status = H5GetStepAttribInfo(this->H5FileId, a, attribName,
                                              128, &attribType, &attribNelem);
        if (status == H5_SUCCESS && !strncmp("TimeValue", attribName, 128) &&
            attribType == H5_FLOAT64_T && attribNelem == 1) {
          status = H5ReadStepAttribFloat64(this->H5FileId, attribName,
                                           &this->TimeStepValues[i]);
          if (status == H5_SUCCESS) {
            validTimes++;
          }
        }
      }
    }
  }
  H5SetStep(this->H5FileId, 0);

  if (this->NumberOfTimeSteps == 0) {
    vtkErrorMacro("No time steps in data");
    return 0;
  }

  if (this->IntegerTimeStepValues ||
      (this->NumberOfTimeSteps > 0 && this->NumberOfTimeSteps != validTimes)) {
    for (int i = 0; i < this->NumberOfTimeSteps; ++i) {
      this->TimeStepValues[i] = i;
    }
  }
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
               this->TimeStepValues.data(),
               static_cast<int>(this->TimeStepValues.size()));
  double timeRange[2];
  timeRange[0] = this->TimeStepValues.front();
  timeRange[1] = this->TimeStepValues.back();
  if (this->TimeStepValues.size() > 1) {
    this->TimeStepTolerance =
        0.01 * (this->TimeStepValues[1] - this->TimeStepValues[0]);
  } else {
    this->TimeStepTolerance = 1E-3;
  }
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

  vtkIdType partitions =
      this->IgnorePartitionBoxes ? 0 : this->ReadBoundingBoxes();
  if (partitions == 0) {
    this->PartitionCount.clear();
    this->PartitionOffset.clear();
    this->PieceId.clear();
    this->PartitionBoundsTable.clear();
    this->PartitionBoundsTableHalo.clear();
  }

  this->CloseFileIntermediate();
  return 1;
}

//----------------------------------------------------------------------------
template <class T1, class T2>
static void CopyIntoTuple(int offset, vtkDataArray *source,
                          vtkDataArray *dest) {
  vtkIdType N = source->GetNumberOfTuples();
  T1 *sptr = static_cast<T1 *>(source->GetVoidPointer(0));
  T2 *dptr = static_cast<T2 *>(dest->WriteVoidPointer(0, N)) + offset;
  for (vtkIdType i = 0; i < N; ++i) {
    *dptr = *sptr++;
    dptr += 3;
  }
}

//----------------------------------------------------------------------------
template <class T2>
void vtkH5PartReaderV2::CopyIntoVector(int offset, vtkDataArray *source,
                                       vtkDataArray *dest) {
  switch (source->GetDataType()) {
  case VTK_CHAR:
    CopyIntoTuple<char, T2>(offset, source, dest);
    break;
  case VTK_SIGNED_CHAR:
    CopyIntoTuple<signed char, T2>(offset, source, dest);
    break;
  case VTK_UNSIGNED_CHAR:
    CopyIntoTuple<unsigned char, T2>(offset, source, dest);
    break;
  case VTK_SHORT:
    CopyIntoTuple<short, T2>(offset, source, dest);
    break;
  case VTK_UNSIGNED_SHORT:
    CopyIntoTuple<unsigned short, T2>(offset, source, dest);
    break;
  case VTK_INT:
    CopyIntoTuple<int, T2>(offset, source, dest);
    break;
  case VTK_UNSIGNED_INT:
    CopyIntoTuple<unsigned int, T2>(offset, source, dest);
    break;
  case VTK_LONG:
    CopyIntoTuple<long, T2>(offset, source, dest);
    break;
  case VTK_UNSIGNED_LONG:
    CopyIntoTuple<unsigned long, T2>(offset, source, dest);
    break;
  case VTK_LONG_LONG:
    CopyIntoTuple<long long, T2>(offset, source, dest);
    break;
  case VTK_UNSIGNED_LONG_LONG:
    CopyIntoTuple<unsigned long long, T2>(offset, source, dest);
    break;
  case VTK_FLOAT:
    CopyIntoTuple<float, T2>(offset, source, dest);
    break;
  case VTK_DOUBLE:
    CopyIntoTuple<double, T2>(offset, source, dest);
    break;
  case VTK_ID_TYPE:
    CopyIntoTuple<vtkIdType, T2>(offset, source, dest);
    break;
  default:
    vtkErrorMacro("Unexpected data type");
  }
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector) {
  using SDDP = vtkStreamingDemandDrivenPipeline;

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::GetData(outInfo);

  const int piece = outInfo->Has(SDDP::UPDATE_PIECE_NUMBER())
                        ? outInfo->Get(SDDP::UPDATE_PIECE_NUMBER())
                        : 0;
  const int numPieces = outInfo->Has(SDDP::UPDATE_NUMBER_OF_PIECES())
                            ? outInfo->Get(SDDP::UPDATE_NUMBER_OF_PIECES())
                            : 1;

  if (this->TimeStepValues.empty()) {
    return 0;
  }

  this->TimeOutOfRange = 0;
  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    double requestedTimeValue =
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    // Prefer an exact match, otherwise clamp to the nearest valid time step.
    auto exact = std::find_if(
        this->TimeStepValues.begin(), this->TimeStepValues.end(),
        [this, &requestedTimeValue](double timeStepValue) {
          return vtkMathUtilities::FuzzyCompare(timeStepValue, requestedTimeValue,
                                                this->TimeStepTolerance);
        });
    if (exact != this->TimeStepValues.end()) {
      this->ActualTimeStep =
          static_cast<int>(exact - this->TimeStepValues.begin());
    } else if (requestedTimeValue <= this->TimeStepValues.front()) {
      this->ActualTimeStep = 0;
    } else if (requestedTimeValue >= this->TimeStepValues.back()) {
      this->ActualTimeStep = static_cast<int>(this->TimeStepValues.size()) - 1;
    } else {
      auto lower = std::lower_bound(this->TimeStepValues.begin(),
                                    this->TimeStepValues.end(),
                                    requestedTimeValue);
      if (lower == this->TimeStepValues.begin()) {
        this->ActualTimeStep = 0;
      } else if (lower == this->TimeStepValues.end()) {
        this->ActualTimeStep =
            static_cast<int>(this->TimeStepValues.size()) - 1;
      } else {
        auto prev = lower - 1;
        if (std::abs(*prev - requestedTimeValue) <=
            std::abs(*lower - requestedTimeValue)) {
          this->ActualTimeStep =
              static_cast<int>(prev - this->TimeStepValues.begin());
        } else {
          this->ActualTimeStep =
              static_cast<int>(lower - this->TimeStepValues.begin());
        }
      }
    }

    if (requestedTimeValue < this->TimeStepValues.front() ||
        requestedTimeValue > this->TimeStepValues.back()) {
      this->TimeOutOfRange = 1;
    }
    output->GetInformation()->Set(
        vtkDataObject::DATA_TIME_STEP(),
        this->TimeStepValues[this->ActualTimeStep]);
  } else {
    double timevalue[1];
    unsigned int index = this->ActualTimeStep;
    if (index < this->TimeStepValues.size()) {
      timevalue[0] = this->TimeStepValues[index];
    } else {
      timevalue[0] = this->TimeStepValues[0];
    }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                                  timevalue[0]);
  }

  this->ActualTimeStep = std::max(
      0, std::min(this->ActualTimeStep,
                  static_cast<int>(this->TimeStepValues.size()) - 1));

  if (this->TimeOutOfRange && this->MaskOutOfTimeRangeOutput) {
    return 1;
  }

  if (!this->OpenFile()) {
    return 0;
  }

  H5SetStep(this->H5FileId, this->ActualTimeStep);

  // Build the field map: arrays to load, combining vector components.
  typedef std::map<std::string, std::vector<std::string>> FieldMap;
  FieldMap scalarFields;

  std::vector<std::string> coordarrays(3, "");
  int N = this->PointDataArraySelection->GetNumberOfArrays();
  for (int i = 0; i < N; ++i) {
    const char *name = this->PointDataArraySelection->GetArrayName(i);
    bool processarray = false;
    if (!vtksys::SystemTools::Strucmp(name, this->Xarray)) {
      processarray = true;
      coordarrays[0] = name;
    }
    if (!vtksys::SystemTools::Strucmp(name, this->Yarray)) {
      processarray = true;
      coordarrays[1] = name;
    }
    if (!vtksys::SystemTools::Strucmp(name, this->Zarray)) {
      processarray = true;
      coordarrays[2] = name;
    }
    if (this->PointDataArraySelection->ArrayIsEnabled(name)) {
      processarray = true;
    }
    if (!processarray) {
      continue;
    }

    int vectorcomponent = this->IndexOfVectorComponent(name);
    if (vectorcomponent > 0) {
      std::string vectorname = this->NameOfVectorComponent(name) + "_v";
      FieldMap::iterator pos = scalarFields.find(vectorname);
      if (pos == scalarFields.end()) {
        std::vector<std::string> arraylist(1, name);
        scalarFields.insert(FieldMap::value_type(vectorname, arraylist));
      } else {
        pos->second.reserve(vectorcomponent);
        pos->second.resize(
            std::max(static_cast<int>(pos->second.size()), vectorcomponent));
        pos->second[vectorcomponent - 1] = name;
      }
    } else {
      std::vector<std::string> arraylist(1, name);
      scalarFields.insert(FieldMap::value_type(name, arraylist));
    }
  }

  FieldMap::iterator coordvector = scalarFields.end();
  for (FieldMap::iterator pos = scalarFields.begin(); pos != scalarFields.end();
       ++pos) {
    if (pos->second.size() == 3 && (pos->second[0] == coordarrays[0]) &&
        (pos->second[1] == coordarrays[1]) &&
        (pos->second[2] == coordarrays[2])) {
      FieldMap::value_type element("Coords", pos->second);
      scalarFields.erase(pos);
      coordvector = scalarFields.insert(element).first;
      break;
    }
  }
  if (coordvector == scalarFields.end()) {
    scalarFields.insert(FieldMap::value_type("Coords", coordarrays));
  }

  // Determine the range of particles to read on this piece.
  H5PartSetView(this->H5FileId, -1, -1);
  vtkIdType totalParticles =
      static_cast<vtkIdType>(H5PartGetNumItems(this->H5FileId));
  vtkIdType particleStart = 0;
  vtkIdType particleEnd = totalParticles - 1;

  std::vector<vtkIdType> minIds, maxIds, ids;
  if (this->PartitionCount.size() > 0 &&
      this->PartitionByBoundingBoxes(numPieces, minIds, maxIds,
                                     this->PieceBounds,
                                     this->PieceBoundsHalo)) {
    particleStart = minIds[piece];
    particleEnd = maxIds[piece];
  } else if (piece < totalParticles) {
    if (this->RandomizePartitionExtents) {
      this->PartitionByExtentsRandomized(totalParticles, piece, numPieces, ids);
      particleStart = ids[0];
      particleEnd = ids[1];
    } else {
      this->PartitionByExtents(totalParticles, piece, numPieces, ids);
      particleStart = ids[0];
      particleEnd = ids[1];
    }
  } else {
    particleStart = 0;
    particleEnd = -1;
  }

  vtkIdType Nt = particleEnd - particleStart + 1;
  bool viewIsIndexed = false;

  // If MaxParticlesPerRank is set, clamp the read to
  // MaxParticlesPerRank particles of this rank's partition.  By
  // default this reads the first N contiguous particles; when
  // UseStridedMaxParticlesPerRank is on, the N particles are spread
  // evenly across the partition, sampling its full spatial extent.
  if (this->MaxParticlesPerRank > 0 && Nt > this->MaxParticlesPerRank)
  {
    if (this->UseStridedMaxParticlesPerRank)
    {
      std::vector<h5_size_t> indices;
      indices.reserve(this->MaxParticlesPerRank);
      vtkIdType stride = Nt / this->MaxParticlesPerRank;
      if (stride < 1)
      {
        stride = 1;
      }
      for (vtkIdType i = 0; i < this->MaxParticlesPerRank; ++i)
      {
        indices.push_back(
            static_cast<h5_size_t>(particleStart + i * stride));
      }
      H5PartSetViewIndices(this->H5FileId, indices.data(), indices.size());
      Nt = static_cast<vtkIdType>(indices.size());
      viewIsIndexed = true;
    }
    else
    {
      particleEnd = particleStart + this->MaxParticlesPerRank - 1;
      Nt = this->MaxParticlesPerRank;
    }
  }

  if (Nt > 0) {
    if (!viewIsIndexed)
    {
      H5PartSetView(this->H5FileId, particleStart, particleEnd);
    }
  } else {
    H5PartSetView(this->H5FileId, -1, -1);
    Nt = 0;
  }

  // Read coordinate/scalar data.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDataArray> coords = nullptr;
  for (const auto &field : scalarFields) {
    const std::vector<std::string> &arraylist = field.second;
    const char *array_name = arraylist[0].c_str();
    std::string rootname = this->NameOfVectorComponent(array_name);
    int Nc = static_cast<int>(arraylist.size());

    h5_int64_t datatype, datatype_comp;
    h5_size_t nelem, nelem_comp;
    if (H5PartGetDatasetInfoByName(this->H5FileId, array_name, &datatype,
                                   &nelem) != H5_SUCCESS) {
      vtkErrorMacro("Could not get dataset info for array " << array_name);
      return 0;
    }
    int vtk_datatype = GetVTKDataType(datatype);

    if (vtk_datatype != VTK_VOID) {
      vtkSmartPointer<vtkDataArray> dataArray;
      dataArray.TakeReference(vtkDataArray::CreateDataArray(vtk_datatype));
      dataArray->SetNumberOfComponents(Nc);
      dataArray->SetNumberOfTuples(Nt);
      dataArray->SetName(rootname.c_str());

      for (int c = 0; c < Nc; ++c) {
        const char *name_comp = arraylist[c].c_str();
        if (H5PartGetDatasetInfoByName(this->H5FileId, name_comp,
                                       &datatype_comp,
                                       &nelem_comp) != H5_SUCCESS) {
          vtkErrorMacro("Could not get dataset info for array " << name_comp);
          return 0;
        }
        if (datatype_comp != datatype) {
          vtkErrorMacro("Inconsistent data types for vector components of "
                        << rootname);
          return 0;
        }

        vtkSmartPointer<vtkDataArray> component;
        component.TakeReference(vtkDataArray::CreateDataArray(vtk_datatype));
        component->SetNumberOfComponents(1);
        component->SetNumberOfTuples(Nt);
        component->SetName(name_comp);
        if (H5hutReadDataArray(this->H5FileId, name_comp, component) !=
            H5_SUCCESS) {
          vtkErrorMacro("Failed to read component " << name_comp);
          return 0;
        }

        if (Nc == 1) {
          dataArray->DeepCopy(component);
        } else {
          switch (vtk_datatype) {
          case VTK_FLOAT:
            this->CopyIntoVector<float>(c, component, dataArray);
            break;
          case VTK_DOUBLE:
            this->CopyIntoVector<double>(c, component, dataArray);
            break;
          case VTK_CHAR:
          case VTK_SIGNED_CHAR:
          case VTK_UNSIGNED_CHAR:
            this->CopyIntoVector<char>(c, component, dataArray);
            break;
          case VTK_SHORT:
            this->CopyIntoVector<short>(c, component, dataArray);
            break;
          case VTK_UNSIGNED_SHORT:
            this->CopyIntoVector<unsigned short>(c, component, dataArray);
            break;
          case VTK_INT:
            this->CopyIntoVector<int>(c, component, dataArray);
            break;
          case VTK_UNSIGNED_INT:
            this->CopyIntoVector<unsigned int>(c, component, dataArray);
            break;
          case VTK_LONG:
            this->CopyIntoVector<long>(c, component, dataArray);
            break;
          case VTK_UNSIGNED_LONG:
            this->CopyIntoVector<unsigned long>(c, component, dataArray);
            break;
          case VTK_LONG_LONG:
            this->CopyIntoVector<long long>(c, component, dataArray);
            break;
          case VTK_UNSIGNED_LONG_LONG:
            this->CopyIntoVector<unsigned long long>(c, component, dataArray);
            break;
          case VTK_ID_TYPE:
            this->CopyIntoVector<vtkIdType>(c, component, dataArray);
            break;
          default:
            vtkErrorMacro("H5Part : Unhandled vector type");
          }
        }
      }

      if (field.first == "Coords") {
        coords = dataArray;
      } else {
        output->GetPointData()->AddArray(dataArray);
        if (!output->GetPointData()->GetScalars()) {
          output->GetPointData()->SetActiveScalars(dataArray->GetName());
        }
      }
    }
  }

  if (this->GenerateVertexCells && Nt > 0) {
    vtkNew<vtkTypeInt64Array> connectivity;
    connectivity->SetNumberOfValues(Nt);
    std::iota(connectivity->Begin(), connectivity->End(), 0);
    vtkNew<vtkCellArray> vertices;
    vertices->SetData(1, connectivity);
    output->SetVerts(vertices);
  }

  if (coords) {
    coords->SetName("Points");
    points->SetData(coords);
    output->SetPoints(points);
  }

  if (!this->IgnorePartitionBoxes &&
      (this->DisplayPartitionBoxes || this->DisplayPieceBoxes) && Nt > 0) {
    this->DisplayBoundingBoxes(coords, output, particleStart, particleEnd);
  }

  this->CloseFileIntermediate();
  return 1;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::GetCoordinateArrayStatus(const char *name) {
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::SetCoordinateArrayStatus(const char *name, int status) {
  if (status) {
    this->PointDataArraySelection->EnableArray(name);
  } else {
    this->PointDataArraySelection->DisableArray(name);
  }
}

//----------------------------------------------------------------------------
const char *vtkH5PartReaderV2::GetPointArrayName(int index) {
  return this->PointDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::GetPointArrayStatus(const char *name) {
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::SetPointArrayStatus(const char *name, int status) {
  if (status != this->GetPointArrayStatus(name)) {
    if (status) {
      this->PointDataArraySelection->EnableArray(name);
    } else {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::Enable(const char *name) {
  this->SetPointArrayStatus(name, 1);
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::Disable(const char *name) {
  this->SetPointArrayStatus(name, 0);
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::EnableAll() {
  this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::DisableAll() {
  this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::GetNumberOfPointArrays() {
  return this->PointDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
bool vtkH5PartReaderV2::HasStep(int step) {
  if (!this->OpenFile()) {
    return false;
  }
  return (step >= 0 && step < this->NumberOfTimeSteps);
}

//----------------------------------------------------------------------------
void vtkH5PartReaderV2::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << (this->FileName ? this->FileName : "(none)")
     << "\n";
  os << indent << "NumberOfSteps: " << this->NumberOfTimeSteps << "\n";
}

//----------------------------------------------------------------------------
vtkIdType vtkH5PartReaderV2::ReadBoundingBoxes() {
  vtkIdType partitions = 0;

  hid_t partitiongroup = -1;
  hid_t dataset_id = -1;
  H5E_BEGIN_TRY {
    hid_t hdf5FileId = vtkH5hutGetHDF5FileId(this->H5FileId);
    partitiongroup = H5Gopen(hdf5FileId, "Partition#0", H5P_DEFAULT);
    if (partitiongroup > 0) {
      dataset_id = H5Dopen(partitiongroup, "Box", H5P_DEFAULT);
    }
  }
  H5E_END_TRY;

  if (partitiongroup > 0 && dataset_id > 0) {
    hid_t space = H5Dget_space(dataset_id);
    hsize_t dims[2], maxdims[2];
    herr_t err = H5Sget_simple_extent_dims(space, dims, maxdims);
    if (err != 1) {
      vtkErrorMacro("Error in H5Part bounding box dimensions read");
    }
    partitions = dims[0] / 13;

    std::vector<double> data(dims[0], 0.0);
    err = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &data[0]);
    if (err < 0) {
      vtkErrorMacro("Error in H5Part bounding box data read");
    }

    this->PartitionCount.resize(partitions);
    this->PartitionOffset.resize(partitions + 1, 0);
    this->PieceId.resize(partitions, 0);
    this->PartitionBoundsTable.resize(partitions * 6);
    this->PartitionBoundsTableHalo.resize(partitions * 6);

    for (vtkIdType i = 0; i < partitions; ++i) {
      vtkIdType offset1 = i * 6;
      vtkIdType offset2 = i * 13;
      this->PartitionCount[i] = static_cast<vtkIdType>(data[0 + offset2]);
      this->PartitionBoundsTable[0 + offset1] = data[1 + offset2];
      this->PartitionBoundsTable[1 + offset1] = data[4 + offset2];
      this->PartitionBoundsTable[2 + offset1] = data[2 + offset2];
      this->PartitionBoundsTable[3 + offset1] = data[5 + offset2];
      this->PartitionBoundsTable[4 + offset1] = data[3 + offset2];
      this->PartitionBoundsTable[5 + offset1] = data[6 + offset2];

      this->PartitionBoundsTableHalo[0 + offset1] = data[1 + 6 + offset2];
      this->PartitionBoundsTableHalo[1 + offset1] = data[4 + 6 + offset2];
      this->PartitionBoundsTableHalo[2 + offset1] = data[2 + 6 + offset2];
      this->PartitionBoundsTableHalo[3 + offset1] = data[5 + 6 + offset2];
      this->PartitionBoundsTableHalo[4 + offset1] = data[3 + 6 + offset2];
      this->PartitionBoundsTableHalo[5 + offset1] = data[6 + 6 + offset2];
    }
    std::partial_sum(this->PartitionCount.begin(), this->PartitionCount.end(),
                     this->PartitionOffset.begin() + 1);

    H5Sclose(space);
  }
  if (dataset_id > 0) {
    H5Dclose(dataset_id);
  }
  if (partitiongroup > 0) {
    H5Gclose(partitiongroup);
  }

  return partitions;
}

//----------------------------------------------------------------------------
vtkIdType vtkH5PartReaderV2::DisplayBoundingBoxes(vtkDataArray *coords,
                                                  vtkPolyData *output,
                                                  vtkIdType extent0,
                                                  vtkIdType extent1) {
  vtkIdType partitions =
      this->DisplayPartitionBoxes ? this->PartitionCount.size() : 0;
  vtkIdType pieces = this->DisplayPieceBoxes ? this->PieceBounds.size() : 0;
  vtkIdType newBoxes = partitions + pieces;

  if (newBoxes == 0) {
    return 0;
  }

  vtkIdType N1 = coords->GetNumberOfTuples();
  vtkIdType N2 = newBoxes * 8 * 2;
  vtkIdType N3 = N1 + N2;

  vtkSmartPointer<vtkIdTypeArray> occupation =
      vtkSmartPointer<vtkIdTypeArray>::New();
  occupation->SetNumberOfTuples(N3);
  occupation->SetName("Occupation");
  vtkSmartPointer<vtkIdTypeArray> boxId =
      vtkSmartPointer<vtkIdTypeArray>::New();
  boxId->SetNumberOfTuples(N3);
  boxId->SetName("Partition");
  vtkSmartPointer<vtkIdTypeArray> piece =
      vtkSmartPointer<vtkIdTypeArray>::New();
  piece->SetNumberOfTuples(N3);
  piece->SetName("PieceId");

  vtkIdType index = 0;
  vtkBoundingBox box;
  vtkSmartPointer<vtkAppendPolyData> polys =
      vtkSmartPointer<vtkAppendPolyData>::New();

  for (size_t i = 0; i < this->PartitionCount.size(); ++i) {
    vtkIdType numParticles = this->PartitionCount[i];
    vtkIdType pieceId = this->PieceId[i];
    if ((index + numParticles) >= extent0 && (index <= extent1)) {
      for (vtkIdType p = 0; p < numParticles; ++p) {
        if (index >= extent0 && index <= extent1) {
          boxId->SetValue(index - extent0, static_cast<vtkIdType>(i));
          occupation->SetValue(index - extent0, numParticles);
          piece->SetValue(index - extent0, pieceId);
        }
        ++index;
      }
    } else {
      index += numParticles;
    }
  }

  index = (extent1 - extent0) + 1;

  for (vtkIdType i = 0; i < partitions; ++i) {
    vtkSmartPointer<vtkOutlineSource> cube1 =
        vtkSmartPointer<vtkOutlineSource>::New();
    cube1->SetBounds(&this->PartitionBoundsTable[i * 6]);
    cube1->Update();
    polys->AddInputData(cube1->GetOutput());

    vtkSmartPointer<vtkOutlineSource> cube2 =
        vtkSmartPointer<vtkOutlineSource>::New();
    cube2->SetBounds(&this->PartitionBoundsTableHalo[i * 6]);
    cube2->Update();
    polys->AddInputData(cube2->GetOutput());
  }

  for (size_t i = 0; i < pieces; ++i) {
    vtkSmartPointer<vtkOutlineSource> cube1 =
        vtkSmartPointer<vtkOutlineSource>::New();
    double bbb[6];
    this->PieceBounds[i].GetBounds(bbb);
    cube1->SetBounds(bbb);
    cube1->Update();
    polys->AddInputData(cube1->GetOutput());

    this->PieceBoundsHalo[i].GetBounds(bbb);
    vtkSmartPointer<vtkOutlineSource> cube2 =
        vtkSmartPointer<vtkOutlineSource>::New();
    cube2->SetBounds(bbb);
    cube2->Update();
    polys->AddInputData(cube2->GetOutput());
  }
  polys->Update();

  vtkPoints *points = polys->GetOutput()->GetPoints();
  coords->Resize(N3);
  coords->SetNumberOfTuples(N3);
  for (vtkIdType P = 0; P < N2; ++P) {
    coords->SetTuple(N1 + P, points->GetPoint(P));
  }

  output->SetLines(polys->GetOutput()->GetLines());
  vtkIdType L = output->GetLines()->GetNumberOfCells();
  vtkIdTypeArray *linedata = output->GetLines()->GetData();
  for (vtkIdType B = 0; B < L; ++B) {
    linedata->SetValue(B * 3 + 1, N1 + linedata->GetValue(B * 3 + 1));
    linedata->SetValue(B * 3 + 2, N1 + linedata->GetValue(B * 3 + 2));
  }

  for (vtkIdType i = 0; i < partitions; ++i) {
    vtkIdType numParticles = this->PartitionCount[i];
    vtkIdType pieceId = this->PieceId[i];
    for (vtkIdType p = 0; p < 8 * 2; ++p) {
      boxId->SetValue(index, i);
      occupation->SetValue(index, numParticles);
      piece->SetValue(index, pieceId);
      ++index;
    }
  }
  for (size_t i = 0; i < pieces; ++i) {
    vtkIdType pieceId = static_cast<vtkIdType>(i);
    for (vtkIdType p = 0; p < 8 * 2; ++p) {
      boxId->SetValue(index, 0);
      occupation->SetValue(index, 0);
      piece->SetValue(index, pieceId);
      ++index;
    }
  }

  output->GetPointData()->AddArray(occupation);
  output->GetPointData()->AddArray(boxId);
  output->GetPointData()->AddArray(piece);

  return N2;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::PartitionByExtents(vtkIdType n, int piece, int numPieces,
                                          std::vector<vtkIdType> &startend) {
  vtkIdType WholeExtent[6] = {0, n, 0, 0, 0, 0};
  this->SplitExtent(piece, numPieces, WholeExtent);
  startend.push_back(WholeExtent[0]);
  startend.push_back(WholeExtent[1] - 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::PartitionByExtentsRandomized(
    vtkIdType n, int piece, int numPieces, std::vector<vtkIdType> &startend) {
  Random r(12345);
  vtkIdType partitionsize = n / numPieces;
  vtkIdType rand_max = partitionsize / 2;
  vtkIdType rand_half = partitionsize / 4;

  vtkIdType pstart = 0;
  vtkIdType pend = 0;
  vtkIdType rnddev;
  for (int i = 0; i < numPieces; ++i) {
    vtkIdType epstart = i * partitionsize;
    vtkIdType epend = (i + 1) * partitionsize;
    rnddev = (r.nextNumberInt() % rand_max) - rand_half;
    pend = epend + rnddev;
    if (pend >= n) {
      pend = n - 1;
    }
    if (i == (numPieces - 1)) {
      pend = n - 1;
    }
    if (i == piece) {
      startend.push_back(pstart);
      startend.push_back(pend);
    }
    pstart = pend + 1;
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::PartitionByBoundingBoxes(
    int numPieces, std::vector<vtkIdType> &minIds,
    std::vector<vtkIdType> &maxIds, std::vector<vtkBoundingBox> &pieceBounds,
    std::vector<vtkBoundingBox> &pieceHaloBounds) {
  if (this->IgnorePartitionBoxes) {
    return 0;
  }
  vtkIdType np = this->PartitionCount.size();
  vtkIdType nr = static_cast<vtkIdType>(numPieces);

  if (nr == 1 || np == 0) {
    return 0;
  }

  pieceBounds.assign(nr, vtkBoundingBox());
  pieceHaloBounds.assign(nr, vtkBoundingBox());
  minIds.assign(nr, VTK_INT_MAX);
  maxIds.assign(nr, VTK_INT_MIN);

  if (this->UseLinearBoxPartitioning) {
    for (vtkIdType pieceId = 0; pieceId < nr; ++pieceId) {
      vtkIdType startPartition = (pieceId * np) / nr;
      vtkIdType endPartition = ((pieceId + 1) * np) / nr - 1;
      if (startPartition > endPartition) {
        minIds[pieceId] = 0;
        maxIds[pieceId] = -1;
        continue;
      }
      for (vtkIdType i = startPartition; i <= endPartition; ++i) {
        this->PieceId[i] = pieceId;
        minIds[pieceId] = std::min(this->PartitionOffset[i], minIds[pieceId]);
        maxIds[pieceId] =
            std::max(this->PartitionOffset[i + 1] - 1, maxIds[pieceId]);
        pieceBounds[pieceId].AddBounds(&this->PartitionBoundsTable[i * 6]);
        pieceHaloBounds[pieceId].AddBounds(
            &this->PartitionBoundsTableHalo[i * 6]);
      }
    }
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkH5PartReaderV2::SplitExtent(int piece, int numPieces, vtkIdType *ext) {
  int numPiecesInFirstHalf;
  unsigned long size[3];
  int splitAxis;
  vtkIdType mid;

  if (piece >= numPieces || piece < 0) {
    return 0;
  }

  while (numPieces > 1) {
    size[0] = ext[1] - ext[0];
    size[1] = ext[3] - ext[2];
    size[2] = ext[5] - ext[4];
    if (size[2] >= size[1] && size[2] >= size[0] && size[2] / 2 >= 1) {
      splitAxis = 2;
    } else if (size[1] >= size[0] && size[1] / 2 >= 1) {
      splitAxis = 1;
    } else if (size[0] / 2 >= 1) {
      splitAxis = 0;
    } else {
      splitAxis = -1;
    }

    if (splitAxis == -1) {
      return 1;
    }

    mid =
        (ext[splitAxis * 2 + 1] - ext[splitAxis * 2]) / 2 + ext[splitAxis * 2];

    numPiecesInFirstHalf = (numPieces / 2);
    if (piece < numPiecesInFirstHalf) {
      ext[splitAxis * 2 + 1] = mid;
    } else {
      ext[splitAxis * 2] = mid + 1;
      piece -= numPiecesInFirstHalf;
      numPieces -= numPiecesInFirstHalf;
      continue;
    }
    numPieces = numPiecesInFirstHalf;
  }

  return 1;
}
