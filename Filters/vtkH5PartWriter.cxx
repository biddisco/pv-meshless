/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5PartWriter.cxx
  Copyright (C) CSCS - Swiss National Supercomputing Centre.

=========================================================================*/
#include "vtkH5PartWriter.h"
#include "vtkH5hutHelper.h"

#include "vtkDataArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
//
#include "vtkCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#endif
#include "vtkDummyController.h"
//
#include <algorithm>
#include <functional>
#include <sstream>
#include <stdlib.h>
//
// vtksys
//
#include <vtksys/SystemTools.hxx>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkH5PartWriter);
vtkCxxSetObjectMacro(vtkH5PartWriter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
namespace {
// Replace any character from @c bad with @c replacement.
std::string ReplaceBadChars(const std::string &input, const char *bad,
                            char replacement) {
  std::string output = input;
  std::string::size_type pos = 0;
  while ((pos = output.find_first_of(bad, pos)) != std::string::npos) {
    output[pos] = replacement;
    ++pos;
  }
  return output;
}
} // namespace
//----------------------------------------------------------------------------
#ifdef JB_DEBUG__
#ifdef WIN32
#define OUTPUTTEXT(a) vtkOutputWindowDisplayText(a);
#else
#define OUTPUTTEXT(a)                                                          \
  std::cout << (a) << "\n";                                                    \
  std::cout.flush();
#endif

#undef vtkDebugMacro
#define vtkDebugMacro(a)                                                       \
  {                                                                            \
    vtkOStreamWrapper::EndlType endl;                                          \
    vtkOStreamWrapper::UseEndl(endl);                                          \
    vtkOStrStreamWrapper vtkmsg;                                               \
    vtkmsg a << "\n";                                                          \
    OUTPUTTEXT(vtkmsg.str());                                                  \
    vtkmsg.rdbuf()->freeze(0);                                                 \
  }

#undef vtkErrorMacro
#define vtkErrorMacro(a) vtkDebugMacro(a)
#endif
//----------------------------------------------------------------------------
vtkH5PartWriter::vtkH5PartWriter() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  //
  this->NumberOfTimeSteps = 0;
  this->TimeStep = 0;
  this->TimeValue = 0.0;
  this->NumberOfParticles = 0;
  this->FileName = NULL;
  this->H5FileId = 0;
  this->FileMode = H5_O_RDWR;
  this->UpdatePiece = -1;
  this->UpdateNumPieces = -1;
  this->VectorsWithStridedWrite = 0;
  this->DisableInformationGather = 0;
  this->StepName = NULL;
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkH5PartWriter::~vtkH5PartWriter() {
  this->CloseFile();
  this->SetController(NULL);
  //
  delete[] this->StepName;
  this->StepName = NULL;
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::FillInputPortInformation(int, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::FillOutputPortInformation(int vtkNotUsed(port),
                                               vtkInformation *info) {
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::RequestInformation(vtkInformation *vtkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
    int NumberOfInputTimeSteps =
        inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    //
    // Get list of input time step values
    this->InputTimeValues.resize(NumberOfInputTimeSteps);
    inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                &this->InputTimeValues[0]);
  }
  return 1;
}
//----------------------------------------------------------------------------
vtkPointSet *vtkH5PartWriter::GetInput() {
  return vtkPointSet::SafeDownCast(this->GetInput(0));
}
//----------------------------------------------------------------------------
vtkPointSet *vtkH5PartWriter::GetInput(int port) {
  return vtkPointSet::SafeDownCast(vtkAbstractParticleWriter::GetInput(port));
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::SetFileModeToWrite() { this->SetFileMode(H5_O_WRONLY); }
//----------------------------------------------------------------------------
void vtkH5PartWriter::SetFileModeToReadWrite() { this->SetFileMode(H5_O_RDWR); }
//----------------------------------------------------------------------------
void vtkH5PartWriter::CloseFile() { H5hutCloseFile(this->H5FileId); }
//----------------------------------------------------------------------------
int vtkH5PartWriter::OpenFile() {
  if (!this->FileName) {
    vtkErrorMacro("FileName must be specified.");
    return 0;
  }
  // if file doesn't exist already, use write mode to force creation
  int actualMode = this->FileMode;
  if (this->FileMode == H5_O_RDWR &&
      !vtksys::SystemTools::FileExists(this->FileName)) {
    actualMode = H5_O_WRONLY;
  }

  if (!this->H5FileId) {
    this->H5FileId =
        H5hutOpenFile(this->FileName, actualMode, this->Controller);
  }
  if (!this->H5FileId) {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
  }
  return 1;
}
//----------------------------------------------------------------------------
struct vtkH5PW_datainfo {
  int datatype;
  int numC;
  char name[64];
  vtkH5PW_datainfo() : datatype(-1), numC(-1) {};
};
//----------------------------------------------------------------------------
bool vtkH5PartWriter::GatherDataArrayInfo(vtkDataArray *data, int &datatype,
                                          std::string &dataname,
                                          int &numComponents) {
#ifdef PARAVIEW_USE_MPI
  std::vector<vtkH5PW_datainfo> datatypes(this->UpdateNumPieces);
  if (data) {
    ((vtkH5PW_datainfo *)&datatypes[this->UpdatePiece])->datatype =
        data->GetDataType();
    ((vtkH5PW_datainfo *)&datatypes[this->UpdatePiece])->numC =
        data->GetNumberOfComponents();
    strncpy(((vtkH5PW_datainfo *)&datatypes[this->UpdatePiece])->name,
            data->GetName(), 64);
  }
  vtkMPICommunicator *com =
      vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int result = com->AllGather((char *)MPI_IN_PLACE, (char *)&datatypes[0],
                              sizeof(vtkH5PW_datainfo));
  for (int i = 0; i < this->UpdateNumPieces; i++) {
    vtkH5PW_datainfo &newdata = datatypes[i];
    if (newdata.datatype != -1) {
      datatype = newdata.datatype;
      numComponents = newdata.numC;
      dataname = newdata.name;
    }
  }
  return (result == 1);
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
bool vtkH5PartWriter::GatherScalarInfo(vtkPointData *pd, int N,
                                       int &numScalar) {
#ifdef PARAVIEW_USE_MPI
  numScalar = N;
  std::vector<int> numScalars(this->UpdateNumPieces, 0);
  if (pd)
    numScalars[this->UpdatePiece] = pd->GetNumberOfArrays();
  vtkMPICommunicator *com =
      vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int result = com->AllGather((const int *)MPI_IN_PLACE, &numScalars[0], 1);
  for (int i = 0; i < this->UpdateNumPieces; i++) {
    if (numScalars[i] > 0)
      numScalar = numScalars[i];
  }
  return (result == 1);
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
template <class T1, class T2>
void CopyFromVector_T(int offset, vtkDataArray *source, vtkDataArray *dest) {
  int N = source->GetNumberOfTuples();
  T1 *sptr = static_cast<T1 *>(source->GetVoidPointer(0)) + offset;
  T2 *dptr = static_cast<T2 *>(dest->WriteVoidPointer(0, N));
  for (int i = 0; i < N; ++i) {
    *dptr++ = *sptr;
    sptr += 3;
  }
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::CopyFromVector(int offset, vtkDataArray *source,
                                     vtkDataArray *dest) {
  switch (source->GetDataType()) {
  case VTK_CHAR:
  case VTK_SIGNED_CHAR:
  case VTK_UNSIGNED_CHAR:
    CopyFromVector_T<char, char>(offset, source, dest);
    break;
  case VTK_SHORT:
    CopyFromVector_T<short int, short int>(offset, source, dest);
    break;
  case VTK_UNSIGNED_SHORT:
    CopyFromVector_T<unsigned short int, unsigned short int>(offset, source,
                                                             dest);
    break;
  case VTK_INT:
    CopyFromVector_T<int, int>(offset, source, dest);
    break;
  case VTK_UNSIGNED_INT:
    CopyFromVector_T<unsigned int, unsigned int>(offset, source, dest);
    break;
  case VTK_LONG:
    CopyFromVector_T<long int, long int>(offset, source, dest);
    break;
  case VTK_UNSIGNED_LONG:
    CopyFromVector_T<unsigned long int, unsigned long int>(offset, source,
                                                           dest);
    break;
  case VTK_LONG_LONG:
    CopyFromVector_T<long long, long long>(offset, source, dest);
    break;
  case VTK_UNSIGNED_LONG_LONG:
    CopyFromVector_T<unsigned long long, unsigned long long>(offset, source,
                                                             dest);
    break;
  case VTK_FLOAT:
    CopyFromVector_T<float, float>(offset, source, dest);
    break;
  case VTK_DOUBLE:
    CopyFromVector_T<double, double>(offset, source, dest);
    break;
  case VTK_ID_TYPE:
    CopyFromVector_T<vtkIdType, vtkIdType>(offset, source, dest);
    break;
  default:
    vtkErrorMacro("Unexpected data type");
  }
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::WriteDataArray(int i, vtkDataArray *indata) {
  vtkSmartPointer<vtkDataArray> data;
  if (this->UpdateNumPieces > 1 && !this->DisableInformationGather) {
    int correctType = -1, numComponents = -1;
    std::string correctName;
    GatherDataArrayInfo(indata, correctType, correctName, numComponents);
    if (!indata) {
      vtkDebugMacro("NULL data found, used MPI_Gather to find :"
                    << " DataType " << correctType << " Name "
                    << correctName.c_str() << " NumComponents "
                    << numComponents);
      data = vtkDataArray::CreateDataArray(correctType);
      data->FastDelete();
      data->SetNumberOfComponents(numComponents);
      data->SetName(correctName.c_str());
    } else {
      data = indata;
    }
  } else
    data = indata;
  //
  int Nt = data->GetNumberOfTuples();
  int Nc = data->GetNumberOfComponents();
  //
  vtkSmartPointer<vtkDataArray> component;
  if (Nc > 1 && !this->VectorsWithStridedWrite) {
    component.TakeReference(data->NewInstance());
    component->SetNumberOfComponents(1);
    component->SetNumberOfTuples(Nt);
    component->WriteVoidPointer(0, Nt);
  }

  char buffer[8];
  char BadChars[] = "/\\:*?\"<> ";
  for (int c = 0; c < Nc; c++) {
    std::string name;
    if (data->GetName()) {
      name = data->GetName();
    } else {
      std::ostringstream oss;
      oss << "Scalars_" << i;
      name = oss.str();
    }
    name = ReplaceBadChars(name, BadChars, '_');

    vtkDataArray *finalData = data;
    if (Nc > 1) {
      std::ostringstream oss;
      oss << name << "_" << c;
      name = oss.str();
      if (!this->VectorsWithStridedWrite) {
        this->CopyFromVector(c, data, component);
        finalData = component;
      }
      // strided write not supported in this simplified port
    }

    if (H5hutWriteDataArray(this->H5FileId, name.c_str(), finalData) !=
        H5_SUCCESS) {
      vtkErrorMacro("Array write failed for name " << name.c_str());
    } else {
      vtkDebugMacro("Wrote " << name.c_str() << " " << Nt << " " << Nc);
    }
  }
}
//----------------------------------------------------------------------------
class H5PartWriterToleranceCheck
    : public std::binary_function<double, double, bool> {
public:
  H5PartWriterToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
  result_type operator()(first_argument_type a, second_argument_type b) const {
    bool result = (fabs(a - b) <= (this->tolerance));
    return (result_type)result;
  }
};
//----------------------------------------------------------------------------
bool vtkH5PartWriter::WriteDataAndReturn() {
  if (this->Controller) {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  } else {
    this->UpdatePiece = 0;
    this->UpdateNumPieces = 1;
  }

  //
  // Make sure file is open
  //
  if (!this->OpenFile()) {
    vtkErrorMacro("Couldn't create file " << this->FileName);
    return false;
  }

  vtkInformation *outInfo = this->GetExecutive()->GetOutputInformation(0);

  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    this->TimeValue =
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    this->ActualTimeStep =
        std::find_if(
            this->InputTimeValues.begin(), this->InputTimeValues.end(),
            std::bind2nd(H5PartWriterToleranceCheck(1E-6), this->TimeValue)) -
        this->InputTimeValues.begin();
    //
  }

  //
  // Set Step. This will create data group for us
  //
  H5SetStep(this->H5FileId, this->ActualTimeStep);
  //
  // Write out a TimeValue attribute for this step
  //
  if (h5_write_iteration_attrib(this->H5FileId, "TimeValue", H5_FLOAT64_T,
                                &this->TimeValue, 1) != H5_SUCCESS) {
    vtkErrorMacro("TimeValue attrib write failed ");
  }
  //
  // Get the input to write and Set Num-Particles
  //
  vtkPointSet *input = this->GetInput();
  this->NumberOfParticles = input->GetNumberOfPoints();
  H5PartSetNumParticles(this->H5FileId, this->NumberOfParticles);
  //
  // Write coordinate data
  //
  vtkSmartPointer<vtkPoints> points = input->GetPoints();
  if (points && points->GetData()) {
    points->GetData()->SetName("Coords");
    this->WriteDataArray(0, points->GetData());
  } else {
    this->WriteDataArray(0, NULL);
  }
  //
  // Write point data
  //
  int numScalars, numFound;
  vtkPointData *pd = input->GetPointData();
  numScalars = pd->GetNumberOfArrays();
  if (this->UpdateNumPieces > 1 && !this->DisableInformationGather) {
    GatherScalarInfo(pd, numScalars, numFound);
    if (numScalars != numFound) {
      vtkDebugMacro("No scalars found, used MPI_Gather to find " << numFound
                                                                 << " Arrays");
    }
  } else
    numFound = numScalars;
  for (int i = 0; i < numFound; i++) {
    vtkDataArray *data = pd->GetArray(i);
    this->WriteDataArray(i, data);
  }
  //
  // We are done.
  //
  vtkDebugMacro("Time Step written " << this->FileName << " Step "
                                     << this->ActualTimeStep << " Time "
                                     << this->TimeValue);
  return true;
}
//----------------------------------------------------------------------------
bool vtkH5PartWriter::IsTimeStepPresent(int timestep) {
  if (!this->OpenFile()) {
    vtkErrorMacro("Couldn't open file " << this->FileName);
    return false;
  }
  return (timestep >= 0 &&
          timestep < static_cast<int>(H5GetNumSteps(this->H5FileId)));
}
//----------------------------------------------------------------------------
// Delete Timestep not working yet
//----------------------------------------------------------------------------
void vtkH5PartWriter::DeleteTimeStep(int vtkNotUsed(timestep)) {
  return; // because this does not work
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "FileName: " << (this->FileName ? this->FileName : "(none)")
     << "\n";

  os << indent << "NumberOfSteps: " << this->NumberOfTimeSteps << "\n";
}
