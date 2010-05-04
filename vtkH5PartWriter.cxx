/*=========================================================================

  Project                 : vtkCSCS
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
#include "vtkH5PartWriter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPointSet.h"
#include "vtkDataArray.h"
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
//
#ifdef VTK_USE_MPI
#include "vtkMPI.h"
#include "vtkMultiProcessController.h"
#include "vtkMPICommunicator.h"
#endif
//
#include "H5Part.h"
#include <stdlib.h>
#include <algorithm>
//
// vtksys
//
#include <vtksys/SystemTools.hxx>
//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkH5PartWriter, "$Revision: 153 $");
vtkStandardNewMacro(vtkH5PartWriter);
#ifdef VTK_USE_MPI
vtkCxxSetObjectMacro(vtkH5PartWriter, Controller, vtkMultiProcessController);
#endif
//----------------------------------------------------------------------------
#ifdef JB_DEBUG__
  #ifdef WIN32
      #define OUTPUTTEXT(a) vtkOutputWindowDisplayText(a);
  #else
      #define OUTPUTTEXT(a) std::cout << (a) << "\n"; std::cout.flush();
  #endif

    #undef vtkDebugMacro
    #define vtkDebugMacro(a)  \
    { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    }

  #undef vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
vtkH5PartWriter::vtkH5PartWriter() 
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  //
  this->NumberOfTimeSteps         = 0;
  this->TimeStep                  = 0;
  this->TimeValue                 = 0.0;
  this->NumberOfParticles         = 0;
  this->FileName                  = NULL;
  this->H5FileId                  = NULL;
  this->FileMode                  = H5PART_APPEND;
  this->UpdatePiece               = -1;
  this->UpdateNumPieces           = -1;
  this->VectorsWithStridedWrite   = 0;
  this->DisableInformationGather  = 0;
#ifdef VTK_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkH5PartWriter::~vtkH5PartWriter()
{ 
  this->CloseFile();
#ifdef VTK_USE_MPI
  this->SetController(NULL);
#endif
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // we might be multiblock, we might be dataset
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);

  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) )
    {
    int NumberOfInputTimeSteps = inInfo->Length( 
      vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
    //
    // Get list of input time step values
    this->InputTimeValues.resize(NumberOfInputTimeSteps);
    inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
      &this->InputTimeValues[0] );
  }
  return 1;
}
//----------------------------------------------------------------------------
vtkPointSet* vtkH5PartWriter::GetInput()
{
  return vtkPointSet::SafeDownCast(this->GetInput(0));
}
//----------------------------------------------------------------------------
vtkPointSet* vtkH5PartWriter::GetInput(int port)
{
  return vtkPointSet::SafeDownCast(vtkAbstractParticleWriter::GetInput(port));
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::SetFileModeToWrite() 
{ 
  this->SetFileMode(H5PART_WRITE); 
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::SetFileModeToReadWrite() 
{ 
  this->SetFileMode(H5PART_APPEND); 
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::CloseFile()
{
  if (this->H5FileId != NULL) {
    if (H5PartCloseFile(this->H5FileId)!=H5PART_SUCCESS) {
      vtkErrorMacro(<<"CloseFile failed");
    }
    this->H5FileId = NULL;
  }
}
//----------------------------------------------------------------------------
int vtkH5PartWriter::OpenFile()
{
  if (!this->FileName) {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
  }
  // if file doesn't exists already, make sure mode is write to force a create
  int actualMode = this->FileMode;
  if (this->FileMode==H5PART_APPEND && !vtksys::SystemTools::FileExists(this->FileName))
  {
    actualMode = H5PART_WRITE;
  }

  // H5FD_MPIO_INDEPENDENT
  // H5FD_MPIO_COLLECTIVE

  if (!this->H5FileId) {
#ifdef PARALLEL_IO
    if (this->Controller && this->UpdateNumPieces>1) {
      vtkMPICommunicator *vtkComm = 
        vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
      MPI_Comm *handle = vtkComm->GetMPIComm()->GetHandle();

//      if (this->CollectiveIO) {
        this->H5FileId = H5PartOpenFileParallel(
          this->FileName, actualMode, *handle);
//      }
//      else {
//        this->H5FileId = H5PartOpenFileParallelIndependent(
//          this->FileName, actualMode, *handle);
//      }
    }
    else {
      this->H5FileId = H5PartOpenFile(this->FileName, actualMode);
      if (this->H5FileId) {
        this->H5FileId->comm = 0;
      }
    }
#else
      this->H5FileId = H5PartOpenFile(this->FileName, this->FileMode);
#endif  
  }
  if (!this->H5FileId) {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
  }
  return 1;
}
//----------------------------------------------------------------------------
struct vtkH5PW_datainfo {
  int  datatype;
  int  numC;
  char name[64];
  vtkH5PW_datainfo() : datatype(-1), numC(-1) {};
};
//----------------------------------------------------------------------------
bool vtkH5PartWriter::GatherDataArrayInfo(vtkDataArray *data, 
  int &datatype, std::string &dataname, int &numComponents)
{
#ifdef VTK_USE_MPI
  std::vector< vtkH5PW_datainfo > datatypes(this->UpdateNumPieces);
  if (data) {
    ((vtkH5PW_datainfo*)&datatypes[this->UpdatePiece])->datatype = data->GetDataType();
    ((vtkH5PW_datainfo*)&datatypes[this->UpdatePiece])->numC     = data->GetNumberOfComponents();
    strncpy(((vtkH5PW_datainfo*)&datatypes[this->UpdatePiece])->name, data->GetName(), 64);
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(
    this->Controller->GetCommunicator()); 
  int result = com->AllGather((char*)&datatypes[this->UpdatePiece], (char*)&datatypes[0], sizeof(vtkH5PW_datainfo));
  for (int i=0; i<this->UpdateNumPieces; i++) {
    vtkH5PW_datainfo &newdata = datatypes[i];
    if (newdata.datatype!=-1) {
      datatype = newdata.datatype;
      numComponents = newdata.numC;
      dataname = newdata.name;
    }
  }
  return (result == 1) ;
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
bool vtkH5PartWriter::GatherScalarInfo(vtkPointData *pd, int N, int &numScalar)
{
#ifdef VTK_USE_MPI
  numScalar = N;
  std::vector<int> numScalars(this->UpdateNumPieces, 0);
  if (pd) numScalars[this->UpdatePiece] = pd->GetNumberOfArrays();
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(
    this->Controller->GetCommunicator()); 
  int result = com->AllGather(&numScalars[this->UpdatePiece], &numScalars[0], 1);
  for (int i=0; i<this->UpdateNumPieces; i++) {
    if (numScalars[i]>0) numScalar = numScalars[i];
  }
  return (result == 1) ;
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
template <class T1, class T2>
void CopyFromVector_T(int offset, vtkDataArray *source, vtkDataArray *dest)
{
  int N = source->GetNumberOfTuples();
  T1 *sptr = static_cast<T1*>(source->GetVoidPointer(0)) + offset;
  T2 *dptr = static_cast<T2*>(dest->WriteVoidPointer(0,N));
  for (int i=0; i<N; ++i) {
    *dptr++ = *sptr;
    sptr += 3;
  }
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::CopyFromVector(int offset, vtkDataArray *source, vtkDataArray *dest)
{
  switch (source->GetDataType()) 
  {
    case VTK_CHAR:
    case VTK_SIGNED_CHAR:
    case VTK_UNSIGNED_CHAR:
      CopyFromVector_T<char,char>(offset, source, dest);
      break;
    case VTK_SHORT:
      CopyFromVector_T<short int,short int>(offset, source, dest);
      break;
    case VTK_UNSIGNED_SHORT:
      CopyFromVector_T<unsigned short int,unsigned short int>(offset, source, dest);
      break;
    case VTK_INT:
      CopyFromVector_T<int,int>(offset, source, dest);
      break;
    case VTK_UNSIGNED_INT:
      CopyFromVector_T<unsigned int,unsigned int>(offset, source, dest);
      break;
    case VTK_LONG:
      CopyFromVector_T<long int,long int>(offset, source, dest);
      break;
    case VTK_UNSIGNED_LONG:
      CopyFromVector_T<unsigned long int,unsigned long int>(offset, source, dest);
      break;
    case VTK_LONG_LONG:
      CopyFromVector_T<long long,long long>(offset, source, dest);
      break;
    case VTK_UNSIGNED_LONG_LONG:
      CopyFromVector_T<unsigned long long,unsigned long long>(offset, source, dest);
      break;
    case VTK_FLOAT:
      CopyFromVector_T<float,float>(offset, source, dest);
      break;
    case VTK_DOUBLE:
      CopyFromVector_T<double,double>(offset, source, dest);
      break;
    case VTK_ID_TYPE:
      CopyFromVector_T<vtkIdType,vtkIdType>(offset, source, dest);
      break;
    default:
      break;
      vtkErrorMacro(<<"Unexpected data type");
  }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// A Convenience Macro which does what we want for any dataset type
// note that we use a dummy null parameter and the ## macro operator to delay expansion
// of the T2 parameter (which causes problems if we don't)
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  #define h_params H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
#else
  #define h_params H5P_DEFAULT
#endif

#define H5PartWriteDataArray(null, T2, f, name, dataarray) \
  sprintf(typestring, "%s", #T2); \
  dataset = H5Dcreate(f->timegroup, name.c_str(), null##T2, f->shape, h_params); \
  if (dataset<0) { \
    vtkErrorMacro(<<"Dataset create failed for " << name.c_str() \
    << " Timestep " << f->timestep \
    << " Shape " << f->shape \
    << " Data Type " << #T2); \
    r = -1; \
  } else { \
    void *dataptr = dataarray->GetVoidPointer(0); \
    dataptr = dataptr ? dataptr : &buffer[0]; \
    r = H5Dwrite(dataset, T2, memshape, diskshape, H5P_DEFAULT, dataptr); \
  } \
  H5Dclose(dataset);  
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void vtkH5PartWriter::WriteDataArray(int i, vtkDataArray *indata)
{
  // if a process has zero points/scalars, then this routine is entered with
  // a null pointer, we must find out what the other processes are sending 
  // and do an empty send with the same type etc.
  vtkSmartPointer<vtkDataArray> data;
  if (this->UpdateNumPieces>1 && !this->DisableInformationGather) {
    int correctType = -1, numComponents = -1;
    std::string correctName;
    GatherDataArrayInfo(indata, correctType, correctName, numComponents);
    if (!indata) {
      vtkDebugMacro(<<"NULL data found, used MPI_Gather to find :" 
        << " DataType " << correctType
        << " Name " << correctName.c_str()
        << " NumComponents " << numComponents);
      data = vtkDataArray::CreateDataArray(correctType);
      data->Delete(); // smartpointer copied it
      data->SetNumberOfComponents(numComponents);
      data->SetName(correctName.c_str());
    }
    else {
      data = indata;
    }
  }
  else data = indata;
  //

  hid_t dataset;
  hid_t &diskshape = H5FileId->diskshape;
  hid_t &memshape  = H5FileId->memshape;
  if (memshape!=H5S_ALL) {
    if (H5Sclose(memshape)<0) vtkErrorMacro(<<"memshape : HANDLE_H5S_CLOSE_ERR");
		memshape = H5S_ALL;
	}
  herr_t r=0;
  int Nt = data->GetNumberOfTuples();
  int Nc = data->GetNumberOfComponents();
  hsize_t     count1_mem[] = { Nt*Nc };
  hsize_t     count2_mem[] = { Nt };
  hsize_t     offset_mem[] = { 0 };
  hsize_t     stride_mem[] = { Nc };
  //
  vtkSmartPointer<vtkDataArray> component;
  if (!this->VectorsWithStridedWrite) {
    // we copy from original to a single component array
    count1_mem[0] = Nt;
    if (Nc>1) {
      component.TakeReference(data->NewInstance());
      component->SetNumberOfComponents(1);
      component->SetNumberOfTuples(Nt);
      component->WriteVoidPointer(0, Nt);
    }
  }
  char buffer[8];
  char BadChars[] = "/\\:*?\"<> ";
  char typestring[128];
  for (int c=0; c<Nc; c++) {
    // set the array name
    sprintf(buffer,"%i", i);
    vtkstd::string name = vtkstd::string("Scalars_").append(buffer);
    if (data->GetName()) name = data->GetName();
    char *tempname = const_cast<char *>(name.c_str());
    name = vtksys::SystemTools::ReplaceChars(tempname, BadChars, '_');
    // shape
    memshape = H5Screate_simple(1, count1_mem, NULL);   
    // single vector write or component by component
    vtkSmartPointer<vtkDataArray> finalData = data;
    if (Nc>1) {
      sprintf(buffer,"_%i", c);
      name = name.append(buffer);
      if (!this->VectorsWithStridedWrite) {
        this->CopyFromVector(c, data, component);
        finalData = component;
      }
      else {
        offset_mem[0] = c;
        r = H5Sselect_hyperslab(
          memshape, H5S_SELECT_SET, offset_mem, stride_mem, count2_mem, NULL);
      }
    }
    else {
      // we don't need a hyperslab here because we're writing 
      // a contiguous block from mem to disk with the same flat shape
    }
    //
    switch (finalData->GetDataType()) {
      case VTK_FLOAT:
        H5PartWriteDataArray(,H5T_NATIVE_FLOAT, this->H5FileId, name, finalData);
        break;
      case VTK_DOUBLE:
        H5PartWriteDataArray(,H5T_NATIVE_DOUBLE, this->H5FileId, name, finalData);
        break;
      case VTK_CHAR:
        if (VTK_TYPE_CHAR_IS_SIGNED) {
          H5PartWriteDataArray(,H5T_NATIVE_SCHAR, this->H5FileId, name, finalData);
        }
        else {
          H5PartWriteDataArray(,H5T_NATIVE_UCHAR, this->H5FileId, name, finalData);
        }
        break;
      case VTK_SIGNED_CHAR:
        H5PartWriteDataArray(,H5T_NATIVE_SCHAR, this->H5FileId, name, finalData);
        break;
      case VTK_UNSIGNED_CHAR:
        H5PartWriteDataArray(,H5T_NATIVE_UCHAR, this->H5FileId, name, finalData);
        break;
      case VTK_SHORT:
        H5PartWriteDataArray(,H5T_NATIVE_SHORT, this->H5FileId, name, finalData);
        break;
      case VTK_UNSIGNED_SHORT:
        H5PartWriteDataArray(,H5T_NATIVE_USHORT, this->H5FileId, name, finalData);
        break;
      case VTK_INT:
        H5PartWriteDataArray(,H5T_NATIVE_INT, this->H5FileId, name, finalData);
        break;
      case VTK_UNSIGNED_INT:
        H5PartWriteDataArray(,H5T_NATIVE_UINT, this->H5FileId, name, finalData);
        break;
      case VTK_LONG:
        H5PartWriteDataArray(,H5T_NATIVE_LONG, this->H5FileId, name, finalData);
        break;
      case VTK_UNSIGNED_LONG:
        H5PartWriteDataArray(,H5T_NATIVE_ULONG, this->H5FileId, name, finalData);
        break;
      case VTK_LONG_LONG:
        H5PartWriteDataArray(,H5T_NATIVE_LLONG, this->H5FileId, name, finalData);
        break;
      case VTK_UNSIGNED_LONG_LONG:
        H5PartWriteDataArray(,H5T_NATIVE_ULLONG, this->H5FileId, name, finalData);
        break;
      case VTK_ID_TYPE:
        if (VTK_SIZEOF_ID_TYPE==8) { 
          H5PartWriteDataArray(,H5T_NATIVE_LLONG, this->H5FileId, name, finalData);
        }
        else if (VTK_SIZEOF_ID_TYPE==4) {
          H5PartWriteDataArray(,H5T_NATIVE_LONG, this->H5FileId, name, finalData);
        }
        break;
      default:
        vtkErrorMacro(<<"Unexpected data type");
    }
    if (memshape!=H5S_ALL) {
      if (H5Sclose(memshape)<0) vtkErrorMacro(<<"memshape : HANDLE_H5S_CLOSE_ERR");
      memshape = H5S_ALL;
    }
    if (dataset>=0 && r<0) {
      vtkErrorMacro(<<"Array write failed for name "
        << name.c_str() << " Timestep " << this->H5FileId->timestep);
    } 
    else if (dataset>=0 && r>=0) {
      vtkDebugMacro(<<"Wrote " << name.c_str() << " " << Nt << " " << Nc << " " << typestring); 
    }
  }
}
//----------------------------------------------------------------------------
class H5PartWriterToleranceCheck: public std::binary_function<double, double, bool>
{
public:
  H5PartWriterToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(this->tolerance));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
void vtkH5PartWriter::WriteData()
{
#ifdef VTK_USE_MPI
  if (this->Controller) {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
  else {
    this->UpdatePiece = 0;
    this->UpdateNumPieces = 1;
  }
#else
  this->UpdatePiece = 0;
  this->UpdateNumPieces = 1;
#endif

  //
  // Make sure file is open
  //
  if (!this->OpenFile()) {
    vtkErrorMacro(<<"Couldn't create file " << this->FileName);
    return;
  }

  vtkInformation *outInfo = this->GetExecutive()->GetOutputInformation(0);

  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    this->TimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())[0];
    this->ActualTimeStep = vtkstd::find_if(
      this->InputTimeValues.begin(), this->InputTimeValues.end(), 
      vtkstd::bind2nd( H5PartWriterToleranceCheck( 1E-6 ), this->TimeValue)) 
      - this->InputTimeValues.begin();
    //
  }

  //
  // If this step exists then skip (not sure if we should keep this)
  //
  if (0 && H5PartHasStep(this->H5FileId, this->ActualTimeStep))
  {
    vtkErrorMacro(<<"Time Step already present in file. Aborting " 
      << this->FileName << " Step " << this->ActualTimeStep);
    return;
  }
  //
  // Set Step. This will create data group for us
  //
  H5PartSetStep(this->H5FileId, this->ActualTimeStep);
  //
  // Write out a TimeValue attribute for this step
  // only need do it on one process
  //
//  if (this->UpdatePiece==0) {
    h5part_int64_t attrib = H5PartWriteStepAttrib(
      this->H5FileId, "TimeValue", H5T_NATIVE_DOUBLE, &this->TimeValue, 1);
    if (attrib!=H5PART_SUCCESS) {
      vtkErrorMacro(<<"TimeValue attrib write failed ");
    }
//  }
  //
  // Get the input to write and Set Num-Particles
  // H5PartSetNumParticles does an MPI_All_Gather so we must do it
  // even if we are not writing anything out on this process
  // (ie if this process has zero particles we still do it)
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
  }
  else {
    this->WriteDataArray(0, NULL);
  }
  //
  // Write point data
  //
  int numScalars, numFound;
  vtkPointData *pd = input->GetPointData();
  numScalars = pd->GetNumberOfArrays();
  if (this->UpdateNumPieces>1 && !this->DisableInformationGather) {
    GatherScalarInfo(pd, numScalars, numFound);
    if (numScalars!=numFound) {
      vtkDebugMacro(<<"No scalars found, used MPI_Gather to find " << numFound << " Arrays");
    }
  }
  else numFound=numScalars;
  for (int i=0; i<numFound; i++) {
    vtkDataArray *data = pd->GetArray(i);
    this->WriteDataArray(i, data);
  }
  //
  // We are done.
  //
  vtkDebugMacro(<<"Time Step written " << this->FileName 
    << " Step " << this->ActualTimeStep
    << " Time " << this->TimeValue);
}
//----------------------------------------------------------------------------
bool vtkH5PartWriter::IsTimeStepPresent(int timestep)
{
  if (!this->OpenFile()) {
    vtkErrorMacro(<<"Couldn't open file " << this->FileName);
    return false;
  }
  return (H5PartHasStep(this->H5FileId, timestep)!=0);
}
//----------------------------------------------------------------------------
// Delete Timestep not working yet
//----------------------------------------------------------------------------
herr_t h5part_unlink(hid_t loc_id, const char *name, void *opdata)
{
  return H5Gunlink(loc_id,name);
}

void vtkH5PartWriter::DeleteTimeStep(int timestep)
{
  return; //because this does not work
  if (!this->OpenFile()) {
    vtkErrorMacro(<<"Couldn't create file " << this->FileName);
    return;
  }
  char name[256];
  hid_t timegroup = -1;
  sprintf(name,"Step#%u",timestep);
  sprintf(name,"Particles#%u",timestep);
  H5E_BEGIN_TRY {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    timegroup = H5Gopen(H5FileId->file, name, H5P_DEFAULT);
#else
    timegroup = H5Gopen(H5FileId->file, name);
#endif
  } H5E_END_TRY;
  if (timegroup>0) {
    H5Giterate(timegroup, "/", 0, h5part_unlink, NULL);
    H5Gunlink(timegroup,name);
  }
}
//----------------------------------------------------------------------------
void vtkH5PartWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " <<
    (this->FileName ? this->FileName : "(none)") << "\n";

  os << indent << "NumberOfSteps: " <<  this->NumberOfTimeSteps << "\n";
}
