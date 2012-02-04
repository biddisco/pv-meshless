/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5SPHReader.h
  Revision of last commit : $Rev: 501 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2008-03-11 20:17:29 +0100 #$

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
#include "vtkH5SPHReader.h"
//
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
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
#include "vtkCellArray.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtkstd/vector>
#include <stdexcept>
#include <algorithm>
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"

#include "FileSeriesFinder.h"

#ifdef VTK_USE_MPI
#include "vtkMultiProcessController.h"
#endif

//#define LIMIT_PARTITIONS 1024

#pragma warning( disable : 4996)

extern "C" hid_t h5tools_get_native_type(hid_t type);

//----------------------------------------------------------------------------
extern hid_t H5PartGetDiskShape(H5PartFile *f, hid_t dataset);
//----------------------------------------------------------------------------
#if 1 // def JB_DEBUG__
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
vtkCxxRevisionMacro(vtkH5SPHReader, "$Revision: 501 $");
vtkStandardNewMacro(vtkH5SPHReader);
//----------------------------------------------------------------------------
int H5DataTypeToVTKType(hid_t dataset_type)
{
  int vtktype = VTK_VOID;
  hid_t native_type = h5tools_get_native_type(dataset_type);
  //
  if (H5Tequal(native_type,H5T_NATIVE_FLOAT)) {
    vtktype = VTK_FLOAT;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_DOUBLE)) {
    vtktype = VTK_DOUBLE;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_SCHAR)) {
    vtktype = VTK_CHAR;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_UCHAR)) {
    vtktype = VTK_UNSIGNED_CHAR;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_SHORT)) {
    vtktype = VTK_SHORT;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_USHORT)) {
    vtktype = VTK_UNSIGNED_SHORT;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_INT)) {
    vtktype = VTK_INT;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_UINT)) {
    vtktype = VTK_UNSIGNED_INT;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_LONG)) {
    vtktype = VTK_LONG;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_ULONG)) {
    vtktype = VTK_UNSIGNED_LONG;
  }
  else if (H5Tequal(native_type,H5T_NATIVE_LLONG)) {
    vtktype = VTK_LONG_LONG; 
  }
  else if (H5Tequal(native_type,H5T_NATIVE_ULLONG)) {
    vtktype = VTK_UNSIGNED_LONG_LONG;
  }
  // H5T_NATIVE_LDOUBLE ?
  else {
    H5Tclose(native_type);
    throw std::runtime_error("An unexpected data type was encountered");
  }
  H5Tclose(native_type);
  return vtktype;
}
//----------------------------------------------------------------------------
vtkH5SPHReader::vtkH5SPHReader() 
{
  this->FileNamePattern         = NULL;
  this->StepNamePrefix          = NULL;
  this->SetStepNameFormat("Step",0);
  //
  this->Finder         = NULL;
  this->SetFileNamePattern("PREFIX TIME TEXT0 EXT");
}
//----------------------------------------------------------------------------
vtkH5SPHReader::~vtkH5SPHReader()
{ 
  this->CloseFile(); 
  if (this->FileNamePattern) {
    delete [] this->FileNamePattern;
    this->FileNamePattern = NULL;
  }
}
//----------------------------------------------------------------------------
bool vtkH5SPHReader::HasStep(int Step)
{
  if (!this->OpenFile()) return false;
  //
  char name[128];
  sprintf( name, "%s#%0*lld", this->StepNamePrefix, this->StepNameWidth, (long long) Step ); 
  // 1.8.0 herr_t herr = H5Oget_info_by_name(this->H5FileId, name, NULL, NULL);
  herr_t herr = H5Gget_objinfo( H5FileId->file, name, 1, NULL );
  return ( herr >= 0 );
}
//----------------------------------------------------------------------------
void vtkH5SPHReader::CloseFile()
{
  if (this->H5FileId) {
    H5PartCloseFile(this->H5FileId);
    this->H5FileId = NULL;
  }
}
//----------------------------------------------------------------------------
int vtkH5SPHReader::OpenFile()
{
  this->FileNameInternal = this->Finder->GenerateFileName(this->ActualTimeStep, 0, 0);
  //
  if (this->FileNameInternal.size()==0) {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
  }

  if (this->FileNameInternal!=this->FileNameLast) {
    this->CloseFile();
    this->FileNameLast = this->FileNameInternal;
  }

  if (!this->H5FileId) {
    this->H5FileId = H5PartOpenFile(this->FileNameInternal.c_str(), H5PART_READ);
    this->FileOpenedTime.Modified();
  }

  if (!this->H5FileId) {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileNameInternal);
    return 0;
  }

  return 1;
}
//----------------------------------------------------------------------------
int vtkH5SPHReader::IndexOfVectorComponent(const char *name)
{
  if (!this->CombineVectorComponents) return 0;
  //
  vtksys::RegularExpression re1(".*_([0-9]+)");
  if (re1.find(name)) {
    int index = atoi(re1.match(1).c_str());
    return index+1;
  }
  return 0;
}
//----------------------------------------------------------------------------
vtkstd::string vtkH5SPHReader::NameOfVectorComponent(const char *name)
{
  if (!this->CombineVectorComponents) return name;
  //
  vtksys::RegularExpression re1("(.*)_[0-9]+");
  if (re1.find(name)) {
    return re1.match(1);
  }
  return name;
}
//----------------------------------------------------------------------------
vtkIdType vtkH5SPHReader::GetNumberOfParticles()
{
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
	hid_t dataset_id = H5Dopen ( this->H5FileId->file, this->CompoundName.c_str(), H5P_DEFAULT );
#else
	hid_t dataset_id = H5Dopen ( this->H5FileId->file, this->CompoundName.c_str());
#endif
  if ( dataset_id < 0 ) return -1;

	hid_t space_id = H5Dget_space(dataset_id);
	if ( space_id < 0 ) return -1;

  vtkIdType nparticles = H5Sget_simple_extent_npoints ( space_id );

  if ( space_id != H5S_ALL ) {
		H5Sclose ( space_id );
	}
	H5Dclose ( dataset_id );
	return nparticles;
}
//----------------------------------------------------------------------------
struct _iter_op_data {
        int           stop_idx;
        int           count;
        int           type;
        char         *name;
        size_t        len;
        char         *pattern;
        vtkH5SPHReader::CompoundInfo *compounddata;
        vtkstd::string               *compundname;
        int                          *compoundsize;
};
//----------------------------------------------------------------------------
int vtkH5SPHReader::ScanCompoundType(hid_t loc_id, const char *name, void *opdata)
{
  H5G_stat_t     statbuf;
  H5T_class_t    classtype;
  hid_t          dataset_id, dataset_type, element_type;
  herr_t         err;
  char          *element_name;
  int            N, vtk_type;
  _iter_op_data *data = (_iter_op_data*)opdata;
  //
  H5Gget_objinfo(loc_id, name, 0, &statbuf);
  if (statbuf.type == H5G_DATASET) {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    dataset_id   = H5Dopen( loc_id, name, H5P_DEFAULT );
#else
    dataset_id   = H5Dopen( loc_id, name);
#endif
    dataset_type = H5Dget_type(dataset_id) ;
    classtype    = H5Tget_class( dataset_type  );
    N = 0;    
    if (classtype == H5T_COMPOUND) {
      data->count++;
      printf(" Object with name %s is compound \n", name);
      N = H5Tget_nmembers(dataset_type);
      *(data->compundname)  = name;
      *(data->compoundsize) = H5Tget_size(dataset_type) ;
      for (int i=0; i<N; i++) {
        element_name  = H5Tget_member_name(dataset_type, i);
        element_type  = H5Tget_member_type(dataset_type, i);
        int offset    = H5Tget_member_offset(dataset_type, i); 
        int size      = H5Tget_size(element_type); 
        vtk_type      = H5DataTypeToVTKType(element_type);
        printf(" \"%s\"; \"%s\";\n", element_name, vtkImageScalarTypeNameMacro(vtk_type));
        hid_t H5NativeType = h5tools_get_native_type(element_type);
        err = H5Tclose(element_type);
        data->compounddata->insert(
          CompoundType(element_name, 
          compound_info(vtk_type, H5NativeType, offset,size)));
        free(element_name);
      }
    }
    err = H5Tclose(dataset_type);
    err = H5Dclose(dataset_id);
    if (N>1) return 1; // stop further iteration
  }
  return 0;
}
//----------------------------------------------------------------------------
int vtkH5SPHReader::FindCompoundDataSet(
  hid_t group_id, const char *group_name, const hid_t type, char * const pattern) 
{
  struct _iter_op_data data;
  memset ( &data, 0, sizeof ( data ) );
  int idx           = 0;
  data.type         = type;
  data.pattern      = pattern;
  data.compounddata = &this->CompoundData;
  data.compundname  = &this->CompoundName;
  data.compoundsize = &this->CompoundSize;
  //
  h5part_int64_t herr = H5Giterate( group_id, group_name, &idx, ScanCompoundType, &data );
  if (herr<0) return herr;
  return data.count;
}
//----------------------------------------------------------------------------
int vtkH5SPHReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

#ifdef VTK_USE_MPI
  if (this->Controller) {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#else 
  this->UpdatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = 1;
#endif

  if (!this->Finder) {   
    this->Finder = new FileSeriesFinder(this->FileNamePattern);
    this->Finder->SetPrefixRegEx("(.*/)");
    this->Finder->SetTimeRegEx("([0-9]+)");
    this->Finder->Scan(this->FileName);
    this->Finder->GetTimeValues(this->TimeStepValues);
    this->NumberOfTimeSteps = this->Finder->GetNumberOfTimeSteps();
  }

  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime || !this->H5FileId);

  if (!this->OpenFile()) {
    return 0;
  }

  if (1 || NeedToReadInformation) {
    int nds = FindCompoundDataSet(this->H5FileId->file, "/", H5G_DATASET, NULL);
    for (CompoundInfo::iterator it=this->CompoundData.begin(); it!=this->CompoundData.end(); ++it) {
      this->PointDataArraySelection->AddArray(it->first.c_str());
    }
      
    // if TIME information was either not present ot not consistent, then 
    // set something so that consumers of this data can iterate sensibly
    if (this->NumberOfTimeSteps>0) {
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
        &this->TimeStepValues[0], this->TimeStepValues.size());
      double timeRange[2];
      timeRange[0] = this->TimeStepValues.front();
      timeRange[1] = this->TimeStepValues.back();
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
  }
  //
  return 1;
}
//----------------------------------------------------------------------------
//int vtkH5SPHReader::ReadArrayFromCompound(
//  hid_t loc_id, vtkDataArray *data, hid_t h5type, char *name)
//{ 
//}
//----------------------------------------------------------------------------
/*
template <class T1, class T2>
void CopyIntoCoords_T(int offset, vtkDataArray *source, vtkDataArray *dest)
{
  int N = source->GetNumberOfTuples();
  T1 *sptr = static_cast<T1*>(source->GetVoidPointer(0));
  T2 *dptr = static_cast<T2*>(dest->WriteVoidPointer(0,N)) + offset;
  for (int i=0; i<N; ++i) {
    *dptr = *sptr++;
    dptr += 3;
  }
}
//----------------------------------------------------------------------------
void vtkH5SPHReader::CopyIntoCoords(int offset, vtkDataArray *source, vtkDataArray *dest)
{
  switch (source->GetDataType()) 
  {
    case VTK_CHAR:
    case VTK_SIGNED_CHAR:
    case VTK_UNSIGNED_CHAR:
      CopyIntoCoords_T<char,float>(offset, source, dest);
      break;
    case VTK_SHORT:
      CopyIntoCoords_T<short int,float>(offset, source, dest);
      break;
    case VTK_UNSIGNED_SHORT:
      CopyIntoCoords_T<unsigned short int,float>(offset, source, dest);
      break;
    case VTK_INT:
      CopyIntoCoords_T<int,float>(offset, source, dest);
      break;
    case VTK_UNSIGNED_INT:
      CopyIntoCoords_T<unsigned int,float>(offset, source, dest);
      break;
    case VTK_LONG:
      CopyIntoCoords_T<long int,float>(offset, source, dest);
      break;
    case VTK_UNSIGNED_LONG:
      CopyIntoCoords_T<unsigned long int,float>(offset, source, dest);
      break;
    case VTK_LONG_LONG:
      CopyIntoCoords_T<long long,float>(offset, source, dest);
      break;
    case VTK_UNSIGNED_LONG_LONG:
      CopyIntoCoords_T<unsigned long long,float>(offset, source, dest);
      break;
    case VTK_FLOAT:
      CopyIntoCoords_T<float,float>(offset, source, dest);
      break;
    case VTK_DOUBLE:
      CopyIntoCoords_T<double,float>(offset, source, dest);
      break;
    case VTK_ID_TYPE:
      CopyIntoCoords_T<vtkIdType,float>(offset, source, dest);
      break;
    default:
      break;
      vtkErrorMacro(<<"Unexpected data type");
  }
}
*/
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class H5SPH_WithinTolerance: public std::binary_function<double, double, bool>
{
public:
  H5SPH_WithinTolerance(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(a*this->tolerance));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
int vtkH5SPHReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  //
  // Get the parallel Id for partitioning
  //
  this->UpdatePiece = -1;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()))
    {
    this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    }
  else
    {
    this->UpdatePiece = 0;
    }
#ifdef VTK_USE_MPI
  if (this->Controller && 
      (this->UpdatePiece != this->Controller->GetLocalProcessId() ||
       this->UpdateNumPieces != this->Controller->GetNumberOfProcesses()))
  {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#endif    

  //
  typedef vtkstd::map< vtkstd::string, vtkstd::vector<vtkstd::string> > FieldMap;
  FieldMap scalarFields;
  //
//  if (this->TimeStepValues.size()==0) return 0;
  //
  // Make sure that the user selected arrays for coordinates are represented
  //
  vtkstd::vector<vtkstd::string> coordarrays(3, "");
  //
  int N = this->PointDataArraySelection->GetNumberOfArrays();
  for (int i=0; i<N; i++) {
    const char *name = this->PointDataArraySelection->GetArrayName(i);
    // Do we want to load this array
    bool processarray = false;
    if (!vtksys::SystemTools::Strucmp(name,this->Xarray)) {
      processarray = true;
      coordarrays[0] = name;
    }
    else if (!vtksys::SystemTools::Strucmp(name,this->Yarray)) {
      processarray = true;
      coordarrays[1] = name;
    }
    else if (!vtksys::SystemTools::Strucmp(name,this->Zarray)) {
      processarray = true;
      coordarrays[2] = name;
    }
    else if (this->PointDataArraySelection->ArrayIsEnabled(name)) {
      processarray = true;
    }
    if (!processarray)  continue;

    // make sure we cater for multi-component vector fields
    int vectorcomponent;
    if ((vectorcomponent=this->IndexOfVectorComponent(name))>0) {
      vtkstd::string vectorname = this->NameOfVectorComponent(name) + "_v";
      FieldMap::iterator pos = scalarFields.find(vectorname);
      if (pos==scalarFields.end()) {
        vtkstd::vector<vtkstd::string> arraylist(1, name);
        FieldMap::value_type element(vectorname, arraylist);
        scalarFields.insert(element);
      }
      else {
        pos->second.reserve(vectorcomponent);
        pos->second.resize(vtkstd::max((int)(pos->second.size()), vectorcomponent));
        pos->second[vectorcomponent-1] = name;
      }
    }
    else {
      vtkstd::vector<vtkstd::string> arraylist(1, name);
      FieldMap::value_type element(name, arraylist);
      scalarFields.insert(element);
    }
  }
  //
  FieldMap::iterator coordvector=scalarFields.end();
  for (FieldMap::iterator pos=scalarFields.begin(); pos!=scalarFields.end(); ++pos)
  {
    if (pos->second.size()==3 &&
       (pos->second[0]==coordarrays[0]) && 
       (pos->second[1]==coordarrays[1]) && 
       (pos->second[2]==coordarrays[2])) 
    {
      // change the name of this entry to "coords" to ensure we use it as such
      FieldMap::value_type element("Coords", pos->second);
      scalarFields.erase(pos);
      coordvector = scalarFields.insert(element).first;
      break;
    }
  }
  if (coordvector==scalarFields.end()) {
    FieldMap::value_type element("Coords", coordarrays);
    scalarFields.insert(element);
  }

  //
  // Get the TimeStep Requested from the information if present
  //
  vtkOStreamWrapper::EndlType endl;
  vtkOStreamWrapper::UseEndl(endl);

  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())[0];
    this->ActualTimeStep = vtkstd::find_if(
      this->TimeStepValues.begin(), 
      this->TimeStepValues.end(), 
      vtkstd::bind2nd( H5SPH_WithinTolerance(0.001), requestedTimeValue )) 
      - this->TimeStepValues.begin();
    this->ActualTimeStep = this->ActualTimeStep;
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), &requestedTimeValue, 1);
  }
  else 
  {
    double timevalue[1];
    unsigned int index = this->ActualTimeStep;
    if (index<this->TimeStepValues.size()) {
      timevalue[0] = this->TimeStepValues[index];
    }
    else {
      timevalue[0] = this->TimeStepValues[0];
    }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), &timevalue[0], 1);
  }

  // timestep may have changed, so reopen file.
  this->OpenFile();

//    this->UpdatePiece = 1;

  //
  // Get the number of particles for this timestep
  //
  vtkIdType Nparticles = this->GetNumberOfParticles();
  //
  // Split particles up per process for parallel load
  //
  vtkExtentTranslator *extTran = vtkExtentTranslator::New();
  extTran->SetSplitModeToBlock();
  int WholeExtent[6] = { 0, Nparticles, 0, 0, 0, 0 };
#if !defined(LIMIT_PARTITIONS)
  extTran->SetNumberOfPieces(this->UpdateNumPieces);
#else 
  extTran->SetNumberOfPieces(LIMIT_PARTITIONS);
#endif
  extTran->SetPiece(this->UpdatePiece);
  extTran->SetWholeExtent(WholeExtent);
  extTran->PieceToExtent();
  int PartitionExtents[6];
  extTran->GetExtent(PartitionExtents);
  extTran->Delete();
  int ParticleStart = PartitionExtents[0];
  int ParticleEnd   = PartitionExtents[1]-1;
  vtkIdType Nt = ParticleEnd - ParticleStart + 1;


  //
  //
  //
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
	hid_t dataset_id = H5Dopen (this->H5FileId->file, this->CompoundName.c_str(), H5P_DEFAULT);
#else
	hid_t dataset_id = H5Dopen (this->H5FileId->file, this->CompoundName.c_str());
#endif

  // Setup arrays for reading data
  vtkSmartPointer<vtkPoints>    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDataArray> coords = NULL;
  int c, datatype = 0;
  for (FieldMap::iterator it=scalarFields.begin(); it!=scalarFields.end(); it++) {
    // use the type of the first array for all if it is a vector field
    vtkstd::vector<vtkstd::string> &arraylist = (*it).second;
    const char *name = arraylist[0].c_str();
    vtkstd::string rootname = this->NameOfVectorComponent(name);
    int Nc = arraylist.size();
    //
    vtkSmartPointer<vtkDataArray> dataarray = NULL;
    int   vtk_type = this->CompoundData[name].vtk_type;
    hid_t datatype = this->CompoundData[name].H5DataType;
    int       size = this->CompoundData[name].size;
    //
    dataarray.TakeReference(vtkDataArray::SafeDownCast(vtkDataArray::CreateArray(vtk_type)));
    dataarray->SetNumberOfComponents(Nc);
    dataarray->SetNumberOfTuples(Nt);
    dataarray->SetName(rootname.c_str());
    hsize_t     count1_mem[] = { Nt };                       
    hsize_t     count2_mem[] = { Nt };                          
    hsize_t     offset_mem[] = { 0 };                           
    hsize_t     stride_mem[] = { 1 };                          
    int structsize = 0;
    int offset     = 0;
    for (c=0; c<Nc; c++) {                                      
      const char *cname = arraylist[c].c_str();
      structsize += this->CompoundData[arraylist[c]].size;
    }
    // Create a data type for this component
    hid_t struct_id = H5Tcreate(H5T_COMPOUND, structsize);
    for (c=0; c<Nc; c++) {                                      
      const char *cname = arraylist[c].c_str();
      H5Tinsert(struct_id, arraylist[c].c_str(), offset, this->CompoundData[arraylist[c]].H5DataType);
      offset += this->CompoundData[arraylist[c]].size;
    }

    // the memory space is easy...we read data into a contiguous chunk of length Nt
    hid_t memspace  = H5Screate_simple(1, count1_mem, NULL);   

    // the disk space is more complicated. We must read Nt from the 
    // part of the file offset by ParticleStart 
    hid_t diskshape = H5PartGetDiskShape(H5FileId, dataset_id);

    offset_mem[0] = ParticleStart;
    herr_t err = H5Sselect_hyperslab(diskshape, H5S_SELECT_SET,
      offset_mem, stride_mem, count2_mem, NULL);  

      err = H5Dread(dataset_id, struct_id, memspace, diskshape, H5P_DEFAULT, dataarray->GetVoidPointer(0));

//      H5Dread(dataset, null##T2, memspace,                    
//        diskshape, H5P_DEFAULT, dataarray->GetVoidPointer(0));  
      if (memspace!=H5S_ALL) H5Sclose(memspace);                
      if (diskshape!=H5S_ALL) H5Sclose(diskshape);              


//    herr_t err = H5Dread(dataset_id, struct_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataarray->GetVoidPointer(0));
    if (err<0) {
      vtkErrorMacro(<<"Failed to read H5T_Compound data " << arraylist[0].c_str()); 
    }
    H5Tclose(struct_id);
/*
      hid_t diskshape = H5PartGetDiskShape(H5FileId,dataset_id);   
      hid_t memspace  = H5Screate_simple(1, count1_mem, NULL);   
      offset_mem[0] = c;                                        
      r = H5Sselect_hyperslab(                                  
        memspace, H5S_SELECT_SET,                               
        offset_mem, stride_mem, count2_mem, NULL);  

//      H5Dread(dataset, null##T2, memspace,                    
//        diskshape, H5P_DEFAULT, dataarray->GetVoidPointer(0));  
      if (memspace!=H5S_ALL) H5Sclose(memspace);                
      if (diskshape!=H5S_ALL) H5Sclose(diskshape);              
    }
*/
    //
    if (dataarray) {
      if ((*it).first=="Coords") coords = dataarray;
      else {
        output->GetPointData()->AddArray(dataarray);
        if (!output->GetPointData()->GetScalars()) {
          output->GetPointData()->SetActiveScalars(dataarray->GetName());
        }
      }
    }
  }
  H5Dclose(dataset_id);
  //
  if (this->GenerateVertexCells) {
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = vertices->WritePointer(Nt, Nt*2);
    for (vtkIdType i=0; i<Nt; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }
    output->SetVerts(vertices);
  }
  //
  points->SetData(coords);
  output->SetPoints(points);
  return 1;
}
//----------------------------------------------------------------------------
void vtkH5SPHReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
