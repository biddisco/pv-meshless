/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5PartReader.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

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
#include "vtkH5PartReader.h"
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
#include "vtkOutlineSource.h"
#include "vtkAppendPolyData.h"
#include "vtkBoundingBox.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtkstd/vector>
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
//
#include "vtkDummyController.h"
//
#include <algorithm>

#include "H5Part.h"
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkH5PartReader, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
/*!
  \ingroup h5part_utility

  This function can be used to query the Type of a dataset
  It is not used by the core H5Part library but is useful when
  reading generic data from the file.
  An example of usage would be (H5Tequal(datatype,H5T_NATIVE_FLOAT))
  any NATIVE type can be used to test.

  \return  \c an hdf5 handle to the native type of the data
*/
static hid_t H5PartGetNativeDatasetType(H5PartFile *f, const char *name)
{
  hid_t dataset, datatype, datatypen;
  if (!f->timegroup)
    {
    H5PartSetStep(f,f->timestep); /* choose current step */
    }
#if (!H5_USE_16_API && ((H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))))
  dataset=H5Dopen(f->timegroup, name, H5P_DEFAULT);
#else
  dataset=H5Dopen(f->timegroup, name);
#endif
  datatype  = H5Dget_type(dataset);
  datatypen = H5Tget_native_type(datatype, H5T_DIR_DEFAULT);
  H5Tclose(datatype);
  H5Dclose(dataset);
  return datatypen;
}

//----------------------------------------------------------------------------
hid_t H5PartGetDiskShape(H5PartFile *f, hid_t dataset)
{
  hid_t space = H5Dget_space(dataset);
  if (H5PartHasView(f))
    {
    int r;
    hsize_t stride, count;
    hsize_t range[2];
    /* so, is this selection inclusive or exclusive? */
    range[0]=f->viewstart;
    range[1]=f->viewend;
    count = range[1]-range[0]; /* to be inclusive */
    stride=1;
    /* now we select a subset */
    if (f->diskshape>0)
      {
      r=H5Sselect_hyperslab(f->diskshape,H5S_SELECT_SET,
        range/* only using first element */,
        &stride,&count,NULL);
      }
    /* now we select a subset */
    r=H5Sselect_hyperslab(space,H5S_SELECT_SET,
      range,&stride,&count,NULL);
    if (r<0)
      {
      fprintf(stderr,"Abort: Selection Failed!\n");
      return space;
      }
    }
  return space;
}
//----------------------------------------------------------------------------
//#define JB_DEBUG__
#ifdef JB_DEBUG__
  #define OUTPUTTEXT(a) std::cout << (a) << std::endl; std::cout.flush();

    #undef vtkDebugMacro
    #define vtkDebugMacro(a)  \
    { \
      if (this->UpdatePiece>=0) { \
        vtkOStreamWrapper::EndlType endl; \
        vtkOStreamWrapper::UseEndl(endl); \
        vtkOStrStreamWrapper vtkmsg; \
        vtkmsg << "P(" << this->UpdatePiece << "): " a << "\n"; \
        OUTPUTTEXT(vtkmsg.str()); \
        vtkmsg.rdbuf()->freeze(0); \
      } \
    }

  #undef vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkH5PartReader);
//----------------------------------------------------------------------------
vtkH5PartReader::vtkH5PartReader()
{
  this->SetNumberOfInputPorts(0);
  //
  this->NumberOfTimeSteps        = 0;
  this->TimeStep                 = 0;
  this->ActualTimeStep           = 0;
  this->TimeStepTolerance        = 1E-6;
  this->CombineVectorComponents  = 1;
  this->GenerateVertexCells      = 0;
  this->FileName                 = NULL;
  this->H5FileId                 = NULL;
  this->Xarray                   = NULL;
  this->Yarray                   = NULL;
  this->Zarray                   = NULL;
  this->UpdatePiece              = 0;
  this->UpdateNumPieces          = 0;
  this->TimeOutOfRange           = 0;
  this->MaskOutOfTimeRangeOutput = 0;
  this->IntegerTimeStepValues    = 0;
  this->ExportPartitionBoxes     = 0;
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  this->SetXarray("Coords_0");
  this->SetYarray("Coords_1");
  this->SetZarray("Coords_2");
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkH5PartReader::~vtkH5PartReader()
{
  this->CloseFile();
  delete [] this->FileName;
  this->FileName = NULL;

  delete [] this->Xarray;
  this->Xarray = NULL;

  delete [] this->Yarray;
  this->Yarray = NULL;

  delete [] this->Zarray;
  this->Zarray = NULL;

  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = 0;

  this->SetController(NULL);
}
//----------------------------------------------------------------------------
bool vtkH5PartReader::HasStep(int Step)
{
  if (!this->OpenFile())
    {
    return false;
    }

  if (H5PartHasStep(this->H5FileId, Step))
    {
    return true;
    }
  return false;
}
//----------------------------------------------------------------------------
void vtkH5PartReader::SetFileName(char *filename)
{
  if (this->FileName == NULL && filename == NULL)
    {
    return;
    }
  if (this->FileName && filename && (!strcmp(this->FileName,filename)))
    {
    return;
    }
  delete [] this->FileName;
  this->FileName = NULL;

  if (filename)
    {
    this->FileName = vtksys::SystemTools::DuplicateString(filename);
    this->SetFileModified();
    }
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkH5PartReader::SetFileModified()
{
  this->FileModifiedTime.Modified();
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkH5PartReader::CloseFile()
{
  
  if (this->H5FileId != NULL)
    {
    H5PartCloseFile(this->H5FileId);
    this->H5FileId = NULL;
    }
}
//----------------------------------------------------------------------------
void vtkH5PartReader::CloseFileIntermediate()
{
}
//----------------------------------------------------------------------------
int vtkH5PartReader::OpenFile()
{
  if (!this->FileName)
    {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
    }

  if (FileModifiedTime>FileOpenedTime)
    {
    this->CloseFile();
    }

  if (!this->H5FileId)
    {
    this->H5FileId = H5PartOpenFile(this->FileName, H5PART_READ);
    this->FileOpenedTime.Modified();
    }

  if (!this->H5FileId)
    {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
    }

  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartReader::IndexOfVectorComponent(const char *name)
{
  if (!this->CombineVectorComponents)
    {
    return 0;
    }
  //
  vtksys::RegularExpression re1(".*_([0-9]+)");
  if (re1.find(name))
    {
    int index = atoi(re1.match(1).c_str());
    return index+1;
    }
  return 0;
}
//----------------------------------------------------------------------------
vtkstd::string vtkH5PartReader::NameOfVectorComponent(const char *name)
{
  if (!this->CombineVectorComponents)
    {
    return name;
    }
  //
  vtksys::RegularExpression re1("(.*)_[0-9]+");
  if (re1.find(name))
    {
    return re1.match(1);
    }
  return name;
}
//----------------------------------------------------------------------------
int vtkH5PartReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime || !this->H5FileId);

  if (NeedToReadInformation)
    {
    if (!this->OpenFile())
      {
      return 0;
      }

    this->NumberOfTimeSteps = H5PartGetNumSteps(this->H5FileId);
    H5PartSetStep(this->H5FileId, 0);
    int nds = H5PartGetNumDatasets(this->H5FileId);
    char name[512];
    for (int i=0; i<nds; i++)
      {
      // return 0 for no, 1,2,3,4,5 etc for index (1 based offset)
      H5PartGetDatasetName(this->H5FileId, i, name, 512);
      this->PointDataArraySelection->AddArray(name);
      }

    this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);
    int validTimes = 0;
    for (int i=0; i<this->NumberOfTimeSteps; ++i)
      {
      H5PartSetStep(this->H5FileId, i);
      // Get the time value if it exists
      h5part_int64_t numAttribs = H5PartGetNumStepAttribs(this->H5FileId);
      if (numAttribs>0)
        {
        char           attribName[128];
        h5part_int64_t attribNameLength = 128;
        h5part_int64_t attribType       = 0;
        h5part_int64_t attribNelem      = 0;
        for (h5part_int64_t a=0; a<numAttribs; a++)
          {
          h5part_int64_t status = H5PartGetStepAttribInfo (
            this->H5FileId, a, attribName, attribNameLength,
            &attribType, &attribNelem);
          if (status==H5PART_SUCCESS && !strncmp("TimeValue",attribName,attribNameLength))
            {
            if (H5Tequal(attribType,H5T_NATIVE_DOUBLE) && attribNelem==1)
              {
              status=H5PartReadStepAttrib(this->H5FileId, attribName, &this->TimeStepValues[i]);
              if (status==H5PART_SUCCESS)
                {
                validTimes++;
                }
              }
            }
          }
        }
      }
    H5PartSetStep(this->H5FileId, 0);

    if (this->NumberOfTimeSteps==0)
      {
      vtkErrorMacro(<<"No time steps in data");
      return 0;
      }

    // if TIME information was either not present ot not consistent, then
    // set something so that consumers of this data can iterate sensibly
    if (this->IntegerTimeStepValues || (this->NumberOfTimeSteps>0 && this->NumberOfTimeSteps!=validTimes))
      {
      for (int i=0; i<this->NumberOfTimeSteps; i++)
        {
        // insert read of Time array here
        this->TimeStepValues[i] = i;
        }
      }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
      &this->TimeStepValues[0],
      static_cast<int>(this->TimeStepValues.size()));
    double timeRange[2];
    timeRange[0] = this->TimeStepValues.front();
    timeRange[1] = this->TimeStepValues.back();
    if (this->TimeStepValues.size()>1)
      {
      this->TimeStepTolerance = 0.01*(this->TimeStepValues[1]-this->TimeStepValues[0]);
      }
    else
      {
      this->TimeStepTolerance = 1E-3;
      }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }

  this->CloseFileIntermediate();

  return 1;
}
//----------------------------------------------------------------------------
int GetVTKDataType(int datatype)
{
  if (H5Tequal(datatype,H5T_NATIVE_FLOAT))
    {
    return VTK_FLOAT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_DOUBLE))
    {
    return VTK_DOUBLE;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_SCHAR))
    {
    return VTK_CHAR;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_UCHAR))
    {
    return VTK_UNSIGNED_CHAR;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_SHORT))
    {
    return VTK_SHORT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_USHORT))
    {
    return VTK_UNSIGNED_SHORT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_INT))
    {
    return VTK_INT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_UINT))
    {
    return VTK_UNSIGNED_INT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_LONG))
    {
    return VTK_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_ULONG))
    {
    return VTK_UNSIGNED_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_LLONG))
    {
    return VTK_LONG_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_ULLONG))
    {
    return VTK_UNSIGNED_LONG_LONG;
    }
  return VTK_VOID;
}

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
void vtkH5PartReader::CopyIntoCoords(int offset, vtkDataArray *source, vtkDataArray *dest)
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
/*
vtkstd::pair<double, double> GetClosest(vtkstd::vector<double> &sortedlist, const double& val) const
{
  std::vector<double>::const_iterator it = std::lower_bound(sortedlist.begin(), sortedlist.end(), val);
  if (it == sortedlist.end())        return std::make_pair(sortedlist.back(), sortedlist.back());
  else if (it == sortedlist.begin()) return std::make_pair(sortedlist.front(), sortedlist.front());
  else return std::make_pair(*(it - 1), *(it));
}
*/
class H5PartToleranceCheck: public std::binary_function<double, double, bool>
{
public:
  H5PartToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(this->tolerance));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
#if (!H5_USE_16_API && ((H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))))
  #define h_params ,H5P_DEFAULT
#else
  #define h_params
#endif
//----------------------------------------------------------------------------
int vtkH5PartReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData     *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  typedef vtkstd::map< vtkstd::string, vtkstd::vector<vtkstd::string> > FieldMap;
  FieldMap scalarFields;
  //
  if (this->TimeStepValues.size()==0) return 0;
  //
  // Make sure that the user selected arrays for coordinates are represented
  //
  vtkstd::vector<vtkstd::string> coordarrays(3, "");
  //
  int N = this->PointDataArraySelection->GetNumberOfArrays();
  for (int i=0; i<N; i++)
    {
    const char *name = this->PointDataArraySelection->GetArrayName(i);
    // Do we want to load this array
    bool processarray = false;
    if (!vtksys::SystemTools::Strucmp(name,this->Xarray))
      {
      processarray = true;
      coordarrays[0] = name;
      }
    if (!vtksys::SystemTools::Strucmp(name,this->Yarray))
      {
      processarray = true;
      coordarrays[1] = name;
      }
    if (!vtksys::SystemTools::Strucmp(name,this->Zarray))
      {
      processarray = true;
      coordarrays[2] = name;
      }
    if (this->PointDataArraySelection->ArrayIsEnabled(name))
      {
      processarray = true;
      }
    if (!processarray)
      {
      continue;
      }

    // make sure we cater for multi-component vector fields
    int vectorcomponent;
    if ((vectorcomponent=this->IndexOfVectorComponent(name))>0)
      {
      vtkstd::string vectorname = this->NameOfVectorComponent(name) + "_v";
      FieldMap::iterator pos = scalarFields.find(vectorname);
      if (pos==scalarFields.end())
        {
        vtkstd::vector<vtkstd::string> arraylist(1, name);
        FieldMap::value_type element(vectorname, arraylist);
        scalarFields.insert(element);
        }
      else
        {
        pos->second.reserve(vectorcomponent);
        pos->second.resize(vtkstd::max((int)(pos->second.size()), vectorcomponent));
        pos->second[vectorcomponent-1] = name;
        }
      }
    else
      {
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

  if (coordvector==scalarFields.end())
    {
    FieldMap::value_type element("Coords", coordarrays);
    scalarFields.insert(element);
    }

  //
  // Get the TimeStep Requested from the information if present
  //
  this->TimeOutOfRange = 0;
  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())[0];
    this->ActualTimeStep = vtkstd::find_if(
      this->TimeStepValues.begin(), this->TimeStepValues.end(),
      vtkstd::bind2nd( H5PartToleranceCheck( 
          this->IntegerTimeStepValues ? 0.5 : this->TimeStepTolerance ), requestedTimeValue ))
      - this->TimeStepValues.begin();
    //
    if (requestedTimeValue<this->TimeStepValues.front() || requestedTimeValue>this->TimeStepValues.back())
      {
      this->TimeOutOfRange = 1;
      }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), &requestedTimeValue, 1);
    }
  else
    {
    double timevalue[1];
    unsigned int index = this->ActualTimeStep;
    if (index<this->TimeStepValues.size())
      {
      timevalue[0] = this->TimeStepValues[index];
      }
    else
      {
      timevalue[0] = this->TimeStepValues[0];
      }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), &timevalue[0], 1);
    }

  if (this->TimeOutOfRange && this->MaskOutOfTimeRangeOutput)
    {
    // don't do anything, just return success
    return 1;
    }

  // open the file if not already done
  if (!this->OpenFile())
    {
    return 0;
    }

  // Set the TimeStep on the H5 file
  H5PartSetStep(this->H5FileId, this->ActualTimeStep);
  //
  // Get the number of particles for this timestep
  //
  vtkIdType Nparticles = H5PartGetNumParticles(this->H5FileId);

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
  vtkDebugMacro(<< "Reading Particles " << ParticleStart << " : " << ParticleEnd << " = " << Nt);

  // Setup arrays for reading data
  vtkSmartPointer<vtkPoints>    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDataArray> coords = NULL;
  for (FieldMap::iterator it=scalarFields.begin(); it!=scalarFields.end(); it++)
    {
    // use the type of the first array for all if it is a vector field
    vtkstd::vector<vtkstd::string> &arraylist = (*it).second;
    const char *array_name = arraylist[0].c_str();
    vtkstd::string rootname = this->NameOfVectorComponent(array_name);
    int Nc = static_cast<int>(arraylist.size());
    //
    vtkSmartPointer<vtkDataArray> dataarray = NULL;
    hid_t datatype = H5PartGetNativeDatasetType(H5FileId,array_name);
    int vtk_datatype = GetVTKDataType(datatype);

    if (vtk_datatype != VTK_VOID)
      {
      dataarray.TakeReference(vtkDataArray::CreateDataArray(vtk_datatype));
      dataarray->SetNumberOfComponents(Nc);
      dataarray->SetNumberOfTuples(Nt);
      dataarray->SetName(rootname.c_str());

      // now read the data components.
      herr_t r;
      hsize_t count1_mem[] = { Nt*Nc };
      hsize_t count2_mem[] = { Nt };
      hsize_t offset_mem[] = { 0 };
      hsize_t stride_mem[] = { Nc };
      hsize_t     count1_dsk[] = { Nt };
      hsize_t     offset_dsk[] = { ParticleStart };
      hsize_t     stride_dsk[] = { 1 };
      //
      for (int c=0; c<Nc; c++)
        {
        const char *name = arraylist[c].c_str();
        hid_t dataset   = H5Dopen(H5FileId->timegroup, name h_params);
        hid_t diskshape = H5PartGetDiskShape(H5FileId,dataset);
        /* parallel read needs hyperslab for disk */
        r = H5Sselect_hyperslab(diskshape, H5S_SELECT_SET,
          offset_dsk, stride_dsk, count1_dsk, NULL);
        hid_t memspace = H5Screate_simple(1, count1_mem, NULL);
        hid_t component_datatype = H5PartGetNativeDatasetType(H5FileId, name);
        /* read x/y/z arrays into strided mem - use hyperslab */
        offset_mem[0] = c;
        r = H5Sselect_hyperslab(
          memspace, H5S_SELECT_SET,
          offset_mem, stride_mem, count2_mem, NULL);
        if (component_datatype == datatype)
          {
          H5Dread(dataset, datatype, memspace,
            diskshape, H5P_DEFAULT, dataarray->GetVoidPointer(0));
          }
        else
          {
          // read data into a temporary array of the right type and then copy it
          // over to the "dataarray".
          // This can be optimized to create a single component array. But I
          // don't understand the stride/offset stuff too well to fix that.
          vtkDataArray* temparray =
            vtkDataArray::CreateDataArray(GetVTKDataType(component_datatype));
          temparray->SetNumberOfComponents(Nc);
          temparray->SetNumberOfTuples(Nt);
          r = H5Sselect_hyperslab(
            memspace, H5S_SELECT_SET,
            offset_mem, stride_mem, count2_mem, NULL);
          H5Dread(dataset, component_datatype, memspace,
            diskshape, H5P_DEFAULT, temparray->GetVoidPointer(0));
          dataarray->CopyComponent(c, temparray, c);
          temparray->Delete();
          }
        if (memspace!=H5S_ALL)
          {
          H5Sclose(memspace);
          }
        if (diskshape!=H5S_ALL)
          {
          H5Sclose(diskshape);
          }
        H5Dclose(dataset);
        }
      }
    else
      {
      H5Tclose(datatype);
      vtkErrorMacro("An unexpected data type was encountered");
      return 0;
      }
    H5Tclose(datatype);
    //
    if (dataarray)
      {
      if ((*it).first=="Coords") {
        coords = dataarray;
        coords->SetName("Coordinates");
      }
      else
        {
        output->GetPointData()->AddArray(dataarray);
        if (!output->GetPointData()->GetScalars())
          {
          output->GetPointData()->SetActiveScalars(dataarray->GetName());
          }
        }
      }
    }

  //
  // generate cells
  //
  if (this->GenerateVertexCells)
    {
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = vertices->WritePointer(Nt, 2*Nt);
    for (vtkIdType i=0; i<Nt; ++i)
      {
      cells[2*i] = 1;
      cells[2*i+1] = i;
      }
    output->SetVerts(vertices);
    }
  //
  //
  //
  if (this->ExportPartitionBoxes) {
    int newpoints = this->ReadBoundingBoxes(coords, output, ParticleStart, ParticleEnd);
  }
  //
  //
  //
  points->SetData(coords);
  output->SetPoints(points);
  //
  //
  // only subclasses actually close the file.
  //
  this->CloseFileIntermediate();

  return 1;
}
//----------------------------------------------------------------------------
int vtkH5PartReader::ReadBoundingBoxes(vtkDataArray *coords, vtkPolyData *output, vtkIdType extent0, vtkIdType extent1)
{
   H5E_auto2_t  errfunc;
   void        *errdata;
   vtkIdType    N2 = 0;

   // Prevent HDF5 to print out handled errors, first save old error handler
   H5Eget_auto(H5E_DEFAULT, &errfunc, &errdata);
   // Turn off error handling
   H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // @TODO, use group/name string passed into reader ...
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t partitiongroup = H5Gopen(H5FileId->file, "Partition#0", H5P_DEFAULT );
	hid_t dataset_id = (partitiongroup>0) ? H5Dopen ( partitiongroup, "Box", H5P_DEFAULT ) : -1;
#else
  hid_t partitiongroup = H5Gopen(H5FileId->file, "Partition#0");
	hid_t dataset_id = (partitiongroup>0) ? H5Dopen ( partitiongroup, "Box") : -1;
#endif
  if (partitiongroup>0 && dataset_id>0) { 
    // the disk space is more complicated. We must read Nt from the 
    // part of the file offset by ParticleStart 
    hid_t diskshape = H5PartGetDiskShape(H5FileId, dataset_id);
    hsize_t dims[2]; 
    hsize_t maxdims[2];
    int val = H5Sget_simple_extent_dims(diskshape, dims, maxdims);
    int partitions = dims[0]/13;
    //
    std::vector<double> data(dims[0], 0.0);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    vtkIdType N1 = coords->GetNumberOfTuples();
    N2 = partitions*8*2;
    vtkIdType N3 = (N1+N2);
    //
    // We will add 2 new data arrays per point, boxId and occupation(count)
    //
    vtkSmartPointer<vtkIdTypeArray> occupation = vtkSmartPointer<vtkIdTypeArray>::New();
    occupation->SetNumberOfTuples(N3);
    occupation->SetName("Occupation");
    vtkSmartPointer<vtkIdTypeArray> boxId = vtkSmartPointer<vtkIdTypeArray>::New();
    boxId->SetNumberOfTuples(N3);
    boxId->SetName("Partition");
    //
    vtkIdType index = 0;
    double bounds[6];
    vtkBoundingBox box;
    vtkSmartPointer<vtkAppendPolyData> polys = vtkSmartPointer<vtkAppendPolyData>::New();
    //
    // set scalars for each particle
    //
    for (int i=0; i<partitions; i++) {
      vtkIdType numParticles = static_cast<vtkIdType>(data[0+i*13]);
      if ((index+numParticles)>=extent0 && (index<=extent1)) {
        for (vtkIdType p=0; p<numParticles; p++) {
          if (index>=extent0 && index<=extent1) {
            boxId->SetValue(index-extent0, i);
            occupation->SetValue(index-extent0, numParticles);
          }
          index++;
        }
      }
      else {
        index +=numParticles;
        continue;
      }
    }
    //
    // We haven't read all particles, so set index to last valid ID
    //
    index = (extent1-extent0)+1;
    //
    // generate boxes, 2 per partition (6*2+1 = 13 vals per partition)
    // count, min{x,y,z}, max{x,y,z}, ghostmin{x,y,z}, ghostmax{x,y,z}
    //
    for (int i=0; i<partitions; i++) {
      double p1[3],p2[3];
      box.SetMinPoint(&data[1+i*13]);
      box.SetMaxPoint(&data[4+i*13]);
      //
      vtkSmartPointer<vtkOutlineSource> cube1 = vtkSmartPointer<vtkOutlineSource>::New();
      cube1->SetBounds(data[0+1+i*13],data[0+4+i*13],data[0+2+i*13],data[0+5+i*13],data[0+3+i*13],data[0+6+i*13]);
      cube1->Update();
      polys->AddInput(cube1->GetOutput());
      //
      vtkSmartPointer<vtkOutlineSource> cube2 = vtkSmartPointer<vtkOutlineSource>::New();
      cube2->SetBounds(data[6+1+i*13],data[6+4+i*13],data[6+2+i*13],data[6+5+i*13],data[6+3+i*13],data[6+6+i*13]);
      cube2->Update();
      polys->AddInput(cube2->GetOutput());
    }
    polys->Update();
    //
    // Add the generated box points to our original list
    //
    vtkPoints *points = polys->GetOutput()->GetPoints();
    coords->Resize(N3); // might allocate more memory than we need.
    coords->SetNumberOfTuples(N3); // set limits to correct count
    for (vtkIdType P=0; P<N2; P++) {
      coords->SetTuple(N1+P, points->GetPoint(P));
    }
    //
    // Copy lines to output, but add N1 to the point ID of each line vertex
    //
    output->SetLines(polys->GetOutput()->GetLines()); 
    vtkIdType L = output->GetLines()->GetNumberOfCells();
    vtkIdTypeArray *linedata = output->GetLines()->GetData();
    // lines stored as : {2, p1, p2}, {2, p1, p2}, ... 
    for (vtkIdType B=0; B<L; B++) {
      linedata->SetValue(B*3+1, N1 + linedata->GetValue(B*3+1));
      linedata->SetValue(B*3+2, N1 + linedata->GetValue(B*3+2));
    }
    //
    // set scalars for line/box
    //
    for (int i=0; i<partitions; i++) {
      vtkIdType numParticles = static_cast<vtkIdType>(data[0+i*13]);
      for (vtkIdType p=0; p<8*2; p++) {
        boxId->SetValue(index, i);
        occupation->SetValue(index, numParticles);
        index++;
      }
    }
    //
    // The existing scalar array must be padded out for the new points we added
    //
    output->GetPointData()->CopyAllocate(output->GetPointData(),N3);
    for (vtkIdType i=0; i<N2; i++) {
      output->GetPointData()->NullPoint(N1+i);
    }
    //
    // And now add the new arrays (correct size already)
    //
    output->GetPointData()->AddArray(boxId);
    output->GetPointData()->AddArray(occupation);
    //
    // cleanup hdf5
    //
    H5Sclose(diskshape);
  }
  if (dataset_id>0) H5Dclose(dataset_id);
  if (partitiongroup>0) H5Gclose(partitiongroup);
  // put error handling back to previous state.
  H5Eset_auto(H5E_DEFAULT, errfunc, errdata); 

  return N2;
}
//----------------------------------------------------------------------------
int vtkH5PartReader::GetCoordinateArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkH5PartReader::SetCoordinateArrayStatus(const char* name, int status)
{
  if (status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

//----------------------------------------------------------------------------
const char* vtkH5PartReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkH5PartReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkH5PartReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
    {
    if (status)
      {
      this->PointDataArraySelection->EnableArray(name);
      }
    else
      {
      this->PointDataArraySelection->DisableArray(name);
      }
    this->Modified();
    }
}
//----------------------------------------------------------------------------
void vtkH5PartReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkH5PartReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkH5PartReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkH5PartReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkH5PartReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkH5PartReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " <<
    (this->FileName ? this->FileName : "(none)") << "\n";

  os << indent << "NumberOfSteps: " <<  this->NumberOfTimeSteps << "\n";
}
