/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkASCIIParticleReader.cpp,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkASCIIParticleReader.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkIdList.h"
#include "vtkFloatArray.h"
#include "vtkObjectFactory.h"
#include "vtkAppendPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUtils.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

/*
K is the particle number
T is the current time
X1(k), Y1(k) is the position of particle K
U1(k), V1(k) is the velocity of particle K
P(K) is the pressure of particle K
DENSTAR(K) is the density of particle K
IFLAG(K) is the particle type (1 = fixed wall; 2 = fluid; 3 = moving wall)
VISEFF(K)/VIS0 is the (normalised) viscosity of particle K.

        WRITE(11,12) T
12      FORMAT(E14.8)
           
        DO K=1,N
        write(11,13) K,X1(K),Y1(K),P(K),U1(K),V1(K),DENSTAR(K),&
                     IFLAG(K),VISEFF(K)/VIS0
13      FORMAT(1X,I5,1X,6(E14.8,1X),1X,I1,1X,E14.8)
        ENDDO
*/

//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkASCIIParticleReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkASCIIParticleReader);
//----------------------------------------------------------------------------
vtkASCIIParticleReader::vtkASCIIParticleReader() {
  this->HeaderLines         = 0;
  this->FieldIndices        = NULL;
  this->FieldNames          = NULL;
  this->TimeExpression      = NULL;
  this->IgnoreExpression    = NULL;
  this->InvalidTime         = 0;
  this->DFormat             = 0;
  this->MultiFileScalarMode = 0;
  this->MultiFileCollectNode= 0;

  // Correct for EDF Data (8 fields, 14 numbers)
  this->SetFieldNames("Position,Masse(kg),Velocity(m/s),rho(kg/m3),P(Pa),k(m2/s2),eps(m2/s3),nut(m2/s),S(s-1),kpar,kfluid,kent,priv");
  this->SetFieldIndices("1 2 3,0,4 5 6,7,8,9,10,11,12,13,14,15,16");
  this->SetTimeExpression("T += +([0-9E+-.]+)");
  this->SetIgnoreExpression("TITLE|VARIABLES|ZONE");
  this->SetHeaderLines(0);

  // Correct for Ben Rodgers Data (8 fields, 14 numbers)
  this->SetFieldNames("Position,Velocity,Density,Pressure,Mass,ObsoleteInteger,Vorticity,ObsoleteReal");
  this->SetFieldIndices("0 1 2,3 4 5,6,7,8,9,10 11 12,13");
  this->SetTimeExpression("");
  this->SetIgnoreExpression("");
  this->SetHeaderLines(1);
  this->SetTimeValue(0.0);

  // Correct for David Grahame Data (7 fields, 9 numbers)
  this->SetFieldNames("Position,K,Pressure,Velocity,Density,Type,Viscosity");
  this->SetFieldIndices("1 2 -1,0,3,4 5 -1,6,7,8");
  this->SetTimeExpression("^([0-9E+-.]+)$");
  this->SetIgnoreExpression("");
  this->SetHeaderLines(0);

}
//---------------------------------------------------------------------------
vtkASCIIParticleReader::~vtkASCIIParticleReader() {
  delete []this->FieldIndices;
  delete []this->FieldNames;
  delete []this->TimeExpression;
  delete []this->IgnoreExpression;
}
//----------------------------------------------------------------------------
int vtkASCIIParticleReader::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//---------------------------------------------------------------------------
vtkPolyData *vtkASCIIParticleReader::GetPolyDataOutput()
{
  return vtkPolyData::SafeDownCast(this->GetOutputDataObject(0));
}
//----------------------------------------------------------------------------
void vtkASCIIParticleReader::GetTimeStepValues(std::vector<double> &values)
{
  if (!this->InvalidTime) {
    values = this->TimeStepValues;
  }
}
//----------------------------------------------------------------------------
bool vtkASCIIParticleReader::ReadMetaInformation()
{
  vtkstd::string metaFile = vtkstd::string(this->FileName) + ".txt";
  if (!vtksys::SystemTools::FileExists(metaFile.c_str())) return false;
  //
  vtkstd::ifstream infile(metaFile.c_str());
  int N;
  infile >> N;
  double t;
  long f;
  int r;
  this->TimeStepValues.clear();
  this->FileMarkers.clear();
  this->RowsPerStep.clear();
  if (!infile.good()) return false;
  for (int i=0; i<N; ++i) {
    infile >> t;
    if (!infile.good()) { 
      this->InvalidTime = 1;
      this->RowsPerStep.push_back(0);
      return false;
    }
    infile >> f >> r;
    TimeStepValues.push_back(t);
    FileMarkers.push_back(f);
    RowsPerStep.push_back(r);
  }
  if (!infile.good()) return false;
  return true;
}
//----------------------------------------------------------------------------
bool vtkASCIIParticleReader::WriteMetaInformation()
{
  vtkstd::string metaFile = vtkstd::string(this->FileName) + ".txt";
  vtkstd::ofstream outfile(metaFile.c_str());
  int N = static_cast<int>(TimeStepValues.size());
  outfile << N << std::endl;
  if (!outfile.good()) return false;
  for (int i=0; i<N; ++i) {
    if (!this->InvalidTime) {
      outfile << TimeStepValues[i] << " " << FileMarkers[i] << " " << RowsPerStep[i] << std::endl;
    }
    else {
      outfile << "0" << " " << FileMarkers[i] << " " << RowsPerStep[i] << std::endl;
    }
  }
  outfile << vtkstd::endl;
  return true;
}
//----------------------------------------------------------------------------
int vtkASCIIParticleReader::RequestInformation(vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  if (this->LastFileRead == this->FileName) return 1;
  //
  ifstream infile(this->FileName, ios::binary);
  if (!infile.good()) {
    vtkErrorMacro(<< "Error opening " << this->FileName);
    return 0;
  }
  //
  // Size of file - rows/columns
  //
  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = 0;
  this->InvalidTime = 0;
  //
  infile.seekg(0);
  const int buf_size = 5192;
  char buf[buf_size+1];
  int NumPoints   = 0;
  int linenum     = 0;
  int headercount = this->HeaderLines;
  unsigned int CurrentStep = 0;
  this->RowsPerStep.clear();
  this->FileMarkers.clear();
  this->TimeStepValues.clear();
  // so far zero data points for first time step
  this->RowsPerStep.push_back(0);
  
  vtksys::RegularExpression re1(this->TimeExpression);
  vtksys::RegularExpression re2(this->IgnoreExpression);
  bool useTimeExpression   = strlen(this->TimeExpression)>0;
  bool useIgnoreExpression = strlen(this->IgnoreExpression)>0;
  //
  double timeval;
  //
  if (!this->ReadMetaInformation()) {
    ifstream::pos_type marker = infile.tellg();
    while (infile.getline(buf, buf_size)) {
      vtkUtils::removeTrailing(buf);
      if (strlen(buf)>0 && headercount==0) {
        if (useTimeExpression && re1.find(buf)) {
          // we have found a time stamp, note the time
          timeval = atof(re1.match(1).c_str());
          this->TimeStepValues.push_back(timeval);
          // if time stamp at end of data
          if (NumPoints>0) {
            this->RowsPerStep[CurrentStep] = NumPoints;
            this->RowsPerStep.push_back(0);
            CurrentStep++;
          }
          NumPoints = 0;
          std::cout << "TimeStep : " << CurrentStep << "\tLineNum " << linenum << "\tTIME : " << timeval  << " (" << re1.match(1).c_str() << ")" << std::endl;
        }
        else if (useIgnoreExpression && re2.find(buf)) {
          // skip this line
          std::cout << "SKIPPED : " << "\tLineNum " << linenum << "\t( " << buf << " )" << std::endl;
        }
        else {
          // we have a data line, increment the number of points
          NumPoints++;
          // if we have not marked the start of the timestepdata, do it now
          if (NumPoints==1) {
            this->FileMarkers.push_back(marker);
          }
          // make sure we track how many points are in this step
          this->RowsPerStep[CurrentStep] = NumPoints;
        }
      }
      if (headercount>0) headercount--;
      linenum++;
      marker = infile.tellg();
    }
    if (this->RowsPerStep.size()>0 && this->RowsPerStep.back()==0) {
      this->RowsPerStep.erase(this->RowsPerStep.end()-1);
    }
    if (this->TimeStepValues.size()==this->FileMarkers.size()-1) {
//      std::cout << "FORCED TIMESTEP to  : " << this->TimeValue << std::endl;    
      this->TimeStepValues.push_back(this->TimeValue);
      this->InvalidTime = 1;
    }
  }
  //
  this->LastFileRead = this->FileName;
  if (this->FileMarkers.size()==this->RowsPerStep.size() && 
      this->FileMarkers.size()==this->TimeStepValues.size()) 
  {
    this->TimeStepRange[1] = this->FileMarkers.size()-1;
  }
  else {
    vtkErrorMacro(<< "Error Finding TimeSteps" << this->FileName);
    this->TimeStepRange[1] = this->RowsPerStep.size()-1;
    return 1;
  }

  // Make sure output is ready for data, the number of Columns may
  // be modified by subclass, so allow a new variable to hold the
  // data in our reading code
  if (!this->ParseFields(outputVector)) {
    return 0;
  }
  //
  this->WriteMetaInformation();
  //
  return Superclass::RequestInformation(request, inputVector, outputVector);
}
//----------------------------------------------------------------------------
bool vtkASCIIParticleReader::ProcessOneLine(vtkInformationVector* outputVector, vtkIdType rownum, const char *line)
{
  return 0;
}
//---------------------------------------------------------------------------
// The function object to determine the average
class SumSizes {
  public:
    SumSizes() : sum(0) {}
    // 
    void operator() (const vtkASCIIParticleReader::ScalarPair &elem) {
      sum += elem.second.size();
    }
    operator int() { return (sum); }
    //
    int sum;
};
//----------------------------------------------------------------------------
bool vtkASCIIParticleReader::ParseFields(vtkInformationVector* outputVector) 
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  // get the output
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  // generate list of scalar names
  vtkstd::vector<vtkstd::string> scalarnameslist;
  this->ScalarNamesList.clear();
  vtksys::SystemTools::Split(this->FieldNames, scalarnameslist, ',');
  // generate list of indices for scalars
  vtkstd::vector<vtkstd::string> templist;
  vtksys::SystemTools::Split(this->FieldIndices, templist, ','); // comma delimeter
  // from indices lists, generate integer lists
  for (unsigned int i=0; i<templist.size(); i++) {
    vtkstd::vector<vtkstd::string> intlist;
    vtksys::SystemTools::Split(templist[i].c_str(), intlist, ' '); // space delimiter
    vtkstd::vector<int> IntList;
    for (unsigned int j=0; j<intlist.size(); j++) {
      int index = atoi(intlist[j].c_str());
//      std::cout << index << ":";
      IntList.push_back(index);
    }
    ScalarPair tempPair(scalarnameslist[i], IntList);
    this->ScalarNamesList.push_back(tempPair);
  }
//  std::cout << "\n";

  int totalSize = for_each(
    this->ScalarNamesList.begin(), 
    this->ScalarNamesList.end(),
    SumSizes()
  );

  this->NumColumns = totalSize;
  return true;
}
//---------------------------------------------------------------------------
bool vtkASCIIParticleReader::AllocateOutputData(vtkInformationVector* outputVector, int numrows, int numcols, int &actualcolumns) 
{
  vtkPolyData *output = vtkPolyData::SafeDownCast(this->GetOutput());
  this->Points  = vtkPoints::New();
  this->Cells   = vtkCellArray::New();
  //
  this->Points->SetNumberOfPoints(numrows);
  this->Cells->Allocate(numrows);
  //
  this->Scalars = new vtkFloatArray*[this->ScalarNamesList.size()];
  for (unsigned int i=0; i<this->ScalarNamesList.size(); i++) {
    this->Scalars[i] = vtkFloatArray::New();
    this->Scalars[i]->SetNumberOfComponents(this->ScalarNamesList[i].second.size());
    this->Scalars[i]->SetNumberOfTuples(numrows);
    this->Scalars[i]->SetName(this->ScalarNamesList[i].first.c_str());
    // the zero-th scalar array is the point coordinates.
    if (i==0) {
      this->Points->SetData(this->Scalars[0]);
    }
    else if (this->MultiFileCollectNode || !this->MultiFileScalarMode) { 
      output->GetPointData()->AddArray(this->Scalars[i]);
    }
    this->Scalars[i]->Delete();
  }
  actualcolumns = this->NumColumns;
  //
  output->SetPoints(this->Points);
  output->SetVerts(this->Cells);
  this->Points->Delete();
  this->Cells->Delete();

  return 1;
}
//---------------------------------------------------------------------------
void vtkASCIIParticleReader::CleanUp(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  if (this->ScalarNamesList.size()>0) {
    for (unsigned int i=0; i<this->ScalarNamesList.size(); i++) {
//      this->Scalars[i]->Delete();
    }
    delete []Scalars;
  }
}
//----------------------------------------------------------------------------
int vtkASCIIParticleReader::RequestData(vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())) {
    this->ActualTimeStep = outInfo->Get(vtkDataObject::DATA_TIME_STEPS())[0];
  }
  //
  const int buf_size = 65536;
  char    buf[buf_size];
  vtkstd::vector<double> onerow;
  //
  // File OK
  //
  ifstream infile(this->FileName);
  infile.seekg(this->FileMarkers[this->ActualTimeStep]);

  //
  // Read all the data and put into appropriate dataset
  //
  onerow.assign(this->NumColumns, 0);
  int rownum = 0;
  int numrows = RowsPerStep[this->ActualTimeStep];
  int actualcolumns = this->NumColumns;
  int numcols = this->NumColumns;
  //
  this->AllocateOutputData(outputVector, numrows, numcols, actualcolumns);
  //
  while (infile.getline(buf, buf_size) && rownum<numrows) {
    if (strlen(buf)>0) {
      if (!ProcessOneLine(outputVector, rownum, buf)) {
        std::string tstr(buf);
        std::replace(tstr.begin(),tstr.end(),this->Delimiter,' ');
        if (this->DFormat) std::replace(tstr.begin(),tstr.end(),'D','E');
        std::istringstream tempstr(tstr.c_str());
        for (int i=0; i<this->NumColumns; i++) onerow[i] = 0.0;
        int col = 0;
        while (col<this->NumColumns && tempstr >> onerow[col++]) { }
        this->ProcessOneRow(outputVector, rownum++, &onerow[0]);
      }
    }
  }
  //
  // Allow subclasses to clean up memory
  //
  if (!this->MultiFileScalarMode) {
    this->CleanUp(NULL, NULL, NULL);
  }
  return 1;
}
//---------------------------------------------------------------------------
void vtkASCIIParticleReader::ProcessOneRow(vtkInformationVector* outputVector, vtkIdType rownum, double *onerow) {
  float point[3] = {0.0, 0.0, 0.0};
  int NumScalars = this->ScalarNamesList.size();
  double temptuple[3];  
  if (this->MultiFileScalarMode) {
    for (int i=0; i<3; i++) temptuple[i] = 0.0;
    temptuple[0] = onerow[0];
    this->Scalars[0]->SetTuple(rownum, temptuple);
  }
  else {
    // clear
    for (int i=0; i<3; i++) temptuple[i] = 0.0;
    //
    for (int i=0; i<NumScalars; i++)
    {
      ScalarPair &temp = this->ScalarNamesList[i];
      vtkstd::vector<int> &indices = temp.second;
      for (unsigned int s=0; s<indices.size(); s++) {
        if (indices[s]>-1) temptuple[s] = onerow[indices[s]];
      }
      this->Scalars[i]->SetTuple(rownum, temptuple);
    }
  }
  Cells->InsertNextCell(1,&rownum);
}
//---------------------------------------------------------------------------
void vtkASCIIParticleReader::CollectMultiFileScalars(ScalarList &scalars) 
{
  if (this->NumColumns!=scalars.size()) {
    vtkErrorMacro(<<"Mismatch in number of columns/scalars");
    return;
  }

  float point[3] = {0.0, 0.0, 0.0};
  int NumScalars = this->ScalarNamesList.size();
  double temptuple[3];  
  // clear
  for (int i=0; i<3; i++) temptuple[i] = 0.0;
  //
  vtkSmartPointer<vtkDataArray> coords;
  coords.TakeReference(this->Scalars[0]->NewInstance());
  coords->DeepCopy(this->Scalars[0]);
  //
  for (int t=0; t<scalars[0]->GetNumberOfTuples(); t++) {
    for (int i=0; i<NumScalars; i++)
    {
      ScalarPair &temp = this->ScalarNamesList[i];
      vtkstd::vector<int> &indices = temp.second;
      for (unsigned int s=0; s<indices.size(); s++) {
        if (indices[s]>-1) temptuple[s] = scalars[indices[s]]->GetTuple(t)[0];
//        if (indices[s]>-1) temptuple[s] = scalars[s]->GetTuple(t)[0];
      }

      // we must not overwrite the array we are contributing
      if (i==0) {
        coords->SetTuple(t, temptuple);
      } 
      else {
        this->Scalars[i]->SetTuple(t, temptuple);
      }
    }
  }
  this->Scalars[0]->DeepCopy(coords);
  this->CleanUp(NULL, NULL, NULL);
}
//---------------------------------------------------------------------------
