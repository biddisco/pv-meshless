/*=========================================================================

  Project                 : XdmfGenerator
  Module                  : FileSeriesFinder.cxx

  Authors:
     John Biddiscombe     Jerome Soumagne
     biddisco@cscs.ch     soumagne@cscs.ch

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This work has received funding from the European Community's Seventh
  Framework Programme (FP7/2007-2013) under grant agreement 225967 “NextMuSE”

=========================================================================*/
#include "FileSeriesFinder.h"

#define VERBOSE_FINDER

//----------------------------------------------------------------------------
FileSeriesFinder::FileSeriesFinder(std::string filenamepattern)
{
    FileNamePattern  = filenamepattern;
    PrefixRegEx      = "(.*)[^0-9]";
    TimeRegEx        = "([0-9]+)";
    BlockRegEx       = "([0-9]+)";
    BlockSubDirRegEx = "";
    VarRegEx         = "(.+)";
    Text0RegEx       = "(.*)";
    ExtRegEx         = "(\\..*)";
}
//----------------------------------------------------------------------------
FileSeriesFinder::FileSeriesFinder()
{
    FileNamePattern  = "PREFIX TIME EXT";
    PrefixRegEx      = "(.*)[^0-9]";
    TimeRegEx        = "([0-9]+)";
    BlockRegEx       = "([0-9]+)";
    BlockSubDirRegEx = "";
    VarRegEx         = "(.+)";
    Text0RegEx       = "(.*)";
    ExtRegEx         = "(\\..*)";
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetFileName(const char *filename)
{
  this->TimeIndex = 0;
  this->Tstrings.clear();
  this->Tstrings.push_back(filename);
  this->NumberOfTimeSteps = 1;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetPrefixRegEx(const char *prefix_regex)
{
  if (prefix_regex) PrefixRegEx = prefix_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetTimeRegEx(const char *time_regex)
{
  if (time_regex) TimeRegEx = time_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetBlockRegEx(const char *block_regex)
{
  if (block_regex) BlockRegEx = block_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetBlockSubDirRegEx(const char *block_regex)
{
  if (block_regex) BlockSubDirRegEx = block_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetVarRegEx(const char *var_regex)
{
  if (var_regex) VarRegEx = var_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetText0RegEx(const char *text0_regex)
{
  if (text0_regex) Text0RegEx = text0_regex;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::SetExtRegEx(const char *ext_regex)
{
  if (ext_regex) ExtRegEx = ext_regex;
}
//----------------------------------------------------------------------------
std::string FileSeriesFinder::GenerateNumericPattern(std::string match)
{
  std::string result;
  std::stringstream t1;
  size_t len = match.length();
  if (len>1 && match[0]=='0') {
    t1 << "%0" << len << "i";
    result = t1.str();
  }
  else {
    result = "%i";
  }
  return result;
}
//----------------------------------------------------------------------------
std::string FileSeriesFinder::GenerateGlobString(unsigned int index, stringlist &patterns, bool regexmode)
{
  std::string globpattern;
  for (unsigned int i=0; i<patterns.size(); i++) {
    if (i==index) {
      if (regexmode) {
        globpattern = globpattern + "(.*)";
      }
      else {
        globpattern = globpattern + "*";
      }
    } 
    else {
      globpattern = globpattern + patterns[i];
    }
  }
  return globpattern;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::FindGlobbedSegments(std::string &pattern, stringlist &files, bool numeric)
{
  vtksys::RegularExpression re(pattern.c_str());
  for (stringlist::iterator it=files.begin(); it!=files.end(); ++it) {
    if (re.find(*it)) (*it) = re.match(1);
  }
  if (numeric) {
    typedef std::pair<int,std::string> pairType;
    typedef std::map< int, std::string > indexMap;
    std::map< int, std::string > newlist;
    for (stringlist::iterator it=files.begin(); it!=files.end(); ++it) {
      int num;
      std::stringstream temp;
      temp << it->c_str();
      temp >> num;
      if (temp.eof()) {
        newlist.insert( pairType(atoi(it->c_str()),*it) );
      }
    }
    files.clear();
    // read all strings back from map - iterator traverses in sorted order
    // so all int-like strings come out as we want them
    for (indexMap::iterator it=newlist.begin(); it!=newlist.end(); ++it) {
      files.push_back(it->second);
    }
  }
  else {
    std::sort(files.begin(), files.end());
  }
}
//----------------------------------------------------------------------------
std::string FileSeriesFinder::GenerateFileName(int T)
{
  return this->GenerateFileName(T,0,0);
}
//----------------------------------------------------------------------------
std::string FileSeriesFinder::GenerateFileName(int T, int B, int V)
{
  std::string result;
  for (unsigned int i=0; i<regexmatches.size(); i++) {
    if      ((int)i==TimeIndex) result += Tstrings[T];
    else if ((int)i==BlockIndex) result += Bstrings[B];
    else if ((int)i==VarIndex) result += Vstrings[V];
    // when prefix is generated, we must fudge the processor path :(
    else if (i==0 && BlockSubDirRegEx.size()>0) {
      result += regexmatches[i];
      vtksys::SystemTools::ReplaceString(result,BlockSubDirRegEx.c_str(), 
        std::string("proc" + Bstrings[B]).c_str());
    }
    else result += regexmatches[i];
  }
  if (regexmatches.size()==0 && this->Tstrings.size()==1) {
    return this->Tstrings[0];
  }
  return result;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::Scan(const char *inputfile)
{
  //
  // initialize
  //
  NumberOfTimeSteps = 0;
  NumberOfBlocks    = 0;
  NumberOfVars      = 0;
  TimeIndex     = -1;
  BlockIndex    = -1;
  VarIndex      = -1; 
  regex.clear();
  //
  //
  //
  std::string cleaninput = inputfile;
  vtksys::SystemTools::ConvertToUnixSlashes(cleaninput);
  //
  // find which parts of search pattern are present
  //
  vtksys::SystemTools::Split(FileNamePattern.c_str(), patterns, ' ');
  if (PrefixRegEx.size()>0)      RegExMap.insert(maptype("PREFIX", PrefixRegEx));
  if (BlockRegEx.size()>0)       RegExMap.insert(maptype("BLOCK",  BlockRegEx));
  if (BlockSubDirRegEx.size()>0) RegExMap.insert(maptype("SUBDIR", BlockSubDirRegEx));
  if (VarRegEx.size()>0)         RegExMap.insert(maptype("VAR",    VarRegEx));
  if (ExtRegEx.size()>0)         RegExMap.insert(maptype("EXT",    ExtRegEx));
  if (TimeRegEx.size()>0)        RegExMap.insert(maptype("TIME",   TimeRegEx));
  if (Text0RegEx.size()>0)       RegExMap.insert(maptype("TEXT0",  Text0RegEx));
  for (stringlist::iterator it=patterns.begin(); it!=patterns.end(); ++it) {
    regex += RegExMap[*it];    
  }
  //
#ifdef VERBOSE_FINDER
  std::cout << "File name search using Regular expression : " << regex.c_str() << "\n";;
#endif
  //
  // Get the order of pattern parts right
  //
  stringlist::iterator pos;
  if ((pos=std::find(patterns.begin(), patterns.end(), "TIME")) != patterns.end()) {
    TimeIndex = pos-patterns.begin();
  }
  if ((pos = std::find(patterns.begin(), patterns.end(), "BLOCK")) != patterns.end()) {
    BlockIndex = pos-patterns.begin();
  }
  if ((pos = std::find(patterns.begin(), patterns.end(), "VAR")) != patterns.end()) {
    VarIndex = pos-patterns.begin();
  }
  if ((pos = std::find(patterns.begin(), patterns.end(), "VAR")) != patterns.end()) {
    VarIndex = pos-patterns.begin();
  }

  //
  // Detect if any numeric parts need zero padding 
  //
  if (regex.size()==0) return;
  vtksys::RegularExpression re1(regex.c_str());
  if (re1.find(cleaninput)) {
#ifdef VERBOSE_FINDER
    for (unsigned int i=0; i<patterns.size(); i++) {
      std::cout << patterns[i] << "=" << re1.match(i+1).c_str() << "  " << "\n";;
    }
#endif
    //
    // Are block and time numbers padded with zeroes?
    //    
    if (TimeIndex!=-1) { 
      timeform = GenerateNumericPattern(re1.match(TimeIndex+1));
    }
    if (BlockIndex!=-1) {
      blockform = GenerateNumericPattern(re1.match(BlockIndex+1));
    }
    if (VarIndex!=-1) {
      varform = "%s";
    }
    //
    // put it all together      
    for (unsigned int i=0; i<patterns.size(); i++) {
      regexmatches.push_back(re1.match(i+1));
      if ((int)i==TimeIndex) {
        filepattern = filepattern + timeform;
      } 
      else if ((int)i==BlockIndex) {
        filepattern = filepattern + blockform;
      }
      else if ((int)i==VarIndex) {
        filepattern = filepattern + varform;
      }
      else {
        filepattern = filepattern + re1.match(i+1);
      }
    }
#ifdef VERBOSE_FINDER
    std::cout << "Complete Pattern = " << filepattern.c_str() << "\n";;
#endif
  }
  //
  // Now scan and see how many files there are matching the patterns
  //
  vtksys::Glob glob;
  // Time 
  if (TimeIndex!=-1) {
    std::string globpattern = GenerateGlobString(TimeIndex, regexmatches, false);
    vtksys::SystemTools::ConvertToUnixSlashes(globpattern);
    if (glob.FindFiles(globpattern)) {
      Tstrings = glob.GetFiles();
      std::string regexpattern = GenerateGlobString(TimeIndex, regexmatches, true);
      FindGlobbedSegments(regexpattern, Tstrings, true);
      NumberOfTimeSteps = (int)Tstrings.size();
      //
#ifdef VERBOSE_FINDER
      std::cout << "#######################################################" << "\n";
      std::cout << "Found " << NumberOfTimeSteps << " Time files \n";
      std::cout << "#######################################################" << "\n";
      for (stringlist::iterator it=Tstrings.begin(); it!=Tstrings.end(); ++it) {
//          std::cout << (*it).c_str() << "\n";
      }
#endif
    }
  }
  // Block
  if (BlockIndex!=-1) {
    std::string globpattern = GenerateGlobString(BlockIndex, regexmatches, false);
    if (BlockSubDirRegEx!="") {
      vtksys::SystemTools::ReplaceString(globpattern,BlockSubDirRegEx.c_str(), "*");
    }
    if (glob.FindFiles(globpattern)) {
      Bstrings = glob.GetFiles();
      NumberOfBlocks = (int)Bstrings.size();
      std::string regexpattern = GenerateGlobString(BlockIndex, regexmatches, true);
      if (BlockSubDirRegEx!="") {
        vtksys::SystemTools::ReplaceString(regexpattern,BlockSubDirRegEx.c_str(), ".*");
      }
      FindGlobbedSegments(regexpattern, Bstrings, true);
#ifdef VERBOSE_FINDER
      std::cout << "#######################################################" << "\n";
      std::cout << "Found " << NumberOfBlocks << " Block files \n";
      std::cout << "#######################################################" << "\n";
      for (stringlist::iterator it=Bstrings.begin(); it!=Bstrings.end(); ++it) {
//          std::cout << (*it).c_str() << "\n";
      }
#endif
    }
  }
  // Var
  if (VarIndex!=-1) {
    std::string globpattern = GenerateGlobString(VarIndex, regexmatches, false);
    if (glob.FindFiles(globpattern)) {
      Vstrings = glob.GetFiles();
      NumberOfVars = (int)Vstrings.size();
      std::string regexpattern = GenerateGlobString(VarIndex, regexmatches, true);
      FindGlobbedSegments(regexpattern, Vstrings, false);
#ifdef VERBOSE_FINDER
      std::cout << "#######################################################" << "\n";
      std::cout << "Found " << NumberOfVars << " Var files \n";
      std::cout << "#######################################################" << "\n";
      for (stringlist::iterator it=Vstrings.begin(); it!=Vstrings.end(); ++it) {
//          std::cout << (*it).c_str() << "\n";
      }
#endif
    }
  }
}
//----------------------------------------------------------------------------
void FileSeriesFinder::TestFilenameGeneration()
{
  //
  //
  //
  std::cout << "#######################################################" << "\n";
  std::cout << "Testing Filename Generation" << "\n";
  std::cout << "#######################################################" << "\n";
  int t=0, b=0, v=0;
  do {
    b=0;
    do {
      v=0;
      do {
        std::cout << this->GenerateFileName(t, b, v).c_str() << std::endl;
      } while (++v<NumberOfVars);
    } while (++b<NumberOfBlocks);
  } while (++t<NumberOfTimeSteps);
}
//----------------------------------------------------------------------------
int FileSeriesFinder::GetNumberOfTimeSteps()
{
  return NumberOfTimeSteps;
}
//----------------------------------------------------------------------------
int FileSeriesFinder::GetNumberOfBlocks()
{
  return NumberOfBlocks;
}
//----------------------------------------------------------------------------
int FileSeriesFinder::GetNumberOfVars()
{
  return NumberOfVars;
}
//----------------------------------------------------------------------------
double FileSeriesFinder::GetTimeValue(int index)
{
  if (index>=0 && index<this->NumberOfTimeSteps) {
    std::stringstream t;
    t << Tstrings[index];
    double v;
    t >> v;
    return v;
  }
  return 0.0;
}
//----------------------------------------------------------------------------
void FileSeriesFinder::GetTimeValues(std::vector<double> &values)
{
  values.clear();
  for (int i=0; i<this->NumberOfTimeSteps; i++) {
    std::stringstream t;
    t << Tstrings[i];
    double v;
    t >> v;
    values.push_back(v);
  }
}
//----------------------------------------------------------------------------
int main_test(int argc, char **argv)
{

  FileSeriesFinder finder("PREFIX TIME EXT");
  finder.Scan("D:/data/csf/CSF16transient-16000.dat");
  finder.TestFilenameGeneration();
  //
  int T = finder.GetNumberOfTimeSteps();
  std::cout << T << " " << finder.GetNumberOfBlocks() << " " << finder.GetNumberOfVars() << " ";
  for (int t=0; t<T; t++) {
    std::cout << finder.GenerateFileName(t,0,0).c_str() << std::endl;
  }
  

//  NumberOfBlocks = 1;

/*
  bool multiblock = (NumberOfBlocks>0);
  bool multivar   = (NumberOfVars>0);
  bool multitime  = (NumberOfTimeSteps>0);

  std::vector<double> time_values;
  typedef std::vector< vtkSmartPointer<vtkASCIIParticleReader> > varReader;
  typedef std::vector< varReader > blockReader;
  blockReader readers;
  t=0;
  b=0;
  do {
    v=0;
    varReader varreaders;
    do {
      vtkSmartPointer<vtkASCIIParticleReader> reader = vtkSmartPointer<vtkASCIIParticleReader>::New();
      varreaders.push_back(reader);     
      std::string filename = GenerateFileName(t, b, v);
      reader->SetFileName(filename.c_str());
      reader->SetFieldNames(FieldNames.c_str());
      reader->SetFieldIndices(FieldIndices.c_str());
      reader->SetTimeExpression(TimeExpr.c_str());
      reader->SetIgnoreExpression(IgnoreExp.c_str());
      reader->SetHeaderLines(headerlines);
      reader->SetDFormat(DFormat);
      reader->SetTimeStep(0);
      reader->SetMultiFileScalarMode(multivar);
      reader->SetMultiFileCollectNode(multivar && v==0);
      if (b==0 && v==0) {
        if (!multitime) {
          std::cout << "Scanning for Time Steps " << asciifile.c_str() << std::endl;
          reader->UpdateInformation();
          reader->GetTimeStepValues(time_values);
          NumberOfTimeSteps = time_values.size();
          std::cout << "Found " << NumberOfTimeSteps << " Time steps in file " << std::endl;
        }
        else {
          for (int i=0; i<NumberOfTimeSteps; i++) {
            time_values.push_back(i);
          }
        }
      }
    } while (++v<NumberOfVars);
    //
    readers.push_back(varreaders);
  } while (++b<NumberOfBlocks);
  //

  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileName((char*)hdf5file.c_str());
  vtkSmartPointer<vtkASCIIParticleReader> reader;
  vtkSmartPointer<vtkPolyData> polys;

  t=0;
  do { // Time
    b=0;
    std::cout << "Reading TimeStep " << t << std::endl;
    //
    vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
    do { // Block    
      v=0;
      ScalarList scalars;
      std::cout << "Reading Block " << b << std::endl;

      do { // Var
        std::cout << "Reading Var " << v << std::endl;
        reader = readers[b][v];
        std::string filename = GenerateFileName(t,b,v);
	      reader->SetFileName(filename.c_str());
        if (!multitime) {
          reader->SetTimeStep(t);
        }
        reader->Update();
        std::vector<double> values;
	      reader->GetTimeStepValues(values);
        if (values.size()>0) time_values[t] = values[0];

        if (multivar) {
          vtkPolyData *p = vtkPolyData::SafeDownCast(reader->GetOutput());
          vtkFloatArray *floats = vtkFloatArray::SafeDownCast(p->GetPoints()->GetData());
          scalars.push_back(floats);
        }
        //
      } while (++v<NumberOfVars);

      if (multivar) {
        readers[b][0]->CollectMultiFileScalars(scalars);
      }
      polys = vtkPolyData::SafeDownCast(readers[b][0]->GetOutput());
      //
      if (multiblock && polys) append->AddInput(polys);
    } while (++b<NumberOfBlocks);
    //
    vtkPolyData *realData = polys;
    if (multiblock) {
      append->Update();
      realData = append->GetOutput();
    }
    if (idMap.size()) {
      vtkSmartPointer<vtkFloatArray> IdScalars = vtkSmartPointer<vtkFloatArray>::New(); 
      IdScalars->SetName("PointIdFlags");
      double scalar = 0.0;
      for (int i=0; i<realData->GetNumberOfPoints(); i++) {
        for (std::map<long int, double>::iterator it=idMap.begin(); it!=idMap.end(); ++it) {
          if (i<=it->first) {
            scalar = it->second;
            IdScalars->InsertNextTuple1(scalar);
            break;
          }
        }
      }
      if (IdScalars->GetNumberOfTuples()==realData->GetNumberOfPoints()) {
        realData->GetPointData()->AddArray(IdScalars);
      }
    }
    writer->SetInput(realData);
    writer->SetTimeStep(t);
    writer->SetTimeValue(time_values[t]);
    if (OverrideTime) {
      double val = OverrideTimeStep*(double)t;
      writer->SetTimeValue(val);
    }
    writer->Write();

    // to make sure reference counting deletes stuff
    reader = NULL;
    polys = NULL;
  } while (++t<NumberOfTimeSteps);

  std::cout << "Closing hdf file " << std::endl;
  std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  std::cout << "Done" << std::endl;
*/

  return 0;
}

