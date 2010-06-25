// DGraham Data
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/DAT_t_6_2.ascii 
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/Graham_wave_data_New.ascii 
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/DATA_2.ascii 
//
// EDF Data
// Ascii2H5Part -c C:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.edf.cfg -f C:/data/ParticleData/EDF/test_case_output95.dat
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.edf.cfg -f H:\ParticleData\EDF\test_case\test_case_output1.dat
//
// UMan Data
// Ascii2H5Part -c C:/cmakebuild/plugins/bin/ASCII2H5Part.manchester.cfg  -f C:/Data/ParticleData/UMan/Case6/PART_0001 -a C:/Data/ParticleData/UMan/Case6/PART_0000_point_attribs.txt   
//
// ECN Data
// Ascii2H5Part -c C:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.ecn.cfg -f C:/Data/ParticleData/ECN/BGL_SPhere0POSX.dat
// Ascii2H5Part -c D:/Code/CSCS/csviz/Shared/vtkCSCS/vtkH5Part/Tools/ASCII2H5Part.ecn.cfg -f D:/data/ecn/MillenISOP0POSX.dat 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib> 
#include <vtksys/SystemTools.hxx>
#include <vtksys/Glob.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>
#include <sstream>

#include "vtkSmartPointer.h"
#include "vtkTesting.h"
#include "Testing/Cxx/vtkTestUtilities.h"
#include "vtkAppendPolyData.h"
#include "vtkPointData.h"
#include "vtkASCIIParticleReader.h"
#include "vtkH5PartWriter.h"
#include "vtkFloatArray.h"
//
#include "FileSeriesFinder.h"
//
typedef vtkstd::vector<vtkstd::string> stringlist;
typedef vtkstd::pair< vtkstd::string, vtkstd::string > maptype;
//
vtkstd::map<vtkstd::string, vtkstd::string> RegExMap;
stringlist patterns, regexmatches, Tstrings, Bstrings, Vstrings;
vtkstd::string regex;
//
int Tindex = -1, Tnum = 0;
int Bindex = -1, Bnum = 0;
int Vindex = -1, Vnum = 0;
vtkstd::string timeform, blockform, varform, filepattern;

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkstd::cout << "Usage : ASCII2H5Part "  
    << "-c full/path/to/.cfg "
    << "-a full/path/to/IdFile for optional point Id flagging "
    << "-t full/path/to/timefile.txt " 
    << "-f full/path/to/input/file " << vtkstd::endl;
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }

  std::string asciipath, asciiname, errormsg, asciifile, hdf5file, ascii2h5config, ascii2h5IdFile;

  // config file
  char *cfg_name = vtkTestUtilities::GetArgOrEnvOrDefault("-c", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::GetFilenamePath(cfg_name).size()==0) {
    std::string cfg_path = vtksys::SystemTools::GetCurrentWorkingDirectory();
    ascii2h5config = cfg_path + "/" + cfg_name;
  }
  else {
    ascii2h5config = cfg_name;
  }
  delete []cfg_name;
  //
  if (!vtksys::SystemTools::FileExists(ascii2h5config.c_str()) ||
       vtksys::SystemTools::FileIsDirectory(ascii2h5config.c_str())) 
  {
    std::cout << "Can't find config file " << ascii2h5config.c_str() << "\n";
    return 0;
  }
  vtkstd::cout << "Config file found    : "  << ascii2h5config.c_str() << "\n";;

  // input file
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault("-f", argc, argv, "DUMMY_ENV_VAR", "");
  asciipath = vtksys::SystemTools::GetFilenamePath(filename);
  asciiname = vtksys::SystemTools::GetFilenameName(filename);
  asciifile = asciipath + "/" + asciiname;
  delete []filename;
  if (!vtksys::SystemTools::FileExists(asciifile.c_str()))
  {
    std::cout << "Can't find input file " << asciifile.c_str() << "\n";
    return 0;
  }
  vtkstd::cout << "Input file found     : "  << asciifile.c_str() << "\n";;

  // Point Id file
  char *idf_name = vtkTestUtilities::GetArgOrEnvOrDefault("-a", argc, argv, "DUMMY_ENV_VAR", "");
  vtkstd::map<long int, double> idMap;
  if (vtksys::SystemTools::FileExists(idf_name)) {
    ascii2h5IdFile = idf_name;
    vtkstd::cout << "Using Id file        : "  << ascii2h5IdFile.c_str() << "\n";
    vtkstd::ifstream idfile(idf_name);
    long int low, high;
    double flag;
    while (idfile.good()) {
      idfile >> low >> high >> flag;
      idMap.insert( vtkstd::pair<long int,double>(high, flag));
    }
  }
  delete []idf_name;

  // time override 
  bool   OverrideTime = false;
  double OverrideTimeStep = 0.0;
  char *time_steps = vtkTestUtilities::GetArgOrEnvOrDefault("-t", argc, argv, "DUMMY_ENV_VAR", "");
  std::string timeoverride = time_steps;
  if (timeoverride.size()>0) {
    OverrideTime = true;
    OverrideTimeStep = atof(timeoverride.c_str());
    vtkstd::cout << "Time override set to : " << OverrideTimeStep << " per step \n";;
  }
  delete []time_steps;

  //
  // generate new h5part file name
  //
  hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(asciifile);
  hdf5file  = asciipath + "/" + hdf5file + ".h5part";
  vtkstd::cout << "Output HDF5 filename : "  << hdf5file.c_str() << "\n";;
  //
  // Now scan the file names to build up our pattern
  //
  vtkstd::string FieldNames;
  vtkstd::string FieldIndices;
  vtkstd::string TimeExpr;
  vtkstd::string IgnoreExp;
  //
  vtkstd::string FileNamePattern;
  vtkstd::string PrefixRegEx, BlockRegEx, VarRegEx, ExtRegEx, TimeRegEx, Text0RegEx;
  //
  int headerlines=0;
  int DFormat=0;
  // open config file
  vtkstd::ifstream cfg(ascii2h5config.c_str());
  char buf[512];
  while (cfg.good()) {
    cfg.getline(buf,512);
    if (strstr(buf, "FieldNames" )!=0)  { cfg.getline(buf,512); FieldNames = buf; }
    if (strstr(buf, "FieldIndices")!=0) { cfg.getline(buf,512); FieldIndices = buf; }
    if (strstr(buf, "TimeExpr")!=0)     { cfg.getline(buf,512); TimeExpr = buf; }
    if (strstr(buf, "IgnoreExp")!=0)    { cfg.getline(buf,512); IgnoreExp = buf; }
    if (strstr(buf, "HeaderLines")!=0)  { cfg.getline(buf,512); sscanf(buf,"%i", &headerlines);}
    if (strstr(buf, "DFormat")!=0)      { cfg.getline(buf,512); sscanf(buf,"%i", &DFormat); }    
    //
    if (strstr(buf, "FileNamePattern")!=0)  { cfg.getline(buf,512); FileNamePattern = buf; }
    if (strstr(buf, "PREFIX_RegEx")!=0)     { cfg.getline(buf,512); PrefixRegEx = buf; }
    if (strstr(buf, "BLOCK_RegEx")!=0)      { cfg.getline(buf,512); BlockRegEx = buf; }
    if (strstr(buf, "VAR_RegEx")!=0)        { cfg.getline(buf,512); VarRegEx = buf; }
    if (strstr(buf, "EXT_RegEx")!=0)        { cfg.getline(buf,512); ExtRegEx = buf; }
    if (strstr(buf, "TIME_RegEx")!=0)       { cfg.getline(buf,512); TimeRegEx = buf; }
    if (strstr(buf, "TEXT0_RegEx")!=0)      { cfg.getline(buf,512); Text0RegEx = buf; }
  }
  //
  FileSeriesFinder *Finder = new FileSeriesFinder(FileNamePattern);
  Finder->SetPrefixRegEx(PrefixRegEx.c_str());
  Finder->SetTimeRegEx(TimeRegEx.c_str());
  Finder->SetBlockRegEx(BlockRegEx.c_str());
  Finder->SetVarRegEx(VarRegEx.c_str());
  Finder->SetExtRegEx(ExtRegEx.c_str());
  Finder->SetText0RegEx(Text0RegEx.c_str());
  Finder->Scan(asciifile.c_str());
//  Finder->GetTimeValues(TimeStepValues);
//  NumberOfTimeSteps = Finder->GetNumberOfTimeSteps();

  bool multiblock = (Finder->GetNumberOfBlocks()>0);
  bool multivar   = (Finder->GetNumberOfVars()>0);
  bool multitime  = (Finder->GetNumberOfTimeSteps()>0);

  std::vector<double> time_values;
  typedef vtkstd::vector< vtkSmartPointer<vtkASCIIParticleReader> > varReader;
  typedef vtkstd::vector< varReader > blockReader;
  blockReader readers;
  int t=0;
  int b=0;
  do {
    int v=0;
    varReader varreaders;
    do {
      vtkSmartPointer<vtkASCIIParticleReader> reader = vtkSmartPointer<vtkASCIIParticleReader>::New();
      varreaders.push_back(reader);     
      vtkstd::string filename = Finder->GenerateFileName(t, b, v);
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
          Tnum = time_values.size();
          std::cout << "Found " << Tnum << " Time steps in file " << std::endl;
        }
        else {
          for (int i=0; i<Tnum; i++) {
            time_values.push_back(i);
          }
        }
      }
    } while (++v<Vnum);
    //
    readers.push_back(varreaders);
  } while (++b<Bnum);
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
      int v=0;
      ScalarList scalars;
      std::cout << "Reading Block " << b << std::endl;

      do { // Var
        std::cout << "Reading Var " << v << std::endl;
        reader = readers[b][v];
        vtkstd::string filename = Finder->GenerateFileName(t,b,v);
	      reader->SetFileName(filename.c_str());
        if (!multitime) {
          reader->SetTimeStep(t);
        }
        reader->Update();
        vtkstd::vector<double> values;
	      reader->GetTimeStepValues(values);
        if (values.size()>0) time_values[t] = values[t];

        if (multivar) {
          vtkPolyData *p = vtkPolyData::SafeDownCast(reader->GetOutput());
          vtkFloatArray *floats = vtkFloatArray::SafeDownCast(p->GetPoints()->GetData());
          scalars.push_back(floats);
        }
        //
      } while (++v<Vnum);

      if (multivar) {
        readers[b][0]->CollectMultiFileScalars(scalars);
      }
      polys = vtkPolyData::SafeDownCast(readers[b][0]->GetOutput());
      //
      if (multiblock && polys) append->AddInput(polys);
    } while (++b<Bnum);
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
        for (vtkstd::map<long int, double>::iterator it=idMap.begin(); it!=idMap.end(); ++it) {
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
    if (time_values.size()>t) {
      writer->SetTimeValue(time_values[t]);
    }
    else {
      writer->SetTimeValue(t);
    }
    if (OverrideTime) {
      double val = OverrideTimeStep*(double)t;
      writer->SetTimeValue(val);
    }
    writer->Write();

    // to make sure reference counting deletes stuff
    reader = NULL;
    polys = NULL;
  } while (++t<Tnum);

  std::cout << "Closing hdf file " << std::endl;
  std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  std::cout << "Done" << std::endl;
}

