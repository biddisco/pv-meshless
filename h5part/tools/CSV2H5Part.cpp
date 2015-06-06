// VTM2H5Part -a -f H:/mhdparticles/parallel-1.00000.vtm -o D:/Data -t /project/s168/simonsch/s15h_64_complete/timetable.dat

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib> 
#include <sstream>
#ifdef WIN32
  #include <Windows.h>
#else
  #include <unistd.h>
#endif

#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>

#include "vtkSmartPointer.h"
#include "vtkAlgorithmOutput.h"

#include "vtkDelimitedTextReader.h"
#include "vtkTableToPolyData.h"
#include "vtkTable.h"
#include "vtkMaskPoints.h"

#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkAppendPolyData.h"
#include "vtkH5PartReader.h"
#include "vtkH5PartWriter.h"
//#include "vtkStreamOutputWindow.h"
#include "vtkTesting.h"
#include "vtkTestUtilities.h"

#include "FileSeriesFinder.h"

//----------------------------------------------------------------------------
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
//----------------------------------------------------------------------------
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}
//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  std::cout << "Usage : CSV2H5Part " << std::endl << "\t"
    << "[-a append to file] " << std::endl << "\t"
    << "[-w wait for new data] " << std::endl << "\t"
    << "[-n just scan input, don't write] " << std::endl << "\t"
    << "[-o specify an output directory] " << std::endl << "\t"
    << "[-t ascii time file name] " << std::endl << "\t"
    << "-f file.csv" << std::endl;

  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }


  // input file
  char *inputfile = vtkTestUtilities::GetArgOrEnvOrDefault("-f", argc, argv, "DUMMY_ENV_VAR", "");
  if (!vtksys::SystemTools::FileExists(inputfile))
  {
    std::cout << "Can't find file " << inputfile << std::endl;
    return 0;
  }
  std::string inputpath = vtksys::SystemTools::GetFilenamePath(inputfile);
  std::string inputname = vtksys::SystemTools::GetFilenameName(inputfile);
  std::string ascii2h5config;
  std::string outputdir = inputpath;
  // output dir
  char *outdir = vtkTestUtilities::GetArgOrEnvOrDefault("-o", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::FileIsDirectory(outdir))
  {
    outputdir = outdir;
    std::cout << "Using output directory " << outdir << std::endl;
  }
  std::string hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(inputfile);
  hdf5file  = outputdir + "/" + hdf5file + ".h5part";
  std::cout << "Output file name is : " << hdf5file.c_str() << std::endl;
  // time file
  bool useTimeFile = false;
  char *timefile = vtkTestUtilities::GetArgOrEnvOrDefault("-t", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::FileExists(timefile))
  {
    std::cout << "Using time from " << timefile << std::endl;
    useTimeFile = true;
  }

  //
  // Is no write flag present
  //
  bool no_write = false;
  if (test->IsFlagSpecified("-n")) {
    std::cout << "No write (-n flag) set" << std::endl;
    no_write = true;
  }

  //
  // Is append flag set for adding to existing data
  //
  bool append = false;
  if (test->IsFlagSpecified("-a")) {
    std::cout << "Append new time steps mode set " << std::endl;
    append = true;
  }
  //
  // Is append flag set for adding to existing data
  //
  bool wait_for_simulation = false;
  if (test->IsFlagSpecified("-w")) {
    std::cout << "Wait mode activited " << std::endl;
    wait_for_simulation = true;
  }

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
  std::cout << "Config file found    : "  << ascii2h5config.c_str() << "\n";;



  //
  // Force the creation of our output window object
  //
  //vtkStreamOutputWindow *outwin = vtkStreamOutputWindow::New();
  //vtkOutputWindow::SetInstance(outwin);
  //outwin->Delete();
  //outwin->SetOutputStream(&std::cout);

  //
  // Now scan the file names to build up our pattern
  //
  std::string FieldNames;
  std::string FieldIndices;
  std::string TimeExpr;
  std::string IgnoreExp;
  //
  std::string FileNamePattern;
  std::string PrefixRegEx, BlockRegEx, VarRegEx, ExtRegEx, TimeRegEx, Text0RegEx;
  //
  int headerlines=0;
  int DFormat=0;
  // open config file
  std::ifstream cfg(ascii2h5config.c_str());
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

  std::cout << "Using field names \n";
  std::vector<std::string> FieldNamesList = split(FieldNames,',');
  std::copy(FieldNamesList.begin(), FieldNamesList.end(), std::ostream_iterator<std::string>(std::cout, ","));
  std::cout << "\n";

  std::cout << "Using indices \n";
  std::vector<std::string> FieldIndicesList = split(FieldIndices,',');
  std::copy(FieldIndicesList.begin(), FieldIndicesList.end(), std::ostream_iterator<std::string>(std::cout, ","));
  std::cout << "\n";

  std::vector< std::vector<std::string> > indices_lists;
  for (std::vector<std::string>::iterator it=FieldIndicesList.begin(); it!=FieldIndicesList.end(); ++it) {
      indices_lists.push_back(split(trim(*it),' '));
      std::cout << "index size " << indices_lists.back().size() << " : ";
      std::copy(indices_lists.back().begin(), indices_lists.back().end(), std::ostream_iterator<std::string>(std::cout, ","));
      std::cout << "\n";
  }
  //
  FileSeriesFinder *Finder = new FileSeriesFinder(FileNamePattern);
  Finder->SetPrefixRegEx(PrefixRegEx.c_str());
  Finder->SetTimeRegEx(TimeRegEx.c_str());
  Finder->SetBlockRegEx(BlockRegEx.c_str());
  Finder->SetVarRegEx(VarRegEx.c_str());
  Finder->SetExtRegEx(ExtRegEx.c_str());
  Finder->SetText0RegEx(Text0RegEx.c_str());
  Finder->Scan(inputfile);

  bool multiblock = (Finder->GetNumberOfBlocks()>0);
  bool multivar   = (Finder->GetNumberOfVars()>0);
  bool multitime  = (Finder->GetNumberOfTimeSteps()>0);

  std::vector<double> time_values;
  int Tnum = 0;
  //
  if (multitime) {
    Tnum = Finder->GetNumberOfTimeSteps();
    Finder->GetTimeValues(time_values);
  }

  //
  // Stop here if not writing results
  //
  if (no_write) return 0;

  std::cout << "Found " << Finder->GetNumberOfTimeSteps() << " Time steps" << std::endl;

  do {

// #define JB_SPECIAL_TIME_OVERRIDE 1
#if defined(JB_SPECIAL_TIME_OVERRIDE) 
    {
      for (int i=0; i<time_values.size(); i++) {
        double tval = 0.213202+ 0.0001*i;
	      std::cout << "overwriting " << i << " with " << tval << std::endl;
        values[i] = tval;
      }
    }
#endif

    if (useTimeFile) {
      std::ifstream timedata(timefile);
      unsigned int tstart=0;
      while (timedata.good()) {
        char timebuff[512];
        timedata.getline(timebuff,512);
        std::stringstream temp;
        temp << timebuff;
        double tval;
        temp >> tval;
        if (tstart>=0 && tstart<time_values.size()) {
            time_values[tstart] = tval;
  	      std::cout << "overwriting " << tstart << " with " << tval << std::endl;
        }
	      tstart++;
      }
    }
    //
    // scan existing file to see which time step to resume from
    //
    int StartStep = 0;
    if (append && vtksys::SystemTools::FileExists((char*)hdf5file.c_str())) {
      vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
      reader->SetFileName((char*)hdf5file.c_str());
      while (reader->HasStep(StartStep)) StartStep++;
      std::cout << "Appending data from Step number " << StartStep << std::endl;
    } 

    //
    // Create writer
    //
    vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
    writer->SetFileName((char*)hdf5file.c_str());
    writer->SetFileModeToReadWrite();

    //
    // Loop over input files and read/write
    //
    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetDetectNumericColumns(1);
    reader->SetTrimWhitespacePriorToNumericConversion(1);
    reader->SetHaveHeaders(0);
    reader->SetFieldDelimiterCharacters(" ");
    reader->SetMergeConsecutiveDelimiters(1);

    std::vector<double>::iterator it2 = time_values.begin();
    for (int t=0; t<Finder->GetNumberOfTimeSteps(); ++t, ++it2)
    {
      std::string name = Finder->GenerateFileName(t);
      int Step = static_cast<int>(it2 - time_values.begin());
      if (append && Step<StartStep) continue;
      //
      std::cout << "Converting " << name.c_str() << std::endl;
      reader->SetFileName(name.c_str());
      reader->Update();

      vtkTable *table = reader->GetOutput();
      for (int i=0; i<table->GetNumberOfColumns(); i++) {
          vtkAbstractArray *arr = table->GetColumn(i);

          std::string name = std::string("col") + std::to_string(i);
          for (int f=0; f<FieldIndicesList.size(); ++f) {
              if (indices_lists[f].size()==1) {
                  int index = std::stoi(indices_lists[f][0]);
                  if (index == i) {
                      name = FieldNamesList[f];
                  }
              }
          }

          arr->SetName(name.c_str());
          std::cout << "Set column name " << name.c_str() << std::endl;
      }
      //
      // convert table data to polydata
      //
      int numParticles = 0;
      if (1) {
        vtkSmartPointer<vtkTableToPolyData> poly = vtkSmartPointer<vtkTableToPolyData>::New();
        vtkSmartPointer<vtkMaskPoints> mask = vtkSmartPointer<vtkMaskPoints>::New();
        poly->SetInputConnection(reader->GetOutputPort());

        std::vector<std::string> list3 = split(FieldIndicesList[0].c_str(),' ');
        std::copy(list3.begin(), list3.end(), std::ostream_iterator<std::string>(std::cout, ","));
        std::cout << "\n";

        poly->SetXColumn( std::string(std::string("col") + list3[0]).c_str() );
        poly->SetYColumn( std::string(std::string("col") + list3[1]).c_str() );
        poly->SetZColumn( std::string(std::string("col") + list3[2]).c_str() );

        mask->SetInputConnection(poly->GetOutputPort());
        mask->SetSingleVertexPerCell(1);
        mask->SetGenerateVertices(1);
        mask->SetOnRatio(1);
        mask->SetMaximumNumberOfPoints(VTK_ID_MAX);
        mask->Update();
        //
        writer->SetInputConnection(mask->GetOutputPort());
        numParticles = mask->GetOutput()->GetNumberOfCells();
      }
      else {    
      }
      //
      // write to hdf5
      //
      writer->SetTimeStep(Step);
      writer->SetTimeValue(*it2);
      std::cout << "Writing step " << Step << " \tParticles : " << numParticles << "\ttime = " << *it2 << std::endl;
      writer->Write();
      writer->CloseFile();
      std::cout << "Writing step " << Step << " done" << std::endl;
    }
    //
    // Clean up
    //
    std::cout << "Closing hdf file " << std::endl;
    std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
    writer->CloseFile();
    std::cout << "Done" << std::endl;
    //
    // wait if told to restart
    //
    if (wait_for_simulation) {
      std::cout << "Waiting 2 minutes before restart" << std::endl;
#ifdef WIN32
      Sleep(1000*120);
#else
      sleep(120);
#endif
    }
  } while (wait_for_simulation);
  //
  // really done now
  //
}


