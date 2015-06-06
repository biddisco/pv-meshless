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
#include "vtkXMLMultiBlockDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkAppendPolyData.h"
#include "vtkH5PartReader.h"
#include "vtkH5PartWriter.h"
//#include "vtkStreamOutputWindow.h"
#include "vtkTesting.h"
#include "vtkTestUtilities.h"

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  std::cout << "Usage : VTM2H5Part " << std::endl << "\t"
    << "[-a append to file] " << std::endl << "\t"
    << "[-w wait for new data] " << std::endl << "\t"
    << "[-n just scan input, don't write] " << std::endl << "\t"
    << "[-o specify an output directory] " << std::endl << "\t"
    << "[-t ascii time file name] " << std::endl << "\t"
    << "-f file.vtm|vtp" << std::endl;

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

  //
  // Force the creation of our output window object
  //
  //vtkStreamOutputWindow *outwin = vtkStreamOutputWindow::New();
  //vtkOutputWindow::SetInstance(outwin);
  //outwin->Delete();
  //outwin->SetOutputStream(&std::cout);

  //
  // From input file name, guess pattern for all time steps
  //
  std::cout << "Scanning directory for Time Steps " << inputpath.c_str() << std::endl;
  bool MultiBlock = false;
  char fpattern[256];
  vtksys::RegularExpression re("(.*[^0-9]+)([0-9]+)\\.(.*)");
  if (re.find(inputfile)) {
    // are we padding with zeroes
    int  len = re.match(2).length();
    bool pad = false;
    if (len>1 && re.match(2)[0]=='0') pad = true;
    // create format string for number
    std::stringstream t1;
    if (pad) t1 << "%0" << len << "i";
    else t1 << "%" << len << "i";
    sprintf(fpattern, "%s%s.%s", re.match(1).c_str(), t1.str().c_str(), re.match(3).c_str());  
    if (re.match(3)=="vtm" || re.match(3)=="vtmb") {
      MultiBlock = true;
      std::cout << "File extension indicates Multi-Block data" << std::endl;
    }
    else if (re.match(3)=="vtp") {
      // ok
    }
    else {
      std::cout << "Expected vtm or vtp, assuming vtp type data" << std::endl;
    }
  }
  else { 
    std::cout << "Pattern match failed" << std::endl;
    return 0;
  }
  //
  // Stop here if not writing results
  //
  if (no_write) return 0;

  do {
    //
    // How many time steps actually exist
    //
    char fname[256];
    int i, c=0;
    std::vector<double> values;
    std::vector<std::string> names;
    for (i=0; i<9999; i++) {
      sprintf(fname, fpattern, i);
      if (vtksys::SystemTools::FileExists(fname)) {
        values.push_back(i);
        names.push_back(fname);
        c = 0;
      }
      else {
        std::cout << "No data for " << i << " : " << fname << std::endl;
        c++;
      }
      // only allow one missing time step, if more (consecutive), stop
      if (c>2) {
        std::cout << "Consecutive steps missing, aborting scan " << i << std::endl;
        break;
      }
    }
    std::cout << "Found " << values.size() << " Time steps" << std::endl;

// #define JB_SPECIAL_TIME_OVERRIDE 1
#if defined(JB_SPECIAL_TIME_OVERRIDE) 
    {
      for (int i=0; i<values.size(); i++) {
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
        if (tstart>=0 && tstart<values.size()) {
          values[tstart] = tval;
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
    vtkSmartPointer<vtkXMLReader> reader;
    std::vector<double>::iterator it2 = values.begin();
    for (std::vector<std::string>::iterator it=names.begin(); it!=names.end(); ++it, ++it2)
    {
      int Step = static_cast<int>(it2 - values.begin());
      if (append && Step<StartStep) continue;
      //
      if (MultiBlock) {
        reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
      } else {
        reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
      }
      reader->SetFileName((*it).c_str());
      reader->Update();
      //
      // convert multiblock data to flat dataset
      //
      int numParticles = 0;
      if (MultiBlock) {
        vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
        vtkSmartPointer<vtkXMLMultiBlockDataReader> Xreader = vtkXMLMultiBlockDataReader::SafeDownCast(reader);
        vtkCompositeDataSet *dobj = vtkCompositeDataSet::SafeDownCast(Xreader->GetOutput());
        if (!dobj) continue;

        vtkSmartPointer<vtkCompositeDataIterator> iter;
        iter.TakeReference(dobj->NewIterator());
        for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
        {
          vtkPolyData *polys = vtkPolyData::SafeDownCast(iter->GetCurrentDataObject());
          if (polys) append->AddInputData(polys);
        }
        append->Update();
        //
        writer->SetInputConnection(append->GetOutputPort());
        numParticles = append->GetOutput()->GetNumberOfCells();
      }
      else {    
        writer->SetInputConnection(reader->GetOutputPort());
        numParticles = reader->GetOutputAsDataSet()->GetNumberOfCells();
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


