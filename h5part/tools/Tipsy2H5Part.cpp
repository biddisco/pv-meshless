// Ensight2H5Part -p1 H:\ParticleData\granulair\RibbonBlender\RibbonBlender.case
// Ensight2H5Part "C:\Data\vatech\RESULTS\fluidSPH.case"
// Ensight2H5Part H:\ParticleData\vevey\deflector\FLUID.case

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib> 
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>

#include "vtkSmartPointer.h"
#include "vtkTipsyReader.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkDataArrayCollection.h"
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkStreamOutputWindow.h"
#include "vtkH5PartWriter.h"
#include "vtkTesting.h"
#include "Testing/Cxx/vtkTestUtilities.h"

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkstd::cout << "Usage : VTM2H5Part " << vtkstd::endl << "\t"
    << "[-a append to file] " << vtkstd::endl << "\t"
    << "[-w wait for new data] " << vtkstd::endl << "\t"
    << "[-n just scan input, don't write] " << vtkstd::endl << "\t"
    << "[-o specify an output directory] " << vtkstd::endl << "\t"
    << "[-t ascii time file name] " << vtkstd::endl << "\t"
    << "-f file.vtm|vtp" << vtkstd::endl;

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
  vtkstd::string inputpath = vtksys::SystemTools::GetFilenamePath(inputfile);
  vtkstd::string inputname = vtksys::SystemTools::GetFilenameName(inputfile);
  vtkstd::string outputdir = inputpath;
  // output dir
  char *outdir = vtkTestUtilities::GetArgOrEnvOrDefault("-o", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::FileIsDirectory(outdir))
  {
    outputdir = outdir;
    std::cout << "Using output directory " << outdir << std::endl;
  }
  vtkstd::string hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(inputfile);
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
  vtkStreamOutputWindow *outwin = vtkStreamOutputWindow::New();
  vtkOutputWindow::SetInstance(outwin);
  outwin->Delete();
  outwin->SetOutputStream(&std::cout);

  int DataIndex = 0;
  int LittleEndian = 0;

  //
  std::cout << "Checking input file for Time Steps " << inputfile << std::endl;

  std::vector<double> timesteps;
  vtkSmartPointer<vtkTipsyReader> reader = vtkSmartPointer<vtkTipsyReader>::New();
  reader->SetFileName(vtksys::SystemTools::ConvertToUnixOutputPath(inputfile).c_str());
  vtkDemandDrivenPipeline* ddp = vtkDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  ddp->UpdateInformation();
  //
  reader->GetTimeStepValues(timesteps);
  int N = timesteps.size();

  std::cout << "Found " << N << " TimeSteps " << std::endl;

  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileName((char*)hdf5file.c_str());
  writer->SetFileModeToReadWrite();
  for (int i=0; i<N; i++) {
    std::cout << "Converting TimeStep " << i << std::endl;
    reader->SetTimeStep(i);
    reader->Update();
    vtkPointSet *output = reader->GetOutput();
    if (output) {
      writer->SetTimeStep(i);
      writer->SetTimeValue(timesteps[i]);
      writer->SetInput(output);
      writer->Write();
      writer->CloseFile();
    }
  }
  std::cout << "Closing hdf file " << std::endl;
  std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  std::cout << "Done" << std::endl;
}

