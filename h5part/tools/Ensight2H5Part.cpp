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
#include "vtkGenericEnSightReader.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkDataArrayCollection.h"
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkStreamOutputWindow.h"
#include "vtkH5PartWriter.h"
#include "vtkTesting.h"

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }

  const char *inputfile = argv[argc-1];
  std::string inputpath, errormsg, casefile, hdf5file;
  if (argc<2 || !inputfile) {
    std::cout << "No input file " << std::endl;
    return 0;
  }
  casefile = inputfile;
  if (!vtksys::SystemTools::FileExists(casefile.c_str()))
  {
    std::cout << "Can't find file " << casefile.c_str() << "\n";
    return 0;
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

  if (test->IsFlagSpecified("-p1")) {
    std::cout << "p1 flag set, using output index 1 for particles" << std::endl;
    DataIndex = 1;
  }

  if (test->IsFlagSpecified("-o")) {
    std::cout << "o flag set, using current dir for output" << std::endl;
    inputpath = vtksys::SystemTools::GetCurrentWorkingDirectory();
  }

  if (test->IsFlagSpecified("-l")) {
    std::cout << "l flag set, setting LittleEndian" << std::endl;
    LittleEndian = true;
  }

//  char *dir = test->GetTempDirectory();

  if (inputpath=="") {
    inputpath = vtksys::SystemTools::GetFilenamePath(casefile.c_str());
  }
  hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(casefile.c_str());
  hdf5file  = inputpath + "/" + hdf5file + ".h5part";
  //
  std::cout << "Checking Ensight file for Time Steps " << casefile.c_str() << std::endl;

  std::vector<double> Values;
  vtkSmartPointer<vtkGenericEnSightReader> reader = vtkSmartPointer<vtkGenericEnSightReader>::New();
  reader->SetCaseFileName(casefile.c_str());
  reader->SetParticleCoordinatesByIndex(1);
  vtkDemandDrivenPipeline* ddp = vtkDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  ddp->UpdateInformation();
  //
  vtkDataArrayCollection  *timesets = reader->GetTimeSets();
  if (timesets) {
    vtkDataArray *timearray = timesets->GetItem(0);
    //
    Values.reserve(timearray->GetNumberOfTuples());
    for (int i=0; i<timearray->GetNumberOfTuples(); i++)
    {
      Values.push_back(timearray->GetTuple1(i));
    }
  }

  int N = Values.size();

  std::cout << "Found " << N << " TimeSteps " << std::endl;

  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileName((char*)hdf5file.c_str());
  writer->SetFileModeToReadWrite();
  for (int i=0; i<N; i++) {
    std::cout << "Converting TimeStep " << i << std::endl;
    if (LittleEndian) {
      reader->SetByteOrderToLittleEndian();
    }
    reader->SetTimeValue(Values[i]);
    reader->Update();
    vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    if (output) {
      vtkSmartPointer<vtkCompositeDataIterator> iter;
      iter.TakeReference(output->NewIterator());
      int p=0;
      for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem(), p++)
      {
        vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
        if (ds && p==DataIndex) {
          writer->SetTimeStep(i);
          writer->SetTimeValue(Values[i]);
          writer->SetInputData(ds);
          writer->Write();
          writer->CloseFile();
        }
      }
    }
  }
  std::cout << "Closing hdf file " << std::endl;
  std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  std::cout << "Done" << std::endl;
}

