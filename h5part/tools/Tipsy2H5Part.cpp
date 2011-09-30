// cd C:\cmakebuild\pv-meshless\bin\Debug
// "C:\Program Files\MPICH2\bin\mpiexec" -localonly -n 2 c:\cmakebuild\plugins\bin\Debug\Tipsy2H5Part.exe -o C:\data\corbett -f C:\data\corbett\b1.00300.d0-1000.std

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib> 
#include "vtkToolkits.h"                // For VTK_USE_MPI 
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>

#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDataArrayCollection.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkStreamOutputWindow.h"
#include "vtkTimerLog.h"
#include "vtkTesting.h"
#include "Testing/Cxx/vtkTestUtilities.h"

#ifdef VTK_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif

#include "vtkTipsyReader.h"
#include "vtkH5PartWriter.h"
#include "vtkParticlePartitionFilter.h"
#include "vtkSPHProbeFilter.h"
#include "vtkSPHManager.h"

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // 
  // Initialize and setup MPI
  //
  MPI_Init(&argc, &argv);
  vtkMPIController* controller = vtkMPIController::New();
  controller->Initialize(&argc, &argv, 1);
  int rank = controller->GetLocalProcessId();
  int numProcs = controller->GetNumberOfProcesses();
  vtkMultiProcessController::SetGlobalController(controller);
  vtkAlgorithm::SetDefaultExecutivePrototype(vtkStreamingDemandDrivenPipeline::New());
  //
  //
  //
  if (rank==0) {
    vtkstd::cout << "Usage : VTM2H5Part " << vtkstd::endl << "\t"
      << "[-a append to file] " << vtkstd::endl << "\t"
      << "[-w wait for new data] " << vtkstd::endl << "\t"
      << "[-n just scan input, don't write] " << vtkstd::endl << "\t"
      << "[-o specify an output directory] " << vtkstd::endl << "\t"
      << "[-t ascii time file name] " << vtkstd::endl << "\t"
      << "-f file.vtm|vtp" << vtkstd::endl;
    //
    char c;
    std::cin >> c;
  }
  MPI_Barrier(*vtkMPICommunicator::SafeDownCast(controller->GetCommunicator())->GetMPIComm()->GetHandle());

  //
  // Command line args
  //
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }

  //
  // input file
  //
  char *inputfile = vtkTestUtilities::GetArgOrEnvOrDefault("-f", argc, argv, "DUMMY_ENV_VAR", "");
  if (!vtksys::SystemTools::FileExists(inputfile))
  {
    std::cout << "Can't find file " << inputfile << std::endl;
    return 0;
  }
  vtkstd::string inputpath = vtksys::SystemTools::GetFilenamePath(inputfile);
  vtkstd::string inputname = vtksys::SystemTools::GetFilenameName(inputfile);
  vtkstd::string outputdir = inputpath;

  //
  // output dir
  //
  char *outdir = vtkTestUtilities::GetArgOrEnvOrDefault("-o", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::FileIsDirectory(outdir))
  {
    outputdir = outdir;
    if (rank==0) std::cout << "Using output directory " << outdir << std::endl;
  }
  vtkstd::string hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(inputfile);
  hdf5file  = outputdir + "/" + hdf5file + ".h5part";
  if (rank==0)  std::cout << "Output file name is : " << hdf5file.c_str() << std::endl;

  //
  // time file. Not used here
  //
  bool useTimeFile = false;
  char *timefile = vtkTestUtilities::GetArgOrEnvOrDefault("-t", argc, argv, "DUMMY_ENV_VAR", "");
  if (vtksys::SystemTools::FileExists(timefile))
  {
    if (rank==0) std::cout << "Using time from " << timefile << std::endl;
    useTimeFile = true;
  }

  //
  // Is no write flag present
  //
  bool no_write = false;
  if (test->IsFlagSpecified("-n")) {
    if (rank==0) std::cout << "No write (-n flag) set" << std::endl;
    no_write = true;
  }

  //
  // Is append flag set for adding to existing data
  // not used here
  //
  bool append = false;
  if (test->IsFlagSpecified("-a")) {
    if (rank==0) std::cout << "Append new time steps mode set " << std::endl;
    append = true;
  }
  //
  // Is append flag set for adding to existing data
  // not used here
  //
  bool wait_for_simulation = false;
  if (test->IsFlagSpecified("-w")) {
    if (rank==0) std::cout << "Wait mode activited " << std::endl;
    wait_for_simulation = true;
  }

  //
  // Force the creation of our output window object
  //
  vtkStreamOutputWindow *outwin = vtkStreamOutputWindow::New();
  vtkOutputWindow::SetInstance(outwin);
  outwin->Delete();
  outwin->SetOutputStream(&std::cout);

  //
  // update reader information to parse files and get time
  //
  if (rank==0) std::cout << "Checking input file for Time Steps " << inputfile << std::endl;
  std::vector<double> timesteps;
  vtkSmartPointer<vtkTipsyReader> reader = vtkSmartPointer<vtkTipsyReader>::New();
  reader->SetFileName(vtksys::SystemTools::ConvertToUnixOutputPath(inputfile).c_str());
  reader->DebugOff();
  vtkDemandDrivenPipeline* ddp = vtkDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  vtkSmartPointer<vtkInformation> execInfo = ddp->GetOutputInformation(0);
  execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numProcs);
  execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), rank);
  if (rank==0) std::cout << "Setting piece " << rank << " of " << numProcs << std::endl;
  ddp->UpdateInformation();

  //
  // get timesteps from reader
  //
  reader->GetTimeStepValues(timesteps);
  unsigned int N = static_cast<unsigned int>(timesteps.size());
  if (rank==0) std::cout << "Found " << N << " TimeSteps " << std::endl;

  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetController(controller);
  partitioner->SetGhostCellOverlap(0.01);
  partitioner->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkSPHManager> manager = vtkSmartPointer<vtkSPHManager>::New();
  manager->SetInterpolationMethodToKernel();  
  manager->SetLimitSearchByNeighbourCount(1);  
  manager->SetMaximumNeighbours(16);
  manager->SetMaximumRadius(0.01);
  manager->SetKernelTypeToCubicSpline();
  manager->SetKernelDimension(3);
  manager->SetDefaultParticleSideLength(0.01);
  manager->SetDefaultDensity(1);
  manager->SetHCoefficient(1.5); 

  vtkSmartPointer<vtkSPHProbeFilter> sph = vtkSmartPointer<vtkSPHProbeFilter>::New();
  sph->SetInputConnection(partitioner->GetOutputPort());
  sph->SetProbeConnection(partitioner->GetOutputPort());
  sph->SetComputeDensityFromNeighbourVolume(1);
  sph->SetSPHManager(manager);
  sph->SetMassScalars("mass");

  //
  // Setup writer
  //
  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileName((char*)hdf5file.c_str());
  writer->SetFileModeToWrite();
  writer->SetController(controller);
  for (unsigned int i=0; i<N; i++) {
    vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
    timer->StartTimer();
    //
    if (rank==0) std::cout << "Converting TimeStep " << i << std::endl;
    reader->SetTimeStep(i);
    //
    vtkDataSet *output = sph->GetOutput();
    if (output) {
      std::cout << "Setting piece information " << rank << " of " << numProcs << std::endl;
      vtkStreamingDemandDrivenPipeline *sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(writer->GetExecutive());
      vtkSmartPointer<vtkInformation> execInfo = sddp->GetOutputInformation(0);    
      //
      writer->SetTimeStep(i);
      writer->SetTimeValue(timesteps[i]);
      writer->SetInputConnection(sph->GetOutputPort());
      //
      // To get parallel writing correct, we need to make sure that piece
      // information is passed upstream. first update information,
      // then set piece update extent, 
      //
      if (!no_write) {
        // no piece info set yet, assumes infor is not piece dependent
        sddp->UpdateInformation();
        sddp->SetUpdateExtent(0, rank, numProcs, 0);
        // this calls writer->Write(); with correct piece
        sddp->Update();
      }
      writer->CloseFile();
    }
    timer->StopTimer();
    controller->Barrier();
    std::cout << "Process Id : " << rank << " time " << timer->GetElapsedTime() << " s" << std::endl;
  }
  if (rank==0) std::cout << "Closing hdf file " << std::endl;
  if (rank==0) std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  //
  MPI_Finalize();
  //
  std::cout << "Done" << std::endl;
}

