// mpiexec : C:\Program Files\MPICH2\bin\mpiexec
// mpiargs : -localonly -n 2 -env PATH D:\cmakebuild\pv-meshless\bin\debug;c:\bin
// appcmd  : $(TargetPath)
// appargs : "-D" "D:/Code/plugins/pv-meshless/Testing/data" "-F" "jet.h5part" "-T" "D:/cmakebuild/plugins/pv-meshless/Testing/Temporary" "-particlesize" "0.001" "-ghost_region" "0.005" "-massScalars" "mass" "-scalar" "ShepardCoeff" "-contour" "0.5" "-densityScalars" "FLUID_density" "-imagetest" "1" "-V" "D:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResampleSerial.png"
// appargs : -D D:/Code/plugins/pv-meshless/Testing/data -F dam-17.h5part -T D:/cmakebuild/plugins/pv-meshless/Testing/Temporary -particlesize 0.005 -ghost_region 0.01 -gridSpacing 0.0025 0.0025 0.0025 -scalar ShepardCoeff -contour 0.1 -densityScalars FLUID_rho -imagetest 1 -V D:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResample-large-parallel-4.png
// mpishim : C:\Program Files\Microsoft Visual Studio 9.0\Common7\IDE\Remote Debugger\x64\mpishim.exe

// //breno
// set path=c:\bin;d:\cmakebuild\pv-meshless\bin\release;d:\cmakebuild\cmake\bin\Debug
// "C:\Program Files\MPICH2\bin\mpiexec.exe" "-localonly" "-n" "4" "C:/cmakebuild/plugins/bin/Debug/TestSPHProbeFilter.exe" "-D" "C:/Code/plugins/pv-meshless/Testing/data" "-F" "dam-17.h5part" "-T" "C:/cmakebuild/plugins/pv-meshless/Testing/Temporary" "-particlesize" "0.005" "-ghost_region" "0.01" "-gridSpacing" "0.0025 0.0025 0.0025" "-scalar" "ShepardCoeff" "-contour" "0.1" "-densityScalars" "FLUID_rho" "-imagetest" "1" "-V" "C:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResample-large-parallel-4.png"
// "C:\Program Files\MPICH2\bin\mpiexec.exe" "-localonly" "-n" "4" "C:/cmakebuild/plugins/bin/Debug/TestSPHProbeFilter.exe" "-D" "C:/Code/plugins/pv-meshless/Testing/data" "-F" "dam-17.h5part" "-T" "C:/cmakebuild/plugins/pv-meshless/Testing/Temporary" "-particlesize" "0.005" "-ghost_region" "0.01" "-gridSpacing" "0.0025 0.0025 0.0025" "-scalar" "ShepardCoeff" "-contour" "0.1" "-densityScalars" "FLUID_rho" "-imagetest" "1" "-V" "C:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResample-large-parallel-4.png"

// Breno Extent debugging
//
// "C:\Program Files\MPICH2\bin\mpiexec.exe" "-localonly" "-n" "8" "D:/cmakebuild/plugins/bin/Debug/TestSPHProbeFilter.exe" 
// -D D:/Code/plugins/pv-meshless/Testing/data -F dam-17.h5part -T D:/cmakebuild/plugins/pv-meshless/Testing/Temporary -particlesize 0.005 -ghost_region 0.01 -gridSpacing 0.01 0.01 0.01 -scalar ShepardCoeff -contour 0.1 -densityScalars FLUID_rho -imagetest 1 -V D:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResample-large-parallel-8-0.005.png -I
#ifdef _WIN32
  #include <windows.h>
#else 
  #include <sys/time.h>
#endif

#include "Testing/Cxx/vtkTestUtilities.h"
#include "Testing/Cxx/vtkRegressionTestImage.h"
//
#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkPointSource.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkMPIController.h"
#include "vtkParallelFactory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkDebugLeaks.h"
#include "vtkElevationFilter.h"
#include "vtkH5PartReader.h"
#include "vtkMaskPoints.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkTimerLog.h"
#include "vtkProcessIdScalars.h"
#include "vtkContourFilter.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataNormals.h"
//
#include "vtkSPHProbeFilter.h"
#include "vtkSPHImageResampler.h"
#include "vtkSPHManager.h"
#include "vtkParticleIdFilter.h"
#include "vtkStreamOutputWindow.h"
//
#ifdef PV_MESHLESS_TRILINOS
  #include "vtkParticlePartitionFilter.h"
#endif
//
#include <vtksys/SystemTools.hxx>
#include <sstream>
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#define _USE_MATH_DEFINES
#include <cmath>

#if 0
  #define OUTPUTTEXT(a) std::cout << (a);

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << myRank << " : " a << endl; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }
#else
  #undef vtkWarningMacro
  #undef vtkDebugMacro
  #define vtkDebugMacro(a) 
#endif

//----------------------------------------------------------------------------
static const int ISO_OUTPUT_TAG=301;
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
template <typename T>
void DisplayParameter(const char *prefix1, const char *prefix2, T *data, int components, int rank)
{
  vtkstd::stringstream temp;
  temp << prefix1 << prefix2 << std::ends;
  std::cout.width(30);
  std::cout << temp.str().c_str() << " : {";
  std::cout.width(0);
  for (int i=0; i<components; i++) {
    std::cout << data[i];
    (i==(components-1)) ? std::cout << "}" : std::cout << ",";
  }
  std::cout << std::endl;
}
//----------------------------------------------------------------------------
template <typename T>
T GetParameter(const char *argstr, const char *message, int argc, char **argv, T defaultvalue, int rank, bool &valueset)
{
  char *tempChar = vtkTestUtilities::GetArgOrEnvOrDefault(argstr, argc, argv, "", "");
  T newValue = defaultvalue;
  valueset = false;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> newValue;
    if (rank==0) {
      DisplayParameter<T>(message, "", &newValue, 1, rank);
    }
    valueset = true;
  }
  delete []tempChar;
  return newValue;
}
//----------------------------------------------------------------------------
template <typename T>
bool GetArrayParameter(const char *argstr, const char *message, T *data, int components, int argc, char **argv, int rank)
{
  char *tempChar = vtkTestUtilities::GetArgOrEnvOrDefault(argstr, argc, argv, "", "");
  bool valueset = false;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    for (int i=0; i<components; i++) temp >> data[i];
    if (rank==0) {
      std::cout.width(30);
      std::cout << message << " : {";
      std::cout.width(0);
      for (int i=0; i<components; i++) {
        std::cout << data[i];
        (i==(components-1)) ? std::cout << "}" : std::cout << ",";
      }
      std::cout << std::endl;
    }
    valueset = true;
  }
  delete []tempChar;
  return valueset;
}
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  int retVal = 1;

  // This is here to avoid false leak messages from vtkDebugLeaks when
  // using mpich. It appears that the root process which spawns all the
  // main processes waits in MPI_Init() and calls exit() when
  // the others are done, causing apparent memory leaks for any objects
  // created before MPI_Init().
  MPI_Init(&argc, &argv);
  vtkSmartPointer<vtkMPIController> controller = vtkSmartPointer<vtkMPIController>::New();
  controller->Initialize(&argc, &argv, 1);
  
  // Obtain the id of the running process and the total
  // number of processes
  vtkTypeInt64 myRank = controller->GetLocalProcessId();
  vtkTypeInt64 numProcs = controller->GetNumberOfProcesses();
  if (myRank==0) {
//    char ch;
//    std::cout << "Attach debugger" << std::endl;
//    std::cin >> ch;
  }
  controller->Barrier();

  //--------------------------------------------------------------
  // command line params : Setup testing utilities/args etc
  //--------------------------------------------------------------
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }
  // Get test filename etc
  char *empty = "";
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-F", argc, argv, "DUMMY_ENV_VAR", "temp.h5");
  char* fullname = vtkTestUtilities::ExpandDataFileName(argc, argv, filename);
  if (myRank==0) {
    DisplayParameter<char *>("====================", "", &empty, 1, myRank);
    DisplayParameter<char *>("FileName", "", &fullname, 1, myRank);
  }

  //
  // Force the creation of our output window object
  //
  vtkSmartPointer<vtkStreamOutputWindow> outwin = vtkSmartPointer<vtkStreamOutputWindow>::New();
  vtkOutputWindow::SetInstance(outwin);
  outwin->SetOutputStream(&std::cout);

  //--------------------------------------------------------------
  // Testing params
  //--------------------------------------------------------------
  bool unused, fixNeighbours, fixRadius, cameraSet;
  double gridSpacing[3] = {0.0, 0.0, 0.0};
  double vminmax[2] = {0.0, 0.0};
  double vpos[3] = {0.0, 0.0, 0.0};
  double cameraPosition[3] = {0.0, 0.0, 0.0};
  double cameraFocus[3] = {0.0, 1.0, 0.0};
  double cameraViewUp[3] = {0.0, 0.0, 1.0};
  int windowSize[2] = {400,400};

  //
  // General test info
  //
  std::string testName = GetParameter<std::string>("-testName", "Test name", argc, argv, "", myRank, unused);
  //
  // SPH kernel or neighbour info
  //
  double particleSize = GetParameter<double>("-particlesize", "Particle Size", argc, argv, 0, myRank, fixRadius);
  double        ghost = GetParameter<double>("-ghost_region", "Ghost Region", argc, argv, 0.0, myRank, unused);
  unused = GetArrayParameter<double>("-gridSpacing", "Grid Spacing", gridSpacing, 3, argc, argv, myRank);
  int            maxN = GetParameter<int>("-neighbours", "Fixed Neighbours", argc, argv, 0, myRank, fixNeighbours);
  std::string    massScalars = GetParameter<std::string>("-massScalars", "Mass Scalar Array", argc, argv, "", myRank, unused);
  std::string densityScalars = GetParameter<std::string>("-densityScalars", "Density Scalar Array", argc, argv, "", myRank, unused);
  //
  // Test/Display of results
  //
  std::string     scalarname = GetParameter<std::string>("-scalar", "Testing Scalar Array", argc, argv, "", myRank, unused);
  double   contourVal = GetParameter<double>("-contour", "Contour Value", argc, argv, 0.0, myRank, unused);
  bool      imageTest = GetParameter<bool>("-imagetest", "imageTest", argc, argv, 0, myRank, unused);
  std::string   imageScalars = GetParameter<std::string>("-imageScalars", "Image Scalar Array", argc, argv, "", myRank, unused);
  unused = GetArrayParameter<double>("-value_range", "Expected Value Range", vminmax, 2, argc, argv, myRank);
  unused = GetArrayParameter<double>("-peak_position", "Expected Peak Position", vpos, 3, argc, argv, myRank);
  //
  // Window/Camera
  //
  cameraSet = GetArrayParameter<double>("-cameraPosition", "Camera Position", cameraPosition, 3, argc, argv, myRank);
  unused = GetArrayParameter<double>("-cameraFocus", "Camera Focus", cameraFocus, 3, argc, argv, myRank);
  unused = GetArrayParameter<double>("-cameraViewUp", "Camera ViewUp", cameraViewUp, 3, argc, argv, myRank);
  unused = GetArrayParameter<int>("-windowSize", "Window Size", windowSize, 2, argc, argv, myRank);
  DisplayParameter<char *>("--------------------", "", &empty, 1, myRank);
  
  // bug fix for cmd line params on windows with debugger (only first read properly)
  gridSpacing[2] = gridSpacing[1] = gridSpacing[0];

  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  if (!fixNeighbours && !fixRadius) {
    std::cerr << "Error : Test requires ParticleSize or Num Neighbours parameter " << std::endl;
    return 1;
  }
  if (scalarname.size()==0) {
    std::cerr << "Error : no scalar name supplied " << std::endl;
    return 1;
  }
  
  //--------------------------------------------------------------
  // Test
  //--------------------------------------------------------------
  vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
  reader->SetFileName(fullname);
  reader->SetController(controller);
  reader->SetGenerateVertexCells(1);

  //--------------------------------------------------------------
  // Update in parallel:
  // To get parallel operation correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------
  //
  vtkSmartPointer<vtkTimerLog> readtimer = vtkSmartPointer<vtkTimerLog>::New();
  readtimer->StartTimer();

  // We can't compute the extents of the image until we know the bounds of the data, 
  // so update reader first.
  vtkStreamingDemandDrivenPipeline *reader_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  reader_sddp->UpdateDataObject();
  reader_sddp->UpdateInformation();
  vtkDebugMacro( "Reader Information Updated " << myRank << " of " << numProcs );
  reader_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  reader_sddp->Update();
  vtkDebugMacro( "Reader Updated " << myRank << " of " << numProcs );
  controller->Barrier();
  readtimer->StopTimer();
  double read_elapsed = readtimer->GetElapsedTime();
  vtkIdType totalParticles = 0;
  vtkIdType localParticles = reader->GetOutput()->GetNumberOfPoints();
  controller->AllReduce(&localParticles, &totalParticles, 1, vtkCommunicator::SUM_OP);

  vtkSmartPointer<vtkAlgorithm> data_algorithm = reader; 

#ifdef PV_MESHLESS_TRILINOS
  vtkSmartPointer<vtkTimerLog> partitiontimer = vtkSmartPointer<vtkTimerLog>::New();
  partitiontimer->StartTimer();
  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  vtkDebugMacro( "Creating Partitioner " << myRank << " of " << numProcs );
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(reader->GetOutputPort());
  partitioner->SetGhostCellOverlap(ghost);
  partitioner->SetController(controller);
  //
  vtkStreamingDemandDrivenPipeline *partition_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(partitioner->GetExecutive());
  partition_sddp->UpdateDataObject();
  vtkDebugMacro( "Partition DataObject Updated " << myRank << " of " << numProcs );
  partition_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  partition_sddp->UpdateInformation();
  vtkDebugMacro( "Partition Information Updated " << myRank << " of " << numProcs );
  vtkDebugMacro( "Partition Update coming " << myRank << " of " << numProcs );
  partition_sddp->Update();
  vtkDebugMacro( "Partition Updated " << myRank << " of " << numProcs );
  
  controller->Barrier();
  partitiontimer->StopTimer();
  double partition_elapsed = partitiontimer->GetElapsedTime();
  //
  data_algorithm = partitioner;
#endif

  vtkSmartPointer<vtkTimerLog> sphtimer = vtkSmartPointer<vtkTimerLog>::New();
  sphtimer->StartTimer();

  vtkDebugMacro( "Creating SPHManager " << myRank << " of " << numProcs );
  controller->Barrier();
  vtkSmartPointer<vtkSPHManager> sphManager = vtkSmartPointer<vtkSPHManager>::New();
  sphManager->SetDefaultDensity(1000.0);
  if (fixNeighbours) {
    sphManager->SetInterpolationMethodToLinear();
    sphManager->SetMaximumNeighbours(maxN);
    sphManager->SetKernelDimension(3);
    sphManager->SetDefaultParticleSideLength(0.0);
    sphManager->SetHCoefficient(0.0);
  }
  else if (fixRadius) {
    sphManager->SetInterpolationMethodToKernel();
    sphManager->SetKernelDimension(3);
    sphManager->SetKernelTypeToCubicSpline();
    sphManager->SetHCoefficient(1.5);
    sphManager->SetDefaultParticleSideLength(particleSize);
    sphManager->SetMaximumSearchRadius(particleSize*1.5*3.0); 
  }

  vtkSmartPointer<vtkAlgorithm> resample_algorithm;
  int wholeExtent[6]={0,-1,0,-1,0,-1};
  if (imageTest) {
    vtkSmartPointer<vtkSPHImageResampler> sphProbe = vtkSmartPointer<vtkSPHImageResampler>::New();
    sphProbe->SetInputConnection(data_algorithm->GetOutputPort());
    if (gridSpacing[0]>0.0) {
      sphProbe->SetSpacing(gridSpacing);
      sphProbe->SetResolution(0,0,0);
    }
    else {
      sphProbe->SetResolution(32,64,32);
    }
    sphProbe->SetDelta(particleSize);
    sphProbe->SetSPHManager(sphManager);
    if (massScalars.size()) {
      sphProbe->SetMassScalars(massScalars.c_str());
    }
    if (densityScalars.size()) {
      sphProbe->SetDensityScalars(densityScalars.c_str());
    }
    sphProbe->SetComputeDensityFromNeighbourVolume(1);
    sphProbe->SetController(controller);
    resample_algorithm = sphProbe;
  }
  else {
    vtkSmartPointer<vtkSPHProbeFilter> sphProbe = vtkSmartPointer<vtkSPHProbeFilter>::New();
    sphProbe->SetInputConnection(data_algorithm->GetOutputPort());
    sphProbe->SetProbeConnection(data_algorithm->GetOutputPort());
    sphProbe->SetSPHManager(sphManager);
    if (massScalars.size()) {
      sphProbe->SetMassScalars(massScalars.c_str());
    }
    if (densityScalars.size()) {
      sphProbe->SetDensityScalars(densityScalars.c_str());
    }
    sphProbe->SetComputeDensityFromNeighbourVolume(1);
    resample_algorithm = sphProbe;
  }

  //
  // Update SPH Pipeline
  //
  vtkStreamingDemandDrivenPipeline *resample_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(resample_algorithm->GetExecutive());
  vtkDebugMacro( "Setting resample piece information " << myRank << " of " << numProcs );
  resample_sddp->UpdateDataObject();
  resample_sddp->UpdateInformation();
  resample_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  resample_sddp->Update();

  controller->Barrier();
  sphtimer->StopTimer();
  double sph_elapsed = sphtimer->GetElapsedTime();
  if (myRank==0) {
    vtkDebugMacro( "Probe completed in " << sph_elapsed << " seconds" );
  }

  vtkSmartPointer<vtkTimerLog> viztimer = vtkSmartPointer<vtkTimerLog>::New();
  viztimer->StartTimer();

  //--------------------------------------------------------------
  // Fetch smoothed output for testing results
  //--------------------------------------------------------------
  vtkDataSet *sph_results = vtkDataSet::SafeDownCast(resample_sddp->GetOutputData(0));
  if (imageTest) {
    vtkImageData::SafeDownCast(sph_results)->GetWholeExtent(wholeExtent);
  }

  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  bool ok = true;
  std::cout.precision(15);
  
  //
  // Testing use scalar values extracted from filters
  //
  if (!imageTest) {
    vtkPointData *sph_pd = sph_results->GetPointData();
    vtkDataArray *scalar_array = sph_pd->GetArray(scalarname.c_str());
    //--------------------------------------------------------------
    // Collect results from parallel filters
    //--------------------------------------------------------------
    double scalar_range_local[2], scalar_range_global[2];
    double scalar_pos[3];
    //
    scalar_array->GetRange(scalar_range_local);
    controller->AllReduce(&scalar_range_local[0], &scalar_range_global[0], 1, vtkCommunicator::MIN_OP);
    controller->AllReduce(&scalar_range_local[1], &scalar_range_global[1], 1, vtkCommunicator::MAX_OP);
    if (myRank==0) DisplayParameter<double>(scalarname.c_str(), " Found Min/Max", scalar_range_global, 2, myRank);
    vtkIdType index = -1;
    for (vtkIdType i=0; i<scalar_array->GetNumberOfTuples(); i++) {
      if (scalar_array->GetTuple1(i)==scalar_range_global[1]) {
        index = i;
        break;
      }
    }
    // Since we computed the density in parallel, only one process will find the peak value
    // the others will have a -1 as the index of the peak.
    if (index!=-1) {
      vtkIdType ids[2] = {myRank, index};
      DisplayParameter<vtkIdType>(scalarname.c_str(), " Peak Proc,Index", ids, 2, myRank);
      sph_results->GetPoint(index,scalar_pos);
      DisplayParameter<double>(scalarname.c_str(), " Peak Position", scalar_pos, 3, myRank);
      //
      double tol_min = std::fabs(vminmax[0]/1000.0);
      double tol_max = std::fabs(vminmax[1]/1000.0);
      if (std::fabs(vminmax[0]-scalar_range_global[0])>tol_min || std::fabs(vminmax[1]-scalar_range_global[1])>tol_max) {
        ok = false;
        DisplayParameter<char *>("++++++++++++++++++++", "", &empty, 1, myRank);
        std::cout << "min/max check failed " << std::endl;
        std::cout << "expected {" << vminmax[0] << ',' << vminmax[1] << "}" << std::endl;
        std::cout << "got {" << scalar_range_global[0] << ',' << scalar_range_global[1] << "}" << std::endl;
        std::cout << "err {" << std::abs(vminmax[0]-scalar_range_global[0]) << ',' << std::abs(vminmax[1]-scalar_range_global[1]) << "}" << std::endl;
      }
      if (std::fabs(vpos[0]-scalar_pos[0])>1E-5 ||
          std::fabs(vpos[1]-scalar_pos[1])>1E-5 ||
          std::fabs(vpos[2]-scalar_pos[2])>1E-5)
      {
        ok = false;
        DisplayParameter<char *>("++++++++++++++++++++", "", &empty, 1, myRank);
        std::cout << "position check failed " << std::endl;
        std::cout << "expected {" << vpos[0] << "," << vpos[1] << ',' << vpos[2] << "}" << std::endl;
        std::cout << "got {" << scalar_pos[0] << "," << scalar_pos[1] << ',' << scalar_pos[2] << "}" << std::endl;
        std::cout << "err {" << std::abs(scalar_pos[0]-vpos[0]) << ',' << std::abs(scalar_pos[1]-vpos[1]) << ',' << std::abs(scalar_pos[2]-vpos[2]) << "}" << std::endl;
      }
    }
  }
  //
  // Testing using the standard VTK Image Regression
  //
  else {
    //
    vtkSmartPointer<vtkContourFilter> iso = vtkSmartPointer<vtkContourFilter>::New();
    iso->SetInputConnection(resample_algorithm->GetOutputPort());
    iso->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,scalarname.c_str());
    iso->SetValue(0, contourVal);
    iso->ComputeScalarsOff();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();
    //
    //
    vtkSmartPointer<vtkAlgorithm> vis_algorithm; 
    vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
    processId->SetInputConnection(iso->GetOutputPort());
    processId->SetController(controller);
    //
    vis_algorithm = processId;
    //
    vtkStreamingDemandDrivenPipeline *vis_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(vis_algorithm->GetExecutive());
    // no piece info set yet, assumes info is not piece dependent
    vtkDebugMacro( "Setting viz piece information " << myRank << " of " << numProcs );
    vis_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
    vtkDebugMacro( "Viz UpdateInformation " );
    vis_sddp->UpdateInformation();
    vtkDebugMacro( "Viz Update " );
    vis_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
    vis_sddp->Update();
    vtkDebugMacro( "Viz Done " );
    //
    vtkSmartPointer<vtkPolyData> contourData = vtkSmartPointer<vtkPolyData>::New();
    contourData->ShallowCopy(vis_sddp->GetOutputData(0));

    //
    // Release all the probing pipeline to free memory so that process zero
    // is less likely to crash when collecting results.
    //
    partitioner->SetInputConnection(NULL);
    resample_algorithm->SetInputConnection(NULL);
    iso->SetInputConnection(NULL);
    processId->SetInputConnection(NULL);
    reader             = NULL;
    partitioner        = NULL;
    sphManager         = NULL;
    resample_algorithm = NULL;
    iso                = NULL;
    processId          = NULL;
    controller->Barrier();

    //
    // Rank >0 Send contour pieces to rank 0
    //
    if (myRank>0) {
      const char *array_name = "ProcessId";
      vtkDataArray *da = contourData->GetPointData()->GetArray(array_name);
      if (da) {      
        double *rg = da->GetRange();
        vtkDebugMacro( "Range (Sending) " << rg[0] << " " << rg[1] );
      }
      controller->Send(contourData, 0, ISO_OUTPUT_TAG);
    }
    //
    // Rank 0 collect all contour pieces from parallel processes
    //
    else if (myRank==0) {
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      vtkSmartPointer<vtkAppendPolyData>       append = vtkSmartPointer<vtkAppendPolyData>::New();
      //
      for (int i=0; i<numProcs; i++) {
        vtkSmartPointer<vtkPolyData> pd;
        if (i==0) {
          pd = contourData;
        }
        else {
          pd = vtkSmartPointer<vtkPolyData>::New();
          controller->Receive(pd, i, ISO_OUTPUT_TAG);
        }
        append->AddInput(pd);
        const char *array_name = "ProcessId";
        vtkDataArray *da = pd->GetPointData()->GetArray(array_name);
        if (da) {
          double *rg = da->GetRange();
          vtkDebugMacro( "Range (Receive) " << rg[0] << " " << rg[1] );
        }
        //
        // Display boxes for each partition
        //
        vtkSmartPointer<vtkOutlineFilter> boxsource = vtkSmartPointer<vtkOutlineFilter>::New();
        boxsource->SetInput(pd);
        vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>          bactor = vtkSmartPointer<vtkActor>::New();
        bmapper->SetInputConnection(boxsource->GetOutputPort());
        bmapper->SetImmediateModeRendering(1);
        bactor->SetMapper(bmapper);
        ren->AddActor(bactor);
      }
      vtkSmartPointer<vtkPolyDataNormals> norm = vtkSmartPointer<vtkPolyDataNormals>::New();
      norm->SetInputConnection(append->GetOutputPort());
      norm->Update();
      //
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize(windowSize);
      mapper->SetInputConnection(norm->GetOutputPort());
      mapper->SetImmediateModeRendering(1);
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointFieldData();
      mapper->SetInterpolateScalarsBeforeMapping(0);
      mapper->SetUseLookupTableScalarRange(0);
      vtkDataArray *da = append->GetOutput()->GetPointData()->GetArray(imageScalars.c_str());
      mapper->SetScalarRange(da->GetRange());
      mapper->SelectColorArray(imageScalars.c_str());
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      if (cameraSet) {
        ren->GetActiveCamera()->SetPosition(cameraPosition);
        ren->GetActiveCamera()->SetFocalPoint(cameraFocus);
        ren->GetActiveCamera()->SetViewUp(cameraViewUp);
        ren->ResetCameraClippingRange();
      }
      else {
        ren->ResetCamera();
      }
      //
      vtkDebugMacro( "Process Id : " << myRank << " About to Render" );
      renWindow->Render();
      //
      retVal = vtkRegressionTester::Test(argc, argv, renWindow, 45);
      if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
        iren->Start();
      }
      vtkDebugMacro( "Process Id : " << myRank << " Rendered" );
      ok = (retVal==vtkRegressionTester::PASSED);
    }
  }
  controller->Barrier();
  viztimer->StopTimer();
  double viz_elapsed = viztimer->GetElapsedTime();

  if (ok && myRank==0) {
    DisplayParameter<char *>("____________________", "", &empty, 1, myRank);
    DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, myRank);
    DisplayParameter<double>("Read Time", "", &read_elapsed, 1, myRank);
    DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, myRank);
    DisplayParameter<double>("SPH Probe Time", "", &sph_elapsed, 1, myRank);
    DisplayParameter<double>("Visualization/Check Time", "", &viz_elapsed, 1, myRank);
    if (imageTest) {
      DisplayParameter<int>("Image Whole Extent", "", wholeExtent, 6, myRank);
      vtkIdType voxels = (1+wholeExtent[1]-wholeExtent[0])*(1+wholeExtent[3]-wholeExtent[2])*(1+wholeExtent[5]-wholeExtent[4]);
      DisplayParameter<vtkIdType>("Voxels", "", &voxels, 1, myRank);
    }
    DisplayParameter<char *>("====================", "", &empty, 1, myRank);
  }

  retVal = ok;
  
  delete []fullname;
  delete []filename;
  controller->Finalize();

  return !retVal;
}
//----------------------------------------------------------------------------

