#define _USE_MATH_DEFINES
#include <cmath>
//
#include "TestUtils.h"
//
#ifdef _WIN32
  #include <windows.h>
  #undef min
  #undef max
#else 
  #include <sys/time.h>
#endif
//
#include <algorithm>
#include <vtksys/SystemTools.hxx>
//
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
 #include "vtkParticlePartitionFilter.h"
#endif
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

//
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#ifdef _WIN32
void known_seed() {
  srand(1);
}
#else
void known_seed() {
  srandom(12345);
}
#endif
//----------------------------------------------------------------------------
#ifdef _WIN32
unsigned long int random_seed()
{
  LARGE_INTEGER lpPerformanceCount;
  QueryPerformanceCounter(&lpPerformanceCount);
  long int seed = lpPerformanceCount.LowPart + lpPerformanceCount.HighPart;
  srand(seed);
  return seed;
}
#else
unsigned long int random_seed()
{
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;
  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } 
  else {
    if (fread(&seed,sizeof(seed),1,devrandom) == 1) {
      fclose(devrandom);
    } 
    else {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
    }
  }
  srandom(seed);
  return seed;
}
#endif
//----------------------------------------------------------------------------
void SpherePoints(int n, float radius, float X[]) {
  double x, y, z, w, t;
  Random r(12345);
  double rmin=1E6, rmax=-1E6;
  for(int i=0; i<n; i++ ) {
    double r1 = r.nextNumber(); // double(rand())/RAND_MAX;
    double r2 = r.nextNumber(); // double(rand())/RAND_MAX;
    rmin = std::min(rmin, r1);
    rmin = std::min(rmin, r2);
    rmax = std::max(rmax, r1);
    rmax = std::max(rmax, r2);
    z = 2.0 * r1 - 1.0;
    t = 2.0 * M_PI * r2;
    w = radius * sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    X[3*i+0] = static_cast<float>(x);
    X[3*i+1] = static_cast<float>(y);
    X[3*i+2] = static_cast<float>(z*radius);
  }
}
//----------------------------------------------------------------------------
void finalizeTest(TestStruct &test)
{
/*
  test.reader->SetController(NULL);
  test.partitioner->SetController(NULL);
  if (test.imageResample) {
    vtkSmartPointer<vtkSPHImageResampler> sphProbe = vtkSPHImageResampler::SafeDownCast(test.sphResampler);
    sphProbe->SetController(NULL);
  }
  else {
    vtkSmartPointer<vtkSPHProbeFilter> sphProbe = vtkSPHProbeFilter::SafeDownCast(test.sphResampler);
    sphProbe->SetController(NULL);
  }
*/
  //
  test.reader       = NULL;
  test.partitioner  = NULL;
  test.sphManager   = NULL;
  test.sphResampler = NULL;
  test.controller   = NULL;
  vtkMultiProcessController::SetGlobalController(NULL);
}
//----------------------------------------------------------------------------
int initTest(int argc, char* argv[], TestStruct &test)
{
  char *empty = "";
#ifdef PARAVIEW_USE_MPI
  int dummy=0;
//  DisplayParameter<char *>("====================", "Init MPI", &empty, 1, 0);
  MPI_Init(&argc,&argv);
  test.controller = vtkSmartPointer<vtkMPIController>::New();
  test.controller->Initialize(&argc, &argv, 1);
  vtkMultiProcessController::SetGlobalController(test.controller);
#else
  test.controller = vtkDummyController::New();
#endif
  
  bool unused;
  // Obtain the id of the running process and the total number of processes
  test.myRank = test.controller->GetLocalProcessId();
  test.numProcs = test.controller->GetNumberOfProcesses();
  //
  test.gridSpacing[0] = test.gridSpacing[1] = test.gridSpacing[2] = 0.0;
  test.gridResolution[0] = test.gridResolution[1] = test.gridResolution[2] = -1;
  test.vminmax[0] = test.vminmax[1] = 0.0;
  test.vpos[0] = test.vpos[1] = test.vpos[2] = 0.0;
  test.cameraPosition[0] = test.cameraPosition[1] = test.cameraPosition[2] = 0.0;
  test.cameraFocus[0] = test.cameraFocus[1] = test.cameraFocus[2] = 0.0;
  test.cameraViewUp[0] = 0.0;
  test.cameraViewUp[1] = 0.0;
  test.cameraViewUp[2] = 1.0;
  test.windowSize[0] = test.windowSize[1] = 400;

  // uncomment this to wait for debugger (needs to be after MPI_Init because it checks rank
  //DEBUG_WAIT
  //
  test.controller->Barrier();

  //--------------------------------------------------------------
  // command line params : Setup testing utilities/args etc
  //--------------------------------------------------------------
  vtkSmartPointer<vtkTesting> vtktest = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    vtktest->AddArgument(argv[c]);
  }

  //
  // Force the creation of our output window object
  //
//  vtkSmartPointer<vtkStreamOutputWindow> outwin = vtkSmartPointer<vtkStreamOutputWindow>::New();
//  vtkOutputWindow::SetInstance(outwin);
//  outwin->SetOutputStream(&std::cout);

  //
  // General test flags/info
  //
  DisplayParameter<char *>("====================", "", &empty, 1, (test.myRank==0)?0:-1);
  test.testName = GetParameter<std::string>("-testName", "Test name", argc, argv, "", test.myRank, unused);
  test.doRender = GetParameter<bool>("-doRender", "Enable Render", argc, argv, 0, test.myRank, unused);
  test.keepTempFiles = GetParameter<bool>("-X", "Keep Temporary Files", argc, argv, 0, test.myRank, unused);

  //
  // ParticleGenerate info
  //
  test.generateN = GetParameter<vtkIdType>("-generateParticles", "Generated particles", argc, argv, 0, test.myRank, unused);
  test.memoryMB  = GetParameter<vtkIdType>("-memoryMB", "Maximum memory (MB)", argc, argv, 0, test.myRank, unused);
  test.pieceValidation = GetParameter<bool>("-pieceValidation", "Validate by pieces", argc, argv, false, test.myRank, unused);
  test.iterations = GetParameter<int>("-iterations", "Number of iterations", argc, argv, 1, test.myRank, unused);

  //
  // H5Part info
  //
  std::string filename = GetParameter<std::string>("-F", "Filename", argc, argv, "", test.myRank, unused);
  std::string filepath = GetParameter<std::string>("-D", "Filepath", argc, argv, "", test.myRank, unused);
  if (filename.size() && filepath.size()) {
    test.fullName = vtksys::SystemTools::ConvertToOutputPath(std::string(filepath+"/"+filename).c_str());
    DisplayParameter<std::string>("FullName", "", &test.fullName, 1, (test.myRank==0)?0:-1);
  }
  test.Xarray = GetParameter<std::string>("-Xarray", "Xarray name", argc, argv, "", test.myRank, unused);
  test.Yarray = GetParameter<std::string>("-Yarray", "Yarray name", argc, argv, "", test.myRank, unused);
  test.Zarray = GetParameter<std::string>("-Zarray", "Zarray name", argc, argv, "", test.myRank, unused);
  test.ignorePartitions = GetParameter<bool>("-ignorePartitions", "Ignore Partitions", argc, argv, 0, test.myRank, unused);
  test.randomizeExtents = GetParameter<bool>("-randomizeExtents", "Randomize Extents", argc, argv, 0, test.myRank, unused);

  //
  // SPH kernel or neighbour info
  //
  test.particleSize = GetParameter<double>("-particlesize", "Particle Size", argc, argv, 0, test.myRank, test.fixRadius);
  test.ghostOverlap = GetParameter<double>("-ghost_region", "Ghost Region", argc, argv, 0.0, test.myRank, unused);
  GetArrayParameter<double>("-gridSpacing", "Grid Spacing", test.gridSpacing, 3, argc, argv, test.myRank);
  GetArrayParameter<int>("-gridResolution", "Grid Resolution", test.gridResolution, 3, argc, argv, test.myRank);
  test.maxN = GetParameter<int>("-neighbours", "Fixed Neighbours", argc, argv, 0, test.myRank, test.fixNeighbours);
  test.massScalars = GetParameter<std::string>("-massScalars", "Mass Scalar Array", argc, argv, "", test.myRank, unused);
  test.densityScalars = GetParameter<std::string>("-densityScalars", "Density Scalar Array", argc, argv, "", test.myRank, unused);
  test.expectedN = GetParameter<vtkIdType>("-expectedparticles", "Simulation Particles", argc, argv, 0, test.myRank, unused);

  //
  // Test/Display of results
  //
  test.scalarname = GetParameter<std::string>("-scalar", "Testing Scalar Array", argc, argv, "", test.myRank, unused);
  test.contourVal = GetParameter<double>("-contour", "Contour Value", argc, argv, 0.0, test.myRank, unused);
  test.imageResample = GetParameter<bool>("-imageResample", "imageResample", argc, argv, 0, test.myRank, unused);
  test.skipImageTest = GetParameter<bool>("-skipImageTest", "skipImageTest", argc, argv, 0, test.myRank, unused);
  test.imageScalars = GetParameter<std::string>("-imageScalars", "Image Scalar Array", argc, argv, "", test.myRank, unused);
  GetArrayParameter<double>("-value_range", "Expected Value Range", test.vminmax, 2, argc, argv, test.myRank);
  GetArrayParameter<double>("-peak_position", "Expected Peak Position", test.vpos, 3, argc, argv, test.myRank);
  test.imageThreshold = GetParameter<int>("-imageThreshold", "Image Threshold Pass/Fail", argc, argv, 1, test.myRank, unused);
  test.benchmarkPartition = GetParameter<bool>("-benchmarkPartition", "benchmarkPartition", argc, argv, 0, test.myRank, unused);
  test.numNodes = GetParameter<int>("-numNodes", "Number of Nodes used (info)", argc, argv, 1, test.myRank, unused);
  test.processesPerNode = GetParameter<int>("-processesPerNode", "Processes per node (info)", argc, argv, 1, test.myRank, unused);

  //
  // Window/Camera
  //
  test.cameraSet = GetArrayParameter<double>("-cameraPosition", "Camera Position", test.cameraPosition, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<double>("-cameraFocus", "Camera Focus", test.cameraFocus, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<double>("-cameraViewUp", "Camera ViewUp", test.cameraViewUp, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<int>("-windowSize", "Window Size", test.windowSize, 2, argc, argv, test.myRank);
  if (unused) { // why have window sizes changed?
//    test.windowSize[0] += 8;
//    test.windowSize[1] += 8;
  }
  // bug fix for cmd line params on windows with debugger (only first read properly)
  test.gridSpacing[2] = test.gridSpacing[1] = test.gridSpacing[0];
  test.gridResolution[2] = test.gridResolution[1] = test.gridResolution[0];
  //
  DisplayParameter<char *>("--------------------", "", &empty, 1, (test.myRank==0)?0:-1);
  //
  return 1;
}  
//----------------------------------------------------------------------------
void TestStruct::CreateReader()
{
  this->reader = vtkSmartPointer<vtkH5PartReader>::New();
  this->reader->SetFileName((char*)(this->fullName.c_str()));
  this->reader->SetController(this->controller);
  this->reader->SetGenerateVertexCells(1);
  if (this->Xarray.size()>0) {
    this->reader->SetXarray(this->Xarray.c_str());
  }
  if (this->Yarray.size()>0) {
    this->reader->SetYarray(this->Yarray.c_str());
  }
  if (this->Zarray.size()>0) {
    this->reader->SetZarray(this->Zarray.c_str());
  }
  this->reader->SetIgnorePartitionBoxes(this->ignorePartitions);
  this->reader->SetRandomizePartitionExtents(this->randomizeExtents);
}
//----------------------------------------------------------------------------
double TestStruct::UpdateReader()
{
  vtkSmartPointer<vtkTimerLog> readtimer = vtkSmartPointer<vtkTimerLog>::New();
  readtimer->StartTimer();

  //--------------------------------------------------------------
  // Update in parallel:
  // To get parallel operation correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------

  // We can't compute the extents of the image until we know the bounds of the data, 
  // so update test.reader first.
  vtkStreamingDemandDrivenPipeline *reader_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(this->reader->GetExecutive());
  reader_sddp->UpdateDataObject();
  reader_sddp->UpdateInformation();
  testDebugMacro( "Reader Information Updated " << this->myRank << " of " << this->numProcs );
  reader_sddp->SetUpdateExtent(0, this->myRank, this->numProcs, 0);
  reader_sddp->Update();
  testDebugMacro( "Reader Updated " << this->myRank << " of " << this->numProcs );
  this->controller->Barrier();
  readtimer->StopTimer();
  double read_elapsed = readtimer->GetElapsedTime();
  testDebugMacro( "Reader completed in " << read_elapsed << " seconds" );
  return read_elapsed;
}
//----------------------------------------------------------------------------
void TestStruct::CreatePartitioner()
{
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
  testDebugMacro( "Creating Partitioner " << this->myRank << " of " << this->numProcs );
  this->partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  if (this->reader) {
    this->partitioner->SetInputConnection(this->reader->GetOutputPort());
  }
  this->partitioner->SetGhostCellOverlap(this->ghostOverlap);
  this->partitioner->SetController(this->controller);
#endif
}
//----------------------------------------------------------------------------
double TestStruct::UpdatePartitioner()
{
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
  vtkSmartPointer<vtkTimerLog> partitiontimer = vtkSmartPointer<vtkTimerLog>::New();
  partitiontimer->StartTimer();
  //
  vtkStreamingDemandDrivenPipeline *partition_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(this->partitioner->GetExecutive());
  partition_sddp->UpdateDataObject();
  testDebugMacro( "Partition DataObject Updated " << this->myRank << " of " << this->numProcs );
  partition_sddp->SetUpdateExtent(0, this->myRank, this->numProcs, 0);
  partition_sddp->UpdateInformation();
  testDebugMacro( "Partition Information Updated " << this->myRank << " of " << this->numProcs );
  testDebugMacro( "Partition Update coming " << this->myRank << " of " << this->numProcs );
  partition_sddp->Update();
  testDebugMacro( "Partition Updated " << this->myRank << " of " << this->numProcs );
  this->controller->Barrier();
  partitiontimer->StopTimer();
  double partition_elapsed = partitiontimer->GetElapsedTime();
  testDebugMacro( "Partition completed in " << partition_elapsed << " seconds" );
#else
  double partition_elapsed = 0.0;
#endif
  return partition_elapsed;
}
//----------------------------------------------------------------------------
void TestStruct::DeletePartitioner()
{
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
  this->partitioner->SetInputConnection(NULL);
  this->partitioner = NULL;
#endif
}
//----------------------------------------------------------------------------
void TestStruct::CreateSPHManager()
{
  testDebugMacro( "Creating SPHManager " << this->myRank << " of " << this->numProcs );
  this->sphManager = vtkSmartPointer<vtkSPHManager>::New();
  this->sphManager->SetDefaultDensity(1000.0);
  if (this->fixNeighbours) {
    this->sphManager->SetInterpolationMethodToLinear();
    this->sphManager->SetMaximumNeighbours(this->maxN);
    this->sphManager->SetKernelDimension(3);
    this->sphManager->SetDefaultParticleSideLength(0.0);
    this->sphManager->SetHCoefficient(0.0);
  }
  else if (this->fixRadius) {
    this->sphManager->SetInterpolationMethodToKernel();
    this->sphManager->SetKernelDimension(3);
    this->sphManager->SetKernelTypeToCubicSpline();
    this->sphManager->SetHCoefficient(1.5);
    this->sphManager->SetDefaultParticleSideLength(this->particleSize);
    this->sphManager->SetMaximumSearchRadius(this->particleSize*1.5*3.0); 
  }
}
//----------------------------------------------------------------------------
void TestStruct::CreateSPHResampler(vtkAlgorithm *input)
{
  if (this->imageResample) {
    vtkSmartPointer<vtkSPHImageResampler> sphProbe = vtkSmartPointer<vtkSPHImageResampler>::New();
    sphProbe->SetInputConnection(input->GetOutputPort());
    if (this->gridSpacing[0]>0.0) {
      sphProbe->SetSpacing(this->gridSpacing);
      sphProbe->SetResolution(0,0,0);
    }
    else if (this->gridResolution[0]>-1) {
      sphProbe->SetResolution(this->gridResolution);
    }
    else {
      sphProbe->SetResolution(32,32,32);
    }
    sphProbe->SetDelta(this->particleSize);
    sphProbe->SetSPHManager(this->sphManager);
    if (this->massScalars.size()) {
      sphProbe->SetMassScalars(this->massScalars.c_str());
      sphProbe->SetComputeDensityFromNeighbourVolume(1);
    }
    if (this->densityScalars.size()) {
      sphProbe->SetDensityScalars(this->densityScalars.c_str());
      sphProbe->SetComputeDensityFromNeighbourVolume(1);
    }
    sphProbe->SetController(this->controller);
    this->sphResampler = sphProbe;
  }
  else {
    vtkSmartPointer<vtkSPHProbeFilter> sphProbe = vtkSmartPointer<vtkSPHProbeFilter>::New();
    sphProbe->SetInputConnection(input->GetOutputPort());
    sphProbe->SetProbeConnection(input->GetOutputPort());
    sphProbe->SetSPHManager(this->sphManager);
    if (this->massScalars.size()) {
      sphProbe->SetMassScalars(this->massScalars.c_str());
      sphProbe->SetComputeDensityFromNeighbourVolume(1);
    }
    if (this->densityScalars.size()) {
      sphProbe->SetDensityScalars(this->densityScalars.c_str());
      sphProbe->SetComputeDensityFromNeighbourVolume(1);
    }
    this->sphResampler = sphProbe;
  }
}
//----------------------------------------------------------------------------
double TestStruct::UpdateSPHResampler()
{
  vtkSmartPointer<vtkTimerLog> sphtimer = vtkSmartPointer<vtkTimerLog>::New();
  sphtimer->StartTimer();
  //
  vtkStreamingDemandDrivenPipeline *resample_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(this->sphResampler->GetExecutive());
  testDebugMacro( "Setting resample piece information " << this->myRank << " of " << this->numProcs );
  // make sure data output is setup
  resample_sddp->UpdateDataObject();
  // make sure data extents and types are declared
  resample_sddp->UpdateInformation();
  
  //
  // set extent to be updated to all with pieces
  //
  resample_sddp->SetUpdateExtent(0, this->myRank, this->numProcs, 0);
//  vtkInformation *info = resample_sddp->GetOutputInformation(0);
//  vtkStreamingDemandDrivenPipeline::SetUpdateExtentToWholeExtent(info);
//  info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), this->myRank);
//  info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), this->numProcs);

  //
  // Update resampler
  //
  this->sphResampler->Update();
  //
  sphtimer->StopTimer();
  double sph_elapsed = sphtimer->GetElapsedTime();
  testDebugMacro( "Probe completed in " << sph_elapsed << " seconds" );
  return sph_elapsed;
}
