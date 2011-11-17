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
#include "vtkH5PartWriter.h"
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

//----------------------------------------------------------------------------
static const int ISO_OUTPUT_TAG=301;
//----------------------------------------------------------------------------

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
      std::cout.width(30);
      std::cout << message << " {" << newValue << "}" << std::endl;
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
      std::cout << message << " {";
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
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-F", argc, argv, "DUMMY_ENV_VAR", "temp.h5");
  char* fullname = vtkTestUtilities::ExpandDataFileName(argc, argv, filename);
  if (myRank==0) {
    std::cout << "Process Id : " << myRank << " FileName : " << fullname << std::endl;
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
  bool unused, fixNeighbours, fixRadius;
  double gridSpacing[3] = {0.0, 0.0, 0.0};
  double vminmax[2] = {0.0, 0.0};
  double vpos[3] = {0.0, 0.0, 0.0};
  double cameraPosition[3] = {0.0, 0.0, 0.0};
  double cameraFocus[3] = {0.0, 1.0, 0.0};
  double cameraViewUp[3] = {0.0, 0.0, 1.0};

  //
  int            maxN = GetParameter<int>("-neighbours", "Fixed Neighbours", argc, argv, 0, myRank, fixNeighbours);
  double particleSize = GetParameter<double>("-particlesize", "Particle Size", argc, argv, 0, myRank, fixRadius);
  double        ghost = GetParameter<double>("-ghost_region", "Ghost Region", argc, argv, 0.0, myRank, unused);
  double   contourVal = GetParameter<double>("-contour", "Contour Value", argc, argv, 0.0, myRank, unused);
  bool      ImageTest = GetParameter<bool>("-imagetest", "ImageTest", argc, argv, 0, myRank, unused);

  std::string     scalarname = GetParameter<std::string>("-scalar", "Testing Scalar array", argc, argv, "", myRank, unused);
  std::string    massScalars = GetParameter<std::string>("-massScalars", "Mass Scalar array", argc, argv, "", myRank, unused);
  std::string densityScalars = GetParameter<std::string>("-densityScalars", "Density Scalar array", argc, argv, "", myRank, unused);

  unused = GetArrayParameter<double>("-gridSpacing", "Grid Spacing", gridSpacing, 3, argc, argv, myRank);
  unused = GetArrayParameter<double>("-value_range", "Test Valid : value_range", vminmax, 2, argc, argv, myRank);
  unused = GetArrayParameter<double>("-peak_position", "Test Valid : scalar_peak_position", vpos, 3, argc, argv, myRank);

  unused = GetArrayParameter<double>("-cameraPosition", "Camera Position", cameraPosition, 3, argc, argv, myRank);
  unused = GetArrayParameter<double>("-cameraFocus", "Camera Focus", cameraFocus, 3, argc, argv, myRank);
  unused = GetArrayParameter<double>("-cameraViewUp", "Camera ViewUp", cameraViewUp, 3, argc, argv, myRank);

  gridSpacing[2] = gridSpacing[1] = gridSpacing[0];

  if (myRank==0) {
    std::cout << "//--------------------------------------------------------------" << std::endl;
    std::cout << "//" << std::endl;
    std::cout << "//--------------------------------------------------------------" << std::endl;
  }
  controller->Barrier();
  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  if (!fixNeighbours && !fixRadius) {
    std::cout << "Error : Test requires ParticleSize or Num Neighbours parameter " << std::endl;
    return 1;
  }
  if (scalarname.size()==0) {
    std::cout << "Error : no scalar name supplied " << std::endl;
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

  // We can't compute the extents of the image until we know the bounds of the data, 
  // so update reader first.
  vtkStreamingDemandDrivenPipeline *reader_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  reader_sddp->UpdateDataObject();
  reader_sddp->UpdateInformation();
  std::cout << "Reader Information Updated " << myRank << " of " << numProcs << std::endl;
  reader_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  reader_sddp->Update();
  std::cout << "Reader Updated " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();

  vtkSmartPointer<vtkAlgorithm> data_algorithm = reader; 
#ifdef PV_MESHLESS_TRILINOS
  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  std::cout << "Creating Partitioner " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(reader->GetOutputPort());
  partitioner->SetGhostCellOverlap(ghost);
  partitioner->SetController(controller);
  //
  vtkStreamingDemandDrivenPipeline *partition_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(partitioner->GetExecutive());
  partition_sddp->UpdateDataObject();
  std::cout << "Partition DataObject Updated " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();
  partition_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  partition_sddp->UpdateInformation();
  std::cout << "Partition Information Updated " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();
  std::cout << "Partition Update coming " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();
  partition_sddp->Update();
  std::cout << "Partition Updated " << myRank << " of " << numProcs << std::endl;
  controller->Barrier();
  //
  data_algorithm = partitioner;
#endif

  std::cout << "Creating SPHManager " << myRank << " of " << numProcs << std::endl;
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
  if (ImageTest) {
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
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  //
  // Update SPH Pipeline
  //
  vtkStreamingDemandDrivenPipeline *resample_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(resample_algorithm->GetExecutive());
  std::cout << "Setting resample piece information " << myRank << " of " << numProcs << std::endl;
  resample_sddp->UpdateDataObject();
  resample_sddp->UpdateInformation();
  resample_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  controller->Barrier();
  resample_sddp->Update();

  //
  timer->StopTimer();
  controller->Barrier();
  if (myRank==0) {
    double elapsed = timer->GetElapsedTime();
    std::cout << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    std::cout << "Probe completed in " << elapsed << " seconds" << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
  }

  //--------------------------------------------------------------
  // Fetch smoothed output for testing results
  //--------------------------------------------------------------
  vtkDataSet *sph_results = vtkDataSet::SafeDownCast(resample_sddp->GetOutputData(0));

  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  bool ok = true;
  std::cout.precision(15);
  
  //
  // Testing use scalar values extracted from filters
  //
  if (!ImageTest) {
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
    if (myRank==0) {
      std::cout << "Min and Max of SmoothedDensity are {" << scalar_range_global[0] << "," << scalar_range_global[1] << "}" << std::endl;
    }
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
      std::cout << "Index of peak SmoothedDensity particle is " << index << " on process " << myRank << std::endl;
      sph_results->GetPoint(index,scalar_pos);
      std::cout << "Position of peak SmoothedDensity particle is {" << scalar_pos[0] << "," << scalar_pos[1] << "," << scalar_pos[2] << "}" << std::endl;
      //
      double tol_min = std::fabs(vminmax[0]/1000.0);
      double tol_max = std::fabs(vminmax[1]/1000.0);
      if (std::fabs(vminmax[0]-scalar_range_global[0])>tol_min || std::fabs(vminmax[1]-scalar_range_global[1])>tol_max) {
        ok = false;
        std::cout << "//--------------------------------------------------------------" << std::endl;
        std::cout << "min/max check failed " << std::endl;
        std::cout << "expected {" << vminmax[0] << ',' << vminmax[1] << "}" << std::endl;
        std::cout << "got {" << scalar_range_global[0] << ',' << scalar_range_global[1] << "}" << std::endl;
        std::cout << "err {" << std::abs(vminmax[0]-scalar_range_global[0]) << ',' << std::abs(vminmax[1]-scalar_range_global[1]) << "}" << std::endl;
        std::cout << "//--------------------------------------------------------------" << std::endl;
      }
      if (std::fabs(vpos[0]-scalar_pos[0])>1E-5 ||
          std::fabs(vpos[1]-scalar_pos[1])>1E-5 ||
          std::fabs(vpos[2]-scalar_pos[2])>1E-5)
      {
        ok = false;
        std::cout << "//--------------------------------------------------------------" << std::endl;
        std::cout << "position check failed " << std::endl;
        std::cout << "expected {" << vpos[0] << "," << vpos[1] << ',' << vpos[2] << "}" << std::endl;
        std::cout << "got {" << scalar_pos[0] << "," << scalar_pos[1] << ',' << scalar_pos[2] << "}" << std::endl;
        std::cout << "err {" << std::abs(scalar_pos[0]-vpos[0]) << ',' << std::abs(scalar_pos[1]-vpos[1]) << ',' << std::abs(scalar_pos[2]-vpos[2]) << "}" << std::endl;
        std::cout << "//--------------------------------------------------------------" << std::endl;
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
    std::cout << "Setting viz piece information " << myRank << " of " << numProcs << std::endl;
    vis_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
    std::cout << "Viz UpdateInformation " << std::endl;
    vis_sddp->UpdateInformation();
    std::cout << "Viz Update " << std::endl;
    vis_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
    vis_sddp->Update();
    std::cout << "Viz Done " << std::endl;
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
      double *rg = da->GetRange();
      std::cout << "Range (Sending) " << rg[0] << " " << rg[1] << std::endl;
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
        double *rg = da->GetRange();
        std::cout << "Range (Receive) " << rg[0] << " " << rg[1] << std::endl;
        //
        // Display boxes for each partition
        //
        vtkSmartPointer<vtkOutlineFilter> boxsource = vtkSmartPointer<vtkOutlineFilter>::New();
        boxsource->SetInput(pd);
        vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>          bactor = vtkSmartPointer<vtkActor>::New();
        bmapper->SetInputConnection(boxsource->GetOutputPort());
        bactor->SetMapper(bmapper);
        ren->AddActor(bactor);
      }
      vtkSmartPointer<vtkPolyDataNormals> norm = vtkSmartPointer<vtkPolyDataNormals>::New();
      norm->SetInputConnection(append->GetOutputPort());
      norm->Update();
      //
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInputConnection(norm->GetOutputPort());
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointFieldData();
      mapper->SetInterpolateScalarsBeforeMapping(0);
      mapper->SetUseLookupTableScalarRange(0);
//      const char *array_name = "FLUID_pressure";
      const char *array_name = "ProcessId";
      vtkDataArray *da = append->GetOutput()->GetPointData()->GetArray(array_name);
      mapper->SetScalarRange(da->GetRange());
      mapper->SelectColorArray(array_name);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      ren->GetActiveCamera()->SetPosition(cameraPosition);
      ren->GetActiveCamera()->SetFocalPoint(cameraFocus);
      ren->GetActiveCamera()->SetViewUp(cameraViewUp);
      ren->ResetCamera();
      //
      std::cout << "Process Id : " << myRank << " About to Render" << std::endl;
      renWindow->Render();
      //
      retVal = vtkRegressionTester::Test(argc, argv, renWindow, 10);
      if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
        iren->Start();
      }
      std::cout << "Process Id : " << myRank << " Rendered" << std::endl;
      ok = (retVal==vtkRegressionTester::PASSED);
    }
  }
  retVal = ok;
  
  delete []fullname;
  delete []filename;
  controller->Finalize();

  return !retVal;
}
//----------------------------------------------------------------------------

