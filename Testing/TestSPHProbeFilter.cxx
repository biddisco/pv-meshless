// mpiexec : C:\Program Files\MPICH2\bin\mpiexec
// mpiargs : -localonly -n 2 -env PATH D:\cmakebuild\pv-meshless\bin\debug;c:\bin
// appcmd  : $(TargetPath)
// appargs : "-D" "D:/Code/plugins/pv-meshless/Testing/data" "-F" "jet.h5part" "-T" "D:/cmakebuild/plugins/pv-meshless/Testing/Temporary" "-particlesize" "0.001" "-ghost_region" "0.005" "-massScalars" "mass" "-scalar" "ShepardCoeff" "-contour" "0.5" "-densityScalars" "FLUID_density" "-imagetest" "1" "-V" "D:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResampleSerial.png"
// appargs : -D D:/Code/plugins/pv-meshless/Testing/data -F dam-17.h5part -T D:/cmakebuild/plugins/pv-meshless/Testing/Temporary -particlesize 0.005 -ghost_region 0.01 -gridSpacing 0.0025 0.0025 0.0025 -scalar ShepardCoeff -contour 0.1 -densityScalars FLUID_rho -imagetest 1 -V D:/Code/plugins/pv-meshless/Testing/baseline/TestSPHImageResample-large-parallel-4.png
// mpishim : C:\Program Files\Microsoft Visual Studio 9.0\Common7\IDE\Remote Debugger\x64\mpishim.exe

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
std::string usage = "\n"\
"\t-D path to use for h5part file \n" \
"\t-F name to use for h5part file \n" \
"\t-R Render. Displays points using vtk renderwindow \n" \
"\t-I Interactive (waits until user closes render window if -R selected) \n";


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
    std::cout << usage.c_str() << std::endl;
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
  // 
  //--------------------------------------------------------------
  char *tempChar;

  //--------------------------------------------------------------
  // File names/dirs
  //--------------------------------------------------------------
  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-neighbours", argc, argv, "", "");
  int maxN;
  bool fixNeighbours = false;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> maxN;
    if (myRank==0) {
      std::cout << "Fixed Neighbours      {" << maxN << "}" << std::endl;
    }
    fixNeighbours = true;
  }
  delete []tempChar;

  //--------------------------------------------------------------
  // Testing params
  //--------------------------------------------------------------
  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-particlesize", argc, argv, "", "");
  double p_size;
  bool fixRadius = false;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> p_size;
    if (myRank==0) {
      std::cout << "Particle Size         {" << p_size << "}" << std::endl;
    }
    fixRadius = true;
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-scalar", argc, argv, "", "");
  std::string scalarname;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> scalarname;
    if (myRank==0) {
      std::cout << "Testing Scalar array  {" << scalarname << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-massScalars", argc, argv, "", "");
  std::string massScalars;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> massScalars;
    if (myRank==0) {
      std::cout << "Mass Scalar array  {" << massScalars << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-densityScalars", argc, argv, "", "");
  std::string densityScalars;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> densityScalars;
    if (myRank==0) {
      std::cout << "Density Scalar array  {" << densityScalars << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-ghost_region", argc, argv, "", "");
  double ghost = 0.0;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> ghost;
    if (myRank==0) {
      std::cout << "Ghost Region {" << ghost << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-contour", argc, argv, "", "");
  double contourVal = 0.0;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> contourVal;
    if (myRank==0) {
      std::cout << "Contour Value {" << contourVal << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-gridSpacing", argc, argv, "", "");
  double gridSpacing[3] = {0.005, 0.005, 0.005};
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
//    temp >> gridSpacing[0] >> gridSpacing[1] >> gridSpacing[2];
//    gridSpacing[2] = gridSpacing[1] = gridSpacing[0];
    if (myRank==0) {
      std::cout << "Grid Spacing {" << gridSpacing[0] << ',' << gridSpacing[1] << "," << gridSpacing[2] << "}" << std::endl;
    }
  }
  delete []tempChar;

  //--------------------------------------------------------------
  // Valid test results
  //--------------------------------------------------------------
  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-value_range", argc, argv, "", "");
  double vmin=0,vmax=0;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> vmin >> vmax;
    if (myRank==0) {
      std::cout << "Test Valid : value_range {" << vmin << ',' << vmax << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-peak_position", argc, argv, "", "");
  double vpos[3];
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> vpos[0] >> vpos[1] >> vpos[2];
    if (myRank==0) {
      std::cout << "Test Valid : scalar_peak_position {" << vpos[0] << ',' << vpos[1] << "," << vpos[2] << "}" << std::endl;
    }
  }
  delete []tempChar;

  tempChar = vtkTestUtilities::GetArgOrEnvOrDefault("-imagetest", argc, argv, "", "");
  bool ImageTest = false;
  if (std::string(tempChar).size()) {
    vtkstd::stringstream temp(tempChar);
    temp >> ImageTest;
    if (myRank==0) {
      std::cout << "ImageTest : " << ImageTest << std::endl;
    }
  }
  delete []tempChar;

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

  vtkSmartPointer<vtkAlgorithm> data_algorithm = reader; 
#ifdef PV_MESHLESS_TRILINOS
  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(reader->GetOutputPort());
  partitioner->SetGhostCellOverlap(ghost);
  partitioner->SetController(controller);
  data_algorithm = partitioner;
#endif

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
    sphManager->SetDefaultParticleSideLength(p_size);
    sphManager->SetMaximumSearchRadius(p_size*1.5*3.0); 
  }

  vtkSmartPointer<vtkAlgorithm> resample_algorithm; 
  if (ImageTest) {
    vtkSmartPointer<vtkSPHImageResampler> sphProbe = vtkSmartPointer<vtkSPHImageResampler>::New();
    sphProbe->SetInputConnection(data_algorithm->GetOutputPort());
    if (gridSpacing[0]>0.0) {
      sphProbe->SetSpacing(gridSpacing);
    }
    else {
      sphProbe->SetResolution(32,64,32);
    }
    sphProbe->SetDelta(p_size);
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
  //--------------------------------------------------------------
  // Update in parallel:
  // To get parallel operation correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------
  std::cout << "Setting resample piece information " << myRank << " of " << numProcs << std::endl;
  vtkStreamingDemandDrivenPipeline *resample_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(resample_algorithm->GetExecutive());
  // no piece info set yet, assumes info is not piece dependent
  if (ImageTest) {
    // We can't compute the extents of the image until we know the bounds of the data, 
    // so update reader first.
    vtkStreamingDemandDrivenPipeline *reader_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
    reader_sddp->UpdateDataObject();
    reader_sddp->UpdateInformation();
    reader_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
    reader_sddp->Update();
    //
    resample_sddp->UpdateDataObject();
    resample_sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  }
  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();
  //
  // Update SPH Pipeline
  //
  resample_sddp->UpdateInformation();
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
      double tol_min = std::fabs(vmin/1000.0);
      double tol_max = std::fabs(vmax/1000.0);
      if (std::fabs(vmin-scalar_range_global[0])>tol_min || std::fabs(vmax-scalar_range_global[1])>tol_max) {
        ok = false;
        std::cout << "//--------------------------------------------------------------" << std::endl;
        std::cout << "min/max check failed " << std::endl;
        std::cout << "expected {" << vmin << ',' << vmax << "}" << std::endl;
        std::cout << "got {" << scalar_range_global[0] << ',' << scalar_range_global[1] << "}" << std::endl;
        std::cout << "err {" << std::abs(vmin-scalar_range_global[0]) << ',' << std::abs(vmax-scalar_range_global[1]) << "}" << std::endl;
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
      const char *array_name = "FLUID_pressure";
//      const char *array_name = "ProcessId";
      vtkDataArray *da = append->GetOutput()->GetPointData()->GetArray(array_name);
      mapper->SetScalarRange(da->GetRange());
      mapper->SelectColorArray(array_name);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      ren->GetActiveCamera()->SetPosition(  -0.160102, 0.106656, -0.212755);
      ren->GetActiveCamera()->SetFocalPoint(-0.160102, 0.106656,  0.023384);
      ren->GetActiveCamera()->SetViewUp(0,1,0);
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

