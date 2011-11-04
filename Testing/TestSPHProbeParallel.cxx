// appcmd  : $(TargetPath)
// appargs : -D C:\share -F sph-test.h5 -N 10000 -I -R
// mpishim : C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\Remote Debugger\x64\mpishim100.exe
// mpiexec : C:\Program Files\MPICH2\bin\mpiexec.exe 
// mpiargs : -localonly -n 2 -env PATH C:\cmakebuild\pv-meshless\bin\debug;c:\bin
//
// windows mpi command line
// cd D:\cmakebuild\csviz\bin\relwithdebinfo
// mpiexec --localonly -n 4 TestSPHProbeParallel -D C:\share -F sph-test.h5 -N 10000 -I -R


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
#include "vtkSPHProbeFilter.h"
#include "vtkSPHManager.h"
#include "vtkParticleIdFilter.h"
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

  //--------------------------------------------------------------
  // File names/dirs
  //--------------------------------------------------------------
  char *n1 = vtkTestUtilities::GetArgOrEnvOrDefault("-neighbours", argc, argv, "", "");
  int maxN;
  bool fixNeighbours = false;
  if (std::string(n1).size()) {
    vtkstd::stringstream temp(n1);
    temp >> maxN;
    if (myRank==0) {
      std::cout << "Fixed Neighbours      {" << maxN << "}" << std::endl;
    }
    fixNeighbours = true;
  }
  delete []n1;

  //--------------------------------------------------------------
  // Testing params
  //--------------------------------------------------------------
  char *n2 = vtkTestUtilities::GetArgOrEnvOrDefault("-particlesize", argc, argv, "", "");
  double p_size;
  bool fixRadius = false;
  if (std::string(n2).size()) {
    vtkstd::stringstream temp(n2);
    temp >> p_size;
    if (myRank==0) {
      std::cout << "Particle Size         {" << p_size << "}" << std::endl;
    }
    fixRadius = true;
  }
  delete []n2;

  char *m1 = vtkTestUtilities::GetArgOrEnvOrDefault("-scalar", argc, argv, "", "");
  std::string scalarname;
  if (std::string(m1).size()) {
    vtkstd::stringstream temp(m1);
    temp >> scalarname;
    if (myRank==0) {
      std::cout << "Testing Scalar array  {" << scalarname << "}" << std::endl;
    }
  }
  delete []m1;

  char *n4 = vtkTestUtilities::GetArgOrEnvOrDefault("-ghost_region", argc, argv, "", "");
  double ghost = 0.0;
  if (std::string(n4).size()) {
    vtkstd::stringstream temp(n4);
    temp >> ghost;
    if (myRank==0) {
      std::cout << "Ghost Region {" << ghost << "}" << std::endl;
    }
  }
  delete []n4;

  //--------------------------------------------------------------
  // Valid test results
  //--------------------------------------------------------------
  char *v1 = vtkTestUtilities::GetArgOrEnvOrDefault("-value_range", argc, argv, "", "");
  double vmin=0,vmax=0;
  if (std::string(v1).size()) {
    vtkstd::stringstream temp(v1);
    temp >> vmin >> vmax;
    if (myRank==0) {
      std::cout << "Test Valid : value_range {" << vmin << ',' << vmax << "}" << std::endl;
    }
  }
  delete []v1;

  char *v2 = vtkTestUtilities::GetArgOrEnvOrDefault("-peak_position", argc, argv, "", "");
  double vpos[3];
  if (std::string(v2).size()) {
    vtkstd::stringstream temp(v2);
    temp >> vpos[0] >> vpos[1] >> vpos[2];
    if (myRank==0) {
      std::cout << "Test Valid : scalar_peak_position {" << vpos[0] << ',' << vpos[1] << "," << vpos[2] << "}" << std::endl;
    }
  }
  delete []v2;

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

  vtkSmartPointer<vtkAlgorithm> algorithm = reader; 
#ifdef PV_MESHLESS_TRILINOS
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(reader->GetOutputPort());
  partitioner->SetGhostCellOverlap(ghost);
  partitioner->SetController(controller);
  algorithm = partitioner;
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
    sphManager->SetDefaultParticleSideLength(p_size);
    sphManager->SetKernelTypeToCubicSpline();
    sphManager->SetKernelDimension(3);
    sphManager->SetHCoefficient(1.5);
    sphManager->SetMaximumSearchRadius(p_size*1.5*5.0); 
  }

  vtkSmartPointer<vtkSPHProbeFilter> sphProbe = vtkSmartPointer<vtkSPHProbeFilter>::New();
  sphProbe->SetInputConnection(algorithm->GetOutputPort());
  sphProbe->SetProbeConnection(algorithm->GetOutputPort());
  sphProbe->SetSPHManager(sphManager);
  sphProbe->SetMassScalars("mass");
  sphProbe->SetComputeDensityFromNeighbourVolume(1);

  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();
  //
  //--------------------------------------------------------------
  // Update in parallel:
  // To get parallel operation correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------
  std::cout << "Setting piece information " << myRank << " of " << numProcs << std::endl;
  vtkStreamingDemandDrivenPipeline *sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(sphProbe->GetExecutive());
  vtkSmartPointer<vtkInformation> execInfo = sddp->GetOutputInformation(0);
  // no piece info set yet, assumes info is not piece dependent
  sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  sddp->UpdateInformation();
  sddp->Update();

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
  // Fetch smoothed output
  //--------------------------------------------------------------
  vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
  polys->ShallowCopy(sphProbe->GetOutput());
  vtkPointData *pd = polys->GetPointData();
  vtkDataArray *scalar_array = pd->GetArray(scalarname.c_str());
  
  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  bool ok = true;
  std::cout.precision(15);
  
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
    polys->GetPoint(index,scalar_pos);
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
 
  retVal = ok;
  
  delete []fullname;
  delete []filename;
  controller->Finalize();

  return !retVal;
}
//----------------------------------------------------------------------------
