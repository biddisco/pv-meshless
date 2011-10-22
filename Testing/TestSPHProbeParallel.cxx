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
#include "Testing/Cxx/vtkTestUtilities.h"
#include "Testing/Cxx/vtkRegressionTestImage.h"
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
#include "vtkBoundingBox.h"
#include "vtkOutlineSource.h"
#include "vtkParticlePartitionFilter.h"
#include "vtkProcessIdScalars.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkSPHProbeFilter.h"
#include "vtkSPHManager.h"
#include "vtkParticleIdFilter.h"
//
#include <vtksys/SystemTools.hxx>
#include <sstream>
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#define _USE_MATH_DEFINES
#include <math.h>

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
  if (numProcs!=4) {
//    std::cout << "This test must be run using 4 processors" << std::endl;
//    controller->Finalize();
//    return 1;
  }
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

  char *v1 = vtkTestUtilities::GetArgOrEnvOrDefault("-range", argc, argv, "", "");
  double vmin=0,vmax=0;
  if (std::string(v1).size() && (myRank==0)) {
    vtkstd::stringstream temp(v1);
    temp >> vmin >> vmax;
    std::cout << "Test Valid : range    {" << vmin << ',' << vmax << "}" << std::endl;
  }
  delete []v1;

  char *v2 = vtkTestUtilities::GetArgOrEnvOrDefault("-index", argc, argv, "", "");
  vtkIdType vindex;
  if (std::string(v2).size() && (myRank==0)) {
    vtkstd::stringstream temp(v2);
    temp >> vindex;
    std::cout << "Test Valid : index     " << vindex << std::endl;
  }
  delete []v2;

  char *v3 = vtkTestUtilities::GetArgOrEnvOrDefault("-coord", argc, argv, "", "");
  double vpos[3];
  if (std::string(v3).size() && (myRank==0)) {
    vtkstd::stringstream temp(v3);
    temp >> vpos[0] >> vpos[1] >> vpos[2];
    std::cout << "Test Valid : position {" << vpos[0] << ',' << vpos[1] << "," << vpos[2] << "}" << std::endl;
  }
  delete []v3;

  vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
  reader->SetFileName(fullname);
  reader->SetController(controller);
  reader->SetGenerateVertexCells(1);

  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(reader->GetOutputPort());
  partitioner->SetGhostCellOverlap(0.01);
  partitioner->SetController(controller);

  vtkSmartPointer<vtkSPHManager> sphManager = vtkSmartPointer<vtkSPHManager>::New();
  sphManager->SetInterpolationMethodToKernel();
  sphManager->SetHCoefficient(1.5);
  sphManager->SetDefaultParticleSideLength(0.005);
  sphManager->SetDefaultDensity(1000.0);
  sphManager->SetKernelTypeToGaussian();
  sphManager->SetKernelDimension(3);
  sphManager->SetLimitSearchByNeighbourCount(1);
  sphManager->SetMaximumNeighbours(32);
  sphManager->SetMaximumRadius(0.01); 

  vtkSmartPointer<vtkSPHProbeFilter> sphProbe = vtkSmartPointer<vtkSPHProbeFilter>::New();
  sphProbe->SetInputConnection(partitioner->GetOutputPort());
  sphProbe->SetProbeConnection(partitioner->GetOutputPort());
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
  // processes 1-N doing nothing for now
  //--------------------------------------------------------------
  double range[2];
  vtkIdType index = -1;
  double pos[3];
  if (myRank == 0) {
    vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
    polys->ShallowCopy(sphProbe->GetOutput());
    vtkDataArray *da = polys->GetPointData()->GetArray("SmoothedDensity");
    da->GetRange(range);
    std::cout << "Min and Max of SmoothedDensity are {" << range[0] << "," << range[1] << "}" << std::endl;
    // 0.567804, 0.586049, 0.559157 
    for (vtkIdType i=0; i<da->GetNumberOfTuples(); i++) {
      if (da->GetTuple1(i)==range[1]) {
        index = i;
      }
    }
    if (index!=-1) {
      std::cout << "Index of peak SmoothedDensity particle is " << index << std::endl;
      polys->GetPoint(index,pos);
      std::cout << "Position of peak SmoothedDensity particle is {" << pos[0] << "," << pos[1] << "," << pos[2] << "}" << std::endl;
    }
    bool ok = true;
    if (index!=vindex) {
      ok = false;
      std::cout << "index check failed " << std::endl;
    }
    if (std::abs(vmin-range[0])>1E-3 || std::abs(vmax-range[1])>1E-3) {
      ok = false;
      std::cout << "min/max check failed " << std::endl;
      std::cout << "expected {" << vmin << ',' << vmax << "}" << std::endl;
      std::cout << "got {" << range[0] << ',' << range[1] << "}" << std::endl;
      std::cout << "err {" << std::abs(vmin-range[0]) << ',' << std::abs(vmax-range[1]) << "}" << std::endl;
    }
    if (std::abs(vpos[0]-pos[0])>1E-5 ||
        std::abs(vpos[1]-pos[1])>1E-5 ||
        std::abs(vpos[2]-pos[2])>1E-5)
    {
      ok = false;
      std::cout << "position check failed " << std::endl;
    }
    if (ok) {
      std::cout << "All checks passed " << std::endl;
    }
    retVal = ok;
  }

  delete []fullname;
  delete []filename;
  controller->Finalize();

  return !retVal;
}
//----------------------------------------------------------------------------
