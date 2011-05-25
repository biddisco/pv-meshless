// -D C:\share -F sph-test.h5 -N 10000 -I -R
// mpishim : C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\Remote Debugger\x64\mpishim100.exe
// mpiexec : C:\Program Files\MPICH2\bin\mpiexec.exe 
// mpiargs : -localonly -n 2 -env PATH C:\cmakebuild\pv-meshless\bin\debug;c:\bin
//
// windows mpi command line
// cd D:\cmakebuild\csviz\bin\relwithdebinfo
// mpiexec --localonly -n 4 TestParticlePartition -D C:\share -F sph-test.h5 -N 10000 -I -R


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
#include "zoltan.h"

//----------------------------------------------------------------------------
std::string usage = "\n"\
"\t-D path to use for temp h5part file \n" \
"\t-F name to use for temp h5part file \n" \
"\t-C collective IO for MPI/HDF write (default independent IO) \n" \
"\t-N Generate exactly N points per process (no cube structure) \n" \
"\t-R Render. Displays points using vtk renderwindow \n" \
"\t-I Interactive (waits until user closes render window if -R selected) \n" \
"\t-X delete h5part file on completion \n" \
"\t\n" \
"\tWIN32 mpi example   : mpiexec -n 4 -localonly TestH5PartParallelWriter.exe -D D:\\ -F sph-test.h5 -X\n"\
"\tLinux mpich example : /project/csvis/biddisco/build/mpich2-1.0.7/bin/mpiexec -n 2 ./TestH5PartParallelWriter -D /project/csvis/hdf-scratch -C -N 1000000 -X\n" \
"\tHorus Slurm example : mpirun -prot -srun mpirun -e DISPLAY=horus11:0 -e MPI_IB_CARD_ORDER=0:0 -e MPI_IB_MTU=1024 -e MPI_IC_ORDER=ibv:vapi:udapl:itapi:TCP -srun -N 10 -n 20 bin/TestH5PartParallelWriter -D /project/csvis/hdf-scratch -C -N 100000 -R -I -X\n";

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
// Just pick a tag which is available
static const int RMI_TAG=300; 
//----------------------------------------------------------------------------
struct ParallelArgs_tmp
{
  int* retVal;
  int    argc;
  char** argv;
};
//----------------------------------------------------------------------------
struct ParallelRMIArgs_tmp
{
  vtkMultiProcessController* Controller;
};
//----------------------------------------------------------------------------
void SetStuffRMI(void *localArg, void* vtkNotUsed(remoteArg), 
                    int vtkNotUsed(remoteArgLen), int vtkNotUsed(id))
{ 
  ParallelRMIArgs_tmp* args = (ParallelRMIArgs_tmp*)localArg;
  vtkMultiProcessController* contrl = args->Controller;
}
//----------------------------------------------------------------------------
void SpherePoints(int n, float radius, float X[]) {
  double x, y, z, w, t;
  for(int i=0; i< n; i++ ) {
    #ifdef WIN32
     double r1 = double(rand())/RAND_MAX;
     double r2 = double(rand())/RAND_MAX;
    #else
     double r1 = drand48();
     double r2 = drand48();
    #endif
    z = 2.0 * r1 - 1.0;
    t = 2.0 * M_PI * r2;
    w = radius * sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    X[3*i+0] = x;
    X[3*i+1] = y;
    X[3*i+2] = z*radius;
  }
}
//----------------------------------------------------------------------------
// This will be called by all processes
void MyMain( vtkMultiProcessController *controller, void *arg )
{
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
  ParallelArgs_tmp* args = reinterpret_cast<ParallelArgs_tmp*>(arg);
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<args->argc; c++ ) {
    test->AddArgument(args->argv[c]);
  }
  // Get test filename etc
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-F", args->argc, args->argv, "DUMMY_ENV_VAR", "temp.h5");
  char* fullname = vtkTestUtilities::ExpandDataFileName(args->argc, args->argv, filename);
  if (myRank==0) {
    std::cout << "Process Id : " << myRank << " FileName : " << fullname << std::endl;
  }

  vtkTypeInt64 numPoints = 100;

  char *number = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-N", args->argc, args->argv, "DUMMY_ENV_VAR", "");
  if (std::string(number)!=std::string("")) {
    vtkstd::stringstream temp;
    temp << number;
    temp >> numPoints;
    if (myRank<2) {
      std::cout << "Process Id : " << myRank << " Requested Particles : " << numPoints << std::endl;
    }
  }

  //--------------------------------------------------------------
  // allocate scalar arrays
  //--------------------------------------------------------------
  vtkSmartPointer<vtkPolyData>  Sprites = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>   verts = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdTypeArray>   Ids = vtkSmartPointer<vtkIdTypeArray>::New();
  vtkSmartPointer<vtkIntArray>    Ranks = vtkSmartPointer<vtkIntArray>::New();
  //
  points->SetNumberOfPoints(numPoints);
  //
  verts->Allocate(numPoints,numPoints);
  Sprites->SetPoints(points);
  Sprites->SetVerts(verts);
  //
  Ids->SetNumberOfTuples(numPoints);
  Ids->SetNumberOfComponents(1);
  Ids->SetName("PointIds");
  Sprites->GetPointData()->AddArray(Ids);  
  //
  Ranks->SetNumberOfTuples(numPoints);
  Ranks->SetNumberOfComponents(1);
  Ranks->SetName("Rank");
  Sprites->GetPointData()->AddArray(Ranks);  
  //
  //--------------------------------------------------------------
  // Create default scalar arrays
  //--------------------------------------------------------------
  double radius  = 500.0;
  const double a = 0.9;
  SpherePoints(numPoints, radius*(1.5+myRank)/(numProcs+0.5), vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0));
  for (vtkIdType Id=0; Id<numPoints; Id++) {
    Ids->SetTuple1(Id, Id + myRank*numPoints);
    Ranks->SetTuple1(Id, myRank);
    verts->InsertNextCell(1,&Id);
  }

  //--------------------------------------------------------------
  // Make up some max kernel size for testing
  //--------------------------------------------------------------
  double KernelMaximum  = radius*0.1;

  //--------------------------------------------------------------
  // Add colour by elevation
  //--------------------------------------------------------------
  vtkSmartPointer<vtkElevationFilter> elev = vtkSmartPointer<vtkElevationFilter>::New();
  elev->SetInput(Sprites);
  elev->SetLowPoint(-radius, -radius, -radius);
  elev->SetHighPoint(radius,  radius,  radius);
  elev->Update();

  bool collective = false;
  if (test->IsFlagSpecified("-C")) {
    if (myRank==0) {
     std::cout << "Process Id : " << myRank << " Collective IO requested" << std::endl;
    }
   collective = true;
  }

  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInputConnection(elev->GetOutputPort());
  partitioner->SetIdChannelArray("PointIds");
  partitioner->SetGhostCellOverlap(KernelMaximum);
  partitioner->Update();

  vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
  processId->SetInputConnection(partitioner->GetOutputPort());
  processId->Update();

  //*****************************************************************
  // Visualize the mesh partitioning
  //*****************************************************************

  //--------------------------------------------------------------
  // Create writer on all processes
  //--------------------------------------------------------------
  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileModeToWrite();
  writer->SetFileName(fullname);
  writer->SetInputConnection(processId->GetOutputPort());
  writer->SetCollectiveIO(collective);
  writer->SetDisableInformationGather(1);
  writer->SetVectorsWithStridedWrite(0);

//
  
  // Randomly give some processes zero points to improve test coverage
  random_seed();
  if (0 && numProcs>1 && rand()%2==3) {
    numPoints = 0;
    Sprites = vtkSmartPointer<vtkPolyData>::New();
    writer->SetInput(Sprites);
  }
  
  controller->Barrier();
  if (myRank==0) {
    std::cout << "Process Id : " << myRank << " Generated N Points : " << numPoints << std::endl;
  }
  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  //if (numProcs>1 && myRank==0) {
  //  char ch;  
  //  std::cin >> ch;
  //}

  //--------------------------------------------------------------
  // Write in parallel
  //--------------------------------------------------------------
  writer->SetTimeStep(0);
  writer->SetTimeValue(0.5);
  writer->Write();

  // 
  // make sure they have all finished writing before going on to the read part
  //
  timer->StopTimer();
  writer->CloseFile();
  controller->Barrier();

  // memory usage - Ids(int) Size(double) Elevation(float) Verts(double*3)
  double bytes = numPoints*(sizeof(int) + sizeof(double) + sizeof(float) + 3*sizeof(float));
  double MBytes = bytes/(1024*1024);
  double elapsed = timer->GetElapsedTime();
  std::cout << "Process Id : " << myRank << " File Written in " << elapsed << " seconds" << std::endl;
  std::cout << "Process Id : " << myRank << " IO-Speed " << MBytes/timer->GetElapsedTime() << " MB/s" << std::endl;
  //
  //--------------------------------------------------------------
  // processes 1-N doing nothing for now
  //--------------------------------------------------------------
  if (myRank != 0)
    {
    // If I am not the root process
    ParallelRMIArgs_tmp args2;
    args2.Controller = controller;

    // We are not using any RMI's yet
    controller->AddRMI(SetStuffRMI, (void *)&args2, RMI_TAG);
    controller->ProcessRMIs();
    
    }
  //--------------------------------------------------------------
  // Read back all particles on process zero
  //--------------------------------------------------------------
  else
    {
    std::cout << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    std::cout << "Process Id : " << myRank << " Expected " << static_cast<vtkTypeInt64>(numPoints*numProcs) << std::endl;

    // Read the file we just wrote on N processes
    vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
    // we want to read all the particles on this node, so don't use MPI/Parallel
    reader->SetController(NULL);
    reader->SetFileName(fullname);
    reader->Update();
    vtkTypeInt64 ReadPoints = reader->GetOutput()->GetNumberOfPoints();
    std::cout << "Process Id : " << myRank << " Read : " << ReadPoints << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;

    bool doRender = false;
    if (test->IsFlagSpecified("-R")) {
     std::cout << "Process Id : " << myRank << " Rendering" << std::endl;
     doRender = true;
    }

    if (doRender) {
      // Generate vertices from the points
      vtkSmartPointer<vtkMaskPoints> verts = vtkSmartPointer<vtkMaskPoints>::New();
      verts->SetGenerateVertices(1);
      verts->SetOnRatio(1);
      verts->SetMaximumNumberOfPoints(2*numPoints*numProcs);
      verts->SetInputConnection(reader->GetOutputPort());
      verts->Update();

      // Display something to see if we got a good output
      vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
      polys->ShallowCopy(verts->GetOutput());
      polys->GetPointData()->SetScalars(polys->GetPointData()->GetArray("ProcessId"));
      std::cout << "Process Id : " << myRank << " Created vertices : " << polys->GetNumberOfPoints() << std::endl;
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInput(polys);
      mapper->SetScalarRange(0,numProcs-1);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      vtkSmartPointer<vtkPolyData> polys2 = vtkSmartPointer<vtkPolyData>::New();
      polys2->ShallowCopy(verts->GetOutput());
      polys2->GetPointData()->SetScalars(polys->GetPointData()->GetArray("GhostLevels"));
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor2 = vtkSmartPointer<vtkActor>::New();
      mapper2->SetInput(polys2);
      mapper2->SetScalarRange(0,1);
      mapper2->SetScalarModeToUsePointData();
      actor2->SetMapper(mapper2);
      actor2->GetProperty()->SetPointSize(2);
      actor2->SetPosition(2.0*radius, 0.0, 0.0);
      ren->AddActor(actor2);

      //
      // Display boxes for each partition
      //
      for (int i=0; i<numProcs; i++) {
        vtkBoundingBox *box = partitioner->GetPartitionBoundingBox(i);
        double bounds[6];
        box->GetBounds(bounds);
        vtkSmartPointer<vtkOutlineSource> boxsource = vtkSmartPointer<vtkOutlineSource>::New();
        boxsource->SetBounds(bounds);
        vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
        mapper->SetInputConnection(boxsource->GetOutputPort());
        actor->SetMapper(mapper);
        ren->AddActor(actor);
      }

      std::cout << "Process Id : " << myRank << " About to Render" << std::endl;
      renWindow->Render();

      *(args->retVal) = 
        vtkRegressionTester::Test(args->argc, args->argv, renWindow, 10);

      if ( *(args->retVal) == vtkRegressionTester::DO_INTERACTOR)
        {
        iren->Start();
        }
      std::cout << "Process Id : " << myRank << " Rendered" << std::endl;
    }

    // Tell the other processors to stop processing RMIs.
    for (int i = 1; i < numProcs; ++i)
    {
      controller->TriggerRMI(i, vtkMultiProcessController::BREAK_RMI_TAG); 
    }
  }

//  writer = NULL;

  if (myRank==0 && test->IsFlagSpecified("-X")) {
   std::cout << "Process Id : " << myRank << " About to Delete file" << std::endl;
   vtksys::SystemTools::RemoveFile(fullname);
  }
  delete []fullname;
  delete []filename;
}
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  // This is here to avoid false leak messages from vtkDebugLeaks when
  // using mpich. It appears that the root process which spawns all the
  // main processes waits in MPI_Init() and calls exit() when
  // the others are done, causing apparent memory leaks for any objects
  // created before MPI_Init().
  MPI_Init(&argc, &argv);

  // Note that this will create a vtkMPIController if MPI
  // is configured, vtkThreadedController otherwise.
  vtkMPIController* controller = vtkMPIController::New();

  controller->Initialize(&argc, &argv, 1);

  vtkParallelFactory* pf = vtkParallelFactory::New();
  vtkObjectFactory::RegisterFactory(pf);
  pf->Delete();
 
  // Added for regression test.
  // ----------------------------------------------
  int retVal = 1;
  ParallelArgs_tmp args;
  args.retVal = &retVal;
  args.argc = argc;
  args.argv = argv;
  // ----------------------------------------------

  controller->SetSingleMethod(MyMain, &args);
  controller->SingleMethodExecute();

  controller->Barrier();
  controller->Finalize();
  controller->Delete();

  return EXIT_SUCCESS; // !retVal;
}
//----------------------------------------------------------------------------
