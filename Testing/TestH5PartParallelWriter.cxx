// TestHDF5PartWriter -I -D D:/ -F test.h5

// windows mpi command line
// cd D:\cmakebuild\csviz\bin\relwithdebinfo
// mpiexec --localonly -n 4 TestH5PartParallelWriter -I -D D:\ -F sph-test.h5 -R -X -N 1000000

// horus slurm command line for auto allocation of nodes
// mpirun -prot -srun -N 6 -n 6 TestH5PartParallelWriter -R -D . -N 1000000

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
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"
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

#include <vtksys/SystemTools.hxx>
#include <sstream>

#include <mpi.h>
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
// call back to set the iso surface value.
void SetStuffRMI(void *localArg, void* vtkNotUsed(remoteArg), 
                    int vtkNotUsed(remoteArgLen), int vtkNotUsed(id))
{ 
  ParallelRMIArgs_tmp* args = (ParallelRMIArgs_tmp*)localArg;
  vtkMultiProcessController* contrl = args->Controller;
}
//----------------------------------------------------------------------------
#define ROWS 50
#define POINTS_PER_PROCESS ROWS*ROWS*ROWS
//----------------------------------------------------------------------------
// This will be called by all processes
void MyMain( vtkMultiProcessController *controller, void *arg )
{
  // Obtain the id of the running process and the total
  // number of processes
  vtkTypeInt64 myId = controller->GetLocalProcessId();
  vtkTypeInt64 numProcs = controller->GetNumberOfProcesses();

  if (myId==0) {
    std::cout << usage.c_str() << std::endl;
  }
  controller->Barrier();

  // Setup testing utilities/args etc
  ParallelArgs_tmp* args = reinterpret_cast<ParallelArgs_tmp*>(arg);
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<args->argc; c++ ) {
    test->AddArgument(args->argv[c]);
  }
  // Get test filename etc
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-F", args->argc, args->argv, "DUMMY_ENV_VAR", "temp.h5");
  char* fullname = vtkTestUtilities::ExpandDataFileName(args->argc, args->argv, filename);
  if (myId==0) {
    std::cout << "Process Id : " << myId << " FileName : " << fullname << std::endl;
  }

  vtkTypeInt64 numPoints = POINTS_PER_PROCESS;
  double       rows      = ROWS;

  char *number = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-N", args->argc, args->argv, "DUMMY_ENV_VAR", "");
  if (std::string(number)!=std::string("")) {
    std::stringstream temp;
    temp << number;
    temp >> numPoints;
    rows = floor(pow(numPoints,1.0/3.0)+0.5);
    numPoints = static_cast<vtkTypeInt64>(pow(rows,3));
    if (myId==0) {
      std::cout << "Process Id : " << myId << " Requested Particles : " << numPoints << std::endl;
    }
  }
  //
  vtkSmartPointer<vtkPolyData>  Sprites = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>   verts = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkDoubleArray> sizes = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray>      Ids = vtkSmartPointer<vtkIntArray>::New();
  //
  points->SetNumberOfPoints(numPoints);
  //
  verts->Allocate(numPoints,numPoints);
  Sprites->SetPoints(points);
  Sprites->SetVerts(verts);
  //
  sizes->SetNumberOfTuples(numPoints);
  sizes->SetNumberOfComponents(1);
  sizes->SetName("PointSizes");
  Sprites->GetPointData()->AddArray(sizes);
  //
  Ids->SetNumberOfTuples(numPoints);
  Ids->SetNumberOfComponents(1);
  Ids->SetName("PointIds");
  Sprites->GetPointData()->AddArray(Ids);  
  //
  double radius  = 0.0001;
  radius  = 0.08;
  double spacing = radius*2.0;
  double offset  = myId*spacing*rows;
  for (vtkTypeInt64 z=0; z<rows; z++) {
    for (vtkTypeInt64 y=0; y<rows; y++) {
      for (vtkTypeInt64 x=0; x<rows; x++) {
        vtkIdType Id = static_cast<vtkIdType>(z*rows*rows + y*rows + x);
        points->SetPoint(Id, x*spacing, y*spacing, z*spacing + offset);
        sizes->SetValue(Id, radius);
        Ids->SetTuple1(Id, Id + myId*numPoints);
        verts->InsertNextCell(1,&Id);
      }
    }
  }

  vtkSmartPointer<vtkElevationFilter> elev = vtkSmartPointer<vtkElevationFilter>::New();
  elev->SetInputData(Sprites);
  elev->SetLowPoint(offset, 0.0, 0.0);
  elev->SetHighPoint(offset+rows*spacing, 0.0, 0.0);
  elev->Update();

  bool collective = false;
  if (test->IsFlagSpecified("-C")) {
    if (myId==0) {
     std::cout << "Process Id : " << myId << " Collective IO requested" << std::endl;
    }
   collective = true;
  }

  // Create writer
  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileModeToWrite();
  writer->SetFileName(fullname);
  writer->SetInputConnection(elev->GetOutputPort());
  writer->SetCollectiveIO(collective);
  writer->SetDisableInformationGather(1);
  writer->SetVectorsWithStridedWrite(0);


/*
  // Randomly give some processes zero points to improve test coverage
  random_seed();
  if (numProcs>1 && rand()%2==3) {
    numPoints = 0;
    Sprites = vtkSmartPointer<vtkPolyData>::New();
    writer->SetInput(Sprites);
  }
  else {
    writer->SetInputConnection(elev->GetOutputPort());
  }
*/

  controller->Barrier();
  if (myId==0) {
    std::cout << "Process Id : " << myId << " Generated N Points : " << numPoints << std::endl;
  }
  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  //if (numProcs>1 && myId==0) {
  //  char ch;  
  //  std::cin >> ch;
  //}

  writer->SetTimeStep(0);
  writer->SetTimeValue(0.5);
  writer->Write();
  // 
  // make sure they have all finished writing before going on to the read part
  //
  timer->StopTimer();
  writer->CloseFile();
  controller->Barrier();
  //
  // memory usage - Ids(int) Size(double) Elevation(float) Verts(double*3)
  double bytes = numPoints*(sizeof(int) + sizeof(double) + sizeof(float) + 3*sizeof(float));
  double MBytes = bytes/(1024*1024);
  double elapsed = timer->GetElapsedTime();
  std::cout << "Process Id : " << myId << " File Written in " << elapsed << " seconds" << std::endl;
  std::cout << "Process Id : " << myId << " IO-Speed " << MBytes/timer->GetElapsedTime() << " MB/s" << std::endl;
  //
  if (myId != 0)
    {
    // If I am not the root process
    ParallelRMIArgs_tmp args2;
    args2.Controller = controller;

    // We are not using any RMI's yet
    controller->AddRMI(SetStuffRMI, (void *)&args2, RMI_TAG);
    controller->ProcessRMIs();
    
    }
  else
    {
    std::cout << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    std::cout << "Process Id : " << myId << " Expected " << static_cast<vtkTypeInt64>(numPoints*numProcs) << std::endl;

    // Read the file we just wrote on N processes
    vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
    reader->SetFileName(fullname);
    reader->Update();
    vtkTypeInt64 ReadPoints = reader->GetOutput()->GetNumberOfPoints();
    std::cout << "Process Id : " << myId << " Read : " << ReadPoints << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    //
    // Validate the point Ids to make sure nothing went wrong in the writing
    //
    vtkIntArray *pointIds = vtkIntArray::SafeDownCast(reader->GetOutput()->GetPointData()->GetArray("PointIds"));
    bool IdsGood = true;
    vtkTypeInt64 valid = 0;
    for (vtkTypeInt64 i=0; IdsGood && pointIds && i<ReadPoints; i++) {
      if (pointIds->GetValue(i)==i) {
        valid++;
      }
      else {
        IdsGood = false;
      }
    }
    if (IdsGood && ReadPoints==numPoints*numProcs) {
      *(args->retVal) = 0;
      std::cout << " " << std::endl;
      std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * "   << std::endl;
      std::cout << " All Points read back and Validated OK "               << std::endl;
      std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * \n" << std::endl;
      //
      if (myId==0) {
        unsigned long size = vtksys::SystemTools::FileLength(fullname);
        double filesize = size/(1024.0*1024.0);
        std::cout << "Process Id : " << myId << " Total IO/Disk-Speed " << filesize/elapsed << " MB/s" << std::endl;
      }
    }
    else {
      *(args->retVal) = 1;
      std::cout << " " << std::endl;
      std::cout << " # # # # # # # # # # # # # # # # # # # # # # # # # " << std::endl;
      std::cout << " FAIL "                                              << std::endl;
      std::cout << " Valid Ids "                                << valid << std::endl;
      std::cout << " # # # # # # # # # # # # # # # # # # # # # # # # # " << std::endl;
      std::cout << " " << std::endl;
    }

    bool doRender = false;
    if (test->IsFlagSpecified("-R")) {
     std::cout << "Process Id : " << myId << " Rendering" << std::endl;
     doRender = true;
    }

    if (doRender) {
      // Generate vertices from the points
      vtkSmartPointer<vtkMaskPoints> verts = vtkSmartPointer<vtkMaskPoints>::New();
      verts->SetGenerateVertices(1);
      verts->SetOnRatio(1);
      verts->SetMaximumNumberOfPoints(numPoints*numProcs);
      verts->SetInputConnection(reader->GetOutputPort());
      verts->Update();

      // Display something to see if we got a good output
      vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
      polys->ShallowCopy(verts->GetOutput());
      polys->GetPointData()->SetScalars(polys->GetPointData()->GetArray("Elevation"));
      std::cout << "Process Id : " << myId << " Created vertices : " << polys->GetNumberOfPoints() << std::endl;
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInputData(polys);
      mapper->SelectColorArray("Elevation");
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
    
      std::cout << "Process Id : " << myId << " About to Render" << std::endl;
      renWindow->Render();

      *(args->retVal) = 
        vtkRegressionTester::Test(args->argc, args->argv, renWindow, 10);

      if ( *(args->retVal) == vtkRegressionTester::DO_INTERACTOR)
        {
        iren->Start();
        }
      std::cout << "Process Id : " << myId << " Rendered" << std::endl;
    }

    // Tell the other processors to stop processing RMIs.
    for (int i = 1; i < numProcs; ++i)
    {
      controller->TriggerRMI(i, vtkMultiProcessController::BREAK_RMI_TAG); 
    }
  }

//  writer = NULL;

  if (myId==0 && test->IsFlagSpecified("-X")) {
   std::cout << "Process Id : " << myId << " About to Delete file" << std::endl;
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
