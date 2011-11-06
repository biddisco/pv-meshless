// -D C:\share -F sph-test.h5 -N 10000 -I -R
// mpishim : C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\Remote Debugger\x64\mpishim100.exe
// mpiexec : C:\Program Files\MPICH2\bin\mpiexec.exe 
// mpiargs : -localonly -n 2 -env PATH C:\cmakebuild\pv-meshless\bin\debug;c:\bin
//
// windows mpi command line
// cd D:\cmakebuild\csviz\bin\relwithdebinfo
// mpiexec --localonly -n 4 TestParticlePartition -D C:\scratch -F sph-test.h5 -N 10000 -I -R


#ifdef _WIN32
  #include <windows.h>
#else 
  #include <sys/time.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
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
#include <algorithm>
//
#include "zoltan.h"
#undef min
#undef max

//----------------------------------------------------------------------------
std::string usage = "\n"\
"\t-D path to use for temp h5part file \n" \
"\t-F name to use for temp h5part file \n" \
"\t-N Generate exactly N points per process (no cube structure) \n" \
"\t-R Render. Displays points using vtk renderwindow \n" \
"\t-I Interactive (waits until user closes render window if -R selected) \n" \
"\t-X delete h5part file on completion \n";

class Random {
  public:
    unsigned int __seed;
    Random(int seed) {
      __seed = seed;
    }

    unsigned int getseed() {
      return __seed;
    }
    void setseed(int seed) {
        __seed = seed;
    }
    double nextNumber()
    {
      __seed = (__seed*9301+49297) % 233280;
      return __seed / 233280.0;
    }
};
//----------------------------------------------------------------------------
#ifdef _WIN32
void known_seed()
{
  srand(1);
}
#else
void known_seed()
{
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
//    std::cout << r1 << " " << r2 << std::endl;
    z = 2.0 * r1 - 1.0;
    t = 2.0 * M_PI * r2;
    w = radius * sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    X[3*i+0] = x;
    X[3*i+1] = y;
    X[3*i+2] = z*radius;
  }
//  std::cout << "min max " << rmin << " " << rmax << std::endl;
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
  if (numProcs!=4) {
    std::cout << "This test must be run using 4 processors" << std::endl;
    controller->Finalize();
    return 1;
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

  vtkTypeInt64 numPoints = 1000;

  char *number = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-N", argc, argv, "DUMMY_ENV_VAR", "");
  if (std::string(number)!=std::string("")) {
    vtkstd::stringstream temp;
    temp << number;
    temp >> numPoints;
    if (myRank<2) {
      std::cout << "Process Id : " << myRank << " Requested Particles : " << numPoints << std::endl;
    }
  }

  //--------------------------------------------------------------
  // allocate points + arrays
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
  

  known_seed();
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
  // Create partitioning filter
  //--------------------------------------------------------------
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  partitioner->SetInput(Sprites);
  partitioner->SetIdChannelArray("PointIds");
  partitioner->SetGhostCellOverlap(KernelMaximum);
  partitioner->SetController(controller);

  //--------------------------------------------------------------
  // Add process Id's
  //--------------------------------------------------------------
  vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
  processId->SetInputConnection(partitioner->GetOutputPort());
  processId->SetController(controller);

  //--------------------------------------------------------------
  // Create parallel writer on all processes
  //--------------------------------------------------------------
  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileModeToWrite();
  writer->SetFileName(fullname);
  writer->SetInputConnection(processId->GetOutputPort());
  writer->SetCollectiveIO(true);
  writer->SetDisableInformationGather(0);
  writer->SetVectorsWithStridedWrite(0);
  writer->SetController(controller);
  writer->SetTimeStep(0);
  writer->SetTimeValue(0.0);

  // Randomly give some processes zero points to improve test coverage
  random_seed();
  if (0 && numProcs>1 && rand()%2==1) {
    numPoints = 0;
    Sprites = vtkSmartPointer<vtkPolyData>::New();
    writer->SetInput(Sprites);
  }
  
  controller->Barrier();
  if (myRank==0) {
    std::cout << "Process Id : " << myRank << " Generated N Points : " << numPoints << std::endl;
  }

  //--------------------------------------------------------------
  // Write in parallel
  //
  // To get parallel writing correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------
  std::cout << "Setting piece information " << myRank << " of " << numProcs << std::endl;
  vtkStreamingDemandDrivenPipeline *sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(writer->GetExecutive());
  vtkSmartPointer<vtkInformation> execInfo = sddp->GetOutputInformation(0);
  // no piece info set yet, assumes info is not piece dependent
  sddp->UpdateInformation();
  sddp->SetUpdateExtent(0, myRank, numProcs, 0);
  // this calls writer->Write(); using correct piece information
  sddp->Update();

  // 
  // make sure they have all finished writing before going on to the read part
  //
  writer->CloseFile();
  controller->Barrier();

  //--------------------------------------------------------------
  // Visualization : Read back all particles on process zero
  //--------------------------------------------------------------
  if (myRank == 0) {
    std::cout << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    std::cout << "Process Id : " << myRank << " Expected " << static_cast<vtkTypeInt64>(numPoints*numProcs) << std::endl;

    // Read the file we just wrote on N processes
    vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
    // we want to read all the particles on this node, so don't use MPI/Parallel
    reader->SetController(NULL);
    reader->SetFileName(fullname);
    reader->SetGenerateVertexCells(1);
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
      // Display something to see if we got a good output
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInputConnection(reader->GetOutputPort());
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointFieldData();
      mapper->SetUseLookupTableScalarRange(0);
      mapper->SetScalarRange(0,numProcs-1);
      mapper->SetInterpolateScalarsBeforeMapping(0);
      mapper->SelectColorArray("ProcessId");
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor2 = vtkSmartPointer<vtkActor>::New();
      mapper2->SetInputConnection(reader->GetOutputPort());
      mapper2->SetColorModeToMapScalars();
      mapper2->SetScalarModeToUsePointFieldData();
      mapper2->SetUseLookupTableScalarRange(0);
      mapper2->SetScalarRange(0,1);
      mapper2->SetInterpolateScalarsBeforeMapping(0);
      mapper2->SelectColorArray("vtkGhostLevels");
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
        vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>          bactor = vtkSmartPointer<vtkActor>::New();
        bmapper->SetInputConnection(boxsource->GetOutputPort());
        bactor->SetMapper(bmapper);
        ren->AddActor(bactor);
      }
      
      ren->GetActiveCamera()->SetPosition(0,4*radius,0);
      ren->GetActiveCamera()->SetFocalPoint(0,0,0);
      ren->GetActiveCamera()->SetViewUp(0,0,-1);
      ren->ResetCamera();
      std::cout << "Process Id : " << myRank << " About to Render" << std::endl;
      renWindow->Render();

      retVal = vtkRegressionTester::Test(argc, argv, renWindow, 10);
      if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
        iren->Start();
      }
      std::cout << "Process Id : " << myRank << " Rendered" << std::endl;
    }
  }

  if (myRank==0 && test->IsFlagSpecified("-X")) {
   std::cout << "Process Id : " << myRank << " About to Delete file" << std::endl;
   vtksys::SystemTools::RemoveFile(fullname);
  }

  controller->Barrier();
  controller->Finalize();

  delete []fullname;
  delete []filename;

  return !retVal;
}
//----------------------------------------------------------------------------
