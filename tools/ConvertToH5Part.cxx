// -D D:/data/cfxtestdata/HSLU -F Vers1_Variable_001.res -O D:/data/xdmf/test
// -D D:/data/cfxtestdata/HSLU -F Vers1_Variable_001.res -O D:/data/xdmf/scratch -I 0.002
// -D D:\data\mhd -F s15t_00000_0000-00.xdmf -O D:/data/xdmf/scratch
// -D D:/data/VectorField -F damBreakMarin.pvd -O D:/data/xdmf/scratch
// -D D:/data/cfxtestdata/HSLU -F Vers1_Variable_001.res -O D:/data/xdmf/scratch
// -D D:/data -F Test_UG_hex_elev.vtu -O D:/data/xdmf/scratch -v
// -D D:/data -F one_cell.vtu -O D:/data/xdmf/scratch -v
// -D H:/paraviewdata/data -F Elements.vtu -O D:/data/xdmf/scratch -v
// -D D:\data\cfx\SmallExample -F TRS_Default_004.res -O D:/data/xdmf/scratch
//
// jfavre data
// -D C:\data\jwalther\case004.dtvis_small -F H_helI00000.pvti -O D:/data/xdmf/scratch 
// -D C:\data\jwalther\ring -F ringRe25k_omegaI02000.pvti -O C:/data/xdmf/scratch 
//
// For Testing DSM
// Machine a=dino) -D C:\data\xdmf\other -F cav_DSM_Test.xmf -Client
// Machine b=agno) -Server
//
// cd D:\cmakebuild\pv-shared\bin\RelWithDebInfo
// mpiexec -localonly -channel mt -n 1 D:\cmakebuild\cscs-plugins\bin\RelWithDebInfo\ConvertToXdmf.exe -D D:\data\xdmf -F damBreakMarin.xmf -Client
//----------------------------------------------------------------------------
#include "vtkH5PartWriter.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkRenderer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCompositeDataIterator.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkInformation.h"
#include "vtkCompositePolyDataMapper2.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTemporalDataSetCache.h"
#include "vtkTemporalInterpolator.h"
#include "vtkTemporalDataSet.h"
#include "vtkExtentTranslator.h"

// vtk MPI
#ifdef VTK_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
// Otherwise
#include "vtkMultiProcessController.h"

#include "vtkParallelFactory.h"

// vtk Testing 
#include "vtkTesting.h"
#include "Testing/Cxx/vtkTestUtilities.h"
#include <vtksys/SystemTools.hxx>

// Readers we support
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLMultiBlockDataReader.h"
#include "vtkXMLCollectionReader.h"
#include "vtkXMLPImageDataReader.h"
#include "vtkPVDReader.h"


// CSCS
#include "vtkStreamOutputWindow.h"

// Sys
#include <sstream>
#include "vtkToolkits.h" // For VTK_USE_MPI
#ifdef VTK_USE_MPI
 #include <mpi.h>
#endif

#ifdef USE_PV_ASTRO
#include "vtkRamsesReader.h"
#endif

#define PARALLEL_PIECES

typedef std::vector<double>::size_type itype;


//----------------------------------------------------------------------------
std::string usage = "\n"\
    "Usage : ConvertToXdmf \n" \
    " Win32  : ConvertToXdmf -D D:/data/xdmf -F elements.vtu -O D:/data/xdmf/test \n" \
    " Win32  : ConvertToXdmf -D D:/data/VectorField -F damBreakMarin.pvd -O D:/data/xdmf/test \n" \
    " Win32  : ConvertToXdmf -D D:/data/cfxtestdata/Example -F TRS_Default_003.res -O D:/data/xdmf/test \n" \
    " Win32  : ConvertToXdmf -D C:/data/xdmf/test -F TRS_Default_003.xmf -O C:/data/xdmf/scratch \n" \
    " laptop : ConvertToXdmf -D C:\\Data\\cfxtestdata\\Example -F TRS_Default_003.res -O C:\\Data\\cfxtestdata\\xml \n" \
    " Linux  : ConvertToXdmf -r -b -s -D /project/csvis/CFXData/Francis -F 1800_001.res -O /project/csvis/CFXData/XML \n" \
    "  -D /project/csvis/CFXData/Francis \n"        \
    "  -F 1800_001.res \n"                          \
    "  -O output directory \n"                      \
    "  -i timestep - interpolate data using timestep between outputs \n" \
    "  -r rotate vectors for stn frame \n"          \
    "  -f flip 180degrees (csf data) \n"            \
    "  -a append to existing set (otherwise overwrite) \n" \
    "  -s Force Static flag \n" \
    "  -v Subdivide cells (only Tetrahedra + Hexahedra supported currently) \n" \
    "  -Client Send data from DSM Client to DSM Server \n" \
    "  -Server Wait for data from DSM Client \n" \
    "  -DSM Write to self contained DSM \n";
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupDataSetReader(const std::string &filename) {
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupMultiBlockReader(const std::string &filename) {
  vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupCollectionReader(const std::string &filename) {
  vtkSmartPointer<vtkPVDReader> reader = vtkSmartPointer<vtkPVDReader>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupPvtiReader(const std::string &filename) {
  vtkSmartPointer<vtkXMLPImageDataReader> reader = vtkSmartPointer<vtkXMLPImageDataReader>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
/*
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupXdmfReader(const std::string &filename) {
  vtkSmartPointer<vtkXdmfReader> reader = vtkSmartPointer<vtkXdmfReader>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkAlgorithm> SetupXdmfReader3(const std::string &filename) {
  vtkSmartPointer<vtkXdmfReader4> reader = vtkSmartPointer<vtkXdmfReader4>::New();
  reader->SetFileName(filename.c_str());
  return reader;
}
*/
//----------------------------------------------------------------------------
#ifdef USE_PV_ASTRO
vtkSmartPointer<vtkAlgorithm> SetupRamsesReader(const std::string &filename) {
  vtkSmartPointer<vtkRamsesReader> reader = vtkSmartPointer<vtkRamsesReader>::New();
  reader->SetFileName(filename.c_str());
  reader->SetHasParticleData(1);
  return reader;
}
#endif
//----------------------------------------------------------------------------
#ifdef HAS_CFX
vtkSmartPointer<vtkAlgorithm> SetupCFXReader(const std::string &filename) {
  vtkSmartPointer<vtkCFXReader> reader = vtkSmartPointer<vtkCFXReader>::New();
  reader->SetFileName(filename.c_str());
  reader->SetLoadBoundaries(0);
  reader->SetRotateRotatingZoneVectors(1);
  reader->SetLoadVolumes(1);
  return reader;
}
#endif
//----------------------------------------------------------------------------
#ifdef HAS_FLUENT
vtkSmartPointer<vtkAlgorithm> SetupFluentReader(const std::string &filename, bool flip) {
  vtkSmartPointer<vtkFluentReaderCSCS> reader = vtkSmartPointer<vtkFluentReaderCSCS>::New();
  reader->SetDataByteOrderToLittleEndian();
  reader->SetFileName(filename.c_str());
  reader->SetTimeStep(0);
  reader->SetFlip180Degrees(flip);
  return reader;
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
// Call back we can use for future control of processes
void SetStuffRMI(void *localArg, void* vtkNotUsed(remoteArg), 
                    int vtkNotUsed(remoteArgLen), int vtkNotUsed(id))
{ 
  // ParallelRMIArgs_tmp* args = (ParallelRMIArgs_tmp*)localArg;
  // vtkMultiProcessController* contrl = args->Controller;
}
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
    std::cout << "Process number " << myId << " of " << numProcs << std::endl;
  }
  else {
    std::cout << "Process number " << myId << " of " << numProcs << std::endl;
  }
  controller->Barrier();

  //
  // begin scope block so auto-deletion of memory can be checked
  //
  { 
    //
    // Force the creation of our output window object
    //
    vtkSmartPointer<vtkStreamOutputWindow> outwin = vtkSmartPointer<vtkStreamOutputWindow>::New();
    vtkOutputWindow::SetInstance(outwin);
    outwin->SetOutputStream(&std::cout);

    //
    // Use utilities in vtkTesting for getting command line params
    //
    ParallelArgs_tmp* args = reinterpret_cast<ParallelArgs_tmp*>(arg);
    vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
    for (int c=1; c<args->argc; c++ ) {
      test->AddArgument(args->argv[c]);
    }

    //
    // Extract filename from test parameters
    //
    char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
      "-F", args->argc, args->argv, "DUMMY_ENV_VAR", 
      // Default file name
      "unused.abc");
    char* inputname = NULL;
    // Client does not load anything from file
    std::cout << "ConvertToH5Part " << std::endl;
    std::cout << "Got Filename " << filename << std::endl;
    inputname = vtkTestUtilities::ExpandDataFileName(args->argc, args->argv, filename);
    std::string infilename = vtksys::SystemTools::ConvertToOutputPath(inputname);
    std::cout << "Got Fullname " << infilename << std::endl << std::endl;
    delete []inputname;

    if (!vtksys::SystemTools::FileExists(infilename.c_str())) {
      std::cout << "Invalid filename " << infilename << std::endl;
      return;
    }

  #ifndef WIN32
    char *out_dir = vtkTestUtilities::GetArgOrEnvOrDefault(
      "-O", args->argc, args->argv, "DUMMY_ENV_VAR", 
      // Default file name
      "/projects/biddisco/MHD");
  #else
    char *out_dir = vtkTestUtilities::GetArgOrEnvOrDefault(
      "-O", args->argc, args->argv, "DUMMY_ENV_VAR", 
      // Default file name
      "D:\\Temp\\");
  #endif
    std::string file_prefix = vtksys::SystemTools::GetFilenameWithoutExtension(filename);
    std::string hdf_name    = std::string(out_dir) + "/" + file_prefix;
    std::string extension   = vtksys::SystemTools::GetFilenameExtension(filename);
    std::string domainname  = file_prefix;
    hdf_name = vtksys::SystemTools::ConvertToOutputPath(hdf_name.c_str());

    //
    // Is static flag present
    //
    bool Static = false;
    if (test->IsFlagSpecified("-s")) {
      std::cout << "Static (-s flag) set" << std::endl;
      Static = true;
    }

    //
    // Is rotate flag present
    //
    bool rotate = false;
    if (test->IsFlagSpecified("-r")) {
      std::cout << "Rotate (-r flag) set" << std::endl;
      rotate = true;
    }

    //
    // Is Subdivide flag present
    //
    bool subdivide = false;
    if (test->IsFlagSpecified("-v")) {
      std::cout << "Subdivide (-v flag) set" << std::endl;
      subdivide = true;
    }

    //
    // is Flip flag specified
    //
    bool flip = false;
    if (test->IsFlagSpecified("-f")) {
      std::cout << "Flipping 180 degrees " << std::endl;
      flip = true;
    }

    //
    // Is append flag present
    //
    bool append = false;
    if (test->IsFlagSpecified("-a")) {
      std::cout << "Append (-a flag) set" << std::endl;
      rotate = true;
    }

    if (!append && vtksys::SystemTools::FileExists((hdf_name+".h5part").c_str())) {
      vtksys::SystemTools::RemoveFile(std::string(hdf_name+".h5part").c_str());
    }

    //
    // Interpolation requested?
    //
    double interpolation = 0;
    char *number = vtkTestUtilities::GetArgOrEnvOrDefault(
      "-i", args->argc, args->argv, "DUMMY_ENV_VAR", "");
    if (std::string(number)!=std::string("")) {
      std::stringstream temp;
      temp << number;
      temp >> interpolation;
      std::cout << "Interpolating with timestep : " << interpolation << std::endl;
    }
    std::cout << std::endl << std::endl;

    //
    // Start the pipeline with a reader
    //
    vtkSmartPointer<vtkAlgorithm> algorithm, reader;
    vtkSmartPointer<vtkRamsesReader> ramses;
    if (extension==".vtu") {
      reader = SetupDataSetReader(infilename);
    }
    else if (extension==".pvd") {
      reader = SetupCollectionReader(infilename);
    }
    else if (extension==".pvti") {
      reader = SetupPvtiReader(infilename);
    }    
    else if (extension==".xmf" || extension==".xdmf") {
//      reader = SetupXdmfReader(infilename);
    }
    else if (extension==".ramses") {
      reader = SetupRamsesReader(infilename);
      ramses = vtkRamsesReader::SafeDownCast(reader);
    }
#ifdef HAS_CFX
    else if (extension==".res") {
      reader = SetupCFXReader(infilename);
    }
#endif
#ifdef HAS_FLUENT
    else if (extension==".cas" || extension==".cas2") {
      reader = SetupFluentReader(infilename, flip);
    }
#endif
    algorithm = reader;

    vtkSmartPointer<vtkCompositeDataPipeline> defaultexecutive = 
      vtkSmartPointer<vtkCompositeDataPipeline>::New();
    vtkAlgorithm::SetDefaultExecutivePrototype(defaultexecutive);      

    //
    // Interpolate data between time steps?
    //
    vtkSmartPointer<vtkTemporalDataSetCache> cache = vtkSmartPointer<vtkTemporalDataSetCache>::New();
    vtkSmartPointer<vtkTemporalInterpolator> inter = vtkSmartPointer<vtkTemporalInterpolator>::New();
    if (interpolation>0) {
      cache->SetInputConnection(algorithm->GetOutputPort());
      cache->SetCacheSize(2);
      inter->SetInputConnection(cache->GetOutputPort());
      inter->SetDiscreteTimeStepInterval(interpolation);
      //
      algorithm = inter;
    }

    //
    // Set Update piece information so pipeline can split correctly
    // actually, it's invalid to do it here, but we're testing
    //
    vtkSmartPointer<vtkInformation> execInfo;
    execInfo = algorithm->GetExecutive()->GetOutputInformation(0);
    execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numProcs);
    execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), myId);

    //
    // update information to set/get TIME_XXX keys and Extents
    //
    algorithm->UpdateInformation();

    //
    // Split into pieces for parallel loading
    //
    // int MaxPieces = execInfo->Get(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES());
    int WholeExtent[6], UpdateExtent[6];
    execInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), WholeExtent);
    vtkSmartPointer<vtkExtentTranslator> extTran = vtkSmartPointer<vtkExtentTranslator>::New();
      extTran->SetSplitModeToBlock();
      extTran->SetNumberOfPieces(numProcs);
      extTran->SetPiece(myId);
      extTran->SetWholeExtent(WholeExtent);
      extTran->PieceToExtent();
      extTran->GetExtent(UpdateExtent);
    execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), UpdateExtent, 6);

    //
    // Get time values
    //    
    std::vector<double> TimeSteps;
    if (execInfo && execInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) 
    {
      int NumberOfTimeSteps = execInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      TimeSteps.assign(NumberOfTimeSteps, 0.0);
      execInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &TimeSteps[0]);
    }
    else {
      TimeSteps.assign(1, 0.0);
    }


    //
    // if we are running a parallel job, sync here
    //
    controller->Barrier();

    //-------------------------------------------------------------
    // Limit time steps for testing to 3 so we don't wait all day
    //-------------------------------------------------------------
    itype t0=0, t1=TimeSteps.size();
//    if (t1<2) {
//      t1=10;
//      TimeSteps.assign(t1, 1);
//    t1 = TimeSteps.size()>3 ? 3 : TimeSteps.size();
//    }

    //
    // Create the Xdmf Writer
    //
    vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();


    //
    // Write out time steps to File or DSM.
    // We loop twice for ramses, first time, using ReadHeaderOnly flag to get the time
    // from each dataset. Then the second time we actually do the read.
    //
    for (int loop=0; loop<(ramses?2:1); ++loop) {
            
      for (itype t=t0; t<t1; ++t) {
        //
        // Update the time we want
        //
        itype TimeStep = t; // was missingsteps[t]
        double current_time = TimeSteps[TimeStep];
        //
        std::cout << "Beginning TimeStep " << setw(5) << setfill('0') << t 
          << " for time " << setw(8) << setfill('0') << current_time << std::endl;
        //
        if (ramses && loop==1) {
          double simplet = t;
          execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS(), &simplet, 1);
        }
        else {
          execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS(), &current_time, 1);
        }
  #ifdef PARALLEL_PIECES
        execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numProcs);
        execInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), myId);
  #endif
        if (ramses) {
          ramses->SetReadHeaderOnly(loop==0);
        }
        algorithm->Update();
        if (ramses && loop==0) {
          vtkDataSet *out = ramses->GetOutput();
          vtkFieldData *fd = out->GetFieldData();
          TimeSteps[TimeStep] = fd->GetArray("time")->GetTuple1(0);
          ramses->SetReadHeaderOnly(loop==0);
        }
        else {
          //
          // get the output temporal/multiblock - depends if we are interpolating or not
          //
          vtkSmartPointer<vtkTemporalDataSet>     temporal = vtkTemporalDataSet::SafeDownCast(algorithm->GetOutputDataObject(0));
          vtkSmartPointer<vtkMultiBlockDataSet> multiblock = vtkMultiBlockDataSet::SafeDownCast(algorithm->GetOutputDataObject(0));
          vtkSmartPointer<vtkDataSet>              dataset = vtkDataSet::SafeDownCast(algorithm->GetOutputDataObject(0));
          if (temporal) {
            multiblock = vtkMultiBlockDataSet::SafeDownCast(temporal->GetTimeStep(0));
          }

          //
          // int d=0;
          if (multiblock || dataset) {

            if (multiblock) writer->SetInput(multiblock);
            if (dataset)    writer->SetInput(dataset);
            //
#ifdef VTK_USE_MPI
            writer->SetController(controller);
#endif
            writer->SetFileName(std::string(hdf_name+".h5part").c_str());
            writer->SetTimeValue(current_time);
            writer->SetTimeStep((int)TimeStep);
            writer->Modified();
            writer->Update();
            writer->CloseFile();
          }
        }
      }
    }

    //
    // After all time steps have been written.
    //
    delete []filename;
    //
    controller->Barrier();
   } // end scope block
}

//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  // This is here to avoid false leak messages from vtkDebugLeaks when
  // using mpich. It appears that the root process which spawns all the
  // main processes waits in MPI_Init() and calls exit() when
  // the others are done, causing apparent memory leaks for any objects
  // created before MPI_Init().
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //
  if (rank == 0) {
  }

  // Note that this will create a vtkMPIController if MPI
  // is configured, vtkThreadedController otherwise.
#ifdef VTK_USE_MPI
  vtkMPIController* controller = vtkMPIController::New();
#else
  vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
#endif
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

  return !retVal;
}
