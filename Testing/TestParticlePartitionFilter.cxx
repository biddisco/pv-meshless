// mpishim : C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\Remote Debugger\x64\mpishim100.exe
// mpiexec : C:\Program Files\MPICH2\bin\mpiexec.exe
// mpiargs : -localonly -n 2 -env PATH C:\cmakebuild\pv-meshless\bin\debug;c:\bin

#ifdef _WIN32
  #include <windows.h>
#else
  #include <sys/time.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
//
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkMultiProcessController.h"
//
#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkPointSource.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkDebugLeaks.h"
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
#include "TestUtils.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#define DATA_SEND_TAG 301
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  int retVal = 1;
  char *empty = "";
  bool ok = true;

  //--------------------------------------------------------------
  // Setup Test Params
  //--------------------------------------------------------------
  TestStruct test;
  initTest(argc, argv, test);

  // if testing partition from file
  double read_elapsed = 0.0;
  double partition_elapsed = 0.0;
  vtkSmartPointer<vtkAlgorithm> data_algorithm;
  vtkIdType totalParticles = 0;
  if (vtksys::SystemTools::FileExists(test.fullName.c_str())) {
    //--------------------------------------------------------------
    // Create H5Part Reader
    //--------------------------------------------------------------
    test.CreateReader();
    read_elapsed = test.UpdateReader();
    data_algorithm = test.reader;
    //
    vtkIdType localParticles = test.reader->GetOutput()->GetNumberOfPoints();
    test.controller->AllReduce(&localParticles, &totalParticles, 1, vtkCommunicator::SUM_OP);
    DisplayParameter<vtkIdType>("Particle Count", "", &totalParticles, 1, (test.myRank==0)?0:-1);

    //--------------------------------------------------------------
    // Parallel partition
    //--------------------------------------------------------------
    test.CreatePartitioner();
    test.partitioner->SetGhostCellOverlap(test.ghostOverlap);
    partition_elapsed = test.UpdatePartitioner();
  }
  else {
    //--------------------------------------------------------------
    // allocate points + arrays
    //--------------------------------------------------------------
    vtkSmartPointer<vtkPolyData>  Sprites = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>   verts = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIdTypeArray>   Ids = vtkSmartPointer<vtkIdTypeArray>::New();
    vtkSmartPointer<vtkIntArray>    Ranks = vtkSmartPointer<vtkIntArray>::New();
    //
    points->SetNumberOfPoints(test.generateN);
    //
    verts->Allocate(test.generateN,test.generateN);
    Sprites->SetPoints(points);
    Sprites->SetVerts(verts);
    //
    Ids->SetNumberOfTuples(test.generateN);
    Ids->SetNumberOfComponents(1);
    Ids->SetName("PointIds");
    Sprites->GetPointData()->AddArray(Ids);
    //
    Ranks->SetNumberOfTuples(test.generateN);
    Ranks->SetNumberOfComponents(1);
    Ranks->SetName("Rank");
    Sprites->GetPointData()->AddArray(Ranks);
    //
    //--------------------------------------------------------------
    // Create default scalar arrays
    //--------------------------------------------------------------
    double radius  = 500.0;
    const double a = 0.9;
    test.ghostOverlap = radius*0.1; // ghost_region

    known_seed();
    SpherePoints(test.generateN, radius*(1.5+test.myRank)/(test.numProcs+0.5), vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0));
    for (vtkIdType Id=0; Id<test.generateN; Id++) {
      Ids->SetTuple1(Id, Id + test.myRank*test.generateN);
      Ranks->SetTuple1(Id, test.myRank);
      verts->InsertNextCell(1,&Id);
    }
/*
    // Randomly give some processes zero points to improve test coverage
    random_seed();
    if (test.numProcs>1 && rand()%2==1) {
      test.generateN = 0;
      Sprites = vtkSmartPointer<vtkPolyData>::New();
    }
*/
    //--------------------------------------------------------------
    // Parallel partition
    //--------------------------------------------------------------
    test.CreatePartitioner();
    test.partitioner->SetInputData(Sprites);
    test.partitioner->SetIdChannelArray("PointIds");
    test.partitioner->SetGhostCellOverlap(test.ghostOverlap);
    partition_elapsed = test.UpdatePartitioner();

    //--------------------------------------------------------------
    // Add process Id's
    //--------------------------------------------------------------
    vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
    processId->SetInputConnection(test.partitioner->GetOutputPort());
    processId->SetController(test.controller);
    //
    test.controller->Barrier();
    if (test.myRank==0) {
      testDebugMacro( "Process Id : " << test.myRank << " Generated N Points : " << test.generateN );
    }

    //--------------------------------------------------------------
    // Update in parallel
    //
    // To get parallel operation correct, we need to make sure that piece
    // information is passed upstream. first update information,
    // then set piece update extent,
    //--------------------------------------------------------------
    testDebugMacro( "Setting piece information " << test.myRank << " of " << test.numProcs );
    vtkStreamingDemandDrivenPipeline *sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(processId->GetExecutive());
    // no piece info set yet, assumes info is not piece dependent
    sddp->UpdateInformation();
    // now set piece info and update
    sddp->SetUpdateExtent(0, test.myRank, test.numProcs, 0);
    sddp->Update();

    if (test.doRender) {
      //
      // Send all the data to process zero for display
      //
      vtkPolyData *OutputData = vtkPolyData::SafeDownCast(sddp->GetOutputData(0)->NewInstance());
      OutputData->ShallowCopy(sddp->GetOutputData(0));
      if (test.myRank>0) {
        test.controller->Send(OutputData, 0, DATA_SEND_TAG);
      }
      //
      // Rank 0 collect all data pieces from parallel processes
      //
      else if (test.myRank==0) {
        //
        vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(renWindow);
        ren->SetBackground(0.1, 0.1, 0.1);
        renWindow->SetSize( 400+8, 400+8);
        renWindow->AddRenderer(ren);
        //
        for (int i=0; i<test.numProcs; i++) {
          vtkSmartPointer<vtkPolyData> pd;
          if (i==0) {
            pd = OutputData;
          }
          else {
            pd = vtkSmartPointer<vtkPolyData>::New();
            test.controller->Receive(pd, i, DATA_SEND_TAG);
          }
          vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
          vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
          mapper->SetInputData(pd);
          mapper->SetColorModeToMapScalars();
          mapper->SetScalarModeToUsePointFieldData();
          mapper->SetUseLookupTableScalarRange(0);
          mapper->SetScalarRange(0,test.numProcs-1);
          mapper->SetInterpolateScalarsBeforeMapping(0);
          mapper->SelectColorArray("ProcessId");
          actor->SetMapper(mapper);
          actor->GetProperty()->SetPointSize(2);
          ren->AddActor(actor);
          //
          vtkSmartPointer<vtkPolyDataMapper>       mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
          vtkSmartPointer<vtkActor>                 actor2 = vtkSmartPointer<vtkActor>::New();
          mapper2->SetInputData(pd);
          mapper2->SetColorModeToMapScalars();
          mapper2->SetScalarModeToUsePointFieldData();
          mapper2->SetUseLookupTableScalarRange(0);
          mapper2->SetScalarRange(0,1);
          mapper2->SetInterpolateScalarsBeforeMapping(0);
          mapper2->SelectColorArray("vtkGhostType");
          actor2->SetMapper(mapper2);
          actor2->GetProperty()->SetPointSize(2);
          actor2->SetPosition(2.0*radius, 0.0, 0.0);
          ren->AddActor(actor2);
        }
        //
        // Display boxes for each partition
        //
        for (int i=0; i<test.numProcs; i++) {
          vtkBoundingBox *box = test.partitioner->GetPartitionBoundingBox(i);
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
        testDebugMacro( "Process Id : " << test.myRank << " About to Render" );
        renWindow->Render();

        retVal = vtkRegressionTester::Test(argc, argv, renWindow, 10);
        if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
          iren->Start();
        }
        ok = (retVal==vtkRegressionTester::PASSED);
        testDebugMacro( "Process Id : " << test.myRank << " Rendered " << (ok?"Pass":"Fail"));
      }
    }
  }

  if (ok && test.myRank==0) {
    DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, test.myRank);
    DisplayParameter<double>("Read Time", "", &read_elapsed, 1, test.myRank);
    DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, test.myRank);
    DisplayParameter<char *>("====================", "", &empty, 1, test.myRank);
  }

  test.controller->Finalize();

  return !retVal;
}
//----------------------------------------------------------------------------
