#include "TestUtils.h"
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
#include "vtkH5PartReader.h"
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
#include "vtkParticleIdFilter.h"
//
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
  #include "vtkParticlePartitionFilter.h"
#endif
//

// #define DEBUG_WAIT 1
//----------------------------------------------------------------------------
static const int CONTOURDATA_TAG=301;
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  int retVal = 1;
  const char *empty = "";

  //--------------------------------------------------------------
  // Setup Test Params
  //--------------------------------------------------------------
  TestStruct test;
  initTest(argc, argv, test);

  //--------------------------------------------------------------
  // Create H5Part Reader
  //--------------------------------------------------------------
  test.CreateReader();
  double read_elapsed = test.UpdateReader();
  //
  vtkIdType totalParticles = 0;
  vtkIdType localParticles = test.reader->GetOutput()->GetNumberOfPoints();
  test.controller->AllReduce(&localParticles, &totalParticles, 1, vtkCommunicator::SUM_OP);
  DisplayParameter<vtkIdType>("Particle Count", "", &totalParticles, 1, (test.myRank==0)?0:-1);

  vtkSmartPointer<vtkAlgorithm> data_algorithm = test.reader; 
  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  test.CreatePartitioner();
  double partition_elapsed = test.UpdatePartitioner();
  //
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
  data_algorithm = test.partitioner;
#endif

  //--------------------------------------------------------------
  // SPH Manager
  //--------------------------------------------------------------
  test.CreateSPHManager();

  //--------------------------------------------------------------
  // SPH Resample
  //--------------------------------------------------------------
  test.CreateSPHResampler(data_algorithm);
  double sph_elapsed = test.UpdateSPHResampler();

  //--------------------------------------------------------------
  // Fetch smoothed output for testing results
  //--------------------------------------------------------------
  vtkDataSet *sph_results = vtkDataSet::SafeDownCast(
    test.sphResampler->GetOutputDataObject(0));

  int wholeExtent[6] = {0,-1,0,-1,0,-1};
  if (test.imageResample) {
    vtkImageData::SafeDownCast(sph_results)->GetInformation()->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),wholeExtent);
  }

  vtkSmartPointer<vtkTimerLog> viztimer = vtkSmartPointer<vtkTimerLog>::New();
  viztimer->StartTimer();

  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------
  bool ok = true;
  std::cout.precision(15);
  
  //
  // Testing use scalar values extracted from filters
  //
  if (!test.imageResample) {
    vtkPointData *sph_pd = sph_results->GetPointData();
    vtkDataArray *scalar_array = sph_pd->GetArray(test.scalarname.c_str());
    //--------------------------------------------------------------
    // Collect results from parallel filters
    //--------------------------------------------------------------
    double scalar_range_local[2], scalar_range_global[2];
    double scalar_pos[3];
    //
    scalar_array->GetRange(scalar_range_local);
    test.controller->AllReduce(&scalar_range_local[0], &scalar_range_global[0], 1, vtkCommunicator::MIN_OP);
    test.controller->AllReduce(&scalar_range_local[1], &scalar_range_global[1], 1, vtkCommunicator::MAX_OP);
    DisplayParameter<double>(test.scalarname.c_str(), " Found Min/Max", scalar_range_global, 2, (test.myRank==0)?0:-1);
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
      vtkIdType ids[2] = {test.myRank, index};
      DisplayParameter<vtkIdType>(test.scalarname.c_str(), " Peak Proc,Index", ids, 2, (test.myRank==0)?0:-1);
      sph_results->GetPoint(index,scalar_pos);
      DisplayParameter<double>(test.scalarname.c_str(), " Peak Position", scalar_pos, 3, (test.myRank==0)?0:-1);
      //
      double tol_min = std::fabs(test.vminmax[0]/1000.0);
      double tol_max = std::fabs(test.vminmax[1]/1000.0);
      if (std::fabs(test.vminmax[0]-scalar_range_global[0])>tol_min || std::fabs(test.vminmax[1]-scalar_range_global[1])>tol_max) {
        ok = false;
        DisplayParameter<const char *>("++++++++++++++++++++", "", &empty, 1, (test.myRank==0)?0:-1);
        std::cout << "min/max check failed " << std::endl;
        std::cout << "expected {" << test.vminmax[0] << ',' << test.vminmax[1] << "}" << std::endl;
        std::cout << "got {" << scalar_range_global[0] << ',' << scalar_range_global[1] << "}" << std::endl;
        std::cout << "err {" << std::abs(test.vminmax[0]-scalar_range_global[0]) << ',' << std::abs(test.vminmax[1]-scalar_range_global[1]) << "}" << std::endl;
      }
      if (std::fabs(test.vpos[0]-scalar_pos[0])>1E-5 ||
          std::fabs(test.vpos[1]-scalar_pos[1])>1E-5 ||
          std::fabs(test.vpos[2]-scalar_pos[2])>1E-5)
      {
        ok = false;
        DisplayParameter<const char *>("++++++++++++++++++++", "", &empty, 1, (test.myRank==0)?0:-1);
        std::cout << "position check failed " << std::endl;
        std::cout << "expected {" << test.vpos[0] << "," << test.vpos[1] << ',' << test.vpos[2] << "}" << std::endl;
        std::cout << "got {" << scalar_pos[0] << "," << scalar_pos[1] << ',' << scalar_pos[2] << "}" << std::endl;
        std::cout << "err {" << std::abs(scalar_pos[0]-test.vpos[0]) << ',' << std::abs(scalar_pos[1]-test.vpos[1]) << ',' << std::abs(scalar_pos[2]-test.vpos[2]) << "}" << std::endl;
      }
    }
  }
  //
  // Testing using the standard VTK Image Regression
  //
  else if (!test.benchmarkPartition) {
    //
    vtkSmartPointer<vtkContourFilter> iso = vtkSmartPointer<vtkContourFilter>::New();
    iso->SetInputConnection(test.sphResampler->GetOutputPort());
    iso->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,test.scalarname.c_str());
    iso->SetValue(0, test.contourVal);
    iso->ComputeScalarsOff();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();
    //
    //
    vtkSmartPointer<vtkAlgorithm> vis_algorithm; 
    vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
    processId->SetInputConnection(iso->GetOutputPort());
    processId->SetController(test.controller);
    //
    vis_algorithm = processId;
    //
    vtkStreamingDemandDrivenPipeline *vis_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(vis_algorithm->GetExecutive());
    // no piece info set yet, assumes info is not piece dependent
    testDebugMacro( "Setting viz piece information " << test.myRank << " of " << test.numProcs );
    vis_sddp->SetUpdateExtent(0, test.myRank, test.numProcs, 0);
//    vis_sddp->SetUpdateResolution(0, 1.0);
    testDebugMacro( "Viz UpdateInformation " );
    vis_sddp->UpdateInformation();
    testDebugMacro( "Viz Update " );
    vis_sddp->SetUpdateExtent(0, test.myRank, test.numProcs, 0);
    vis_sddp->Update();
    testDebugMacro( "Viz Done " );
    //
    vtkSmartPointer<vtkPolyData> contourData = vtkSmartPointer<vtkPolyData>::New();
    contourData->ShallowCopy(vis_sddp->GetOutputData(0));

    //
    // Release all the probing pipeline to free memory so that process zero
    // is less likely to crash when collecting results.
    //
    test.DeletePartitioner();

    test.sphResampler->SetInputConnection(NULL);
    iso->SetInputConnection(NULL);
    processId->SetInputConnection(NULL);
    test.reader       = NULL;
    test.sphManager   = NULL;
    test.sphResampler = NULL;
    iso               = NULL;
    processId         = NULL;
    test.controller->Barrier();

    //
    // Rank >0 Send contour pieces to rank 0
    //
    if (test.myRank>0) {
      const char *array_name = "ProcessId";
      vtkDataArray *da = contourData->GetPointData()->GetArray(array_name);
      if (da) {      
        double *rg = da->GetRange();
        testDebugMacro( "Range (Sending) " << rg[0] << " " << rg[1] );
      }
      test.controller->Send(contourData, 0, CONTOURDATA_TAG);
    }
    //
    // Rank 0 collect all contour pieces from parallel processes
    //
    else if (test.myRank==0) {
      //
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      vtkSmartPointer<vtkAppendPolyData>       append = vtkSmartPointer<vtkAppendPolyData>::New();
      //
      for (int i=0; i<test.numProcs; i++) {
        vtkSmartPointer<vtkPolyData> pd;
        if (i==0) {
          pd = contourData;
        }
        else {
          pd = vtkSmartPointer<vtkPolyData>::New();
          test.controller->Receive(pd, i, CONTOURDATA_TAG);
        }
        append->AddInputData(pd);
        const char *array_name = "ProcessId";
        vtkDataArray *da = pd->GetPointData()->GetArray(array_name);
        if (da) {
          double *rg = da->GetRange();
          testDebugMacro( "Range (Receive) " << rg[0] << " " << rg[1] );
        }
        //
        // Display boxes for each partition
        //
        vtkSmartPointer<vtkOutlineFilter> boxsource = vtkSmartPointer<vtkOutlineFilter>::New();
        boxsource->SetInputData(pd);
        vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>          bactor = vtkSmartPointer<vtkActor>::New();
        bmapper->SetInputConnection(boxsource->GetOutputPort());
        bmapper->SetImmediateModeRendering(1);
        bactor->SetMapper(bmapper);
        ren->AddActor(bactor);
      }
      vtkSmartPointer<vtkPolyDataNormals> norm = vtkSmartPointer<vtkPolyDataNormals>::New();
      norm->SetInputConnection(append->GetOutputPort());
      norm->Update();
      //
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize(test.windowSize);
      mapper->SetInputConnection(norm->GetOutputPort());
      mapper->SetImmediateModeRendering(1);
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointFieldData();
      mapper->SetInterpolateScalarsBeforeMapping(0);
      mapper->SetUseLookupTableScalarRange(0);
      vtkDataArray *da = append->GetOutput()->GetPointData()->GetArray(test.imageScalars.c_str());
      mapper->SetScalarRange(da->GetRange());
      mapper->SelectColorArray(test.imageScalars.c_str());
      actor->SetMapper(mapper);
      actor->GetProperty()->SetPointSize(2);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);
      //
      if (test.cameraSet) {
        ren->GetActiveCamera()->SetPosition(test.cameraPosition);
        ren->GetActiveCamera()->SetFocalPoint(test.cameraFocus);
        ren->GetActiveCamera()->SetViewUp(test.cameraViewUp);
        ren->ResetCameraClippingRange();
      }
      else {
        ren->ResetCamera();
      }
      //
      testDebugMacro( "Process Id : " << test.myRank << " About to Render" );
      renWindow->Render();
      //
      if (test.skipImageTest) {
        retVal=vtkRegressionTester::PASSED;
        ok = true;
      }
      else {
        retVal = vtkRegressionTester::Test(argc, argv, renWindow, test.imageThreshold);
        if (retVal == vtkRegressionTester::DO_INTERACTOR) {
          iren->Start();
        }
        ok = (retVal==vtkRegressionTester::PASSED);
      }
      testDebugMacro( "Process Id : " << test.myRank << " Rendered " << (ok?"Pass":"Fail"));
    }
  }

  test.controller->Barrier();
  viztimer->StopTimer();
  double viz_elapsed = viztimer->GetElapsedTime();

  if (ok && test.myRank==0) {
    DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, test.myRank);
    DisplayParameter<double>("Read Time", "", &read_elapsed, 1, test.myRank);
#ifdef PV_MESHLESS_ZOLTAN_SUPPORT
//    DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, test.myRank);
#endif
    DisplayParameter<double>("SPH Probe Time", "", &sph_elapsed, 1, test.myRank);
    DisplayParameter<double>("Visualization/Check Time", "", &viz_elapsed, 1, test.myRank);
    if (test.imageResample) {
      DisplayParameter<int>("Image Whole Extent", "", wholeExtent, 6, test.myRank);
      vtkIdType voxels = (1+wholeExtent[1]-wholeExtent[0])*(1+wholeExtent[3]-wholeExtent[2])*(1+wholeExtent[5]-wholeExtent[4]);
      DisplayParameter<vtkIdType>("Voxels", "", &voxels, 1, test.myRank);
    }
    DisplayParameter<const char *>("====================", "", &empty, 1, test.myRank);
  }
  
  test.controller->Finalize();
  finalizeTest(test);

  retVal = ok;

  return !retVal;
}
//----------------------------------------------------------------------------

