//
#include <iostream>
#include <sstream>
//
// For VTK_USE_MPI 
#include "vtkToolkits.h"     
#ifdef VTK_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
// Otherwise
#include "vtkDummyController.h"
//
#include "Testing/Cxx/vtkTestUtilities.h"
#include "Testing/Cxx/vtkRegressionTestImage.h"
//
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkStreamingDemandDrivenPipeline.h"

// CSCS
#include "vtkStreamOutputWindow.h"
#include "vtkH5PartReader.h"
#include "vtkSPHManager.h"
#include "vtkSPHProbeFilter.h"
#include "vtkSPHImageResampler.h"
// only defined if trilinos used
class vtkParticlePartitionFilter;
//----------------------------------------------------------------------------
#if 0
  #define OUTPUTTEXT(a) std::cout << (a);
  #define testDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << a << endl; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }
#else
  #define testDebugMacro(a) 
#endif
//----------------------------------------------------------------------------
#define DEBUG_WAIT \
  if (test.myRank==0) { \
    char ch;       \
    std::cout << "Attach debugger" << std::endl; \
    std::cin >> ch; \
  }
//----------------------------------------------------------------------------
class TestStruct {
 public:
  //
  vtkSmartPointer<vtkMultiProcessController>  controller;
  vtkSmartPointer<vtkH5PartReader>            reader;
  vtkSmartPointer<vtkParticlePartitionFilter> partitioner;
  vtkSmartPointer<vtkSPHManager>              sphManager;
  vtkSmartPointer<vtkAlgorithm>               sphResampler;
  //
  vtkTypeInt64 myRank;
  vtkTypeInt64 numProcs;
  //--------------------------------------------------------------
  // Testing params
  //--------------------------------------------------------------
  bool   unused, fixNeighbours, fixRadius, cameraSet;
  double gridSpacing[3];
  int    gridResolution[3];
  double vminmax[2];
  double vpos[3];
  double cameraPosition[3];
  double cameraFocus[3];
  double cameraViewUp[3];
  int    windowSize[2];
  //
  std::string testName;

  bool        doRender;
  bool        keepTempFiles;
  //
  // (Random) Particle Generation
  //
  vtkIdType   generateN;

  //
  // H5Part Reader 
  //
  bool        ReadData;
  std::string fullName;
  std::string Xarray;
  std::string Yarray;
  std::string Zarray;
  bool        ignorePartitions;
  bool        randomizeExtents;

  //
  // SPH kernel or neighbour info
  //
  double      particleSize;
  double      ghostOverlap;
  int         maxN;
  std::string massScalars;
  std::string densityScalars;
  vtkIdType   expectedN;

  //
  // Test/Display of results
  //
  std::string scalarname;
  double      contourVal;
  bool        imageResample;
  bool        skipImageTest;
  std::string imageScalars;
  int         imageThreshold;
  bool        benchmarkPartition;
  //
  void    CreateReader();
  double  UpdateReader();
  void    CreatePartitioner();
  double  UpdatePartitioner();
  void    DeletePartitioner();
  void    CreateSPHManager();
  void    CreateSPHResampler(vtkAlgorithm *input);
  double  UpdateSPHResampler();
  //
};
//----------------------------------------------------------------------------
int initTest(int argc, char* argv[], TestStruct &test);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
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
    double nextNumber() {
      __seed = (__seed*9301+49297) % 233280;
      return __seed / 233280.0;
    }
    int nextNumberInt() {
      __seed = (__seed*9301+49297) % 233280;
      return __seed;
    }
};
//----------------------------------------------------------------------------
unsigned long int random_seed();
void known_seed();
void SpherePoints(int n, float radius, float X[]);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
template <typename T>
void DisplayParameter(const char *prefix1, const char *prefix2, T *data, int components, int rank)
{
  if (rank==-1) {
    return;
  }
  vtkstd::stringstream temp;
  temp << prefix1 << prefix2 << std::ends;
  std::cout.width(30);
  std::cout << temp.str().c_str() << " : {";
  std::cout.width(0);
  for (int i=0; i<components; i++) {
    std::cout << data[i];
    (i==(components-1)) ? std::cout << "}" : std::cout << ",";
  }
  std::cout << std::endl;
}
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
      DisplayParameter<T>(message, "", &newValue, 1, rank);
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
      std::cout << message << " : {";
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
