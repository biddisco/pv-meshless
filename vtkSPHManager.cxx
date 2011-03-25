/*=========================================================================

  Project                 : vtkCSCSMeshless
  Module                  : vtkSPHManager.cxx

  Authors:
     John Biddiscombe     Jerome Soumagne
     biddisco@cscs.ch     soumagne@cscs.ch

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkSPHProbeFilter2 - SPH Convenience Manager Class
// .SECTION Description
#include "vtkSPHManager.h"
//
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
//
#include "KernelGaussian.h"
#include "KernelWendland.h"
#include "KernelQuadratic.h"
#include "KernelSpline3rdOrder.h"
#include "KernelSpline5thOrder.h"
//
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
extern vtkObject* vtkInstantiatorvtkSPHManagerNew(); 
vtkSPHManager *vtkSPHManager::SPHManagerSingleton = NULL;
//----------------------------------------------------------------------------
#undef ErrorMacro
#define ErrorMacro(x)                                           \
   {                                                            \
   if (vtkObject::GetGlobalWarningDisplay())                    \
     {                                                          \
     vtkOStreamWrapper::EndlType endl;                          \
     vtkOStreamWrapper::UseEndl(endl);                          \
     vtkOStrStreamWrapper vtkmsg;                               \
     vtkmsg << "ERROR: In " __FILE__ ", line " << __LINE__      \
            << "\n" x << "\n\n";                                \
     vtkOutputWindowDisplayErrorText(vtkmsg.str());             \
     vtkmsg.rdbuf()->freeze(0); vtkObject::BreakOnError();      \
     }                                                          \
   }
//----------------------------------------------------------------------------
vtkSPHManager *vtkSPHManager::New() 
{
  if (vtkSPHManager::SPHManagerSingleton) {
    // increase the reference count
    vtkSPHManager::SPHManagerSingleton->Register(NULL);
    return vtkSPHManager::SPHManagerSingleton;
  }
  //
  vtkSPHManager *temp = new vtkSPHManager();
  if (temp!=vtkSPHManager::SPHManagerSingleton) {
    ErrorMacro(<<"A serious (probably thread related error) has occured in vtkSPHManager singleton creation")
  }
  return vtkSPHManager::SPHManagerSingleton;
}
//----------------------------------------------------------------------------
vtkObject *vtkInstantiatorvtkSPHManagerNew()
{
  return vtkSPHManager::New();
}
//----------------------------------------------------------------------------
// Singleton creation, not really thread-safe, but unlikely to ever be tested
vtkSPHManager::vtkSPHManager()
{
  if (vtkSPHManager::SPHManagerSingleton==NULL) {
    vtkSPHManager::SPHManagerSingleton = this;
  }

  // switch shepard/sph
  this->InterpolationMethod = POINT_INTERPOLATION_KERNEL;

  // Shepard Mode
  this->LimitSearchByNeighbourCount = 1;
  this->MaximumNeighbours           = 128;
  this->MaximumRadius               = 0.1;

  // SPH Mode
  this->KernelType          = SPH_KERNEL_SPLINE_3RD;
  this->KernelDimension     = 3;
  this->HScalarsRegex       = NULL;
  this->VolumeScalarsRegex  = NULL;
  this->MassScalarsRegex    = NULL;
  this->DensityScalarsRegex = NULL;

  // SPH Variables
  this->DefaultParticleSideLength = 0.18333;
  this->DefaultDensity            = 1000.0;
  this->HCoefficient              = 1.5; 
}
//----------------------------------------------------------------------------
vtkSPHManager::~vtkSPHManager()
{
  this->SetHScalarsRegex(NULL);
  this->SetVolumeScalarsRegex(NULL);
  this->SetMassScalarsRegex(NULL);
  this->SetDensityScalarsRegex(NULL);
  vtkSPHManager::SPHManagerSingleton=NULL;
}
//----------------------------------------------------------------------------
