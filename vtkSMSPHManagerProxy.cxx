#include "vtkSMSPHManagerProxy.h"
#include "vtkObjectFactory.h"
#include "vtkSMProxyManager.h"
#include "vtkProcessModuleConnectionManager.h"
#include "vtkProcessModule.h"
#include "vtkSPHManager.h"
#include "vtkCallbackCommand.h"
//-----------------------------------------------------------------------------
vtkSMSPHManagerProxy *vtkSMSPHManagerProxy::ReferenceProxy = NULL;
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSMSPHManagerProxy);
//----------------------------------------------------------------------------
vtkSMSPHManagerProxy::vtkSMSPHManagerProxy() : vtkSMProxy()
{
/*
  this->SPHSingleton = vtkSPHManager::New();
  this->SPHObserver = vtkCallbackCommand::New();
  this->SPHObserver->SetCallback(
    &vtkSMSPHManagerProxy::SPHModifiedCallback);
  this->SPHObserver->SetClientData(this);
//  this->SPHSingleton->AddObserver(vtkCommand::ModifiedEvent, this->SPHObserver);
  //
  //
  //
*/
}
//----------------------------------------------------------------------------
vtkSMSPHManagerProxy::~vtkSMSPHManagerProxy()
{
  vtkSMSPHManagerProxy::ReferenceProxy = NULL;
//  this->SPHSingleton->RemoveObserver(this->SPHObserver);
//  this->SPHSingleton->Delete();
}
//----------------------------------------------------------------------------
void vtkSMSPHManagerProxy::CreateVTKObjects()
{
  if (!vtkSMSPHManagerProxy::ReferenceProxy) {
    vtkSMSPHManagerProxy::ReferenceProxy = this;
    vtkSMProxy::CreateVTKObjects();
  }
  else {
    this->InitializeAndCopyFromProxy(vtkSMSPHManagerProxy::ReferenceProxy);
  }
}
//----------------------------------------------------------------------------
void vtkSMSPHManagerProxy::SPHModifiedCallback(
  vtkObject* caller, unsigned long eid, void* clientdata, void* calldata)
{
  static bool incallback = true;
  if (!incallback) {
    static_cast<vtkSMSPHManagerProxy*>(clientdata)->Modified();
    static_cast<vtkSMSPHManagerProxy*>(clientdata)->UpdatePropertyInformation();
  }
  incallback = false;
}
//----------------------------------------------------------------------------
