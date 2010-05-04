#include "pqSPHManagerPanel.h"

// Qt includes
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QVariant>
#include <QLabel>
#include <QComboBox>
#include <QTableWidget>
#include <QMessageBox>
#include <QProgressDialog>
#include <QTimer>
#include <QInputDialog>
#include <QFileDialog>
#include <QUrl>
#include <QDesktopServices>
#include <QThread>

// VTK includes

// ParaView Server Manager includes
#include "vtkSMInputProperty.h"
#include "vtkSMProxyManager.h"
#include "vtkSMSourceProxy.h"
#include "vtkSMStringVectorProperty.h"
#include "vtkSMIntVectorProperty.h"
#include "vtkSMArraySelectionDomain.h"
#include "vtkSMProxyProperty.h"
#include "vtkSMViewProxy.h"
#include "vtkSMRepresentationProxy.h"
#include "vtkProcessModule.h"
#include "vtkProcessModuleConnectionManager.h"

// ParaView includes
#include "pqActiveServer.h"
#include "pqApplicationCore.h"
#include "pqSettings.h"
#include "pqOutputPort.h"
#include "pqPipelineSource.h"
#include "pqPropertyLinks.h"
#include "pqProxy.h"
#include "pqServer.h"
#include "pqServerManagerSelectionModel.h"
#include "pqServerManagerModelItem.h"
#include "pqServerManagerModel.h"
#include "pqSMAdaptor.h"
#include "pqTreeWidgetCheckHelper.h"
#include "pqTreeWidgetItemObject.h"
#include "pqTreeWidget.h"
#include "pqTreeWidgetItem.h"
#include "pqView.h"
#include "pqRenderView.h"
#include "pqActiveView.h"
#include "pqDataRepresentation.h"
#include "pqActiveObjects.h"
#include "pqDisplayPolicy.h"
#include "pqAnimationScene.h"

//
#include "ui_pqSPHManagerPanel.h"

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
class pqSPHManagerPanel::pqUI : public QObject, public Ui::SPHManagerPanel 
{
public:
  pqUI(pqSPHManagerPanel* p) : QObject(p)
  {
    this->Links = new pqPropertyLinks;
    this->SPHInitialized   = 0;
    this->ActiveSourcePort = 0;
    this->ActiveServer     = 0;
  }
  //
  ~pqUI() {
    delete this->Links;
  }

  void CreateProxy() {
    vtkSMProxyManager *pm = vtkSMProxy::GetProxyManager();
    SPHProxy = pm->NewProxy("meshless_helpers", "SPHManager");
    this->SPHProxy->UpdatePropertyInformation();
  }

  //
  bool ProxyCreated() { return this->SPHProxy!=NULL; }
  pqPropertyLinks            *Links;
  int                         SPHInitialized;
  vtkSmartPointer<vtkSMProxy> SPHProxy;
  vtkSmartPointer<vtkSMProxy> ActiveSourceProxy;
  int                         ActiveSourcePort;
  pqServer                   *ActiveServer;
};
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
pqSPHManagerPanel::pqSPHManagerPanel(QWidget* p) :
QDockWidget("SPH Manager", p)
{
  this->UI = new pqUI(this);
  this->UI->setupUi(this);


  //
  // Link paraview events to callbacks
  //
  pqServerManagerModel* smModel =
    pqApplicationCore::instance()->getServerManagerModel();

  this->connect(smModel, SIGNAL(serverAdded(pqServer*)),
    this, SLOT(onServerAdded(pqServer*)));

  this->connect(&pqActiveObjects::instance(),
    SIGNAL(serverChanged(pqServer*)),
    this, SLOT(onActiveServerChanged(pqServer*)));

  this->connect(smModel, SIGNAL(aboutToRemoveServer(pqServer *)),
    this, SLOT(StartRemovingServer(pqServer *)));

/*
  //
  // Button Group for Standalone/Server/Client buttons
  //
  this->SPHServerGroup = new QButtonGroup(this->UI);
  this->SPHServerGroup->addButton(this->UI->dsmIsStandalone);
  this->SPHServerGroup->addButton(this->UI->dsmIsServer);
  this->SPHServerGroup->addButton(this->UI->dsmIsClient);
  this->SPHServerGroup->setId(this->UI->dsmIsStandalone,0);
  this->SPHServerGroup->setId(this->UI->dsmIsServer,1);
  this->SPHServerGroup->setId(this->UI->dsmIsClient,2);
*/
  //
  this->LoadSettings();
}
//----------------------------------------------------------------------------
pqSPHManagerPanel::~pqSPHManagerPanel()
{
  this->SaveSettings();

//  if (this->SPHServerGroup) delete this->SPHServerGroup;
//  this->SPHServerGroup = NULL;
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::LoadSettings()
{
  pqSettings *settings = pqApplicationCore::instance()->settings();
/*
  settings->beginGroup("SPHManager");
  int size = settings->beginReadArray("Servers");
  if (size>0) {
    for (int i=0; i<size; ++i) {
      settings->setArrayIndex(i);
      QString server = settings->value("server").toString();
      if (this->UI->dsmServerName->findText(server) < 0) {
        this->UI->dsmServerName->addItem(server);
      }
    }
  }
  settings->endArray();
  // Active server
  this->UI->dsmServerName->setCurrentIndex(settings->value("Selected", 0).toInt());
  // Method
  this->UI->xdmfCommTypeComboBox->setCurrentIndex(settings->value("Communication", 0).toInt());
  // Port
  this->UI->xdmfCommPort->setValue(settings->value("Port", 0).toInt());
  // Client/Server mode
  int index = settings->value("ClientServerMode", 0).toInt();
  if (index!=-1) this->SPHServerGroup->buttons()[index]->click();
  // Description file type
  this->UI->xdmfFileTypeComboBox->setCurrentIndex(settings->value("DescriptionFileType", 0).toInt());
  // Description file path
  QString descFilePath = settings->value("DescriptionFilePath").toString();
  if(!descFilePath.isEmpty()) {
    this->UI->xdmfFilePathLineEdit->insert(descFilePath);
  }
  //
  settings->endGroup();
*/
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::SaveSettings()
{
  pqSettings *settings = pqApplicationCore::instance()->settings();
/*
  settings->beginGroup("SPHManager");
  // servers
  settings->beginWriteArray("Servers");
  for (int i=0; i<this->UI->dsmServerName->model()->rowCount(); i++) {
    settings->setArrayIndex(i);
    settings->setValue("server", this->UI->dsmServerName->itemText(i));
  }
  settings->endArray();
  // Active server
  settings->setValue("Selected", this->UI->dsmServerName->currentIndex());
  // Method
  settings->setValue("Communication", this->UI->xdmfCommTypeComboBox->currentIndex());
  // Port
  settings->setValue("Port", this->UI->xdmfCommPort->value());
  // Client/Server mode
  int index = this->SPHServerGroup->checkedId();
  settings->setValue("ClientServerMode", index);
  // Description file type
  settings->setValue("DescriptionFileType", this->UI->xdmfFileTypeComboBox->currentIndex());
  // Description file path
  settings->setValue("DescriptionFilePath", this->UI->xdmfFilePathLineEdit->text());
  //
  settings->endGroup();
*/
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::onServerAdded(pqServer *server)
{
  this->UI->ActiveServer = server;
}
//-----------------------------------------------------------------------------
void pqSPHManagerPanel::onActiveServerChanged(pqServer* server)
{
  this->UI->ActiveServer = server;
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::StartRemovingServer(pqServer *server)
{
  if (this->UI->ProxyCreated()) {
    this->UI->SPHProxy->InvokeCommand("DestroySPH");
    this->UI->SPHProxy = NULL;
    this->UI->SPHInitialized = 0;
  }
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::LinkServerManagerProperties()
{
  // TBD
}
//---------------------------------------------------------------------------
bool pqSPHManagerPanel::ProxyReady()
{
  if (!this->UI->ProxyCreated()) {
    this->UI->CreateProxy();
    this->LinkServerManagerProperties();
    return this->UI->ProxyCreated();
  }
  return true;
}
//---------------------------------------------------------------------------
bool pqSPHManagerPanel::SPHReady()
{
  if (!this->ProxyReady()) return 0;
  //
  if (!this->UI->SPHInitialized) {
    //
//    bool server = (this->UI->dsmIsServer->isChecked() || this->UI->dsmIsStandalone->isChecked());
//    pqSMAdaptor::setElementProperty(
//      this->UI->SPHProxy->GetProperty("DsmIsServer"), server);
    //
    //
    this->UI->SPHProxy->UpdateVTKObjects();
//    this->UI->SPHProxy->InvokeCommand("CreateSPH");
    this->UI->SPHInitialized = 1;
  }
  return this->UI->SPHInitialized;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
