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
#include "vtkSMPropertyHelper.h"
#include "vtkProcessModule.h"

// ParaView includes
#include "pqActiveServer.h"
#include "pqApplicationCore.h"
#include "pqAutoGeneratedObjectPanel.h"
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
#include "pqPropertyManager.h"
#include "pqNamedWidgets.h"
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
  }
  //
  ~pqUI() {
    delete this->Links;
  }
  //
  pqPropertyLinks            *Links;
};
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
pqSPHManagerPanel::pqSPHManagerPanel(pqProxy* proxy, QWidget* p) : 
  pqNamedObjectPanel(proxy, p) 
{
  this->UI = new pqUI(this);
  this->UI->setupUi(this);

  //
  // Link Gui events
  //
  connect(this->UI->InterpolationMethod,
    SIGNAL(currentIndexChanged(int)), this, SLOT(onPointMethodChanged(int)));

  // inherited from pqNamedObjectPanel
  this->linkServerManagerProperties();
  
  // Link this so that when reset is clicked, the combo goes back to its old state
  pqNamedWidgets::linkObject(this->UI->InterpolationMethod, proxy->getProxy(), 
      "InterpolationMethod", this->propertyManager());

  // after changes are accepted, we must release hold of the DSM
  this->connect(this, SIGNAL(onaccept()), this, SLOT(onAccept()));
}
//----------------------------------------------------------------------------
pqSPHManagerPanel::~pqSPHManagerPanel()
{
}
//-----------------------------------------------------------------------------
void pqSPHManagerPanel::onPointMethodChanged(int mode) 
{
  if (mode==0) {
    this->UI->SPHBox->setEnabled(1);
  }
  else {
    this->UI->SPHBox->setEnabled(0);
  }
}
//-----------------------------------------------------------------------------
void pqSPHManagerPanel::onLimitChanged(bool checked)
{
  if (checked) {
    this->UI->MaximumNeighbours->setEnabled(1);
    this->UI->MaximumSearchRadius->setEnabled(0);
  }
  else {
    this->UI->MaximumNeighbours->setEnabled(0);
    if (this->UI->InterpolationMethod->currentIndex()==1) {
      this->UI->MaximumSearchRadius->setEnabled(1);
    }
  }
}
//-----------------------------------------------------------------------------
void pqSPHManagerPanel::LoadSettings()
{
  pqSettings *settings = pqApplicationCore::instance()->settings();
  settings->beginGroup("SPHManager");
  // InterpolationMethod
  this->UI->InterpolationMethod->setCurrentIndex(settings->value("InterpolationMethod", 0).toInt());
  // KernelType
  this->UI->KernelType->setCurrentIndex(settings->value("KernelType", 0).toInt());
  // KernelDimension
  this->UI->KernelDimension->setCurrentIndex(settings->value("KernelDimension", 0).toInt());
  // HCoefficient
  this->UI->HCoefficient->setText(settings->value("HCoefficient", 1.5).toString());
  // DefaultDensity
  this->UI->DefaultDensity->setText(settings->value("DefaultDensity", 1000.0).toString());
  // DefaultParticleSideLength
  this->UI->DefaultParticleSideLength->setText(settings->value("DefaultParticleSideLength", 0.18333).toString());
  // Max Neighbours
  this->UI->MaximumNeighbours->setText(settings->value("MaximumNeighbours", 64).toString()); 
  // MaxSearchRadius
  this->UI->MaximumSearchRadius->setText(settings->value("MaximumSearchRadius", 0.0).toString());
  settings->endGroup();
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::SaveSettings()
{
  pqSettings *settings = pqApplicationCore::instance()->settings();
  settings->beginGroup("SPHManager");
  // InterpolationMethod
  settings->setValue("InterpolationMethod", this->UI->InterpolationMethod->currentIndex());
  // KernelType
  settings->setValue("KernelType", this->UI->KernelType->currentIndex());
  // KernelDimension
  settings->setValue("KernelDimension", this->UI->KernelDimension->currentIndex());
  // HCoefficient
  settings->setValue("HCoefficient", this->UI->HCoefficient->text());
  // DefaultDensity
  settings->setValue("DefaultDensity", this->UI->DefaultDensity->text());
  // DefaultParticleSideLength
  settings->setValue("DefaultParticleSideLength", this->UI->DefaultParticleSideLength->text());
  // Max Neighbours
  settings->setValue("MaximumNeighbours", this->UI->MaximumNeighbours->text());
  // MaxSearchRadius
  settings->setValue("MaximumSearchRadius", this->UI->MaximumSearchRadius->text());
  settings->endGroup();
}
//-----------------------------------------------------------------------------
void pqSPHManagerPanel::onAccept()
{
  this->ModifyGUIFilters();
}
//----------------------------------------------------------------------------
void pqSPHManagerPanel::ModifyGUIFilters()
{
    //
    // Set status of registered pipeline source to unmodified 
    //

    QList<pqPipelineSource*> sources = 
      pqApplicationCore::instance()->getServerManagerModel()->findItems<pqPipelineSource*>((pqServer*)(0));

    for (QList<pqPipelineSource*>::Iterator it=sources.begin(); it!=sources.end(); it++) {
      pqPipelineSource *source = *it;
      std::string xmlName = source->getProxy()->GetXMLName();
      if (xmlName.find("vtkSPH")!=xmlName.npos) {
//        source->getProxy()->MarkAllPropertiesAsModified();
        vtkSMPropertyHelper modified(source->getProxy(), "ModifiedNumber");
        modified.UpdateValueFromServer();
        modified.Set(modified.GetAsInt()+1);
        source->setModifiedState(pqProxy::MODIFIED);
       //          InvokeCommand("Modified");

//        source->getProxy()->GetProperty("SPHManager")->Modified();
      }
    }
}
