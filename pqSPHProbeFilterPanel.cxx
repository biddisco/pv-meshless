#include "ui_pqProbeCustomFilterPanel.h"

// Qt includes
#include <QTreeWidget>
#include <QVariant>
#include <QLabel>
#include <QComboBox>
#include <QTableWidget>
#include <QLayout>

// VTK includes

// ParaView Server Manager includes
#include "vtkSMProxyManager.h"
#include "vtkSMSourceProxy.h"
#include "vtkSMStringVectorProperty.h"
#include "vtkSMArraySelectionDomain.h"
#include "vtkSMProxy.h"
#include "vtkSMPropertyHelper.h"

// ParaView includes
#include "pqProxy.h"
#include "pqSMAdaptor.h"
#include "pqTreeWidgetCheckHelper.h"
#include "pqTreeWidgetItemObject.h"
#include "pqPropertyLinks.h"
#include "pqRegularGridSourceWidget.h"

#include "pqProbeCustomFilterPanel.h"
//-----------------------------------------------------------------------------
class pqProbeCustomFilterPanel::pqImplementation  : public QObject, public Ui::pqProbeCustomFilterPanel
{
public:
  pqImplementation() : RegularGridSourceWidget(NULL)
  {
  }

  ~pqImplementation()
  {
    delete this->RegularGridSourceWidget;
  }

  /// Provides a UI for managing a RegularGridSource
  pqRegularGridSourceWidget *RegularGridSourceWidget;
  /// Provides the remaining Qt controls for the panel
//  Ui::pqProbeCustomFilterPanel UI;
};
//-----------------------------------------------------------------------------
pqProbeCustomFilterPanel::pqProbeCustomFilterPanel(pqProxy* pxy, QWidget* p)
  : Superclass(pxy, p), Implementation(new pqImplementation())
{
  // Setup the Qt control inside this panel
  this->Implementation->setupUi(this);

  // bind server-manager properties to GUI widgets
  // names must be same for widgets as properties in XML
  this->linkServerManagerProperties();

  connect(this->Implementation->InterpolationMode,
    SIGNAL(currentIndexChanged(int)), this, SLOT(onModeChanged(int)));

  connect(this->Implementation->InterpolationMethod,
    SIGNAL(currentIndexChanged(int)), this, SLOT(onPointMethodChanged(int)));

  connect(this->Implementation->limit1,
    SIGNAL(toggled(bool)), this, SLOT(onLimitChanged(bool)));

  int meshstatus = 0;
  vtkSMPropertyHelper(this->proxy(), "MeshStatus").Get(&meshstatus, 1);
  if (meshstatus) {
    this->Implementation->InterpolationMode->setCurrentIndex(0);
    this->onModeChanged(0);
  }
  else {
    this->Implementation->InterpolationMode->setCurrentIndex(1);
    this->onModeChanged(1);
  }
  this->onLimitChanged(this->Implementation->limit1->isChecked());
}
//-----------------------------------------------------------------------------
pqProbeCustomFilterPanel::~pqProbeCustomFilterPanel()
{
}
//-----------------------------------------------------------------------------
void pqProbeCustomFilterPanel::onModeChanged(int mode) 
{
  if (mode==0) {
    this->Implementation->modestack->setCurrentIndex(0);
  }
  else {
    this->Implementation->modestack->setCurrentIndex(1);
  }
}
//-----------------------------------------------------------------------------
void pqProbeCustomFilterPanel::onPointMethodChanged(int mode) 
{
  if (mode==0) {
    this->Implementation->methodstack->setCurrentIndex(0);
  }
  else {
    this->Implementation->methodstack->setCurrentIndex(1);
  }
}
//-----------------------------------------------------------------------------
void pqProbeCustomFilterPanel::onLimitChanged(bool checked)
{
  if (checked) {
    this->Implementation->MaximumNeighbours->setEnabled(1);
    this->Implementation->MaximumRadius->setEnabled(0);
    vtkSMPropertyHelper(this->proxy(), "LimitSearchByNeighbourCount").Set(1);
  }
  else {
    this->Implementation->MaximumNeighbours->setEnabled(0);
    this->Implementation->MaximumRadius->setEnabled(1);
    vtkSMPropertyHelper(this->proxy(), "LimitSearchByNeighbourCount").Set(0);
  }
  this->setModified();
}

