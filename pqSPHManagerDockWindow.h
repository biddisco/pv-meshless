#ifndef _pqSPHManagerDockWindow_h
#define _pqSPHManagerDockWindow_h

#include <QDockWidget>
//#include "pqProxyPanel.h"

#include <vtkSmartPointer.h>

class QCheckBox;
class QComboBox;
class QPushButton;
class pqServer;
class pqView;
class QTreeWidgetItem;
class QProgressDialog;
class QButtonGroup;
class pqObjectPanel;

class vtkSMSourceProxy;
class vtkSMRepresentationProxy;
class pqSPHManagerPanel;

class pqSPHManagerDockWindow : public QDockWidget
{
  Q_OBJECT

public:
  /// constructor
  pqSPHManagerDockWindow(QWidget* p = NULL);
 ~pqSPHManagerDockWindow();

  bool ProxyReady();
  bool SPHReady();

  ///////////////////////////////////
  // From pqObjectInspectorWidget
  ///////////////////////////////////
signals:
  /// emitted before accept.
  void preaccept();
  /// emitted on accept() after preaccept() but before postaccept()/
  void accepted();
  ///emitted after accept;
  void postaccept();

  /// emitted before reject.
  void prereject();
  /// emitted after reject.
  void postreject();

  /// emitted when render module is changed
  void viewChanged(pqView*);

  void helpRequested(const QString& proxyType);

  void canAccept();

public slots:

  void updateAcceptState();

  /// accept the changes made to the properties
  /// changes will be propogated down to the server manager
  void accept();

  /// reset the changes made
  /// editor will query properties from the server manager
  void reset();

  /// Updates the accept/reset button state.
  void canAccept(bool status);


  ///////////////////////////////////
  // From pqObjectPanel
  ///////////////////////////////////
  /// Fires modified
  virtual void setModified();
 
private slots:

signals:
  void modified();

protected:
  /// populate widgets with properties from the server manager
  virtual void showEvent( QShowEvent * event );

  class pqUI;
  pqUI* UI;

protected slots:

protected:
  QPushButton   *AcceptButton;
  QPushButton   *ResetButton;
  pqSPHManagerPanel *SPHManagerPanel;

};

#endif

