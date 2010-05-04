#ifndef _pqSPHManagerPanel_h
#define _pqSPHManagerPanel_h

#include <QDockWidget>

#include <vtkSmartPointer.h>

class QCheckBox;
class QComboBox;
class QPushButton;
class pqServer;
class pqView;
class QTreeWidgetItem;
class QProgressDialog;
class QButtonGroup;

class vtkSMSourceProxy;
class vtkSMRepresentationProxy;

class pqSPHManagerPanel : public QDockWidget
{
  Q_OBJECT

public:
  /// constructor
  pqSPHManagerPanel(QWidget* p = NULL);
 ~pqSPHManagerPanel();

  bool ProxyReady();
  bool SPHReady();

signals:

public slots:
  void onServerAdded(pqServer *server);
  void onActiveServerChanged(pqServer *server);
  void StartRemovingServer(pqServer *server);

  void LoadSettings();
  void SaveSettings();

private slots:

protected:
  /// populate widgets with properties from the server manager
  virtual void LinkServerManagerProperties();

  class pqUI;
  pqUI* UI;

protected slots:

};

#endif

