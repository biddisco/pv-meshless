#ifndef _pqSPHManagerPanel_h
#define _pqSPHManagerPanel_h

#include "pqProxy.h"
#include "pqNamedObjectPanel.h"

class pqSPHManagerPanel : public pqNamedObjectPanel 
{
  Q_OBJECT

public:
  /// constructor
  pqSPHManagerPanel(pqProxy* proxy, QWidget* p);
 ~pqSPHManagerPanel();

  void LoadSettings();
  void SaveSettings();

private slots:
  void onPointMethodChanged(int mode);
  void onLimitChanged(bool checked);

protected:

  class pqUI;
  pqUI* UI;

protected slots:

};

#endif

