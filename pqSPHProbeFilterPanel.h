
#ifndef _pqProbeCustomFilterPanel_h
#define _pqProbeCustomFilterPanel_h

#include "pqNamedObjectPanel.h"
#include "pqObjectPanelInterface.h"

class pqProbeCustomFilterPanel : public pqNamedObjectPanel
{
  Q_OBJECT
  typedef pqNamedObjectPanel Superclass;
public:
  /// constructor
   pqProbeCustomFilterPanel(pqProxy* proxy, QWidget* p);
  /// destructor
  ~pqProbeCustomFilterPanel();

  
private slots:
  void onModeChanged(int mode);
  void onPointMethodChanged(int mode);
  void onLimitChanged(bool checked);

private:
  class pqImplementation;
  pqImplementation* const Implementation;

};

#endif

