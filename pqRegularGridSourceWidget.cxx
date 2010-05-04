/*=========================================================================

   Program:   ParaQ
   Module:    $RCSfile: pqRegularGridSourceWidget.cxx,v $

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaQ is a free software; you can redistribute it and/or modify it
   under the terms of the ParaQ license version 1.1. 

   See License_v1.1.txt for the full ParaQ license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include "pqApplicationCore.h"
#include "pqRegularGridSourceWidget.h"
#include "pqPropertyLinks.h"
#include "pqSMSignalAdaptors.h"
#include "pqRenderView.h"
#include "pq3DWidgetFactory.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqSMAdaptor.h"

#include <QDoubleValidator>
#include "ui_pqRegularGridSourceWidget.h"

#include "vtkSMPropertyLink.h"
#include <vtkSMProxyProperty.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMNewWidgetRepresentationProxy.h>
#include <vtkSMDoubleVectorProperty.h>
#include <vtkSMRenderViewProxy.h>
#include "vtkSMPropertyHelper.h"
#include <vtkMath.h>
#include <vtkCamera.h>
#include <vtkBoundingBox.h>
#include <vtkBox.h>
#include <vtkTransform.h>

#define PVBOXWIDGET_LINK(ui, smproperty, index)\
{\
  this->Implementation->Links.addPropertyLink(\
    this->Implementation->ui, "text2",\
    SIGNAL(textChanged(const QString&)),\
    widget, widget->GetProperty(smproperty), index);\
}

#define PVBOXWIDGET_LINK2(ui, smproperty)\
{\
  this->Implementation->Links.addPropertyLink(\
    this->Implementation->ui, "value",\
    SIGNAL(valueChanged(int)), \
    widget, widget->GetProperty(smproperty),0);\
}

#define PVBOXWIDGET_TRIGGER_UPDATE(ui)  \
  QObject::connect(this->Implementation->ui,\
    SIGNAL(editingFinished()),\
    this, SLOT(UpdatePlanePoints()), Qt::QueuedConnection);

/////////////////////////////////////////////////////////////////////////
// pqRegularGridSourceWidget::pqImplementation

class pqRegularGridSourceWidget::pqImplementation : public Ui::pqRegularGridSourceWidget
{
public:
  pqImplementation() : Ui::pqRegularGridSourceWidget()
  {
//    this->Links.setUseUncheckedProperties(false);
//    this->Links.setAutoUpdateVTKObjects(true);
  } 

  ~pqImplementation()
  {
    this->Links.removeAllPropertyLinks();
  }
    
  /// Maps Qt widgets to the 3D widget
  pqPropertyLinks Links;
  /// links between SM properties
  QList<vtkSmartPointer<vtkSMPropertyLink> > PropertyLinks;

};

/////////////////////////////////////////////////////////////////////////
// pqRegularGridSourceWidget

//-----------------------------------------------------------------------------
pqRegularGridSourceWidget::pqRegularGridSourceWidget(vtkSMProxy* refProxy, vtkSMProxy* pxy, QWidget* _parent) :
  Superclass(refProxy, pxy, _parent)
{
  this->Implementation = new pqImplementation();
  this->Implementation->setupUi(this);
  this->Implementation->show3DWidget->setChecked(this->widgetVisible());  

  // Setup validators for all line edits.
  QDoubleValidator* validator = new QDoubleValidator(this);
  this->Implementation->positionX->setValidator(validator);
  this->Implementation->positionY->setValidator(validator);
  this->Implementation->positionZ->setValidator(validator);
  this->Implementation->scaleX->setValidator(validator);
  this->Implementation->scaleY->setValidator(validator);
  this->Implementation->scaleZ->setValidator(validator);
  this->Implementation->rotationX->setValidator(validator);
  this->Implementation->rotationY->setValidator(validator);
  this->Implementation->rotationZ->setValidator(validator);
  this->Implementation->normalX->setValidator(validator);
  this->Implementation->normalY->setValidator(validator);
  this->Implementation->normalZ->setValidator(validator);

  PVBOXWIDGET_TRIGGER_UPDATE(positionX);
  PVBOXWIDGET_TRIGGER_UPDATE(positionY);
  PVBOXWIDGET_TRIGGER_UPDATE(positionZ);

  PVBOXWIDGET_TRIGGER_UPDATE(scaleX);
  PVBOXWIDGET_TRIGGER_UPDATE(scaleY);
  PVBOXWIDGET_TRIGGER_UPDATE(scaleZ);

  PVBOXWIDGET_TRIGGER_UPDATE(rotationX);
  PVBOXWIDGET_TRIGGER_UPDATE(rotationY);
  PVBOXWIDGET_TRIGGER_UPDATE(rotationZ);

  PVBOXWIDGET_TRIGGER_UPDATE(normalX);
  PVBOXWIDGET_TRIGGER_UPDATE(normalY);
  PVBOXWIDGET_TRIGGER_UPDATE(normalZ);

  PVBOXWIDGET_TRIGGER_UPDATE(resolutionX);
  PVBOXWIDGET_TRIGGER_UPDATE(resolutionY);
  PVBOXWIDGET_TRIGGER_UPDATE(resolutionZ);

  QObject::connect(this->Implementation->show3DWidget,
    SIGNAL(toggled(bool)), this, SLOT(setWidgetVisible(bool)));

  QObject::connect(this, SIGNAL(widgetVisibilityChanged(bool)),
    this, SLOT(onWidgetVisibilityChanged(bool)));

  QObject::connect(this->Implementation->OBB,
    SIGNAL(clicked()), this, SLOT(onComputeOBB()));

  QObject::connect(&this->Implementation->Links, 
    SIGNAL(qtWidgetChanged()), this, SLOT(setModified()));

  connect(this->Implementation->useXNormal,
    SIGNAL(clicked()), this, SLOT(onUseXNormal()));
  connect(this->Implementation->useYNormal,
    SIGNAL(clicked()), this, SLOT(onUseYNormal()));
  connect(this->Implementation->useZNormal,
    SIGNAL(clicked()), this, SLOT(onUseZNormal()));
  connect(this->Implementation->useCameraNormal,
    SIGNAL(clicked()), this, SLOT(onUseCameraNormal()));
  connect(this->Implementation->resetBounds,
    SIGNAL(clicked()), this, SLOT(resetBounds()));
  connect(this->Implementation->WidgetMode,
    SIGNAL(currentIndexChanged(int)), this, SLOT(onWidgetModeChanged(int)));  

  pqServerManagerModel* smmodel =
    pqApplicationCore::instance()->getServerManagerModel();
  //
  this->createWidget(smmodel->findServer(refProxy->GetConnectionID()));

  //
  //
  //
  vtkSMNewWidgetRepresentationProxy *widget = this->getWidgetProxy();
  vtkSMProxy *reg_grid_source = this->getControlledProxy();

  // Server Manager Property Links
  vtkSMPropertyLink* link = vtkSMPropertyLink::New();
  link->AddLinkedProperty(reg_grid_source, "OriginInfo", vtkSMPropertyLink::INPUT);
  link->AddLinkedProperty(widget, "OriginInfo", vtkSMPropertyLink::OUTPUT);
  this->Implementation->PropertyLinks.push_back(link);
  link->Delete();

}
//-----------------------------------------------------------------------------
pqRegularGridSourceWidget::~pqRegularGridSourceWidget()
{
  delete this->Implementation;
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::createWidget(pqServer* server)
{
  vtkSMNewWidgetRepresentationProxy* widget = 
    pqApplicationCore::instance()->get3DWidgetFactory()->
    get3DWidget("CustomBoxWidgetRepresentation", server);
  this->setWidgetProxy(widget);
  widget->UpdateVTKObjects();
  widget->UpdatePropertyInformation();

  // Now bind the GUI controls to the 3D widget.

  PVBOXWIDGET_LINK(positionX, "Position", 0);
  PVBOXWIDGET_LINK(positionY, "Position", 1);
  PVBOXWIDGET_LINK(positionZ, "Position", 2);

  PVBOXWIDGET_LINK(rotationX, "Rotation", 0);
  PVBOXWIDGET_LINK(rotationY, "Rotation", 1);
  PVBOXWIDGET_LINK(rotationZ, "Rotation", 2);

  PVBOXWIDGET_LINK(scaleX, "Scale", 0);
  PVBOXWIDGET_LINK(scaleY, "Scale", 1);
  PVBOXWIDGET_LINK(scaleZ, "Scale", 2);

  PVBOXWIDGET_LINK(resolutionX, "Resolution", 0);
  PVBOXWIDGET_LINK(resolutionY, "Resolution", 1);
  PVBOXWIDGET_LINK(resolutionZ, "Resolution", 2);

//  PVBOXWIDGET_LINK(normalX, "NormalInfo", 0);
//  PVBOXWIDGET_LINK(normalY, "NormalInfo", 1);
//  PVBOXWIDGET_LINK(normalZ, "NormalInfo", 2);

}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::accept()
{
  this->Superclass::accept();
  this->ComputePlanePoints();
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onWidgetVisibilityChanged(bool visible)
{
  this->Implementation->show3DWidget->blockSignals(true);
  this->Implementation->show3DWidget->setChecked(visible);
  this->Implementation->show3DWidget->blockSignals(false);
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::resetBounds(double bounds[6])
{
  vtkBoundingBox box(bounds);
  double c[3], o[3];
  box.GetCenter(c);
//  box.GetMinPoint(o[0], o[1], o[2]);
  for (int i=0; i<3; i++) { o[i] = 0; } // /*o[i]-*/c[i]; }

  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  vtkSMPropertyHelper(widget, "PlaceWidget").Set(bounds, 6);
  widget->UpdateVTKObjects();
  widget->UpdatePropertyInformation();
  //
  // reset the Transform 
  double rotation[3] = { 0.0, 0.0, 0.0 };
  double scale[3] = { 1.0, 1.0, 1.0 };
  vtkSMPropertyHelper(widget, "Rotation").Set(rotation, 3);
  vtkSMPropertyHelper(widget, "Scale").Set(scale, 3);
  vtkSMPropertyHelper(widget, "Position").Set(o, 3);
  widget->UpdateVTKObjects();
  widget->UpdatePropertyInformation();
  this->setModified();
  this->render();
}
//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::resetBounds()
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  this->SizeFactor = 1.2;
  if (!widget)
    {
      return;
    }
  vtkSMPropertyHelper(widget, "PlaceFactor").Get(&this->SizeFactor, 1);
  if (!this->getReferenceInputBounds(this->InputBounds))
    {
    this->InputBounds[0] = -0.5*this->SizeFactor;
    this->InputBounds[1] =  0.5*this->SizeFactor;
    this->InputBounds[2] = -0.5*this->SizeFactor;
    this->InputBounds[3] =  0.5*this->SizeFactor;
    this->InputBounds[4] = -0.5*this->SizeFactor;
    this->InputBounds[5] =  0.5*this->SizeFactor;
    }
  this->resetBounds(this->InputBounds);
}

//
enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::UpdatePlanePoints()
{
  this->ComputePlanePoints();
  this->render();
}
//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::ComputePlanePoints()
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if (!widget)
  {
    return;
  }
    
  double plane_origin[3];
  double point1[3], point2[3], point3[3];

  vtkSMPropertyHelper(widget, "OriginInfo").Get(plane_origin, 3);
//  vtkSMPropertyHelper(widget, "NormalInfo").Get(plane_normal, 3);   
  vtkSMPropertyHelper(widget, "Point1Info").Get(point1, 3);
  vtkSMPropertyHelper(widget, "Point2Info").Get(point2, 3);
  vtkSMPropertyHelper(widget, "Point3Info").Get(point3, 3);
  
  vtkSMDoubleVectorProperty* const o1 = vtkSMDoubleVectorProperty::SafeDownCast(
      this->getControlledProxy()->GetProperty("Origin"));
  if(o1)
    {
    o1->SetElements(plane_origin);
    }
  this->getControlledProxy()->UpdateProperty("Origin");

  vtkSMDoubleVectorProperty* const p1 = vtkSMDoubleVectorProperty::SafeDownCast(
      this->getControlledProxy()->GetProperty("Point1"));
  if(p1)
    {
    p1->SetElements(point1);
    }
  this->getControlledProxy()->UpdateProperty("Point1");

  vtkSMDoubleVectorProperty* const p2 = vtkSMDoubleVectorProperty::SafeDownCast(
      this->getControlledProxy()->GetProperty("Point2"));
  if(p2)
    {
    p2->SetElements(point2);
    }
  this->getControlledProxy()->UpdateProperty("Point2");

  int planemode = 0;
  vtkSMPropertyHelper(widget, "PlaneMode").Get(&planemode, 1);
  //
  vtkSMDoubleVectorProperty* const p3 = 
    vtkSMDoubleVectorProperty::SafeDownCast(
      this->getControlledProxy()->GetProperty("Point3"));
  if(p3)
    {
    if (planemode) 
      {
      p3->SetElements(plane_origin);
      }
    else
      {
      p3->SetElements(point3);
      }
    this->getControlledProxy()->UpdateProperty("Point3");
    }
  widget->UpdateVTKObjects();

}
//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onUseXNormal()
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if(widget)
    {
      widget->InvokeCommand("OrientX");
      widget->UpdatePropertyInformation();
      widget->UpdateVTKObjects();
      this->render();
      this->setModified();
    }
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onUseYNormal()
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if(widget)
    {
      widget->InvokeCommand("OrientY");
      widget->UpdatePropertyInformation();
      widget->UpdateVTKObjects();
      this->render();
      this->setModified();
    }
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onUseZNormal()
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if(widget)
    {
      widget->InvokeCommand("OrientZ");
      widget->UpdatePropertyInformation();
      widget->UpdateVTKObjects();
      this->render();
      this->setModified();
    }
}

//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onUseCameraNormal()
{
  return; // disabled until fixed
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if (widget)
    {
    pqRenderView* renView = qobject_cast<pqRenderView*>(this->renderView());
    if (vtkCamera* const camera = renView? 
      renView->getRenderViewProxy()->GetActiveCamera() : 0)
      {
      double camera_normal[3];
      camera->GetViewPlaneNormal(camera_normal);
      camera_normal[0] = -camera_normal[0];
      camera_normal[1] = -camera_normal[1];
      camera_normal[2] = -camera_normal[2];
      vtkSMPropertyHelper(widget, "Normal").Set(camera_normal, 3);
      widget->UpdateVTKObjects();
      this->render();
      this->setModified();
      }
    }
}
//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onComputeOBB()
{
  double boxdata[12];

  vtkSMNewWidgetRepresentationProxy *widget = this->getWidgetProxy();
  vtkSMProxy *reg_grid_source = this->getControlledProxy();

  // update source object to compute OBB
  reg_grid_source->InvokeCommand("ComputeOBB");
  // properties which are information must be extracted again
  reg_grid_source->UpdatePropertyInformation();
  // Get the params of the new box from the source where they were computed
  vtkSMPropertyHelper(reg_grid_source, "OriginInfo").Get(&boxdata[0], 3);
  vtkSMPropertyHelper(reg_grid_source, "Point1Info").Get(&boxdata[3], 3);
  vtkSMPropertyHelper(reg_grid_source, "Point2Info").Get(&boxdata[6], 3);
  vtkSMPropertyHelper(reg_grid_source, "Point3Info").Get(&boxdata[9], 3);

  // copy OBB to the widget representation which we see on screen
  vtkSMDoubleVectorProperty* const obb = 
    vtkSMDoubleVectorProperty::SafeDownCast(widget->GetProperty("SetOBB"));
  obb->SetElements(boxdata);
  obb->Modified();
  widget->UpdateProperty("SetOBB");
  //
  // make sure new data is propagated to all widget properties
  //
  widget->UpdatePropertyInformation();

  double rotation[3] = { 0.0, 0.0, 0.0 };
  double scale[3] = { 1.0, 1.0, 1.0 };
  double o[3] = { 0.0, 0.0, 0.0 };
  vtkSMPropertyHelper(widget, "Rotation").Get(rotation, 3);
  vtkSMPropertyHelper(widget, "Scale").Get(scale, 3);
  vtkSMPropertyHelper(widget, "Position").Get(o, 3);

  //
  this->render();
  this->setModified();
}
//-----------------------------------------------------------------------------
void pqRegularGridSourceWidget::onWidgetModeChanged(int index)
{
  vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  if (index==0) {
    int planemode = 1;
    vtkSMPropertyHelper(widget, "PlaneMode").Set(&planemode, 1);
  } 
  else {
    int planemode = 0;
    vtkSMPropertyHelper(widget, "PlaneMode").Set(&planemode, 1);
  }
  widget->UpdateProperty("PlaneMode");
  widget->UpdatePropertyInformation();
  this->render();
  this->setModified();
}
