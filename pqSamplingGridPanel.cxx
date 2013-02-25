/*=========================================================================

   Program:   ParaQ
   Module:    $RCSfile: pqSamplingGridPanel.cxx,v $

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
#include "pqSamplingGridPanel.h"
#include "pqPropertyLinks.h"
#include "pqSMSignalAdaptors.h"
#include "pqRenderView.h"
#include "pq3DWidgetFactory.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqActiveObjects.h"
#include "pqSMAdaptor.h"
#include "pqBoxWidget.h"
#include "pqImplicitPlaneWidget.h"

#include <QDoubleValidator>
#include <QSpacerItem>
#include "ui_pqSamplingGridPanel.h"

#include "vtkSMPropertyLink.h"
#include "vtkSMProxyManager.h"
#include "vtkSMInputProperty.h"
#include <vtkSMProxyProperty.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMNewWidgetRepresentationProxy.h>
#include <vtkPVDataInformation.h>
#include <vtkSMDoubleVectorProperty.h>
#include <vtkSMRenderViewProxy.h>
#include "vtkSMPropertyHelper.h"
#include <vtkMath.h>
#include <vtkCamera.h>
#include <vtkBoundingBox.h>
#include <vtkBox.h>
#include <vtkTransform.h>
#include <vtkCutter.h>
#include <vtkImageData.h>
#include <vtkPlane.h>

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
// pqSamplingGridPanel::pqImplementation

class pqSamplingGridPanel::pqImplementation : public Ui::pqSamplingGridPanel
{
public:
  pqImplementation() : Ui::pqSamplingGridPanel()
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
// pqSamplingGridPanel

//-----------------------------------------------------------------------------
pqSamplingGridPanel::pqSamplingGridPanel(pqProxy* proxy, QWidget* _parent) :
  Superclass(proxy, _parent)
{
  this->Implementation = new pqImplementation();
  
  // Create a frame and copy our custom gui into it
  QFrame *frame = new QFrame();
  frame->setObjectName(QString::fromUtf8("frame"));
  frame->setFrameShape(QFrame::StyledPanel);
  frame->setFrameShadow(QFrame::Raised);
  frame->setFrameStyle(QFrame::NoFrame);
  this->Implementation->setupUi(frame);

  // add custom panel to the existing auto-generated stuff
  int rows = this->PanelLayout->rowCount();
  int cols = this->PanelLayout->columnCount();
  QVBoxLayout* subLayout = new QVBoxLayout();
  subLayout->addWidget(frame, 1);
  subLayout->setMargin(0);
  subLayout->setSpacing(0);
  this->PanelLayout->addLayout(subLayout, rows-1, 0, -1, -1);

  connect(this->Implementation->WidgetMode,
    SIGNAL(currentIndexChanged(int)), this, SLOT(onWidgetModeChanged(int)));  

/*
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
*/

  QObject::connect(&this->Implementation->Links, 
    SIGNAL(qtWidgetChanged()), this, SLOT(setModified()));

  //connect(this->Implementation->useCameraNormal,
  //  SIGNAL(clicked()), this, SLOT(onUseCameraNormal()));

  pqServerManagerModel* smmodel =
    pqApplicationCore::instance()->getServerManagerModel();
  //
  this->currentMode = 0;
  this->createWidget(pqActiveObjects::instance().activeServer());

  connect(this->boxWidget,
    SIGNAL(modified()), this, SLOT(onWidgetModified()));  

  connect(this->planeWidget,
    SIGNAL(modified()), this, SLOT(onWidgetModified()));  

/*
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
  link->FastDelete();
*/
  //
  // Start in Box mode by default
  //
  this->Implementation->stackedWidget->blockSignals(true);
  this->Implementation->stackedWidget->setCurrentIndex(1);
  this->onWidgetModeChanged(1);
  this->Implementation->stackedWidget->blockSignals(false);
}
//-----------------------------------------------------------------------------
pqSamplingGridPanel::~pqSamplingGridPanel()
{
  delete this->Implementation;
}

//-----------------------------------------------------------------------------
void pqSamplingGridPanel::setView(pqView* pqview)
{ 
  Superclass::setView(pqview);
  this->boxWidget->setView(pqview);
  this->planeWidget->setView(pqview);
}

//-----------------------------------------------------------------------------
void pqSamplingGridPanel::createWidget(pqServer* server)
{
  vtkSMProxy *proxy = this->referenceProxy()->getProxy();
  //
  this->boxWidget = new pqBoxWidget(proxy,proxy);
  this->Implementation->boxwidgetholder->addWidget(this->boxWidget);
  QSpacerItem *vSpacer1 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);
  this->Implementation->boxwidgetholder->addItem(vSpacer1);
  this->boxWidget->resetBounds();
  //
  this->planeWidget = new pqImplicitPlaneWidget(proxy,proxy);
  this->Implementation->planewidgetholder->addWidget(this->planeWidget);
  QSpacerItem *vSpacer2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);
  this->Implementation->planewidgetholder->addItem(vSpacer2);
  this->planeWidget->resetBounds();

  // Now bind the GUI controls to the 3D widget.
/*
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
*/
//  PVBOXWIDGET_LINK(normalX, "NormalInfo", 0);
//  PVBOXWIDGET_LINK(normalY, "NormalInfo", 1);
//  PVBOXWIDGET_LINK(normalZ, "NormalInfo", 2);

}

//-----------------------------------------------------------------------------
void pqSamplingGridPanel::accept()
{
  this->onWidgetModified();
  //
  this->Superclass::accept();
  this->ComputePlanePoints();
}

//-----------------------------------------------------------------------------
void pqSamplingGridPanel::onWidgetModified()
{
  double Pos[3],Rot[3],Sca[3],Nor[3];
  double bounds[6],lengths[3];
  if (!this->getReferenceInputBounds(bounds)) {
    bounds[0] = bounds[2] = bounds[4] = -0.25;
    bounds[1] = bounds[3] = bounds[5] =  0.25;
  }
  vtkBoundingBox box(bounds);
  box.GetLengths(lengths);

  // Box mode
  if (this->currentMode == 1) {
    this->boxWidget->getWidgetProxy()->UpdatePropertyInformation();
    //
    vtkSMPropertyHelper pos(this->boxWidget->getWidgetProxy(), "PositionInfo");
    pos.UpdateValueFromServer();
    pos.Get(Pos,3);
    //
    vtkSMPropertyHelper rot(this->boxWidget->getWidgetProxy(), "RotationInfo");
    rot.UpdateValueFromServer();
    rot.Get(Rot,3);
    //
    vtkSMPropertyHelper sca(this->boxWidget->getWidgetProxy(), "ScaleInfo");
    sca.UpdateValueFromServer();
    sca.Get(Sca,3);
    //
    vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
    trans->Identity();
    trans->Translate(Pos);
    trans->RotateZ(Rot[2]);
    trans->RotateX(Rot[0]);
    trans->RotateY(Rot[1]);
    trans->Scale(Sca);
    trans->Translate(box.GetMinPoint());
    //
    double *p0 = trans->TransformDoublePoint(0,0,0);
    vtkSMPropertyHelper pp0(this->ReferenceProxy->getProxy(), "Origin");
    pp0.Set(p0,3);
    //
    double *p1 = trans->TransformDoublePoint(lengths[0],0,0);
    vtkSMPropertyHelper pp1(this->ReferenceProxy->getProxy(), "Point1");
    pp1.Set(p1,3);
    //
    double *p2 = trans->TransformDoublePoint(0,lengths[1],0);
    vtkSMPropertyHelper pp2(this->ReferenceProxy->getProxy(), "Point2");
    pp2.Set(p2,3);
    //
    double *p3 = trans->TransformDoublePoint(0,0,lengths[2]);
    vtkSMPropertyHelper pp3(this->ReferenceProxy->getProxy(), "Point3");
    pp3.Set(p3,3);
  }
  // Plane mode
  else {
    this->planeWidget->getWidgetProxy()->UpdatePropertyInformation();
    //
    vtkSMPropertyHelper pos(this->planeWidget->getWidgetProxy(), "OriginInfo");
    pos.UpdateValueFromServer();
    pos.Get(Pos,3);
    //
    vtkSMPropertyHelper nor(this->planeWidget->getWidgetProxy(), "NormalInfo");
    nor.UpdateValueFromServer();
    nor.Get(Nor,3);
    //
//    vtkSMPropertyHelper pf(this->planeWidget->getWidgetProxy(), "PlaceFactor");
//    pf.UpdateValueFromServer();
//    pf.Get(Sca,1);   
    //
    vtkSmartPointer<vtkImageData> Box = vtkSmartPointer<vtkImageData>::New();
    Box->SetDimensions(2,2,2);
    vtkSmartPointer<vtkCutter> Cutter = vtkSmartPointer<vtkCutter>::New();
    Cutter->SetInputData(Box);
    Box->SetOrigin(bounds[0], bounds[2], bounds[4]);
    Box->SetSpacing(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]);
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(Pos);
    plane->SetNormal(Nor);

    Cutter->SetCutFunction(plane);
    Cutter->SetInputData(Box);
    Cutter->Update();
    vtkSmartPointer<vtkPolyData>  polys = Cutter->GetOutput();
    vtkSmartPointer<vtkPoints> pts = polys->GetPoints();
    double p1[3], p2[3], origin[3];
    if (pts->GetNumberOfPoints()>0) pts->GetPoint(0,origin);
    else Box->GetOrigin(origin);
    if (pts->GetNumberOfPoints()>1) pts->GetPoint(1,p1);
//    else memcpy(p1, this->Point1, 3*sizeof(double));
    if (pts->GetNumberOfPoints()>2) pts->GetPoint(2,p2);
//    else memcpy(p2, this->Point2, 3*sizeof(double));

    vtkSMPropertyHelper pp0(this->ReferenceProxy->getProxy(), "Origin");
    pp0.Set(origin,3);
    //
    vtkSMPropertyHelper pp1(this->ReferenceProxy->getProxy(), "Point1");
    pp1.Set(p1,3);
    //
    vtkSMPropertyHelper pp2(this->ReferenceProxy->getProxy(), "Point2");
    pp2.Set(p2,3);

    vtkSMPropertyHelper pp3(this->ReferenceProxy->getProxy(), "Point3");
    pp3.Set(origin,3);
  }
  this->setModified();
//  this->ReferenceProxy->getProxy()->UpdateVTKObjects();
}
//-----------------------------------------------------------------------------
void pqSamplingGridPanel::UpdatePlanePoints()
{
  this->ComputePlanePoints();
//  this->render();
}
//-----------------------------------------------------------------------------
void pqSamplingGridPanel::ComputePlanePoints()
{
/*
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
*/
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void pqSamplingGridPanel::onUseCameraNormal()
{
  //return; // disabled until fixed
  //vtkSMNewWidgetRepresentationProxy* widget = this->getWidgetProxy();
  //if (widget)
  //  {
  //  pqRenderView* renView = qobject_cast<pqRenderView*>(this->renderView());
  //  if (vtkCamera* const camera = renView? 
  //    renView->getRenderViewProxy()->GetActiveCamera() : 0)
  //    {
  //    double camera_normal[3];
  //    camera->GetViewPlaneNormal(camera_normal);
  //    camera_normal[0] = -camera_normal[0];
  //    camera_normal[1] = -camera_normal[1];
  //    camera_normal[2] = -camera_normal[2];
  //    vtkSMPropertyHelper(widget, "Normal").Set(camera_normal, 3);
  //    widget->UpdateVTKObjects();
  //    this->render();
  //    this->setModified();
  //    }
  //  }
}
//-----------------------------------------------------------------------------
void pqSamplingGridPanel::onWidgetModeChanged(int index)
{
  this->Implementation->stackedWidget->setCurrentIndex(index);
  if (index==0) {
    this->boxWidget->hideWidget();
    this->planeWidget->showWidget();
    this->planeWidget->select();
  } 
  else {
    this->planeWidget->hideWidget();
    this->boxWidget->showWidget();
    this->boxWidget->select();
  }
  this->currentMode = index;
  this->onWidgetModified();
}
//-----------------------------------------------------------------------------
int pqSamplingGridPanel::getReferenceInputBounds(double bounds[6]) const
{
  vtkSMSourceProxy* input = NULL;
  vtkSMInputProperty* ivp = vtkSMInputProperty::SafeDownCast(
    this->ReferenceProxy->getProxy()->GetProperty("Input"));
  int output_port = 0;
  if (ivp && ivp->GetNumberOfProxies())
    {
    vtkSMProxy* pxy = ivp->GetProxy(0);
    input = vtkSMSourceProxy::SafeDownCast(pxy);
    output_port =ivp->GetOutputPortForConnection(0);
    }
  else
    {
    // reference proxy has no input. This generally happens when the widget is
    // controlling properties of a source. In that case, if the source has been
    // "created", simply use the source's bounds.
    input = vtkSMSourceProxy::SafeDownCast(this->ReferenceProxy->getProxy());
    }

  if(input)
    {
    input->GetDataInformation(output_port)->GetBounds(bounds);
    return (bounds[1] >= bounds[0] && bounds[3] >= bounds[2] && bounds[5] >=
      bounds[4]) ? 1 : 0;
    }
  return 0;
}
