/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCustomBoxRepresentation.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkCustomBoxRepresentation.h"
#include "vtkActor.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkCallbackCommand.h"
#include "vtkBox.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkInteractorObserver.h"
#include "vtkMath.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkTransform.h"
#include "vtkDoubleArray.h"
#include "vtkBox.h"
#include "vtkPlanes.h"
#include "vtkCamera.h"
#include "vtkAssemblyPath.h"
#include "vtkWindow.h"
#include "vtkObjectFactory.h"
#include "vtkConeSource.h"
#include "vtkLineSource.h"
#include "vtkBoundingBox.h"

vtkCxxRevisionMacro(vtkCustomBoxRepresentation, "$Revision: 1.9 $");
vtkStandardNewMacro(vtkCustomBoxRepresentation);

//----------------------------------------------------------------------------
vtkCustomBoxRepresentation::vtkCustomBoxRepresentation()
{
  // in C++ constructors don't call virtual methods, so call our function
  // to add to the ones created in the vtkBoxRepresentation constructor
  vtkCustomBoxRepresentation::CreateDefaultProperties();

  //default to off
  this->PlaneMode = 0;

  // Create the plane normal
  this->ConeSource = vtkConeSource::New();
  this->ConeSource->SetResolution(12);
  this->ConeSource->SetAngle(25.0);
  this->ConeMapper = vtkPolyDataMapper::New();
  this->ConeMapper->SetInput(this->ConeSource->GetOutput());
  this->ConeActor = vtkActor::New();
  this->ConeActor->SetMapper(this->ConeMapper);
  this->ConeActor->SetProperty(this->NormalProperty);
  //
  this->GenerateConnectedCells = 1;
  this->SetPlaneMode(0);
}

//----------------------------------------------------------------------------
vtkCustomBoxRepresentation::~vtkCustomBoxRepresentation()
{  
  this->ConeSource->Delete();
  this->ConeMapper->Delete();
  this->ConeActor->Delete();

  this->NormalProperty->Delete();
  this->SelectedNormalProperty->Delete();
}
//----------------------------------------------------------------------
void vtkCustomBoxRepresentation::SetPlaneMode(int planemode)
{
  if (this->PlaneMode == planemode)
  {
    return;
  }
  this->PlaneMode = planemode;
  if (this->PlaneMode) {
    this->HandlePicker->DeletePickList(this->Handle[4]);
    this->HandlePicker->DeletePickList(this->Handle[5]);
    this->Handle[4]->SetVisibility(0);
    this->Handle[5]->SetVisibility(0);
    this->ConeActor->SetVisibility(1);
    this->FlattenToPlane();
    this->PositionHandles();    
  }
  else {
    this->HandlePicker->AddPickList(this->Handle[4]);
    this->HandlePicker->AddPickList(this->Handle[5]);
    this->Handle[4]->SetVisibility(1);
    this->Handle[5]->SetVisibility(1);
    this->ConeActor->SetVisibility(0);
  }
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::ComputeNormals()
{
  vtkBoxRepresentation::ComputeNormals();
  this->ConeSource->SetDirection(this->N[5]);
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::CreateDefaultProperties()
{
  // Normal properties
  this->NormalProperty = vtkProperty::New();
  this->NormalProperty->SetColor(1,1,1);
  this->NormalProperty->SetLineWidth(2);

  this->SelectedNormalProperty = vtkProperty::New();
  this->SelectedNormalProperty->SetColor(1,0,0);
  this->NormalProperty->SetLineWidth(2);
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::PlaceWidget(double bds[6])
{
  int i;
  double bounds[6], center[3];
  
  this->AdjustBounds(bds,bounds,center);
  
  this->Points->SetPoint(0, bounds[0], bounds[2], bounds[4]);
  this->Points->SetPoint(1, bounds[1], bounds[2], bounds[4]);
  this->Points->SetPoint(2, bounds[1], bounds[3], bounds[4]);
  this->Points->SetPoint(3, bounds[0], bounds[3], bounds[4]);
  this->Points->SetPoint(4, bounds[0], bounds[2], bounds[5]);
  this->Points->SetPoint(5, bounds[1], bounds[2], bounds[5]);
  this->Points->SetPoint(6, bounds[1], bounds[3], bounds[5]);
  this->Points->SetPoint(7, bounds[0], bounds[3], bounds[5]);

  if (this->PlaneMode) {
    this->FlattenToPlane();
    bounds[4] = bounds[5] = (bounds[4] + bounds[5])/2.0;
    double epsilion = ((bounds[3]-bounds[2])+(bounds[1]-bounds[0]))/1000.0;
    bounds[4] -= epsilion;
    bounds[5] += epsilion;
  }

  for (i=0; i<6; i++)
    {
    this->InitialBounds[i] = bounds[i];
    }
  this->InitialLength = sqrt((bounds[1]-bounds[0])*(bounds[1]-bounds[0]) +
                             (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) +
                             (bounds[5]-bounds[4])*(bounds[5]-bounds[4]));

  this->PositionHandles();
  this->ComputeNormals();
  this->ValidPick = 1; //since we have set up widget
  this->SizeHandles();
}

//----------------------------------------------------------------------------
#define BOXWIDGET_XMY(x,y,z) \
  z[0] = x[0] - y[0]; \
  z[1] = x[1] - y[1]; \
  z[2] = x[2] - y[2]; 
#define BOXWIDGET_AXPY(a,x,y,z) \
  z[0] = a*x[0] + y[0]; \
  z[1] = a*x[1] + y[1]; \
  z[2] = a*x[2] + y[2]; 
#define BOXWIDGET_PAX(a,x,z) \
  z[0] += a*x[0]; \
  z[1] += a*x[1]; \
  z[2] += a*x[2]; 
#define BOXWIDGET_AXPBY(a,x,b,y,z) \
  z[0] = a*x[0] + b*y[0]; \
  z[1] = a*x[1] + b*y[1]; \
  z[2] = a*x[2] + b*y[2]; 
#define BOXWIDGET_AXPBYPCZ(a,x,b,y,c,z,r) \
  r[0] = a*x[0] + b*y[0] + c*z[0]; \
  r[1] = a*x[1] + b*y[1] + c*z[1]; \
  r[2] = a*x[2] + b*y[2] + c*z[2]; 
#define VTK_AVERAGE(a,b,c) \
  c[0] = (a[0] + b[0])/2.0; \
  c[1] = (a[1] + b[1])/2.0; \
  c[2] = (a[2] + b[2])/2.0;

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::FlattenToPlane()
{
  vtkBoundingBox box(this->GetBounds());
  double epsilon = box.GetDiagonalLength()/1000.0;
  double *pts =
    static_cast<vtkDoubleArray *>(this->Points->GetData())->GetPointer(0);
  double *hzm = pts + 3*12;
  double *hzp = pts + 3*13;
  double *x0 = pts + 3*0;
  double *x1 = pts + 3*1;
  double *x2 = pts + 3*2;
  double *x3 = pts + 3*3;
  double *x4 = pts + 3*4;
  double *x5 = pts + 3*5;
  double *x6 = pts + 3*6;
  double *x7 = pts + 3*7;
  //
  double vp0p4[3];
  BOXWIDGET_XMY(x4,x0,vp0p4);
  vtkMath::Normalize(vp0p4);
  for (int i=0; i<3; i++) {
    x0[i] = x4[i] = (x0[i]+x4[i])/2.0;
    x1[i] = x5[i] = (x1[i]+x5[i])/2.0;
    x2[i] = x6[i] = (x2[i]+x6[i])/2.0;
    x3[i] = x7[i] = (x3[i]+x7[i])/2.0;
    x0[i] -= vp0p4[i]*epsilon;
    x4[i] += vp0p4[i]*epsilon;
    x1[i] -= vp0p4[i]*epsilon;
    x5[i] += vp0p4[i]*epsilon;
    x2[i] -= vp0p4[i]*epsilon;
    x6[i] += vp0p4[i]*epsilon;
    x3[i] -= vp0p4[i]*epsilon;
    x7[i] += vp0p4[i]*epsilon;
    hzm[i] = hzp[i] = (hzm[i]+hzp[i])/2.0;
  }
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::ReOrient(double norm[3])
{
  double vec1[3], vec2[3];
  vtkMath::Perpendiculars(norm, vec1, vec2, 0.0);
  //
  double sx,sx2,sy,sy2,sz,sz2;
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x0 = pts + 3*0;
  double *x1 = pts + 3*1;
  double *x2 = pts + 3*2;
  double *x3 = pts + 3*3;
  double *x4 = pts + 3*4;
  double *x5 = pts + 3*5;
  double *x6 = pts + 3*6;
  double *x7 = pts + 3*7;
  double *bc = pts + 3*14;
  //
  sx2 = sqrt(vtkMath::Distance2BetweenPoints(x0, x1));
  sx  = sx2/2.0;
  sy2 = sqrt(vtkMath::Distance2BetweenPoints(x0, x3));
  sy  = sy2/2.0;
  sz2 = sqrt(vtkMath::Distance2BetweenPoints(x0, x4));
  sz  = sz2/2.0;
  // compute corner p0
  BOXWIDGET_AXPBY(-sx, vec1, -sy, vec2, x0);
  BOXWIDGET_PAX  (-sz, norm, x0);
  BOXWIDGET_PAX  (1.0, bc,   x0);
  // compute corner p1
  BOXWIDGET_AXPBY( sx, vec1, -sy, vec2, x1);
  BOXWIDGET_PAX  (-sz, norm, x1);
  BOXWIDGET_PAX  (1.0, bc,   x1);
  // make {p4, p5} by adding z to {p0, p1}
  BOXWIDGET_AXPY(sz2, norm, x0, x4);
  BOXWIDGET_AXPY(sz2, norm, x1, x5);
  // make remaining 4 corners by offsetting from first 4
  BOXWIDGET_AXPY(sy2, vec2, x0, x3);
  BOXWIDGET_AXPY(sy2, vec2, x1, x2);
  BOXWIDGET_AXPY(sy2, vec2, x5, x6);
  BOXWIDGET_AXPY(sy2, vec2, x4, x7);
  //
  this->PositionHandles();
  if (this->PlaneMode) {
    this->FlattenToPlane();
  }
  this->ComputeNormals();
  this->SizeHandles();
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::OrientX(void)
{
  double norm[3] = { 1, 0, 0 };
  this->ReOrient(norm);
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::OrientY(void)
{
  double norm[3] = { 0, 1, 0 };
  this->ReOrient(norm);
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::OrientZ(void)
{
  double norm[3] = { 0, 0, 1 };
  this->ReOrient(norm);
}
//----------------------------------------------------------------------------
double *vtkCustomBoxRepresentation::GetNormal()
{
  this->ComputeNormals();
  return this->N[0];
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::GetNormal(double xyz[3])
{
  this->ComputeNormals();
  xyz[0] = this->N[0][0];
  xyz[1] = this->N[0][1];
  xyz[2] = this->N[0][2];
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::SetOBB(double o[12])
{
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x0 = pts + 3*0;
  double *x1 = pts + 3*1;
  double *x2 = pts + 3*2;
  double *x3 = pts + 3*3;
  double *x4 = pts + 3*4;
  double *x5 = pts + 3*5;
  double *x6 = pts + 3*6;
  double *x7 = pts + 3*7;
  double *bc = pts + 3*14;
  x0[0] = o[0]; x0[1] = o[1];  x0[2] = o[2];
  x1[0] = o[3]; x1[1] = o[4];  x1[2] = o[5];
  x3[0] = o[6]; x3[1] = o[7];  x3[2] = o[8];
  x4[0] = o[9]; x4[1] = o[10]; x4[2] = o[11];
  // axes vectors
  double v0v1[3], v0v3[3], v0v4[3];
  BOXWIDGET_XMY(x1,x0,v0v1);
  BOXWIDGET_XMY(x3,x0,v0v3);
  BOXWIDGET_XMY(x4,x0,v0v4);
  // make {p2} by adding y to {p1}
  BOXWIDGET_AXPY(1.0, v0v3, x1, x2);
  // make {p5} by adding z to {p1}
  BOXWIDGET_AXPY(1.0, v0v4, x1, x5);
  // make {p6} by adding z to {p2}
  BOXWIDGET_AXPY(1.0, v0v4, x2, x6);
  // make {p7} by adding z to {p3}
  BOXWIDGET_AXPY(1.0, v0v4, x3, x7);

  // compute new centre by midpoint of p0,p6
  VTK_AVERAGE(x0,x6,bc);
  //
  this->PositionHandles();
  if (this->PlaneMode) {
    this->FlattenToPlane();
  }
  this->ComputeNormals();
  this->SizeHandles();
}
//----------------------------------------------------------------------------
double *vtkCustomBoxRepresentation::GetOrigin()
{
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x0 = pts + 3*0;
  return x0;
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::GetOrigin(double xyz[3])
{
  double *p = this->GetOrigin();
  xyz[0] = p[0]; xyz[1] = p[1]; xyz[2] = p[2];
}
//----------------------------------------------------------------------------
double* vtkCustomBoxRepresentation::GetPoint1()
{
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x1 = pts + 3*1;
  return x1;
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::GetPoint1(double xyz[3])
{
  double *p = this->GetPoint1();
  xyz[0] = p[0]; xyz[1] = p[1]; xyz[2] = p[2];
}
//----------------------------------------------------------------------------
double* vtkCustomBoxRepresentation::GetPoint2()
{
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x3 = pts + 3*3;
  return x3;
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::GetPoint2(double xyz[3])
{
  double *p = this->GetPoint2();
  xyz[0] = p[0]; xyz[1] = p[1]; xyz[2] = p[2];
}
//----------------------------------------------------------------------------
double* vtkCustomBoxRepresentation::GetPoint3()
{
  double *pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double *x4 = pts + 3*4;
  return x4;
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::GetPoint3(double xyz[3])
{
  double *p = this->GetPoint3();
  xyz[0] = p[0]; xyz[1] = p[1]; xyz[2] = p[2];
}
//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  vtkBoxRepresentation::ReleaseGraphicsResources(w);
  this->ConeActor->ReleaseGraphicsResources(w);
}

//----------------------------------------------------------------------------
int vtkCustomBoxRepresentation::RenderOpaqueGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();
  
  count += this->ConeActor->RenderOpaqueGeometry(v);
  count += this->HexActor->RenderOpaqueGeometry(v);
  count += this->HexOutline->RenderOpaqueGeometry(v);
  count += this->HexFace->RenderOpaqueGeometry(v);
  // render the handles
  for (int i=0; i<7; i++)
    {
    if (this->PlaneMode && (i==4 || i==5)) 
      {
        continue;
      }
    count += this->Handle[i]->RenderOpaqueGeometry(v);
    }

  return count;
}

//----------------------------------------------------------------------------
int vtkCustomBoxRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();
  
  count += this->ConeActor->RenderTranslucentPolygonalGeometry(v);
  count += this->HexActor->RenderTranslucentPolygonalGeometry(v);
  count += this->HexOutline->RenderTranslucentPolygonalGeometry(v);
  count += this->HexFace->RenderTranslucentPolygonalGeometry(v);
  // render the handles
  for (int i=0; i<7; i++)
    {
    if (this->PlaneMode && (i==4 || i==5)) 
      {
        continue;
      }
    count += this->Handle[i]->RenderTranslucentPolygonalGeometry(v);
    }

  return count;
}

//----------------------------------------------------------------------------
int vtkCustomBoxRepresentation::HasTranslucentPolygonalGeometry()
{
  int result=0;
  this->BuildRepresentation();

  result |= this->ConeActor->HasTranslucentPolygonalGeometry();
  result |= this->HexActor->HasTranslucentPolygonalGeometry();
  result |= this->HexOutline->HasTranslucentPolygonalGeometry();

  // If the face is not selected, we are not really rendering translucent faces,
  // hence don't bother taking it's opacity into consideration.
  // Look at BUG #7301.
  if (this->HexFace->GetProperty() == this->SelectedFaceProperty)
    {
    result |= this->HexFace->HasTranslucentPolygonalGeometry();
    }

  // render the handles
  for (int i=0; i<7; i++)
    {
    if (this->PlaneMode && (i==4 || i==5)) 
      {
        continue;
      }
    result |= this->Handle[i]->HasTranslucentPolygonalGeometry();
    }

  return result;
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::PositionHandles()
{
  vtkBoxRepresentation::PositionHandles();
  this->ConeSource->SetCenter(this->Points->GetPoint(13));
}
#undef VTK_AVERAGE

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::HandlesOn()
{
  for (int i=0; i<7; i++)
    {
    if (this->PlaneMode && (i==4 || i==5)) 
      {
        continue;
      }
    this->Handle[i]->VisibilityOn();
    }
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::SizeHandles()
{
  vtkBoxRepresentation::SizeHandles();
  double *center 
    = static_cast<vtkDoubleArray *>(this->Points->GetData())->GetPointer(3*14);
  double radius =
      this->vtkWidgetRepresentation::SizeHandlesInPixels(1.5,center);
  this->ConeSource->SetHeight(2.0*radius);
  this->ConeSource->SetRadius(2.0*radius);
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
