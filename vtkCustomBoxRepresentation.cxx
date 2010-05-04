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
  // The initial state
  this->InteractionState = vtkCustomBoxRepresentation::Outside;

  // Handle size is in pixels for this widget
  this->HandleSize = 5.0;

  // Control orientation of normals
  this->InsideOut = 0;
  this->OutlineFaceWires = 0;
  this->OutlineCursorWires = 1;

  //default to off
  this->PlaneMode = 0;

  // Set up the initial properties
  this->CreateDefaultProperties();

  // Create the plane normal
  this->ConeSource = vtkConeSource::New();
  this->ConeSource->SetResolution(12);
  this->ConeSource->SetAngle(25.0);
  this->ConeMapper = vtkPolyDataMapper::New();
  this->ConeMapper->SetInput(this->ConeSource->GetOutput());
  this->ConeActor = vtkActor::New();
  this->ConeActor->SetMapper(this->ConeMapper);
  this->ConeActor->SetProperty(this->NormalProperty);

  // Construct the poly data representing the hex
  this->HexPolyData = vtkPolyData::New();
  this->HexMapper = vtkPolyDataMapper::New();
  this->HexMapper->SetInput(HexPolyData);
  this->HexActor = vtkActor::New();
  this->HexActor->SetMapper(this->HexMapper);
  this->HexActor->SetProperty(this->OutlineProperty);

  // Construct initial points
  this->Points = vtkPoints::New(VTK_DOUBLE);
  this->Points->SetNumberOfPoints(15);//8 corners; 6 faces; 1 center
  this->HexPolyData->SetPoints(this->Points);
  
  // Construct connectivity for the faces. These are used to perform
  // the picking.
  int i;
  vtkIdType pts[4];
  vtkCellArray *cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(6,4));
  pts[0] = 3; pts[1] = 0; pts[2] = 4; pts[3] = 7;
  cells->InsertNextCell(4,pts);
  pts[0] = 1; pts[1] = 2; pts[2] = 6; pts[3] = 5;
  cells->InsertNextCell(4,pts);
  pts[0] = 0; pts[1] = 1; pts[2] = 5; pts[3] = 4;
  cells->InsertNextCell(4,pts);
  pts[0] = 2; pts[1] = 3; pts[2] = 7; pts[3] = 6;
  cells->InsertNextCell(4,pts);
  pts[0] = 0; pts[1] = 3; pts[2] = 2; pts[3] = 1;
  cells->InsertNextCell(4,pts);
  pts[0] = 4; pts[1] = 5; pts[2] = 6; pts[3] = 7;
  cells->InsertNextCell(4,pts);
  this->HexPolyData->SetPolys(cells);
  cells->Delete();
  this->HexPolyData->BuildCells();
  
  // The face of the hexahedra
  cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(1,4));
  cells->InsertNextCell(4,pts); //temporary, replaced later
  this->HexFacePolyData = vtkPolyData::New();
  this->HexFacePolyData->SetPoints(this->Points);
  this->HexFacePolyData->SetPolys(cells);
  this->HexFaceMapper = vtkPolyDataMapper::New();
  this->HexFaceMapper->SetInput(HexFacePolyData);
  this->HexFace = vtkActor::New();
  this->HexFace->SetMapper(this->HexFaceMapper);
  this->HexFace->SetProperty(this->FaceProperty);
  cells->Delete();

  // Create the outline for the hex
  this->OutlinePolyData = vtkPolyData::New();
  this->OutlinePolyData->SetPoints(this->Points);
  this->OutlineMapper = vtkPolyDataMapper::New();
  this->OutlineMapper->SetInput(this->OutlinePolyData);
  this->HexOutline = vtkActor::New();
  this->HexOutline->SetMapper(this->OutlineMapper);
  this->HexOutline->SetProperty(this->OutlineProperty);
  cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(15,2));
  this->OutlinePolyData->SetLines(cells);
  cells->Delete();

  // Create the outline
  this->GenerateOutline();

  // Create the handles
  this->Handle = new vtkActor* [7];
  this->HandleMapper = new vtkPolyDataMapper* [7];
  this->HandleGeometry = new vtkSphereSource* [7];
  for (i=0; i<7; i++)
    {
    this->HandleGeometry[i] = vtkSphereSource::New();
    this->HandleGeometry[i]->SetThetaResolution(16);
    this->HandleGeometry[i]->SetPhiResolution(8);
    this->HandleMapper[i] = vtkPolyDataMapper::New();
    this->HandleMapper[i]->SetInput(this->HandleGeometry[i]->GetOutput());
    this->Handle[i] = vtkActor::New();
    this->Handle[i]->SetMapper(this->HandleMapper[i]);
    }
  
  // Define the point coordinates
  double bounds[6];
  bounds[0] = -0.5;
  bounds[1] = 0.5;
  bounds[2] = -0.5;
  bounds[3] = 0.5;
  bounds[4] = -0.5;
  bounds[5] = 0.5;
  // Points 8-14 are down by PositionHandles();
  this->BoundingBox = vtkBox::New();
  this->PlaceWidget(bounds);

  //Manage the picking stuff
  this->HandlePicker = vtkCellPicker::New();
  this->HandlePicker->SetTolerance(0.001);
  for (i=0; i<7; i++)
    {
    this->HandlePicker->AddPickList(this->Handle[i]);
    }
  this->HandlePicker->PickFromListOn();

  this->HexPicker = vtkCellPicker::New();
  this->HexPicker->SetTolerance(0.001);
  this->HexPicker->AddPickList(HexActor);
  this->HexPicker->PickFromListOn();
  
  this->CurrentHandle = NULL;

  // Internal data memebers for performance
  this->Transform = vtkTransform::New();
  this->PlanePoints = vtkPoints::New(VTK_DOUBLE);
  this->PlanePoints->SetNumberOfPoints(6);
  this->PlaneNormals = vtkDoubleArray::New();
  this->PlaneNormals->SetNumberOfComponents(3);
  this->PlaneNormals->SetNumberOfTuples(6);
  this->Matrix = vtkMatrix4x4::New();

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
  double *pts =
     static_cast<vtkDoubleArray *>(this->Points->GetData())->GetPointer(0);
  double *p0 = pts;
  double *px = pts + 3*1;
  double *py = pts + 3*3;
  double *pz = pts + 3*4;
  int i;
  
  for (i=0; i<3; i++)
    {
    this->N[0][i] = p0[i] - px[i];
    this->N[2][i] = p0[i] - py[i];
    this->N[4][i] = p0[i] - pz[i];
    }
  vtkMath::Normalize(this->N[0]);
  vtkMath::Normalize(this->N[2]);
  vtkMath::Normalize(this->N[4]);
  for (i=0; i<3; i++)
    {
    this->N[1][i] = -this->N[0][i];
    this->N[3][i] = -this->N[2][i];
    this->N[5][i] = -this->N[4][i];
    }
  this->ConeSource->SetDirection(this->N[5]);

}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::CreateDefaultProperties()
{
  // Handle properties
  this->HandleProperty = vtkProperty::New();
  this->HandleProperty->SetColor(1,1,1);

  this->SelectedHandleProperty = vtkProperty::New();
  this->SelectedHandleProperty->SetColor(1,0,0);

  // Face properties
  this->FaceProperty = vtkProperty::New();
  this->FaceProperty->SetColor(1,1,1);
  this->FaceProperty->SetOpacity(0.0);

  this->SelectedFaceProperty = vtkProperty::New();
  this->SelectedFaceProperty->SetColor(1,1,0);
  this->SelectedFaceProperty->SetOpacity(0.25);
  
  // Outline properties
  this->OutlineProperty = vtkProperty::New();
  this->OutlineProperty->SetRepresentationToWireframe();
  this->OutlineProperty->SetAmbient(1.0);
  this->OutlineProperty->SetAmbientColor(1.0,1.0,1.0);
  this->OutlineProperty->SetLineWidth(2.0);

  this->SelectedOutlineProperty = vtkProperty::New();
  this->SelectedOutlineProperty->SetRepresentationToWireframe();
  this->SelectedOutlineProperty->SetAmbient(1.0);
  this->SelectedOutlineProperty->SetAmbientColor(0.0,1.0,0.0);
  this->SelectedOutlineProperty->SetLineWidth(2.0);

  // Normal properties
  this->NormalProperty = vtkProperty::New();
  this->NormalProperty->SetColor(1,1,1);
  this->NormalProperty->SetLineWidth(2);

  this->SelectedNormalProperty = vtkProperty::New();
  this->SelectedNormalProperty->SetColor(1,0,0);
  this->NormalProperty->SetLineWidth(2);

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
#define BOXWIDGET_AVERAGE(a,b,c) \
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
  BOXWIDGET_AVERAGE(x0,x6,bc);
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
void vtkCustomBoxRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  this->ConeActor->ReleaseGraphicsResources(w);
  this->HexActor->ReleaseGraphicsResources(w);
  this->HexOutline->ReleaseGraphicsResources(w);
  this->HexFace->ReleaseGraphicsResources(w);
  // render the handles
  for (int j=0; j<7; j++)
    {
    this->Handle[j]->ReleaseGraphicsResources(w);
    }

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
  double *pts =
     static_cast<vtkDoubleArray *>(this->Points->GetData())->GetPointer(0);
  double *p0 = pts;
  double *p1 = pts + 3*1;
  double *p2 = pts + 3*2;
  double *p3 = pts + 3*3;
  //double *p4 = pts + 3*4;
  double *p5 = pts + 3*5;
  double *p6 = pts + 3*6;
  double *p7 = pts + 3*7;
  double x[3];

  BOXWIDGET_AVERAGE(p0,p7,x);
  this->Points->SetPoint(8, x);
  BOXWIDGET_AVERAGE(p1,p6,x);
  this->Points->SetPoint(9, x);
  BOXWIDGET_AVERAGE(p0,p5,x);
  this->Points->SetPoint(10, x);
  BOXWIDGET_AVERAGE(p2,p7,x);
  this->Points->SetPoint(11, x);
  BOXWIDGET_AVERAGE(p1,p3,x);
  this->Points->SetPoint(12, x);
  BOXWIDGET_AVERAGE(p5,p7,x);
  this->Points->SetPoint(13, x);
  BOXWIDGET_AVERAGE(p0,p6,x);
  this->Points->SetPoint(14, x);

  int i;
  for (i = 0; i < 7; ++i)
    {
    this->HandleGeometry[i]->SetCenter(this->Points->GetPoint(8+i));
    }

  this->ConeSource->SetCenter(this->Points->GetPoint(13));

  this->Points->GetData()->Modified();
  this->HexFacePolyData->Modified();
  this->HexPolyData->Modified();
  this->GenerateOutline();
}
#undef BOXWIDGET_AVERAGE

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
  double *center 
    = static_cast<vtkDoubleArray *>(this->Points->GetData())->GetPointer(3*14);
  double radius =
      this->vtkWidgetRepresentation::SizeHandlesInPixels(1.5,center);
  for(int i=0; i<7; i++)
    {
    this->HandleGeometry[i]->SetRadius(radius);
    }
  this->ConeSource->SetHeight(2.0*radius);
  this->ConeSource->SetRadius(2.0*radius);
}

//----------------------------------------------------------------------------
void vtkCustomBoxRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  double *bounds=this->InitialBounds;
  os << indent << "Initial Bounds: "
     << "(" << bounds[0] << "," << bounds[1] << ") "
     << "(" << bounds[2] << "," << bounds[3] << ") " 
     << "(" << bounds[4] << "," << bounds[5] << ")\n";

  if ( this->HandleProperty )
    {
    os << indent << "Handle Property: " << this->HandleProperty << "\n";
    }
  else
    {
    os << indent << "Handle Property: (none)\n";
    }
  if ( this->SelectedHandleProperty )
    {
    os << indent << "Selected Handle Property: " 
       << this->SelectedHandleProperty << "\n";
    }
  else
    {
    os << indent << "SelectedHandle Property: (none)\n";
    }

  if ( this->FaceProperty )
    {
    os << indent << "Face Property: " << this->FaceProperty << "\n";
    }
  else
    {
    os << indent << "Face Property: (none)\n";
    }
  if ( this->SelectedFaceProperty )
    {
    os << indent << "Selected Face Property: " 
       << this->SelectedFaceProperty << "\n";
    }
  else
    {
    os << indent << "Selected Face Property: (none)\n";
    }

  if ( this->OutlineProperty )
    {
    os << indent << "Outline Property: " << this->OutlineProperty << "\n";
    }
  else
    {
    os << indent << "Outline Property: (none)\n";
    }
  if ( this->SelectedOutlineProperty )
    {
    os << indent << "Selected Outline Property: " 
       << this->SelectedOutlineProperty << "\n";
    }
  else
    {
    os << indent << "Selected Outline Property: (none)\n";
    }

  os << indent << "Outline Face Wires: "
     << (this->OutlineFaceWires ? "On\n" : "Off\n");
  os << indent << "Outline Cursor Wires: "
     << (this->OutlineCursorWires ? "On\n" : "Off\n");
  os << indent << "Inside Out: " << (this->InsideOut ? "On\n" : "Off\n");
}
