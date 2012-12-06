/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCustomBoxRepresentation.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCustomBoxRepresentation - a class defining the representation for the vtkCustomBoxWidget
// .SECTION Description
// This class is a concrete representation for the vtkCustomBoxWidget. It
// represents a box with seven handles: one on each of the six faces, plus a
// center handle. Through interaction with the widget, the box
// representation can be arbitrarily positioned in the 3D space.
//
// To use this representation, you normally use the PlaceWidget() method
// to position the widget at a specified region in space.
//
// .SECTION Caveats
// This class, and vtkvtkCustomBoxRepresentationBoxWidget, are second generation VTK
// widgets. An earlier version of this functionality was defined in the
// class vtkBoxWidget.

// .SECTION See Also
// vtkvtkCustomBoxRepresentationBoxWidget vtkBoxWidget


#ifndef __vtkCustomBoxRepresentation_h
#define __vtkCustomBoxRepresentation_h

#include "vtkBoxRepresentation.h"

class vtkConeSource;

class VTK_EXPORT vtkCustomBoxRepresentation : public vtkBoxRepresentation
{
public:
  // Description:
  // Instantiate the class.
  static vtkCustomBoxRepresentation *New();

  // Description:
  // Standard methods for the class.
  vtkTypeMacro(vtkCustomBoxRepresentation,vtkBoxRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The box widget/representation can also be used to manipulate a plane
  // which is bounded by the box, but flat (2D). When plane mode is set
  // one side of the box (the z side) is collapsed to an xy plane.
  void SetPlaneMode(int planemode);
  vtkGetMacro(PlaneMode,int);
  vtkBooleanMacro(PlaneMode,int);

  // Description:
  // Switches handles (the spheres) on or off by manipulating the underlying
  // actor visibility.
  virtual void HandlesOn();
  
  // Description:
  // These are methods that satisfy vtkWidgetRepresentation's API.
  virtual void PlaceWidget(double bounds[6]);
  
  // Description
  // These functions are used interactively to snap the widget to X/Y/Z planes
  // leaving the centre unchanged.
  virtual void OrientX(void);
  virtual void OrientY(void);
  virtual void OrientZ(void);

  // Description:
  // Get the normal to the plane.
  virtual void SetOBB(double o[12]);
  
  // Description:
  // Get the normal to the plane.
  virtual double* GetNormal();
  virtual void GetNormal(double xyz[3]);
  
  // Description:
  // Get the origin of the box
  virtual double* GetOrigin();
  virtual void GetOrigin(double xyz[3]);

  // Description:
  // Get the point representing the X axis line relative to origin
  virtual double* GetPoint1();
  virtual void GetPoint1(double xyz[3]);

  // Description:
  // Get the point representing the y axis line relative to origin
  virtual double* GetPoint2();
  virtual void GetPoint2(double xyz[3]);

  // Description:
  // Get the point representing the z axis line relative to origin
  virtual double* GetPoint3();
  virtual void GetPoint3(double xyz[3]);

  // Description
  // Used when controlling a regular grid source. This parameter is not used, but is
  // required o tht when the property is set, it can be linked to the RegularGridSource
  // controlled by the widget
  vtkSetVector3Macro(Resolution, int);
  vtkGetVector3Macro(Resolution, int);

  // Description
  // Used when controlling a regular grid source. This parameter is not used, but is
  // required o tht when the property is set, it can be linked to the RegularGridSource
  // controlled by the widget
  vtkSetMacro(GenerateConnectedCells,int);
  vtkGetMacro(GenerateConnectedCells,int);
  vtkBooleanMacro(GenerateConnectedCells,int);

  // Description:
  // Methods supporting, and required by, the rendering process.
  virtual void ReleaseGraphicsResources(vtkWindow*);
  virtual int  RenderOpaqueGeometry(vtkViewport*);
  virtual int  RenderTranslucentPolygonalGeometry(vtkViewport*);
  virtual int  HasTranslucentPolygonalGeometry();

protected:
   vtkCustomBoxRepresentation();
  ~vtkCustomBoxRepresentation();

  // The +direction normal cone
  vtkConeSource     *ConeSource;
  vtkPolyDataMapper *ConeMapper;
  vtkActor          *ConeActor;

  // The +direction normal line
  vtkLineSource     *LineSource;
  vtkPolyDataMapper *LineMapper;
  vtkActor          *LineActor;

  // The +direction normal property
  vtkProperty *NormalProperty;
  vtkProperty *SelectedNormalProperty;

  // Options to control RegularGridSource
  int PlaneMode;
  int Resolution[3];
  int GenerateConnectedCells;

  // Functions which are overriden by our widget
  virtual void CreateDefaultProperties();
  virtual void ComputeNormals();
  virtual void PositionHandles();
  virtual void SizeHandles();
  virtual void FlattenToPlane();
  virtual void ReOrient(double norm[3]);

private:
  vtkCustomBoxRepresentation(const vtkCustomBoxRepresentation&);  //Not implemented
  void operator=(const vtkCustomBoxRepresentation&);  //Not implemented
};

#endif
