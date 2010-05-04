/*=========================================================================

   Program:   ParaQ
   Module:    $RCSfile: pqRegularGridSourceWidget.h,v $

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

#ifndef _pqRegularGridSourceWidget_h
#define _pqRegularGridSourceWidget_h

#include "pq3DWidget.h"

class pqServer;

/// Provides a complete Qt UI for working with a 3D box/plane widget

/// BUG: Currently we work directly with the source, when compuing OBB
/// or snap to X/Y/Z axis, this means properties are set directly on
/// the objects and we can't "undo" by pressing reset before accept.

class pqRegularGridSourceWidget :
  public pq3DWidget
{
  typedef pq3DWidget Superclass;
  
  Q_OBJECT
  
public:
   pqRegularGridSourceWidget(vtkSMProxy* o, vtkSMProxy* pxy, QWidget* p = 0);
  ~pqRegularGridSourceWidget();

  virtual void createWidget(pqServer* server);

  /// Resets the bounds of the 3D widget to the reference proxy bounds.
  virtual void resetBounds(double bounds[6]);
  virtual void resetBounds();

  /// We need to generate the correct source settings
  /// (origin, point1, point2, point3) when the widget is accepted.
  virtual void accept();

  /// Do the correct (origin, point1, point2, point3) calculation
  void ComputePlanePoints();

  /// these are used internally and stored so that the 
  /// pqRegularGridSourceWidget subclass can use them to
  /// compute plane p1 and p2 points
  double InputBounds[6];
  double SizeFactor;


private slots:
  /// Called to set the widget normal to the X axis
  void onUseXNormal();
  /// Called to set the widget normal to the Y axis
  void onUseYNormal();
  /// Reset box to OBB of data
  void onComputeOBB();
  /// Called to set the widget normal to the camera direction
  void onUseZNormal();
  /// Called to set the widget normal to the camera direction
  void onUseCameraNormal();
  /// called when switching between plane/box mode
  void onWidgetModeChanged(int index);
  /// Called whenever the 3D widget visibility is modified
  void onWidgetVisibilityChanged(bool visible);
  /// When a plane parater changes, we must recompute
  void UpdatePlanePoints();

private:
  class pqImplementation;
  pqImplementation *Implementation;
};

#endif
