/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParticleIdFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkParticleIdFilter - generate scalars or field data from point and cell ids
// .SECTION Description
// vtkParticleIdFilter is a filter to that generates scalars or field data
// using cell and point ids. That is, the point attribute data scalars
// or field data are generated from the point ids, and the cell
// attribute data scalars or field data are generated from the the
// cell ids.
//
// Typically this filter is used with vtkLabeledDataMapper (and possibly
// vtkSelectVisiblePoints) to create labels for points and cells, or labels
// for the point or cell data scalar values.

#ifndef __vtkParticleIdFilter_h
#define __vtkParticleIdFilter_h

#include "vtkIdFilter.h"
class vtkMultiProcessController;

class VTK_EXPORT vtkParticleIdFilter : public vtkIdFilter
{
public:
  vtkTypeMacro(vtkParticleIdFilter,vtkIdFilter);

  // Description:
  // Construct object with PointIds and CellIds on; and ids being generated
  // as scalars.
  static vtkParticleIdFilter *New();

  // Description:
  // By default this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

protected:
  vtkParticleIdFilter();
  ~vtkParticleIdFilter();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  //
  vtkMultiProcessController *Controller;
  int UpdatePiece;
  int UpdateNumPieces;

private:
  vtkParticleIdFilter(const vtkParticleIdFilter&);  // Not implemented.
  void operator=(const vtkParticleIdFilter&);  // Not implemented.
};

#endif


