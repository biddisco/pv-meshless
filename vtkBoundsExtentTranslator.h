/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBoundsExtentTranslator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBoundsExtentTranslator - Extent translation through lookup table.
// .SECTION Description
// vtkBoundsExtentTranslator provides a vtkExtentTranslator that is
// programmed with a specific extent corresponding to each piece
// number.  Readers can provide this to an application to allow the
// pipeline to execute using the same piece breakdown that is provided
// in the input file.

#ifndef __vtkBoundsExtentTranslator_h
#define __vtkBoundsExtentTranslator_h

#include "vtkExtentTranslator.h"
#include <vector>

class VTK_EXPORT vtkBoundsExtentTranslator : public vtkExtentTranslator
{
public:
  vtkTypeMacro(vtkBoundsExtentTranslator,vtkExtentTranslator);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  static vtkBoundsExtentTranslator* New();
  
  // Description:

  // Set the number of pieces into which the whole extent will be
  // split.  If this is 1 then the whole extent will be returned.  If
  // this is more than the number of pieces in the table then the
  // extra pieces will be empty data.  If this is more than one but
  // less than the number of pieces in the table then only this many
  // pieces will be returned (FIXME).
  void SetNumberOfPieces(int pieces);

  // Description:
  // Not supported by this subclass of vtkExtentTranslator.
  int PieceToExtent();  
  int PieceToExtentByPoints();
  int PieceToExtentThreadSafe(int piece, int numPieces, 
                              int ghostLevel, int *wholeExtent, 
                              int *resultExtent, int splitMode, 
                              int byPoints);
  
  // Description:
  // Set the extent to be used for a piece.  This sets the extent table
  // entry for the piece.
  virtual void SetBoundsForPiece(int piece, double* bounds);
  
  // Description:  
  // Get the extent table entry for the given piece.  This is only for
  // code that is setting up the table.  Extent translation should
  // always be done through the PieceToExtent method.
  virtual void GetBoundsForPiece(int piece, double* bounds);
  virtual double *GetBoundsForPiece(int piece);
  
  // Description:
  // Set the maximum ghost overlap region that is required 
  vtkSetMacro(MaximumGhostDistance, double);
  vtkGetMacro(MaximumGhostDistance, double);
  
protected:
  vtkBoundsExtentTranslator();
  ~vtkBoundsExtentTranslator();
  
  // Store the extent table in a single array.  Every 6 values form an extent.
  std::vector<double> BoundsTable;
  double MaximumGhostDistance;
   
private:
  vtkBoundsExtentTranslator(const vtkBoundsExtentTranslator&);  // Not implemented.
  void operator=(const vtkBoundsExtentTranslator&);  // Not implemented.
};

#endif
