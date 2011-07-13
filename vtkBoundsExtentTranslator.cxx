/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBoundsExtentTranslator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkBoundsExtentTranslator.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkBoundsExtentTranslator);

//----------------------------------------------------------------------------
vtkBoundsExtentTranslator::vtkBoundsExtentTranslator()
{
  this->MaximumGhostDistance = 0;
}

//----------------------------------------------------------------------------
vtkBoundsExtentTranslator::~vtkBoundsExtentTranslator()
{
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  if(!this->BoundsTable.empty())
    {
    vtkIndent nextIndent = indent.GetNextIndent();
    double* bounds = &this->BoundsTable[0];
    int i;
    
    os << indent << "BoundsTable: 0: "
       << bounds[0] << " " << bounds[1] << " "
       << bounds[2] << " " << bounds[3] << " "
       << bounds[4] << " " << bounds[5] << "\n";
    for(i=1; i<this->GetNumberOfPieces();++i)
      {
      bounds += 6;
      os << nextIndent << "             " << i << ": "
         << bounds[0] << " " << bounds[1] << " "
         << bounds[2] << " " << bounds[3] << " "
         << bounds[4] << " " << bounds[5] << "\n";
      }
    }
  else
    {
    os << indent << "BoundsTable: (none)\n";
    }
  os << indent << "MaximumGhostDistance: " << this->MaximumGhostDistance << "\n";
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetNumberOfPieces(int pieces)
{
  // Allocate a table for this number of pieces.
  this->BoundsTable.resize(pieces*6);
  this->Superclass::SetNumberOfPieces(pieces);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetBoundsForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(&this->BoundsTable[piece*6], bounds, sizeof(double)*6);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::GetBoundsForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(bounds, &this->BoundsTable[piece*6], sizeof(double)*6);
}

//----------------------------------------------------------------------------
double* vtkBoundsExtentTranslator::GetBoundsForPiece(int piece)
{
  static double emptyBounds[6] = {0,-1,0,-1,0,-1};
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return emptyBounds;
    }
  return &this->BoundsTable[piece*6];
}

//----------------------------------------------------------------------------
// Make sure these inherited methods report an error is anyone is calling them
//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtentByPoints()
{
  vtkErrorMacro("PieceToExtentByPoints not supported.");
  return 0;
}

//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtent()
{
  return 
    this->PieceToExtentThreadSafe(this->Piece, this->NumberOfPieces,
                                  this->GhostLevel, this->WholeExtent,
                                  this->Extent, this->SplitMode, 0);
}

//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtentThreadSafe(int piece, int numPieces, 
                                                 int ghostLevel, 
                                                 int *wholeExtent, 
                                                 int *resultExtent, 
                                                 int splitMode, 
                                                 int byPoints)
{
  memcpy(resultExtent, wholeExtent, sizeof(int)*6);   
  return 1;
}

