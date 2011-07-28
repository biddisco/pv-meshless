/*=========================================================================

  Project                 : vtkCSCS
  Module                  : vtkParticleBoxTreeRepresentation.h
  Revision of last commit : $Rev: 155 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2006-07-13 10:23:31 +0200 #$

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This code is derived from an earlier work and is distributed
  with permission from, and thanks to

  ------------------------------------------
  Copyright (C) 2000-2004 John Biddiscombe
  Skipping Mouse Software Ltd,
  Blewbury, England
  ------------------------------------------

=========================================================================*/
// .NAME vtkParticleBoxTreeRepresentation - Create Polygonal Bars from an image histogram
//
// .SECTION Description
//
// .SECTION See Also
// vtkImageAccumulate

#ifndef _vtkBarPolyGenerator_h
#define _vtkBarPolyGenerator_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkParticleBoxTree;

class VTK_EXPORT vtkParticleBoxTreeRepresentation : public vtkPolyDataAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeRevisionMacro(vtkParticleBoxTreeRepresentation,vtkPolyDataAlgorithm);

    // Description:
    // Create an instance of vtkParticleBoxTreeRepresentation
    static vtkParticleBoxTreeRepresentation *New();

    vtkSetMacro(Level, int);
    vtkGetMacro(Level, int);

    vtkSetMacro(MaxDepth, int);
    vtkGetMacro(MaxDepth, int);

    vtkSetMacro(MaxCellsPerNode, int);
    vtkGetMacro(MaxCellsPerNode, int);
    
    vtkSetMacro(ParticleSize, double);
    vtkGetMacro(ParticleSize, double);    

  protected:
     vtkParticleBoxTreeRepresentation(void);
    ~vtkParticleBoxTreeRepresentation();
    //
    virtual int RequestInformation (vtkInformation *,
                                    vtkInformationVector **,
                                    vtkInformationVector *);
    //
    virtual int RequestData(vtkInformation *request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);
    //
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    int    Level;
    int    MaxDepth;
    int    MaxCellsPerNode;
    double ParticleSize;
  
    vtkSmartPointer<vtkParticleBoxTree> BSPTree;

private:
  // Not implemented.
  vtkParticleBoxTreeRepresentation(const vtkParticleBoxTreeRepresentation&);
  void operator=(const vtkParticleBoxTreeRepresentation&);
};

#endif

