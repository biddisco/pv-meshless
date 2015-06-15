/*=========================================================================

  Project                 : pv-meshless
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

#include "vtkParticleBoxTreeBSP.h";
#include "vtkParticleBoxTreeCell.h";
class vtkDataArray;

class VTK_EXPORT vtkParticleBoxTreeRepresentation : public vtkPolyDataAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeMacro(vtkParticleBoxTreeRepresentation,vtkPolyDataAlgorithm);

    // Description:
    // Create an instance of vtkParticleBoxTreeRepresentation
    static vtkParticleBoxTreeRepresentation *New();

    // Description:
    // Set BSP Tree type -
    // 0 : vtkCellTreeLocator
    // 1 : vtkModifiedBSPTree
    vtkSetMacro(TreeType,int);
    vtkGetMacro(TreeType,int);

    vtkSetMacro(Level, int);
    vtkGetMacro(Level, int);

    vtkSetMacro(MaxDepth, int);
    vtkGetMacro(MaxDepth, int);

    vtkSetMacro(MaxCellsPerNode, int);
    vtkGetMacro(MaxCellsPerNode, int);
    
    vtkSetMacro(ParticleSize, double);
    vtkGetMacro(ParticleSize, double);    

    // Description:
    // Supply an array to be used for particle sizes.
    // one usually specifies an 'H' smoothing length array
    vtkSetStringMacro(ParticleSizeArray);
    vtkGetStringMacro(ParticleSizeArray);

    // Description:
    // Supply an array of bounding boxes, one per particle
    // this allows asymmetrical boxes to be handled.
    // The array should have N*6 components of bounds[] entries
    vtkSetStringMacro(ParticleBoundsArray);
    vtkGetStringMacro(ParticleBoundsArray);

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

    int    TreeType;
    int    Level;
    int    MaxDepth;
    int    MaxCellsPerNode;
    double ParticleSize;
    char  *ParticleSizeArray;
    char  *ParticleBoundsArray;
  
    vtkSmartPointer<vtkParticleBoxTreeCell> BSPTree1;
    vtkSmartPointer<vtkParticleBoxTreeBSP>  BSPTree2;

private:
  // Not implemented.
  vtkParticleBoxTreeRepresentation(const vtkParticleBoxTreeRepresentation&);
  void operator=(const vtkParticleBoxTreeRepresentation&);
};

#endif

