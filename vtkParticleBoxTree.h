/*=========================================================================

  Project         : vtkCSCS
  Module          : $RCSfile: vtkParticleBoxTree.h,v $
  Revision of last commit : $Rev: 155 $
  Author of last commit   : $Author: biddisco $
  Date of last commit   : $Date:: 2006-07-13 10:23:31 +0200 #$

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkParticleBoxTree - Special BSP tree for dimensionless particles
//
// .SECTION Description
// vtkParticleBoxTree is a BSP tree which takes particle data and
// assumes a box size for each particle. The box size is either fixed
// for all particles, or if supplied, a 'support radius' array with one entry
// per particle can be used. Raytracing of simple points would 
// normally produce no output, using this tree enables a finite 
// volume to be assigned to each particle and permits raycasting
// operations to be performed on them.

#ifndef _vtkParticleBoxTree_h
#define _vtkParticleBoxTree_h

#include "vtkCellTreeLocator.h"
class vtkDataArray;

class VTK_EXPORT vtkParticleBoxTree : public vtkCellTreeLocator {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeRevisionMacro(vtkParticleBoxTree,vtkCellTreeLocator);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Construct with maximum 32 cells per node. (average 16->31)
    static vtkParticleBoxTree *New();

//BTX
    using vtkAbstractCellLocator::IntersectWithLine;
    using vtkAbstractCellLocator::FindClosestPoint;
    using vtkAbstractCellLocator::FindClosestPointWithinRadius;
//ETX

    // Description:
    // State the particle size as a fixed Box dimension - 
    // the size is one axis of a cube shaped particle
    // if a ParticleSize array is provided, this is ignored
    vtkSetMacro(ParticleSize,double);
    vtkGetMacro(ParticleSize,double);

    // Description:
    // Supply an array to be used for particle sizes.
    // one usually specifies an 'H' smoothing length array
    virtual void SetParticleSizeArray(vtkDataArray*);
    vtkGetObjectMacro(ParticleSizeArray, vtkDataArray);

    void GenerateRepresentation(int level, vtkPolyData *pd);
    
  protected:
   vtkParticleBoxTree();
  ~vtkParticleBoxTree();

  // override the cell bounds caching so that we can add a real size to the box
  virtual bool StoreCellBounds();

//BTX
  // override the ray/cell test to give the test result from our fake
  // cell based on particle size instead of the real cell which is a vtkVertex
  // with zero dimension/size
  virtual int IntersectCellInternal(vtkIdType cell_ID, double p1[3], double p2[3], 
    double tol, double &t, double ipt[3], double pcoords[3], int &subId);
//ETX

  double ParticleSize;
  vtkDataArray *ParticleSizeArray;

private:
  vtkParticleBoxTree(const vtkParticleBoxTree&);  // Not implemented.
  void operator=(const vtkParticleBoxTree&);      // Not implemented.
};

#endif


