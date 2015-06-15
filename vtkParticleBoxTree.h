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
// .NAME vtkParticleBoxTree - Special BSP tree for particles with no size
//
// .SECTION Description
// vtkParticleBoxTree is a BSP tree which takes particle data and
// assumes a box size for each particle. The box size is either fixed
// for all particles, or if supplied, a support radius with one entry
// per particle can be used. Raytracing of simple points would 
// normally produce no output, using this BSP tree enables a finite 
// volume to be assigned to each particle and permits raycasting 
// operations to be performed on them.

//#ifndef _vtkParticleBoxTree_h
//#define _vtkParticleBoxTree_h

class vtkDataArray;

class VTK_EXPORT PARTICLE_BOX_TREE_NAME : public PARTICLE_BOX_TREE_BASE {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeMacro(PARTICLE_BOX_TREE_NAME, PARTICLE_BOX_TREE_BASE);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Construct with default maximum 32 cells per node. (average 16->31)
    static PARTICLE_BOX_TREE_NAME *New();

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
    
    // Description:
    // Supply an array of bounding boxes, one per particle
    // this allows asymmetrical boxes to be handled.
    // The array should have N*6 components of bounds[] entries
    virtual void SetParticleBoundsArray(vtkDataArray*);
    vtkGetObjectMacro(ParticleBoundsArray, vtkDataArray);
    
    // Description:
    // Generate BBox representation of Nth level
    virtual void GenerateRepresentation(int level, vtkPolyData *pd);

    // Description:
    // Generate BBox representation of all leaf nodes
    virtual void GenerateRepresentationLeafs(vtkPolyData *pd);

    // Description:
    // Generate BBox representation of all particles stored in the tree
    virtual void GenerateRepresentationParticles(vtkPolyData *pd);

  protected:
   PARTICLE_BOX_TREE_NAME();
  ~PARTICLE_BOX_TREE_NAME();

  // override the cell bounds caching so that we can add a real size to the box
  virtual bool StoreCellBounds();

//BTX
  // override the ray/cell test to give the test result from our fake
  // cell based on particle size instead of the real cell which is a vtkVertex
  // with zero dimension/size
  virtual int IntersectCellInternal(vtkIdType cell_ID, double p1[3], double p2[3], 
    double tol, double &t, double ipt[3], double pcoords[3], int &subId);
//ETX

  // Used internally when generating representation of nodes/particles
  void AddBox(vtkPolyData *pd, double *bounds, int level);

  double ParticleSize;
  vtkDataArray *ParticleSizeArray;
  vtkDataArray *ParticleBoundsArray;

private:
  PARTICLE_BOX_TREE_NAME(const PARTICLE_BOX_TREE_NAME&);  // Not implemented.
  void operator=(const PARTICLE_BOX_TREE_NAME&);      // Not implemented.
};

//#endif


