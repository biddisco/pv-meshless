/*=========================================================================

  Project                 : pv-meshless
  Module                  : $RCSfile: vtkParticleBoxTree.cpp,v $
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
//
#include <stack>
#include <vector>
//
#include "vtkObjectFactory.h"
#include "vtkParticleBoxTree.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkGenericCell.h"
#include "vtkMath.h"
//
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticleBoxTree);
vtkCxxSetObjectMacro(vtkParticleBoxTree, ParticleSizeArray, vtkDataArray);
vtkCxxSetObjectMacro(vtkParticleBoxTree, ParticleBoundsArray, vtkDataArray);
//----------------------------------------------------------------------------
#define JB_DEBUG__
#if defined JB_DEBUG__
  #ifdef WIN32
    #define OUTPUTTEXT(a) vtkOutputWindowDisplayText(a);
  #else
    #define OUTPUTTEXT(a) std::cout << (a); std::cout.flush();
  #endif

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << "" a << "\n"; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }

  #undef vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
vtkParticleBoxTree::vtkParticleBoxTree(void) {
  this->ParticleSize = 0.05;
  this->ParticleSizeArray = NULL;
  this->ParticleBoundsArray = NULL;
}
//---------------------------------------------------------------------------
vtkParticleBoxTree::~vtkParticleBoxTree(void) {
  this->SetParticleSizeArray(NULL);
  this->SetParticleBoundsArray(NULL);
}
//----------------------------------------------------------------------------
bool vtkParticleBoxTree::StoreCellBounds()
{
  if (this->CellBounds) return false;
  if (!this->DataSet) return false;
  // Allocate space for cell bounds storage, then fill
  vtkIdType numCells = this->DataSet->GetNumberOfCells();
  this->CellBounds = new double [numCells][6];
  //
  double size = this->ParticleSize;
  for (vtkIdType j=0; j<numCells; j++) 
  { 
    if (this->ParticleBoundsArray) {
      this->ParticleBoundsArray->GetTuple(j, &this->CellBounds[j][0]);
    }
    else {
      if (this->ParticleSizeArray) {
        size = this->ParticleSize*sqrt(this->ParticleSizeArray->GetTuple1(j));
      }
      this->DataSet->GetCellBounds(j, CellBounds[j]);
      for (int i=0; i<3; i++) {
        this->CellBounds[j][i*2+0] -= size/2.0;    
        this->CellBounds[j][i*2+1] += size/2.0;    
      }
    }
  }
  return true;
}
//----------------------------------------------------------------------------
// this routine assumes the sphere centred at origin, 
// so transform  ray origin before entering here.
bool intersect_sphere(double o[3], double d[3], double r, double &t)
{
  //Compute A, B and C coefficients
  float a = vtkMath::Dot(d,d);
  float b = 2 * vtkMath::Dot(d, o);
  float c = vtkMath::Dot(o,o) - (r * r);

  //Find discriminant
  float disc = b * b - 4 * a * c;
    
  // if discriminant is negative there are no real roots, so return 
  // false as ray misses sphere
  if (disc < 0)
      return false;

  // compute q as described above
  float distSqrt = sqrtf(disc);
  float q;
  if (b < 0)
      q = (-b - distSqrt)/2.0;
  else
      q = (-b + distSqrt)/2.0;

  // compute t0 and t1
  float t0 = q / a;
  float t1 = c / q;

  // make sure t0 is smaller than t1
  if (t0 > t1)
  {
    // if t0 is bigger than t1 swap them around
    float temp = t0;
    t0 = t1;
    t1 = temp;
  }

  // if t1 is less than zero, the object is in the ray's negative direction
  // and consequently the ray misses the sphere
  if (t1 < 0)
    return false;

  // if t0 is less than zero, the intersection point is at t1
  if (t0 < 0)
  {
    t = t1;
    return true;
  }
  // else the intersection point is at t0
  else
  {
    t = t0;
    return true;
  }
}
//---------------------------------------------------------------------------
int vtkParticleBoxTree::IntersectCellInternal(
  vtkIdType cell_ID, double p1[3], double p2[3], 
  double tol, double &t, double ipt[3], double pcoords[3], int &subId)
{
  double ctmin=0, ctmax=1;
  double ray_vec[3] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
  // radius of sphere
  double r = CellBounds[cell_ID][1] - CellBounds[cell_ID][0];
  // shift ray because sphere eqn is expecting centre on origin
  double mid[3] = { (CellBounds[cell_ID][0] + CellBounds[cell_ID][1])/2.0, 
                    (CellBounds[cell_ID][2] + CellBounds[cell_ID][3])/2.0, 
                    (CellBounds[cell_ID][4] + CellBounds[cell_ID][5])/2.0 };
  double ray_origin[3] = { p1[0]-mid[0], p1[1]-mid[1], p1[2]-mid[2] };
  if (intersect_sphere(ray_origin, ray_vec, r/2.0, t)) {
    t = vtkMath::Distance2BetweenPoints(p1, mid)/vtkMath::Distance2BetweenPoints(p1, p2);
    t = sqrt(t);
    ipt[0] = p1[0]+t*ray_vec[0];
    ipt[1] = p1[1]+t*ray_vec[1];
    ipt[2] = p1[2]+t*ray_vec[2];
    return 1;
  }
  return 0;
/*
  // we already did this in the IntersectWithLine fn, but we will
  // now set the tmin and intersection values. We use Tmin = centre of box
  // to avoid problems of overlapping cubes
  if (BSPNode::RayMinMaxT(CellBounds[cell_ID], p1, ray_vec, ctmin, ctmax)) {
    double mid[3] = { (CellBounds[cell_ID][0] + CellBounds[cell_ID][1])/2.0, 
                      (CellBounds[cell_ID][2] + CellBounds[cell_ID][3])/2.0, 
                      (CellBounds[cell_ID][4] + CellBounds[cell_ID][5])/2.0 };
    double tt[3] = {  (mid[0] - p1[0]) / ray_vec[0],
                      (mid[1] - p1[1]) / ray_vec[1],
                      (mid[2] - p1[2]) / ray_vec[2] };

    t = vtkMath::Distance2BetweenPoints(p1, mid)/vtkMath::Distance2BetweenPoints(p1, p2);
    t = sqrt(t);
    ipt[0] = p1[0]+t*ray_vec[0];
    ipt[1] = p1[1]+t*ray_vec[1];
    ipt[2] = p1[2]+t*ray_vec[2];
    return 1;
  }
  return 0;
  */

}
//---------------------------------------------------------------------------
class _box {
  public:
  double bounds[6];
  _box(double *b) { for (int i=0; i<6; i++) { bounds[i] = b[i]; } };
};

typedef std::vector<_box> boxlist;
typedef std::stack<BSPNode*, std::vector<BSPNode*> > nodestack;
//---------------------------------------------------------------------------
void vtkParticleBoxTree::GenerateRepresentation(int level, vtkPolyData *pd)
{
  nodestack ns;
  boxlist   bl;
  BSPNode   *node;
  this->BuildLocatorIfNeeded();
  ns.push(mRoot);
  // lets walk the tree and get all the level n node boxes
  while (!ns.empty())  {
    node = ns.top();
    ns.pop();
    if (node->depth==level) bl.push_back(_box(node->bounds));
    else {
      if (node->mChild[0]) {
        ns.push(node->mChild[0]);
        if (node->mChild[1]) ns.push(node->mChild[1]);
        ns.push(node->mChild[2]);
      }
      else if (level==-1) bl.push_back(_box(node->bounds));
    }
  }
  //
  //
  //
  // For each node, add the bbox to our polydata
  size_t s = bl.size();
  for (size_t i=0; i<s; i++) {
    this->AddBox(pd, bl[i].bounds);
  }
}
//---------------------------------------------------------------------------
void vtkParticleBoxTree::GenerateRepresentationLeafs(vtkPolyData *pd)
{
  GenerateRepresentation(-1,pd);
}
//---------------------------------------------------------------------------
void vtkParticleBoxTree::GenerateRepresentationParticles(vtkPolyData *pd)
{
  /*
  //
  this->BuildLocatorIfNeeded();
  //
  nodestack ns;
  BSPNode   *node;
  ns.push(mRoot);
  double closestPoint[3], dist2;
  int subId;
  //
  while (!ns.empty())  {
    node = ns.top();
    ns.pop();
    if (node->mChild[0]) { // this must be a parent node    
      if (node->mChild[0]->Inside(x)) ns.push(node->mChild[0]);
      if (node->mChild[1] && node->mChild[1]->Inside(x)) ns.push(node->mChild[1]);
      if (node->mChild[2]->Inside(x)) ns.push(node->mChild[2]);
    }
    else { // a leaf, so test the cells
      for (int i=0; i<node->num_cells; i++) {
        int cell_ID = node->sorted_cell_lists[0][i];
        //
        if (vtkModifiedBSPTree_Inside(CellBounds[cell_ID], x)) {
          this->DataSet->GetCell(cell_ID, cell);
          if (cell->EvaluatePosition(x, closestPoint, subId, pcoords, dist2, weights)==1) {
            return cell_ID;
          }
//          if (dist2<tol2) return cell_ID;
        }
      }
    }
  }
  return -1;
*/
}
//---------------------------------------------------------------------------
void vtkParticleBoxTree::AddBox(vtkPolyData *pd, double *bounds)
{
  vtkPoints      *pts = pd->GetPoints();
  vtkCellArray *lines = pd->GetLines();
  double x[3];
  vtkIdType cells[8], ids[2];
  //
  x[0] = bounds[0]; x[1] = bounds[2]; x[2] = bounds[4];
  cells[0] = pts->InsertNextPoint(x);
  x[0] = bounds[1]; x[1] = bounds[2]; x[2] = bounds[4];
  cells[1] = pts->InsertNextPoint(x);
  x[0] = bounds[0]; x[1] = bounds[3]; x[2] = bounds[4];
  cells[2] = pts->InsertNextPoint(x);
  x[0] = bounds[1]; x[1] = bounds[3]; x[2] = bounds[4];
  cells[3] = pts->InsertNextPoint(x);
  x[0] = bounds[0]; x[1] = bounds[2]; x[2] = bounds[5];
  cells[4] = pts->InsertNextPoint(x);
  x[0] = bounds[1]; x[1] = bounds[2]; x[2] = bounds[5];
  cells[5] = pts->InsertNextPoint(x);
  x[0] = bounds[0]; x[1] = bounds[3]; x[2] = bounds[5];
  cells[6] = pts->InsertNextPoint(x);
  x[0] = bounds[1]; x[1] = bounds[3]; x[2] = bounds[5];
  cells[7] = pts->InsertNextPoint(x);
  //
  ids[0] = cells[0]; ids[1] = cells[1];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[2]; ids[1] = cells[3];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[4]; ids[1] = cells[5];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[6]; ids[1] = cells[7];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[0]; ids[1] = cells[2];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[1]; ids[1] = cells[3];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[4]; ids[1] = cells[6];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[5]; ids[1] = cells[7];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[0]; ids[1] = cells[4];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[1]; ids[1] = cells[5];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[2]; ids[1] = cells[6];
  lines->InsertNextCell(2,ids);
  ids[0] = cells[3]; ids[1] = cells[7];
  lines->InsertNextCell(2,ids);
}
//---------------------------------------------------------------------------
void vtkParticleBoxTree::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
