#include "KernelQuadratic.h"

#define _USE_MATH_DEFINES
#include <math.h> 
//----------------------------------------------------------------------------
KernelQuadratic::KernelQuadratic(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  Hinverse = 1.0/smoothingLength;
  maxRadius = 2.0*smoothingLength;
  
  // set the dimension dependent auxiliary factors
  switch (dim)
  {
    case 2:
      norm = (21.0/20.0) * 10.0/(7.0*M_PI);          
      break;
    case 3:
      norm = (15.0/16.0) * 1.0/M_PI; 
      break;
    default:
      // Error: dim not in [2..3]
      cerr << "ERROR: Kernel gets the wrong dimension: " << dim << "!" << endl;
      exit(1);
  }
  factorW     = norm * pow(Hinverse, static_cast<int>(dim));
  factorGradW = norm * pow(Hinverse, static_cast<int>(dim+1));
}
//----------------------------------------------------------------------------
double KernelQuadratic::w(double distance) const
{
  double Q = distance * Hinverse;
  if (Q<2.0) {
    return this->factorW * ((1.0/4.0)*Q*Q - Q + 1.0);
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
double
KernelQuadratic::w(double h, double distance)
{
  this->Hinverse = 1.0/h;
  this->factorW = norm * pow(this->Hinverse, this->dim);
  double Q = distance * this->Hinverse;

  if (Q<2.0) {
    return this->factorW * ((1.0/4.0)*Q*Q - Q + 1.0);
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector 
KernelQuadratic::gradW(double distance, const Vector& distanceVector) const
{
  double Q = distance * Hinverse;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<2.0) {
    return this->factorGradW * ((1.0/2.0)*Q - 1.0) * distanceVector;
  }
  return Vector(0.0);
}
//----------------------------------------------------------------------------
Vector 
KernelQuadratic::gradW(double h, double distance, const Vector& distanceVector)
{
  double Q = distance / h;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<2.0) {
    return this->factorGradW * ((1.0/2.0)*Q - 1.0) * distanceVector;
  }
  return Vector(0.0);
}
//----------------------------------------------------------------------------
double KernelQuadratic::maxDistance() const
{
  return maxRadius;
}
//----------------------------------------------------------------------------


