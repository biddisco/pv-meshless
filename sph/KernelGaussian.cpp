#include "KernelGaussian.h"

#define _USE_MATH_DEFINES
#include <math.h> 

// Constructor initializes auxiliary values, so the calculation is faster
KernelGaussian::KernelGaussian(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  this->Hinverse = 1.0/smoothingLength;
  this->maxRadius = 5.0*smoothingLength;
  
  // set the dimension dependent auxiliary factors
  norm        = 1.0 / pow(M_PI, static_cast<double>(this->dim) / 2.0);
  factorW     = norm * pow(this->Hinverse, this->dim);
  factorGradW = - 2.0 * pow(this->Hinverse, this->dim+1) / pow(M_PI, static_cast<double>(dim) / 2.0);
  // a comment on the gradW sign:
  // In Speith (3.21) gradW has a minus in front!
  // With the one from the derivation of exp(-...), we get a plus
}
//----------------------------------------------------------------------------
double
KernelGaussian::w(double distance) const
{
  // dist/smoothingLength is often needed
  double normedDist = distance * this->Hinverse;
  //!!! fit this, so that the kernel is steady (use an h dependant offset)
  return factorW * exp(-normedDist * normedDist);
}
//----------------------------------------------------------------------------
double
KernelGaussian::w(double h, double distance)
{
  this->Hinverse   = 1.0/h;
  double normedDist = distance*this->Hinverse;
  this->factorW = norm * pow(this->Hinverse, this->dim);

  //!!! fit this, so that the kernel is steady (use an h dependant offset)
  return factorW * exp(-normedDist * normedDist);
}
//----------------------------------------------------------------------------
Vector
KernelGaussian::gradW(double distance, const Vector& distanceVector) const
{
  // dist/smoothingLength is often needed
  double normedDist = distance * this->Hinverse;
  //!!! check this due to the fitting offset
  if (distance != 0.0)
  {
      return factorGradW * exp(-normedDist * normedDist) * distanceVector;
  }
  else
  {
      // for distance == 0, the distanceVector is a null vector
      return distanceVector;
  }
}
//----------------------------------------------------------------------------
Vector KernelGaussian::gradW(double h, double distance, const Vector& distanceVector)
{
  this->Hinverse   = 1.0/h;
  double normedDist = distance*this->Hinverse;
  this->factorGradW = - 2.0 * pow(this->Hinverse, this->dim+1) / pow(M_PI, static_cast<double>(dim) / 2.0);

  //!!! check this due to the fitting offset
  if (distance != 0.0)
  {
      return factorGradW * exp(-normedDist * normedDist) * distanceVector;
  }
  else
  {
      // for distance == 0, the distanceVector is a null vector
      return distanceVector;
  }
}
//----------------------------------------------------------------------------
double KernelGaussian::maxDistance() const
{
  return maxRadius;
}
