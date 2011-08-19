#include "KernelGaussian.h"

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
KernelGaussian::w(double h, double distance) const
{
  double Hinverse = 1.0/h;
  double factorW = norm * pow(Hinverse, this->dim);
  double Q = distance * Hinverse;

  return factorW * exp(-Q*Q);
}
//----------------------------------------------------------------------------
Vector
KernelGaussian::gradW(double distance, const Vector& distanceVector) const
{
  double Q = distance * this->Hinverse;
  if (Q != 0.0)
  {
    return factorGradW * exp(-Q * Q) * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
Vector KernelGaussian::gradW(double h, double distance, const Vector& distanceVector) const
{
  double Hinverse = 1.0/h;
  double factorGradW = - 2.0 * pow(Hinverse, this->dim+1) / pow(M_PI, static_cast<double>(dim) / 2.0);
  double Q = distance * Hinverse;

  //!!! check this due to the fitting offset
  if (distance != 0.0)
  {
    return factorGradW * exp(-Q * Q) * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
double KernelGaussian::maxDistance() const
{
  return maxRadius;
}
