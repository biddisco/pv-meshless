#include "KernelSpline3rdOrder.h"

#define _USE_MATH_DEFINES
#include <math.h> 
//----------------------------------------------------------------------------
KernelSpline3rdOrder::KernelSpline3rdOrder(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  this->Hinverse = 1.0/smoothingLength;
  this->maxRadius = 2.0*smoothingLength;
  
  // set the dimension dependent auxiliary factors
  switch (this->dim)
  {
    case 2:
      this->norm = 10.0/(7.0*M_PI);
      break;
    case 3:
    default:
      this->norm = 1.0/M_PI;
      break;
  }
  factorW     = norm * pow(Hinverse, this->dim);
  factorGradW = norm * pow(Hinverse, this->dim+1);
}
//----------------------------------------------------------------------------
double KernelSpline3rdOrder::w(double distance) const
{
  double Q = distance * Hinverse;
  if (Q<1.0) {
    return this->factorW *(1.0 - (3.0/2.0)*Q*Q + (3.0/4.0)*Q*Q*Q);
  }
  else if (Q<2.0) {
    double q2 = (2.0-Q);
    return this->factorW * (1.0/4.0) * (q2*q2*q2);;
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
double KernelSpline3rdOrder::w(double h, double distance)
{
  this->Hinverse = 1.0/h;
  this->factorW = norm * pow(this->Hinverse, this->dim);
  double Q = distance * this->Hinverse;

  if (Q<1.0) {
    return this->factorW *(1.0 - (3.0/2.0)*Q*Q + (3.0/4.0)*Q*Q*Q);
  }
  else if (Q<2.0) {
    double q2 = (2.0-Q);
    return this->factorW * (1.0/4.0) * (q2*q2*q2);;
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector 
KernelSpline3rdOrder::gradW(double distance, const Vector& distanceVector) const
{
  double Q = distance * Hinverse;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<1.0) {
    return this->factorGradW * (-3.0*Q + (9.0/4.0)*Q*Q) * distanceVector;
  }
  else if (Q<2.0) {
    double q2 = (2.0-Q);
    return this->factorGradW * (-3.0/4.0)*q2*q2 * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
Vector 
KernelSpline3rdOrder::gradW(double h, double distance, const Vector& distanceVector)
{
  this->Hinverse = 1.0/h;
  this->factorGradW = norm * pow(this->Hinverse, this->dim+1);
  double Q = distance * this->Hinverse;

  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<1.0) {
    return this->factorGradW * (-3.0*Q + (9.0/4.0)*Q*Q) * distanceVector;
  }
  else if (Q<2.0) {
    double q2 = (2.0-Q);
    return this->factorGradW * (-3.0/4.0)*q2*q2 * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
double KernelSpline3rdOrder::maxDistance() const
{
  return maxRadius;
}
//----------------------------------------------------------------------------


