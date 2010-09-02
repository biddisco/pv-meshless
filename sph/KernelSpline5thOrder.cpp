#include "KernelSpline5thOrder.h"

//----------------------------------------------------------------------------
KernelSpline5thOrder::KernelSpline5thOrder(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  this->Hinverse = 1.0/smoothingLength;
  this->maxRadius = 3.0*smoothingLength;
  
  // set the dimension dependent auxiliary factors
  switch (dim)
  {
    case 2:
      norm = 7.0/(478.0*M_PI);
      break;
    case 3:
      norm = 1.0/(120.0*M_PI);
      break;
    default:
      // Error: dim not in [2..3]
      cerr << "ERROR: Kernel gets the wrong dimension: " << dim << "!" << endl;
      exit(1);
  }
  factorW     = norm * pow(Hinverse, this->dim);
  factorGradW = norm * pow(Hinverse, this->dim+1);
}
//----------------------------------------------------------------------------
double KernelSpline5thOrder::w(double distance) const
{
  double Q = distance * Hinverse;
  if (Q<1.0) {
    return this->factorW *(pow(3.0-Q,5) - 6.0*pow(2.0-Q,5) + 15.0*pow(1.0-Q,5));
  }
  else if (Q<2.0) {
    return this->factorW *(pow(3.0-Q,5) - 6.0*pow(2.0-Q,5));
  }
  else if (Q<3.0) {
    return this->factorW *(pow(3.0-Q,5));
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
double KernelSpline5thOrder::w(double h, double distance)
{
  this->Hinverse = 1.0/h;
  this->factorW = norm * pow(this->Hinverse, this->dim);
  double Q = distance * this->Hinverse;

  if (Q<1.0) {
    return this->factorW *(pow(3.0-Q,5) - 6.0*pow(2.0-Q,5) + 15.0*pow(1.0-Q,5));
  }
  else if (Q<2.0) {
    return this->factorW *(pow(3.0-Q,5) - 6.0*pow(2.0-Q,5));
  }
  else if (Q<3.0) {
    return this->factorW *(pow(3.0-Q,5));
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector 
KernelSpline5thOrder::gradW(double distance, const Vector& distanceVector) const
{
  double Q = distance * Hinverse;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<1.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4) + 30.0*pow(2.0-Q,4) + 75.0*pow(1.0-Q,4)) * distanceVector;
  }
  else if (Q<2.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4) + 30.0*pow(2.0-Q,4)) * distanceVector;
  }
  else if (Q<3.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4)) * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
Vector 
KernelSpline5thOrder::gradW(double h, double distance, const Vector& distanceVector)
{
  this->Hinverse = 1.0/h;
  this->factorGradW = norm * pow(this->Hinverse, this->dim+1);
  double Q = distance * this->Hinverse;

  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<1.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4) + 30.0*pow(2.0-Q,4) + 75.0*pow(1.0-Q,4)) * distanceVector;
  }
  else if (Q<2.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4) + 30.0*pow(2.0-Q,4)) * distanceVector;
  }
  else if (Q<3.0) {
    return this->factorGradW * (-5.0*pow(3.0-Q,4)) * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
double KernelSpline5thOrder::maxDistance() const
{
  return maxRadius;
}
//----------------------------------------------------------------------------


