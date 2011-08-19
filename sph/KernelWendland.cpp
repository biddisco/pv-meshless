#include "KernelWendland.h"

//----------------------------------------------------------------------------
KernelWendland::KernelWendland(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  this->Hinverse = 1.0/smoothingLength;
  this->maxRadius = 2.0*smoothingLength;
  
  // set the dimension dependent auxiliary factors
  switch (dim)
  {
    case 2:
      norm = 7.0/(4.0*M_PI);
      break;
    case 3:
      norm = 21.0/(16.0*M_PI);
      break;
    default:
      // Error: dim not in [2..3]
      cerr << "ERROR: Kernel gets the wrong dimension: " << dim << "!" << endl;
      exit(1);
  }
  factorW     = norm * pow(Hinverse, this->dim);
  factorGradW = norm * pow(Hinverse, this->dim+2);
}
//----------------------------------------------------------------------------
double KernelWendland::w(double distance) const
{
  double Q = distance * Hinverse;
  if (Q<2.0) {
    return this->factorW *(pow(1.0-0.50*Q,4.0)*(2.0*Q+1.0));
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
double KernelWendland::w(double h, double distance) const
{
  double Hinverse = 1.0/h;
  double factorW = norm * pow(Hinverse, this->dim);
  double Q = distance * Hinverse;
  if (Q<2.0) {
    return factorW *(pow(1.0-0.50*Q,4.0)*(2.0*Q+1.0));
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector 
KernelWendland::gradW(double distance, const Vector& distanceVector) const
{
  double Q = distance * Hinverse;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<2.0) {
    return this->factorGradW * -5.0*pow(1.0-0.50*Q,3.0);
  }
  return Vector(0.0);
}
//----------------------------------------------------------------------------
Vector KernelWendland::gradW(double h, double distance, const Vector& distanceVector) const
{
  double Hinverse = 1.0/h;
  double factorGradW = norm * pow(Hinverse, this->dim+1);
  double Q = distance * Hinverse;
  if (Q==0.0) {
    return Vector(0.0);
  }
  else if (Q<2.0) {
    return factorGradW * -5.0*pow(1.0-0.50*Q,3.0);
  }
  return Vector(0.0);
}
//----------------------------------------------------------------------------
double KernelWendland::maxDistance() const
{
  return maxRadius;
}
//----------------------------------------------------------------------------


