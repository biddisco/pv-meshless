#include "KernelCusp.h"

//----------------------------------------------------------------------------
KernelCusp::KernelCusp(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
  // initialize the auxiliary factors
  Hinverse = 1.0/smoothingLength;
  maxRadius = 1.0*smoothingLength;

  // set the dimension dependent auxiliary factors
  switch (dim) {
    case 2:
      norm = 6.0 / M_PI;
      break;
    case 3:
      norm = 7.5 / M_PI;
      break;
    default:
      // Error: dim not in [1..3]
      cerr << "FAILURE: Kernel gets the wrong dimension: " << dim << "!" << endl;
      exit(1);
  }
  factorW = norm * pow(Hinverse, int(dim));
  factorGradW = - 2.0 * norm / pow(smoothingLength, int(dim+1));
  // a comment on the gradW sign:
  // In Speith (3.21) gradW has a minus in front!
  // This is the one before the 2.0 * ...
}
//----------------------------------------------------------------------------
double KernelCusp::w(double distance) const
{
    // dist/smoothingLength is often needed
    double Q = distance * Hinverse;

    // determine, were we are
    if (Q <= 1.0)
    {
        // we are inside the kernel
        double aux = 1.0 - Q;
        return factorW * aux * aux;
    }
    else {
        return 0.0;
    }
}
//----------------------------------------------------------------------------
double KernelCusp::w(double h, double distance) const
{
  // dist/smoothingLength is often needed
  double Q = distance * 1.0/h;
  // determine, were we are
  if (Q <= 1.0)
  {
    double aux = 1.0 - Q;
    return factorW * aux * aux;
  }
  else {
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector
KernelCusp::gradW(double distance, const Vector& distanceVector) const
{
    // dist/smoothingLength is often needed
    double normedDist = distance * Hinverse;

    // determine, were we are
    if ((normedDist <= 1.0) && (distance != 0.0))
    {
        // we are inside the kernel
        return (factorGradW * (Hinverse - 1.0 / distance)) * distanceVector;
    }
    // dist is bigger then the kernel.
    // Normaly, we should have tested this before, 
    // so the kernel will only be called for interacting particle pairs!
    // (that's the reason we put this condition past the others)
    //
    // and handle the case distance == 0 extra to avoid infty values!
    else
    {
        return Vector(0.0);
    }
}
//----------------------------------------------------------------------------
Vector KernelCusp::gradW(double h, double distance, const Vector& distanceVector) const
{
  double normedDist = distance * 1.0/h;
  if ((normedDist <= 1.0) && (distance != 0.0))
  {
    return (factorGradW * (Hinverse - 1.0 / distance)) * distanceVector;
  }
  else {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
double KernelCusp::maxDistance() const
{
  return maxRadius;
}
