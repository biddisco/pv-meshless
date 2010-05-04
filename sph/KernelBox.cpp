#include "KernelBox.h"

#define _USE_MATH_DEFINES
#include <math.h> 

//----------------------------------------------------------------------------
Box::Box(int dim, double smoothingLength)
  : Kernel(dim, smoothingLength)
{
  Hinverse = 1.0 / smoothingLength;
  
  switch (dim) {
    case 2:
      this->ScaleFactor     = 1.0/M_PI;
      this->GradientScaleFactor = -3.0/M_PI;
      break;
    case 3:
      this->ScaleFactor     = (3.0/4.0)/M_PI;
      this->GradientScaleFactor = -3.0/M_PI;
      break;
    default:
      // Error: dim not in [1..3]
      cerr << "FAILURE: Kernel gets the wrong dimension: " << dim << "!" << endl;
      exit(1);
  }
  factorW     = this->ScaleFactor * pow(Hinverse, int(dim));
  factorGradW = this->GradientScaleFactor * pow(Hinverse, int(dim+1));
}
//----------------------------------------------------------------------------
double Box::w(double distance) const
{
  if (distance <= smoothingLength) 
  // we are in the region of the kernel
  {
    return factorW;
  }
  else
  // dist is bigger then the kernel.
  // Normaly, we should have tested this before, 
  // so the kernel will only be called for interacting particle pairs!
  // (that's the reason we put this condition past the others)
  {
#     if RUN >= RUNALL
      cout << "WARNING: Called kernel for distance > smoothing length!\n";
#     endif
    return 0.0;
  }
}
//----------------------------------------------------------------------------
Vector
Box::gradW(double distance, const Vector& distanceVector) const
{
  // determine, whether we are inside or ouside the kernel region
  if ((distance <= smoothingLength) && (distance != 0.0)) 
  // we are in the region of the kernel
  {
    return (factorGradW / distance) * distanceVector;
  }
  else
  // dist is bigger then the kernel.
  // Normaly, we should have tested this before, 
  // so the kernel will only be called for interacting particle pairs!
  // (that's the reason we put this condition past the others)
  //
  // and the special case distance == 0, which must be handled extra, not to get infty values!
  {
    return Vector(0.0);
  }
}
//----------------------------------------------------------------------------
