// filename:
//    BetaSpline.cpp
// description:
//    A concrete kernel class.
// authors:
//    Sven Ganzenmueller <ganzenmu@informatik.uni-tuebingen.de>
//    Frank Heuser <heuserf@informatik.uni-tuebingen.de>
// last modified:
//    $Date: 2005/03/08 08:28:41 $
// project:
//    sph2000
// filetype:
//    c++-source
//

#include "BetaSpline.hpp"

//----------------------------------------------------------------------------
// Constructor initializes auxiliary values, so the calculation is faster
BetaSpline::BetaSpline(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
    // initialize the auxiliary factors
    reciprocH = 1.0 / smoothingLength;
    
    // set the dimension dependent auxiliary factors
    switch (dim)
    {
/*        case 1:
            norm = 4.0 / 3.0;
            break;
*/
        case 2:
            norm = 40.0 / 7.0 / pi;
            break;
        case 3:
            norm = 8.0 / pi;
            break;
        default:
            // Error: dim not in [1..3]
            cerr << "FAILURE: Kernel gets the wrong dimension: " << dim << "!" << endl;
            exit(1);
    }
      
    factorW     = norm * pow(reciprocH, int(dim));
    factorGradW = - 6.0 * norm * pow(reciprocH, int(dim+1));
    
    // a comment on the gradW sign:
    // In Speith (3.21) gradW has a minus in front!
    // This is the one before the 6.0 * ...
}
//----------------------------------------------------------------------------
// Calculates the kernel value for the given distance of two particles
// We take this from Speith, p.78, who references Monaghan & Lattenzio (1985)
// but used a doubled smoothing length for the definition of the interaction radius.
double
BetaSpline::w(double distance) const
{
    // dist/smoothingLength is often needed
    double normedDist = distance * reciprocH;

    // the beta-spline is composed of three functions, so we must determine, were we are
    if (normedDist < 0.5) 
    // we are in the inner region of the kernel
    {
        return factorW * (6.0 * normedDist * normedDist * normedDist - 
                            6.0 * normedDist * normedDist + 1.0);
    }
    else if (normedDist <= 1.0)
    // we are in the outer region of the kernel (not outside!)
    {
        double aux = 1.0 - normedDist;
        return factorW * ( 2.0 * aux * aux * aux );
    }
    else
    // dist is bigger then the kernel.
    // Normaly, we should have tested this before, 
    // so the kernel will only be called for interacting particle pairs!
    // (that's the reason we put this condition past the others)
    {
#       if RUN >= RUNALL
            cout << "WARNING: Called kernel for distance > smoothing length!\n";
#       endif
        return 0.0;
    }
}
//----------------------------------------------------------------------------
// Calculates the kernel derivation for the given distance of two particles
// We take this from Speith, p.78, who references Monaghan & Lattenzio (1985)
// but used a doubled smoothing length for the definition of the interaction radius.
Vector<double>
BetaSpline::gradW(double distance, const Vector<double>& distanceVector) const
{
    // dist/smoothingLength is often needed
    double normedDist = distance * reciprocH;

    // the beta-spline is composed of three functions (so the derivate is also), 
    // we must determine, were we are
    if (normedDist < 0.5) 
    // we are in the inner region of the kernel
    {
        return (factorGradW * (3.0 * normedDist * reciprocH - 2.0 * reciprocH)) * distanceVector;
    }
    else if (normedDist <= 1.0)
    // we are in the outer region of the kernel (not outside!)
    {
        return (- factorGradW * (1.0 / distance - 2.0 * reciprocH + normedDist * reciprocH)) *
                   distanceVector;
    }
    else
    // dist is bigger then the kernel.
    // Normaly, we should have tested this before, 
    // so the kernel will only be called for interacting particle pairs!
    // (that's the reason we put this condition past the others)
    {
#       if RUN >= RUNALL
            cout << "WARNING: Called kernel gradient for distance > smoothing length!\n";
#       endif
        return Vector<double>(0.0);
    }
}
//----------------------------------------------------------------------------
