// filename:
//    Spline35.cpp
// description:
//    A concrete kernel class.
// authors:
//    Sven Ganzenmueller <ganzenmu@informatik.uni-tuebingen.de>
//    Frank Heuser <heuserf@informatik.uni-tuebingen.de>
// last modified:
//    $Date: 2005/03/08 08:28:42 $
// project:
//    sph2000
// filetype:
//    c++-source
//



// ***** system includes *****
#include <iostream>



// ***** project includes *****
#include "Spline35.hpp"



// ***** forward references *****



// ***** class Spline35 : public Kernel *****

// Constructor initializes auxiliary values, so the calculation is faster
Spline35::Spline35(int dim, double smoothingLength)
    : Kernel(dim, smoothingLength)
{
    // initialize the auxiliary factors
    reciprocH = 1.0 / smoothingLength;
    
    // set the dimension dependent auxiliary factors
    switch (dim)
    {
        case 1:
            norm = 243.0 / 40.0;
            break;
        case 2:
            norm = 15309.0 / 478.0 / pi;
            break;
        case 3:
            norm = 2187.0 / 40.0 / pi;
            break;
        default:
            // Error: dim not in [1..3]
            cerr << "ERROR: Kernel gets the wrong dimension: " << dim << "!" << endl;
            exit(1);
    }
    factorW     = norm * pow(reciprocH, static_cast<int>(dim));
    factorGradW = + 5.0 * norm * pow(reciprocH, static_cast<int>(dim+1));
    // a comment on the gradW sign:
    // In Speith (3.21) gradW has a minus in front!
    // With the one from the derivation of the Spline, we get a plus for the gradient
}



// Calculates the kernel value for the given distance of two particles
// We take this from Speith, p.78, who references Monaghan & Lattenzio (1985)
// but used a doubled smoothing length for the definition of the interaction radius.
double
Spline35::w(double distance) const
{
    // dist/smoothingLength is often needed
    double normedDist = distance * reciprocH;

    // the Spline35 is composed of three functions, so we must determine, were we are
    if (normedDist < 1.0/3.0) 
    // we are in the inner region of the kernel
    {
        return factorW * (pow(1.0 - normedDist, 5) - 6.0 * pow(2.0/3.0 - normedDist, 5) +
                          15.0 * pow(1.0/3.0 - normedDist, 5));
    }
    else if (normedDist < 2.0/3.0)
    // we are in an outer inner region of the kernel
    {
        return factorW * (pow(1.0 - normedDist, 5) - 6.0 * pow(2.0/3.0 - normedDist, 5));
    }
    else if (normedDist <= 1.0)
    // we are in the outer region of the kernel (not outside!)
    {
        return factorW * pow(1.0 - normedDist, 5);
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



double
Spline35::w(double distance, double h) const
{
    // determine, where we are
    if (distance <= h)
    {
        double hpow = pow(h,dim);
        return norm / hpow * pow((1.0 - distance/h), 2);
    }
    else
    {
#       if RUN >= RUNALL
            cout << "WARNING: Called kernel for distance > smoothing length!\n";
#       endif
        return 0.0;
    }
}



// Calculates the kernel derivation for the given distance of two particles
// We take this from Speith, p.78, who references Monaghan & Lattenzio (1985)
// but used a doubled smoothing length for the definition of the interaction radius.
Vector<double>
Spline35::gradW(double distance, const Vector<double>& distanceVector) const
{
    // dist/smoothingLength is often needed
    double normedDist = distance * reciprocH;

    // the Spline35 is composed of three functions, so we must determine, were we are
    if (distance == 0.0)
    {
        return distanceVector;
    }
    else if (normedDist < 1.0/3.0)
    // we are in the inner region of the kernel
    {
        return (factorGradW *
                (pow(1.0 - normedDist, 4) - 6.0 * pow(2.0/3.0 - normedDist, 4) +
                15.0 * pow(1.0/3.0 - normedDist, 4)) / distance) * distanceVector;
    }
    else if (normedDist < 2.0/3.0)
    // we are in an outer inner region of the kernel
    {
        return (factorGradW *
                (pow(1.0 - normedDist, 4) - 6.0 * pow(2.0/3.0 - normedDist, 4)) / distance) * distanceVector;
    }
    else if (normedDist <= 1.0)
    // we are in the outer region of the kernel (not outside!)
    {
        return (factorGradW * pow(1.0 - normedDist, 4) / distance) * distanceVector;
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



Vector<double>
Spline35::gradW(double distance, const Vector<double>& distanceVector, double h) const
{
    // determine, were we are
    double normedDist = distance / h;
    double hpow = pow(h, dim + 1);
    if (distance == 0.0)
    {
        return distanceVector;
    }
    else if (normedDist < 1.0/3.0)
    {
        return (-5.0 * norm /hpow *
                (pow(1.0 - normedDist, 4) - 6.0 * pow(2.0/3.0 - normedDist, 4) +
                15.0 * pow(1.0/3.0 - normedDist, 4)) / distance) * distanceVector;
    }
    else if (normedDist < 2.0/3.0)
    {
        return (-5.0 * norm /hpow *
                (pow(1.0 - normedDist, 4) - 6.0 * pow(2.0/3.0 - normedDist, 4)) / distance) * distanceVector;
    }
    else if (normedDist <= 1.0)
    {
        return (-5.0 * norm /hpow * pow(1.0 - normedDist, 4) / distance) * distanceVector;
    }
    else
    {
#       if RUN >= RUNALL
            cout << "WARNING: Called kernel gradient for distance > smoothing length!\n";
#       endif
        return Vector<double>(0.0);
    }
}



