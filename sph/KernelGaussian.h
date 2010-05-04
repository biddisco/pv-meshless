// filename:
//    KernelGaussian.hpp
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
//    c++-class
//



// ***** global macros (prevent multiple inclusion) *****
#ifndef GAUSS_HPP
#define GAUSS_HPP

// ***** project includes *****
#include "Kernel.h"
#include "Vector.hpp"

/**
 * A gaussian function as kernel.
 *  See Speith (1998), p.77.
 * The gauss kernel is non-compact, which means it is larger than 0 until infinity.
 * A calculation with this kernel would make all particles interact with all others
 * (a pow(n, 2) algorithm for the interaction search). This is something we want to prevent!
 * A possibility would be to cut the kernel for ordinates larger than the smoothing length;
 * to get a continuous kernel, an offset must be defined, fitting the cutted gauss to the
 * zero line. But the gradient would still be unsteady!
 * (This implementation only cuts the kernel without fitting).
 */
class KernelGaussian : public Kernel
{
  public:

    /**
     * Constructor to initialize the data members and
     * pre-calculate some factors to save time when calling w() and gradW().
     */
    KernelGaussian(int dim, double smoothingLength);

    /**
     * Calculates the kernel value for the given distance of two particles. 
     * The used formula is Speith (3.126).
     */
    virtual double w(double distance) const;
    virtual double w(double h, double distance);

    /**
     * Calculates the kernel derivation for the given distance of two particles. 
     * The used formula is the derivation of Speith (3.126) for the value
     * with (3.21) for the direction of the gradient vector.
     * Be careful: grad W is antisymmetric in r (3.25)!.
     */
    virtual Vector gradW(double distance, const Vector& distanceVector) const;
    virtual Vector gradW(double h, double distance, const Vector& distanceVector);

    /** return the maximum distance at which this kernel is non zero */
    virtual double maxDistance() const;

private:

    /** The assignment operator.
     * This class does not need an assignment operator.
     * By declaring it privat, without defining it
     * (no implementation in the cpp source file),
     * we prevent the compiler from generating one
     * (See Scott Meyers, Effective C++, Lections 11+27 for more details!)
     */
    KernelGaussian& operator=(const KernelGaussian& non);

    /** Auxiliary factors for intermediate results: a pi-factor */
    double norm;

    /** Auxiliary factors for intermediate results: The inverse smoothing length */
    double Hinverse;
    double maxRadius;

    /** Auxiliary factors for intermediate results: A pre-factor for w */
    double factorW;

    /** Auxiliary factors for intermediate results: A pre-factor for grad w */
    double factorGradW;
};
#endif
// GAUSS_HPP

