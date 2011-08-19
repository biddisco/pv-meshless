// filename:
//    KernelCusp.hpp
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

#ifndef CUSP_HPP
#define CUSP_HPP

#include "Kernel.h"
#include "Vector.hpp"

/**
 * A parabolic sector as kernel named "cusp".
 * See Speith (1998), p.81.
 * <br>
 * This Kernel has a well-developed peak in the center. 
 */
class VTK_EXPORT KernelCusp : public Kernel
{
  public:

    /**
     * Constructor to initialize the data members and
     * pre-calculate some factors to save time when calling w() and gradW().
     */
    KernelCusp(int dim, double smoothingLength);

    /**
     * Calculates the kernel value for the given distance of two particles. 
     * The used formula is Speith (3.133).
     */
    virtual double w(double distance) const;
    virtual double w(double h, double distance) const;

    /**
     * Calculates the kernel derivation for the given distance of two particles. 
     * The used formula is Speith (3.134) for the value with (3.21) for the direction of the
     * gradient vector.
     * Be careful: grad W is antisymmetric in r (3.25)!.
     */
    virtual Vector gradW(double distance, const Vector& distanceVector) const;
    virtual Vector gradW(double h, double distance, const Vector& distanceVector) const;

    /** return the maximum distance at which this kernel is non zero */
    virtual double maxDistance() const;

    /** return the multiplier between smoothing length and max cutoff distance */
    virtual double getDilationFactor() const { return 1.0; }

  private:

    /** The assignment operator.
     * This class does not need an assignment operator.
     * By declaring it privat, without defining it
     * (no implementation in the cpp source file),
     * we prevent the compiler from generating one
     * (See Scott Meyers, Effective C++, Lections 11+27 for more details!)
     */
    KernelCusp& operator=(const KernelCusp& non);

    /** Normalization factor */
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

