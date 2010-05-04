// filename:
//    Spline35.hpp
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
//    c++-class
//



// ***** global macros (prevent multiple inclusion) *****
#ifndef SPLINE35_HPP
#define SPLINE35_HPP



// ***** system includes *****



// ***** project includes *****
#include "Kernel.hpp"
#include "Vector.hpp"



// ***** forward references *****



/**
 * 3 splines of 5.th order fitted as a kernel.
 * See Speith (1998), p.81.
 */
class Spline35 : public Kernel
{

public:

    /**
     * Constructor to initialize the data members and
     * pre-calculate some factors to save time when calling w() and gradW().
     */
    Spline35(int dim, double smoothingLength);

    /**
     * Calculates the kernel value for the given distance of two particles. 
     * The used formula is Speith (3.135).
     */
    virtual double w(double distance) const;

    /**
     * Calculates the kernel value for the given distance of two particles
     * with a variable smoothing length h.
     */
    virtual double w(double distance, double h) const;

    /**
     * Calculates the kernel derivation for the given distance of two particles. 
     * The used formula is the derivation of Speith (3.135) for the value 
     * with (3.21) for the direction of the gradient vector.
     * Be careful: grad W is antisymmetric in r (3.25)!.
     */
    virtual Vector<double> gradW(double distance, const Vector<double>& distanceVector) const;

    /**
     * Calculates the kernel derivation for the given distance of two particles
     * with a variable smoothing length h.
     */
    virtual Vector<double> gradW(double distance, const Vector<double>& distanceVector, double h) const;

private:

    /** The assignment operator.
     * This class does not need an assignment operator.
     * By declaring it privat, without defining it
     * (no implementation in the cpp source file),
     * we prevent the compiler from generating one
     * (See Scott Meyers, Effective C++, Lections 11+27 for more details!)
     */
    Spline35& operator=(const Spline35& non);

    /** Normalization factor */
    double norm;

    /** Auxiliary factors for intermediate results: The inverse smoothing length */
    double reciprocH;

    /** Auxiliary factors for intermediate results: A pre-factor for w */
    double factorW;

    /** Auxiliary factors for intermediate results: A pre-factor for grad w */
    double factorGradW;
};
#endif
// SPLINE35_HPP

