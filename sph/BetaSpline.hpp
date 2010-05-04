// filename:
//    BetaSpline.hpp
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
#ifndef BETASPLINE_HPP
#define BETASPLINE_HPP



// ***** system includes *****



// ***** project includes *****
#include "Kernel.hpp"
#include "Vector.hpp"



// ***** forward references *****



/**
 * Cubic beta-spline as kernel.
 * See Speith (1998), p.78 or
 * Monaghan & Lattanzio (1985).
 *  
 * Most often used kernel.
 */
class BetaSpline : public Kernel
{

public:

    /**
     * Constructor to initialize the data members and
     * pre-calculate some factors to save time when calling w() and gradW().
     */
    BetaSpline(int dim, double smoothingLength);

    /**
     * Calculates the kernel value for the given distance of two particles. 
     * The used formula is Speith (3.127).
     */
    virtual double w(double distance) const;

    /**
     * Calculates the kernel derivation for the given distance of two particles. 
     * The used formula is Speith (3.129) for the value with (3.21) for the direction of the
     * gradient vector.
     * Be careful: grad W is antisymmetric in r (3.25)!.
     */
    virtual Vector<double> gradW(double distance, const Vector<double>& distanceVector) const;

private:

    /** The assignment operator.
     * This class does not need an assignment operator.
     * By declaring it privat, without defining it
     * (no implementation in the cpp source file),
     * we prevent the compiler from generating one
     * (See Scott Meyers, Effective C++, Lections 11+27 for more details!)
     */
    BetaSpline& operator=(const BetaSpline& non);

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
// BETASPLINE_HPP

