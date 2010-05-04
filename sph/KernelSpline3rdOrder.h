#ifndef SPLINE3RDORDER_HPP
#define SPLINE3RDORDER_HPP

#include "Kernel.h"
#include "Vector.hpp"

class KernelSpline3rdOrder : public Kernel
{

public:
    /**
     * Constructor to initialize the data members and
     * pre-calculate some factors to save time when calling w() and gradW().
     */
    KernelSpline3rdOrder(int dim, double smoothingLength);

    /**
     * Calculates the kernel value for the given distance of two particles. 
     */
    virtual double w(double distance) const;
    virtual double w(double h, double distance);

    /**
     * Calculates the kernel derivation for the given distance of two particles. 
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
    KernelSpline3rdOrder& operator=(const KernelSpline3rdOrder& non);

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


