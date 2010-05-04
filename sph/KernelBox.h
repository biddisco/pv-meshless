
#ifndef _KERNELBOX_H
#define _KERNELBOX_H

#include "Kernel.h"
#include "Vector.hpp"

class Box : public Kernel {
  public:
    /**
     * Constructor : Initializes constants
     */
    Box(int dim, double smoothingLength);

    /**
     * Calculates kernel value given a particle separation
     */
    virtual double w(double distance) const;

    /**
     * Calculates kernel gradient given a particle separation
     */
    virtual Vector gradW(double distance, const Vector& distanceVector) const;

private:
    Box& operator=(const Box& non); // Not implemented.
    Box(const Box&);                // Not implemented.

    double Hinverse;
    double ScaleFactor;

    /** The Box' gradient norm (the same gradient magnitude for all inner points) */
    double GradientScaleFactor;

    /** The Box' norm multiplied with the dimensioned h */
    double factorW;

    /** The Box' gradient norm  multiplied with the dimensioned h*/
    double factorGradW;
};
#endif
// _KERNELBOX_H

