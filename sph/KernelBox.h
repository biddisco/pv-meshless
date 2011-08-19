
#ifndef _KERNELBOX_H
#define _KERNELBOX_H

#include "Kernel.h"
#include "Vector.hpp"

class VTK_EXPORT KernelBox : public Kernel {
  public:
    /**
     * Constructor : Initializes constants
     */
    KernelBox(int dim, double smoothingLength);

    /**
     * Calculates kernel value given a particle separation
     */
    virtual double w(double distance) const;
    virtual double w(double h, double distance) const;

    /**
     * Calculates kernel gradient given a particle separation
     */
    virtual Vector gradW(double distance, const Vector& distanceVector) const;
    virtual Vector gradW(double h, double distance, const Vector& distanceVector) const;

    /** return the maximum distance at which this kernel is non zero */
    virtual double maxDistance() const;

private:
    KernelBox& operator=(const KernelBox& non); // Not implemented.
    KernelBox(const KernelBox&);                // Not implemented.

    double Hinverse;
    double ScaleFactor;
    double maxRadius;

    /** The Box' gradient norm (the same gradient magnitude for all inner points) */
    double GradientScaleFactor;

    /** The Box' norm multiplied with the dimensioned h */
    double factorW;

    /** The Box' gradient norm  multiplied with the dimensioned h*/
    double factorGradW;
};
#endif
// _KERNELBOX_H

