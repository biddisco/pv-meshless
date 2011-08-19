#ifndef KERNEL_HPP
#define KERNEL_HPP 

#include <string>
#define _USE_MATH_DEFINES
#include <math.h> 

#include "vtkSystemIncludes.h"
#include "Vector.hpp"

class VTK_EXPORT Kernel
{
  public:

    /**
     * Conctructor to set the data members.
     * Concrete kernels will pre-calculate some factors in their constructor
     * to save time when calling w() and gradW().
     */
    Kernel(int dim, double smoothingLength)     
      : dim(dim), smoothingLength(smoothingLength) {};

    /** Base classes with virtual member functions should have a virtual destructor. */
    virtual ~Kernel();

    /** Calculates the kernel value at the given distance. */
    virtual double w(double distance) const = 0;

    /** Calculates the kernel value at the given distance using a variable h value */
    virtual double w(double h, double distance) const = 0;

    /** Calculates the kernel derivative at the given distance using a variable h value */
    virtual Vector gradW(double distance, const Vector& distanceVector) const = 0;
    virtual Vector gradW(double h, double distance, const Vector& distanceVector) const = 0;

    /** return the multiplier between smoothing length and max cutoff distance */
    virtual double getDilationFactor() const { return 2.0; }

    /** return the maximum distance at which this kernel is non zero */
    virtual double maxDistance() const = 0;

  protected:

    /**
     * A local copy of the dimension.
     */
    const int dim;

    /**
     * A local copy of the smoothing length.
     */
    const double smoothingLength;

  private:
    Kernel& operator=(const Kernel& non);

};
#endif


