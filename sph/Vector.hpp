// filename:
//    Vector.hpp
// description:
//    The geometric vector class.
// authors:
//    Sven Ganzenmueller <ganzenmu@informatik.uni-tuebingen.de>
//    Frank Heuser <heuserf@informatik.uni-tuebingen.de>
// last modified:
//    $Date: 2005/06/06 06:16:06 $
// project:
//    sph2000
// filetype:
//    c++-class
//



// ***** global macros (prevent multiple inclusion) *****
#ifndef VECTOR_HPP
#define VECTOR_HPP



// ***** system includes *****
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
//#include <tpo++.H>


// ***** forward references *****

class Vector;

using namespace std;

#define X 0
#define Y 1
#define Z 2

/**
 * Addition of two Vectors.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 */
Vector
operator+(const Vector& v, const Vector& w);

/**
 * Addition of a scalar and a Vector.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 * <br>
 * Be careful: This operation isn't mathematical correct, but for physicians it is useful
 * (so long as you know, that there have to be written in mind a '1-Vector' behind the double value).
 */
Vector
operator+(double d, const Vector& v);

/**
 * Addition of a Vector and a scalar.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 * <br>
 * Be careful: This operation isn't mathematical correct, but for physicians it is useful
 * (so long as you know, that there have to be written in mind a '1-Vector' behind the double value).
 */
Vector
operator+(const Vector& v, double d);

/**
 * Subtraction of two Vectors.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 */
Vector
operator-(const Vector& v, const Vector& w);

/**
 * Subtraction of a scalar from a Vector.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 * <br>
 * Be careful: This operation isn't mathematical correct, but for physicians it is useful
 * (so long as you know, that there have to be written in mind a '1-Vector' behind the double value).
 */
Vector
operator-(const Vector& v, double d);

/**
 * Subtraction of a Vector from a scalar.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 * <br>
 * Be careful: This operation isn't mathematical correct, but for physicians it is useful
 * (so long as you know, that there have to be written in mind a '1-Vector' behind the double value).
 */
Vector
operator-(double d, const Vector& v);

/**
 * Multiplication of a scalar with a Vector.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 */
Vector
operator*(double d, const Vector& v);

/**
 * Multiplication of a Vector with a scalar.
 * <br>
 * Because a member function always uses the parameter order: *this as 1.st operator,
 * only non-member function can be used for mixed operand types with mixed operand orders.
 */
Vector
operator*(const Vector& v, double d);

/** Division of a Vector by a scalar. */
Vector
operator/(const Vector& v, double d);

/** 
 * Calculates the magnitude of the given Vector if the Vector 
 * can't be adressed by v.magnitude().
 */
double
magnitude(const Vector& v);

/**
 * A geometric vector for physical vector quantities like the position, the velocity,
 * the velocity derivative or the kernel gradient.
 * <br>
 * The Initializer must set the dimension, before the first using of a Vector!
 * <br>
 * Implementation hints: <br>
 * One way would be to use the valarray<double> as base class, because most of the needed functionality is given there.
 * With the valarray::operator[ ] there would be an array-like access. <br>
 * (For an intuitive access we define index constants X, Y and Z (as enum in auxiliary.hpp).)
 * <br>
 * Another way (the first thought) is to use the design pattern &gt;&gt; Strategy &lt;&lt;:<br>
 * An abstract base Vector class with concret 2D- and 3D-Vectors.
 * But that would mean to use pointers to the base class, which doesn't fit in numeric particle simulations:<br>
 * Pointers need extra memory and extra time (because of the indirection) and that with a very big
 * amount of needed Vectors...)
 * (by the way, it doesn't looks intuitive, that the e.g. position of a particle is a 
 * pointer to a vector).
 * <br>
 * So the implementation uses a simple C-array of doubles.
 * (A possibility to be dimension independent would be to store a pointer to the array.)
 */

class Vector
{

public:

    class VectorIterator : public std::iterator <std::input_iterator_tag, double, std::ptrdiff_t, double*, double&>
    {

    protected:

        /** The Vector, the iterator refers to. */
        Vector& Vec;

        /** The current position pointer into the referenced Vector. */
        double* pid;



    public:

        typedef VectorIterator self;
        typedef input_iterator_tag iterator_category;
        typedef iterator_traits<VectorIterator>::value_type value_type;
        typedef iterator_traits<VectorIterator>::reference reference;
        typedef iterator_traits<VectorIterator>::pointer pointer;
        typedef iterator_traits<VectorIterator>::difference_type difference_type;



        /**
         * Constructor.
         * The new iterator points to vec.begin().
         */
        explicit VectorIterator(Vector& vec);

        /**
         * Constructor.
         * The new iterator points to the element of vec, adressed by pid.
         */
        VectorIterator(Vector& vec, double* pid);

        /**
         * *iter
         * provides read access to the actual element.
         */
        double& operator*();

        /**
         * iter1 = iter2
         * assignment operator.
         */
        self& operator=(const self& it);

        /**
         * ++iter (prefix)
         * steps forward,
         * returning the new position.
         */
        self& operator++();

        /**
         * iter1 == iter2
         * return whether the two iterators are equal.
         */
        bool operator==(const self& it);

        /**
         * iter1!= iter2
         * returns whether the two iterators are unequal.
         */
        bool operator!=(const self& it);

        /**
         * TYPE(iter)
         * copy constructor.
         */

    }; // end class VectorIterator



    /**
     * A constructor to create a Vector in valarray style: <br>
     * All elements are initialized with 'value'.
     * Used as a standard constructor (no argument) it sets all coordinates to 0.0.
     */
    Vector(double value = 0);

    inline Vector(double v1, double v2, double v3) {
      values[0] = v1; values[1] = v2; values[2] = v3;
    }

    /** The copy constructor to create (clone) a Vector from a given Vector. */
    Vector(const Vector& v);

    /** The copy constructor to create (clone) a Vector from a given Vector. */
/*    Vector<int>(const Vector& v);
*/

    VectorIterator begin();
    VectorIterator end();

    /** The access-operator for non-const Vectors. */
    double& operator[](int id);

    /** The access-operator for const Vectors. */
    const double& operator[](int id) const;

    /** The assignment operator to assign a Vector to the Vector object. */
    Vector& operator=(const Vector& v);

    /**
     * The assignment operator to assign a number to all components.
     * We need this, because we often need to reset a Vector to 0, but to assign a nullvector,
     * we need the dimension to generate the nullvector.
     * So this member sets all coordinates of a given Vector to the given value.
     * Be careful: This operation isn't mathematical correct, but for physicians it is useful
     * (so long as you know, that there have to be written in mind a '1-Vector' behind the double value).
     */
    Vector& operator=(double d);

    /** Add a Vector to the given one. */
    Vector& operator+=(const Vector& v);

    /** Subtract a Vector from the given one. */
    Vector& operator-=(const Vector& v);

    /** Divide a Vector with a double. */
    Vector& operator/=(double d);

    /** Negate a Vector. */
    Vector operator-() const;



    /**
     * The comparison operator for our Vector.
     * It returns true, if the two Vectors are equal, false, if they are unequal.
     */
    bool operator==(const Vector& v) const;



    /** The mathematical dot product of two vectors.
     * <br>
     * To distinguish between the cross and the dot product, this operation does not overload
     * C++ operators, but define this member function.
     * (Who could intuitive know, what (v * w) means? Dot- or Crossproduct?)
     */
    double scalarProduct(const Vector& v) const;



    /** Returns the smallest element of the Vector. */
    double min() const;

    /** Returns the largest element of the Vector. */
    double max() const;

    /** Calculates the distance to another Vector. */
    double distance(const Vector& v) const;

    /** Calculates the magnitude of the given Vector. */
    double magnitude() const;

    /**
     * Calculates the magnitude of the given Vector in the x-plane
     * due to the given origin.
     * This is needed to calculate the velocity profile of the jets.
     */
    double xPlaneMagnitude(const Vector& origin) const;



    /**
     * The dimension of the simulation.
     * To be initialized from the Initializer.
     */
    static int dim;

    /**
     * The compiled dimension of the Vector.
     */
    static const int constdim = 3;

    /**
     * The coordinates of the Vector.
     * We use a C array here to save time and space.
     * (Every other thing needs more space, so we make it fast by spending one extra 'double'
     * for 2D. If you really need memory for 2D than you have to change constdim from 3 to 2 and to recompile.
     * But be carefull: the resulting exe is only for 2D. The *.config has to define Dimension : 2.
     * Otherwise we must use an extra pointer to an variable array!)
     */
    double values[constdim];



};

/** TPO++: make the Vector sendable. */
//TPO_TRIVIAL(Vector);





inline
Vector::VectorIterator::VectorIterator(Vector& vec)
    : Vec(vec), pid(vec.values) // point to the beginning
{
}
    



inline
Vector::VectorIterator::VectorIterator(Vector& vec, double* pid)
    : Vec(vec), pid(pid)
{
}




inline
double&
Vector::VectorIterator::operator*()
{
    return *pid;
}




inline
 Vector::VectorIterator::self&
Vector::VectorIterator::operator=(const  Vector::VectorIterator::self& it)
{
    Vec = it.Vec;
    pid = it.pid;
    return *this;
}




inline
 Vector::VectorIterator::self&
Vector::VectorIterator::operator++()
{
    if (pid != Vec.values + Vector::dim) // test end()
    {
        ++pid;
    }
    return *this;
}



/*

inline
Vector::VectorIterator<double>::self
Vector::VectorIterator<double>::operator++(int)
{
    self tmp = *this;
    ++(*this);
    //!!! test end()
    return tmp;
}
*/




inline
bool
Vector::VectorIterator::operator==(const  Vector::VectorIterator::self& it)
{
    return pid == it.pid;
}




inline
bool
Vector::VectorIterator::operator!=(const  Vector::VectorIterator::self& it)
{
    return pid != it.pid;
}



// ***** Auxiliary functions and types for Vector *****



// ***** inline functions *****

// we must include the implementation of inline member functions in the hpp-file, because the compiler 
// needs them in all files, which include the class;
// but to write the code in the class definition above wouldn't be nice to read, so we place it here.
// Private member functions could be placed in the cpp file.



inline
Vector::Vector(double value)
{
    // Stroustrup 18.6.6
    fill_n(values, dim, value);
}




inline
Vector::Vector(const Vector& v)
{
    // Stroustrup 18.6.1
    copy(v.values, v.values + dim, values);
}



/*template <>
inline
Vector<int>::Vector(const Vector& v)
{
    for (int i = 0; i < Vector::dim; ++i)
    {
        values[i] = static_cast<int>(v.values[i]);
    }
}
*/



inline
 Vector::VectorIterator
Vector::begin()
{
    return VectorIterator(*this, values);
}




inline
 Vector::VectorIterator
Vector::end()
{
    return VectorIterator(*this, values + dim);
}




inline
double&
Vector::operator[](int id)
{
    return values[id];
}




inline
const double&
Vector::operator[](int id) const
{
    return values[id];
}




inline
Vector&
Vector::operator=(const Vector& v)
{
    // Stroustrup 18.6.1
    copy(v.values, v.values + dim, values);
    return *this;
}




inline
Vector&
Vector::operator=(double d)
{
    // Stroustrup 18.6.6
    fill_n(values, dim, d);
    return *this;
}




inline
Vector&
Vector::operator+=(const Vector& v)
{
    // Stroustrup 18.4.3 and 18.6.2
    transform(values, values + Vector::dim, v.values, values, plus<double>());
    return *this;
}




inline
Vector&
Vector::operator-=(const Vector& v)
{
    // Stroustrup 18.4.3 and 18.6.2
    transform(values, values + Vector::dim, v.values, values, minus<double>());
    return *this;
}




inline
Vector&
Vector::operator/=(double d)
{
    Vector divisor(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(values, values + Vector::dim, divisor.values, divisor.values, divides<double>());
    return *this;
}




inline
Vector
Vector::operator-() const
{
    Vector res(*this);
    for (VectorIterator i = res.begin(); i != res.end(); ++i)
    {
        *i = -(*i);
    }
//!!! does not work:    for_each(res.values, res.values + Vector::dim, negate<double>());
//    cerr << "negate " << *this << " to " << res << endl;
//    exit(0);
    return res;
}




inline
bool
Vector::operator==(const Vector& v) const
{
    // Stroustrup 18.5.4
    return equal(values, values + dim, v.values);
}




inline
double
Vector::scalarProduct(const Vector& v) const
{
    // Stroustrup 22.6.2
    return inner_product(values, values + dim, v.values, 0.0);
}



//!!! not for Vector<int>

inline
double
Vector::min() const
{
    double elem = fabs(values[X]);
    for (int i = Y; i < Vector::dim; ++i)
    {
        elem = std::min(elem, fabs(values[i]));
    }
    return elem;
}



//!!! not for Vector<int>

inline
double
Vector::max() const
{
    double elem = fabs(values[X]);
    for (int i = Y; i < Vector::dim; ++i)
    {
        elem = std::max(elem, fabs(values[i]));
    }
    return elem;
}




inline
double
Vector::distance(const Vector& vec) const
{
    Vector dist(*this - vec);
    return sqrt(dist.scalarProduct(dist));
}




inline
double
Vector::magnitude() const
{
    // Stroustrup 22.6.2
    return sqrt(inner_product(values, values + dim, values, 0.0));
}

inline
Vector
operator+(const Vector& v, const Vector& w)
{
    Vector res(w);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, plus<double>());
    return res;
}




inline
Vector
operator+(double d, const Vector& v)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, plus<double>());
    return res;
}




inline
Vector
operator+(const Vector& v, double d)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, plus<double>());
    return res;
}




inline
Vector
operator-(const Vector& v, const Vector& w)
{
    Vector res(w);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, minus<double>());
    return res;
}




inline
Vector
operator-(double d, const Vector& v)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(res.values, res.values + Vector::dim, v.values, res.values, minus<double>());
    return res;
}




inline
Vector
operator-(const Vector& v, double d)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, minus<double>());
    return res;
}




inline
Vector
operator*(double d, const Vector& v)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, multiplies<double>());
    return res;
}




inline
Vector
operator*(const Vector& v, double d)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, multiplies<double>());
    return res;
}




inline
Vector
operator/(const Vector& v, double d)
{
    Vector res(d);
    // Stroustrup 18.4.3 and 18.6.2
    transform(v.values, v.values + Vector::dim, res.values, res.values, divides<double>());
    return res;
}




inline
double
magnitude(const Vector& v)
{
    return v.magnitude();
}

#undef X
#undef Y
#undef Z

#endif
// VECTOR_HPP

