/** @file
 *  @brief In memory vector of Tipsy particles
 *  @author Doug Potter
 */

/** @page vtipsy_page The vTipsy class
 *
 *  Most file operations can be done by reading particles one at a time,
 *  directly from the file (or writing directly to a file).  This is the
 *  best way to accomplish what you want because with large input files,
 *  the particles may not fit in memory.
 *
 *  There are, however, cases where it is necessary to keep the particles
 *  in memory.  The vTipsy class was created for just this purpose.  It
 *  allows you to read an entire Tipsy file into memory, or to create all
 *  of the particles in memory, then write them to a file.
 *
 *  Reading and writing to a file are optional, so it is possible to use
 *  this class as a simple particle container.
 *
 *  To get a feel for how this class works, take a look at the following
 *  test program.  It reads a Tipsy file into memory, then copies the
 *  particles, one at a time to another memory array.  Finally, the new
 *  array is written to disk.
 *
 *  The example actually creates two copies of the file in memory which is
 *  a bad thing normally.  The example was used to test the vTipsy class as
 *  the output file should be identical to the input file.
 *
 *  @include vtest.cpp
 *
 */



#ifndef VTIPSY_H
#define VTIPSY_H

#include <vector>
#include "ftipsy.hpp"

/** @brief A vector of Tipsy particles
 *  This class allows an entire Tipsy file to be read in memory, or constructed
 *  in memory and later written to disk.
 */
class vTipsy {
public:
    //! Array type of gas particles.
    typedef std::vector<TipsyGasParticle>  gas_type;
    //! Array type of dark particles.
    typedef std::vector<TipsyDarkParticle> dark_type;
    //! Array type of star particles.
    typedef std::vector<TipsyStarParticle> star_type;

    double    expFactor; //!< Expansion factor (simulation time)
    gas_type  gas;       //!< The gas particles
    dark_type dark;      //!< The dark particles
    star_type star;      //!< The star particles

public:
    //! @brief Construct a vTipsy particle array.
    //! The expansion factor (simulation time) is initialized to 0.0.
    //! If writing, and a different value is required, the expFactor member
    //! must be set appropriately.
    vTipsy() : expFactor(0.0) {}

    //! @brief Remove all elements from the vector.
    void clear(void);

    //! @brief Read an entire Tipsy file into memory.
    //! @param fin The iTipsy or ifTipsy object from which to read.
    void operator<<( iTipsy &fin );

    //! @brief Write all particles to a Tipsy file.
    //! @param fout The iTipsy or ifTipsy object into which to write.
    void operator>>( oTipsy &fout );

    //! @brief Add a gas particle to the array.
    //! @param val The gas particle to add.
    void operator<<( const TipsyGasParticle &val )
	{ gas.push_back(val); }

    //! @brief Add a dark particle to the array.
    //! @param val The dark particle to add.
    void operator<<( const TipsyDarkParticle &val )
	{ dark.push_back(val); }

    //! @brief Add a star particle to the array.
    //! @param val The star particle to add.
    void operator<<( const TipsyStarParticle &val )
	{ star.push_back(val); }
};

#endif
