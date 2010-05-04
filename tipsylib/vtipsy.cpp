/** @file
 *  @brief In memory vector of Tipsy particles
 *  @author Doug Potter
 */

#include "vtipsy.hpp"

//! @brief Remove all elements from the vector.
void vTipsy::clear(void)
{
    gas.clear();
    dark.clear();
    star.clear();
}

//! @brief Read an entire Tipsy file into memory.
//! @param fin The iTipsy or ifTipsy object from which to read.
void vTipsy::operator<<( iTipsy &fin )
{
    TipsyHeader h; // The header structure
    gas_type::iterator  gi;
    dark_type::iterator di;
    star_type::iterator si;

    //! Read in the header.
    fin >> h;

    //! Save the expansion factor (simulation time).
    expFactor = h.h_time;

    //! Resize the arrays to hold the number of particles of each type that
    //! were specified in the header.  Note that any particles that were
    //! already loaded into the array are removed.
    gas.resize(h.h_nSph);
    dark.resize(h.h_nDark);
    star.resize(h.h_nStar);

    for( gi=gas.begin();  gi!=gas.end();  gi++ ) { fin >> *gi; }
    for( di=dark.begin(); di!=dark.end(); di++ ) { fin >> *di; }
    for( si=star.begin(); si!=star.end(); si++ ) { fin >> *si; }
}

//! @brief Write all particles to a Tipsy file.
//! @param fout The iTipsy or ifTipsy object into which to write.
void vTipsy::operator>>( oTipsy &fout )
{
    TipsyHeader h; // The header structure
    gas_type::iterator  gi;
    dark_type::iterator di;
    star_type::iterator si;

    //! Write a header.
    h.h_time  = expFactor;
    h.h_nDims = 3;
    h.h_nSph  = gas.size();
    h.h_nDark = dark.size();
    h.h_nStar = star.size();
    h.h_nBodies = h.h_nSph + h.h_nDark + h.h_nStar;
    fout << h;

    //! Write every particle of every type.  If there are not particles of
    //! a given type, then nothing is written.
    for( gi=gas.begin();  gi!=gas.end();  gi++ ) { fout << *gi; }
    for( di=dark.begin(); di!=dark.end(); di++ ) { fout << *di; }
    for( si=star.begin(); si!=star.end(); si++ ) { fout << *si; }
}
