/**
 *  @file
 *  @brief Adapter for Tipsy native format files
 *  @author Doug Potter
 *
 *  @section Tipsy File Format
 *
 *  A Tipsy file consists of four different sections.  Except for the header,
 *  all of the sections are optional and are only present if there are a
 *  non-zero number of particles specified in the header.
 *
 *
\verbatim
  +----------------+
  | Header         |
  +----------------+
  | Gas Particles  |
  +----------------+
  | Dark Particles |
  +----------------+
  | Star Particles |
  +----------------+
 \endverbatim
 *
 *  @subsection hdr Header
\verbatim
                     4               8               12
     +---------------+---------------+---------------+---------------+
   0 | time (expansion factor)       | total count   | dimensions    |
     +---------------+---------------+---------------+---------------+
  16 | gas count     | dark count    | star count    | 0 (zero)      |
     +---------------+---------------+---------------+---------------+
                     20              24              28              32
\endverbatim
 *
 *  @subsection gasp Gas Particle
\verbatim
     0               4               8               12
     +---------------+---------------+---------------+---------------+
   0 | mass          | x             | y             | z             |
     +---------------+---------------+---------------+---------------+
  16 | Vx            | Vy            | Vz            | rho           |
     +---------------+---------------+---------------+---------------+
  32 | temp          | hsmooth       | metals        | phi           |
     +---------------+---------------+---------------+---------------+
\endverbatim
 *
 *  @subsection darkp Dark Particle
\verbatim
     0               4               8               12
     +---------------+---------------+---------------+---------------+
   0 | mass          | x             | y             | z             |
     +---------------+---------------+---------------+---------------+
  16 | Vx            | Vy            | Vz            | eps           |
     +---------------+---------------+---------------+---------------+
  32 | phi           |
     +---------------+
\endverbatim
 *
 *  @subsection starp Star Particle
\verbatim
     0               4               8               12
     +---------------+---------------+---------------+---------------+
   0 | mass          | x             | y             | z             |
     +---------------+---------------+---------------+---------------+
  16 | Vx            | Vy            | Vz            | metals        |
     +---------------+---------------+---------------+---------------+
  32 | tform         | eps           | phi           |
     +---------------+---------------+---------------+
\endverbatim
 *
 */

#include <iostream> // For debugging
#include <assert.h>
#include "native.h"

typedef uint32_t disk_uint32_t; //!< An integer as stored in the file.
typedef float    disk_float;    //!< A float as stored in the file.
typedef double   disk_double;   //!< A double as stored in the file.

#include "tipsyrec.h"

//! Construct a Tipsy Native Adapter.
TipsyNativeAdapter::TipsyNativeAdapter( streambuf_type *sb )
{
    //! Save the stream buffer and perform initialization of the base class.
    init(sb);

    //! Reserve enough space for the largest particle in the values vector.
    m_values.resize(idcountid,0.0f);

    //! Setup the field map (maps field name to the position in the vector)
    m_idmap["mass"]    = idmass;
    m_idmap["x"]       = idx;
    m_idmap["y"]       = idy;
    m_idmap["z"]       = idz;
    m_idmap["vx"]      = idvx;
    m_idmap["vy"]      = idvy;
    m_idmap["vz"]      = idvz;
    m_idmap["rho"]     = idrho;
    m_idmap["temp"]    = idtemp;    // Star(tform)
    m_idmap["tform"]   = idtform;   // shared
    m_idmap["eps"]     = ideps;     // Dark(eps), Gas(hsmooth)
    m_idmap["hsmooth"] = idhsmooth; // Star(metals)
    m_idmap["metals"]  = idmetals;  // Star(metals)
    m_idmap["phi"]     = idphi;     // Dark(phi)
    m_idmap["density"] = iddensity;
}

//! Initialize our object with the given stream buffer.
void TipsyNativeAdapter::init( streambuf_type *sb )
{
    //! All we have to do is save the stream buffer and indicate that the
    //! buffer points to the header record (normally the start of the file).
    m_sb = sb;
    m_position = tipsypos(tipsypos::header,0);
}

//! Advance the file position forward by one particle.
void TipsyNativeAdapter::forward(void)
{
    //! Advance to the next particle
    m_position.offset()++;

    //! If we have switched sections, advance to the next section
    switch( m_position.section() ) {
    case tipsypos::header:
	m_position = tipsypos(tipsypos::gas, 0);
    case tipsypos::gas:
	if ( m_position.offset() >= m_nSph ) 
	    m_position = tipsypos(tipsypos::dark, 0);
	else
	    break;
    case tipsypos::dark:
	if ( m_position.offset() >= m_nDark ) 
	    m_position = tipsypos(tipsypos::star, 0);
	else
	    break;
    case tipsypos::star:
	if ( m_position.offset() >= m_nStar ) 
	    m_position = tipsypos(tipsypos::eof, 0);
	break;
    default:
	assert( "Attempt to move past EOF" == 0 );
    }
}

//! Read the next particle (at the current position).
void TipsyNativeAdapter::getNext(void)
{
    tipsyrec_header hdr;
    tipsyrec_gas    gas;
    tipsyrec_dark   dark;
    tipsyrec_star   star;

    //! Depending on the section, read the appropriate particle type.
    switch( m_position.section() ) {
    case tipsypos::header:
	assert( m_sb->sgetn( (char *)(&hdr), sizeof(hdr) ) == sizeof(hdr));
	m_dTime  = hdr.h_dTime;
	m_nBodies= hdr.h_nBodies;
	m_nDims  = hdr.h_nDims;
	m_nSph   = hdr.h_nSph;
	m_nDark  = hdr.h_nDark;
	m_nStar  = hdr.h_nStar;
	break;
    case tipsypos::gas:
	assert( m_sb->sgetn( (char *)(&gas), sizeof(gas) ) == sizeof(gas));
	m_values[idmass]   = gas.mass;
	m_values[idx]      = gas.pos[0];
	m_values[idy]      = gas.pos[1];
	m_values[idz]      = gas.pos[2];
	m_values[idvx]     = gas.vel[0];
	m_values[idvy]     = gas.vel[1];
	m_values[idvz]     = gas.vel[2];
	m_values[idrho]    = gas.rho;
	m_values[idtemp]   = gas.temp;
	m_values[idhsmooth]= gas.hsmooth;
	m_values[idmetals] = gas.metals;
	m_values[idphi]    = gas.phi;
	break;
    case tipsypos::dark:
	assert( m_sb->sgetn( (char *)(&dark), sizeof(dark) ) == sizeof(dark) );
	m_values[idmass]   = dark.mass;
	m_values[idx]      = dark.pos[0];
	m_values[idy]      = dark.pos[1];
	m_values[idz]      = dark.pos[2];
	m_values[idvx]     = dark.vel[0];
	m_values[idvy]     = dark.vel[1];
	m_values[idvz]     = dark.vel[2];
	m_values[ideps]    = dark.eps;
	m_values[idphi]    = dark.phi;
	break;
    case tipsypos::star:
	assert( m_sb->sgetn( (char *)(&star), sizeof(star) ) == sizeof(star) );
	m_values[idmass]   = star.mass;
	m_values[idx]      = star.pos[0];
	m_values[idy]      = star.pos[1];
	m_values[idz]      = star.pos[2];
	m_values[idvx]     = star.vel[0];
	m_values[idvy]     = star.vel[1];
	m_values[idvz]     = star.vel[2];
	m_values[idmetals] = star.metals;
	m_values[idtform]  = star.tform;
	m_values[ideps]    = star.eps;
	m_values[idphi]    = star.phi;
	break;
    default:
	assert( "Invalid read position" == 0 );
    }

    //! Move forward one particle
    forward();
}

//! Write the current particle to the stream.
void TipsyNativeAdapter::putNext(void)
{
    tipsyrec_header hdr;
    tipsyrec_gas    gas;
    tipsyrec_dark   dark;
    tipsyrec_star   star;

    //! Depending on the section, write the appropriate particle type.
    switch( m_position.section() ) {
    case tipsypos::header:
	hdr.h_dTime  = m_dTime;
	hdr.h_nBodies= m_nBodies;
	hdr.h_nDims  = m_nDims;
	hdr.h_nSph   = m_nSph;
	hdr.h_nDark  = m_nDark;
	hdr.h_nStar  = m_nStar;
	hdr.h_MBZ    = 0;
	m_sb->sputn( (char *)(&hdr), sizeof(hdr) );
	break;
    case tipsypos::gas:
	gas.mass   = m_values[idmass];
	gas.pos[0] = m_values[idx];
	gas.pos[1] = m_values[idy];
	gas.pos[2] = m_values[idz];
	gas.vel[0] = m_values[idvx];
	gas.vel[1] = m_values[idvy];
	gas.vel[2] = m_values[idvz];
	gas.rho    = m_values[idrho];
	gas.temp   = m_values[idtemp];
	gas.hsmooth= m_values[idhsmooth];
	gas.metals = m_values[idmetals];
	gas.phi    = m_values[idphi];
	m_sb->sputn( (char *)(&gas), sizeof(gas) );
	break;
    case tipsypos::dark:
	dark.mass   = m_values[idmass];
	dark.pos[0] = m_values[idx];
	dark.pos[1] = m_values[idy];
	dark.pos[2] = m_values[idz];
	dark.vel[0] = m_values[idvx];
	dark.vel[1] = m_values[idvy];
	dark.vel[2] = m_values[idvz];
	dark.eps    = m_values[ideps];
	dark.phi    = m_values[idphi];
	m_sb->sputn( (char *)(&dark), sizeof(dark) );
	break;
    case tipsypos::star:
	star.mass   = m_values[idmass];
	star.pos[0] = m_values[idx];
	star.pos[1] = m_values[idy];
	star.pos[2] = m_values[idz];
	star.vel[0] = m_values[idvx];
	star.vel[1] = m_values[idvy];
	star.vel[2] = m_values[idvz];
	star.metals = m_values[idmetals];
	star.tform  = m_values[idtform];
	star.eps    = m_values[ideps];
	star.phi    = m_values[idphi];
	m_sb->sputn( (char *)(&star), sizeof(star) );
	break;
    default:
	assert( "Invalid read position" == 0 );
    }

    //! Move forward one particle
    forward();
}

//! Move the get pointer (seek) to the specified particle.
tipsypos TipsyNativeAdapter::seekg( tipsypos & pos )
{
    streambuf_type::pos_type where = 0;
    tipsypos::offset_type off;

    //! Find the base location for the section
    switch( pos.section() ) {
    case tipsypos::star:
	where += sizeof(tipsyrec_dark) * m_nDark;
    case tipsypos::dark:
	where += sizeof(tipsyrec_gas) * m_nSph;
    case tipsypos::gas:
    case tipsypos::particle:
	where += sizeof(tipsyrec_header);
    case tipsypos::header:
	break;
    default:
	assert( "Invalid seek section" == 0 );
	break;
    }

    //! Now add in the offset
    switch( pos.section() ) {
    case tipsypos::star:
	assert( pos.offset() < m_nStar );
	where += sizeof(tipsyrec_star) * pos.offset();
	break;
    case tipsypos::dark:
	assert( pos.offset() < m_nDark );
	where += sizeof(tipsyrec_dark) * pos.offset();
	break;
    case tipsypos::gas:
	assert( pos.offset() < m_nSph );
	where += sizeof(tipsyrec_gas) * pos.offset();
	break;
    case tipsypos::header:
	assert( pos.offset() == 0 );
	break;
    case tipsypos::particle:
	off = pos.offset();
	if ( off < m_nSph ) {
	    where += sizeof(tipsyrec_gas) * off;
	    pos = tipsypos(tipsypos::gas,off);
	}
	else {
	    off -= m_nSph;
	    where += sizeof(tipsyrec_gas) * m_nSph;
	    if ( off < m_nDark ) {
		where += sizeof(tipsyrec_dark) * off;
		pos = tipsypos(tipsypos::dark,off);
	    }
	    else {
		off -= m_nDark;
		where += sizeof(tipsyrec_dark) * m_nDark;

		assert( off < m_nStar );
		where += sizeof(tipsyrec_star) * off;
		pos = tipsypos(tipsypos::star,off);
	    }
	}
	break;
    default:
	assert( "Invalid seek section" == 0 );
	break;
    }

    //! Try to seek and update the position if successful
    if ( m_sb->pubseekpos(where) == where ) {
	m_position = pos;
    }
    return m_position;
}

//! Seek the put pointer.
tipsypos TipsyNativeAdapter::seekp( tipsypos & pos )
{
    //! We have only a single pointer, so we just seek the get pointer.
    return seekg( pos );
}
