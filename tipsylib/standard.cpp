/**
 *  @file
 *  @brief Adapter for Tipsy standard format files
 *  @author Doug Potter
 */

#include <iostream> // for debugging
#ifdef _WIN32 // ntohl
 #include <winsock2.h>
#else
 #include <netinet/in.h> 
#endif
#include <assert.h>
#include "standard.h"

//! A representation of an integer in a tipsy file.
class disk_uint32_t
{
    uint32_t m_value; //!< Internal representation of an integer.
public:
    //! @brief Assign an integer.
    //! @param rhs The integer to assign.
    inline disk_uint32_t &operator=(const uint32_t &rhs) {
        m_value = ntohl(rhs);
        return *this;
    }
    //! @brief Cast to an integer.
    //! @return The integer value.
    inline operator uint32_t() const
        { return ntohl(m_value); }
};

//! A representation of a float in a tipsy file.
class disk_float
{
    uint32_t m_value; //!< Internal representation of a float.
public:
    //! @brief Assign a float.
    //! @param rhs The float to assign.
    inline disk_float &operator=(const float &rhs) {
	union { float rhs; uint32_t m_value; } x;
	x.rhs = rhs;
        m_value = ntohl(x.m_value);
        return *this;
    }
    //! @brief Cast to a float.
    //! @return The float value.
    inline operator float() const {
	union { float rhs; uint32_t m_value; } ret;
	ret.m_value = ntohl(m_value);
	return ret.rhs;
    }
};

//! A respresentation of a double in a tipsy file.
class disk_double
{
    uint32_t m_value[2]; //!< Internal representation of a double.
public:
    //! @brief Assign a double.
    //! @param rhs The double to assign.
    inline disk_double &operator=(const double &rhs) {
	union { double rhs; uint32_t m_value[2]; } x;
	x.rhs = rhs;
        m_value[1] = ntohl(x.m_value[0]);
        m_value[0] = ntohl(x.m_value[1]);
        return *this;
    }
    //! @brief Cast to a double.
    //! @return The double value.
    inline operator double() const {
	// This union is important.  Just doing a cast can result in the
	// compiler optimizing away the operation.
	union { double rhs; uint32_t m_value[2]; } ret;
	ret.m_value[0] = ntohl(m_value[1]);
	ret.m_value[1] = ntohl(m_value[0]);
	return ret.rhs;
    }
};

#include "tipsyrec.h"

TipsyStandardAdapter::TipsyStandardAdapter( streambuf_type *sb )
    : TipsyNativeAdapter(sb)
{
}

void TipsyStandardAdapter::getNext(void)
{
    tipsyrec_header hdr;
    tipsyrec_gas    gas;
    tipsyrec_dark   dark;
    tipsyrec_star   star;

    //! Depending on the section, read the appropriate particle type.
    switch( m_position.section() ) {
    case tipsypos::header:
	assert( m_sb->sgetn( (char *)(&hdr), sizeof(hdr) ) == sizeof(hdr) );
	m_dTime  = hdr.h_dTime;
	m_nBodies= hdr.h_nBodies;
	m_nDims  = hdr.h_nDims;
	m_nSph   = hdr.h_nSph;
	m_nDark  = hdr.h_nDark;
	m_nStar  = hdr.h_nStar;
	break;
    case tipsypos::gas:
	assert( m_sb->sgetn( (char *)(&gas), sizeof(gas) ) == sizeof(gas) );
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


void TipsyStandardAdapter::putNext(void)
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
