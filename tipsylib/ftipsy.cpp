/**
 *  @file
 *  @brief Generic Tipsy file class
 *  @author Doug Potter
 *
 */

#include <cstring>
#include <assert.h>
#include <iostream>
#ifdef _WIN32 // ntohl
 #include <winsock2.h>
#else
 #include <netinet/in.h>
#endif
#include "ftipsy.hpp"
#include "native.h"
#include "standard.h"

//****************************************************************************
//  TipsyIOS
//****************************************************************************

void TipsyIOS::setAdapter( const char *type ) {
    if ( adapter ) {
	delete adapter;
	adapter = 0;
    }
    if ( strcmp(type,"standard") == 0 )
	adapter = new TipsyStandardAdapter(m_sb);
    else if ( strcmp(type,"native") == 0 )
	adapter = new TipsyNativeAdapter(m_sb);
    else {
	assert( 0==1 /* Unknown adapter type */ );
	// What we really want to do at this point is attempt to dynamically
	// load the correct adapter.  This is a future project.
    }

    // Save field ids for fast lookup
    idmass    = adapter->getFieldID("mass");
    idx       = adapter->getFieldID("x");
    idy       = adapter->getFieldID("y");
    idz       = adapter->getFieldID("z");
    idvx      = adapter->getFieldID("vx");
    idvy      = adapter->getFieldID("vy");
    idvz      = adapter->getFieldID("vz");
    idrho     = adapter->getFieldID("rho");
    idtemp    = adapter->getFieldID("temp");
    idtform   = adapter->getFieldID("tform");
    ideps     = adapter->getFieldID("eps");
    idhsmooth = adapter->getFieldID("hsmooth");
    idmetals  = adapter->getFieldID("metals");
    idphi     = adapter->getFieldID("phi");
    iddensity = adapter->getFieldID("density");
}


void TipsyIOS::init( streambuf_type *sb, const char * type ) {
    m_sb = sb;
    if ( type ) setAdapter(type);
}

TipsyIOS::~TipsyIOS() {
    if ( adapter ) delete adapter;
}

//****************************************************************************
//  iTipsy
//****************************************************************************

inline void iTipsy::read(TipsyHeader &val)
{
    val.h_time    = adapter->m_dTime;
    val.h_nBodies = adapter->m_nBodies;
    val.h_nDims   = adapter->m_nDims;
    val.h_nSph    = adapter->m_nSph;
    val.h_nDark   = adapter->m_nDark;
    val.h_nStar   = adapter->m_nStar;
}

inline void iTipsy::read(TipsyGasParticle &val)
{
    val.mass   = adapter->Field(idmass);
    val.pos[0] = adapter->Field(idx);
    val.pos[1] = adapter->Field(idy);
    val.pos[2] = adapter->Field(idz);
    val.vel[0] = adapter->Field(idvx);
    val.vel[1] = adapter->Field(idvy);
    val.vel[2] = adapter->Field(idvz);
    val.phi    = adapter->Field(idphi);
    val.density= adapter->Field(iddensity);
    val.rho    = adapter->Field(idrho);
    val.temp   = adapter->Field(idtemp);
    val.hsmooth= adapter->Field(idhsmooth);
    val.metals = adapter->Field(idmetals);
}

inline void iTipsy::read(TipsyDarkParticle &val)
{
    val.mass   = adapter->Field(idmass);
    val.pos[0] = adapter->Field(idx);
    val.pos[1] = adapter->Field(idy);
    val.pos[2] = adapter->Field(idz);
    val.vel[0] = adapter->Field(idvx);
    val.vel[1] = adapter->Field(idvy);
    val.vel[2] = adapter->Field(idvz);
    val.phi    = adapter->Field(idphi);
    val.density= adapter->Field(iddensity);
    val.eps    = adapter->Field(ideps);
}

inline void iTipsy::read(TipsyStarParticle &val)
{
    val.mass   = adapter->Field(idmass);
    val.pos[0] = adapter->Field(idx);
    val.pos[1] = adapter->Field(idy);
    val.pos[2] = adapter->Field(idz);
    val.vel[0] = adapter->Field(idvx);
    val.vel[1] = adapter->Field(idvy);
    val.vel[2] = adapter->Field(idvz);
    val.phi    = adapter->Field(idphi);
    val.density= adapter->Field(iddensity);
    val.eps    = adapter->Field(ideps);
    val.metals = adapter->Field(idmetals);
    val.tform  = adapter->Field(idtform);
}


iTipsy& iTipsy::operator>>(TipsyHeader &val)
{
    uint32_t nDensity;
    assert( adapter->tellg() == tipsypos::header );
    adapter->getNext();
    read(val);
    if ( m_fDensity ) {
	*m_fDensity >> nDensity;
    	assert( nDensity == val.h_nBodies );
    }
    return *this;
}

iTipsy& iTipsy::operator>>(TipsyGasParticle &val)
{
    assert( adapter->tellg() == tipsypos::gas );
    adapter->getNext();
    read(val);
    if ( m_fDensity ) *m_fDensity >> val.density;
    return *this;
}

iTipsy& iTipsy::operator>>(TipsyDarkParticle &val)
{
    assert( adapter->tellg() == tipsypos::dark );
    adapter->getNext();
    read(val);
    if ( m_fDensity ) *m_fDensity >> val.density;
    return *this;
}

iTipsy& iTipsy::operator>>(TipsyStarParticle &val)
{
    assert( adapter->tellg() == tipsypos::star );
    adapter->getNext();
    read(val);
    if ( m_fDensity ) *m_fDensity >> val.density;
    return *this;
}

iTipsy &iTipsy::seekg( tipsypos pos )
{
    // It is possible to rewind the density file
    if ( pos == tipsypos::header ) {
	if ( m_fDensity != 0 ) m_fDensity->seekg(0);
    }
    else {
	assert( m_fDensity == 0 );
    }
    adapter->seekg(pos);
    return *this;
}

tipsypos iTipsy::tellg()
{
    return adapter->tellg();
}

//****************************************************************************
//  oTipsy
//****************************************************************************

inline void oTipsy::stow(TipsyHeader &val)
{
    adapter->m_dTime   = val.h_time;
    adapter->m_nBodies = val.h_nBodies;
    adapter->m_nDims   = val.h_nDims;
    adapter->m_nSph    = val.h_nSph;
    adapter->m_nDark   = val.h_nDark;
    adapter->m_nStar   = val.h_nStar;
}

inline void oTipsy::stow(TipsyGasParticle &val)
{
    adapter->Field(idmass)   = val.mass;
    adapter->Field(idx)      = val.pos[0];
    adapter->Field(idy)      = val.pos[1];
    adapter->Field(idz)      = val.pos[2];
    adapter->Field(idvx)     = val.vel[0];
    adapter->Field(idvy)     = val.vel[1];
    adapter->Field(idvz)     = val.vel[2];
    adapter->Field(idphi)    = val.phi;
    adapter->Field(iddensity)= val.density;
    adapter->Field(idrho)    = val.rho;
    adapter->Field(idtemp)   = val.temp;
    adapter->Field(idhsmooth)= val.hsmooth;
    adapter->Field(idmetals) = val.metals;
}

inline void oTipsy::stow(TipsyDarkParticle &val)
{
    adapter->Field(idmass)   = val.mass;
    adapter->Field(idx)      = val.pos[0];
    adapter->Field(idy)      = val.pos[1];
    adapter->Field(idz)      = val.pos[2];
    adapter->Field(idvx)     = val.vel[0];
    adapter->Field(idvy)     = val.vel[1];
    adapter->Field(idvz)     = val.vel[2];
    adapter->Field(idphi)    = val.phi;
    adapter->Field(iddensity)= val.density;
    adapter->Field(ideps)    = val.eps;
}

inline void oTipsy::stow(TipsyStarParticle &val)
{
    adapter->Field(idmass)   = val.mass;
    adapter->Field(idx)      = val.pos[0];
    adapter->Field(idy)      = val.pos[1];
    adapter->Field(idz)      = val.pos[2];
    adapter->Field(idvx)     = val.vel[0];
    adapter->Field(idvy)     = val.vel[1];
    adapter->Field(idvz)     = val.vel[2];
    adapter->Field(idphi)    = val.phi;
    adapter->Field(iddensity)= val.density;
    adapter->Field(ideps)    = val.eps;
    adapter->Field(idmetals) = val.metals;
    adapter->Field(idtform)  = val.tform;
}

oTipsy& oTipsy::operator<<(TipsyHeader &val)
{
    assert( adapter->tellp() == tipsypos::header );
    stow(val);
    adapter->putNext();
    if ( m_fDensity ) *m_fDensity << val.h_nBodies;
    return *this;
}

oTipsy& oTipsy::operator<<(TipsyGasParticle &val)
{
    assert( adapter->tellp() == tipsypos::gas );
    stow(val);
    adapter->putNext();
    if ( m_fDensity ) *m_fDensity << val.density;
    return *this;
}

oTipsy& oTipsy::operator<<(TipsyDarkParticle &val)
{
    assert( adapter->tellp() == tipsypos::dark );
    stow(val);
    adapter->putNext();
    if ( m_fDensity ) *m_fDensity << val.density;
    return *this;
}

oTipsy& oTipsy::operator<<(TipsyStarParticle &val)
{
    assert( adapter->tellp() == tipsypos::star );
    stow(val);
    adapter->putNext();
    if ( m_fDensity ) *m_fDensity << val.density;
    return *this;
}

oTipsy &oTipsy::seekp( tipsypos pos )
{
    // It is possible to rewind the density file
    if ( pos == tipsypos::header ) {
	if ( m_fDensity != 0 ) m_fDensity->seekp(0);
    }
    else {
	assert( m_fDensity == 0 );
    }
    adapter->seekp(pos);
    return *this;
}

tipsypos oTipsy::tellp()
{
    return adapter->tellp();
}

//****************************************************************************
//  ifTipsy
//****************************************************************************

void ifTipsy::open(const char* iname, const char *type) {
    M_filebuf.open(iname, std::ios_base::in|std::ios_base::binary);
    setAdapter(type);
}

//****************************************************************************
//  ofTipsy
//****************************************************************************

void ofTipsy::open(const char* iname, const char *type) {
    M_filebuf.open(iname, std::ios_base::out|std::ios_base::binary);
    setAdapter(type);
}
