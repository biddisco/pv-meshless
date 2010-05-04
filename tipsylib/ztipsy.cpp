/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
/* Doug Potter - 06/06/06 */
#include <math.h>
#include "ztipsy.hpp"

/* ********** ozTipsy ********** */

ozTipsy::ozTipsy( bTipsy::streambuf_type *sb )
#ifdef USE_ZLIB
    : zTipsyHelper(sb), obTipsy(&buf), oxdrstream(&buf)
#else
    : zTipsyHelper(sb), obTipsy(sb), oxdrstream(sb)
#endif
{
#ifdef USE_ZLIB
    if ( buf.DeflateInit() != Z_OK )
        this->setstate(ios_base::failbit);
#endif
    int i;

    // This is great for cosmological simulations.  This allows a range of values
    // for position of (-0.5,0.5).
    for( i=1; i<EXP_NUM; i++ ) m_exp[i] = 0;
    m_exp[INDX_POSITION] = -1;
}

ozTipsy::~ozTipsy()
{
}

/* ********** izTipsy ********** */

izTipsy::izTipsy( bTipsy::streambuf_type *sb )
#ifdef USE_ZLIB
    : zTipsyHelper(sb), ibTipsy(&buf), ixdrstream(&buf)
#else
    : zTipsyHelper(sb), ibTipsy(sb), ixdrstream(sb)
#endif
{
#ifdef USE_ZLIB
    if ( buf.InflateInit() != Z_OK )
        this->setstate(ios_base::failbit);
#endif
    int i;
    for( i=1; i<EXP_NUM; i++ ) m_exp[i] = 0;
    // This is great for cosmological simulations.  This allows a range of values
    // for position of (-0.5,0.5).
    m_exp[INDX_POSITION] = -1;
}

izTipsy::~izTipsy()
{
}

/* ********** ozTipsy ********** */

ozTipsy& ozTipsy::operator<<(TipsyHeader &val)
{
    oxdrstream *s = this;
    int pad = 0;
    header = val;

    if ( sputn( "tzip", 4 ) != 4 ) {
        this->setstate(ios_base::badbit);
        return *this;
    }
    *s << m_flags;

    *s << val.h_time << val.h_nBodies << val.h_nDims
       << val.h_nSph << val.h_nDark << val.h_nStar << pad;
    return *this;
}

zTipsyHelper::zTipsyHelper( bTipsy::streambuf_type *sb )
#ifdef USE_ZLIB
    : buf(sb)
#endif
{
}


#define POS_BITS 21
#define POS_MASK 0x001fffff
#define POS_REST 0xffe00000
#define POS_SIGN 0x00100000

#define DEN_BITS 17
#define DEN_MASK 0x0001ffff

#define MASK_WHERE 21 // 4 bits
#define EXP_WHERE  25 // 7 bits (0 -> 127)

// Here is the magic.  If our current particle limits are incorrect, then we
// adjust them and write out a sentinal node.  A sentinal node is "minus zero"
// in the first 21 bits, followed by the index of the exponent type, followed
// by the actual exponent.
void ozTipsy::AdjustLimit( float v, int index, bool bLog )
{
    oxdrstream *s = this;
    int e;
    long sentinal;
    // Zero can be stored with any exponent
    if ( bLog && v < 0.0 ) v = 0.0;
    if ( /*v != NAN &&*/ v != 0.0 && (m_flags&(1<<index))) { 
        if ( bLog ) v = log10f( v + 1.0 );
        frexpf( v, &e );
        if ( e > m_exp[index] ) {
            m_exp[index] = e;
            sentinal = (m_exp[index]<<EXP_WHERE) | (index<<MASK_WHERE) | POS_SIGN;
            *s << sentinal;
        }
    }
}


// Here is the magic.  If our current particle limits are incorrect, then we
// adjust them and write out a sentinal node.  A sentinal node is "minus zero"
// in the first 21 bits, and the exponent following.
void ozTipsy::AdjustLimits( TipsyBaseParticle &b )
{
    //oxdrstream *s = this;
    //bool bChanged = false;
    //int i, e, max_e;
    int i;
    float v;

    /* All three positions have a single 'max' exponent */
    v = fabs(b.pos[0]);
    for( i=1; i<header.h_nDims; i++ )
        if ( fabs(b.pos[i]) > v ) v = fabs(b.pos[i]);
    AdjustLimit( v, INDX_POSITION );

    AdjustLimit( b.density, INDX_DENSITY, true ); //JDP
    AdjustLimit( b.mass, INDX_MASS );
    AdjustLimit( b.phi,  INDX_PHI );
}

// The goal here is to extract 'n' bits of mantissa, preserving absolute error,
// NOT preserving a certain number of bits of accuracy.  The highest order bit
// will be the sign bit.
long zTipsyHelper::ExtractMantissa( float v, int bits, int exp )
{
    unsigned long mask = (unsigned long)-1 >> (32-bits);
    return (long)(floorf((v+0.5f)*(1<<(bits-exp-1)))) & mask;
}

float zTipsyHelper::BuildFloat( long v, int bits, int exp )
{
    long sign = 1 << (bits-1);
    long rest = (unsigned long)-1 << bits;
    if ( v & sign ) v |= rest;
    return (float)v / (1<<(bits-exp-1));
}

void ozTipsy::WriteFloat16( float v, int index, bool bLog )
{
    char B[2];
    long L;

    // Only save this if it has been requested
    if ( m_flags & (1<<index) ) {
        if ( bLog ) {
            if ( v < 0.0 ) v = 0.0;
            else v = log10f(v + 1.0);
        }
        L = ExtractMantissa(v,16,m_exp[index]) & 0xffff;
        B[0] = L >> 8;
        B[1] = L & 0xff;
        if ( sputn( B, sizeof(B) ) != sizeof(B) )
            this->setstate(ios_base::badbit);
    }
}

void ozTipsy::WriteParticle( TipsyBaseParticle &b )
{
    oxdrstream *s = this;

    //  0         10          21
    // +--------------------------------+
    // |       pos0          |   pos1   |
    // +--------------------------------+
    // +    pos1  |        pos3        ||
    // +--------------------------------+
    // +--------------------------------+

    unsigned long L[3];
    unsigned long P[2];
    //short S;
    //char B[2];

    L[0] = ExtractMantissa(b.pos[0],POS_BITS,m_exp[INDX_POSITION]);
    L[1] = ExtractMantissa(b.pos[1],POS_BITS,m_exp[INDX_POSITION]);
    L[2] = ExtractMantissa(b.pos[2],POS_BITS,m_exp[INDX_POSITION]);
    P[0] = L[0] | (L[1] << 21);
    P[1] = (L[1] >> 11) | (L[2]<<10);
    *s << P[0] << P[1];

    WriteFloat16( b.density, INDX_DENSITY, true ); // JDP
    WriteFloat16( b.mass, INDX_MASS, false );
    WriteFloat16( b.phi, INDX_PHI, false );

}

ozTipsy& ozTipsy::operator<<(TipsyGasParticle &val)
{
    AdjustLimits(val);
    AdjustLimit( val.rho,    INDX_RHO );
    AdjustLimit( val.temp,   INDX_TEMP );
    AdjustLimit( val.hsmooth,INDX_HSMOOTH );
    AdjustLimit( val.metals, INDX_METALS );
    WriteParticle(val);
    WriteFloat16( val.rho,    INDX_RHO,    false );
    WriteFloat16( val.temp,   INDX_TEMP,   false );
    WriteFloat16( val.hsmooth,INDX_HSMOOTH,false );
    WriteFloat16( val.metals, INDX_METALS, false );
    return *this;
}

ozTipsy& ozTipsy::operator<<(TipsyStarParticle &val)
{
    AdjustLimits(val);
    AdjustLimit( val.eps,    INDX_EPS );
    AdjustLimit( val.metals, INDX_METALS );
    AdjustLimit( val.tform,  INDX_TFORM );
    WriteParticle(val);
    WriteFloat16( val.eps,    INDX_EPS,    false );
    WriteFloat16( val.metals, INDX_METALS, false );
    WriteFloat16( val.tform,  INDX_TFORM,  false );
    return *this;
}

ozTipsy& ozTipsy::operator<<(TipsyDarkParticle &val)
{
    AdjustLimits(val);
    AdjustLimit( val.eps,    INDX_EPS );
    WriteParticle(val);
    WriteFloat16( val.eps,   INDX_EPS,     false );
    return *this;
}

/* ********** izTipsy ********** */

float izTipsy::ReadFloat16( int index, bool bLog )
{
    long L;
    char B[2];
    float v;

    if ( m_flags & (1<<index) ) {
        if ( sgetn(B,sizeof(B)) != sizeof(B) ) {
            this->setstate(ios_base::badbit);
            return 0.0;
        }
        L = (B[0]<<8) | (B[1]&0xff);
        v = BuildFloat(L,16,m_exp[index]);
        if ( bLog ) v = expf(v)/expf(10.0) - 1.0;
        return v;
    }
    return 0.0;
}


void izTipsy::ReadParticle(TipsyBaseParticle &val)
{
    ixdrstream *s = this;
    long L[3];
    unsigned long P[2];

    // Handle sentinal nodes
    *s >> P[0];
    while( (P[0]&POS_MASK) == POS_SIGN ) {
        int index = (P[0] >> MASK_WHERE) & 0xf;
        int exp = (P[0]>>EXP_WHERE);
        m_exp[index] = exp;
        *s >> P[0];
    }

    *s >> P[1];

    L[0] = P[0] & POS_MASK;
    L[1] = (((P[0] >> 21) & 0x7ff) | ((P[1]&0x3ff) << 11))  & POS_MASK;
    L[2] = (P[1] >> 10) & POS_MASK;

    val.pos[0] = BuildFloat(L[0],POS_BITS,m_exp[INDX_POSITION]);
    val.pos[1] = BuildFloat(L[1],POS_BITS,m_exp[INDX_POSITION]);
    val.pos[2] = BuildFloat(L[2],POS_BITS,m_exp[INDX_POSITION]);

//     printf( "%08lx %08lx -> %08lx %08lx %08lx <-> % f % f % f\n",
//             P[0], P[1], L[0], L[1], L[2],
//             val.pos[0], val.pos[1], val.pos[2] );
    val.density = ReadFloat16( INDX_DENSITY, true ); //JDP

    val.mass = ReadFloat16( INDX_MASS, false );
    val.phi  = ReadFloat16( INDX_PHI, false );
}

izTipsy& izTipsy::operator>>(TipsyHeader &val)
{
    ixdrstream *s = this;
    int pad = 0;
    char tzip[4];

    if ( sgetn(tzip,sizeof(tzip)) != sizeof(tzip)
         || memcmp(tzip,"tzip",sizeof(tzip)) != 0 ) {
        this->setstate(ios_base::badbit);
        return *this;
    }
    *s >> m_flags;

    *s >> val.h_time >> val.h_nBodies >> val.h_nDims
       >> val.h_nSph >> val.h_nDark >> val.h_nStar >> pad;
    header = val;
    return *this;
}

izTipsy& izTipsy::operator>>(TipsyGasParticle &val)
{
    //val.rho = val.temp = val.hsmooth = val.metals = 0.0;

    ReadParticle(val);
    val.rho     = ReadFloat16( INDX_RHO,    false );
    val.temp    = ReadFloat16( INDX_TEMP,   false );
    val.hsmooth = ReadFloat16( INDX_HSMOOTH,false );
    val.metals  = ReadFloat16( INDX_METALS, false );

    return *this;
}

izTipsy& izTipsy::operator>>(TipsyStarParticle &val)
{
    //val.eps = val.metals = val.tform = 0.0;
    ReadParticle(val);

    val.eps    = ReadFloat16( INDX_EPS,    false );
    val.metals = ReadFloat16( INDX_METALS, false );
    val.tform  = ReadFloat16( INDX_TFORM,  false );

    return *this;
}

izTipsy& izTipsy::operator>>(TipsyDarkParticle &val)
{
    //val.eps = 0.0;
    ReadParticle(val);

    val.eps = ReadFloat16( INDX_EPS,     false );

    return *this;
}

/* ********** ozfTipsy ********** */

ozfTipsy::ozfTipsy()
    : ozTipsy(this)
{
}

ozfTipsy::ozfTipsy( const char *filename, openmode mode )
    : ozTipsy(this)
{
    open(filename,mode);
}

/* ********** izfTipsy ********** */

izfTipsy::izfTipsy()
    : izTipsy(this)
{
}

izfTipsy::izfTipsy( const char *filename, openmode mode )
    : izTipsy(this)
{
    open(filename,mode);
}

