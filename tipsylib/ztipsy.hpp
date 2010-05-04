/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
/* Doug Potter - 06/06/06 */
#ifndef ZTIPSY_H
#define ZTIPSY_H

#include "ftipsy.hpp"
#ifdef USE_ZLIB
#include "zstream.hpp"
#endif

/* ********** zTipsyHelper - Common to compression classes ********** */


class zTipsyHelper
{
protected:
#ifdef USE_ZLIB
    zstreambuf buf;
#endif
    int m_exp[EXP_NUM];

protected:
    static long ExtractMantissa( float v, int bits, int exp );
    static float BuildFloat( long v, int bits, int exp );

    zTipsyHelper( bTipsy::streambuf_type *sb );

};

/* ********** ozTipsy - Output in COMPRESSED format ********** */

class ozTipsy : protected zTipsyHelper, public obTipsy, public oxdrstream
{
protected:
    // Returns true if all particles are "in the middle"
    void AdjustLimits( TipsyBaseParticle &b );
    void AdjustLimit( float v, int index, bool bLog=false );

    void WriteParticle( TipsyBaseParticle &b );
    void WriteFloat16( float v, int index, bool bLog );

public:
    ozTipsy( bTipsy::streambuf_type *sb );
    virtual ~ozTipsy();

    void Store( long mask ) { m_flags |= mask; }

    virtual ozTipsy& operator<<(TipsyHeader &val);
    virtual ozTipsy& operator<<(TipsyGasParticle &val);
    virtual ozTipsy& operator<<(TipsyStarParticle &val);
    virtual ozTipsy& operator<<(TipsyDarkParticle &val);
};

class ozfTipsy : protected std::filebuf, // Must be before ozTipsy (uses it)
                 public ozTipsy
{
public:
    ozfTipsy();
    ozfTipsy( const char *filename, openmode mode = out );

    bool is_open()
    { return std::filebuf::is_open(); }

    void open(const char* __s, 
              ios_base::openmode __mode = ios_base::out | ios_base::trunc)
    {
        if (!std::filebuf::open(__s, __mode | ios_base::out))
            this->setstate(ios_base::failbit);
    }

    void
    close()
    {
        if (!std::filebuf::close())
            this->setstate(ios_base::failbit);
    }

};
/* ********** izTipsy - Output in COMPRESSED format ********** */

class izTipsy : protected zTipsyHelper, public ibTipsy, public ixdrstream
{
protected:
    void ReadParticle( TipsyBaseParticle &b );
    float ReadFloat16( int index, bool bLog );

public:
    izTipsy( bTipsy::streambuf_type *sb );
    virtual ~izTipsy();

    virtual izTipsy& operator>>(TipsyHeader &val);
    virtual izTipsy& operator>>(TipsyGasParticle &val);
    virtual izTipsy& operator>>(TipsyStarParticle &val);
    virtual izTipsy& operator>>(TipsyDarkParticle &val);
};

class izfTipsy : protected std::filebuf, // Must be before izTipsy (uses it)
                 public izTipsy
{
public:
    izfTipsy();
    izfTipsy( const char *filename, openmode mode = in );

    bool is_open()
    { return std::filebuf::is_open(); }

    void open(const char* __s, 
              ios_base::openmode __mode = ios_base::in)
    {
        if (!std::filebuf::open(__s, __mode | ios_base::in))
            this->setstate(ios_base::failbit);
    }

    void
    close()
    {
        if (!std::filebuf::close())
            this->setstate(ios_base::failbit);
    }

};

#endif
