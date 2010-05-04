/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
/* Doug Potter - 06/06/06 */
#include <cstring>
#include <assert.h>
#include "zstream.hpp"
#include <iostream>

zstreambuf::zstreambuf( std::streambuf *sb )
{ 
    theSB = sb;

    m_state = STATE_NONE;

    m_z.zalloc  = Z_NULL;
    m_z.zfree   = Z_NULL;
    m_z.opaque  = Z_NULL;
    m_z.avail_in= 0;
    m_z.next_in = Z_NULL;
    m_z.avail_out = 0;
    m_z.next_out  = Z_NULL;

}

zstreambuf::~zstreambuf()
{
    switch( m_state ) {
    case STATE_DEFLATING: DeflateEnd(); break;
    case STATE_INFLATING: InflateEnd(); break;
    default: break;
    }
}

int zstreambuf::DeflateInit(int level)
{
    assert( m_state == STATE_NONE );
    m_state = STATE_DEFLATING;
    return deflateInit2(&m_z,level,Z_DEFLATED,31,8,Z_DEFAULT_STRATEGY);
}

int zstreambuf::DeflateEnd()
{
    consume(Z_FINISH);
    theSB->pubsync();
    m_state = STATE_NONE;
    return deflateEnd(&m_z);
}

int zstreambuf::InflateInit()
{
    int rc;
    assert( m_state == STATE_NONE );
    m_state = STATE_INFLATING;
    m_z.avail_in = 0;
    m_z.next_in = Z_NULL;
    m_z.avail_out = sizeof(m_out);
    m_z.next_out = Z_NULL;
    rc = inflateInit2(&m_z,31);
    m_have = 0;
    assert( rc == Z_OK );
    return rc;
}

int zstreambuf::InflateEnd()
{
    m_state = STATE_NONE;
    return inflateEnd(&m_z);
}

/* Consume characters in the output buffer */
int zstreambuf::consume(int flush)
{
    int rc;
    m_z.next_in = (Bytef *)m_inp;
    do {
        m_z.avail_out = sizeof(m_out);
        m_z.next_out  = (Bytef *)m_out;
        rc = deflate(&m_z,flush);
        if ( rc == Z_STREAM_ERROR )
            return -1;
        assert(rc != Z_STREAM_ERROR);
        theSB->sputn(m_out, sizeof(m_out)-m_z.avail_out);
    } while (m_z.avail_out == 0);
    assert(m_z.avail_in == 0);     /* all input will be used */
    return 0;
}

int zstreambuf::sync()
{
    if ( consume(Z_SYNC_FLUSH) != 0 )
        return -1;
    return theSB->pubsync();
}

std::streamsize zstreambuf::xsgetn(char_type* __s, std::streamsize __n)
{
    int need = __n;
    int rc;

    do {
        if ( m_have >= need ) {
            memcpy( __s, m_ptr, need );
            m_ptr  += need;
            m_have -= need;
            return __n;
        }
        else if ( m_have ) {
            memcpy( __s, m_ptr, m_have );
            __s  += m_have;
            need -= m_have;
            m_have = 0;
        }
        if ( m_z.avail_in ) {
            m_z.next_out = (Bytef *)m_out;
            m_z.avail_out= sizeof(m_out);
            rc = inflate(&m_z,Z_NO_FLUSH);
            //assert( rc == Z_OK || rc==Z_STREAM_END );
            if ( rc != Z_OK && rc!=Z_STREAM_END)
                return 0;
            m_have = sizeof(m_out) - m_z.avail_out;
            m_ptr = m_out;
        }
        else {
            m_z.avail_in = theSB->sgetn(m_inp,sizeof(m_inp));
            m_z.next_in  = (Bytef *)m_inp;
            m_z.avail_out = sizeof(m_out);
            m_z.next_out  = (Bytef *)m_out;
        }
    } while( need );

    return __n;
}

std::streamsize zstreambuf::xsputn(const char_type* __s, std::streamsize __n)
{
    int avail, left = __n;

    avail = sizeof(m_inp) - m_z.avail_in;
    //cout << m_z.avail_in << " " << __n << " " << avail << endl;
    while( left > avail ) {
        memcpy( m_inp + m_z.avail_in, __s, avail );
        __s += avail;
        left -= avail;
        m_z.avail_in += avail;
        assert( m_z.avail_in == sizeof(m_inp) );
        if ( consume() != 0 )
            return 0;
        avail = sizeof(m_inp) - m_z.avail_in;
    }

    // Buffer remaining characters
    memcpy( m_inp + m_z.avail_in, __s, left );
    m_z.avail_in += left;

    return __n;
}

//zstreambuf::__streambuf_type *
//zstreambuf::setbuf(char_type* __s, std::streamsize __n)
//{
//    return this;
//}

zstreambuf::pos_type
zstreambuf::seekoff(off_type __off, std::ios_base::seekdir __way,
                    std::ios_base::openmode __mode)
{
    return 0;
}

zstreambuf::pos_type
zstreambuf::seekpos(pos_type __pos,
                    std::ios_base::openmode __mode)
{
    return 0;
}
