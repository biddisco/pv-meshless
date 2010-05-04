/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
/* Doug Potter - 06/06/06 */
#ifndef ZSTREAM_H
#define ZSTREAM_H

#include <ios>
#include <streambuf>
#include <ostream>
#include <zlib.h>

/*
**  The idea here is to provide a stream buffer that "compresses" any characters
**  that are written, and decompresses any characters that are read.  This is
**  really just a wrapper for a real stream buffer.
*/
class zstreambuf : public std::streambuf
{
    z_stream m_z;

    enum STATE { STATE_NONE, STATE_DEFLATING, STATE_INFLATING } m_state;

    char     m_inp[32768];
    char     m_out[32768];
    int      m_have;
    char     *m_ptr;


    std::streambuf *theSB;
    int consume( int flush=Z_NO_FLUSH );

public:
    zstreambuf( std::streambuf *sb );
    virtual ~zstreambuf();

    int DeflateInit(int level = Z_DEFAULT_COMPRESSION);
    int DeflateEnd();

    int InflateInit();
    int InflateEnd();

protected:
    virtual int sync();
    //virtual std::streamsize xsgetn(char_type* __s, std::streamsize __n);
    virtual std::streamsize xsgetn(char * __s, std::streamsize __n);
    virtual std::streamsize xsputn(const char_type* __s, std::streamsize __n);

//    virtual __streambuf_type*
//    setbuf(char_type* __s, std::streamsize __n);

    virtual pos_type
    seekoff(off_type __off, std::ios_base::seekdir __way,
            std::ios_base::openmode __mode = std::ios_base::in|std::ios_base::out);

    virtual pos_type
    seekpos(pos_type __pos,
            std::ios_base::openmode __mode = std::ios_base::in|std::ios_base::out);

};

#endif
