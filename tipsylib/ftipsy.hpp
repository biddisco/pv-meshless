/**
 *  @file
 *  @brief Generic Tipsy file class
 *  @author Doug Potter
 */
/** @mainpage Tipsy File Library (tipsylib)
 *
 *  This library provides a set of classes for reading and writing Tipsy
 *  native and standard files.
 *
 *  If you need in-memory Tipsy files, see \ref vtipsy_page.
 *
 *  @section ex1 A simple example
 *
 *  The easiest way to see how these classes work is to look at an example
 *  program.  The following example will convert a tipsy file in "native"
 *  format into a tipsy file in "standard" format.
 *
 *  @include tostd.cpp
 */
#ifndef FTIPSY_H
#define FTIPSY_H

//This is likely not required as the stream classes use 64 bit offsets.
//#define _LARGEFILE_SOURCE 1
//#define _FILE_OFFSET_BITS 64

#include <ios>
#include <map>
#include <vector>
#include <streambuf>
#include <fstream>
#include <istream>

#include "tipsypos.h"

class TipsyAdapter;

//! @brief Description of the contents of a Tipsy file
class TipsyHeader
{
public:
    double  h_time;     //!< Expansion factor
    uint32_t h_nBodies; //!< Total number of particles
    uint32_t h_nDims;   //!< Number of dimensions (3)
    uint32_t h_nSph;    //!< Number of gas particles
    uint32_t h_nDark;   //!< Number of dark particles
    uint32_t h_nStar;   //!< Number of star particles
};
//! @brief Information common to all particle types.
class TipsyBaseParticle
{
public:
    float mass;       //!< Mass of the particle
    float pos[3];     //!< Positions (x,y,z)
    float vel[3];     //!< Velocity (Vx,Vy,Vz)
    float phi;        //!< Potential
    float density;    //!< Density
    // Density is NOT part of the Tipsy file (from an external density file)
};

//! @brief Gas Particle.
class TipsyGasParticle : public TipsyBaseParticle
{
public:
    float rho;        //!< rho
    float temp;       //!< Temperature
    float hsmooth;    //!< hsmooth
    float metals;     //!< Metals.
};

//! @brief Dark Particle.
class TipsyDarkParticle : public TipsyBaseParticle {
public:
    float eps;        //!< Softening
};

//! @brief Star Particle.
class TipsyStarParticle : public TipsyBaseParticle {
public:
    float eps;        //!< Softening
    float metals;     //!< Metals
    float tform;      //!< tform
};

//! @brief Common base class for iTipsy and oTipsy
class TipsyIOS {
public:
    //! Stream input type
    typedef std::basic_streambuf<char> streambuf_type;

protected:
    streambuf_type *m_sb;   //!< The stream buffer object.
    TipsyAdapter *adapter;  //!< The adapter object.

    std::size_t idmass;     //!< Index of "mass" field.
    std::size_t idx;        //!< Index of "x" field.
    std::size_t idy;        //!< Index of "y" field.
    std::size_t idz;        //!< Index of "z" field.
    std::size_t idvx;       //!< Index of "vx" field.
    std::size_t idvy;       //!< Index of "vy" field.
    std::size_t idvz;       //!< Index of "vz" field.
    std::size_t idrho;      //!< Index of "rho" field.
    std::size_t idtemp;     //!< Index of "temp" field.
    std::size_t idtform;    //!< Index of "tform" field.
    std::size_t ideps;      //!< Index of "eps" field.
    std::size_t idhsmooth;  //!< Index of "hsmooth" field.
    std::size_t idmetals;   //!< Index of "metals" field.
    std::size_t idphi;      //!< Index of "phi" field.
    std::size_t iddensity;  //!< Index of "density" field.

protected:
    //! @brief Initialize our state.
    void init( streambuf_type *b, const char * adaptertype );

    //! @brief Default constructor.
    TipsyIOS() { adapter=0; m_sb=0; }

public:
    //! @brief Destructor.
    virtual ~TipsyIOS();

    //! @brief Return the current streambuf object.
    streambuf_type * rdbuf() const
	{ return m_sb; }

    //! @brief Set the adapter type.
    void setAdapter( const char *adaptertype );
};



//! @brief Read a Tipsy formatted stream.
class iTipsy : virtual public TipsyIOS {
protected:
    //! Read a header from the current stream.
    inline void read(TipsyHeader &val);

    //! Read a gas particle from the current stream.
    inline void read(TipsyGasParticle &val);

    //! Read a dark particle from the current stream.
    inline void read(TipsyDarkParticle &val);

    //! Read a star particle from the current stream.
    inline void read(TipsyStarParticle &val);

    //! The density input stream.
    std::istream *m_fDensity;

public:
    //! @brief Construct an iTipsy object.
    explicit iTipsy( std::basic_streambuf<char> *sb,
		     const char * adaptertype = 0 )
	: m_fDensity(0)
	{ this->init(sb,adaptertype); }

    //! @brief Destry an iTipsy object.
    virtual ~iTipsy() {}

    //! @brief Set the density input file.
    void DensityIn( std::istream &fDensity )
	{ m_fDensity = &fDensity; }

    //! @brief Read the Tipsy header.
    //! @param val A TipsyHeader structure.
    virtual iTipsy& operator>>(TipsyHeader &val);

    //! @brief Read a gas particle
    //! @param val A TipsyGasParticle Structure.
    virtual iTipsy& operator>>(TipsyGasParticle &val);

    //! @brief Read a star particle
    //! @param val A TipsyStarParticle Structure.
    virtual iTipsy& operator>>(TipsyStarParticle &val);

    //! @brief Read a dark particle
    //! @param val A TipsyDarkParticle Structure.
    virtual iTipsy& operator>>(TipsyDarkParticle &val);

    //! @brief Seek to a specific particle in the stream.
    //! @param pos The new position.
    virtual iTipsy &seekg( tipsypos pos );

    //! @brief Return the current particle position.
    //! @return The current particle position.
    virtual tipsypos tellg();
};

//! @brief Read a Tipsy formatted stream.
class oTipsy : virtual public TipsyIOS {
protected:
    //! Write a header to the current stream.
    inline void stow(TipsyHeader &val);

    //! Write a gas particle to the current stream.
    inline void stow(TipsyGasParticle &val);

    //! Write a dark particle to the current stream.
    inline void stow(TipsyDarkParticle &val);

    //! Write a star particle to the current stream.
    inline void stow(TipsyStarParticle &val);

    //! The density output stream.
    std::ostream *m_fDensity;

public:
    //! @brief Construct an oTipsy object.
    //! @param sb The output stream buffer.
    //! @param adaptertype The adapter type ("native" or "standard").
    explicit oTipsy( streambuf_type * sb,
		     const char * adaptertype = 0 )
	: m_fDensity(0)
	{ this->init(sb,adaptertype); }

    //! @brief Destroy an oTipsy object.
    virtual ~oTipsy() {}

    //! @brief Set the density output stream.
    void DensityOut( std::ostream &fDensity )
	{ m_fDensity = &fDensity; }

    //! @brief Write the Tipsy header.
    //! @param val A TipsyHeader structure.
    virtual oTipsy& operator<<(TipsyHeader &val);

    //! @brief Write a gas particle.
    //! @param val A TipsyGasParticle structure.
    virtual oTipsy& operator<<(TipsyGasParticle &val);

    //! @brief Write a star particle.
    //! @param val A TipsyStarParticle structure.
    virtual oTipsy& operator<<(TipsyStarParticle &val);

    //! @brief Write a dark particle.
    //! @param val A TipsyDarkParticle structure.
    virtual oTipsy& operator<<(TipsyDarkParticle &val);

    //! @brief Seek to the specified particle in the stream.
    //! @param pos The new position.
    virtual oTipsy &seekp( tipsypos pos );

    //! @brief Return the current particle position.
    //! @return The current particle position.
    virtual tipsypos tellp();
};

/** @brief Tipsy input file.
 */

//! The input file buffer type.
typedef std::basic_filebuf<char> ifilebuf_type;

class ifTipsy : public iTipsy {
protected:

    //! The private file buffer object.
    ifilebuf_type M_filebuf;

public:
    //! Construct an ifTipsy class.
    ifTipsy() : iTipsy( &M_filebuf ) {}

    //! @brief Construct an ifTipsy class and open the file.
    //! @param iname Input file name
    //! @param adaptertype Type of file (standard or native)
    ifTipsy( const char *iname,
	     const char * adaptertype = "standard" )
	: iTipsy( &M_filebuf ) {
	this->open(iname,adaptertype);
    }

    //! Get the internal buffer.
    //! @return A pointer to the file buffer.
    ifilebuf_type *rdbuf() const
      { return const_cast<ifilebuf_type*>(&M_filebuf); }

    //! Check if the file is open.
    //! @return true if the file is open, false if it is not.
    bool is_open()
	{ return M_filebuf.is_open(); }

    //! Open a Tipsy file.
    //! @param iname Input file name
    //! @param adaptertype Type of file (standard or native)
    void open(const char* iname, const char *adaptertype="standard" );

    //! Close the current file.
    void close() {
        M_filebuf.close();
    }

};

//! The output file buffer type.
typedef std::basic_filebuf<char> ofilebuf_type;

/** @brief A tipsy output file.
 */
class ofTipsy : public oTipsy {
protected:

    //! The private file buffer object.
    ofilebuf_type M_filebuf;
public:
    //! Construct an ofTipsy class and open the file.
    ofTipsy() : oTipsy( &M_filebuf ) {}

    //! Construct an ofTipsy class and open the file.
    //! @param oname Output file name
    //! @param adaptertype Type of file (standard or native)
    ofTipsy( const char * oname,
	     const char * adaptertype = "standard" )
	: oTipsy( &M_filebuf ) {
	this->open(oname,adaptertype);
    }

    //! Get the internal buffer.
    //! @return A pointer to the file buffer.
    ofilebuf_type *rdbuf() const
      { return const_cast<ofilebuf_type*>(&M_filebuf); }

    //! Check if the file is open.
    //! @return true if the file is open, false if it is not.
    bool is_open()
	{ return M_filebuf.is_open(); }

    //! Open a Tipsy file.
    //! @param oname Output file name
    //! @param adaptertype Type of file (standard or native)
    void open(const char* oname, const char *adaptertype="standard" );

    //! Close the current file.
    void close() {
        M_filebuf.close();
    }
};


#endif
