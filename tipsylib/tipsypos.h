/**
 *  @file
 *  @brief A position within a cosmological file
 *  @author Doug Potter
 */

#ifndef TIPSYPOS_H
#define TIPSYPOS_H

// I need a 64-bit integer
#ifdef _WIN32
   typedef __int64 int64_t;	// Define it from MSVC's internal type
   typedef unsigned __int64 uint64_t;
   typedef unsigned __int32 uint32_t;
   typedef unsigned    char uint8_t;
#else
   #include <stdint.h>		// Use the C99 official header
#endif

/** @brief A position in a Tipsy file
 */
class tipsypos {
    //! Private enumeration of section types.
    enum _tipsypos_section {
	_tipsypos_invalid  = 0, // The position is invalid
	_tipsypos_header   = 1, // Seek (or at) the header
	_tipsypos_gas      = 2, // Seek (or at) gas particles
	_tipsypos_dark     = 3, // Seek (or at) dark particles
	_tipsypos_star     = 4, // Seek (or at) star particles
	_tipsypos_particle = 5, // Seek (only) to absolute particle number
	_tipsypos_eof      = 6  // At EOF
    };
public:
    //! Type of a section.
    typedef _tipsypos_section section_type;

    //! Type of an offset.
    typedef uint64_t          offset_type;

    //! Section is inknown.
    static const section_type invalid = section_type(_tipsypos_invalid);
    //! Header section.
    static const section_type header  = section_type(_tipsypos_header);
    //! Gas particles section.
    static const section_type gas     = section_type(_tipsypos_gas);
    //! Dark particles section.
    static const section_type dark    = section_type(_tipsypos_dark);
    //! Star particles section.
    static const section_type star    = section_type(_tipsypos_star);
    //! Seek to an absolute particle number.
    static const section_type particle= section_type(_tipsypos_particle);
    //! At the end of file.
    static const section_type eof     = section_type(_tipsypos_eof);
protected:
    section_type m_section; //!< The current section.
    offset_type  m_offset;  //!< Offset within the section.

public:
    //! Assign another tipsypos to this tipsypos.
    tipsypos &assign( const tipsypos &rhs )
	{ m_section = rhs.m_section; m_offset = rhs.m_offset; return *this; }

    //! Construct a Tipsypos.
    tipsypos()
	{ m_section = invalid; m_offset = 0; }
    //! Construct a tipsypos at the given section and offset.
    tipsypos( section_type s, offset_type o = 0 )
	{ m_section = s; m_offset = o; }
    //! Copy constructor.
    tipsypos( const tipsypos &rhs )
	{ assign(rhs); }

    //! @brief Current section.
    //! @return The current section.
    inline section_type & section(void) { return m_section; }
    //! @brief Current offset.
    //! @return The current offset.
    inline offset_type  & offset(void)  { return m_offset; }
public:
    /**
     *  @brief Assign one tipsypos to another.
     *  @param rhs  The source tipsypos
     *  @return the destination position
     */
    inline tipsypos &operator=(const tipsypos &rhs)
	{ return assign(rhs); }

    /**
     *  @brief Compare one tipsypos to another.
     *  @param rhs  A tipsypos
     *  @return True if both are equal
     */
    inline bool operator==(const tipsypos &rhs)
	{ return m_section == rhs.m_section && m_offset==rhs.m_offset; }

    /**
     *  @brief Compare one tipsypos to a section
     *  @param s  A section
     *  @return True if the position is in the section.
     */
    inline bool operator==(section_type s)
	{ return m_section == s; }

};

#endif
