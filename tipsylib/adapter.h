/**
 *  @file
 *  @brief Generic cosmological file format adapter.
 *  @author Doug Potter
 */

#ifndef TIPSY_ADAPTER_H
#define TIPSY_ADAPTER_H

#include <string>
#include <map>
#include <vector>
#include "tipsypos.h"

/**
 *  @brief  Abstract implementation of a file adapter.
 *
 *  All supported file formats will use this class as the base class.
 */
class TipsyAdapter {
protected:
    //! The type of our field name to id map.
    typedef std::map<std::string,std::size_t> idmap_t;

    //! The values vector type.
    typedef std::vector<float> values_t;

    //! Maps between a field name (e.g., "mass") and the corresponding id.
    idmap_t  m_idmap;

    //! A list of the values of each field.
    values_t m_values;

    //! The current position in the file
    tipsypos m_position;

    //FIXME: should not be public.  should be a map like m_values.
public:
    double                m_dTime;   //!< Expansion factor (simulation time).
    int                   m_nDims;   //!< Number of dimensions (2 or 3).
    tipsypos::offset_type m_nBodies; //!< Total number of particles.
    tipsypos::offset_type m_nSph;    //!< Number of gas particles.
    tipsypos::offset_type m_nDark;   //!< Number of dark particles.
    tipsypos::offset_type m_nStar;   //!< Number of star particles.

public:
    /** @brief Returns the id of the specified field
     *  @param field The name of the field (e.g., "mass").
     *  @return The field id (index) or zero if not present.
     */
    std::size_t getFieldID( const std::string &field );

    /** @brief Get or set a specified field.
     *  @param fid Field id (returned from getFieldID).
     *  @return A reference to the field.
     */
    inline values_t::reference Field( std::size_t fid )
	{ return m_values[fid]; }

public:
    /** @brief Construct an adapter.
     */
    TipsyAdapter();

    /** @brief Destroy an adapter.
     */
    virtual ~TipsyAdapter();

    /** @brief Read the next particle from the current position
     */
    virtual void getNext(void) = 0;

    /** @brief Write the next particle at the current position
     */
    virtual void putNext(void) = 0;

    /** @brief Seek to the specified particle
     *  @return The new position (if successful), else the old position.
     */
    virtual tipsypos seekg( tipsypos &pos ) = 0;

    /** @brief Return the current get position
     *  @return The current file position.
     */
    virtual tipsypos tellg( void );

    /** @brief Seek to the specified particle
     *  @return The new position (if successful), else the old position.
     */
    virtual tipsypos seekp( tipsypos &pos ) = 0;

    /** @brief Return the current pu position
     *  @return The current file position.
     */
    virtual tipsypos tellp( void );
};

#endif
