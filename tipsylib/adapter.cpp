/**
 *  @file
 *  @brief Generic cosmological file format adapter.
 *  @author Doug Potter
 */

#include "adapter.h"

TipsyAdapter::TipsyAdapter()
{
    // The first value is always zero because ID 0 means that the field is
    // not supported.
    m_values.resize(1);
    m_values[0] = 0.0f;
}

TipsyAdapter::~TipsyAdapter()
{
}

std::size_t TipsyAdapter::getFieldID( const std::string &field )
{
    idmap_t::iterator i;
    i=m_idmap.find(field);
    if ( i == m_idmap.end() ) return 0;
    return i->second;
}

/** @brief Return the current get position
 *  @return The current file position.
 */
tipsypos TipsyAdapter::tellg( void )
{
    //! By default we just return the position as we know it.  This function
    //! can be overridden, but probably won't be.
    return m_position;
}

/** @brief Return the current put position
 *  @return The current file position.
 */
tipsypos TipsyAdapter::tellp( void )
{
    //! By default we just return the position as we know it.  This function
    //! can be overridden, but probably won't be.
    return m_position;
}
