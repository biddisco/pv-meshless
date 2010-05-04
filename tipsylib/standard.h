/**
 *  @file
 *  @brief Adapter for Tipsy standard format files
 *  @author Doug Potter
 */
#ifndef TIPSY_STANDARD_H
#define TIPSY_STANDARD_H

#include <istream>
#include <ostream>
#include "native.h"

//! Read a Tipsy file in "standard" format
class TipsyStandardAdapter : public TipsyNativeAdapter {
public:
    //! @brief Construct a Tipsy standard adapter.
    //! @param sb The streambuf object to use.
    TipsyStandardAdapter( streambuf_type *sb );

    //! @brief Read the next particle from the stream.
    virtual void getNext();

    //! @brief Write the next particle to the stream.
    virtual void putNext();
};

#endif
