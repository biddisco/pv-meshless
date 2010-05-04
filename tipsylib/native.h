/**
 *  @file
 *  @brief Adapter for Tipsy native format files
 *  @author Doug Potter
 */
#ifndef TIPSY_NATIVE_H
#define TIPSY_NATIVE_H

#include <istream>
#include <ostream>
#include "adapter.h"

//! @brief Read a Tipsy file in "native" format
class TipsyNativeAdapter : public TipsyAdapter {
    //! Private enumeration of field ids.
    enum _tipsyid {
	id_mass    = 1,
	id_x       = 2,
	id_y       = 3,
	id_z       = 4,
	id_vx      = 5,
	id_vy      = 6,
	id_vz      = 7,
	id_rho     = 8,
	id_temp    = 9,  // also id_tform
	id_eps     = 10,
	id_metals  = 11,
	id_phi     = 12,
	id_density = 13,
	count_id   = 14
    };

public:
    static const _tipsyid idmass    = _tipsyid(id_mass);   //!< mass id
    static const _tipsyid idx       = _tipsyid(id_x);      //!< x id
    static const _tipsyid idy       = _tipsyid(id_y);      //!< y id
    static const _tipsyid idz       = _tipsyid(id_z);      //!< z id
    static const _tipsyid idvx      = _tipsyid(id_vx);     //!< Vx id
    static const _tipsyid idvy      = _tipsyid(id_vy);     //!< Vy id
    static const _tipsyid idvz      = _tipsyid(id_vz);     //!< Vz id
    static const _tipsyid idrho     = _tipsyid(id_rho);    //!< rho id
    static const _tipsyid idtemp    = _tipsyid(id_temp);   //!< temp id
    static const _tipsyid idtform   = _tipsyid(id_temp);   //!< tform id
    static const _tipsyid ideps     = _tipsyid(id_eps);    //!< eps id
    static const _tipsyid idhsmooth = _tipsyid(id_eps);    //!< smooth id
    static const _tipsyid idmetals  = _tipsyid(id_metals); //!< metals id
    static const _tipsyid idphi     = _tipsyid(id_phi);    //!< phi id
    static const _tipsyid iddensity = _tipsyid(id_density);//!< density id
    static const _tipsyid idcountid = _tipsyid(count_id);  //!< Greatest id

public:
    //! Stream buffer type
    typedef std::basic_streambuf<char> streambuf_type;

protected:
    //! Reference to the stream buffer object.
    streambuf_type *m_sb;

protected:
    //! @brief Move forward a single particle.
    void forward(void);

    //! @brief Initialize this object.
    //! @param sb The streambuf object to use.
    void init( streambuf_type *sb );

public:
    //! @brief Construct a native adapter.
    //! @param sb The streambuf object to use.
    explicit TipsyNativeAdapter( streambuf_type *sb );

    //! @brief Read the next record (header or particle).
    virtual void getNext();

    //! @brief Write the next record (header or particle).
    virtual void putNext();

    //! @brief Seek the get pointer to a specific particle.
    //! @param pos The particle to seek to.
    virtual tipsypos seekg( tipsypos &pos );

    //! @brief Seek the put pointer to a specific particle.
    //! @param pos The particle to seek to.
    virtual tipsypos seekp( tipsypos &pos );
};

#endif
