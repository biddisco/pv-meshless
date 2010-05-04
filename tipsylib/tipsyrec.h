/**
 *  @file
 *  @brief Disk format of the Tipsy file records
 *  @author Doug Potter
 *
 *  You must define the following types before including this file.
 *    disk_uint32_t
 *    disk_float
 *    disk_double
 *
 *  These structures define the physical layout of each record as they
 *  appear in the file.
 */

#ifndef TIPSYREC_H
#define TIPSYREC_H

//!  This structure is the exact layout of the header on disk.
typedef struct {
    disk_double   h_dTime;   //!< Expansion factor
    disk_uint32_t h_nBodies; //!< Total number of particles
    disk_uint32_t h_nDims;   //!< Number of dimensions (3)
    disk_uint32_t h_nSph;    //!< Number of gas particles
    disk_uint32_t h_nDark;   //!< Number of dark particles
    disk_uint32_t h_nStar;   //!< Number of star particles
    disk_uint32_t h_MBZ;     //!< A zero to pad the structure to 32 bytes.
    } tipsyrec_header;

//! This structure is the exact layout of a gas particle on disk.
typedef struct {
    disk_float mass;    //!< mass
    disk_float pos[3];  //!< position (x,y,z)
    disk_float vel[3];  //!< velocity (Vx,Vy,Vz)
    disk_float rho;     //!< rho
    disk_float temp;    //!< temp
    disk_float hsmooth; //!< hsmooth
    disk_float metals;  //!< metals
    disk_float phi;     //!< Potential
} tipsyrec_gas;

//! This structure is the exact layout of a dark particle on disk.
typedef struct {
    disk_float mass;    //!< mass
    disk_float pos[3];  //!< position (x,y,z)
    disk_float vel[3];  //!< velocity (Vx,Vy,Vz)
    disk_float eps;     //!< epsilon
    disk_float phi;     //!< Potential
} tipsyrec_dark;

//! This structure is the exact layout of a star particle on disk.
typedef struct {
    disk_float mass;    //!< mass
    disk_float pos[3];  //!< position (x,y,z)
    disk_float vel[3];  //!< velocity (Vx,Vy,Vz)
    disk_float metals;  //!< metals
    disk_float tform;   //!< tform
    disk_float eps;     //!< eps
    disk_float phi;     //!< Potential
} tipsyrec_star;

#endif
