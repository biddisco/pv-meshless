/**
 *  @file
 *  @brief Bin Particles.
 *  @author Doug Potter
 */

#ifndef TIPSY_BINNER_H
#define TIPSY_BINNER_H
#include "ftipsy.hpp"
namespace Tipsy {
    // Base class to perform particle binning
    class Binner {
    protected:
	typedef double Real;

    protected:

	class BaseParticle {
	public:
	    Real mass;       //!< Mass of the particle
	    Real pos[3];     //!< Positions (x,y,z)
	    Real vel[3];     //!< Velocity (Vx,Vy,Vz)
	    Real phi;        //!< Potential
	    Real density;    //!< Density
	    Real eps;        //!< Softening (hsmooth for gas??)
	    long N;          //!< Number of particles binned here
	};


	class ParticleVector : public std::vector<BaseParticle> {};

	ParticleVector m_dark;
	ParticleVector m_gas;

    protected:
	float m_cx,m_cy,m_cz; //!< Center coordinates
	float m_sx,m_sy,m_sz; //!< Size of the box

	unsigned int m_gasparticle;
	unsigned int m_darkparticle;

	void Recenter(
	    float &x, float &y, float &z,
	    float cx, float cy, float cz ) const;
	void applyBoundary(
	    float &x, float &y, float &z ) const;
	void applyBoundary(
	    float &x, float &y, float &z,
	    float &rx, float &ry, float &rz) const;

	void Bin( ParticleVector &v,int n,TipsyBaseParticle &d,float eps=0.0 );
	//void Bin( int n, TipsyDarkParticle &d );
	//void Bin( int n, TipsyGasParticle &d );

	int Adjust(TipsyBaseParticle &d, float &x, float &y, float &z) const;


    public:
	// Number of dark particles that have been binned.
	unsigned int numDark() { return m_dark.size(); }

	// Number of gas particles that have been binned.
	unsigned int numGas() { return m_dark.size(); }

	// Returns the bin number of the given coordinates or -1.
	virtual int getBin(float x,float y,float z) const = 0;

    public:
	// Constructor.  Either specify the box information in the constructor,
	// or call setBox after the object is constructed.
	Binner();
	Binner( float cx, float cy, float cz,
		float sx=1.0, float sy=1.0, float sz=1.0 );
	virtual ~Binner();

	// Set the box information (binning center and box size)
	void setBox( float cx, float cy, float cz,
		     float sx=1.0, float sy=1.0, float sz=1.0 );

	// Bin the given particle
	bool Bin( TipsyDarkParticle &d );
	bool Bin( TipsyGasParticle &d );

	// As above, but put it into the bin based on x,y,z
	bool Bin( TipsyDarkParticle &d, float x, float y, float z );
	bool Bin( TipsyGasParticle &d, float x, float y, float z );

	// Returns true (and a single dark particle) or false
	bool getParticle( TipsyDarkParticle &d );
	bool getParticle( TipsyGasParticle &d );
    };

    class RadialBinner : public Binner {
    protected:
	class rFact {
	public:
	    Real m_Ax;
	    Real m_Ay;
	    Real m_Az;
	    rFact( Real Ax, Real Ay, Real Az )
		{ m_Ax = Ax; m_Ay = Ay; m_Az = Az; }
	    rFact( const rFact &rhs )
		{ m_Ax = rhs.m_Ax; m_Ay = rhs.m_Ay; m_Az = rhs.m_Az; }
	};

	class rFactList : public std::vector<rFact> {};

	friend class Cell;
	class Cell {
	protected:
	    int m_Nx;
	    int m_Ny;
	    int m_Nz;
	    Cell *m_Children;
	    int getIndex(int xo,int yo,int zo) const;
	public:
	    Cell();
	    void setSize( int Nx, int Ny=0, int Nz=0 );
	    void makeCells( RadialBinner &binner,
		int level, int &BinNo, const rFactList &rFactors,
		float x1, float y1, float z1,
		float x2, float y2, float z2 );
	    int getBin( float x, float y, float z ) const;
	};

	Cell m_TopCell;
	uint32_t  m_GridSize, m_GridX;
	uint8_t *m_GridData;

	float getOffset(
	    int nLevels,      //!< Number of refinements
	    int & N );        //!< Upper grid spacing

	float m_Rx, m_Ry, m_Rz; // Binning Radius
	// Offsets for special upper bin
	float m_Ox, m_Oy, m_Oz;
	int m_nBins;           // Number of radial bins

	// Precalculated
	rFactList m_rFactors;


    public:
	virtual int getBin(float x,float y,float z) const;

	RadialBinner();
	// Obsolete
	//void setBinning( int nTop, int nLevels,
	//		 float Rx, float Ry, float Rz );
	virtual void setBinning( int nLevels, int Bx, int By, int Bz,
			 float Rx, float Ry, float Rz );

	void loadGrid( const char *Name );

    };

}
#endif
