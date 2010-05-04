/**
 *  @file
 *  @brief Bin Particles
 *  @author Doug Potter
 */
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include "binner.hpp"

using namespace Tipsy;

void Binner::setBox(
    float cx, float cy, float cz,
    float sx, float sy, float sz )
{
    m_cx = cx;
    m_cy = cy;
    m_cz = cz;
    m_sx = sx;
    m_sy = sy;
    m_sz = sz;
}

Binner::Binner(
    float cx, float cy, float cz,
    float sx, float sy, float sz )
{
    setBox(cx,cy,cz,sx,sy,sz);
    m_gasparticle = m_darkparticle = 0;
}

Binner::Binner()
{
    m_cx = m_cy = m_cz = 0.0;
    m_sx = m_sy = m_sz = 0.0;
    m_gasparticle = m_darkparticle = 0;
}


Binner::~Binner()
{
}

void Binner::applyBoundary(
    float &x, float &y, float &z) const // Particle position
{
    float dx=x, dy=y, dz=z;
    applyBoundary(x,y,z,dx,dy,dz);
}

void Binner::applyBoundary(
    float &x, float &y, float &z,  // Particle position
    float &rx, float &ry, float &rz ) const
{
    assert( m_sx > 0.0 );
    assert( m_sy > 0.0 ); 
    assert( m_sz > 0.0 );

    // Apply the periodic boundary
    if ( x >= m_sx*0.5 )      { rx -= m_sx; x -= m_sx; }
    else if ( x < -m_sx*0.5 ) { rx += m_sx; x += m_sx; }
    if ( y >= m_sy*0.5 )      { ry -= m_sy; y -= m_sy; }
    else if ( y < -m_sy*0.5 ) { ry += m_sy; y += m_sy; }
    if ( z >= m_sz*0.5 )      { rz -= m_sz; z -= m_sz; }
    else if ( z < -m_sz*0.5 ) { rz += m_sz; z += m_sz; }

    assert( x >= -m_sx*0.5 && x < m_sx*0.5 );
    assert( y >= -m_sy*0.5 && y < m_sy*0.5 );
    assert( z >= -m_sz*0.5 && z < m_sz*0.5 );
}

void Binner::Recenter(
    float &x, float &y, float &z,        // Particle position
    float cx, float cy, float cz ) const // Center adjustment (m_cx,m_cy,m_cz)
{
    // Recenter the particle
    x -= cx;
    y -= cy;
    z -= cz;
}

void Binner::Bin( ParticleVector &v,int n,TipsyBaseParticle &d,float eps )
{
    // First, let's check to see that we have at least this many bins.
    if ( n >= (int)v.size() ) {
	BaseParticle empty;

	// Reserve space if the derived class wasn't kind enough to do so.
	if ( v.capacity() < (unsigned)(n+1) ) v.reserve(n+1);

	empty.mass = 0.0;
	empty.pos[0] = empty.pos[1] = empty.pos[2] = 0.0;
	empty.vel[0] = empty.vel[1] = empty.vel[2] = 0.0;
	empty.phi = empty.density = empty.eps = 0.0;
	empty.N = 0;

	while( n >= (int)v.size() )
	    v.push_back(empty);
    }

    // We can remove this later, but check to see we have treated the
    // periodic boundary condition properly.
    if ( v[n].N >= 1 ) {
	assert( fabs(v[n].pos[0]/v[n].mass-d.pos[0]) < 0.5*m_sx );
	assert( fabs(v[n].pos[1]/v[n].mass-d.pos[1]) < 0.5*m_sy );
	assert( fabs(v[n].pos[2]/v[n].mass-d.pos[2]) < 0.5*m_sz );
    }

    // Now add the particle to this bin
    v.at(n).mass+= d.mass;
    v[n].eps    += eps;
    v[n].pos[0] += d.mass * d.pos[0];
    v[n].pos[1] += d.mass * d.pos[1];
    v[n].pos[2] += d.mass * d.pos[2];
    v[n].vel[0] += d.mass * d.vel[0];
    v[n].vel[1] += d.mass * d.vel[1];
    v[n].vel[2] += d.mass * d.vel[2];
    v[n].N++;
}


//void Binner::Bin( int n, TipsyDarkParticle &d )
//{
//    Bin( m_dark, n, d, d.eps );
//}

//void Binner::Bin( int n, TipsyGasParticle &g )
//{
//    Bin( m_gas, n, g, g.hsmooth );
//}

int Binner::Adjust(TipsyBaseParticle &d, float &x, float &y, float &z) const
{
    int n;

    Recenter(x,y,z,m_cx,m_cy,m_cz);
    Recenter(d.pos[0],d.pos[1],d.pos[2],m_cx,m_cy,m_cz);
    applyBoundary(x,y,z,d.pos[0],d.pos[1],d.pos[2]);
    applyBoundary(d.pos[0],d.pos[1],d.pos[2]);

    // Get the bin number and return false if this should not be binned.
    n = getBin(x,y,z);
    if ( n < 0 ) {
	return n;
    }

    // Remove periodic boundary.  We add it back later
    if      ( d.pos[0] - x >  0.5 * m_sx ) d.pos[0] -= m_sx;
    else if ( d.pos[0] - x < -0.5 * m_sx ) d.pos[0] += m_sx;
    if      ( d.pos[1] - y >  0.5 * m_sy ) d.pos[1] -= m_sy;
    else if ( d.pos[1] - y < -0.5 * m_sy ) d.pos[1] += m_sy;
    if      ( d.pos[2] - z >  0.5 * m_sz ) d.pos[2] -= m_sz;
    else if ( d.pos[2] - z < -0.5 * m_sz ) d.pos[2] += m_sz;

    return n;
}


bool Binner::Bin( TipsyDarkParticle &d, float x, float y, float z )
{
    int n;
    n = Adjust(d,x,y,z);
    if ( n < 0 ) return false;
    Bin( m_dark, n, d, d.eps );
    return true;
}

bool Binner::Bin( TipsyGasParticle &d, float x, float y, float z )
{
    int n;
    n = Adjust(d,x,y,z);
    if ( n < 0 ) return false;
    Bin( m_gas, n, d, d.hsmooth );
    return true;
}

bool Binner::Bin( TipsyDarkParticle &d )
{
    return Bin(d,d.pos[0],d.pos[1],d.pos[2]);
}

bool Binner::Bin( TipsyGasParticle &d )
{
    return Bin(d,d.pos[0],d.pos[1],d.pos[2]);
}

bool Binner::getParticle( TipsyDarkParticle &d )
{
    while( m_darkparticle < m_dark.size() ) {

	// Only return bins that actually have particles
	if ( m_dark.at(m_darkparticle).mass > 0.0 ) {
	    d.mass   = m_dark[m_darkparticle].mass;
	    d.pos[0] = m_dark[m_darkparticle].pos[0] / d.mass;
	    d.pos[1] = m_dark[m_darkparticle].pos[1] / d.mass;
	    d.pos[2] = m_dark[m_darkparticle].pos[2] / d.mass;
	    d.vel[0] = m_dark[m_darkparticle].vel[0] / d.mass;
	    d.vel[1] = m_dark[m_darkparticle].vel[1] / d.mass;
	    d.vel[2] = m_dark[m_darkparticle].vel[2] / d.mass;
	    d.phi    = m_dark[m_darkparticle].phi;
	    d.density= m_dark[m_darkparticle].density;
	    d.eps    = m_dark[m_darkparticle].eps
		/ m_dark[m_darkparticle].N
		* logf(m_dark[m_darkparticle].N)
		/ logf(8);
	    //TEMP: Centering should be done outside
	    applyBoundary(d.pos[0],d.pos[1],d.pos[2]);
	    m_darkparticle++;
	    return true;
	}
	m_darkparticle++;
    }

    return false;
}

bool Binner::getParticle( TipsyGasParticle &d )
{
    while( m_gasparticle < m_gas.size() ) {

	// Only return bins that actually have particles
	if ( m_gas.at(m_gasparticle).mass > 0.0 ) {
	    d.mass   = m_gas[m_gasparticle].mass;
	    d.pos[0] = m_gas[m_gasparticle].pos[0] / d.mass;
	    d.pos[1] = m_gas[m_gasparticle].pos[1] / d.mass;
	    d.pos[2] = m_gas[m_gasparticle].pos[2] / d.mass;
	    d.vel[0] = m_gas[m_gasparticle].vel[0] / d.mass;
	    d.vel[1] = m_gas[m_gasparticle].vel[1] / d.mass;
	    d.vel[2] = m_gas[m_gasparticle].vel[2] / d.mass;
	    d.phi    = m_gas[m_gasparticle].phi;
	    d.density= m_gas[m_gasparticle].density;
	    d.hsmooth= m_gas[m_gasparticle].eps
		/ m_gas[m_gasparticle].N
		* logf(m_gas[m_gasparticle].N)
		/ logf(8);
	    //TEMP: Centering should be done outside
	    applyBoundary(d.pos[0],d.pos[1],d.pos[2]);
	    m_gasparticle++;
	    return true;
	}
	m_gasparticle++;
    }

    return false;
}

// Set the number of bins in each dimension.  They should result in a more
// or less perfect cube for each bin.  For a cosmological simulation,
// Sx == Sy == Sz because the top level box is a cube.
void RadialBinner::Cell::setSize( int Nx, int Ny, int Nz )
{
    m_Nx = Nx;
    m_Ny = Ny == 0 ? Nx : Ny;
    m_Nz = Nz == 0 ? Nx : Nz;
}

RadialBinner::Cell::Cell()
{
    m_Nx = m_Ny = m_Nz = 0;
    m_Children = 0;
}

int RadialBinner::Cell::getIndex(int xo,int yo,int zo) const
{
    return xo + m_Nx*yo + m_Nx*m_Ny*zo;
}

void RadialBinner::Cell::makeCells(
    RadialBinner &binner,
    int level, int &BinNo, const rFactList &rFactors,
    float x1, float y1, float z1,
    float x2, float y2, float z2 )
{
    int i, j, k, cell;
    float xc, yc, zc;
    float xmax, ymax, zmax;
    float xmin, ymin, zmin;
    float dx = (x2-x1) / m_Nx;
    float dy = (y2-y1) / m_Ny;
    float dz = (z2-z1) / m_Nz;

    m_Children = new Cell[m_Nx*m_Ny*m_Nz];

    for( i=0; i<m_Nx; i++ ) {
	xc = x1 + i*dx + 0.5*dx;
	xmax = std::max(fabs(x1 + i*dx ),fabs(x1 + (i+1)*dx));
	xmin = std::min(fabs(x1 + i*dx ),fabs(x1 + (i+1)*dx));

	for( j=0; j<m_Ny; j++ ) {
	    yc = y1 + j*dy + 0.5*dy;
	    ymax = std::max(fabs(y1 + j*dy ),fabs(y1 + (j+1)*dy));
	    ymin = std::min(fabs(y1 + j*dy ),fabs(y1 + (j+1)*dy));
	    for( k=0; k<m_Nz; k++ ) {
		zc = z1 + k*dz + 0.5*dz;
		zmax = std::max(fabs(z1 + k*dz ),fabs(z1 + (k+1)*dz));
		zmin = std::min(fabs(z1 + k*dz ),fabs(z1 + (k+1)*dz));

		cell = getIndex(i,j,k);
		//cell = i + m_Nx*j + m_Nx*m_Ny*k;

		// If the box is entirely outside Rn, then we bin it.
		if ( (xmin*xmin)*rFactors[level].m_Ax +
			  (ymin*ymin)*rFactors[level].m_Ay +
			  (zmin*zmin)*rFactors[level].m_Az >= 1.0 ) {
		    m_Children[cell].setSize(0,BinNo++,0);
		}

		// We cannot refine any further
		else if ( level == 0 ) {
		    // If this bin is entirely contained within R0, then ignore
		    if ( (xmax*xmax)*rFactors[0].m_Ax +
			 (ymax*ymax)*rFactors[0].m_Ay +
			 (zmax*zmax)*rFactors[0].m_Az <= 1.0 ) {

			// Check for the peanut.
			if ( binner.m_GridData ) {
			    int gx = static_cast<int>(
				floor(((xmax+xmin)/2 + 0.5) * binner.m_GridSize));
			    int gy = static_cast<int>(
				floor(((ymax+ymin)/2 + 0.5) * binner.m_GridSize));
			    int gz = static_cast<int>(
				floor(((zmax+zmin)/2 + 0.5) * binner.m_GridSize));
			    int O = (gx>>3) + 
				binner.m_GridX*(gy + binner.m_GridSize*gz);

			    if ( binner.m_GridData[O] & (1<<(gx&7)) ) {
				m_Children[cell].setSize(-1);
			    }
			    else {
				m_Children[cell].setSize(0,BinNo++,0);
			    }
			}
			else {
			    m_Children[cell].setSize(-1);
			}
		    }
		    // Okay, bin it
		    else {
			m_Children[cell].setSize(0,BinNo++,0);
		    }
		}

		// Okay, we can refine this further
		else {
		    m_Children[cell].setSize(2);
		    m_Children[cell].makeCells(
			binner,
			level-1, BinNo, rFactors,
			x1 +  i    * dx, y1 +  j    * dy, z1 +  k    * dz,
			x1 + (i+1) * dx, y1 + (j+1) * dy, z1 + (k+1) * dz );
		}
	    }
	}
    }

}



RadialBinner::RadialBinner()
    : m_GridSize(0), m_GridX(0), m_GridData(0)
{
    m_Ox = m_Oy = m_Oz = 0.0;
}

void RadialBinner::loadGrid( const char *Name )
{
    std::ifstream f(Name);
    if ( !f.is_open() ) {
	std::cerr << "Unable to open " << Name << std::endl;
	exit(11);
    }

    f.read(reinterpret_cast<char *>(&m_GridSize), sizeof(m_GridSize) );
    m_GridX = (m_GridSize+7) >> 3;
    m_GridData = new uint8_t [m_GridX * m_GridSize * m_GridSize];

    f.read(reinterpret_cast<char *>(m_GridData),
	   m_GridX * m_GridSize * m_GridSize );


}


float RadialBinner::getOffset(
    int nLevels, // Number of refinements
    int &N)      // Upper grid spacing
{
    int Ntotal = N;
    int Nticks;

    // Number of ticks at the greatest refinement level
    // e.g., nLevels = 4, N=488 => 488,244,122,61,30
    N = Ntotal >> nLevels;

    // How may will we use for refinement
    // e.g., nLevels = 4, N=488 => 30, 60, 120, 240, 480
    Nticks = N << nLevels;

    // The offset is just half the size of the difference
    // e.g., 0.5 * 8/488
    return 0.5 * (Ntotal-Nticks) / Ntotal;
}

void RadialBinner::setBinning( int nLevels, int Bx, int By, int Bz,
			       float Rx, float Ry, float Rz )
{
    double dx, dy, dz;
    int i, BinNo;

    // Calculate the offset for each dimension (normalized)
    dx = getOffset( nLevels, Bx ); // Bx changes
    dy = getOffset( nLevels, By ); // By changes
    dz = getOffset( nLevels, Bz ); // Bz changes

    // The size can now be slightly smaller than the original Bx,By,Bz
    m_TopCell.setSize(Bx,By,Bz);

    // Store the offsets in real units
    m_Ox = m_sx * dx;
    m_Oy = m_sy * dy;
    m_Oz = m_sz * dz;

    // Save the binning radius
    m_Rx = Rx;
    m_Ry = Ry;
    m_Rz = Rz;

    dx = ( m_sx*(0.5-1.0/Bx) - m_Ox - Rx) / (nLevels-1);
    dy = ( m_sy*(0.5-1.0/By) - m_Oy - Ry) / (nLevels-1);
    dz = ( m_sz*(0.5-1.0/Bz) - m_Oz - Rz) / (nLevels-1);

    // Calculate the binning radius factors, as in
    // x^2/rx^2 + ... where f=1/rx^2
    for( i=0; i< nLevels; i++ ) {
	m_rFactors.push_back(
	    rFact(1.0 / (Rx * Rx),
		  1.0 / (Ry * Ry),
		  1.0 / (Rz * Rz) )
	    );
	Rx += dx;
	Ry += dy;
	Rz += dz;
    }

    BinNo = 0;
    m_TopCell.makeCells(
	*this,
	nLevels-1, BinNo,
	m_rFactors,
	-m_sx*0.5+m_Ox, m_sx*0.5-m_Ox,
	-m_sy*0.5+m_Oy, m_sy*0.5-m_Oy,
	-m_sz*0.5+m_Oz, m_sz*0.5-m_Oz );

    // Well, we know how many bins we have, so prereserve the space.
    m_dark.reserve(BinNo);
    m_gas.reserve(BinNo);
}


// x, y and z are between 0 and 1, relative to this bin
int RadialBinner::Cell::getBin( float x, float y, float z ) const
{
    int xo, yo, zo, cell, ind;
    float dx, dy, dz;

    dx = 1.0 / m_Nx;
    dy = 1.0 / m_Ny;
    dz = 1.0 / m_Nz;

    xo = (int)floorf(x * m_Nx); if ( xo == m_Nx ) xo = m_Nx - 1;
    yo = (int)floorf(y * m_Ny); if ( yo == m_Ny ) yo = m_Ny - 1;
    zo = (int)floorf(z * m_Nz); if ( zo == m_Nz ) zo = m_Nz - 1;

    assert( xo >= 0 && xo < m_Nx );
    assert( yo >= 0 && yo < m_Ny );
    assert( zo >= 0 && zo < m_Nz );

    cell = getIndex(xo,yo,zo);
    //cell = xo + m_Nx*yo + m_Nx*m_Ny*zo;
    ind = m_Children[cell].m_Nx;
    if ( ind  == -1 ) return -1;
    else if ( ind == 0 ) return m_Children[cell].m_Ny;
    else {
	x = (x - xo*dx) * m_Nx;
	y = (y - yo*dy) * m_Ny;
	z = (z - zo*dz) * m_Nz;
	return m_Children[cell].getBin(x,y,z);
    }
}


int RadialBinner::getBin(float x,float y,float z) const
{
    // Convert coordinate to be from [0,size)
    x += 0.5*m_sx;
    y += 0.5*m_sy;
    z += 0.5*m_sz;
    assert( x >= 0 && x < m_sx );
    assert( y >= 0 && y < m_sy );
    assert( z >= 0 && z < m_sz );

    // Adjust in the upper box offset
    x -= m_Ox;
    y -= m_Oy;
    z -= m_Oz;

    // Normalize coordinates to between 0 and 1
    x /= m_sx - 2 * m_Ox;
    y /= m_sy - 2 * m_Oy;
    z /= m_sz - 2 * m_Oz;

    // If we had an offset, then we can just bin the outer particles.
    if ( x < 0.0 ) x = 0.0;
    else if ( x >= 1 ) x = 1;
    if ( y < 0.0 ) y = 0.0;
    else if ( y >= 1 ) y = 1;
    if ( z < 0.0 ) z = 0.0;
    else if ( z >= 1 ) z = 1;

    return m_TopCell.getBin(x,y,z);
}
