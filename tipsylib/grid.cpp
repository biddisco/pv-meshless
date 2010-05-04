#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>
#include <math.h>
#include "grid.hpp"
#include "tiplim.hpp"

Grid::Grid()
    : m_GridX(0), m_GridSize(0), m_GridData(0),
      m_cx(0), m_cy(0), m_cz(0)
{
}

Grid::~Grid()
{
    if ( m_GridData ) delete [] m_GridData;
}

void Grid::init( int32_t GridSize )
{
    if ( m_GridData ) delete [] m_GridData;
    m_GridSize = GridSize;
    m_GridX = (m_GridSize+7) >> 3;
    m_GridData = new uint8_t [m_GridX * m_GridSize * m_GridSize];
}

void Grid::zero(void)
{
    memset(m_GridData,0,m_GridX * m_GridSize * m_GridSize);
}

static const uint8_t bitsSet[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
};

int32_t Grid::getGridCount() const
{
    int32_t i, count;

    count = 0;
    for( i=0; i<m_GridX*m_GridSize*m_GridSize; i++ )
	count += bitsSet[m_GridData[i]];

    return count;
}


void Grid::operator=(const Grid &rhs)
{
    init(rhs.m_GridSize);
    memcpy( m_GridData, rhs.m_GridData, m_GridX * m_GridSize * m_GridSize );
    //return *this;
}

void Grid::operator-=(const Grid &rhs)
{
    int32_t N = m_GridX * m_GridSize * m_GridSize;
    int32_t i;

    assert( m_GridSize == rhs.m_GridSize );

    for( i=0; i<N; i++ )
	m_GridData[i] &= ~rhs.m_GridData[i];

    //return *this;
}

void Grid::operator+=(const Grid &rhs)
{
    int32_t N = m_GridX * m_GridSize * m_GridSize;
    int32_t i;

    assert( m_GridSize == rhs.m_GridSize );

    for( i=0; i<N; i++ )
	m_GridData[i] |= rhs.m_GridData[i];

    //return *this;
}

void Grid::invert()
{
    int32_t N = m_GridX * m_GridSize * m_GridSize;
    int32_t i;

    for( i=0; i<N; i++ )
	m_GridData[i] = ~m_GridData[i];
}

bool Grid::saveGrid( std::ostream &f ) const
{
    uint8_t gridBlob[4];

    gridBlob[0] = m_GridSize & 0xff;
    gridBlob[1] = (m_GridSize>>8) & 0xff;
    gridBlob[2] = (m_GridSize>>16) & 0xff;
    gridBlob[3] = (m_GridSize>>24) & 0xff;

    f.write(reinterpret_cast<const char *>(&gridBlob), sizeof(gridBlob) );
    f.write(reinterpret_cast<const char *>(m_GridData),
	    m_GridX * m_GridSize * m_GridSize );
    return true;
}

bool Grid::saveGrid( const char *Name ) const
{
    std::ofstream f(Name);
    if ( !f.is_open() ) {
	return false;
    }
    return saveGrid(f);
}

bool Grid::loadGrid( std::istream &f )
{
    uint8_t gridBlob[4];

    f.read(reinterpret_cast<char *>(&gridBlob), sizeof(gridBlob) );
    m_GridSize = gridBlob[0] | (gridBlob[1]<<8) 
	| (gridBlob[2]<<16) | (gridBlob[3]<<24);
    m_GridX = (m_GridSize+7) >> 3;
    m_GridData = new uint8_t [m_GridX * m_GridSize * m_GridSize];

    f.read(reinterpret_cast<char *>(m_GridData),
	   m_GridX * m_GridSize * m_GridSize );
    return true;
}

bool Grid::loadGrid( const char *Name )
{
    std::ifstream f(Name);
    if ( !f.is_open() )
	return false;
    loadGrid(f);
    return true;
}

bool Grid::inGrid( float x, float y, float z )
{
    int32_t i, j, k;

    i = static_cast<int32_t>(floor( (x+0.5) * m_GridSize ));
    j = static_cast<int32_t>(floor( (y+0.5) * m_GridSize ));
    k = static_cast<int32_t>(floor( (z+0.5) * m_GridSize ));

    return getGrid(i,j,k);
}

// Make this a courser grid.  The idea is that if even one cell in the fine
// grid is set in a course grid, then set all of the fine grid cells.
void Grid::reduce( int32_t Factor )
{
    int32_t H, N, F, L;
    int32_t i, j, k;
    int32_t i0,i1,i2,j0,j1,j2,k0,k1,k2;

    H = m_GridSize / 2; // Half the grid size
    N = H / Factor;     // Number of course grid cells on each half
    F = H - Factor*N;   // First grid cell (ignoring rounding)
    L = H + Factor*N;   // Last
    assert( L-F == Factor*N*2);


    //std::cerr << "Reduce " << Factor << " F=" << F << " L=" << L << std::endl;

    for( k=F; k<L; k+=Factor ) {
	k1 = k==F ? 0 : k;
	k2 = (k+Factor) == L ? m_GridSize : (k+Factor);
	for( j=F; j<L; j+=Factor ) {
	    j1 = j==F ? 0 : j;
	    j2 = (j+Factor) == L ? m_GridSize : (j+Factor);
	    for( i=F; i<L; i+=Factor ) {
		i1 = i==F ? 0 : i;
		i2 = (i+Factor) == L ? m_GridSize : (i+Factor);
		bool bHave = false;

		for( k0=k1; k0<k2; k0++ )
		    for( j0=j1; j0<j2; j0++ )
			for( i0=i1; i0<i2; i0++ )
			    if ( getGrid(i0,j0,k0) ) bHave = true;

		if ( bHave ) {
		    for( k0=k1; k0<k2; k0++ )
			for( j0=j1; j0<j2; j0++ )
			    for( i0=i1; i0<i2; i0++ )
				setGrid(i0,j0,k0);
		}
	    }
	}
    }
}

//- Usenet echoes (21:200/1) -------------------------- COMP.GRAPHICS.ALGORITHMS -
// Msg  : 12 of 38                                                                
// From : eberly@cs.unc.edu                   2:5030/144.99   29 Jul 94  17:49:38 
// To   : All                                                 31 Jul 94  01:06:30 
// Subj : (1) Re: 3D Breshenhame lines?                                           
//--------------------------------------------------------------------------------
void Grid::drawLine (int x0, int y0, int z0, int x1, int y1, int z1)
{
    // starting point of line
    register int x = x0, y = y0, z = z0;

    // direction of line
    int dx = x1-x0, dy = y1-y0, dz = z1-z0;

    // increment or decrement depending on direction of line
    register int sx = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
    register int sy = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));
    register int sz = (dz > 0 ? 1 : (dz < 0 ? -1 : 0));

    // decision parameters for voxel selection
    dx = abs(dx); dy = abs(dy); dz = abs(dz);
    register int ax = 2*dx, ay = 2*dy, az = 2*dz;
    register int decx, decy, decz;

    // determine largest direction component, single-step related variable
    int max = dx, var = 0;
    if ( dy > max ) { max = dy; var = 1; }
    if ( dz > max ) { max = dz; var = 2; }

    // traverse Bresenham line
    switch ( var ) {
    case 0:  // single-step in x-direction
        for (decy=ay-dx, decz=az-dx; ; x += sx, decy += ay, decz += az) {
            // routine for displaying (x,y,z) goes here
	    setGrid(x,y,z);

            // take Bresenham step
            if ( x == x1 ) break;
            if ( decy >= 0 ) { decy -= ax; y += sy; }
            if ( decz >= 0 ) { decz -= ax; z += sz; }
        }
        break;
    case 1:  // single-step in y-direction
        for (decx=ax-dy, decz=az-dy; ; y += sy, decx += ax, decz += az) {
            // routine for displaying (x,y,z) goes here
	    setGrid(x,y,z);

            // take Bresenham step
            if ( y == y1 ) break;
            if ( decx >= 0 ) { decx -= ay; x += sx; }
            if ( decz >= 0 ) { decz -= ay; z += sz; }
        }
        break;
    case 2:  // single-step in z-direction
        for (decx=ax-dz, decy=ay-dz; ; z += sz, decx += ax, decy += ay) {
            // routine for displaying (x,y,z) goes here
	    setGrid(x,y,z);

            // take Bresenham step
            if ( z == z1 ) break;
            if ( decx >= 0 ) { decx -= az; x += sx; }
            if ( decy >= 0 ) { decy -= az; y += sy; }
        }
        break;
    }
}






#if 1
void Grid::fill( const Grid &S, int Fill )
{
    int32_t F2 = Fill*Fill;
    int32_t i,j,k;
    int32_t i2,j2,k2;

    init( S.m_GridSize );
    zero();

    // Yes: this is crazy
    for(k=0; k<m_GridSize; k++)
	for(j=0; j<m_GridSize; j++)
	    for(i=0; i<m_GridSize; i++)
		if ( S.getGrid(i,j,k) )
		    for( k2=k-Fill; k2<=k+Fill; k2++ )
			for( j2=j-Fill; j2<=j+Fill; j2++ )
			    for( i2=i-Fill; i2<=i+Fill; i2++ )
				if ((i2-i)*(i2-i)+(j2-j)*(j2-j)+(k2-k)*(k2-k)<=F2)
				    if ( S.getGrid(i2,j2,k2) )
					drawLine(i,j,k,i2,j2,k2);
}
#else
// This tries to connect any region
void Grid::fill( const Grid &S, int Fill )
{
    int32_t i,j,k,t;

    init( S.m_GridSize );
    zero();

    for(k=0; k<m_GridSize; k++)
	for(j=0; j<m_GridSize; j++)
	    for(i=0; i<m_GridSize; i++)
		if ( S.getGrid(i,j,k) ) {
		    setGrid(i,j,k);

		    for(t=i-Fill;t<i;t++) if ( S.getGrid(t,j,k) ) break;
		    for(;t!=i;t++) setGrid(t,j,k);
		    for(t=i+Fill;t>i;t--) if ( S.getGrid(t,j,k) ) break;
		    for(;t!=i;t--) setGrid(t,j,k);

		    for(t=j-Fill;t<j;t++) if ( S.getGrid(i,t,k) ) break;
		    for(;t!=j;t++) setGrid(i,t,k);
		    for(t=j+Fill;t>j;t--) if ( S.getGrid(i,t,k) ) break;
		    for(;t!=j;t--) setGrid(i,t,k);

		    for(t=k-Fill;t<k;t++) if ( S.getGrid(i,j,t) ) break;
		    for(;t!=k;t++) setGrid(i,j,t);
		    for(t=k+Fill;t>k;t--) if ( S.getGrid(i,j,t) ) break;
		    for(;t!=k;t--) setGrid(i,j,t);
		}
}
#endif

// This functions adds a border of "N" cells around any marked cell.
void Grid::addBorder( const Grid &S, int Boundary )
{
    int32_t i, j, k;
    int32_t i0,i1,i2,j0,j1,j2,k0,k1,k2;
    int32_t D2;

    init( S.m_GridSize );
    zero();

    D2 = (Boundary+1) * (Boundary+1);
    for(k=0; k<m_GridSize; k++)
	for(j=0; j<m_GridSize; j++)
	    for(i=0; i<m_GridSize; i++)
		if ( S.getGrid(i,j,k) ) {
		    i1 = i-Boundary; i2 = i+Boundary;
		    j1 = j-Boundary; j2 = j+Boundary;
		    k1 = k-Boundary; k2 = k+Boundary;

		    // The following is an important optimization
		    if ( S.getGrid(i-1,j,k) ) i1 = i;
		    if ( S.getGrid(i+1,j,k) ) i2 = i;
		    if ( S.getGrid(i,j-1,k) ) j1 = j;
		    if ( S.getGrid(i,j+1,k) ) j2 = j;
		    if ( S.getGrid(i,j,k-1) ) k1 = k;
		    if ( S.getGrid(i,j,k+1) ) k2 = k;

		    for( k0=k1; k0<=k2; k0++ )
			for( j0=j1; j0<=j2; j0++ )
			    for( i0=i1; i0<=i2; i0++ ) {
				int32_t di2 = (i0-i) * (i0-i);
				int32_t dj2 = (j0-j) * (j0-j);
				int32_t dk2 = (k0-k) * (k0-k);
				if ( di2 + dj2 + dk2 < D2 )
				    setGrid(i0,j0,k0);
			    }
		}
}

void Grid::findCenter()
{
    int i, j, k;

    Tipsy::limitTracker limits(0,m_GridSize);
    for( k=0; k<m_GridSize; k++ )
	for( j=0; j<m_GridSize; j++ )
	    for( i=0; i<m_GridSize; i++ )
		if ( getGrid(i,j,k) )
		    limits.Adjust(i,j,k);
    limits.Normalize();

    m_cx = m_GridSize/2 - static_cast<int32_t>(round(limits.X.getCenter()));
    m_cy = m_GridSize/2 - static_cast<int32_t>(round(limits.Y.getCenter()));
    m_cz = m_GridSize/2 - static_cast<int32_t>(round(limits.Z.getCenter()));
}

void GridBinner::init( int32_t GridSize )
{
    m_TopNode.splitNode(GridSize);
}

GridBinner::GridBinner( int GridSize )
    : m_cx(0), m_cy(0), m_cz(0)
{
    if ( GridSize )
	init(GridSize);
}

void GridBinner::create( int GridSize )
{
    init( GridSize );
}


GridBinner::~GridBinner()
{
}

GridBinner::Bin::Bin()
    : m_U(0)
{
    m_r[0] = m_r[1] = m_r[2] = 0.0;
    m_v[0] = m_v[1] = m_v[2] = 0.0;
}

GridBinner::Node::~Node() {
    if ( m_GridSize > 1 ) delete [] m_GridNode;
}


void GridBinner::Node::refine( int F, std::list<Bin> &Bins )
{
    int32_t i,j,k,N;
    if ( F == 1 ) {
	if ( !isBinned() ) {
	    Bins.push_back(Bin());
	    setBin(&Bins.back());
	}
    }
    else if ( m_GridSize == 0 ) {
	splitNode(F);
    }
    else {
	assert( F > 1 );
	assert( m_GridSize > 1 );
	N = F / m_GridSize;

	assert( N * m_GridSize == F );
	for( k=0; k<m_GridSize; k++ )
	    for( j=0; j<m_GridSize; j++ )
		for( i=0; i<m_GridSize; i++ )
		    (*this)(i,j,k).refine(N,Bins);
    }
}


bool GridBinner::Node::bin(Real *minR, Real *maxR,
			   Real *p, Real *r, Real *v, int32_t iMass)
{
    Real F[3];
    int32_t i[3],n;

    // Correct for periodic wrap (we put it back later).  The problem here is
    // that we are binning multiple particles.  If some of the particles have
    // wrapped, but others have not, then the center of mass becomes garbage.
    r[0] = fixWrap(r[0],p[0]);
    r[1] = fixWrap(r[1],p[1]);
    r[2] = fixWrap(r[2],p[2]);

    assert( m_GridSize != 0 );

    if ( m_GridSize == 1 ) {
	if ( m_GridBin == 0 ) return false;
	m_GridBin->add( r, v, iMass );
	return true;
    }

    //std::clog << "OKAY: " << minR[0] << " " << maxR[0] << std::endl;
    for( n=0; n<3; n++ ) {
	F[n] = (maxR[n] - minR[n]);
	// Convert vertex position into a grid position
	i[n] = static_cast<int32_t>(floorf(((p[n]-minR[n])/F[n])*getSize()));
	maxR[n] = minR[n] + (Real)(i[n]+1)*F[n]/getSize();
	minR[n] = minR[n] + (Real)i[n]*F[n]/getSize();
    }
#if 0
    if ( m_GridSize != 488 ) {
    std::cout << "Okay: sub bin x: " << i[0]
	      << " : " << minR[0] << " " << maxR[0] << " : " << p[0]
	      << " : " << F[0] << " " << getSize()
	      << std::endl;
    std::cout << "Okay: sub bin y: " << i[1]
	      << " : " << minR[1] << " " << maxR[1] << " : " << p[1]
	      << std::endl;
    std::cout << "Okay: sub bin z: " << i[2]
	      << " : " << minR[2] << " " << maxR[2] << " : " << p[2]
	      << std::endl;
    }
#endif
    return (*this)(i[0],i[1],i[2]).bin(minR,maxR,p,r,v,iMass);
}


void GridBinner::Node::fillInner( std::list<Bin> &Bins )
{
    int32_t N;
    int32_t i, j, k;
    N = getSize();
    if ( N == 0 ) {
	Bins.push_back(Bin());
	setBin(&Bins.back());
    }
    else if ( N > 1 ){
	for(k=0;k<N;k++ )
	    for( j=0; j<N; j++ )
		for( i=0; i<N; i++ )
		    (*this)(i,j,k).fillInner(Bins);
    }
}


void GridBinner::addSubRegion( const Grid *g, int32_t Factor )
{
    int32_t N;
    int32_t i, j, k;
    N = m_TopNode.getSize();

    for(k=0;k<N;k++ )
	for( j=0; j<N; j++ )
	    for( i=0; i<N; i++ )
		if ( g->getGrid(i,j,k) ) {
		    // This grid cell needs to be refined further
		    m_TopNode(i,j,k).refine(Factor,m_Bins);
		}
}


// Adds the specified grid and bins it according to Factor
void GridBinner::addRegionBIN( const Grid *g, int32_t Factor )
{
    int32_t H, N, F, L;
    int32_t i, j, k;
    int32_t i0,i1,i2,j0,j1,j2,k0,k1,k2;

    assert( Factor >= 1 );

    H = m_TopNode.getSize() / 2; // Half the grid size
    N = H / Factor;              // Number of course grid cells on each half
    F = H - Factor*N;            // First grid cell (ignoring rounding)
    L = H + Factor*N;            // Last
    assert( L-F == Factor*N*2);

    //std::cerr << "Region " << Factor << " F=" << F << " L=" << L << std::endl;


    for( k=F; k<L; k+=Factor ) {
	k1 = k==F ? 0 : k;
	k2 = (k+Factor) == L ? m_TopNode.getSize() : (k+Factor);

	for( j=F; j<L; j+=Factor ) {
	    j1 = j==F ? 0 : j;
	    j2 = (j+Factor) == L ? m_TopNode.getSize() : (j+Factor);

	    for( i=F; i<L; i+=Factor ) {
		i1 = i==F ? 0 : i;
		i2 = (i+Factor) == L ? m_TopNode.getSize() : (i+Factor);

		if ( m_TopNode(i1,j1,k1).isBinned() ) /*NADA*/;
		else if ( (g!=NULL && g->getGrid(i1,j1,k1))
		     || (g==NULL && !m_TopNode(i1,j1,k1).isBinned()) ) {
		    m_Bins.push_back(Bin());
		    for( k0=k1; k0<k2; k0++ )
			for( j0=j1; j0<j2; j0++ )
			    for( i0=i1; i0<i2; i0++ ) {
				m_TopNode(i0,j0,k0).setBin(&m_Bins.back());
			    }


		}
#if 0
		int bHave = 0;
		int bDont = 0;
		for( k0=k1; k0<k2; k0++ )
		    for( j0=j1; j0<j2; j0++ )
			for( i0=i1; i0<i2; i0++ ) {
			    if ( m_TopNode(i0,j0,k0).isBinned() ) bHave++;
			    else bDont++;
			}
		if ( bHave && bDont) {
		    std::clog << "OOPS: "
			      << "i1=" << i1 << " i2=" << i2 << " "
			      << "j1=" << j1 << " j2=" << j2 << " "
			      << "k1=" << k1 << " k2=" << k2 << " "
			      << "bHave=" << bHave << " bDont=" << bDont
			      << std::endl;

		}


		assert( !(bHave && bDont) );
#endif
	    }
	}
    }
}

void GridBinner::addRegionHR( const Grid &g )
{
    int32_t i, j, k;
    int32_t N = m_TopNode.getSize();

    for(k=0;k<N;k++)
	for(j=0;j<N;j++)
	    for(i=0;i<N;i++)
		if ( g.getGrid(i,j,k) ) {
		    m_TopNode(i,j,k).setBin(NULL);
		}
}

void GridBinner::addRegion( int32_t Factor )
{
    //std::cerr << "Factor " << Factor << " (final)" << std::endl;
    addRegionBIN( 0, Factor );
    std::clog << m_Bins.size() << " bins before filling the inner region" << std::endl;
    m_TopNode.fillInner( m_Bins );
    std::clog << m_Bins.size() << " bins after filling the inner region" << std::endl;
}


void GridBinner::addRegion( const Grid &g, int32_t Factor )
{
    //std::cerr << "Factor " << Factor << std::endl;
    if ( abs(Factor) < 1 )
	addRegionHR(g);
    else if ( Factor <= -1 )
	addSubRegion(&g,-Factor);
    else
	addRegionBIN(&g,abs(Factor));
}

bool GridBinner::bin( Real *p, Real *r, Real *v, int iMass )
{
    //int32_t i, j, k;
    Real minR[3], maxR[3];

    minR[0] = minR[1] = minR[2] = -0.5;
    maxR[0] = maxR[1] = maxR[2] =  0.5;

    // Recenter on the box
    p[0] = Recenter<Real>(p[0], m_cx, minR[0], maxR[0]);
    p[1] = Recenter<Real>(p[1], m_cy, minR[1], maxR[1]);
    p[2] = Recenter<Real>(p[2], m_cz, minR[2], maxR[2]);
    r[0] = Recenter<Real>(r[0], m_cx, minR[0], maxR[0]);
    r[1] = Recenter<Real>(r[1], m_cy, minR[1], maxR[1]);
    r[2] = Recenter<Real>(r[2], m_cz, minR[2], maxR[2]);

    return m_TopNode.bin(minR,maxR,p,r,v,iMass);

}

void GridBinner::sort()
{
    m_Bins.sort();
}

void GridBinner::getStart()
{
    m_bi = m_Bins.begin();
}
bool GridBinner::getNext( Real *r, Real *v, int32_t &iMass )
{
    for(;;) {
	if ( m_bi == m_Bins.end() ) return false;
	if ( m_bi->get(r,v,iMass) ) {
	    m_bi++;
	    periodic(r[0],r[1],r[2]);
	    return true;
	}
	m_bi++;
    }
}
