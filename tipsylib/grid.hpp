#include <fstream>
#include <stdint.h>
#include <assert.h>
#include <list>

template<typename T>
T Recenter( T r, T c, T b, T e ) {
    r -= c;
    if ( r < b ) r += e-b;
    else if ( r >= e ) r -= e-b;
    return r;
}

class Grid {
protected:
    int32_t m_GridX, m_GridSize;
    uint8_t *m_GridData;
    int32_t m_cx, m_cy, m_cz;

    inline void periodic(int32_t &i, int32_t &j, int32_t &k) const {
	i = Recenter<int32_t>(i, m_cx, 0, m_GridSize);
	j = Recenter<int32_t>(j, m_cy, 0, m_GridSize);
	k = Recenter<int32_t>(k, m_cz, 0, m_GridSize);
    }

    void drawLine (int x0, int y0, int z0, int x1, int y1, int z1);

public:

    Grid();
    virtual ~Grid();

    void operator=(const Grid &rhs);
    void operator-=(const Grid &rhs);
    void operator+=(const Grid &rhs);

    void invert(void);

    virtual void init( int32_t GridSize );
    virtual void zero(void);

    void findCenter();
    uint32_t getCX() const { return m_cx; }
    uint32_t getCY() const { return m_cy; }
    uint32_t getCZ() const { return m_cz; }
    void setCenter( int32_t cx, int32_t cy, int32_t cz )
	{ m_cx=cx; m_cy=cy; m_cz=cz; }

    int32_t getGridSize() const { return m_GridSize; }
    int32_t getGridCount() const;

    bool saveGrid( std::ostream &f ) const;
    bool saveGrid( const char *Name ) const;

    bool loadGrid( std::istream &f );
    bool loadGrid( const char *Name );

    inline void setGrid( int32_t i, int32_t j, int32_t k ) {
	periodic(i,j,k);
	assert(i>=0&&i<m_GridSize&&j>=0&&j<m_GridSize&&k>=0&&k<m_GridSize);
	int32_t O = i/8 + m_GridX * ( j + m_GridSize * k );
	m_GridData[O] |= (1 << (i&7));
    }

    inline bool getGrid( int32_t i, int32_t j, int32_t k) const {
	periodic(i,j,k);
	assert(i>=0&&i<m_GridSize&&j>=0&&j<m_GridSize&&k>=0&&k<m_GridSize);
	int32_t O = i/8 + m_GridX * ( j + m_GridSize * k );
	return (m_GridData[O] & (1 << (i&7))) != 0;
    }

    bool inGrid( float x, float y, float z );

    void addBorder( const Grid &S, int Boundary = 1 );
    void fill( const Grid &S, int Fill = 5 );
    void reduce( int32_t Factor );

};

class GridBinner {
public:
    typedef float Real;

protected:
    class Bin {
    protected:
	Real m_r[3];     //!< Positions (x,y,z)
	Real m_v[3];     //!< Velocity (Vx,Vy,Vz)
	int64_t m_U;     //!< Number of mass resolution units
    public:
	Bin();
	inline void add( Real *r, Real *v, int64_t U ) {
	    m_r[0] += r[0]*U; m_r[1] += r[1]*U; m_r[2] += r[2]*U;
	    m_v[0] += v[0]*U; m_v[1] += v[1]*U; m_v[2] += v[2]*U;
	    m_U += U;
	}

	bool get( Real *r, Real *v, int32_t &iMass ) const {
	    if ( m_U == 0 ) return false;
	    r[0] = m_r[0] / m_U;
	    r[1] = m_r[1] / m_U;
	    r[2] = m_r[2] / m_U;
	    v[0] = m_v[0] / m_U;
	    v[1] = m_v[1] / m_U;
	    v[2] = m_v[2] / m_U;
	    iMass = m_U;
	    return true;
	}


	bool operator<(const Bin &rhs) const {
	    return m_U < rhs.m_U;
	}
    };


    class Node {
    protected:
	union {               // Neither is valid with m_GridSize == 0
	    Bin  *m_GridBin;  // When m_GridSize == 1
	    Node *m_GridNode; // When m_GridSize > 1
	};
	int32_t m_GridSize;

	inline Real fixWrap( Real particle, Real grid ) {
	    Real delta = grid - particle;
	    if ( delta > 0.5 ) particle += 1.0;
	    else if ( delta < -0.5 ) particle -= 1.0;
	    return particle;
	}

    public:
	Node() : m_GridNode(0), m_GridSize(0)  {}
	~Node();

	inline int32_t getSize() const { return m_GridSize; }
	inline int32_t isBinned() { return getSize() != 0; }

	Node &operator()( int32_t i, int32_t j, int32_t k ) {
	    assert( m_GridSize > 1 );

	    int32_t N = getSize();
	    if (i < 0) i += N; else if (i >= N) i -= N;
	    if (j < 0) j += N; else if (j >= N) j -= N;
	    if (k < 0) k += N; else if (k >= N) k -= N;
	    assert( i>=0 && i<N && j>=0 && j<N && k>=0 && k<N );
	    return m_GridNode[i + m_GridSize * (j + m_GridSize*k)];
	}

	void setBin( Bin *b ) {
	    if ( m_GridSize != 0 ) {
		std::clog << "ASSERT: m_GridSize = " << m_GridSize << std::endl;
	    }

	    assert( m_GridSize == 0 );
	    m_GridSize = 1;
	    m_GridBin = b;
	}
	void splitNode( int iSize ) {
	    assert( m_GridSize == 0 );
	    m_GridSize = iSize;
	    assert( m_GridSize > 1 );
	    m_GridNode = new Node[iSize*iSize*iSize];
	}
	void refine( int F, std::list<Bin> &Bins );
	bool bin(Real *minR, Real *maxR, Real *p, Real *r, Real *v, int32_t iMass);

	void fillInner( std::list<Bin> &Bins );
    };

    Node m_TopNode;

    Real m_cx, m_cy, m_cz;
    Bin m_SentinalHR;
    std::list<Bin> m_Bins;
    std::list<Bin>::iterator m_bi;

protected:
    inline void periodic(int32_t &i, int32_t &j, int32_t &k) const {
	int32_t N = m_TopNode.getSize();
	if (i < 0) i += N; else if (i >= N) i -= N;
	if (j < 0) j += N; else if (j >= N) j -= N;
	if (k < 0) k += N; else if (k >= N) k -= N;
    }
    inline void periodic(Real &x, Real &y, Real &z) const {
	if (x < -0.5) x += 1.0; else if (x >= 0.5) x -= 1.0;
	if (y < -0.5) y += 1.0; else if (y >= 0.5) y -= 1.0;
	if (z < -0.5) z += 1.0; else if (z >= 0.5) z -= 1.0;
    }

    virtual void init( int32_t GridSize );

    void addSubRegion( const Grid *g, int32_t Factor );
    void addRegionHR(  const Grid &g );
    void addRegionBIN( const Grid *g, int32_t Factor );

public:
    GridBinner( int GridSize=0 );
    virtual ~GridBinner();

    void create( int GridSize );

    void addRegion( const Grid &g, int32_t Factor );
    void addRegion( int32_t Factor );

    void setCenter( int32_t cx, int32_t cy, int32_t cz ) {
	assert( m_TopNode.getSize() > 1 );

	m_cx = -(Real)cx/m_TopNode.getSize();
	m_cy = -(Real)cy/m_TopNode.getSize();
	m_cz = -(Real)cz/m_TopNode.getSize();
    }
    inline Real getCX() { return m_cx; }
    inline Real getCY() { return m_cy; }
    inline Real getCZ() { return m_cz; }

    bool bin( Real *p, Real *r, Real *v, int32_t iMass );
    int32_t binCount() { return m_Bins.size(); }

    void sort();
    void getStart();
    bool getNext( Real *r, Real *v, int32_t &iMass );

};
