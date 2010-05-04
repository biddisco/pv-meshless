/**
 *  @file
 *  @brief Box Statistics (ala Tipsy)
 *  @author Doug Potter
 */

#include <iostream>
#include <limits>
#include <map>
#include <list>
#include <getopt.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include "ftipsy.hpp"
#include "tiplim.hpp"
#include "grid.hpp"

#define OPT_HELP        'h'
#define OPT_PERIODIC    'p'
#define OPT_NONPERIODIC 'n'
#define OPT_MARK        'm'
#define OPT_BOXUNIT     '9'
#define OPT_DOUBLE      'd'

#define OPT_GRID_SIZE   'g'
#define OPT_GRID_FILE   'o'
#define OPT_GRID_BOUND  'b'
#define OPT_GRID_SLOP   's'
#define OPT_CENTER      'c'


class ZeroList : public std::list<unsigned int>
{
};

class MassMap : public std::map<float,unsigned int>
{
public:
    inline void addMass( float mass );
};

void MassMap::addMass( float mass ) {
    iterator i = find(mass);

    if ( i != end() )
	i->second++;
    else
	(*this)[mass] = 1;
}


class Info : public Tipsy::limitTracker {
public:
    double Mass;
    double Angx, Angy, Angz;
    double Mx, My, Mz;
    double Vx, Vy, Vz;

    Info();

    void Accumulate( TipsyBaseParticle &p );
    void Normalize(void);
};

class AllInfo {
public:
    Info info, ginfo, dinfo, sinfo, binfo;
    MassMap massinfo;
    ZeroList zeros;

    int32_t GridSize;
    int32_t GridX;
    int32_t GridBound;
    float   GridSlop;
    Grid    grid;
    //uint8_t *GridData;
    float m_cx, m_cy, m_cz;
    //std::ofstream GridFile;

public:
    void setLimits( float left, float right );

    void CommonAcc( TipsyBaseParticle &p, unsigned int i );
    void Accumulate( TipsyGasParticle &p, unsigned int i );
    void Accumulate( TipsyDarkParticle &p, unsigned int i );
    void Accumulate( TipsyStarParticle &p, unsigned int i );

    void ReadNext( iTipsy &in, unsigned int i );

    void Normalize(void);

    AllInfo();
    void setGrid( const char *Name, int32_t Size,
		  int32_t Boundary, float slop,
		  float cx, float cy, float cz );
    void adjGrid( TipsyBaseParticle &p );
    void writeGrid( const char *Name);
    int32_t countGrid();

};

AllInfo::AllInfo()
    : GridSize(0), GridX(0), GridSlop(0)
{
}

void AllInfo::setGrid( const char *Name, int32_t Size,
		       int32_t Boundary, float slop,
		       float cx, float cy, float cz )
{
    GridSize = Size;
    GridSlop = slop / GridSize;
    GridBound = Boundary;
    GridX = (GridSize+7) >> 3;
    if ( GridSize == 0 ) return;

    // Round to the best grid position
    std::clog << "OLD: cx=" << cx << " cy=" << cy << " cz=" << cz << std::endl;
    m_cx = roundf( cx * GridSize ) / GridSize;
    m_cy = roundf( cy * GridSize ) / GridSize;
    m_cz = roundf( cz * GridSize ) / GridSize;
    std::clog << "NEW: cx=" << m_cx << " cy=" << m_cy << " cz=" << m_cz << std::endl;

    grid.init(GridSize);
    grid.zero();
}

void AllInfo::writeGrid( const char *Name )
{
    if ( GridSize == 0 ) return;
    grid.saveGrid( Name );
//    GridFile.write( reinterpret_cast<char *>(GridData),
//		    GridX * GridSize * GridSize );
}

void AllInfo::adjGrid( TipsyBaseParticle &p )
{
    int32_t i, j, k;
    int32_t i1, i2, j1, j2, k1, k2;

    if ( GridSize == 0 )
	return;

    // Find the start and ending grid postions
    i1=static_cast<int32_t>(floor((p.pos[0]-m_cx+0.5-GridSlop)*GridSize))-GridBound;
    i2=static_cast<int32_t>(floor((p.pos[0]-m_cx+0.5+GridSlop)*GridSize))+GridBound;
    j1=static_cast<int32_t>(floor((p.pos[1]-m_cx+0.5-GridSlop)*GridSize))-GridBound;
    j2=static_cast<int32_t>(floor((p.pos[1]-m_cx+0.5+GridSlop)*GridSize))+GridBound;
    k1=static_cast<int32_t>(floor((p.pos[2]-m_cx+0.5-GridSlop)*GridSize))-GridBound;
    k2=static_cast<int32_t>(floor((p.pos[2]-m_cx+0.5+GridSlop)*GridSize))+GridBound;

    for( i=i1; i<=i2; i++ ) {
	for( j=j1; j<=j2; j++ ) {
	    for( k=k1; k<=k2; k++ ) {
		grid.setGrid(i,j,k);
	    }
	}
    }
}

int32_t AllInfo::countGrid()
{
    int i,j,k;
    int S;

    S = 0;
//    for( i=0; i<GridX*GridSize*GridSize; i++ )
//	for( j=0; j<8; j++ )
//	    if ( GridData[i]&(1<<j) ) S++;
    for( k=0; k<GridSize; k++ )
	for( j=0; j<GridSize; j++ )
	    for( i=0; i<GridSize; i++ )
		if ( grid.getGrid(i,j,k) ) S++;
    return S;
}


void AllInfo::setLimits( float left, float right )
{
    info.setLimits(left,right);
    ginfo.setLimits(left,right);
    dinfo.setLimits(left,right);
    sinfo.setLimits(left,right);
    binfo.setLimits(left,right);
}



void AllInfo::CommonAcc( TipsyBaseParticle &p, unsigned int i )
{
    if ( p.mass < std::numeric_limits<float>::min() ) {
	zeros.push_back(i);
    }
    massinfo.addMass(p.mass);
    info.Accumulate(p);
    adjGrid(p);
}

void AllInfo::Accumulate( TipsyGasParticle &p, unsigned int i )
{
    CommonAcc(p,i);
    ginfo.Accumulate(p);
    binfo.Accumulate(p);
}

void AllInfo::Accumulate( TipsyDarkParticle &p, unsigned int i )
{
    CommonAcc(p,i);
    dinfo.Accumulate(p);
}

void AllInfo::Accumulate( TipsyStarParticle &p, unsigned int i )
{
    CommonAcc(p,i);
    sinfo.Accumulate(p);
    binfo.Accumulate(p);
}

void AllInfo::ReadNext( iTipsy &in, unsigned int i )
{
    TipsyGasParticle  g;
    TipsyDarkParticle d;
    TipsyStarParticle s;

    switch( in.tellg().section() ) {
    case tipsypos::gas:
	in >> g;
	Accumulate(g,i);
	break;
    case tipsypos::dark:
	in >> d;
	Accumulate(d,i);
	break;
    case tipsypos::star:
	in >> s;
	Accumulate(s,i);
	break;
    default:
	assert(0);
	break;
    }
}

void AllInfo::Normalize(void)
{
    info.Normalize();
    ginfo.Normalize();
    dinfo.Normalize();
    sinfo.Normalize();
    binfo.Normalize();
}



Info::Info()
{
    Mass = Mx = My = Mz = Vx = Vy = Vz = Angx = Angy = Angz = 0.0;
}

void Info::Accumulate( TipsyBaseParticle &p )
{
    Mass += p.mass;

    Mx += p.mass * p.pos[0];
    My += p.mass * p.pos[1];
    Mz += p.mass * p.pos[2];
    Vx += p.mass * p.vel[0];
    Vy += p.mass * p.vel[1];
    Vz += p.mass * p.vel[2];

    // Angular Momentum = ( pos X vel ) * mass
    Angx += p.mass * (p.pos[1]*p.vel[2] - p.pos[2]*p.vel[1]);
    Angy += p.mass * (p.pos[2]*p.vel[0] - p.pos[0]*p.vel[2]);
    Angz += p.mass * (p.pos[0]*p.vel[1] - p.pos[1]*p.vel[0]);

    Adjust(p.pos);
}

void Info::Normalize(void)
{
    if ( Mass < std::numeric_limits<float>::min() /* Mass == 0.0 */  ) {
	Mx = My = Mz = Vx = Vy = Vz = 0.0;
    }
    else {
	Mx /= Mass;  My /= Mass;  Mz /= Mass;
	Vx /= Mass;  Vy /= Mass;  Vz /= Mass;
	Angx /= Mass;  Angy /= Mass;  Angz /= Mass;
    }

    Tipsy::limitTracker::Normalize();
}

void Usage( const char *name ) {
    std::cerr
	<< "Usage: " << name
	<< " [-hmnp] TIPSY [all|dark|baryon|gas|star]" << std::endl
	<< "    -p,--periodic=1.0\t\tBox is periodic (default with MARK)"
	<< std::endl
	<< "    -n,--nonperiodic\t\tBox is nonperiodic (default without)"
	<< std::endl
	<< "    -m,--mark=MARKFILE\t\tReport only on marked particles"
	<< std::endl
	<< "    --boxunit=n\t\t\tMultiply units by the box unit"
	<< std::endl
	<< "    -g,--grid-size=n\t\tOutput a GRID of size n"
	<< std::endl
	<< "    -b,--boundary=1\t\tSurround the GRID with a boundary"
	<< std::endl
	<< "    -o,--output=grid.dat\tOut name (default grid.dat)"
	<< std::endl
	<< "    -c,--center=x,y,z\t\tCorrect for recentering"
	<< std::endl
	<< "    -d,--doublei\t\t\tDouble precision input"
	<< std::endl;
}

static void BadCenter(void) {
    std::clog << "Specify center as --center=x,y,z" << std::endl;
    exit(1);
}

int main( int argc, char *argv[] ) {
    bool bHelp     = false; //!< Help was requested on the command line
    bool bError    = false; //!< An error occurred on the command line
    bool bPeriodSet= false; //!< True if the periodicity has been set
    float BoxUnit  = 1.0;   //!< Box units to scale output display
    float cx = 0.0;
    float cy = 0.0;
    float cz = 0.0;
    int32_t GridSize = 0;   //!< Size of grid box
    int32_t GridBound = 1;  //!< Boundary around grid
    float   GridSlop = 0.0; //!< How close to a boundary before counting it
    std::string nameMark;   //!< Name of the MARK file (optional)
    std::string nameTipsy;  //!< Tipsy input file name
    const char *boxtype = "all";
    const char *GridFile = "grid.dat";
    const char *Adapter = "standard";

    ifTipsy in;

    TipsyHeader       h;

    unsigned int i;
    AllInfo info;
    Info *ip;

    //! Parse command line arguments (flags).
    for(;;) {
        int c, option_index=0;
	char *p;

        static struct option long_options[] = {
            { "help",       0, 0, OPT_HELP },
            { "mark",       1, 0, OPT_MARK },
            { "periodic",   2, 0, OPT_PERIODIC },
            { "nonperiodic",0, 0, OPT_NONPERIODIC },
            { "boxunit",    1, 0, OPT_BOXUNIT },
            { "grid-size",  1, 0, OPT_GRID_SIZE },
            { "grid-slop",  1, 0, OPT_GRID_SLOP },
            { "output",     1, 0, OPT_GRID_FILE },
            { "boundary",   1, 0, OPT_GRID_BOUND },
            { "center",     1, 0, OPT_CENTER },
            { "centre",     1, 0, OPT_CENTER },
            { "recenter",   1, 0, OPT_CENTER },
            { "recentre",   1, 0, OPT_CENTER },
            { "double",     0, 0, OPT_DOUBLE },
	    { NULL,         0, 0, 0 },
        };

        c = getopt_long( argc, argv, "hp:nm:g:o:b:c:s:d",
                         long_options, &option_index );
        if ( c == -1 ) break;

        switch(c) {
        case OPT_HELP:
            bHelp = true;
            break;
	case OPT_PERIODIC:
	    float f;
	    if ( optarg )
		f = atof(optarg) * 0.5;
	    else
		f = 0.5;
	    info.setLimits( -f, f );
	    bPeriodSet = true;
	    break;
	case OPT_NONPERIODIC:
	    info.setLimits( -std::numeric_limits<float>::max(),
			    std::numeric_limits<float>::max() );
	    bPeriodSet = true;
	    break;
	case OPT_MARK:
	    assert( optarg != 0 );
	    nameMark = optarg;
	    break;
	case OPT_BOXUNIT:
	    assert( optarg != 0 );
	    BoxUnit = atof(optarg);
	    break;
	case OPT_GRID_SIZE:
	    assert( optarg != 0 );
	    GridSize = atoi(optarg);
	    break;
	case OPT_GRID_BOUND:
	    assert( optarg != 0 );
	    GridBound = atoi(optarg);
	    break;
	case OPT_GRID_FILE:
	    assert( optarg != 0 );
	    GridFile = optarg;
	    break;
	case OPT_GRID_SLOP:
	    assert( optarg != 0 );
	    GridSlop = atof( optarg );
	    break;
	case OPT_CENTER:
	    assert( optarg != 0 );

            p = strtok(optarg,",");
	    if ( p == NULL ) BadCenter();
            cx = atof(p);
            p = strtok(NULL,",");
	    if ( p == NULL ) BadCenter();
            cy = atof(p);
            p = strtok(NULL,",");
	    if ( p == NULL ) BadCenter();
            cz = atof(p);
            p = strtok(NULL,",");
	    if ( p != NULL ) BadCenter();
	    break;
	case OPT_DOUBLE:
	    Adapter = "double";
	    break;
        default:
            bError = true;
            break;
        }


    }

    //! If an error occurred, or help was explicitly requested, print the
    //! usage and exit.
    if ( bError || bHelp ) {
        Usage(argv[0]);
        exit(1);
    }

    if ( optind < argc ) {
	nameTipsy = argv[optind++];
    }
    else {
	std::cerr << "Missing Tipsy input file" << std::endl;
	Usage(argv[0]);
	exit(1);
    }

    if ( optind < argc ) {
	boxtype = argv[optind++];
    }

    if ( strcmp(boxtype,"all") == 0 ) {
	ip = &info.info;
    }
    else if ( strcmp(boxtype,"dark") == 0 ) {
	ip = &info.dinfo;
    }
    else if ( strcmp(boxtype,"baryon") == 0 ) {
	ip = &info.binfo;
    }
    else if ( strcmp(boxtype,"gas") == 0 ) {
	ip = &info.ginfo;
    }
    else if ( strcmp(boxtype,"star") == 0 ) {
	ip = &info.sinfo;
    }
    else {
	Usage(argv[0]);
	exit(1);
    }

    in.open(nameTipsy.c_str(),Adapter);
    if ( ! in.is_open() ) {
	std::cerr << "Unable to open Tipsy binary " << nameTipsy << std::endl;
        exit(2);
    }

    // If we have been asked to generate a grid map, then set it up
    info.setGrid(GridFile,GridSize,GridBound, GridSlop, cx, cy, cz );

    in >> h;

//number of dark, gas and star particles = 16777216, 0, and 0
    std::cout << "number of dark, gas and star particles = "
	      << h.h_nDark << ", " << h.h_nSph << ", and "
	      << h.h_nStar << std::endl;

    if ( nameMark.empty() ) {
	// By default, don't use a periodic box
	if ( !bPeriodSet )
	    info.setLimits( -std::numeric_limits<float>::max(),
			    std::numeric_limits<float>::max() );
	for( i=0; i<h.h_nBodies; i++ ) {
	    info.ReadNext(in,i);
	}
    }
    else {
	std::ifstream mark(nameMark.c_str());
	uint32_t ng, nd, ns;

	if ( !mark.is_open() ) {
	    std::cerr << "Unable to open file " << nameMark << std::endl;
	    exit(2);
	}
	mark.exceptions( std::fstream::badbit |
			 std::fstream::eofbit |
			 std::fstream::failbit );


	mark >> nd >> ng >> ns;
	if ( h.h_nSph != ng || h.h_nDark != nd || h.h_nStar != ns ) {
	    std::cerr << "Mark file does not match Tipsy file:" << std::endl
		      << "Tipsy: "
		      << h.h_nDark << ", " << h.h_nSph << ", and "
		      << h.h_nStar << std::endl
		      << "Mark:  " << nd << ", " << ng << ", " << ns
		      << std::endl;
	    exit(2);
	}

	while( !mark.eof() ) {
	    try {
		mark >> i;
	    }
	    catch(...) {
		break;
	    }
	    in.seekg( tipsypos(tipsypos::particle,i-1) );
	    info.ReadNext(in, i-1);
	}

	mark.close();


    }
    info.Normalize();
    info.writeGrid(GridFile);

    if ( GridSize ) {
	std::clog << "GRID: marked " << info.countGrid() << " cells" << std::endl;
    }



    float Volume =
	(info.info.X.getMax() - info.info.X.getMin())   * BoxUnit
	* (info.info.Y.getMax() - info.info.Y.getMin()) * BoxUnit
	* (info.info.Z.getMax() - info.info.Z.getMin()) * BoxUnit;
    float Density = info.info.Mass / Volume;


    if ( info.zeros.begin() != info.zeros.end() ) {
	ZeroList::iterator zi;
	unsigned int j;

	std::cout << "Particles with zero mass:" << std::endl;

	zi = info.zeros.begin();
	i = j = *zi;
	while( ++zi != info.zeros.end() ) {
	    if ( *zi != (j+1) ) {
		if ( i == j )
		    std::cout << i << std::endl;
		else
		    std::cout << i << " through " << j 
			      << " (" << j-i+1 << ")" << std::endl;
		i = *zi;
	    }
	    j = *zi;
	}
	if ( i == j )
	    std::cout << i << std::endl;
	else
	    std::cout << i << " through " << j 
		      << " (" << j-i+1 << ")" << std::endl;
    }

    std::cout.flags(std::ios::right|std::ios::scientific);
    std::cout.precision(6);

    std::cout << "mass = " << ip->Mass
	      << ", density = " << Density
	      << ", volume = " << Volume
	      << std::endl;

    std::cout << "center coordinates = {"
	      << ip->X.getCenter() * BoxUnit << " "
	      << ip->Y.getCenter() * BoxUnit << " "
	      << ip->Z.getCenter() * BoxUnit
	      << "}" << std::endl;
    std::cout << "size = {"
	      << ip->X.getSize() * BoxUnit << " "
	      << ip->Y.getSize() * BoxUnit << " "
	      << ip->Z.getSize() * BoxUnit
	      << "}" << std::endl;

    std::cout << "center of mass coordinates = {"
	      << ip->Mx * BoxUnit
	      << " " << ip->My * BoxUnit
	      << " " << ip->Mz * BoxUnit
	      << "}" << std::endl;

    std::cout << "center of mass velocity = {"
	      << ip->Vx * BoxUnit 
	      << " " << ip->Vy * BoxUnit
	      << " " << ip->Vz * BoxUnit
	      << "}" << std::endl;

    std::cout << "angular momentum vector = {"
	      << ip->Angx << " " << ip->Angy << " " << ip->Angz << "}"
	      << std::endl;


    std::cout << "Particle Mass Counts:" << std::endl;

    MassMap::iterator mi;
    for( mi=info.massinfo.begin(); mi != info.massinfo.end(); mi++ ) {
	std::cout << "    " << mi->first << " * " << mi->second 
		  << " = " << mi->first * mi->second << std::endl;
    }

    in.close();
}
