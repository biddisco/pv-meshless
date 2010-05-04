/**
 *  @file
 *  @brief Dump particles in a Tipsy file.
 *  @author Doug Potter
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <getopt.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include "regex.hpp"

#include "ftipsy.hpp"

#define OPT_HELP 'h'
#define OPT_SINGLE '1'
#define OPT_STANDARD 's'
#define OPT_NATIVE 'n'
#define OPT_MARKFILE 'm'
#define OPT_TIPSYOUT 't'
#define OPT_PRECISION 'p'

class range {
public:
    tipsypos::section_type section;
    tipsypos::offset_type  begin;
    tipsypos::offset_type  end;

    range(
	tipsypos::section_type s,
	tipsypos::offset_type  b,
	tipsypos::offset_type  e)
	{ section=s; begin=b; end=e; }

    range() { section=tipsypos::eof; begin=end=0; }

};

class rangelist : public std::vector<range>
{
};


static void Display( TipsyHeader &h ) {
    double red = 1.0 / h.h_time - 1;

    std::cout << "TIME: " << h.h_time << " (redshift " << red << ")"
              << std::endl
              << "#: " << h.h_nBodies 
              << " total (" << h.h_nDims << " dimensions), "
              << h.h_nDark << " dark, "
              << h.h_nSph  << " gas, "
              << h.h_nStar << " star." << std::endl;
}

static void Display( TipsyDarkParticle &p, bool bSplit ) {
    std::cout 
	<< p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
	<< p.vel[0] << " " << p.vel[1] << " " << p.vel[2] << " ";
    if ( bSplit ) std::cout << std::endl;
    std::cout 
	<< p.mass   << " " << p.phi    << " " << p.eps    << std::endl;
}

static void Display( TipsyGasParticle &p, bool bSplit ) {
    std::cout 
	<< p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
	<< p.vel[0] << " " << p.vel[1] << " " << p.vel[2] << " ";
    if ( bSplit ) std::cout << std::endl;
    std::cout
	<< p.mass   << " " << p.phi    << " " << p.hsmooth<< " "
	<< p.metals << " " << p.temp   << " " << p.rho    << std::endl;
}

static void Display( TipsyStarParticle &p, bool bSplit ) {
    std::cout 
	<< p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
	<< p.vel[0] << " " << p.vel[1] << " " << p.vel[2] << " ";
    if ( bSplit) std::cout << std::endl;
    std::cout 
	<< p.mass   << " " << p.phi    << " " << p.eps    << std::endl;
}

static void Usage( const char *name ) {
    std::clog << "Usage: " << name 
	      << " [-h1] <tipsy> [g|d|s]#[-[+]#] ..." 
	      << std::endl
	      << "  -h,--help        Show help (this text)" << std::endl 
	      << "  --standard,--std Tipsy standard file (default)"<<std::endl 
	      << "  --native,--nat   Tipsy native file"<<std::endl 
	      << "  -p,--precision=5 Set output precision"<<std::endl 
	      << "  -1,--single      Display each particle on a single line"
	      << std::endl 
	      << "  --markfile=file  Dump only particles from the mark file"
	      << std::endl 
	      << "  --tipsyout=file  Write a Tipsy file instead"
	      << std::endl 
	      << std::endl
	      << "  You must specify at least one "
	      << "particle to dump.  Particle numbers" << std::endl
	      << "  start at zero.  Specify 'g' for gas, 'd' for dark,"
	      << " or 's' for star" << std::endl
	      << "  particles.  If you omit g, d or s, the "
	      << "absolute particle number is" << std::endl
	      << "  is used.  For example, " << std::endl << std::endl
	      << "    tdump in.std 0           "
	      << "Dump the first particle (particle 0)" << std::endl
	      << "    tdump in.std d0          "
	      << "Dump the first dark particle (particle 0)" << std::endl
	      << "    tdump in.std g0-10       "
	      << "Dump the first 11 gas particles" << std::endl
	      << "    tdump in.std d0 g0 s0    "
	      << "Dump the first dark, gas and star" << std:: endl
	      << "    tdump in.std d100-101    "
	      << "Dump dark particles 100 and 101" << std::endl
	      << "    tdump in.std d100000-+10 "
	      << "Dump dark particles 100,000 to 100,010" << std::endl;
}

int main( int argc, char *argv[] ) {
    bool bHelp  = false;    //!< Help was requested on the command line
    bool bError = false;    //!< An error occurred on the command line
    bool bSplit = true;     //!< Split lines
    int  iPrecision=5;
    std::string ftype="standard";
    std::string filename;
    std::string markname;
    std::string tipsyname;

    ifTipsy in;
    ofTipsy out;
    rangelist worklist;
    range workunit;

    TipsyHeader       h;
    TipsyGasParticle  g;
    TipsyDarkParticle d;
    TipsyStarParticle s;

    unsigned int i;


    //! Parse command line
    for(;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { "help",        0, 0, OPT_HELP },
            { "single",      0, 0, OPT_SINGLE },
            { "standard",    0, 0, OPT_STANDARD },
            { "std",         0, 0, OPT_STANDARD },
            { "native",      0, 0, OPT_NATIVE },
            { "nat",         0, 0, OPT_NATIVE },
            { "markfile",    1, 0, OPT_MARKFILE },
            { "tipsyout",    1, 0, OPT_TIPSYOUT },
            { "precision",   1, 0, OPT_PRECISION },
            { 0,             0, 0, 0 }
        };

        c = getopt_long( argc, argv, "h1p:",
                         long_options, &option_index );
        if ( c == -1 ) break;
        switch(c) {
        case OPT_HELP:
            bHelp = true;
            break;
	case OPT_SINGLE:
	    bSplit = false;
	    break;
	case OPT_STANDARD:
	    ftype = "standard";
	    break;
	case OPT_NATIVE:
	    ftype = "native";
	    break;
	case OPT_MARKFILE:
	    assert(optarg!=NULL);
	    markname = optarg;
	    break;
	case OPT_TIPSYOUT:
	    assert(optarg!=NULL);
	    tipsyname = optarg;
	    break;
	case OPT_PRECISION:
	    assert(optarg!=NULL);
	    iPrecision = atoi(optarg);
	    break;
        default:
            bError = true;
            break;
        }
    }

    if ( bHelp || bError ) {
	Usage(argv[0]);
	exit(1);
    }

    if ( optind < argc ) {
	std::string::size_type sep;
        filename = argv[optind++];
	if ( (sep=filename.find(':')) != std::string::npos ) {
	    ftype = filename.substr(0,sep);
	    filename = filename.substr(sep+1);
	}
    }
    else {
        std::cerr << "Missing tipsy input file" << std::endl
		  << "Try " << argv[0] << " --help" << std::endl;
        exit(2);
    }

    in.open(filename.c_str(),ftype.c_str());
    if ( ! in.is_open() ) {
	std::cerr << "Unable to open Tipsy binary " << filename << std::endl;
        exit(2);
    }

    in >> h;

    if ( h.h_nBodies != h.h_nSph + h.h_nDark + h.h_nStar 
	|| h.h_nDims < 1 || h.h_nDims > 3 ) {
	std::cerr << filename << " has crazy dimensions" << std::endl
		  << "(have you tried --std or --nat?)" << std::endl;
	exit(2);
    }


    Regex re("^([gsd])?([[:digit:]]+)(-(\\+?)([[:digit:]]+))?$");

    while( optind < argc ) {
	if ( re.Search(argv[optind]) >= 0 ) {
	    tipsypos::section_type p;
	    tipsypos::offset_type  b;
	    tipsypos::offset_type  e;

	    if ( re[0].size() == 0 ) p = tipsypos::particle;
	    else switch( re[0][0] ) {
	    case 'g': p = tipsypos::gas;  break;
	    case 'd': p = tipsypos::dark; break;
	    case 's': p = tipsypos::star; break;
	    default: assert( 0==1 /* Invalid particle type */ );
	    }

	    b = atoi(re[1].c_str());
	    if ( re[4].size() == 0 ) e = b;
	    else {
		e = atoi( re[4].c_str() );
		if ( re[3].size() != 0 ) e += b;
	    }
	    switch(p) {
	    case tipsypos::gas:
		bError = b>=h.h_nSph || e>=h.h_nSph;
		break;
	    case tipsypos::dark:
		bError = b>=h.h_nDark || e>=h.h_nDark;
		break;
	    case tipsypos::star:
		bError = b>=h.h_nStar || e>=h.h_nStar;
		break;
	    case tipsypos::particle:
		bError = b>=h.h_nBodies || e>=h.h_nBodies;
		break;
	    default:
		assert( 0==1 /* Invalid file section */ );
	    }
	    if ( bError ) {
		std::cerr << "Range " << argv[optind]
			  << " is outside the file" << std::endl;
		exit(2);
	    }

	    worklist.push_back(range(p,b,e));
	}
	else {
	    std::cerr << "Invalid particle selection: "
		      << argv[optind] << std::endl
		      << "Try " << argv[0] << " --help" << std::endl;
	    exit(2);
	}
	optind++;
    }


    if ( !markname.empty() ) {
	std::ifstream mf(markname.c_str());
	uint32_t nd, ng, ns;
	if ( !mf.is_open() ) {
	    std::cerr << "Unable to open " << markname << std::endl;
	    exit(2);
	}
	mf >> nd >> ng >> ns;
	if ( nd != h.h_nDark || ng != h.h_nSph || ns != h.h_nStar ) {
	    std::cerr << markname << " does not match the tipsy file"
		      << std::endl;
	    exit(2);
	}
	while( !mf.eof() ) {
	    mf >> i;
	    worklist.push_back(range(tipsypos::particle,i-1,i-1));
	}
    }

    if ( worklist.size() == 0 ) {
	worklist.push_back(range(tipsypos::particle,0,h.h_nBodies-1));
    }

    if ( !tipsyname.empty() ) {
	struct stat sbuf;
	int rc = stat( tipsyname.c_str(), &sbuf );

	if ( rc != -1 || errno!=ENOENT ) {
	    std::cerr << "Output file exists, delete it first" << std::endl;
	    std::cerr << rc << std::endl;
	    exit(3);
	}

	out.open( tipsyname.c_str() );
	if ( !out.is_open() ) {
	    std::cerr << "Unable to create " << tipsyname << std::endl;
	    exit(2);
	}
    }

    if ( out.is_open() ) {
	out << h;
    }
    else {
	Display(h);
	std::cout << "     x            y             z      "
		  << "     Vx           Vy            Vz     ";
	if ( bSplit ) std::cout << std::endl;
	std::cout << "    Mass         phi       eps/Smooth  " 
		  << "   Metals        Temp           Rho"<< std::endl;
	std::cout.flags(std::ios::right|std::ios::scientific|std::ios::showpos);
	std::cout.precision(iPrecision);
    }

    h.h_nSph = h.h_nDark = h.h_nStar = h.h_nBodies = 0;

    for( i=0; i<worklist.size(); i++ ) {
	tipsypos::offset_type o;
	in.seekg( tipsypos(worklist[i].section, worklist[i].begin) );

	for( o=worklist[i].begin; o<=worklist[i].end; o++ ) {
	    switch( in.tellg().section() ) {
	    case tipsypos::gas:
		in >> g;
		h.h_nSph++;
		if ( out.is_open() )
		    out << g;
		else
		    Display(g,bSplit);
		break;
	    case tipsypos::dark:
		in >> d;
		h.h_nDark++;
		if ( out.is_open() )
		    out << d;
		else
		    Display(d,bSplit);
		break;
	    case tipsypos::star:
		in >> s;
		h.h_nStar++;
		if ( out.is_open() )
		    out << s;
		else
		    Display(s,bSplit);
		break;
	    default:
		break;
	    }
	}
    }

    if ( out.is_open() ) {
	// Rewrite the header
	h.h_nBodies = h.h_nSph + h.h_nDark + h.h_nStar;
	out.seekp( tipsypos( tipsypos::header, 0 ) );
	out << h;
	out.close();
    }

    in.close();
}
